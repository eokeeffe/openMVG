#include "visual_features.hpp"

#include <opencv2/features2d.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>

namespace AtomicClustering {

cv::Mat extract_orb_descriptors(const std::string &image_path, int resize_max_dim, int n_features) {
    cv::Mat img = cv::imread(image_path, cv::IMREAD_GRAYSCALE);
    if (img.empty()) {
        return cv::Mat();
    }

    int h = img.rows, w = img.cols;
    double scale = static_cast<double>(resize_max_dim) / static_cast<double>(std::max(h, w));
    if (scale < 1.0) {
        cv::resize(img, img,
                   cv::Size(std::max(1, static_cast<int>(w * scale)), std::max(1, static_cast<int>(h * scale))));
    }

    cv::Ptr<cv::ORB> orb = cv::ORB::create(n_features);
    std::vector<cv::KeyPoint> keypoints;
    cv::Mat descriptors;
    orb->detectAndCompute(img, cv::noArray(), keypoints, descriptors);
    return descriptors; // CV_8U, [n x 32], empty if none found
}

cv::Mat build_visual_vocabulary(const std::map<ImageId, cv::Mat> &descriptors, int vocab_size,
                                 int sample_per_image) {
    std::mt19937 rng(0);
    std::vector<cv::Mat> sampled;

    for (const auto &kv : descriptors) {
        const cv::Mat &desc = kv.second;
        if (desc.empty()) {
            continue;
        }

        cv::Mat d = desc;
        if (d.rows > sample_per_image) {
            std::vector<int> idxs(d.rows);
            std::iota(idxs.begin(), idxs.end(), 0);
            std::shuffle(idxs.begin(), idxs.end(), rng);
            idxs.resize(sample_per_image);

            cv::Mat picked(sample_per_image, d.cols, d.type());
            for (int i = 0; i < sample_per_image; ++i) {
                d.row(idxs[i]).copyTo(picked.row(i));
            }
            d = picked;
        }

        cv::Mat d_float;
        d.convertTo(d_float, CV_32F);
        sampled.push_back(d_float);
    }

    if (sampled.empty()) {
        throw std::runtime_error("No ORB descriptors could be extracted from any image.");
    }

    cv::Mat all_desc;
    cv::vconcat(sampled, all_desc);

    int actual_vocab_size = std::max(2, std::min(vocab_size, all_desc.rows));

    cv::Mat labels, centers;
    cv::TermCriteria criteria(cv::TermCriteria::EPS + cv::TermCriteria::MAX_ITER, 20, 1.0);
    cv::kmeans(all_desc, actual_vocab_size, labels, criteria, 3, cv::KMEANS_PP_CENTERS, centers);

    return centers; // CV_32F, [actual_vocab_size x 32]
}

std::map<ImageId, std::vector<double>> compute_bow_vectors(const std::map<ImageId, cv::Mat> &descriptors,
                                                             const cv::Mat &vocabulary) {
    int n_words = vocabulary.rows;
    std::map<ImageId, std::vector<double>> raw;

    for (const auto &kv : descriptors) {
        const ImageId &id = kv.first;
        const cv::Mat &desc = kv.second;

        std::vector<double> vec(static_cast<size_t>(n_words), 0.0);
        if (!desc.empty()) {
            cv::Mat desc_f;
            desc.convertTo(desc_f, CV_32F);

            for (int r = 0; r < desc_f.rows; ++r) {
                double best_d = std::numeric_limits<double>::max();
                int best_w = 0;
                for (int w = 0; w < n_words; ++w) {
                    double d = cv::norm(desc_f.row(r), vocabulary.row(w), cv::NORM_L2SQR);
                    if (d < best_d) {
                        best_d = d;
                        best_w = w;
                    }
                }
                vec[static_cast<size_t>(best_w)] += 1.0;
            }
        }
        raw[id] = vec;
    }

    auto n_images = static_cast<double>(raw.size());
    std::vector<double> doc_freq(static_cast<size_t>(n_words), 0.0);
    for (const auto &kv : raw) {
        for (int w = 0; w < n_words; ++w) {
            if (kv.second[static_cast<size_t>(w)] > 0) {
                doc_freq[static_cast<size_t>(w)] += 1.0;
            }
        }
    }

    std::vector<double> idf(static_cast<size_t>(n_words));
    for (int w = 0; w < n_words; ++w) {
        idf[static_cast<size_t>(w)] = std::log(n_images / std::max(1.0, doc_freq[static_cast<size_t>(w)]));
    }

    std::map<ImageId, std::vector<double>> weighted;
    for (const auto &kv : raw) {
        std::vector<double> wv(static_cast<size_t>(n_words));
        for (int w = 0; w < n_words; ++w) {
            wv[static_cast<size_t>(w)] = kv.second[static_cast<size_t>(w)] * idf[static_cast<size_t>(w)];
        }
        weighted[kv.first] = wv;
    }
    return weighted;
}

std::map<ImageId, std::vector<double>> extract_visual_vectors(const std::vector<ImageId> &image_ids,
                                                                int resize_max_dim, int n_features,
                                                                int vocab_size) {
    std::cout << "Extracting ORB features from " << image_ids.size() << " image(s)...\n";

    std::map<ImageId, cv::Mat> descriptor_map;
    int n_ok = 0;
    for (const auto &id : image_ids) {
        cv::Mat desc = extract_orb_descriptors(id, resize_max_dim, n_features);
        if (!desc.empty()) {
            ++n_ok;
        }
        descriptor_map[id] = desc;
    }
    std::cout << "  " << n_ok << "/" << image_ids.size() << " image(s) yielded usable ORB features\n";
    if (n_ok == 0) {
        throw std::runtime_error(
            "No ORB features could be extracted from any image; cannot build a visual vocabulary.");
    }

    std::cout << "Building visual vocabulary (target size " << vocab_size << ")...\n";
    cv::Mat vocabulary = build_visual_vocabulary(descriptor_map, vocab_size, /*sample_per_image=*/100);
    std::cout << "  vocabulary size: " << vocabulary.rows << " visual words\n";

    std::cout << "Computing IDF-weighted bag-of-visual-words vectors...\n";
    return compute_bow_vectors(descriptor_map, vocabulary);
}

DistanceFn make_visual_distance_fn(const std::map<ImageId, std::vector<double>> &weighted_vectors) {
    auto norms = std::make_shared<std::map<ImageId, std::vector<double>>>();

    for (const auto &kv : weighted_vectors) {
        const std::vector<double> &vec = kv.second;
        double norm = 0.0;
        for (double v : vec) {
            norm += v * v;
        }
        norm = std::sqrt(norm);

        std::vector<double> normalized(vec.size());
        if (norm > 0.0) {
            for (size_t i = 0; i < vec.size(); ++i) {
                normalized[i] = vec[i] / norm;
            }
        } else {
            normalized = vec;
        }
        (*norms)[kv.first] = normalized;
    }

    return [norms](const ImageId &a, const ImageId &b) -> double {
        const std::vector<double> &va = norms->at(a);
        const std::vector<double> &vb = norms->at(b);
        double sum = 0.0;
        for (size_t i = 0; i < va.size(); ++i) {
            double d = va[i] - vb[i];
            sum += d * d;
        }
        return std::sqrt(sum);
    };
}

DistanceFn make_hybrid_distance_fn(const std::map<ImageId, GpsCoord> &points,
                                    const std::map<ImageId, std::vector<double>> &weighted_vectors,
                                    double alpha) {
    auto gps_fn = std::make_shared<DistanceFn>(make_gps_distance_fn(points));
    auto visual_fn = std::make_shared<DistanceFn>(make_visual_distance_fn(weighted_vectors));

    std::vector<ImageId> ids;
    ids.reserve(points.size());
    for (const auto &kv : points) {
        ids.push_back(kv.first);
    }

    std::vector<ImageId> sample = ids;
    if (sample.size() > 200) {
        std::mt19937 rng(42);
        std::shuffle(sample.begin(), sample.end(), rng);
        sample.resize(200);
    }

    double gps_scale = 1.0;
    for (size_t i = 0; i < sample.size(); ++i) {
        for (size_t j = i + 1; j < sample.size(); ++j) {
            gps_scale = std::max(gps_scale, (*gps_fn)(sample[i], sample[j]));
        }
    }

    double visual_scale = std::sqrt(2.0);

    return [gps_fn, visual_fn, gps_scale, visual_scale, alpha](const ImageId &a, const ImageId &b) -> double {
        double norm_gps = (*gps_fn)(a, b) / gps_scale;
        double norm_visual = (*visual_fn)(a, b) / visual_scale;
        return alpha * norm_gps + (1.0 - alpha) * norm_visual;
    };
}

} // namespace AtomicClustering
