#ifndef __VISUAL_FEATURES_HPP__
#define __VISUAL_FEATURES_HPP__

#pragma once

#include "distance.hpp"
#include "types.hpp"

#include <opencv2/core.hpp>

#include <map>
#include <string>
#include <vector>

namespace AtomicClustering {
    // A lightweight stand-in for the paper's precomputed hierarchical
    // vocabulary tree: rather than a tree, we build a flat "bag of visual
    // words" vocabulary via k-means over ORB descriptors sampled from the
    // dataset itself, then weight each image's word-count histogram by IDF --
    // the same weighting scheme as the paper's Equations 1-4
    // (q_i = n_i * w_i, w_i = ln(N / N_i)). This is a simplification for
    // practicality: it captures the same "how visually similar are these two
    // images" signal without needing an externally-supplied vocabulary tree.

    // Extracts ORB descriptors for one image, downscaled so its longer side
    // is at most resize_max_dim pixels (mirrors the paper's own resize-before-
    // feature-extraction step). Returns an empty Mat if the image can't be
    // read or no features are found.
    cv::Mat extract_orb_descriptors(const std::string &image_path, int resize_max_dim, int n_features);

    // K-means over sampled ORB descriptors -> a flat visual-word vocabulary.
    // Returns a [vocab_size x 32] CV_32F matrix of cluster centers.
    cv::Mat build_visual_vocabulary(const std::map<ImageId, cv::Mat> &descriptors, int vocab_size,
                                    int sample_per_image);

    // For each image: a nearest-vocabulary-word histogram, weighted by IDF
    // across the whole dataset.
    std::map<ImageId, std::vector<double>> compute_bow_vectors(
        const std::map<ImageId, cv::Mat> &descriptors, const cv::Mat &vocabulary);

    // Full pipeline: ORB extraction -> vocabulary -> IDF-weighted BoW vectors.
    // Throws std::runtime_error if no image yields usable ORB features.
    std::map<ImageId, std::vector<double>> extract_visual_vectors(const std::vector<ImageId> &image_ids,
                                                                    int resize_max_dim, int n_features,
                                                                    int vocab_size);

    // dist(a, b) = normalized L2 distance between two images' IDF-weighted BoW
    // vectors -- the paper's Equation 4: || q/|q| - d/|d| ||. Range [0, sqrt(2)].
    DistanceFn make_visual_distance_fn(const std::map<ImageId, std::vector<double>> &weighted_vectors);

    // dist(a, b) = alpha * normalized_gps_distance + (1 - alpha) * normalized_visual_distance.
    // alpha=1.0 is equivalent to pure GPS distance; alpha=0.0 is equivalent to
    // pure visual distance. GPS distance is normalized by the dataset's own
    // (sampled) max pairwise distance; visual distance by its theoretical max
    // of sqrt(2).
    DistanceFn make_hybrid_distance_fn(const std::map<ImageId, GpsCoord> &points,
                                        const std::map<ImageId, std::vector<double>> &weighted_vectors,
                                        double alpha);

} // namespace AtomicClustering

#endif /* __VISUAL_FEATURES_HPP__ */
