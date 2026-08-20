#include "distance.hpp"

#include <cmath>
#include <iostream>
#include <memory>
#include <unordered_map>

namespace AtomicClustering {

namespace {
constexpr double EARTH_RADIUS_M = 6371000.0;
constexpr double PI = 3.14159265358979323846;
} // namespace

double haversine_m(const GpsCoord &a, const GpsCoord &b) {
    double lat1 = a.lat * PI / 180.0;
    double lat2 = b.lat * PI / 180.0;
    double dlat = (b.lat - a.lat) * PI / 180.0;
    double dlon = (b.lon - a.lon) * PI / 180.0;

    double h = std::sin(dlat / 2) * std::sin(dlat / 2) +
               std::cos(lat1) * std::cos(lat2) * std::sin(dlon / 2) * std::sin(dlon / 2);
    return 2 * EARTH_RADIUS_M * std::asin(std::sqrt(h));
}

DistanceFn make_gps_distance_fn(const std::map<ImageId, GpsCoord> &points) {
    // Store a copy behind a shared_ptr so the returned lambda doesn't hold
    // a dangling reference if the caller's `points` ever goes out of
    // scope before the lambda does.
    auto points_ptr = std::make_shared<std::map<ImageId, GpsCoord>>(points);
    return [points_ptr](const ImageId &a, const ImageId &b) -> double {
        return haversine_m(points_ptr->at(a), points_ptr->at(b));
    };
}

DistanceFn build_distance_matrix(const std::vector<ImageId> &ids, const DistanceFn &raw_fn) {
    auto idx = std::make_shared<std::unordered_map<ImageId, size_t>>();
    for (size_t i = 0; i < ids.size(); ++i) {
        (*idx)[ids[i]] = i;
    }

    size_t n = ids.size();
    auto mat = std::make_shared<std::vector<std::vector<double>>>(n, std::vector<double>(n, 0.0));

    size_t total_pairs = (n > 0) ? n * (n - 1) / 2 : 0;
    size_t done = 0;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            double d = raw_fn(ids[i], ids[j]);
            (*mat)[i][j] = d;
            (*mat)[j][i] = d;
            ++done;
        }
        if (total_pairs > 20000 && (done % 20000) < n) {
            std::cerr << "  ... computed " << done << "/" << total_pairs << " pairwise distances\n";
        }
    }

    return [idx, mat](const ImageId &a, const ImageId &b) -> double {
        if (a == b) {
            return 0.0;
        }
        return (*mat)[idx->at(a)][idx->at(b)];
    };
}

} // namespace AtomicClustering
