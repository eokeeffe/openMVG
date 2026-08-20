#ifndef __DSITANCE_HPP__
#define __DSITANCE_HPP__

#pragma once

#include "types.hpp"

#include <functional>
#include <map>
#include <vector>

namespace AtomicClustering {

    // dist(a, b) -> a non-negative distance between two images. The atomic
    // clustering pipeline (see atomic_clustering.hpp) is written entirely
    // against this interface, so it doesn't matter whether the underlying
    // metric is GPS distance, visual dissimilarity, or a blend of both.
    using DistanceFn = std::function<double(const ImageId &, const ImageId &)>;

    // Great-circle distance between two GPS coordinates, in meters.
    double haversine_m(const GpsCoord &a, const GpsCoord &b);

    // dist(a, b) = haversine distance between two images' GPS coordinates.
    DistanceFn make_gps_distance_fn(const std::map<ImageId, GpsCoord> &points);

    // Precomputes the full pairwise distance matrix for `raw_fn` once and
    // returns a fast memoized lookup function. O(n^2) time and memory -- fine
    // for up to a few thousand images. For much larger datasets, pre-bucket
    // (e.g. by geohash/grid cell) and cluster per-bucket instead.
    DistanceFn build_distance_matrix(const std::vector<ImageId> &ids, const DistanceFn &raw_fn);

} // namespace AtomicClustering

#endif /* __DSITANCE_HPP__ */
