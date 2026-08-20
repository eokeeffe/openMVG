#ifndef __JSON_EXPORT_HPP__
#define __JSON_EXPORT_HPP__

#pragma once

#include "types.hpp"

#include <map>
#include <string>
#include <vector>

namespace AtomicClustering {

    // Writes `clusters` as a JSON array to `output_path`, one object per
    // cluster with its id and member images' path + GPS coordinate.
    void export_json(const std::vector<Cluster> &clusters, const std::map<ImageId, GpsCoord> &points,
                      const std::string &output_path);

    // Prints a human-readable report to stdout: per-cluster member list and
    // GPS centroid, plus overlap/unclustered counts across the whole run.
    void summarize(const std::vector<Cluster> &clusters, const std::map<ImageId, GpsCoord> &points);

} // namespace AtomicClustering

#endif /* __JSON_EXPORT_HPP__ */
