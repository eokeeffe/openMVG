#ifndef __TYPES_HPP__
#define __TYPES_HPP__

#pragma once

#include <string>
#include <vector>

namespace AtomicClustering {

    // GPS coordinate extracted from an image's EXIF data.
    struct GpsCoord {
        double lat = 0.0;
        double lon = 0.0;
        bool valid = false;

        // Default member initializers make this a non-aggregate under the
        // C++11 rules this project builds with, so brace-init like
        // GpsCoord{lat, lon, true} needs this constructor spelled out.
        GpsCoord() = default;
        GpsCoord(double lat_, double lon_, bool valid_) : lat(lat_), lon(lon_), valid(valid_) {}
    };

    // An image is identified by its filesystem path (must be unique per image).
    using ImageId = std::string;

    // One cluster: an unordered set of image ids, stored sorted for determinism.
    struct Cluster {
        std::vector<ImageId> members;
    };

} // namespace AtomicClustering

#endif /* __TYPES_HPP__ */
