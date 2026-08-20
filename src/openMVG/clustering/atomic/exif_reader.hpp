#ifndef __EXIF_READER_HPP__
#define __EXIF_READER_HPP__

#pragma once

#include "types.hpp"

#include <map>
#include <string>

namespace AtomicClustering {

    // Reads the GPS coordinate embedded in one image's EXIF data, if any.
    // Returns a GpsCoord with valid=false if the file has no EXIF GPS IFD,
    // isn't a format easyexif can parse, or is missing/unreadable.
    GpsCoord read_gps(const std::string &image_path);

    // Recursively walks `folder` for supported image extensions and reads
    // each one's GPS coordinate via read_gps(), keeping only the geotagged
    // ones. Keys are the full filesystem path of each image, matching the
    // ImageId convention expected by the rest of the atomic clustering
    // pipeline (see distance.hpp / visual_features.hpp, which open the id
    // directly as a file path).
    std::map<ImageId, GpsCoord> load_gps_points(const std::string &folder);

} // namespace AtomicClustering

#endif /* __EXIF_READER_HPP__ */
