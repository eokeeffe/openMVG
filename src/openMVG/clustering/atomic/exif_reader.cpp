#include "exif_reader.hpp"

#include "exif.h" // easyexif -- https://github.com/mayanklahiri/easyexif

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <vector>

// The rest of OpenMVG targets C++11 and uses stlplus for filesystem access
// rather than C++17's <filesystem>, so this does the same instead of
// requiring every consumer of this library to opt into a newer standard.
namespace AtomicClustering {

namespace {

bool has_image_extension(const std::string &path) {
    static const std::vector<std::string> exts = {"jpg", "jpeg", "tif", "tiff", "png", "heic"};
    std::string ext = stlplus::extension_part(path);
    std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c) { return std::tolower(c); });
    return std::find(exts.begin(), exts.end(), ext) != exts.end();
}

// Recursively collects every file under `folder` with a supported image
// extension, appending their full paths to `out`.
void collect_image_files(const std::string &folder, std::vector<std::string> &out) {
    for (const auto &filename : stlplus::folder_files(folder)) {
        const std::string full_path = stlplus::create_filespec(folder, filename);
        if (has_image_extension(full_path)) {
            out.push_back(full_path);
        }
    }
    for (const auto &subfolder : stlplus::folder_subdirectories(folder)) {
        collect_image_files(stlplus::create_filespec(folder, subfolder), out);
    }
}

} // namespace

GpsCoord read_gps(const std::string &image_path) {
    GpsCoord coord;

    std::ifstream file(image_path, std::ios::binary | std::ios::ate);
    if (!file) {
        return coord;
    }

    std::streamsize size = file.tellg();
    if (size <= 0) {
        return coord;
    }
    file.seekg(0, std::ios::beg);

    std::vector<unsigned char> buf(static_cast<size_t>(size));
    if (!file.read(reinterpret_cast<char *>(buf.data()), size)) {
        return coord;
    }

    easyexif::EXIFInfo info;
    int parse_code = info.parseFrom(buf.data(), static_cast<unsigned int>(buf.size()));
    if (parse_code != 0) {
        return coord; // not a JPEG, no EXIF segment, or a parse error
    }

    // easyexif leaves GeoLocation at (0, 0) when there is no GPS IFD at
    // all. A genuine (0, 0) GPS reading is astronomically unlikely for a
    // real photo, so treat it the same as "no GPS data".
    if (info.GeoLocation.Latitude == 0.0 && info.GeoLocation.Longitude == 0.0) {
        return coord;
    }

    coord.lat = info.GeoLocation.Latitude;
    coord.lon = info.GeoLocation.Longitude;
    coord.valid = true;
    return coord;
}

std::map<ImageId, GpsCoord> load_gps_points(const std::string &folder) {
    std::map<ImageId, GpsCoord> points;

    std::vector<std::string> files;
    collect_image_files(folder, files);
    std::sort(files.begin(), files.end());

    std::cout << "Found " << files.size() << " candidate image(s) in " << folder << "\n";

    for (const auto &f : files) {
        GpsCoord coord = read_gps(f);
        if (coord.valid) {
            points[f] = coord;
        } else {
            std::cerr << "  [skip] no GPS data: " << stlplus::filename_part(f) << "\n";
        }
    }

    std::cout << "Extracted GPS coordinates for " << points.size() << " image(s)\n";
    return points;
}

} // namespace AtomicClustering
