// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/cameras.hpp"
#include "openMVG/exif/exif_IO_EasyExif.hpp"
#include "openMVG/exif/sensor_width_database/ParseDatabase.hpp"
#include "openMVG/geodesy/geodesy.hpp"
#include "openMVG/image/image_io.hpp"
#include "openMVG/numeric/eigen_alias_definition.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_data_utils.hpp"
#include "openMVG/sfm/sfm_view.hpp"
#include "openMVG/sfm/sfm_view_priors.hpp"
#include "openMVG/system/logger.hpp"
#include "openMVG/system/loggerprogress.hpp"
#include "openMVG/types.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include "openMVG/clustering/atomic/atomic_clustering.hpp"
#include "openMVG/clustering/atomic/distance.hpp"
#include "openMVG/clustering/atomic/types.hpp"
#ifdef ATOMIC_CLUSTERING_HAS_VISUAL
#include "openMVG/clustering/atomic/visual_features.hpp"
#endif

#include <algorithm>
#include <exception>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::exif;
using namespace openMVG::geodesy;
using namespace openMVG::image;
using namespace openMVG::sfm;

// -- Small option parsers, mirrored from main_SfMInit_ImageListing.cpp --

/// Check that Kmatrix is a string like "f;0;ppx;0;f;ppy;0;0;1"
bool checkIntrinsicStringValidity(const std::string &Kmatrix, double &focal, double &ppx, double &ppy)
{
  std::vector<std::string> vec_str;
  stl::split(Kmatrix, ';', vec_str);
  if (vec_str.size() != 9)
  {
    OPENMVG_LOG_ERROR << "\n Missing ';' character";
    return false;
  }
  for (size_t i = 0; i < vec_str.size(); ++i)
  {
    double readvalue = 0.0;
    std::stringstream ss;
    ss.str(vec_str[i]);
    if (!(ss >> readvalue))
    {
      OPENMVG_LOG_ERROR << "\n Used an invalid not a number character";
      return false;
    }
    if (i == 0) focal = readvalue;
    if (i == 2) ppx = readvalue;
    if (i == 5) ppy = readvalue;
  }
  return true;
}

/// Check string of prior weights "x;y;z"
std::pair<bool, Vec3> checkPriorWeightsString(const std::string &sWeights)
{
  std::pair<bool, Vec3> val(true, Vec3::Zero());
  std::vector<std::string> vec_str;
  stl::split(sWeights, ';', vec_str);
  if (vec_str.size() != 3)
  {
    OPENMVG_LOG_ERROR << "Missing ';' character in prior weights";
    val.first = false;
  }
  for (size_t i = 0; i < vec_str.size(); ++i)
  {
    double readvalue = 0.0;
    std::stringstream ss;
    ss.str(vec_str[i]);
    if (!(ss >> readvalue))
    {
      OPENMVG_LOG_ERROR << "Used an invalid not a number character in local frame origin";
      val.first = false;
    }
    val.second[i] = readvalue;
  }
  return val;
}

/// Escape a string for embedding as a JS/JSON double-quoted string literal.
std::string JsEscape(const std::string &s)
{
  std::string out;
  out.reserve(s.size());
  for (char c : s)
  {
    switch (c)
    {
      case '"': out += "\\\""; break;
      case '\\': out += "\\\\"; break;
      case '\n': out += "\\n"; break;
      default: out += c;
    }
  }
  return out;
}

/// Fixed, visually-distinct palette; cycled by cluster index.
const std::vector<std::string> &ClusterPalette()
{
  static const std::vector<std::string> palette = {
    "#e6194b", "#3cb44b", "#4363d8", "#f58231", "#911eb4",
    "#46f0f0", "#f032e6", "#bcf60c", "#fabebe", "#008080",
    "#e6beff", "#9a6324", "#800000", "#808000", "#000075", "#a9a9a9"
  };
  return palette;
}

/// Export the subset of `sfm_data` (views + their intrinsics) referenced by
/// `members` (full image paths) into its own SfM_Data file. There is no
/// reconstruction yet at this stage, so unlike main_ComputeClusters.cpp's
/// exportData() there is no landmark/structure copying to do.
bool ExportClusterSfMData
(
  const SfM_Data &sfm_data,
  const std::string &outFilename,
  const std::vector<AtomicClustering::ImageId> &members,
  const std::map<AtomicClustering::ImageId, IndexT> &pathToViewId
)
{
  SfM_Data cl_sfm_data;
  cl_sfm_data.s_root_path = sfm_data.s_root_path;

  for (const auto &path : members)
  {
    const auto it_id = pathToViewId.find(path);
    if (it_id == pathToViewId.end())
      continue;

    const IndexT view_id = it_id->second;
    const auto it_view = sfm_data.views.find(view_id);
    if (it_view == sfm_data.views.end())
      continue;

    cl_sfm_data.views[view_id] = it_view->second;

    const IndexT id_intrinsic = it_view->second->id_intrinsic;
    if (id_intrinsic != UndefinedIndexT && cl_sfm_data.intrinsics.count(id_intrinsic) == 0)
    {
      const auto it_intrinsic = sfm_data.intrinsics.find(id_intrinsic);
      if (it_intrinsic != sfm_data.intrinsics.end())
        cl_sfm_data.intrinsics[id_intrinsic] = it_intrinsic->second;
    }
  }

  return Save(cl_sfm_data, outFilename, ESfM_Data(VIEWS | INTRINSICS));
}

/// Writes a self-contained HTML map of the clusters to `out_html`, one color
/// per cluster. `basemap` == "osm" renders real OpenStreetMap tiles via
/// Leaflet (loaded from a CDN -- needs internet to view); "none" renders an
/// inline SVG scatter plot with no external requests at all (fully offline,
/// no geographic basemap). Only images with a valid GPS coordinate are
/// plotted -- present regardless of clustering method, since GPS is read
/// for every image up front.
bool WriteLeafletMap
(
  const std::string &out_html,
  const std::vector<AtomicClustering::Cluster> &clusters,
  const std::map<AtomicClustering::ImageId, AtomicClustering::GpsCoord> &points,
  const std::string &basemap
)
{
  struct MapPoint { std::string name; double lat, lon; };
  std::vector<std::vector<MapPoint>> per_cluster_points(clusters.size());

  size_t total_mapped = 0;
  for (size_t i = 0; i < clusters.size(); ++i)
  {
    for (const auto &member : clusters[i].members)
    {
      const auto it = points.find(member);
      if (it == points.end() || !it->second.valid)
        continue;
      per_cluster_points[i].push_back({stlplus::filename_part(member), it->second.lat, it->second.lon});
      ++total_mapped;
    }
  }

  if (total_mapped == 0)
  {
    OPENMVG_LOG_WARNING << "No geotagged images among the clustered images; skipping Leaflet map.";
    return false;
  }

  std::ofstream out(out_html);
  if (!out)
  {
    OPENMVG_LOG_ERROR << "Cannot write: " << out_html;
    return false;
  }
  out << std::fixed << std::setprecision(7);

  const std::vector<std::string> &palette = ClusterPalette();

  out << "<!DOCTYPE html>\n<html>\n<head>\n<meta charset=\"utf-8\">\n"
      << "<title>Atomic Clustering - GPS map</title>\n";

  if (basemap == "osm")
  {
    out << "<link rel=\"stylesheet\" href=\"https://unpkg.com/leaflet@1.9.4/dist/leaflet.css\">\n"
        << "<script src=\"https://unpkg.com/leaflet@1.9.4/dist/leaflet.js\"></script>\n";
  }

  out << "<style>\n"
      << "  html, body { margin:0; padding:0; height:100%; font-family: sans-serif; }\n"
      << "  #map { position:absolute; top:0; left:0; right:260px; bottom:0; }\n"
      << "  #legend { position:absolute; top:0; right:0; width:260px; bottom:0; overflow:auto;\n"
      << "            border-left:1px solid #ccc; box-sizing:border-box; padding:10px; }\n"
      << "  .legend-row { display:flex; align-items:center; margin-bottom:6px; }\n"
      << "  .swatch { width:14px; height:14px; border-radius:50%; margin-right:8px; flex-shrink:0; }\n"
      << "</style>\n</head>\n<body>\n"
      << "<div id=\"map\"></div>\n"
      << "<div id=\"legend\"><h3>Clusters</h3><div id=\"legend-rows\"></div></div>\n";

  // Cluster data, shared by both rendering paths below.
  out << "<script>\nconst clusters = [\n";
  for (size_t i = 0; i < per_cluster_points.size(); ++i)
  {
    out << "  {\"id\": " << i << ", \"color\": \"" << palette[i % palette.size()]
        << "\", \"total\": " << clusters[i].members.size() << ", \"points\": [\n";
    const auto &pts = per_cluster_points[i];
    for (size_t j = 0; j < pts.size(); ++j)
    {
      out << "    {\"name\": \"" << JsEscape(pts[j].name) << "\", \"lat\": " << pts[j].lat
          << ", \"lon\": " << pts[j].lon << "}" << (j + 1 < pts.size() ? ",\n" : "\n");
    }
    out << "  ]}" << (i + 1 < per_cluster_points.size() ? ",\n" : "\n");
  }
  out << "];\n";

  // Legend, common to both rendering paths.
  out << R"JS(
const legendRows = document.getElementById('legend-rows');
clusters.forEach(c => {
  const row = document.createElement('div');
  row.className = 'legend-row';
  const mapped = c.points.length;
  row.innerHTML = '<span class="swatch" style="background:' + c.color + '"></span>' +
    'Cluster ' + c.id + ' (' + mapped + '/' + c.total + ' mapped)';
  legendRows.appendChild(row);
});
)JS";

  if (basemap == "osm")
  {
    out << R"JS(
const map = L.map('map');
L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
  maxZoom: 19,
  attribution: '&copy; OpenStreetMap contributors'
}).addTo(map);

const overlays = {};
const allMarkers = [];
clusters.forEach(c => {
  const layer = L.layerGroup();
  c.points.forEach(p => {
    const marker = L.circleMarker([p.lat, p.lon], {
      radius: 6, color: c.color, weight: 1, fillColor: c.color, fillOpacity: 0.85
    }).bindPopup('<b>Cluster ' + c.id + '</b><br>' + p.name + '<br>' + p.lat.toFixed(6) + ', ' + p.lon.toFixed(6));
    marker.addTo(layer);
    allMarkers.push(marker);
  });
  layer.addTo(map);
  overlays['Cluster ' + c.id] = layer;
});
L.control.layers(null, overlays).addTo(map);

if (allMarkers.length > 0) {
  map.fitBounds(L.featureGroup(allMarkers).getBounds().pad(0.1));
} else {
  map.setView([0, 0], 2);
}
)JS";
  }
  else
  {
    // Fully offline fallback: a plain inline-SVG scatter plot, no external
    // requests, no basemap/tiles -- just relative point positions.
    out << R"JS(
let minLat = Infinity, maxLat = -Infinity, minLon = Infinity, maxLon = -Infinity;
clusters.forEach(c => c.points.forEach(p => {
  minLat = Math.min(minLat, p.lat); maxLat = Math.max(maxLat, p.lat);
  minLon = Math.min(minLon, p.lon); maxLon = Math.max(maxLon, p.lon);
}));
const padLat = Math.max((maxLat - minLat) * 0.05, 1e-5);
const padLon = Math.max((maxLon - minLon) * 0.05, 1e-5);
minLat -= padLat; maxLat += padLat; minLon -= padLon; maxLon += padLon;

const mapDiv = document.getElementById('map');
const W = 1000, H = 1000;
const svgNS = 'http://www.w3.org/2000/svg';
const svg = document.createElementNS(svgNS, 'svg');
svg.setAttribute('viewBox', '0 0 ' + W + ' ' + H);
svg.setAttribute('width', '100%');
svg.setAttribute('height', '100%');
svg.style.background = '#f0f0f0';
mapDiv.appendChild(svg);

function project(lat, lon) {
  const x = (lon - minLon) / (maxLon - minLon) * W;
  const y = H - (lat - minLat) / (maxLat - minLat) * H; // north-up
  return [x, y];
}

const groups = {};
clusters.forEach(c => {
  const g = document.createElementNS(svgNS, 'g');
  g.setAttribute('data-cluster', c.id);
  c.points.forEach(p => {
    const [x, y] = project(p.lat, p.lon);
    const circle = document.createElementNS(svgNS, 'circle');
    circle.setAttribute('cx', x);
    circle.setAttribute('cy', y);
    circle.setAttribute('r', 6);
    circle.setAttribute('fill', c.color);
    circle.setAttribute('fill-opacity', '0.85');
    circle.setAttribute('stroke', '#333');
    circle.setAttribute('stroke-width', '0.5');
    const title = document.createElementNS(svgNS, 'title');
    title.textContent = 'Cluster ' + c.id + ': ' + p.name + ' (' + p.lat.toFixed(6) + ', ' + p.lon.toFixed(6) + ')';
    circle.appendChild(title);
    g.appendChild(circle);
  });
  svg.appendChild(g);
  groups[c.id] = g;
});

legendRows.querySelectorAll('.legend-row').forEach((row, i) => {
  row.style.cursor = 'pointer';
  row.addEventListener('click', () => {
    const g = groups[i];
    g.style.display = (g.style.display === 'none') ? '' : 'none';
  });
});
)JS";
  }

  out << "</script>\n</body>\n</html>\n";
  return true;
}

//
// Create a description of an image dataset for OpenMVG, split into
// overlapping clusters by GPS / visual-similarity / a hybrid of both, and
// export a sfm_data file + a colored Leaflet map per run.
//
int main(int argc, char **argv)
{
  CmdLine cmd;

  std::string sImageDir, sfileDatabase = "", sOutputDir = "", sKmatrix;
  std::string sPriorWeights = "1.0;1.0;1.0";
  std::pair<bool, Vec3> prior_w_info(false, Vec3());

  int i_User_camera_model = PINHOLE_CAMERA_RADIAL3;
  bool b_Group_camera_model = true;
  int i_GPS_XYZ_method = 0;
  double focal_pixels = -1.0;

  std::string sMethod = "gps"; // gps | visual | hybrid
  unsigned int min_atom_size = 15;
  unsigned int max_atom_size = 50;
  unsigned int knn_k = 6;
  unsigned int atom_neighbors = 2;
  unsigned int vocab_size = 300;
  unsigned int orb_features = 500;
  unsigned int resize_dim = 640;
  double hybrid_alpha = 0.5;

  std::string sOutputFormat = "BIN"; // BIN | JSON
  std::string sMapFilename = "clusters_map.html";
  std::string sBasemap = "osm"; // osm | none

  cmd.add(make_option('i', sImageDir, "imageDirectory"));
  cmd.add(make_option('d', sfileDatabase, "sensorWidthDatabase"));
  cmd.add(make_option('o', sOutputDir, "outputDirectory"));
  cmd.add(make_option('f', focal_pixels, "focal"));
  cmd.add(make_option('k', sKmatrix, "intrinsics"));
  cmd.add(make_option('c', i_User_camera_model, "camera_model"));
  cmd.add(make_option('g', b_Group_camera_model, "group_camera_model"));
  cmd.add(make_switch('P', "use_pose_prior"));
  cmd.add(make_option('W', sPriorWeights, "prior_weights"));
  cmd.add(make_option('m', i_GPS_XYZ_method, "gps_to_xyz_method"));

  cmd.add(make_option('M', sMethod, "method"));
  cmd.add(make_option('l', min_atom_size, "min_atom_size"));
  cmd.add(make_option('u', max_atom_size, "max_atom_size"));
  cmd.add(make_option('n', knn_k, "knn_k"));
  cmd.add(make_option('N', atom_neighbors, "atom_neighbors"));
  cmd.add(make_option('V', vocab_size, "vocab_size"));
  cmd.add(make_option('O', orb_features, "orb_features"));
  cmd.add(make_option('R', resize_dim, "resize_dim"));
  cmd.add(make_option('H', hybrid_alpha, "hybrid_alpha"));

  cmd.add(make_option('F', sOutputFormat, "output_format"));
  cmd.add(make_option('p', sMapFilename, "map_filename"));
  cmd.add(make_option('b', sBasemap, "basemap"));

  try
  {
    if (argc == 1) throw std::string("Invalid command line parameter.");
    cmd.process(argc, argv);
  }
  catch (const std::string &s)
  {
    OPENMVG_LOG_INFO << "Usage: " << argv[0] << '\n'
      << "[-i|--imageDirectory]\n"
      << "[-d|--sensorWidthDatabase]\n"
      << "[-o|--outputDirectory]\n"
      << "[-f|--focal] (pixels)\n"
      << "[-k|--intrinsics] Kmatrix: \"f;0;ppx;0;f;ppy;0;0;1\"\n"
      << "[-c|--camera_model] Camera model type (default: " << static_cast<int>(PINHOLE_CAMERA_RADIAL3) << ", Pinhole radial 3)\n"
      << "[-g|--group_camera_model] 0/1 (default: 1)\n"
      << "[-P|--use_pose_prior] Use pose prior if GPS EXIF pose is available\n"
      << "[-W|--prior_weights] \"x;y;z\" weights for the prior (default: 1.0;1.0;1.0)\n"
      << "[-m|--gps_to_xyz_method] 0: ECEF (default), 1: UTM\n"
      << "\n"
      << "[-M|--method] gps (default) | visual | hybrid\n"
#if !defined(ATOMIC_CLUSTERING_HAS_VISUAL)
      << "    (this build has no OpenCV support: only 'gps' is usable --\n"
      << "     reconfigure with -DOpenMVG_USE_OPENCV=ON for visual/hybrid)\n"
#endif
      << "[-l|--min_atom_size] target minimum images per atom (default: 15)\n"
      << "[-u|--max_atom_size] max images per atom (default: 50)\n"
      << "[-n|--knn_k] neighbors per node in the path-search graph (default: 6)\n"
      << "[-N|--atom_neighbors] nearest atoms each atom links to for the molecule (default: 2)\n"
      << "[-V|--vocab_size] visual vocabulary size, visual/hybrid only (default: 300)\n"
      << "[-O|--orb_features] max ORB features per image, visual/hybrid only (default: 500)\n"
      << "[-R|--resize_dim] downscale longer side to N px before ORB (default: 640)\n"
      << "[-H|--hybrid_alpha] GPS vs visual weight in [0,1], hybrid only (default: 0.5)\n"
      << "\n"
      << "[-F|--output_format] BIN (default) | JSON, for the per-cluster sfm_data files\n"
      << "[-p|--map_filename] Leaflet map file name (default: clusters_map.html)\n"
      << "[-b|--basemap] osm (default, needs internet to view) | none (offline scatter plot)";

    OPENMVG_LOG_ERROR << s;
    return EXIT_FAILURE;
  }

  const bool b_Use_pose_prior = cmd.used('P');

  if (sMethod != "gps" && sMethod != "visual" && sMethod != "hybrid")
  {
    OPENMVG_LOG_ERROR << "--method must be one of: gps, visual, hybrid";
    return EXIT_FAILURE;
  }
  if (sMethod == "visual" || sMethod == "hybrid")
  {
#if !defined(ATOMIC_CLUSTERING_HAS_VISUAL)
    OPENMVG_LOG_ERROR
      << "This build was compiled without OpenCV support (OpenMVG_USE_OPENCV=OFF or OpenCV not found); "
      << "--method visual/hybrid are unavailable. Only --method gps can be used. "
      << "Reconfigure with -DOpenMVG_USE_OPENCV=ON to enable them.";
    return EXIT_FAILURE;
#endif
  }
  if (min_atom_size > max_atom_size)
  {
    OPENMVG_LOG_ERROR << "--min_atom_size must be <= --max_atom_size";
    return EXIT_FAILURE;
  }
  if (sMethod == "hybrid" && (hybrid_alpha < 0.0 || hybrid_alpha > 1.0))
  {
    OPENMVG_LOG_ERROR << "--hybrid_alpha must be between 0.0 and 1.0";
    return EXIT_FAILURE;
  }
  std::transform(sOutputFormat.begin(), sOutputFormat.end(), sOutputFormat.begin(), ::toupper);
  if (sOutputFormat != "BIN" && sOutputFormat != "JSON")
  {
    OPENMVG_LOG_ERROR << "--output_format must be one of: BIN, JSON";
    return EXIT_FAILURE;
  }
  if (sBasemap != "osm" && sBasemap != "none")
  {
    OPENMVG_LOG_ERROR << "--basemap must be one of: osm, none";
    return EXIT_FAILURE;
  }

  double width = -1, height = -1, focal = -1, ppx = -1, ppy = -1;
  const EINTRINSIC e_User_camera_model = EINTRINSIC(i_User_camera_model);

  if (!stlplus::folder_exists(sImageDir))
  {
    OPENMVG_LOG_ERROR << "The input directory doesn't exist";
    return EXIT_FAILURE;
  }
  if (sOutputDir.empty())
  {
    OPENMVG_LOG_ERROR << "Invalid output directory";
    return EXIT_FAILURE;
  }
  if (!stlplus::folder_exists(sOutputDir) && !stlplus::folder_create(sOutputDir))
  {
    OPENMVG_LOG_ERROR << "Cannot create output directory";
    return EXIT_FAILURE;
  }
  if (sKmatrix.size() > 0 && !checkIntrinsicStringValidity(sKmatrix, focal, ppx, ppy))
  {
    OPENMVG_LOG_ERROR << "Invalid K matrix input";
    return EXIT_FAILURE;
  }
  if (sKmatrix.size() > 0 && focal_pixels != -1.0)
  {
    OPENMVG_LOG_ERROR << "Cannot combine -f and -k options";
    return EXIT_FAILURE;
  }

  std::vector<Datasheet> vec_database;
  if (!sfileDatabase.empty() && !parseDatabase(sfileDatabase, vec_database))
  {
    OPENMVG_LOG_ERROR << "Invalid input database: " << sfileDatabase << ", please specify a valid file.";
    return EXIT_FAILURE;
  }

  if (b_Use_pose_prior)
    prior_w_info = checkPriorWeightsString(sPriorWeights);

  std::vector<std::string> vec_image = stlplus::folder_files(sImageDir);
  std::sort(vec_image.begin(), vec_image.end());

  // Configure an empty scene with Views and their corresponding cameras.
  SfM_Data sfm_data;
  sfm_data.s_root_path = sImageDir;
  Views &views = sfm_data.views;
  Intrinsics &intrinsics = sfm_data.intrinsics;

  std::map<AtomicClustering::ImageId, IndexT> pathToViewId;
  std::map<AtomicClustering::ImageId, AtomicClustering::GpsCoord> gpsPoints;
  std::vector<AtomicClustering::ImageId> allImageIds;

  system::LoggerProgress my_progress_bar(vec_image.size(), "- Listing images -");
  std::ostringstream error_report_stream;
  for (auto iter_image = vec_image.begin(); iter_image != vec_image.end(); ++iter_image, ++my_progress_bar)
  {
    width = height = ppx = ppy = focal = -1.0;

    const std::string sImageFilename = stlplus::create_filespec(sImageDir, *iter_image);
    const std::string sImFilenamePart = stlplus::filename_part(sImageFilename);

    if (openMVG::image::GetFormat(sImageFilename.c_str()) == openMVG::image::Unknown)
    {
      error_report_stream << sImFilenamePart << ": Unknown image file format." << "\n";
      continue;
    }
    if (sImFilenamePart.find("mask.png") != std::string::npos ||
        sImFilenamePart.find("_mask.png") != std::string::npos)
    {
      error_report_stream << sImFilenamePart << " is a mask image" << "\n";
      continue;
    }

    ImageHeader imgHeader;
    if (!openMVG::image::ReadImageHeader(sImageFilename.c_str(), &imgHeader))
      continue;

    width = imgHeader.width;
    height = imgHeader.height;
    ppx = width / 2.0;
    ppy = height / 2.0;

    // Single EXIF parse per image, reused for both focal-length resolution
    // and GPS extraction (rather than opening the file twice).
    Exif_IO_EasyExif exifReader;
    const bool bHaveValidExifMetadata =
      exifReader.open(sImageFilename) && exifReader.doesHaveExifInfo();

    if (sKmatrix.size() > 0)
    {
      if (!checkIntrinsicStringValidity(sKmatrix, focal, ppx, ppy))
        focal = -1.0;
    }
    else if (focal_pixels != -1)
      focal = focal_pixels;

    if (focal == -1)
    {
      if (bHaveValidExifMetadata && !exifReader.getModel().empty() && !exifReader.getBrand().empty())
      {
        if (exifReader.getFocal() == 0.0f)
        {
          error_report_stream << stlplus::basename_part(sImageFilename) << ": Focal length is missing." << "\n";
          focal = -1.0;
        }
        else
        {
          const std::string sCamModel = exifReader.getBrand() + " " + exifReader.getModel();
          Datasheet datasheet;
          if (getInfo(sCamModel, vec_database, datasheet))
          {
            const double ccdw = datasheet.sensorSize_;
            focal = std::max(width, height) * exifReader.getFocal() / ccdw;
          }
          else
          {
            error_report_stream
              << stlplus::basename_part(sImageFilename) << "\" model \"" << sCamModel
              << "\" doesn't exist in the database" << "\n"
              << "Please consider add your camera model and sensor width in the database." << "\n";
          }
        }
      }
    }

    std::shared_ptr<IntrinsicBase> intrinsic;
    if (focal > 0 && ppx > 0 && ppy > 0 && width > 0 && height > 0)
    {
      switch (e_User_camera_model)
      {
        case PINHOLE_CAMERA:
          intrinsic = std::make_shared<Pinhole_Intrinsic>(width, height, focal, ppx, ppy);
          break;
        case PINHOLE_CAMERA_RADIAL1:
          intrinsic = std::make_shared<Pinhole_Intrinsic_Radial_K1>(width, height, focal, ppx, ppy, 0.0);
          break;
        case PINHOLE_CAMERA_RADIAL3:
          intrinsic = std::make_shared<Pinhole_Intrinsic_Radial_K3>(width, height, focal, ppx, ppy, 0.0, 0.0, 0.0);
          break;
        case PINHOLE_CAMERA_BROWN:
          intrinsic = std::make_shared<Pinhole_Intrinsic_Brown_T2>(width, height, focal, ppx, ppy, 0.0, 0.0, 0.0, 0.0, 0.0);
          break;
        case PINHOLE_CAMERA_FISHEYE:
          intrinsic = std::make_shared<Pinhole_Intrinsic_Fisheye>(width, height, focal, ppx, ppy, 0.0, 0.0, 0.0, 0.0);
          break;
        case CAMERA_SPHERICAL:
          intrinsic = std::make_shared<Intrinsic_Spherical>(width, height);
          break;
        default:
          OPENMVG_LOG_ERROR << "Error: unknown camera model: " << static_cast<int>(e_User_camera_model);
          return EXIT_FAILURE;
      }
    }

    // GPS: read regardless of --method, since the Leaflet map is generated
    // from whatever GPS is available no matter which clustering method ran.
    AtomicClustering::GpsCoord gps;
    double lat = 0.0, lon = 0.0, alt = 0.0;
    if (bHaveValidExifMetadata && exifReader.GPSLatitude(&lat) && exifReader.GPSLongitude(&lon))
    {
      exifReader.GPSAltitude(&alt); // best-effort; altitude defaults to 0 if absent
      gps.lat = lat;
      gps.lon = lon;
      gps.valid = true;
    }

    const IndexT view_id = static_cast<IndexT>(views.size());
    if (b_Use_pose_prior && gps.valid)
    {
      ViewPriors v(*iter_image, view_id, view_id, view_id, width, height);
      if (!intrinsic)
        v.id_intrinsic = UndefinedIndexT;
      else
        intrinsics[v.id_intrinsic] = intrinsic;

      const Vec3 pose_center = (i_GPS_XYZ_method == 1)
        ? lla_to_utm(gps.lat, gps.lon, alt)
        : lla_to_ecef(gps.lat, gps.lon, alt);
      v.SetPoseCenterPrior(pose_center, prior_w_info.first ? prior_w_info.second : Vec3::Constant(1.0));

      views[v.id_view] = std::make_shared<ViewPriors>(v);
    }
    else
    {
      View v(*iter_image, view_id, view_id, view_id, width, height);
      if (!intrinsic)
        v.id_intrinsic = UndefinedIndexT;
      else
        intrinsics[v.id_intrinsic] = intrinsic;

      views[v.id_view] = std::make_shared<View>(v);
    }

    pathToViewId[sImageFilename] = view_id;
    allImageIds.push_back(sImageFilename);
    if (gps.valid)
      gpsPoints[sImageFilename] = gps;
  }

  if (!error_report_stream.str().empty())
  {
    OPENMVG_LOG_WARNING << "Warning & Error messages:\n" << error_report_stream.str();
  }

  if (b_Group_camera_model)
    GroupSharedIntrinsics(sfm_data);

  OPENMVG_LOG_INFO
    << "Listed #File(s): " << vec_image.size() << "\n"
    << "Usable #View(s): " << sfm_data.GetViews().size() << "\n"
    << "Usable #Intrinsic(s): " << sfm_data.GetIntrinsics().size() << "\n"
    << "Geotagged #View(s): " << gpsPoints.size();

  // --- Clustering ---
  std::vector<AtomicClustering::ImageId> geotaggedIds;
  for (const auto &kv : gpsPoints)
    geotaggedIds.push_back(kv.first);
  std::sort(geotaggedIds.begin(), geotaggedIds.end());

  AtomicClustering::DistanceFn dist;
  std::vector<AtomicClustering::ImageId> clusterIds;

  try
  {
    if (sMethod == "gps")
    {
      if (geotaggedIds.size() < 2)
      {
        OPENMVG_LOG_ERROR << "Only " << geotaggedIds.size() << " geotagged image(s) found; need at least 2 for --method gps";
        return EXIT_FAILURE;
      }
      clusterIds = geotaggedIds;
      dist = AtomicClustering::build_distance_matrix(clusterIds, AtomicClustering::make_gps_distance_fn(gpsPoints));
    }
#ifdef ATOMIC_CLUSTERING_HAS_VISUAL
    else if (sMethod == "visual")
    {
      if (allImageIds.size() < 2)
      {
        OPENMVG_LOG_ERROR << "Only " << allImageIds.size() << " usable image(s) found; need at least 2 for --method visual";
        return EXIT_FAILURE;
      }
      clusterIds = allImageIds;
      const auto vectors = AtomicClustering::extract_visual_vectors(clusterIds, resize_dim, orb_features, vocab_size);
      dist = AtomicClustering::build_distance_matrix(clusterIds, AtomicClustering::make_visual_distance_fn(vectors));
    }
    else // hybrid
    {
      if (geotaggedIds.size() < 2)
      {
        OPENMVG_LOG_ERROR << "Only " << geotaggedIds.size() << " geotagged image(s) found; need at least 2 for --method hybrid";
        return EXIT_FAILURE;
      }
      clusterIds = geotaggedIds;
      const auto vectors = AtomicClustering::extract_visual_vectors(clusterIds, resize_dim, orb_features, vocab_size);
      dist = AtomicClustering::build_distance_matrix(clusterIds, AtomicClustering::make_hybrid_distance_fn(gpsPoints, vectors, hybrid_alpha));
    }
#endif
  }
  catch (const std::exception &e)
  {
    OPENMVG_LOG_ERROR << "Error while building the distance function: " << e.what();
    return EXIT_FAILURE;
  }

  const std::vector<AtomicClustering::Cluster> clusters = AtomicClustering::run_atomic_clustering(
    clusterIds, dist, static_cast<int>(min_atom_size), static_cast<int>(max_atom_size),
    static_cast<int>(knn_k), static_cast<int>(atom_neighbors));

  OPENMVG_LOG_INFO << "Number of clusters = " << clusters.size();

  const std::string sExt = (sOutputFormat == "JSON") ? "json" : "bin";
  for (size_t i = 0; i < clusters.size(); ++i)
  {
    std::ostringstream filename;
    filename << "sfm_data_" << std::setw(4) << std::setfill('0') << i << "." << sExt;
    const std::string outPath = stlplus::create_filespec(sOutputDir, filename.str());

    OPENMVG_LOG_INFO << "Writing cluster " << i << " (" << clusters[i].members.size() << " image(s)) to " << outPath;
    if (!ExportClusterSfMData(sfm_data, outPath, clusters[i].members, pathToViewId))
    {
      OPENMVG_LOG_ERROR << "Could not write cluster: " << outPath;
    }
  }

  const std::string sMapPath = stlplus::create_filespec(sOutputDir, sMapFilename);
  WriteLeafletMap(sMapPath, clusters, gpsPoints, sBasemap);

  return EXIT_SUCCESS;
}
