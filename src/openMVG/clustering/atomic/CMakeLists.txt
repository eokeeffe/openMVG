# Atomic Clustering: GPS / visual / hybrid clustering with overlap, see
# README.md. The core (GPS-only) pipeline only needs LEMON (for the
# molecule's shortest-path search) and easyexif (for GPS EXIF) -- both are
# already vendored and always built by OpenMVG, so it never needs OpenCV.
#
# Visual and hybrid clustering additionally need OpenCV (for ORB features +
# k-means), so visual_features.cpp is only compiled in -- and the
# ATOMIC_CLUSTERING_HAS_VISUAL define only set -- when OpenMVG_USE_OPENCV is
# enabled and OpenCV is actually found.

if (OpenMVG_USE_OPENCV)
  # Only request the modules actually used (ORB + imread/resize): an
  # unscoped find_package(OpenCV) pulls in OpenCV_LIBS for every module in
  # the install (opencv_sfm, opencv_hdf, opencv_videoio, ...), several of
  # which commonly drag in unrelated, easily version-mismatched system
  # deps (ffmpeg, HDF5, GDAL, glog) that this library never needs to link.
  find_package(OpenCV QUIET COMPONENTS core imgproc features2d imgcodecs)
  if (NOT OpenCV_FOUND)
    # OpenCV 5.x renamed the "features2d" module to "features" (the
    # opencv2/features2d.hpp header path is unchanged).
    find_package(OpenCV QUIET COMPONENTS core imgproc features imgcodecs)
  endif()
endif()

set(ATOMIC_CLUSTERING_HAS_VISUAL OFF)
if (OpenMVG_USE_OPENCV AND OpenCV_FOUND)
  set(ATOMIC_CLUSTERING_HAS_VISUAL ON)
endif()
# Let callers of add_subdirectory(atomic) (openMVG/clustering/CMakeLists.txt,
# software/SfM/clustering/CMakeLists.txt, and its own test target below) know
# whether --method visual/hybrid will be usable.
set(ATOMIC_CLUSTERING_HAS_VISUAL ${ATOMIC_CLUSTERING_HAS_VISUAL} PARENT_SCOPE)

add_library(openMVG_clustering STATIC
    distance.cpp
    atomic_clustering.cpp
    exif_reader.cpp
    json_export.cpp

    distance.hpp
    types.hpp
    atomic_clustering.hpp
    exif_reader.hpp
    json_export.hpp
)

target_include_directories(openMVG_clustering
    PRIVATE
        # PRIVATE: only exif_reader.cpp's "exif.h" needs this, no public
        # header of this library exposes easyexif types. LEMON's own
        # include dirs come transitively via linking ${LEMON_LIBRARY}
        # below (openMVG_lemon already exports them, generator-expression
        # wrapped, same as openMVG_graph does) -- adding them again here
        # as raw absolute paths would break `install(EXPORT ...)` below.
        "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/third_party/easyexif>"
)
target_link_libraries(openMVG_clustering
    PRIVATE
        ${LEMON_LIBRARY}
        openMVG_easyexif
    PUBLIC
        # PUBLIC: pulls in ${CMAKE_SOURCE_DIR} as a (BUILD_INTERFACE) include
        # dir, so exif_reader.cpp's quoted "third_party/..." includes
        # resolve, and so it's still there for anything that links
        # openMVG_clustering (main_ComputeAtomicClustering, and
        # atomic_clustering_test.cpp's "testing/testing.h", which itself
        # needs Eigen via openMVG_numeric). Both already wrap their own
        # include dirs in BUILD_INTERFACE/INSTALL_INTERFACE themselves.
        openMVG_system
        openMVG_numeric
)
target_compile_options(openMVG_clustering PRIVATE -Wall -Wextra)
set_target_properties(openMVG_clustering PROPERTIES SOVERSION ${OPENMVG_VERSION_MAJOR} VERSION "${OPENMVG_VERSION_MAJOR}.${OPENMVG_VERSION_MINOR}")
set_property(TARGET openMVG_clustering PROPERTY FOLDER OpenMVG/OpenMVG)
install(TARGETS openMVG_clustering DESTINATION ${CMAKE_INSTALL_LIBDIR} EXPORT openMVG-targets)

if (ATOMIC_CLUSTERING_HAS_VISUAL)
  target_sources(openMVG_clustering PRIVATE visual_features.cpp visual_features.hpp)
  # PUBLIC + BUILD_INTERFACE-wrapped: visual_features.hpp itself includes
  # <opencv2/core.hpp>, so consumers (main_ComputeAtomicClustering, the
  # unit test) need this on their include path too when it's enabled.
  # NOTE: the generator expression MUST be quoted -- OpenCV_INCLUDE_DIRS is
  # a multi-value (semicolon-separated) list, and an unquoted
  # $<BUILD_INTERFACE:${OpenCV_INCLUDE_DIRS}> gets split into separate
  # command arguments on those semicolons before CMake ever sees it as one
  # generator expression, corrupting the first/last entries.
  target_include_directories(openMVG_clustering PUBLIC "$<BUILD_INTERFACE:${OpenCV_INCLUDE_DIRS}>")
  target_link_libraries(openMVG_clustering PRIVATE ${OpenCV_LIBS})
  # PUBLIC: anything linking openMVG_clustering (the main tool, the unit
  # test) needs to see this too, since it gates whether they can even
  # reference visual_features.hpp / --method visual|hybrid.
  target_compile_definitions(openMVG_clustering PUBLIC ATOMIC_CLUSTERING_HAS_VISUAL=1)
else()
  message(STATUS "openMVG_clustering: OpenCV not enabled/found -- building GPS-only "
                  "atomic clustering (visual & hybrid methods will be unavailable)")
endif()

UNIT_TEST(openMVG atomic_clustering "openMVG_clustering")
