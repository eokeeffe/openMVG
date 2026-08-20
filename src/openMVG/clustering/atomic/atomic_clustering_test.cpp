// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "atomic_clustering.hpp"
#include "distance.hpp"
#include "types.hpp"
#ifdef ATOMIC_CLUSTERING_HAS_VISUAL
#include "visual_features.hpp"
#endif

#include "CppUnitLite/TestHarness.h"
#include "testing/testing.h"

#include <set>
#include <sstream>

using namespace AtomicClustering;

namespace {

// Every id shows up in at least one cluster -- atoms are a hard partition of
// the whole input, so coverage should always be total.
std::set<ImageId> CoveredIds(const std::vector<Cluster> &clusters)
{
  std::set<ImageId> covered;
  for (const auto &c : clusters)
    for (const auto &m : c.members)
      covered.insert(m);
  return covered;
}

} // namespace

// -- GPS-only: no images, no OpenCV -- always built and run. --

TEST(AtomicClustering, GpsOnly_TwoWellSeparatedGroupsSplitIntoAtoms)
{
  // Two tight groups of points, ~157 km apart (1 degree of lat/lon at the
  // equator) -- far larger than the jitter within each group, so the
  // atoms should cleanly separate them.
  std::map<ImageId, GpsCoord> points;
  for (int i = 0; i < 10; ++i)
  {
    std::ostringstream name;
    name << "groupA_" << i;
    points[name.str()] = GpsCoord{0.0 + i * 0.0001, 0.0 + i * 0.0001, true};
  }
  for (int i = 0; i < 10; ++i)
  {
    std::ostringstream name;
    name << "groupB_" << i;
    points[name.str()] = GpsCoord{1.0 + i * 0.0001, 1.0 + i * 0.0001, true};
  }

  std::vector<ImageId> ids;
  for (const auto &kv : points)
    ids.push_back(kv.first);

  const DistanceFn dist = build_distance_matrix(ids, make_gps_distance_fn(points));
  const std::vector<Cluster> clusters =
    run_atomic_clustering(ids, dist, /*min_atom_size=*/3, /*max_atom_size=*/8, /*knn_k=*/4, /*atom_neighbors=*/2);

  CHECK(!clusters.empty());
  EXPECT_EQ(ids.size(), CoveredIds(clusters).size());

  // groupA and groupB members should not be mixed into the same atom
  // (the molecule, if any, is expected to bridge them -- that's fine).
  bool found_pure_group = false;
  for (const auto &c : clusters)
  {
    bool has_a = false, has_b = false;
    for (const auto &m : c.members)
    {
      if (m.rfind("groupA", 0) == 0) has_a = true;
      if (m.rfind("groupB", 0) == 0) has_b = true;
    }
    if (has_a && !has_b) found_pure_group = true;
  }
  CHECK(found_pure_group);
}

TEST(AtomicClustering, GpsOnly_TooFewImagesReturnsSingleCluster)
{
  std::map<ImageId, GpsCoord> points = {
    {"a", GpsCoord{0.0, 0.0, true}},
    {"b", GpsCoord{0.001, 0.001, true}},
  };
  std::vector<ImageId> ids = {"a", "b"};

  const DistanceFn dist = build_distance_matrix(ids, make_gps_distance_fn(points));
  const std::vector<Cluster> clusters = run_atomic_clustering(ids, dist, 15, 50, 6, 2);

  CHECK_EQUAL(1, static_cast<int>(clusters.size()));
  EXPECT_EQ(size_t(2), clusters[0].members.size());
}

// -- Visual-only / hybrid: need real images + OpenCV, only built when the
// library itself was compiled with ATOMIC_CLUSTERING_HAS_VISUAL. --

#ifdef ATOMIC_CLUSTERING_HAS_VISUAL

namespace {
std::string ExampleImage(const char *filename)
{
  return std::string(THIS_SOURCE_DIR) + "/example_images/" + filename;
}
} // namespace

TEST(AtomicClustering, VisualOnly_RealImagesSplitIntoAtomsAndMolecule)
{
  const std::vector<ImageId> ids = {ExampleImage("fraumunster00.png"), ExampleImage("fraumunster01.png")};

  const auto vectors = extract_visual_vectors(ids, /*resize_max_dim=*/640, /*n_features=*/500, /*vocab_size=*/50);
  const DistanceFn dist = build_distance_matrix(ids, make_visual_distance_fn(vectors));
  // max_atom_size=1 forces a split (each image its own atom) so the
  // knn-graph + molecule-bridging code path actually runs, not just the
  // single-cluster short-circuit.
  const std::vector<Cluster> clusters =
    run_atomic_clustering(ids, dist, /*min_atom_size=*/1, /*max_atom_size=*/1, /*knn_k=*/1, /*atom_neighbors=*/1);

  CHECK(!clusters.empty());
  EXPECT_EQ(ids.size(), CoveredIds(clusters).size());
}

TEST(AtomicClustering, Hybrid_RealImagesWithSyntheticGps)
{
  const std::vector<ImageId> ids = {ExampleImage("fraumunster00.png"), ExampleImage("fraumunster01.png")};
  const std::map<ImageId, GpsCoord> points = {
    {ids[0], GpsCoord{47.3700, 8.5400, true}},
    {ids[1], GpsCoord{47.3701, 8.5401, true}},
  };

  const auto vectors = extract_visual_vectors(ids, 640, 500, 50);
  const DistanceFn dist = build_distance_matrix(ids, make_hybrid_distance_fn(points, vectors, /*alpha=*/0.5));
  const std::vector<Cluster> clusters = run_atomic_clustering(ids, dist, 1, 1, 1, 1);

  CHECK(!clusters.empty());
  EXPECT_EQ(ids.size(), CoveredIds(clusters).size());
}

#endif // ATOMIC_CLUSTERING_HAS_VISUAL

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr); }
/* ************************************************************************* */
