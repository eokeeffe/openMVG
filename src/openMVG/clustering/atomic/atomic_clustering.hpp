#ifndef __ATOMIC_CLUSTERING_HPP__
#define __ATOMIC_CLUSTERING_HPP__

#include "distance.hpp"
#include "types.hpp"

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace AtomicClustering {

    // Adapted from Xie et al., "Hierarchical Clustering-Aligning Framework
    // Based Fast Large-Scale 3D Reconstruction Using Aerial Imagery" (Remote
    // Sensing, 2019). Their original method clusters by image-feature
    // similarity (a precomputed vocabulary tree + ORB) for SfM reconstruction.
    // Here the same two-stage algorithmic idea is generalized to run over any
    // `DistanceFn`, so it can be driven by GPS distance, visual dissimilarity,
    // or a blend of both:
    //
    //   1. "Atoms": a hard partition built by growing regions outward from
    //      well-spaced seed points (farthest-point sampling), each capped to
    //      a size range.
    //   2. "Molecule": adjacent atoms are connected via shortest paths in a
    //      k-nearest-neighbor graph; every image on one of those paths is
    //      added to a shared "molecule" cluster. Because the molecule's
    //      images are also still members of their original atom, this is
    //      where the overlap comes from.

    // Greedily picks n_sources ids that are maximally spread out (by `dist`).
    // These become the seeds atoms grow outward from.
    std::vector<ImageId> farthest_point_sample(const std::vector<ImageId> &ids, const DistanceFn &dist,
                                                int n_sources);

    // Grows one atom per source by repeatedly assigning the nearest unassigned
    // id to whichever atom currently has a member closest to it, capped at
    // max_atom_size. Returns {source -> members}; a source's own membership
    // always includes itself.
    std::map<ImageId, std::vector<ImageId>> generate_atoms(const std::vector<ImageId> &ids,
                                                            const std::vector<ImageId> &sources,
                                                            const DistanceFn &dist, int max_atom_size);

    // Opaque handle around a LEMON ListGraph + node/edge maps (defined in the
    // .cpp so LEMON headers don't leak into every translation unit that
    // includes this header).
    struct KnnGraph;

    // Builds an unweighted (weight=1 per edge) k-nearest-neighbor graph, used
    // for shortest-path search when generating the molecule. Edge weight is
    // fixed at 1 so Dijkstra finds the path with the fewest hops (fewest
    // images), not the geometrically/visually shortest one -- matching the
    // paper's own choice for this step.
    std::shared_ptr<KnnGraph> build_knn_graph(const std::vector<ImageId> &ids, const DistanceFn &dist, int k);

    // Connects each atom to its atom_neighbor_count nearest atoms (by
    // source-to-source `dist`), then finds the shortest path between each
    // pair of atom sources in the knn graph (via LEMON's Dijkstra). Every id
    // on any of those paths joins the returned molecule.
    std::vector<ImageId> generate_molecule(const std::map<ImageId, std::vector<ImageId>> &atom_members,
                                            const DistanceFn &dist, const std::shared_ptr<KnnGraph> &graph,
                                            int atom_neighbor_count);

    // Full pipeline: generate atoms, then a molecule, driven entirely by
    // `dist`. Returns atoms + the molecule (if any) as a flat list of
    // clusters.
    std::vector<Cluster> run_atomic_clustering(const std::vector<ImageId> &ids, const DistanceFn &dist,
                                                int min_atom_size, int max_atom_size, int knn_k,
                                                int atom_neighbor_count);

} // namespace AtomicClustering

#endif /* __ATOMIC_CLUSTERING_HPP__ */
