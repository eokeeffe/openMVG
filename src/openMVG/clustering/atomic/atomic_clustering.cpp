#include "atomic_clustering.hpp"

#include <lemon/core.h>
#include <lemon/dijkstra.h>
#include <lemon/list_graph.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <set>
#include <unordered_map>
#include <unordered_set>

namespace AtomicClustering {

    struct KnnGraph {
        lemon::ListGraph g;
        lemon::ListGraph::EdgeMap<int> weight{g};
        std::unordered_map<ImageId, lemon::ListGraph::Node> node_of;
        std::unordered_map<int, ImageId> id_of_node; // keyed by g.id(node)
    };

    std::vector<ImageId> farthest_point_sample(const std::vector<ImageId> &ids, const DistanceFn &dist,
                                                int n_sources) {
        if (n_sources >= static_cast<int>(ids.size())) {
            return ids;
        }

        std::mt19937 rng(42);
        std::uniform_int_distribution<size_t> pick(0, ids.size() - 1);

        std::vector<ImageId> selected = {ids[pick(rng)]};
        std::unordered_set<ImageId> selected_set(selected.begin(), selected.end());

        while (static_cast<int>(selected.size()) < n_sources) {
            ImageId best_key;
            double best_dist = -1.0;

            for (const auto &k : ids) {
                if (selected_set.count(k)) {
                    continue;
                }
                double d = std::numeric_limits<double>::max();
                for (const auto &s : selected) {
                    d = std::min(d, dist(k, s));
                }
                if (d > best_dist) {
                    best_dist = d;
                    best_key = k;
                }
            }

            selected.push_back(best_key);
            selected_set.insert(best_key);
        }

        return selected;
    }

    std::map<ImageId, std::vector<ImageId>> generate_atoms(const std::vector<ImageId> &ids,
                                                            const std::vector<ImageId> &sources,
                                                            const DistanceFn &dist, int max_atom_size) {
        std::map<ImageId, std::vector<ImageId>> members;
        for (const auto &s : sources) {
            members[s] = {s};
        }

        std::unordered_set<ImageId> unassigned(ids.begin(), ids.end());
        for (const auto &s : sources) {
            unassigned.erase(s);
        }

        while (!unassigned.empty()) {
            bool progressed = false;

            for (const auto &s : sources) {
                if (unassigned.empty()) {
                    break;
                }
                if (static_cast<int>(members[s].size()) >= max_atom_size) {
                    continue;
                }

                ImageId nearest;
                double nearest_d = std::numeric_limits<double>::max();
                for (const auto &u : unassigned) {
                    double d = std::numeric_limits<double>::max();
                    for (const auto &m : members[s]) {
                        d = std::min(d, dist(u, m));
                    }
                    if (d < nearest_d) {
                        nearest_d = d;
                        nearest = u;
                    }
                }

                members[s].push_back(nearest);
                unassigned.erase(nearest);
                progressed = true;
            }

            if (!progressed) {
                break; // every atom is at max_atom_size but ids remain
            }
        }

        // Leftover ids (all atoms capped): assign to nearest atom regardless of cap.
        for (auto it = unassigned.begin(); it != unassigned.end();) {
            const ImageId &u = *it;

            ImageId best_s;
            double best_d = std::numeric_limits<double>::max();
            for (const auto &s : sources) {
                double d = std::numeric_limits<double>::max();
                for (const auto &m : members[s]) {
                    d = std::min(d, dist(u, m));
                }
                if (d < best_d) {
                    best_d = d;
                    best_s = s;
                }
            }

            members[best_s].push_back(u);
            it = unassigned.erase(it);
        }

        return members;
    }

    std::shared_ptr<KnnGraph> build_knn_graph(const std::vector<ImageId> &ids, const DistanceFn &dist, int k) {
        auto kg = std::make_shared<KnnGraph>();

        for (const auto &id : ids) {
            auto node = kg->g.addNode();
            kg->node_of[id] = node;
            kg->id_of_node[kg->g.id(node)] = id;
        }

        for (const auto &a : ids) {
            std::vector<ImageId> others;
            others.reserve(ids.size());
            for (const auto &b : ids) {
                if (b != a) {
                    others.push_back(b);
                }
            }

            int limit = std::min<int>(k, static_cast<int>(others.size()));
            std::partial_sort(others.begin(), others.begin() + limit, others.end(),
                            [&](const ImageId &x, const ImageId &y) { return dist(a, x) < dist(a, y); });

            auto na = kg->node_of[a];
            for (int i = 0; i < limit; ++i) {
                auto nb = kg->node_of[others[static_cast<size_t>(i)]];
                if (lemon::findEdge(kg->g, na, nb) == lemon::INVALID) {
                    auto e = kg->g.addEdge(na, nb);
                    kg->weight[e] = 1;
                }
            }
        }

        return kg;
    }

    std::vector<ImageId> generate_molecule(const std::map<ImageId, std::vector<ImageId>> &atom_members,
                                            const DistanceFn &dist, const std::shared_ptr<KnnGraph> &graph,
                                            int atom_neighbor_count) {
        std::vector<ImageId> sources;
        sources.reserve(atom_members.size());
        for (const auto &kv : atom_members) {
            sources.push_back(kv.first);
        }
        if (sources.size() < 2) {
            return {};
        }

        std::set<std::pair<ImageId, ImageId>> pairs;
        for (const auto &s : sources) {
            std::vector<ImageId> others;
            others.reserve(sources.size());
            for (const auto &o : sources) {
                if (o != s) {
                    others.push_back(o);
                }
            }

            int limit = std::min<int>(atom_neighbor_count, static_cast<int>(others.size()));
            std::partial_sort(others.begin(), others.begin() + limit, others.end(),
                            [&](const ImageId &x, const ImageId &y) { return dist(s, x) < dist(s, y); });

            for (int i = 0; i < limit; ++i) {
                ImageId a = s, b = others[static_cast<size_t>(i)];
                if (b < a) {
                    std::swap(a, b);
                }
                pairs.insert({a, b});
            }
        }

        std::set<ImageId> molecule;
        lemon::Dijkstra<lemon::ListGraph, lemon::ListGraph::EdgeMap<int>> dijkstra(graph->g, graph->weight);

        for (const auto &pr : pairs) {
            auto na = graph->node_of.at(pr.first);
            auto nb = graph->node_of.at(pr.second);

            dijkstra.run(na, nb);
            if (!dijkstra.reached(nb)) {
                continue; // no path between these two atom sources in the knn graph
            }

            for (auto v = nb; v != na; v = dijkstra.predNode(v)) {
                molecule.insert(graph->id_of_node.at(graph->g.id(v)));
            }
            molecule.insert(graph->id_of_node.at(graph->g.id(na)));
        }

        return std::vector<ImageId>(molecule.begin(), molecule.end());
    }

    std::vector<Cluster> run_atomic_clustering(const std::vector<ImageId> &ids, const DistanceFn &dist,
                                                int min_atom_size, int max_atom_size, int knn_k,
                                                int atom_neighbor_count) {
        if (static_cast<int>(ids.size()) <= max_atom_size) {
            std::cerr << "[atomic] Only " << ids.size() << " image(s), which is <= --max-atom-size ("
                    << max_atom_size
                    << "), so no splitting is needed -- returning everything as a single cluster.\n"
                    << "[atomic] If you wanted multiple clusters, lower --max-atom-size "
                    << "(and --min-atom-size below it) to force a split.\n";

            std::vector<ImageId> sorted_ids = ids;
            std::sort(sorted_ids.begin(), sorted_ids.end());
            return {Cluster{sorted_ids}};
        }

        double target_size = (min_atom_size + max_atom_size) / 2.0;
        int n_sources = std::max(2, static_cast<int>(std::lround(ids.size() / target_size)));
        auto sources = farthest_point_sample(ids, dist, n_sources);

        auto atom_members = generate_atoms(ids, sources, dist, max_atom_size);
        auto knn_graph = build_knn_graph(ids, dist, knn_k);
        auto molecule = generate_molecule(atom_members, dist, knn_graph, atom_neighbor_count);

        std::cout << "Generated " << atom_members.size() << " atom(s) (target size ~"
                << static_cast<int>(target_size) << ", range [" << min_atom_size << ", " << max_atom_size
                << "])\n";
        if (!molecule.empty()) {
            std::cout << "Generated 1 molecule bridging adjacent atoms: " << molecule.size() << " image(s)\n";
        }

        std::vector<Cluster> clusters;
        clusters.reserve(atom_members.size() + 1);
        for (const auto &kv : atom_members) {
            std::vector<ImageId> sorted_mem = kv.second;
            std::sort(sorted_mem.begin(), sorted_mem.end());
            clusters.push_back(Cluster{sorted_mem});
        }
        if (!molecule.empty()) {
            clusters.push_back(Cluster{molecule});
        }

        return clusters;
}

} // namespace AtomicClustering
