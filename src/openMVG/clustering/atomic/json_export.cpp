#include "json_export.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>

namespace AtomicClustering {

namespace {

std::string json_escape(const std::string &s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        switch (c) {
            case '"':
                out += "\\\"";
                break;
            case '\\':
                out += "\\\\";
                break;
            case '\n':
                out += "\\n";
                break;
            default:
                out += c;
        }
    }
    return out;
}

} // namespace

void export_json(const std::vector<Cluster> &clusters, const std::map<ImageId, GpsCoord> &points,
                  const std::string &output_path) {
    std::ofstream out(output_path);
    out << std::fixed << std::setprecision(7);
    out << "[\n";

    for (size_t i = 0; i < clusters.size(); ++i) {
        const auto &members = clusters[i].members;

        out << "  {\n";
        out << "    \"cluster_id\": " << (i + 1) << ",\n";
        out << "    \"images\": [\n";
        for (size_t j = 0; j < members.size(); ++j) {
            const GpsCoord &p = points.at(members[j]);
            out << "      {\"path\": \"" << json_escape(members[j]) << "\", "
                << "\"lat\": " << p.lat << ", \"lon\": " << p.lon << "}";
            out << (j + 1 < members.size() ? ",\n" : "\n");
        }
        out << "    ]\n";
        out << "  }" << (i + 1 < clusters.size() ? ",\n" : "\n");
    }

    out << "]\n";
    std::cout << "\nWrote " << output_path << "\n";
}

void summarize(const std::vector<Cluster> &clusters, const std::map<ImageId, GpsCoord> &points) {
    std::map<ImageId, int> membership_count;
    std::set<ImageId> all_clustered;

    for (size_t i = 0; i < clusters.size(); ++i) {
        const auto &members = clusters[i].members;

        double clat = 0.0, clon = 0.0;
        for (const auto &p : members) {
            const GpsCoord &c = points.at(p);
            clat += c.lat;
            clon += c.lon;
            all_clustered.insert(p);
            membership_count[p]++;
        }
        clat /= static_cast<double>(members.size());
        clon /= static_cast<double>(members.size());

        std::cout << "\nCluster " << (i + 1) << ": " << members.size() << " image(s), "
                  << "centroid=(" << clat << ", " << clon << ")\n";
        for (const auto &p : members) {
            std::cout << "    " << p << "\n";
        }
    }

    int overlapping = 0;
    for (const auto &kv : membership_count) {
        if (kv.second > 1) {
            ++overlapping;
        }
    }

    int unclustered = static_cast<int>(points.size()) - static_cast<int>(all_clustered.size());

    std::cout << "\n" << clusters.size() << " cluster(s) total.\n";
    std::cout << overlapping << " image(s) appear in more than one cluster.\n";
    std::cout << unclustered << " image(s) did not join any cluster.\n";
}

} // namespace gpsatomic
