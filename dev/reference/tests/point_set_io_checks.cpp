#include "point_set.h"
#include "grid_io.h"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>

namespace fs = std::filesystem;
using namespace without_sgems;

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}
template<class F> void rejects(F operation, const std::string& message) {
    try { operation(); } catch (const std::exception&) { return; }
    throw std::runtime_error(message);
}
void put(const fs::path& path, const std::string& value) { std::ofstream out(path); out << value; }

int main() {
    try {
        const fs::path dir = fs::current_path() / "test_output" / "point_set_io";
        fs::create_directories(dir);
        const auto file = dir / "fractional.dat";
        put(file, "point set\n5\nx\ny\nz\nphi\nother\n-100.125 0.25 12.75 1 100\n-100.125 0.25 12.75 2 200\n0.75 0.25 12.75 -999 300\n1.25 0.25 12.75 nan 400\n");
        auto data = load_point_set(file, "PHI", 3, -999);
        require(data.nodes.size() == 4 && data.nodes[0].x == -100.125 && data.nodes[0].z == 12.75, "Physical coordinates changed");
        require(data.nodes[0].informed && data.nodes[1].informed && !data.nodes[2].informed && !data.nodes[3].informed, "Missing-value selection failed");
        require(data.nodes[0].index != data.nodes[1].index, "Coincident observations were collapsed");
        const auto other = load_point_set(file, "other", 0);
        require(other.dimensions == 3 && other.nodes[3].value == 400, "Named property selection failed");
        rejects([&] { load_point_set(file, "absent"); }, "Missing property accepted");
        rejects([&] { load_hosc_grid(file, "phi", 3); }, "Fractional coordinates silently rounded by regular loader");
        PointSpatialIndex index(data);
        require(index.query_radius({-100.125,0.25,12.75},0) == std::vector<size_t>({0,1}), "Duplicate coordinates or missing flags mishandled");
        write_point_set(dir / "roundtrip.dat", "roundtrip", "phi", data.nodes);
        const auto again = load_point_set(dir / "roundtrip.dat", "phi", 3);
        require(again.nodes[0].x == data.nodes[0].x && !again.nodes[2].informed, "Point output did not preserve coordinates/missing values");
        put(dir / "named.csv", "x,y,phi,other\n-0.2,2.3,4,8\n0.1,2.5,6,9\n");
        const auto csv = load_point_set(dir / "named.csv", "phi");
        require(csv.dimensions == 2 && csv.nodes[0].x == -0.2 && csv.nodes[1].value == 6, "Named-column CSV inference failed");
        put(dir / "integer.csv", "x,y,phi\n1,2,3\n4,5,6\n");
        const auto integer_csv = load_point_set(dir / "integer.csv", "phi");
        require(integer_csv.nodes.size() == 2 && integer_csv.nodes[1].value == 6, "Integer CSV misclassified as legacy");
        put(dir / "integer1d.csv", "x,phi\n1,2\n4,5\n");
        require(load_point_set(dir / "integer1d.csv", "phi").dimensions == 1, "One-dimensional CSV misclassified");
        rejects([&] { load_point_set(file, "phi", 2); }, "Nonzero Z coordinate silently flattened");
        put(dir / "flat.dat", "flat\n4\nx\ny\nz\nphi\n1.5 2.5 0 7\n");
        require(load_point_set(dir / "flat.dat", "phi", 2).nodes[0].value == 7, "Explicit flat 2-D input rejected");
        put(dir / "legacy.dat", "legacy\n2 4\n# columns\n# data\n7 -0.25 2.5 0 9\n11 0.75 2.5 0 12\n");
        const auto legacy = load_point_set(dir / "legacy.dat", "phi", 2);
        require(legacy.nodes[0].index == 7 && legacy.nodes[1].x == .75, "Legacy point IDs or fractional coordinates lost");
        put(dir / "invalid.dat", "x y value\nnan 0 4\n");
        rejects([&] { load_point_set(dir / "invalid.dat", "value", 2); }, "Nonfinite coordinate accepted");
        put(dir / "invalid.dat", "x y value\n1oops 0 4\n");
        rejects([&] { load_point_set(dir / "invalid.dat", "value", 2); }, "Malformed numeric coordinate accepted");
        put(dir / "type.ini", "[general]\nsource_type=point_set\nreg_irreg=0\n");
        require(uses_point_set(IniConfig::load(dir / "type.ini")), "Point-set selection failed");
        put(dir / "type.ini", "[general]\nsource_type=cartesian\nreg_irreg=0\n");
        rejects([&] { uses_point_set(IniConfig::load(dir / "type.ini")); }, "Conflicting input types accepted");

        // Independent brute-force oracle checks the index on a non-lattice 3-D cloud.
        PointSetData cloud;
        std::mt19937 random(9271);
        std::uniform_real_distribution<double> coordinate(-100,100);
        for (int i = 0; i < 1000; ++i) cloud.nodes.push_back({i,coordinate(random),coordinate(random),coordinate(random),static_cast<double>(i),i % 7 != 0});
        PointSpatialIndex tree(cloud);
        for (int trial = 0; trial < 100; ++trial) {
            std::array<double,3> center = {coordinate(random),coordinate(random),coordinate(random)};
            const double radius = 1.0 + trial;
            std::vector<size_t> expected;
            for (size_t i = 0; i < cloud.nodes.size(); ++i) {
                const auto& p = cloud.nodes[i];
                if (p.informed && std::hypot(p.x-center[0],p.y-center[1],p.z-center[2]) <= radius) expected.push_back(i);
            }
            require(tree.query_radius(center,radius) == expected, "KD radius query disagrees with brute force");
        }
        rejects([&] { tree.query_radius({0,0,0},-1); }, "Negative radius accepted");
        std::cout << "Point-set input, coordinate preservation, validation, and spatial-index checks passed.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
