#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#endif

namespace fs = std::filesystem;
namespace {
using Row = std::array<double, 4>;
const double nan_value = std::numeric_limits<double>::quiet_NaN();

void require(bool condition, const std::string& description) {
    if (!condition) throw std::runtime_error(description);
}

void close(double actual, double expected, const std::string& description) {
    require(std::isfinite(actual) && std::abs(actual - expected) <
        1e-10 * std::max(1.0, std::abs(expected)), description + ": " + std::to_string(actual));
}

void write_points(const fs::path& file, const std::vector<Row>& rows, int dimensions = 3) {
    std::ofstream out(file);
    out << "point set fixture\n" << (dimensions + 1) << "\nx\ny\n";
    if (dimensions == 3) out << "z\n";
    out << "phi\n" << std::setprecision(17);
    for (const auto& row : rows) {
        out << row[0] << ' ' << row[1] << ' ';
        if (dimensions == 3) out << row[2] << ' ';
        if (std::isnan(row[3])) out << "nan\n";
        else out << row[3] << '\n';
    }
}

std::vector<Row> read_points(const fs::path& file) {
    std::ifstream in(file);
    require(static_cast<bool>(in), "Missing point-set result: " + file.string());
    std::string text;
    for (int i = 0; i < 6; ++i) require(static_cast<bool>(std::getline(in, text)), "Truncated GSLIB output");
    std::vector<Row> rows;
    while (std::getline(in, text)) {
        if (text.empty()) continue;
        std::istringstream stream(text);
        Row row;
        std::string field;
        for (double& value : row) {
            require(static_cast<bool>(stream >> field), "Incomplete GSLIB result row");
            value = std::stod(field);
        }
        rows.push_back(row);
    }
    return rows;
}

std::string contents(const fs::path& file) {
    std::ifstream in(file);
    return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

int execute(const fs::path& executable, const fs::path& config) {
#ifdef _WIN32
    const std::wstring command = L"\"" + executable.wstring() + L"\" \"" + config.wstring() + L"\"";
    std::vector<wchar_t> args(command.begin(), command.end());
    args.push_back(0);
    STARTUPINFOW startup = {};
    startup.cb = sizeof(startup);
    PROCESS_INFORMATION process = {};
    require(CreateProcessW(executable.c_str(), args.data(), NULL, NULL, FALSE,
        CREATE_NO_WINDOW, NULL, NULL, &startup, &process) != 0, "Cannot run movewindow_cli");
    const DWORD status = WaitForSingleObject(process.hProcess, 30000);
    DWORD code = 1;
    GetExitCodeProcess(process.hProcess, &code);
    CloseHandle(process.hThread);
    CloseHandle(process.hProcess);
    require(status == WAIT_OBJECT_0, "movewindow_cli exceeded process timeout");
    return static_cast<int>(code);
#else
    const std::string command = "\"" + executable.string() + "\" \"" + config.string() + "\"";
    return std::system(command.c_str());
#endif
}

fs::path config_file(const fs::path& work, const std::string& name, const fs::path& source,
    const std::string& algorithm, int order, const std::string& sliding,
    const std::string& additional = "", const std::string& target = "Cumulant",
    const std::string& source_type = "source_type=point_set") {
    const fs::path config = work / (name + ".ini");
    std::ofstream out(config);
    out << "[general]\n" << source_type << "\nsource_grid=" << source.string()
        << "\nsource_property=phi\ntarget_property=result\ncalculate_target=" << target
        << "\ncalculate_algorithm=" << algorithm << "\nmissing_value=-999\n"
        << "[sliding_window]\nmoment_cumulant_order=" << order << '\n' << sliding
        << "\n[runtime]\noutput_dir=" << (work / name).string() << '\n' << additional;
    return config;
}

std::vector<Row> run(const fs::path& cli, const fs::path& work, const std::string& name,
    const fs::path& source, const std::string& algorithm, int order, const std::string& sliding,
    const std::string& additional = "", const std::string& target = "Cumulant",
    const std::string& source_type = "source_type=point_set") {
    require(execute(cli, config_file(work, name, source, algorithm, order, sliding, additional,
        target, source_type)) == 0, "Production CLI failed for " + name);
    return read_points(work / name / "result.txt");
}

double expected_cumulant(int order) {
    const std::array<double, 4> values = {1, 2, 4, 8};
    const double mean = 3.75;
    if (order == 1) return mean;
    double second = 0.0, central = 0.0;
    for (double value : values) {
        const double d = value - mean;
        second += d * d / values.size();
        central += std::pow(d, order) / values.size();
    }
    return order == 4 ? central - 3 * second * second : central;
}

void check_coordinates(const std::vector<Row>& input, const std::vector<Row>& actual) {
    require(input.size() == actual.size(), "Every input point must retain an output row");
    for (size_t i = 0; i < input.size(); ++i)
        for (size_t axis = 0; axis < 3; ++axis)
            require(input[i][axis] == actual[i][axis], "Physical coordinate was altered");
}
} // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 2, "Usage: window_pointset_checks <movewindow_cli>");
        const fs::path cli = fs::absolute(argv[1]);
        const fs::path work = fs::absolute("test_output/window_pointset_checks");
        fs::create_directories(work);
        const std::vector<Row> points = {
            {-1.25, 2.5, 3.75, 1}, {-0.5, 2.8, 3.3, 2}, {0.25, 1.9, 3.1, 4},
            {1.0, 2.4, 4.3, 8}, {999.125, 2.5, 3.75, 11}, {0, 2.2, 3.4, nan_value},
            {2000.125, 2.1, 3.2, -999}};
        const fs::path source = work / "points.gslib";
        write_points(source, points);
        for (int order = 1; order <= 4; ++order) {
            for (const std::string algorithm : {"Bell", "AVL Based algorithm"}) {
                for (bool centered : {false, true}) {
                    const std::string name = (algorithm == "Bell" ? "bell" : "avl") +
                        std::to_string(order) + (centered ? "centered" : "raw");
                    const std::string settings = "neighborhood_radii=10,10,10\nmin_samples=2\n"
                        "center_cumulant_input_for_order_ge_2=" + std::string(centered ? "true" : "false");
                    const auto rows = run(cli, work, name, source, algorithm, order, settings);
                    check_coordinates(points, rows);
                    for (size_t i : {size_t(0), size_t(1), size_t(2), size_t(3), size_t(5)})
                        close(rows[i][3], expected_cumulant(order), "Nonzero-mean cumulant at " + name);
                    require(std::isnan(rows[4][3]) && std::isnan(rows[6][3]), "Insufficient neighborhoods must stay NaN");
                    const auto counts = read_points(work / name / "result.sample_count.txt");
                    close(counts[0][3], 4, "Finite and missing sample support");
                    close(counts[4][3], 1, "One-point support");
                    close(counts[6][3], 0, "Empty support");
                    const auto flags = read_points(work / name / "result.valid.txt");
                    close(flags[5][3], 1, "Uninformed center must still be evaluated");
                    close(flags[6][3], 0, "Invalid output flag");
                    const std::string meta = contents(work / name / "result.meta.txt");
                    if (algorithm != "Bell") require(meta.find("evaluator=AVL\n") != std::string::npos,
                        "AVL selection must use the production AVL evaluator");
                    require(contents(work / name / "result.window_support.csv").find("empty_neighborhood") != std::string::npos,
                        "Support diagnostics must explain empty outputs");
                }
            }
        }
        const auto raw_moments = run(cli, work, "raw_moment", source, "AVL", 2,
            "neighborhood_radii=10,10,10\nminimum_samples=2", "", "Moment", "reg_irreg=0");
        close(raw_moments[0][3], 21.25, "Raw moment must retain nonzero mean");
        const auto weighted_avl = run(cli, work, "weighted_avl", source, "Weighted AVL Based algorithm", 4,
            "neighborhood_radii=10,10,10\nminimum_samples=2");
        close(weighted_avl[0][3], expected_cumulant(4), "Uniform weighted AVL result");
        require(contents(work / "weighted_avl" / "result.meta.txt").find("evaluator=AVL\n") != std::string::npos,
            "The weighted qualifier must not replace an explicit AVL evaluator");
        const auto weighted_recurrence = run(cli, work, "weighted_recurrence", source, "Weighted", 4,
            "neighborhood_radii=10,10,10\nminimum_samples=2");
        close(weighted_recurrence[0][3], expected_cumulant(4), "Historical weighted recurrence alias");
        require(contents(work / "weighted_recurrence" / "result.meta.txt").find("evaluator=univariate_recurrence\n") != std::string::npos,
            "The weighted-only evaluator alias must be explicit in metadata");

        // A whole-lattice window is the same empirical sample in both paths.
        // Compare the production executables, including their respective loaders.
        for (int dimensions : {2, 3}) {
            std::vector<Row> lattice;
            const int layers = dimensions == 2 ? 1 : 3;
            for (int z = 0; z < layers; ++z)
                for (int y = 0; y < 3; ++y)
                    for (int x = 0; x < 3; ++x)
                        lattice.push_back({static_cast<double>(x), static_cast<double>(y), static_cast<double>(z),
                            1.0 + x + 2.0 * y * y + 0.5 * z * x + z});
            const fs::path lattice_source = work / ("lattice" + std::to_string(dimensions) + ".gslib");
            write_points(lattice_source, lattice, dimensions);
            for (const std::string algorithm : {"Bell", "AVL"}) {
                const std::string id = algorithm + std::to_string(dimensions);
                const std::string grid_settings = "window_size=3,3," + std::string(dimensions == 2 ? "1" : "3") +
                    "\nmove_step=1,1," + std::string(dimensions == 2 ? "0" : "1") +
                    "\ncenter_cumulant_input_for_order_ge_2=true";
                const auto regular = run(cli, work, "regular_" + id, lattice_source, algorithm, 4,
                    grid_settings, "", "Cumulant", "source_type=cartesian");
                const auto irregular = run(cli, work, "point_lattice_" + id, lattice_source, algorithm, 4,
                    "neighborhood_radii=2,2," + std::string(dimensions == 2 ? "0" : "2") +
                    "\ncenter_cumulant_input_for_order_ge_2=true");
                require(regular.size() == 1, "Full lattice regular window must have one output");
                const size_t center_index = dimensions == 2 ? 4 : 13;
                close(irregular[center_index][3], regular[0][3], "Regular and point-set lattice agreement " + id);
            }
        }

        const fs::path anisotropic = work / "anisotropic.gslib";
        const std::vector<Row> plane = {{-2.25,-3.5,0,nan_value}, {-0.75,-3.5,0,2},
            {-3.75,-3.5,0,6}, {-2.25,-2,0,10}, {-2.25,-5,0,14}, {-0.35,-1.6,0,100}};
        write_points(anisotropic, plane, 2);
        const auto east = run(cli, work, "ellipse_east", anisotropic, "Bell", 1,
            "neighborhood_radii=2,0.5,0\nmin_samples=2");
        const auto north = run(cli, work, "ellipse_north", anisotropic, "AVL", 1,
            "neigh_ellipsoid=2,0.5,0,0,0,0\nmin_samples=2");
        check_coordinates(plane, east);
        close(east[0][3], 4, "X-oriented ellipse membership");
        close(north[0][3], 12, "Azimuth rotation must change physical membership");

        const fs::path vertical = work / "vertical.gslib";
        write_points(vertical, {{0,0,0,nan_value}, {1.5,0,0,2}, {0,0,1.5,10}, {0,1.5,0,50}});
        const auto dipped = run(cli, work, "ellipse_dip", vertical, "Bell", 1,
            "neighborhood_radii=2,0.3,0.3\nneighborhood_angles=90,90,0");
        close(dipped[0][3], 10, "Dip must rotate major axis into Z");

        const fs::path duplicate = work / "duplicate.gslib";
        write_points(duplicate, {{-0.125,0.75,0,1}, {-0.125,0.75,0,9}, {0.875,0.75,0,5}}, 2);
        const auto included = run(cli, work, "center_included", duplicate, "Bell", 1,
            "neighborhood_radii=1,1,0");
        const auto excluded = run(cli, work, "center_excluded", duplicate, "Bell", 1,
            "neighborhood_radii=1,1,0\nincludes_center=false");
        close(included[0][3], 5, "Boundary point must be included");
        close(excluded[0][3], 7, "Exclude center identity but retain coincident observation");
        close(excluded[1][3], 3, "Coincident samples retain independent identities");

        const fs::path uninformed = work / "uninformed.gslib";
        write_points(uninformed, {{-1.5,2.25,3.75,nan_value}, {0.5,2.25,3.75,-999}});
        const auto empty = run(cli, work, "all_uninformed", uninformed, "AVL", 4,
            "neighborhood_radii=10,10,10\nenable_percentile_clipping=true");
        require(empty.size() == 2 && std::isnan(empty[0][3]) && std::isnan(empty[1][3]),
            "An all-uninformed point set must produce NaN rows without clipping crashes");

        for (const std::string mode : {"coordinate_gaussian", "value_gaussian", "soft_high_porosity", "adaptive_upper_tail"}) {
            for (const std::string type : {"point_set", "cartesian"}) {
                const fs::path config = config_file(work, "reject_" + mode + "_" + type, source, "Bell", 4,
                    "neighborhood_radii=10,10,10", "[weighted_moment]\nenabled=false\nweight_mode=" + mode,
                    "Cumulant", "source_type=" + type);
                require(execute(cli, config) != 0, "Removed scheme must not silently fall back to uniform: " + mode);
            }
        }
        for (const std::string settings : {"window_size=3,3,3\nmove_step=1,1,1",
                "neighborhood_radii=0,1,1", "neighborhood_radii=1,1,1\nmin_samples=0",
                "neighborhood_radii=2junk,1,1", "neighborhood_radii=nan,1,1",
                "neighborhood_radii=1,,1", "neighborhood_radii=1,1,1,",
                "neighborhood_radii=1,1,1\nneighborhood_angles=90,zero,0",
                "neighborhood_radii=1,1,1\nneighborhood_angles=90,0,1bad",
                "neigh_ellipsoid=1,1,1,90,0,inf",
                "neighborhood_radii=1,1,1\npercentile_vmin=2percent",
                "neighborhood_radii=1,1,1\npercentile_vmax=NaN"}) {
            require(execute(cli, config_file(work, "reject_geometry", source, "Bell", 4, settings)) != 0,
                "Invalid physical neighborhood settings must be rejected");
        }
        for (const std::string name : {"incorrect", "Bell_typo", "Bell AVL", "algorithm"}) {
            require(execute(cli, config_file(work, "reject_algorithm", source, name, 4,
                "neighborhood_radii=1,1,1")) != 0, "Unknown or ambiguous evaluator must be rejected: " + name);
        }
        require(execute(cli, config_file(work, "reject_target", source, "Bell", 4,
            "neighborhood_radii=1,1,1", "", "Cumulant_typo")) != 0,
            "Unknown statistic must not silently become a cumulant");
        for (const std::string option : {"source_origin=0,0,0", "source_spacing=1,1,1",
                "target_origin=0,0,0", "missing_value=-999oops", "missing_value=nan",
                "coordinate_dimensions=2junk", "coordinate_dimensions=2.5"}) {
            require(execute(cli, config_file(work, "reject_general_option", source, "Bell", 4,
                "neighborhood_radii=1,1,1", "[general]\n" + option)) != 0,
                "Unsupported transform or malformed missing sentinel must be rejected: " + option);
        }
        require(execute(cli, config_file(work, "reject_empty_clipping", uninformed, "Bell", 4,
            "neighborhood_radii=1,1,1\nenable_percentile_clipping=true\npercentile_vmin=2bad")) != 0,
            "Clipping parameters must be validated even when every neighborhood is empty");
        std::cout << "Point-set windows: geometry, support, missing values, Bell/AVL and weight removal checks passed.\n";
        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "window_pointset_checks failed: " << error.what() << '\n';
        return 1;
    }
}
