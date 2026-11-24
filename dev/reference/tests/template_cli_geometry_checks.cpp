#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

namespace fs = std::filesystem;

namespace {

constexpr int kGridSide = 7;
constexpr double kTolerance = 1.0e-11;

struct DirectionPair {
    const char* id;
    std::array<std::array<int, 2>, 2> original;
    std::array<std::array<int, 2>, 2> rotated;
};

struct SupportRow {
    std::array<int, 2> positions;
    std::array<std::array<int, 3>, 2> offsets;
    int valid_translations;
    int support_count;
    double cumulant;
};

int fail(const std::string& message) {
    std::cerr << "FAILED: " << message << "\n";
    return 1;
}

std::string quote_for_shell(const fs::path& value) {
    return std::string("\"") + value.string() + "\"";
}

double asymmetric_value(int x, int y) {
    // This field deliberately has no 90-degree rotational symmetry.
    const int residue = (3 * x + 5 * y + 2 * x * y) % 7;
    return 0.011 * static_cast<double>(x * x * x) +
        0.023 * static_cast<double>(y * y) +
        0.071 * static_cast<double>(x * y) +
        0.137 * static_cast<double>(x) -
        0.193 * static_cast<double>(y) +
        0.017 * static_cast<double>(residue);
}

void write_grid(const fs::path& path, bool rotate_counterclockwise) {
    std::ofstream output(path.c_str(), std::ios::out | std::ios::trunc);
    if (!output) {
        throw std::runtime_error("Unable to write grid fixture: " + path.string());
    }
    output << "template_cli geometry fixture\n";
    output << std::setw(11) << kGridSide * kGridSide << std::setw(12) << 4 << "\n";
    output << " # FOR EACH OF THE NODES :  \"INDEX\" \"X\" \"Y\" \"Z\" \"VAL\"\n";
    output << " # IF THE DATA ARE IN 2D : Z = 0\n";
    output << std::setprecision(17);

    int index = 1;
    for (int y = 0; y < kGridSide; ++y) {
        for (int x = 0; x < kGridSide; ++x) {
            double value = 0.0;
            if (!rotate_counterclockwise) {
                value = asymmetric_value(x, y);
            }
            else {
                // T(x, y) = (side - 1 - y, x), whose vector action is
                // (dx, dy) -> (-dy, dx).  Therefore G(T(p)) = F(p).
                value = asymmetric_value(y, kGridSide - 1 - x);
            }
            output << std::setw(10) << index++
                   << std::setw(11) << x
                   << std::setw(11) << y
                   << std::setw(11) << 0
                   << " " << std::setw(21) << value << "\n";
        }
    }
}

std::string vector_text(const std::array<int, 2>& vector) {
    return std::to_string(vector[0]) + "," + std::to_string(vector[1]);
}

void write_config(
    const fs::path& path,
    const fs::path& grid_path,
    const fs::path& output_dir,
    const std::array<std::array<int, 2>, 2>& directions,
    const std::string& algorithm,
    const std::string& target_property) {
    std::ofstream output(path.c_str(), std::ios::out | std::ios::trunc);
    if (!output) {
        throw std::runtime_error("Unable to write template CLI configuration: " + path.string());
    }
    output << "[general]\n";
    output << "calculate_target=Cumulant\n";
    output << "calculate_algorithm=" << algorithm << "\n";
    output << "source_grid=" << grid_path.string() << "\n";
    output << "source_property=value\n";
    output << "coordinate_dimensions=2\n";
    output << "target_property=" << target_property << "\n\n";
    output << "[template_shape]\n";
    output << "template_sampling=exhaustive_integer_offsets\n";
    output << "NOD=2,1\n";
    output << "NLAG=2,2\n";
    output << "LAG=1,1\n";
    output << "LTOL=0,0\n";
    output << "AZM=0,45\n";
    output << "ATOL=0,0\n";
    output << "BANDWH=0,0\n";
    output << "DIP=0,0\n";
    output << "DTOL=0,0\n";
    output << "BANDWD=0,0\n";
    output << "integer_direction_vectors=" << vector_text(directions[0]) << ";"
           << vector_text(directions[1]) << "\n";
    output << "source_spacing=50,50,2\n\n";
    output << "[target_cross_sections]\n";
    output << "displayed_dirs=1,2\n";
    output << "fixed_positions=\n\n";
    output << "[runtime]\n";
    output << "output_dir=" << output_dir.string() << "\n";
}

std::vector<std::array<int, 3>> parse_offsets(const std::string& value) {
    static const std::regex offset_pattern("\\((-?[0-9]+),(-?[0-9]+),(-?[0-9]+)\\)");
    std::vector<std::array<int, 3>> offsets;
    for (std::sregex_iterator match(value.begin(), value.end(), offset_pattern), end;
         match != end; ++match) {
        offsets.push_back({{
            std::stoi((*match)[1].str()),
            std::stoi((*match)[2].str()),
            std::stoi((*match)[3].str()),
        }});
    }
    return offsets;
}

std::vector<int> parse_positions(const std::string& value) {
    std::vector<int> values;
    std::istringstream input(value);
    std::string item;
    while (std::getline(input, item, ';')) {
        values.push_back(std::stoi(item));
    }
    return values;
}

std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> fields;
    std::string field;
    bool quoted = false;
    for (size_t index = 0; index < line.size(); ++index) {
        const char character = line[index];
        if (character == '\"') {
            quoted = !quoted;
        }
        else if (character == ',' && !quoted) {
            fields.push_back(field);
            field.clear();
        }
        else {
            field.push_back(character);
        }
    }
    fields.push_back(field);
    return fields;
}

std::vector<SupportRow> read_support_rows(const fs::path& support_path) {
    std::ifstream input(support_path.c_str());
    if (!input) {
        throw std::runtime_error("Missing production support output: " + support_path.string());
    }
    std::string header;
    std::getline(input, header);
    if (header != "lag_cell,lag_positions,integer_offsets,physical_offsets,valid_translations,support_count,cumulant") {
        throw std::runtime_error("Unexpected support CSV header in: " + support_path.string());
    }

    std::vector<SupportRow> rows;
    std::string line;
    while (std::getline(input, line)) {
        const std::vector<std::string> fields = split_csv_line(line);
        if (fields.size() != 7) {
            throw std::runtime_error("Malformed support CSV row in: " + support_path.string());
        }
        const std::vector<int> positions = parse_positions(fields[1]);
        const std::vector<std::array<int, 3>> offsets = parse_offsets(fields[2]);
        if (positions.size() != 2 || offsets.size() != 2) {
            throw std::runtime_error("Support CSV did not retain both signed template offsets.");
        }
        SupportRow row;
        row.positions = {{positions[0], positions[1]}};
        row.offsets = {{offsets[0], offsets[1]}};
        row.valid_translations = std::stoi(fields[4]);
        row.support_count = std::stoi(fields[5]);
        row.cumulant = std::stod(fields[6]);
        rows.push_back(row);
    }
    if (rows.size() != 4) {
        throw std::runtime_error("Expected exactly four lag cells from NLAG=(2,2).");
    }
    return rows;
}

int expected_translation_count(const std::array<std::array<int, 3>, 2>& offsets) {
    int minimum_x = 0;
    int maximum_x = 0;
    int minimum_y = 0;
    int maximum_y = 0;
    for (const std::array<int, 3>& offset : offsets) {
        minimum_x = std::min(minimum_x, offset[0]);
        maximum_x = std::max(maximum_x, offset[0]);
        minimum_y = std::min(minimum_y, offset[1]);
        maximum_y = std::max(maximum_y, offset[1]);
    }
    return (kGridSide - maximum_x + minimum_x) * (kGridSide - maximum_y + minimum_y);
}

bool nearly_equal(double left, double right) {
    return std::abs(left - right) <= kTolerance * std::max({1.0, std::abs(left), std::abs(right)});
}

int run_template_cli(const fs::path& executable, const fs::path& config_path) {
#ifdef _WIN32
    const std::wstring command = L"\"" + executable.wstring() + L"\" \"" + config_path.wstring() + L"\"";
    std::vector<wchar_t> mutable_command(command.begin(), command.end());
    mutable_command.push_back(L'\0');
    STARTUPINFOW startup = {};
    startup.cb = sizeof(startup);
    PROCESS_INFORMATION process = {};
    if (!CreateProcessW(
            executable.c_str(), mutable_command.data(), NULL, NULL, FALSE, 0, NULL, NULL, &startup, &process)) {
        throw std::runtime_error(
            "CreateProcessW failed for template_cli with Windows error " + std::to_string(GetLastError()) + ".");
    }
    WaitForSingleObject(process.hProcess, INFINITE);
    DWORD exit_code = 1;
    const BOOL got_exit_code = GetExitCodeProcess(process.hProcess, &exit_code);
    CloseHandle(process.hThread);
    CloseHandle(process.hProcess);
    if (!got_exit_code || exit_code != 0) {
        throw std::runtime_error(
            "template_cli failed for " + config_path.string() + " with exit code " + std::to_string(exit_code) + ".");
    }
    return static_cast<int>(exit_code);
#else
    const std::string command = quote_for_shell(executable) + " " + quote_for_shell(config_path);
    const int result = std::system(command.c_str());
    if (result != 0) {
        throw std::runtime_error(
            "template_cli failed for " + config_path.string() + " using command " + command + ".");
    }
    return result;
#endif
}

std::array<int, 2> rotate_vector_90_counterclockwise(const std::array<int, 2>& vector) {
    return {{-vector[1], vector[0]}};
}

DirectionPair make_direction_pair(
    const char* id,
    int first_x,
    int first_y,
    int second_x,
    int second_y) {
    DirectionPair result = {};
    result.id = id;
    result.original[0] = {{first_x, first_y}};
    result.original[1] = {{second_x, second_y}};
    result.rotated[0] = rotate_vector_90_counterclockwise(result.original[0]);
    result.rotated[1] = rotate_vector_90_counterclockwise(result.original[1]);
    return result;
}

int validate_bell_avl_rows(
    const std::vector<SupportRow>& avl_rows,
    const std::vector<SupportRow>& bell_rows,
    const std::string& id) {
    if (avl_rows.size() != bell_rows.size()) {
        return fail(id + ": Bell and AVL returned different numbers of lag cells.");
    }
    for (size_t row_index = 0; row_index < avl_rows.size(); ++row_index) {
        const SupportRow& avl = avl_rows[row_index];
        const SupportRow& bell = bell_rows[row_index];
        if (avl.positions != bell.positions ||
            avl.offsets != bell.offsets ||
            avl.valid_translations != bell.valid_translations ||
            avl.support_count != bell.support_count) {
            return fail(id + ": Bell and AVL did not use the same exhaustive template support.");
        }
        if (!nearly_equal(avl.cumulant, bell.cumulant)) {
            std::ostringstream message;
            message << id << ": Bell and AVL disagree at lag cell " << row_index
                    << " (AVL=" << std::setprecision(17) << avl.cumulant
                    << ", Bell=" << bell.cumulant << ").";
            return fail(message.str());
        }
    }
    return 0;
}

int validate_pair(
    const fs::path& executable,
    const fs::path& root,
    const fs::path& original_grid,
    const fs::path& rotated_grid,
    const DirectionPair& direction_pair) {
    const std::string base_id(direction_pair.id);
    const fs::path original_avl_output = root / (base_id + "_original_avl");
    const fs::path original_bell_output = root / (base_id + "_original_bell");
    const fs::path rotated_avl_output = root / (base_id + "_rotated_avl");
    const fs::path original_avl_config = root / (base_id + "_original_avl.ini");
    const fs::path original_bell_config = root / (base_id + "_original_bell.ini");
    const fs::path rotated_avl_config = root / (base_id + "_rotated_avl.ini");
    const std::string original_avl_target = base_id + "_original_avl";
    const std::string original_bell_target = base_id + "_original_bell";
    const std::string rotated_avl_target = base_id + "_rotated_avl";

    write_config(
        original_avl_config, original_grid, original_avl_output, direction_pair.original,
        "AVL Based algorithm", original_avl_target);
    write_config(
        original_bell_config, original_grid, original_bell_output, direction_pair.original,
        "Bell algorithm", original_bell_target);
    write_config(
        rotated_avl_config, rotated_grid, rotated_avl_output, direction_pair.rotated,
        "AVL Based algorithm", rotated_avl_target);
    run_template_cli(executable, original_avl_config);
    run_template_cli(executable, original_bell_config);
    run_template_cli(executable, rotated_avl_config);

    const std::vector<SupportRow> original_rows = read_support_rows(
        original_avl_output / (original_avl_target + ".support.csv"));
    const std::vector<SupportRow> bell_rows = read_support_rows(
        original_bell_output / (original_bell_target + ".support.csv"));
    const std::vector<SupportRow> rotated_rows = read_support_rows(
        rotated_avl_output / (rotated_avl_target + ".support.csv"));
    const int bell_avl_result = validate_bell_avl_rows(original_rows, bell_rows, base_id);
    if (bell_avl_result != 0) {
        return bell_avl_result;
    }

    for (size_t row_index = 0; row_index < original_rows.size(); ++row_index) {
        const SupportRow& original_row = original_rows[row_index];
        const SupportRow& rotated_row = rotated_rows[row_index];
        if (original_row.positions != rotated_row.positions) {
            return fail(std::string(direction_pair.id) + ": lag-cell positions changed after rotation.");
        }

        std::array<std::array<int, 3>, 2> expected_original_offsets = {{
            {{direction_pair.original[0][0] * original_row.positions[0], direction_pair.original[0][1] * original_row.positions[0], 0}},
            {{direction_pair.original[1][0] * original_row.positions[1], direction_pair.original[1][1] * original_row.positions[1], 0}},
        }};
        if (original_row.offsets != expected_original_offsets) {
            return fail(std::string(direction_pair.id) + ": production support output did not preserve the configured signed integer offsets.");
        }
        const int expected_support = expected_translation_count(expected_original_offsets);
        if (original_row.valid_translations != expected_support ||
            original_row.support_count != expected_support ||
            rotated_row.valid_translations != expected_support ||
            rotated_row.support_count != expected_support) {
            return fail(std::string(direction_pair.id) + ": exhaustive sampler did not retain every complete translation.");
        }

        for (size_t direction = 0; direction < original_row.offsets.size(); ++direction) {
            const std::array<int, 2> original_vector = {{
                original_row.offsets[direction][0], original_row.offsets[direction][1],
            }};
            const std::array<int, 2> expected_rotated = rotate_vector_90_counterclockwise(original_vector);
            if (rotated_row.offsets[direction][0] != expected_rotated[0] ||
                rotated_row.offsets[direction][1] != expected_rotated[1] ||
                rotated_row.offsets[direction][2] != 0) {
                return fail(std::string(direction_pair.id) + ": support offsets do not follow (dx,dy)->(-dy,dx).");
            }
        }
        if (!nearly_equal(original_row.cumulant, rotated_row.cumulant)) {
            std::ostringstream message;
            message << direction_pair.id << ": 90-degree rotation changed the production AVL cumulant at lag cell "
                    << row_index << " from " << std::setprecision(17) << original_row.cumulant
                    << " to " << rotated_row.cumulant;
            return fail(message.str());
        }
    }
    return 0;
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        return fail("Usage: template_cli_geometry_checks <template_cli path>");
    }

    try {
        const fs::path executable = fs::absolute(argv[1]);
        if (!fs::exists(executable)) {
            return fail("template_cli executable does not exist: " + executable.string());
        }

        const fs::path root = fs::current_path() / "test_output" / "template_cli_geometry_checks";
        fs::remove_all(root);
        fs::create_directories(root);
        const fs::path original_grid = root / "asymmetric_original.grid";
        const fs::path rotated_grid = root / "asymmetric_rotated_90ccw.grid";
        write_grid(original_grid, false);
        write_grid(rotated_grid, true);

        const std::array<DirectionPair, 4> direction_pairs = {{
            make_direction_pair("directions_0_45", 1, 0, 1, 1),
            make_direction_pair("directions_90_135", 0, 1, -1, 1),
            make_direction_pair("directions_180_225", -1, 0, -1, -1),
            make_direction_pair("directions_270_315", 0, -1, 1, -1),
        }};

        for (const DirectionPair& direction_pair : direction_pairs) {
            const int result = validate_pair(executable, root, original_grid, rotated_grid, direction_pair);
            if (result != 0) {
                return result;
            }
        }
        std::cout << "template_cli exhaustive signed-direction and 90-degree rotation checks: PASS\n";
        return 0;
    }
    catch (const std::exception& error) {
        return fail(error.what());
    }
}
