#include "text.h"

#include <cctype>
#include <cstdlib>
#include <stdexcept>

namespace without_sgems {

std::string trim(const std::string& value) {
    size_t start = 0;
    while (start < value.size() && std::isspace(static_cast<unsigned char>(value[start])) != 0) {
        start++;
    }

    size_t end = value.size();
    while (end > start && std::isspace(static_cast<unsigned char>(value[end - 1])) != 0) {
        end--;
    }

    return value.substr(start, end - start);
}

std::string to_lower(std::string value) {
    for (size_t i = 0; i < value.size(); ++i) {
        value[i] = static_cast<char>(std::tolower(static_cast<unsigned char>(value[i])));
    }
    return value;
}

std::vector<std::string> split(const std::string& value, char delimiter) {
    std::vector<std::string> parts;
    std::string current;

    for (size_t i = 0; i < value.size(); ++i) {
        if (value[i] == delimiter) {
            parts.push_back(trim(current));
            current.clear();
            continue;
        }
        current.push_back(value[i]);
    }
    parts.push_back(trim(current));
    return parts;
}

std::vector<int> parse_int_list(const std::string& value) {
    std::vector<std::string> parts = split(value, ',');
    std::vector<int> result;
    result.reserve(parts.size());

    for (size_t i = 0; i < parts.size(); ++i) {
        if (parts[i].empty()) {
            throw std::runtime_error("Encountered an empty integer item in list: " + value);
        }
        result.push_back(std::atoi(parts[i].c_str()));
    }
    return result;
}

std::vector<double> parse_double_list(const std::string& value) {
    std::vector<std::string> parts = split(value, ',');
    std::vector<double> result;
    result.reserve(parts.size());

    for (size_t i = 0; i < parts.size(); ++i) {
        if (parts[i].empty()) {
            throw std::runtime_error("Encountered an empty floating-point item in list: " + value);
        }
        result.push_back(std::atof(parts[i].c_str()));
    }
    return result;
}

}
