#include "config.h"

#include "text.h"

#include <fstream>
#include <stdexcept>

namespace without_sgems {

std::string IniConfig::make_storage_key(const std::string& section, const std::string& key) {
    return to_lower(trim(section)) + "." + to_lower(trim(key));
}

IniConfig IniConfig::load(const std::filesystem::path& path) {
    std::ifstream input(path.c_str(), std::ios::in);
    if (!input) {
        throw std::runtime_error("Unable to open configuration file: " + path.string());
    }

    IniConfig config;
    std::string section = "default";
    std::string line;

    while (std::getline(input, line)) {
        std::string normalized = trim(line);
        if (normalized.empty()) {
            continue;
        }
        if (normalized[0] == '#' || normalized[0] == ';') {
            continue;
        }
        if (normalized.size() >= 2 && normalized[0] == '[' && normalized[normalized.size() - 1] == ']') {
            section = trim(normalized.substr(1, normalized.size() - 2));
            continue;
        }

        size_t split_pos = normalized.find('=');
        if (split_pos == std::string::npos) {
            split_pos = normalized.find(':');
        }
        if (split_pos == std::string::npos) {
            throw std::runtime_error("Invalid configuration line: " + line);
        }

        std::string key = trim(normalized.substr(0, split_pos));
        std::string value = trim(normalized.substr(split_pos + 1));
        config.values_[make_storage_key(section, key)] = value;
    }

    return config;
}

bool IniConfig::has(const std::string& section, const std::string& key) const {
    return values_.find(make_storage_key(section, key)) != values_.end();
}

std::string IniConfig::get(const std::string& section, const std::string& key) const {
    std::map<std::string, std::string>::const_iterator it = values_.find(make_storage_key(section, key));
    if (it == values_.end()) {
        throw std::runtime_error("Missing configuration value [" + section + "] " + key);
    }
    return it->second;
}

std::string IniConfig::get_or(const std::string& section, const std::string& key, const std::string& fallback) const {
    std::map<std::string, std::string>::const_iterator it = values_.find(make_storage_key(section, key));
    if (it == values_.end()) {
        return fallback;
    }
    return it->second;
}

}
