#ifndef WITHOUT_SGEMS_CONFIG_H
#define WITHOUT_SGEMS_CONFIG_H

#include <filesystem>
#include <map>
#include <string>

namespace without_sgems {

class IniConfig {
public:
    static IniConfig load(const std::filesystem::path& path);

    bool has(const std::string& section, const std::string& key) const;
    std::string get(const std::string& section, const std::string& key) const;
    std::string get_or(const std::string& section, const std::string& key, const std::string& fallback) const;

private:
    std::map<std::string, std::string> values_;

    static std::string make_storage_key(const std::string& section, const std::string& key);
};

}

#endif
