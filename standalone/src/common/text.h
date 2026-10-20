#ifndef WITHOUT_SGEMS_TEXT_H
#define WITHOUT_SGEMS_TEXT_H

#include <string>
#include <vector>

namespace without_sgems {

std::string trim(const std::string& value);
std::string to_lower(std::string value);
std::vector<std::string> split(const std::string& value, char delimiter);
std::vector<int> parse_int_list(const std::string& value);
std::vector<double> parse_double_list(const std::string& value);

}

#endif
