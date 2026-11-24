#ifndef WITHOUT_SGEMS_EXTERNAL_WEIGHT_IO_H
#define WITHOUT_SGEMS_EXTERNAL_WEIGHT_IO_H
#include "text.h"
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace without_sgems {
inline int external_weight_column(const std::vector<std::string>& labels,
    const std::string& property, int first_property) {
    if (property.empty()) return -1;
    int found = -1;
    for (int i = first_property; i < static_cast<int>(labels.size()); ++i) {
        if (to_lower(trim(labels[i])) != to_lower(trim(property))) continue;
        if (found >= 0) throw std::runtime_error("Duplicate weight_property column: " + property);
        found = i;
    }
    if (found < 0) throw std::runtime_error("weight_property column not found: " + property);
    return found;
}

// Defer validity to the source-weight factory, which knows which observations
// are informed (including a configured missing-value sentinel).
inline double external_weight_token(const std::string& token) {
    try {
        size_t used = 0;
        const double result = std::stod(token, &used);
        if (used == token.size()) return result;
    } catch (const std::exception&) {}
    return std::numeric_limits<double>::quiet_NaN();
}
}
#endif
