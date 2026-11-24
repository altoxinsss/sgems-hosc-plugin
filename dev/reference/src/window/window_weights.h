#ifndef WITHOUT_SGEMS_WINDOW_WEIGHTS_H
#define WITHOUT_SGEMS_WINDOW_WEIGHTS_H

#include "template_weight_config.h"

namespace without_sgems {

// Reuse the template module's source-point weight factory. Window samples
// carry these weights directly, without the template geometric-mean step.
struct WindowWeights {
    bool declustered;
    std::vector<double> values;
    std::string metadata;

    WindowWeights(const IniConfig& config, const PointSetData& data)
        : declustered(uses_source_weights(config)), values(data.nodes.size(), 0.0) {
        if (declustered) {
            const auto options = parse_checked_weight_options(config, true);
            if (options.zero_first_order_moment)
                throw std::runtime_error("Weighted window moments use raw values; zero_first_order_moment must be false.");
            auto result = calculate_template_declustering(config, data);
            values = std::move(result.weights);
            metadata = std::move(result.metadata);
            // The shared factory's historical tuple version label is not a
            // description of the window estimator.
            const auto end = metadata.find('\n');
            metadata.erase(0, end + 1);
            metadata = std::string(template_weight_mode(config)=="external" ? "weighting_scheme_version=" : "declustering_scheme_version=") + template_weight_mode(config) +
                ".source_points.window_normalization.v1\n" + metadata;
        } else {
            for (size_t i = 0; i < values.size(); ++i)
                if (data.nodes[i].informed) values[i] = 1.0;
        }
    }

    void write(const std::filesystem::path& directory, const std::string& target,
               const PointSetData& data) const {
        if (!declustered) return;
        std::ofstream out(directory / (target + ".point_weights.csv"));
        out << "node_id,x,y,z,value,informed,weight\n" << std::setprecision(17);
        for (size_t i = 0; i < values.size(); ++i) {
            const auto& p = data.nodes[i];
            out << p.index << ',' << p.x << ',' << p.y << ',' << p.z << ','
                << p.value << ',' << (p.informed ? 1 : 0) << ',' << values[i] << '\n';
        }
        if (!out) throw std::runtime_error("Cannot write window source point weights.");
    }
};

struct WindowWeightSupport {
    long double sum = 0, squares = 0;
    void add(double weight) { sum += weight; squares += static_cast<long double>(weight) * weight; }
    double ess() const { return squares > 0 ? static_cast<double>(sum * sum / squares) : 0; }
};

} // namespace without_sgems
#endif
