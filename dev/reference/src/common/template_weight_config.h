#ifndef WITHOUT_SGEMS_TEMPLATE_WEIGHT_CONFIG_H
#define WITHOUT_SGEMS_TEMPLATE_WEIGHT_CONFIG_H

#include "cell_declustering.h"
#include "spatial_declustering.h"
#include "weight_config.h"

#include <cmath>
#include <limits>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>

namespace without_sgems {

inline std::string template_weight_mode(const IniConfig& config) {
    return to_lower(trim(config.get_or("weighted_moment", "weight_mode", "uniform")));
}

inline bool uses_source_weights(const IniConfig& config) {
    const auto mode=template_weight_mode(config);
    return mode=="cell_declustering" || mode=="polygonal_declustering" || mode=="kriging_declustering" || mode=="external";
}

inline std::string external_weight_property(const IniConfig& config) {
    return template_weight_mode(config)=="external" ? trim(config.get_or("weighted_moment", "weight_property", "")) : "";
}

inline weighted_moment_core::WeightedMomentOptions parse_template_weight_options(const IniConfig& config) {
    return parse_checked_weight_options(config, true);
}

inline std::vector<double> parse_declustering_vector(const std::string& text, const std::string& key) {
    const auto parts = split(text, ',');
    std::vector<double> values;
    for (const auto& part : parts) {
        std::size_t used = 0;
        double value = 0;
        try {
            value = std::stod(part, &used);
        } catch (const std::exception&) {
            throw std::runtime_error("Invalid numeric value in declustering." + key + ": " + part);
        }
        if (used != part.size() || !std::isfinite(value)) {
            throw std::runtime_error("declustering." + key + " requires finite numeric values.");
        }
        values.push_back(value);
    }
    return values;
}

inline CellDeclusteringOptions parse_cell_declustering_options(const IniConfig& config, int dimensions) {
    if (dimensions != 2 && dimensions != 3) {
        throw std::runtime_error("Template cell declustering requires two- or three-dimensional coordinates.");
    }
    CellDeclusteringOptions options;
    if (config.has("declustering", "dimensions")) {
        const auto parsed = parse_declustering_vector(config.get("declustering", "dimensions"), "dimensions");
        if (parsed.size() != 1 || (parsed[0] != 2 && parsed[0] != 3) || parsed[0] > dimensions)
            throw std::runtime_error("declustering.dimensions must be 2 or 3, no larger than input dimensions.");
        options.dimensions = static_cast<int>(parsed[0]);
        dimensions = options.dimensions;
    }
    const auto sizes = parse_declustering_vector(config.get("declustering", "cell_size"), "cell_size");
    if (sizes.size() != static_cast<std::size_t>(dimensions) && !(dimensions == 2 && sizes.size() == 3)) {
        throw std::runtime_error("declustering.cell_size must specify each active axis (X,Y or X,Y,Z).");
    }
    for (std::size_t axis = 0; axis < sizes.size(); ++axis) {
        if (!(sizes[axis] > 0)) throw std::runtime_error("declustering.cell_size must be positive.");
        options.cell_size[axis] = sizes[axis];
    }
    if (config.has("declustering", "origin")) {
        const auto origin = parse_declustering_vector(config.get("declustering", "origin"), "origin");
        if (origin.size() != static_cast<std::size_t>(dimensions) && !(dimensions == 2 && origin.size() == 3)) {
            throw std::runtime_error("declustering.origin must specify each active axis (X,Y or X,Y,Z).");
        }
        if (dimensions == 2 && origin.size() == 3 && origin[2] != 0) {
            throw std::runtime_error("Two-dimensional declustering.origin requires Z=0 when supplied.");
        }
        for (std::size_t axis = 0; axis < origin.size(); ++axis) options.origin[axis] = origin[axis];
        options.origin_at_data_min = false;
    }
    const std::string shifts = trim(config.get_or("declustering", "shifts_per_axis", "3"));
    std::size_t used = 0;
    long count = 0;
    try {
        count = std::stol(shifts, &used);
    } catch (const std::exception&) {
        throw std::runtime_error("declustering.shifts_per_axis must be a positive integer.");
    }
    if (used != shifts.size() || count < 1 || count > std::numeric_limits<int>::max()) {
        throw std::runtime_error("declustering.shifts_per_axis must be a positive integer.");
    }
    options.shifts_per_axis = static_cast<int>(count);
    if (config.has("declustering", "origin_offsets_file")) {
        const std::filesystem::path path(config.get("declustering", "origin_offsets_file"));
        if (!path.is_absolute()) throw std::runtime_error("declustering.origin_offsets_file requires an absolute path.");
        std::ifstream input(path);
        if (!input) throw std::runtime_error("Cannot open declustering.origin_offsets_file.");
        std::string line;
        while (std::getline(input, line)) {
            line = trim(line);
            if (line.empty() || line[0] == '#') continue;
            const auto values = parse_declustering_vector(line, "origin_offsets_file");
            if (values.size() != static_cast<std::size_t>(dimensions))
                throw std::runtime_error("Each origin offset row must have one fraction per declustering dimension.");
            std::array<double,3> offset{{0,0,0}};
            for (int axis=0; axis<dimensions; ++axis) offset[axis]=values[axis];
            options.origin_offsets.push_back(offset);
        }
        if (options.origin_offsets.empty()) throw std::runtime_error("Origin offsets file is empty.");
    }
    return options;
}

inline DeclusteringDomain parse_declustering_domain(const IniConfig& config, int dimensions) {
    DeclusteringDomain domain;
    domain.dimensions=dimensions;
    if (config.has("declustering","dimensions")) {
        const auto value=parse_declustering_vector(config.get("declustering","dimensions"),"dimensions");
        if (value.size()!=1 || (value[0]!=2 && value[0]!=3) || value[0]>dimensions)
            throw std::runtime_error("declustering.dimensions must be 2 or 3 and no greater than input dimensions.");
        domain.dimensions=static_cast<int>(value[0]);
    }
    for (const std::string key:{"domain_min","domain_max"}) {
        const auto values=parse_declustering_vector(config.get("declustering",key),key);
        if (values.size()!=static_cast<std::size_t>(domain.dimensions))
            throw std::runtime_error("declustering."+key+" must specify exactly one bound per active dimension.");
        auto& target=key=="domain_min"?domain.lower:domain.upper;
        for (std::size_t k=0;k<values.size();++k) target[k]=values[k];
    }
    return domain;
}

// Both template entry points use this factory and identical metadata. This
// selects the point weights only; tuple collection/evaluation is unchanged.
struct TemplateDeclusteringResult {
    std::vector<double> weights;
    std::string metadata;
};

inline TemplateDeclusteringResult calculate_template_declustering(const IniConfig& config, const PointSetData& data) {
    const std::string mode=template_weight_mode(config);
    TemplateDeclusteringResult output;
    std::ostringstream meta;
    meta<<std::setprecision(17)<<(mode=="external"?"weighting_scheme_version=":"declustering_scheme_version=")<<mode<<".distinct_tuple.v1\n";
    if (mode=="external") {
        output.weights.assign(data.nodes.size(), 0.0);
        size_t count=0;
        for (size_t i=0; i<data.nodes.size(); ++i) {
            const auto& p=data.nodes[i];
            if (!p.informed) continue;
            if (!std::isfinite(p.external_weight) || !(p.external_weight>0))
                throw std::runtime_error("External weight must be finite and strictly positive at source record " + std::to_string(p.index));
            output.weights[i]=p.external_weight;
            ++count;
        }
        meta<<"weight_source=external_observation_property\nweight_property="<<external_weight_property(config)
            <<"\nweight_source_file="<<config.get_or("general","source_grid","")
            <<"\nweight_record_mapping=source_row\nweight_valid_count="<<count<<'\n';
    } else if (mode=="cell_declustering") {
        const auto options=parse_cell_declustering_options(config,data.dimensions);
        auto result=calculate_cell_declustering(data,options);
        meta<<"declustering_cell_size="<<options.cell_size[0]<<','<<options.cell_size[1]<<','<<options.cell_size[2]
            <<"\ndeclustering_origin="<<result.origin[0]<<','<<result.origin[1]<<','<<result.origin[2]
            <<"\ndeclustering_origin_at_data_min="<<(options.origin_at_data_min?"true":"false")
            <<"\ndeclustering_shifts_per_axis="<<options.shifts_per_axis
            <<"\ndeclustering_origin_count="<<result.origin_count
            <<"\ndeclustering_dimensions="<<(options.dimensions?options.dimensions:data.dimensions)
            <<"\ndeclustering_shift_schedule="<<(options.origin_offsets.empty()?"regular_cartesian_product":"explicit_fractional_origin_list")
            <<"\ndeclustering_origin_offsets_file="<<config.get_or("declustering","origin_offsets_file","")
            <<"\ndeclustering_valid_count="<<result.valid_count
            <<"\ndeclustering_min_occupied_cells="<<result.min_occupied_cells
            <<"\ndeclustering_max_occupied_cells="<<result.max_occupied_cells
            <<"\ndeclustering_weighted_mean="<<result.weighted_mean<<'\n';
        output.weights=std::move(result.weights);
    } else if (mode=="polygonal_declustering" || mode=="kriging_declustering") {
        const auto domain=parse_declustering_domain(config,data.dimensions);
        meta<<"declustering_dimensions="<<domain.dimensions<<"\ndeclustering_domain_min=";
        for (int k=0;k<domain.dimensions;++k) {if(k)meta<<',';meta<<domain.lower[k];}
        meta<<"\ndeclustering_domain_max=";
        for (int k=0;k<domain.dimensions;++k) {if(k)meta<<',';meta<<domain.upper[k];}
        meta<<"\ncoincident_site_policy=equal_share_of_site_mass_over_informed_records_at_identical_active_coordinates\n";
        SpatialDeclusteringResult result;
        if (mode=="polygonal_declustering") {
            result=calculate_polygonal_declustering(data,domain);
            meta<<"polygonal_method=euclidean_voronoi_exact_halfspace_clipping_to_axis_aligned_domain\n"
                <<"polygonal_measure="<<(domain.dimensions==2?"area":"volume")
                <<"\npolygonal_partition_relative_error="<<result.partition_relative_error<<'\n';
        } else {
            KrigingDeclusteringOptions options;
            options.domain=domain;
            options.covariance_model=to_lower(trim(config.get("declustering","kriging_model")));
            const auto ranges=parse_declustering_vector(config.get("declustering","kriging_ranges"),"kriging_ranges");
            const auto nugget=parse_declustering_vector(config.get("declustering","kriging_nugget_fraction"),"kriging_nugget_fraction");
            const auto cells=parse_declustering_vector(config.get("declustering","kriging_integration_cells"),"kriging_integration_cells");
            if (ranges.size()!=static_cast<std::size_t>(domain.dimensions) || cells.size()!=ranges.size() || nugget.size()!=1)
                throw std::runtime_error("Kriging requires one range/integration count per active axis and one nugget fraction.");
            options.nugget_fraction=nugget[0];
            for (std::size_t k=0;k<ranges.size();++k) {
                if (cells[k]<1 || cells[k]!=std::floor(cells[k]) || cells[k]>std::numeric_limits<int>::max())
                    throw std::runtime_error("Kriging integration counts must be positive integers representable by the integration index type.");
                options.ranges[k]=ranges[k];
                options.integration_cells[k]=static_cast<int>(cells[k]);
            }
            result=calculate_kriging_declustering(data,options);
            meta<<"kriging_target=ordinary_block_mean_over_specified_domain\nkriging_equations=Cw_plus_lambda_one_equals_domain_average_covariance_sum_w_equals_one\n"
                <<"kriging_model="<<options.covariance_model<<"\nkriging_ranges=";
            for (int k=0;k<domain.dimensions;++k) {if(k)meta<<',';meta<<options.ranges[k];}
            meta<<"\nkriging_nugget_fraction="<<options.nugget_fraction<<"\nkriging_total_sill=1\nkriging_integration_cells=";
            for (int k=0;k<domain.dimensions;++k) {if(k)meta<<',';meta<<options.integration_cells[k];}
            meta<<"\nkriging_integration=equal_volume_cell_midpoints\nkriging_integration_point_count="<<result.integration_point_count
                <<"\nkriging_range_definition=spherical_support_exponential_exp_minus_3h_gaussian_exp_minus_3h_squared"
                <<"\nkriging_anisotropy=axis_aligned_ranges\nkriging_nugget_integration=zero_continuous_domain_measure"
                <<"\nkriging_lagrange_multiplier="<<result.kriging_lagrange_multiplier
                <<"\nkriging_max_residual="<<result.kriging_max_residual
                <<"\nkriging_min_cholesky_pivot="<<result.kriging_min_cholesky_pivot
                <<"\nkriging_nonpositive_policy=error_no_clipping_no_absolute_value\nkriging_diagonal_jitter=none\n";
        }
        meta<<"declustering_valid_count="<<result.valid_count
            <<"\ndeclustering_unique_site_count="<<result.unique_site_count
            <<"\ndeclustering_min_site_weight="<<result.min_site_weight
            <<"\ndeclustering_max_site_weight="<<result.max_site_weight
            <<"\ndeclustering_weighted_mean="<<result.weighted_mean<<'\n';
        output.weights=std::move(result.weights);
    } else throw std::runtime_error("No declustering method selected.");
    output.metadata=meta.str();
    return output;
}

} // namespace without_sgems

#endif
