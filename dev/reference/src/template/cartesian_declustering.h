// Included by main.cpp inside its private namespace: reuses the existing
// Cartesian geometry, output layout, and vendor sampler without redirecting
// regular grids through the point-set candidate search.
struct CartesianTupleCollection {
    const GridData* grid = nullptr;
    std::vector<int> original_slot;
    std::vector<std::vector<double>> values;
    std::vector<std::vector<size_t>> ids;
};

void capture_cartesian_tuple(int variables, const int* x, const int* y,
    const int* z, const double* values, void* context) {
    auto& samples = *static_cast<CartesianTupleCollection*>(context);
    std::vector<double> row(static_cast<size_t>(variables));
    std::vector<size_t> ids(static_cast<size_t>(variables));
    for (int j = 0; j < variables; ++j) {
        const size_t slot = samples.original_slot.empty() ? static_cast<size_t>(j)
            : static_cast<size_t>(samples.original_slot[static_cast<size_t>(j)]);
        ids[slot] = (static_cast<size_t>(z[j]) * samples.grid->ny + y[j]) * samples.grid->nx + x[j];
        // GetGrid2 already subtracts a global mean. Recover original samples
        // by the traced source coordinates, never by undoing that subtraction.
        row[slot] = samples.grid->nodes.at(ids[slot]).value;
    }
    samples.values.push_back(std::move(row));
    samples.ids.push_back(std::move(ids));
}

std::vector<double> cartesian_physical_vector(const IniConfig& config,
    const std::string& key, const std::vector<double>& fallback) {
    const std::string value = get_alias_or(config, "template_shape", {key},
        get_alias_or(config, "general", {key}, ""));
    if (value.empty()) return fallback;
    const auto result = without_sgems::parse_declustering_vector(value, key);
    if (result.size() != 3) throw std::runtime_error(key + " must contain three finite coordinates.");
    return result;
}

std::string cartesian_sequence(const std::vector<int>& values) {
    std::ostringstream out;
    for (size_t j = 0; j < values.size(); ++j) { if (j) out << ';'; out << values[j]; }
    return out.str();
}

int run_cartesian_declustering(const IniConfig& config, const fs::path& config_path,
    const GridData& source, const fs::path& source_path,
    const std::string& target, const std::string& algorithm, const std::string& property,
    const std::vector<double>& target_origin, const std::vector<int>& nod,
    const std::vector<int>& nlag, const std::vector<int>& lag,
    const std::vector<int>& ltol, const std::vector<int>& azm, const std::vector<int>& atol,
    const std::vector<int>& bandwh, const std::vector<int>& dip,
    const std::vector<int>& dtol, const std::vector<int>& bandwd) {
    const auto options = without_sgems::parse_template_weight_options(config);
    if (options.zero_first_order_moment)
        throw std::runtime_error("zero_first_order_moment=true is not supported: Cartesian declustered moments use raw values and cumulants use weighted centering.");
    if (nod.empty() || nod.size() > 8 || nod.size() != nlag.size() || nod.size() != lag.size())
        throw std::runtime_error("Cartesian declustering requires one to eight matching NOD, NLAG and LAG entries.");
    for (const auto* list : {&ltol,&azm,&atol,&bandwh,&dip,&dtol,&bandwd})
        if (list->size()!=nod.size()) throw std::runtime_error("Directional parameter lists must all have the same length.");
    int total_order = 1;
    size_t lag_count = 1;
    for (size_t j = 0; j < nod.size(); ++j) {
        if (nod[j] < 1 || nod[j] > 8 || nlag[j] < 1 || lag[j] < 1)
            throw std::runtime_error("Cartesian NOD, NLAG and LAG must be positive integers.");
        total_order += nod[j];
        if (lag_count > 1000000 / static_cast<size_t>(nlag[j]))
            throw std::runtime_error("Cartesian lag lattice exceeds one million combinations.");
        lag_count *= static_cast<size_t>(nlag[j]);
        if (static_cast<long long>(nlag[j]-1) * lag[j] > std::numeric_limits<int>::max())
            throw std::runtime_error("Cartesian lag distance exceeds the integer range.");
    }
    const bool moment = is_moment_target(target), bell = is_bell_algorithm(algorithm);
    const std::string sampling = normalized_template_sampling(config);
    const bool exhaustive = is_exhaustive_integer_sampling(sampling);
    if (!exhaustive && sampling != "legacy_search" && sampling != "legacy" && sampling != "search")
        throw std::runtime_error("template_sampling must be legacy_search or exhaustive_integer_offsets.");
    for (size_t j=0;j<nod.size();++j) {
        if (ltol[j]<0 || atol[j]<0 || bandwh[j]<0 || dtol[j]<0 || dtol[j]>180 || bandwd[j]<0 ||
            std::abs(static_cast<long long>(dip[j]))>90)
            throw std::runtime_error("Cartesian template tolerances and bandwidths must be nonnegative, DTOL<=180 and |DIP|<=90.");
        if (!exhaustive && atol[j]>=90)
            throw std::runtime_error("legacy_search requires ATOL<90 degrees for its finite rectangular tolerance construction.");
    }
    GridData grid = build_dense_grid_for_template(source, true);
    const int dimensions = grid.nz <= 1 ? 2 : 3;
    auto spacing = cartesian_physical_vector(config, "source_spacing", {1,1,1});
    // Preserve the input lattice coordinate scale unless the caller supplies a
    // physical origin. With spacing alone, the original coordinates are scaled.
    auto origin = cartesian_physical_vector(config, "source_origin",
        {source.raw_min_x * spacing[0], source.raw_min_y * spacing[1], source.raw_min_z * spacing[2]});
    for (size_t j = 0; j < 3; ++j)
        if (!(spacing[j] > 0) || !std::isfinite(spacing[j]) || !std::isfinite(origin[j]))
            throw std::runtime_error("source_spacing must be positive and source_origin must be finite.");
    double missing = std::numeric_limits<double>::quiet_NaN();
    if (config.has("general", "missing_value")) {
        const auto parsed = without_sgems::parse_declustering_vector(config.get("general", "missing_value"), "missing_value");
        if (parsed.size() != 1) throw std::runtime_error("missing_value requires one finite number.");
        missing = parsed[0];
    }
    without_sgems::PointSetData physical;
    physical.title = source.title; physical.value_label = source.value_label; physical.dimensions = dimensions;
    std::vector<unsigned char> informed(grid.nodes.size(), 0);
    std::vector<int> source_ids(grid.nodes.size(), 0);
    for (const auto& node : source.nodes) {
        const size_t id = (static_cast<size_t>(node.z) * grid.ny + node.y) * grid.nx + node.x;
        source_ids[id] = node.index;
    }
    for (size_t id = 0; id < grid.nodes.size(); ++id) {
        auto& node = grid.nodes[id];
        if (!std::isfinite(node.value) || node.value == missing) node.value = std::numeric_limits<double>::quiet_NaN();
        without_sgems::PointSample point;
        point.index = static_cast<int>(id) + 1;
        point.x = origin[0] + node.x * spacing[0];
        point.y = origin[1] + node.y * spacing[1];
        point.z = origin[2] + node.z * spacing[2];
        point.value = node.value; point.informed = std::isfinite(node.value);
        point.external_weight = node.external_weight;
        if (!std::isfinite(point.x) || !std::isfinite(point.y) || !std::isfinite(point.z))
            throw std::runtime_error("Cartesian physical coordinate overflow.");
        informed[id] = point.informed ? 1 : 0;
        physical.nodes.push_back(point);
    }
    // Declustering ignores the inactive Z dimension but the audit retains
    // its actual physical elevation.
    const auto weights = without_sgems::calculate_template_declustering(config, physical);
    const auto layout = build_exhaustive_layout(config, nod.size(), nlag, !exhaustive);
    auto directions = exhaustive ? parse_integer_direction_vectors(config, dimensions, nod.size())
        : std::vector<std::array<int,3>>(nod.size(), {{0,0,0}});
    if (exhaustive) for (size_t j = 0; j < directions.size(); ++j)
        for (int axis = 0; axis < 3; ++axis)
            if (std::abs(static_cast<long double>(directions[j][axis])) * (nlag[j]-1) * lag[j] > std::numeric_limits<int>::max())
                throw std::runtime_error("Integer template offset exceeds the supported range.");
    std::vector<ExhaustiveLagCell> cells;
    std::vector<int> positions(nod.size(), 0);
    enumerate_exhaustive_lag_cells(nlag, lag, directions, layout, 0, &positions, &cells);
    std::sort(cells.begin(), cells.end(), [&](const ExhaustiveLagCell& a, const ExhaustiveLagCell& b) {
        for (int direction : layout.displayed_dirs) {
            const size_t j = static_cast<size_t>(direction-1);
            if (a.positions[j] != b.positions[j]) return a.positions[j] < b.positions[j];
        }
        return false;
    });
    const bool audit = without_sgems::parse_weight_boolean(config.get_or("runtime", "write_template_tuples", "false"), "write_template_tuples");
    const auto minimum_values = without_sgems::parse_declustering_vector(config.get_or("template_shape", "minimum_template_count", "1"), "minimum_template_count");
    if (minimum_values.size() != 1 || minimum_values[0] < 1 || minimum_values[0] != std::floor(minimum_values[0]))
        throw std::runtime_error("minimum_template_count must be a positive integer.");
    const double minimum = minimum_values[0];
    const fs::path output = resolve_path(config_path, get_alias_or(config,"runtime",{"output_dir"},"output/hosc_avl"));
    fs::create_directories(output);
    std::ofstream support(output/(property+".support.csv"));
    std::ofstream tuples;
    support << "lag_cell,lag_positions,lag_distances,anchor_count,support_count,value,status,weight_sum,ESS\n" << std::setprecision(17);
    if (audit) {
        tuples.open(output/(property+".tuples.csv"));
        tuples << "lag_cell,tuple,node_ids,values,tuple_weight\n" << std::setprecision(17);
    }
    Keys* key = nullptr;
    bool initialized = false;
    std::vector<int> reordered_dirs;
    std::vector<without_sgems::GridNode> output_nodes;
    try {
        if (!exhaustive) {
            key = build_keys(grid,nod,nlag,lag,ltol,azm,atol,bandwh,dip,dtol,bandwd);
            std::vector<int> fixed;
            for (int value : layout.fixed_positions_zero_based) fixed.push_back(value+1);
            const auto legacy_layout = apply_cross_section_layout(key, layout.displayed_dirs, fixed);
            reordered_dirs = legacy_layout.hidden_dirs;
            reordered_dirs.insert(reordered_dirs.end(), legacy_layout.displayed_dirs.begin(), legacy_layout.displayed_dirs.end());
            // The compatibility layout helper does not reorder LTOL. Keep
            // each tolerance attached to its configured direction here,
            // without changing the historical unweighted entry point.
            for (size_t j=0;j<reordered_dirs.size();++j)
                key->LTOL[j]=ltol[static_cast<size_t>(reordered_dirs[j]-1)];
            GridData vendor_source = grid;
            // Missing cells are blocked by the explicit informed mask. The
            // compatibility initializer computes a global mean, so finite
            // placeholders keep that internal pass well-defined. These
            // placeholders never enter the trace or the weighted moments.
            for (auto& node : vendor_source.nodes) if (!std::isfinite(node.value)) node.value=0;
            double** vendor_grid = without_sgems::create_vendor_grid_array(vendor_source);
            // GetGrid2 consumes the input matrix; Destruction frees its nodes.
            if (InitialStruct4(key, vendor_grid) != 0) throw std::runtime_error("Cannot initialize the legacy Cartesian template sampler.");
            initialized = true;
        }
        for (size_t row = 0; row < cells.size(); ++row) {
            CartesianTupleCollection samples;
            samples.grid = &grid;
            if (exhaustive) {
                for (size_t anchor = 0; anchor < grid.nodes.size(); ++anchor) {
                    if (!informed[anchor]) continue;
                    const auto& a = grid.nodes[anchor];
                    std::vector<size_t> ids(1,anchor);
                    std::vector<double> values(1,a.value);
                    bool complete = true;
                    for (const auto& offset : cells[row].offsets) {
                        const long long x = static_cast<long long>(a.x)+offset[0], y = static_cast<long long>(a.y)+offset[1], z = static_cast<long long>(a.z)+offset[2];
                        if (x<0 || x>=grid.nx || y<0 || y>=grid.ny || z<0 || z>=grid.nz) {complete=false;break;}
                        const size_t id = (static_cast<size_t>(z)*grid.ny+static_cast<size_t>(y))*grid.nx+static_cast<size_t>(x);
                        if (!informed[id]) {complete=false;break;}
                        ids.push_back(id); values.push_back(grid.nodes[id].value);
                    }
                    if (complete) { samples.ids.push_back(std::move(ids)); samples.values.push_back(std::move(values)); }
                }
            } else {
                samples.original_slot.push_back(0);
                std::vector<int> distances;
                for (int direction : reordered_dirs) {
                    samples.original_slot.push_back(direction);
                    distances.push_back(cells[row].positions[static_cast<size_t>(direction-1)]*lag[static_cast<size_t>(direction-1)]);
                }
                SetTemplateRawTrace(capture_cartesian_tuple, &samples, informed.data(), grid.nx, grid.ny);
                Template* templates = nullptr;
                int count = 0;
                const int status = dimensions==2
                    ? TemplateConstruct_2D(distances.data(),key->ORDER,&templates,&count,0,pt_Grids->ND)
                    : TemplateConstruct_3D(distances.data(),key->ORDER,&templates,&count,0,pt_Grids->ND);
                SetTemplateRawTrace(nullptr,nullptr,nullptr,0,0);
                if (templates) {
                    std::free(templates[0].index);
                    templates[0].index=nullptr;
                    FreeTemplate(&templates);
                    G_Template=nullptr;G_TemplateNum=0;
                }
                if (status < 0 || count != static_cast<int>(samples.values.size()))
                    throw std::runtime_error("Legacy Cartesian raw-template trace is inconsistent.");
            }
            std::vector<double> tuple_weights;
            long double weight_sum=0, weight_square_sum=0;
            for (size_t t=0;t<samples.ids.size();++t) {
                auto distinct = samples.ids[t];
                std::sort(distinct.begin(),distinct.end());
                distinct.erase(std::unique(distinct.begin(),distinct.end()),distinct.end());
                long double log_sum=0;
                for (size_t id:distinct) {
                    const double weight=weights.weights[id];
                    if (!(weight>0) || !std::isfinite(weight)) throw std::runtime_error("Invalid informed point weight in Cartesian template.");
                    log_sum+=std::log(static_cast<long double>(weight));
                }
                const double weight=static_cast<double>(std::exp(log_sum/distinct.size()));
                tuple_weights.push_back(weight); weight_sum+=weight; weight_square_sum+=static_cast<long double>(weight)*weight;
                if (audit) {
                    tuples<<row<<','<<t<<",\"";
                    for (size_t j=0;j<samples.ids[t].size();++j) {if(j)tuples<<';';tuples<<samples.ids[t][j]+1;}
                    tuples<<"\",\"";
                    for (size_t j=0;j<samples.values[t].size();++j) {if(j)tuples<<';';tuples<<samples.values[t][j];}
                    tuples<<"\","<<weight<<'\n';
                }
            }
            const double value = samples.values.size()<minimum ? std::numeric_limits<double>::quiet_NaN()
                : without_sgems::evaluate_complete_template_samples(samples.values,tuple_weights,nod,moment,bell,options);
            if (samples.values.size()>=minimum && !std::isfinite(value)) throw std::runtime_error("Cartesian weighted statistic is not finite.");
            without_sgems::GridNode node;
            node.index=static_cast<int>(row)+1;
            node.x=cells[row].positions[static_cast<size_t>(layout.displayed_dirs[0]-1)]*lag[static_cast<size_t>(layout.displayed_dirs[0]-1)];
            node.y=layout.displayed_dirs.size()>1?cells[row].positions[static_cast<size_t>(layout.displayed_dirs[1]-1)]*lag[static_cast<size_t>(layout.displayed_dirs[1]-1)]:0;
            node.z=layout.displayed_dirs.size()>2?cells[row].positions[static_cast<size_t>(layout.displayed_dirs[2]-1)]*lag[static_cast<size_t>(layout.displayed_dirs[2]-1)]:0;
            node.value=value;output_nodes.push_back(node);
            support<<row<<",\""<<cartesian_sequence(cells[row].positions)<<"\",\"";
            for(size_t j=0;j<nod.size();++j){if(j)support<<';';support<<cells[row].positions[j]*lag[j];}
            support<<"\","<<samples.values.size()<<','<<samples.values.size()<<','<<value<<','<<(std::isfinite(value)?"ok":"insufficient_support")
                <<','<<static_cast<double>(weight_sum)<<','<<(weight_square_sum>0?static_cast<double>(weight_sum*weight_sum/weight_square_sum):0)<<'\n';
        }
        if (initialized) {
            if (dimensions==2) {Free2DVec(&ValueMap_2D,MaxX);ValueMap_2D=nullptr;}
            else {Free3DVec(&ValueMap_3D,MaxX,MaxY);ValueMap_3D=nullptr;}
            Destruction();initialized=false;
        }
        free_keys(key);key=nullptr;
    } catch (...) {
        SetTemplateRawTrace(nullptr,nullptr,nullptr,0,0);
        if (initialized) Destruction();
        free_keys(key);
        throw;
    }
    without_sgems::write_text_grid(output/(property+".txt"),source,source.title,property,output_nodes);
    std::ofstream point_audit(output/(property+".point_weights.csv"));
    point_audit<<"node_id,x,y,z,value,informed,weight,source_id\n"<<std::setprecision(17);
    for(size_t i=0;i<physical.nodes.size();++i){const auto& p=physical.nodes[i];point_audit<<p.index<<','<<p.x<<','<<p.y<<','<<p.z<<','<<p.value<<','<<(p.informed?1:0)<<','<<weights.weights[i]<<','<<source_ids[i]<<'\n';}
    std::ofstream meta(output/(property+".meta.txt"));
    meta<<std::setprecision(17)<<"source_type=cartesian\nsource_grid="<<source_path.string()<<"\nsource_property="<<source.value_label
        <<"\ncalculate_target="<<target<<"\ncalculate_algorithm="<<algorithm<<"\ncoordinate_dimensions="<<dimensions
        <<"\ntemplate_sampling="<<(exhaustive?"exhaustive_integer_offsets":"legacy_search")<<"\nweight_mode="<<without_sgems::template_weight_mode(config)<<"\nweighted_moment_enabled=true"
        <<"\nsource_origin="<<origin[0]<<','<<origin[1]<<','<<origin[2]<<"\nsource_spacing="<<spacing[0]<<','<<spacing[1]<<','<<spacing[2]
        <<"\nsource_coordinate_mapping=source_origin_plus_zero_based_normalized_index_times_source_spacing"
        <<"\nsource_origin_default=raw_input_minimum_index_times_source_spacing\nlag_units=integer_grid_steps"
        <<"\nsource_raw_min="<<source.raw_min_x<<','<<source.raw_min_y<<','<<source.raw_min_z
        <<"\ntarget_origin="<<target_origin[0]<<','<<target_origin[1]<<','<<target_origin[2]
        <<"\noutput_coordinates=selected_direction_lag_distances_in_grid_steps"
        <<"\ntuple_policy="<<(exhaustive?"all_complete_integer_template_translations":"legacy_TemplateConstruct_selection_and_global_exclusion_with_explicit_informed_mask")
        <<"\nmissing_policy=omitted_nonfinite_and_configured_missing_value_excluded\nnode_id=dense_zero_based_linear_index_plus_one"
        <<"\nlegacy_validity_override=explicit_informed_mask_instead_of_attribute_sign_or_hardcoded_minus_one_sentinel"
        <<'\n'<<weights.metadata
        <<"point_weight_normalization=sum_one_over_informed_selected_property_points"
        <<"\ntuple_weight=geometric_mean_over_distinct_source_record_ids_in_complete_target_tuple"
        <<"\nsubmoment_support=same_complete_tuple_collection\nsubmoment_weights=same_fixed_target_tuple_weights_for_all_dependency_moments"
        <<"\nmoment_normalization=divide_by_sum_of_complete_tuple_weights\nmoment_target=raw_product_average"
        <<"\ncumulant_weighted_centering=per_variable_using_fixed_complete_tuple_weights"
        <<"\nevaluator="<<(moment?"weighted_raw_product_average":(bell?"weighted_set_partition_reference":"legacy_AVL_recurrence_with_complete_weighted_moment_cache"))
        <<"\nminimum_template_count="<<minimum<<"\ndisplayed_dirs="<<cartesian_sequence(layout.displayed_dirs)
        <<"\nfixed_lag_positions_zero_based="<<cartesian_sequence(layout.fixed_positions_zero_based)
        <<"\noutput_shape="<<nlag[static_cast<size_t>(layout.displayed_dirs[0]-1)]<<','
        <<(layout.displayed_dirs.size()>1?nlag[static_cast<size_t>(layout.displayed_dirs[1]-1)]:1)<<','
        <<(layout.displayed_dirs.size()>2?nlag[static_cast<size_t>(layout.displayed_dirs[2]-1)]:1)<<'\n';
    for(size_t j=0;j<nod.size();++j) meta<<"direction_"<<j+1<<'='<<nod[j]<<','<<nlag[j]<<','<<lag[j]<<','<<ltol[j]<<','<<azm[j]<<','<<atol[j]<<','<<bandwh[j]<<','<<dip[j]<<','<<dtol[j]<<','<<bandwd[j]<<'\n';
    if (exhaustive) {
        meta<<"integer_direction_vectors=";
        for(size_t j=0;j<directions.size();++j) {
            if(j) meta<<';';
            for(int axis=0;axis<dimensions;++axis) {if(axis)meta<<',';meta<<directions[j][axis];}
        }
        meta<<"\ntemplate_geometry=integer_direction_vectors_times_zero_based_lag_position_times_LAG\n";
    }
    if (!support || !point_audit || !meta || (audit&&!tuples)) throw std::runtime_error("Cannot write Cartesian declustering audit output.");
    std::cout<<"Wrote Cartesian "<<without_sgems::template_weight_mode(config)<<' '<<(exhaustive?"integer-offset":"legacy-search")<<" templates to "<<output.string()<<"\n";
    return 0;
}
