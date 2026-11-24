#include "point_set_template.h"
#include "point_set.h"
#include "text.h"
#include "template_weight_config.h"
#include "FunctionCall.h"
#include "WeightedMomentCore.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <new>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace without_sgems {
namespace {
namespace fs = std::filesystem;
constexpr double pi = 3.1415926535897932384626433832795;
using Vec3 = std::array<double, 3>;
namespace wmc = weighted_moment_core;

std::string option(const IniConfig& config, const std::string& section,
    std::initializer_list<const char*> aliases, const std::string& fallback = "") {
    for (const char* key : aliases) if (config.has(section, key)) return config.get(section, key);
    return fallback;
}

double number(const std::string& text, const std::string& name) {
    const std::string clean = trim(text);
    size_t used = 0;
    double value = 0;
    try { value = std::stod(clean, &used); }
    catch (...) { throw std::runtime_error(name + " must contain a finite number."); }
    if (used != clean.size() || !std::isfinite(value))
        throw std::runtime_error(name + " must contain a finite number.");
    return value;
}

int integer(const std::string& text, const std::string& name, int minimum, int maximum) {
    const double value = number(text, name);
    if (value != std::floor(value) || value < minimum || value > maximum)
        throw std::runtime_error(name + " is outside its supported integer range.");
    return static_cast<int>(value);
}

size_t optional_count_limit(const std::string& text, const std::string& name) {
    const std::string clean=trim(text);
    if(clean.empty() || clean.find_first_not_of("0123456789")!=std::string::npos)
        throw std::runtime_error(name+" must be a non-negative integer; zero means unlimited.");
    unsigned long long value=0;
    try { value=std::stoull(clean); }
    catch(const std::exception&) { throw std::runtime_error(name+" exceeds the addressable integer range."); }
    if(value>std::numeric_limits<size_t>::max())
        throw std::runtime_error(name+" exceeds the addressable integer range.");
    return static_cast<size_t>(value);
}

std::vector<double> numbers(const std::string& text, const std::string& name) {
    std::vector<double> values;
    for (const auto& item : split(text, ',')) values.push_back(number(item, name));
    return values;
}

fs::path relative_to(const fs::path& config, const std::string& value) {
    fs::path result(value);
    return result.is_absolute() ? result : config.parent_path() / result;
}

struct Direction {
    int occurrences = 1;
    int count = 1;
    double lag = 1;
    double radial_tolerance = 0;
    double azimuth = 0;
    double azimuth_tolerance = 0;
    double horizontal_band = 0;
    double dip = 0;
    double dip_tolerance = 0;
    double vertical_band = 0;
};

struct Shape {
    std::vector<Direction> directions;
    int dimensions = 0;
    std::string keys_file;
};

std::vector<double> key_line(const std::string& line) {
    std::istringstream stream(line);
    std::string item;
    std::vector<double> values;
    while (stream >> item) values.push_back(number(item, "KEYS"));
    return values;
}

Shape read_shape(const IniConfig& config, const fs::path& config_path) {
    Shape shape;
    const std::string keys = option(config, "template_shape", {"keys_file"});
    if (!keys.empty()) {
        shape.keys_file = relative_to(config_path, keys).string();
        std::ifstream input(shape.keys_file);
        if (!input) throw std::runtime_error("Cannot open template keys_file: " + shape.keys_file);
        std::vector<std::vector<double>> lines;
        std::string line;
        while (std::getline(input, line)) {
            line = trim(line.substr(0, line.find('#')));
            if (!line.empty()) lines.push_back(key_line(line));
        }
        if (lines.size() < 9 || lines[0].size() != 1 || lines[1].size() != 1 ||
            lines[3].size() != 1 || lines[4].size() != 1)
            throw std::runtime_error("Invalid legacy KEYS layout.");
        shape.dimensions = integer(std::to_string(lines[0][0]), "KEYS DIM", 2, 3);
        if (lines[1][0] != 0)
            throw std::runtime_error("Point-set templates require REG_IRREG=0 in keys_file; a regular-grid KEYS file is not silently reinterpreted.");
        const int count = integer(std::to_string(lines[4][0]), "KEYS NDIR", 1, 8);
        if (lines.size() != static_cast<size_t>(8 + count) || lines[5].size() != static_cast<size_t>(count) ||
            lines[6].size() != static_cast<size_t>(count) || lines[7].size() != static_cast<size_t>(count))
            throw std::runtime_error("KEYS direction count does not match parameter rows.");
        if (lines[3][0] != count + 1)
            throw std::runtime_error("Legacy KEYS ORDER must equal NDIR+1; set repeated-variable NOD separately in the INI file.");
        shape.directions.resize(count);
        for (int j = 0; j < count; ++j) {
            if (lines[8 + j].size() != 6) throw std::runtime_error("KEYS direction rows require AZM ATOL BANDWH DIP DTOL BANDWD.");
            Direction& d = shape.directions[j];
            d.count = integer(std::to_string(lines[5][j]), "KEYS NLAG", 1, 100000);
            d.lag = lines[6][j]; d.radial_tolerance = lines[7][j];
            d.azimuth = lines[8+j][0]; d.azimuth_tolerance = lines[8+j][1];
            d.horizontal_band = lines[8+j][2]; d.dip = lines[8+j][3];
            d.dip_tolerance = lines[8+j][4]; d.vertical_band = lines[8+j][5];
        }
    }
    const std::string nod = option(config, "template_shape", {"nod", "number_of_occurrences_in_each_direction"});
    if (keys.empty()) {
        if (nod.empty()) throw std::runtime_error("Point-set template_shape requires NOD or keys_file.");
        shape.directions.resize(numbers(nod, "NOD").size());
    }
    const size_t count = shape.directions.size();
    if (count < 1 || count > 8) throw std::runtime_error("Point-set templates support one to eight directions.");
    auto assign = [&](std::initializer_list<const char*> aliases, const std::string& name,
                      const std::function<void(Direction&, double)>& setter, bool required) {
        const std::string value = option(config, "template_shape", aliases);
        if (value.empty()) {
            if (required && keys.empty()) throw std::runtime_error("Missing point-set template parameter " + name);
            return;
        }
        const auto values = numbers(value, name);
        if (values.size() != count) throw std::runtime_error(name + " must have one entry per direction.");
        for (size_t j = 0; j < count; ++j) setter(shape.directions[j], values[j]);
    };
    assign({"nod", "number_of_occurrences_in_each_direction"}, "NOD", [](Direction& d,double v){
        if (v != std::floor(v) || v < 1 || v > 8) throw std::runtime_error("NOD must be a positive integer <=8."); d.occurrences=static_cast<int>(v);}, false);
    assign({"nlag", "number_of_lags_in_each_direction"}, "NLAG", [](Direction& d,double v){
        if (v != std::floor(v) || v < 1 || v > 100000) throw std::runtime_error("NLAG must be a positive integer <=100000."); d.count=static_cast<int>(v);}, true);
    assign({"lag", "lag_separation_in_each_direction"}, "LAG", [](Direction& d,double v){d.lag=v;}, true);
    assign({"ltol", "lag_tolerance_in_each_direction"}, "LTOL", [](Direction& d,double v){d.radial_tolerance=v;}, false);
    assign({"azm", "azimuth"}, "AZM", [](Direction& d,double v){d.azimuth=v;}, false);
    assign({"atol", "tolerance_of_azimuth"}, "ATOL", [](Direction& d,double v){d.azimuth_tolerance=v;}, false);
    assign({"bandwh", "bandwidth_of_azimuth"}, "BANDWH", [](Direction& d,double v){d.horizontal_band=v;}, false);
    assign({"dip", "dip_angle"}, "DIP", [](Direction& d,double v){d.dip=v;}, false);
    assign({"dtol", "tolerance_of_dip_angle"}, "DTOL", [](Direction& d,double v){d.dip_tolerance=v;}, false);
    assign({"bandwd", "dip_bandwidth"}, "BANDWD", [](Direction& d,double v){d.vertical_band=v;}, false);
    int order = 1;
    for (const auto& d : shape.directions) {
        order += d.occurrences;
        if (d.lag == 0 || d.radial_tolerance < 0 || d.azimuth_tolerance < 0 || d.azimuth_tolerance > 180 ||
            d.horizontal_band < 0 || std::abs(d.dip) > 90 || d.dip_tolerance < 0 || d.dip_tolerance > 180 || d.vertical_band < 0)
            throw std::runtime_error("Require nonzero LAG, nonnegative tolerances/bandwidths, ATOL/DTOL<=180 and |DIP|<=90.");
    }
    return shape;
}

Vec3 displacement(const Direction& d, int position) {
    const double az = d.azimuth*pi/180, dip = d.dip*pi/180, lag = position*d.lag;
    return {{lag*std::cos(dip)*std::cos(az), lag*std::cos(dip)*std::sin(az), lag*std::sin(dip)}};
}

double angle_difference(double a, double b) {
    return std::abs(std::remainder(a-b, 360.0));
}

bool in_direction(const PointSample& anchor, const PointSample& point,
    const Direction& d, int position, int dimensions) {
    const double dx = point.x-anchor.x, dy = point.y-anchor.y, dz = point.z-anchor.z;
    const double distance = std::hypot(std::hypot(dx,dy),dz);
    const double signed_lag = position*d.lag, lag = std::abs(signed_lag);
    const double eps = 64*std::numeric_limits<double>::epsilon()*std::max(1.0, std::max(lag,distance));
    if (distance <= eps || std::abs(distance-lag) > d.radial_tolerance+eps) return false;
    const double az = d.azimuth + (signed_lag < 0 ? 180 : 0);
    const double dip = signed_lag < 0 ? -d.dip : d.dip;
    const double horizontal = std::hypot(dx,dy);
    if (horizontal > eps && std::abs(std::cos(dip*pi/180)) > 1e-12 &&
        angle_difference(std::atan2(dy,dx)*180/pi, az) > d.azimuth_tolerance+1e-10) return false;
    if (dimensions == 3 && std::abs(std::atan2(dz,horizontal)*180/pi-dip) > d.dip_tolerance+1e-10) return false;
    const double azrad=az*pi/180, diprad=dip*pi/180;
    const double across = -dx*std::sin(azrad)+dy*std::cos(azrad);
    const double vertical = -std::sin(diprad)*(dx*std::cos(azrad)+dy*std::sin(azrad))+std::cos(diprad)*dz;
    // Zero bandwidth means no additional band cap; angular tolerance still applies.
    return (d.horizontal_band == 0 || std::abs(across) <= d.horizontal_band+eps) &&
        (d.vertical_band == 0 || std::abs(vertical) <= d.vertical_band+eps);
}

std::vector<size_t> candidate_query(const PointSetData& data, const PointSpatialIndex& spatial,
    const PointSample& anchor, const Direction& direction, int position, const Vec3& target) {
    const double radius=std::abs(position*direction.lag)+direction.radial_tolerance;
    const double eps=64*std::numeric_limits<double>::epsilon()*std::max(1.0,radius);
    const Vec3 center{{anchor.x,anchor.y,anchor.z}};
    // With exact radial and angular targets the accepted region is close to
    // one displacement, so avoid visiting the entire anchor-centered sphere.
    // Huge unnormalized angles retain the general path to avoid assuming an
    // argument-reduction error bound for trigonometric functions.
    const bool exact=direction.radial_tolerance==0 && direction.azimuth_tolerance==0 &&
        (data.dimensions==2 || direction.dip_tolerance==0) && std::abs(direction.azimuth)<=360;
    if(!exact || !std::isfinite(radius+eps)) return spatial.query_radius(center,radius+eps);
    Vec3 lower,upper;
    for(size_t axis=0;axis<3;++axis) {
        const double expected=center[axis]+target[axis];
        if(!std::isfinite(expected)) return spatial.query_radius(center,radius+eps);
        // Cover the unchanged filter's 1e-10 degree angular allowance, its
        // near-pole azimuth exemption, radial epsilon, and coordinate addition
        // / subtraction rounding. This box is only a conservative prefilter.
        const double angular_pad=radius*(8*(1e-10*pi/180+1e-12));
        const double coordinate_pad=256*std::numeric_limits<double>::epsilon()*
            std::max({1.0,radius,std::abs(center[axis]),std::abs(expected)});
        lower[axis]=expected-angular_pad-coordinate_pad;
        upper[axis]=expected+angular_pad+coordinate_pad;
    }
    auto ids=spatial.query_box(lower,upper);
    // Preserve query_radius's exact three-argument hypot cutoff as well as
    // the subsequent in_direction test; neither support nor order changes.
    ids.erase(std::remove_if(ids.begin(),ids.end(),[&](size_t id){
        const auto& point=data.nodes[id];
        if(point.x<center[0]-(radius+eps) || point.x>center[0]+(radius+eps) ||
           point.y<center[1]-(radius+eps) || point.y>center[1]+(radius+eps) ||
           point.z<center[2]-(radius+eps) || point.z>center[2]+(radius+eps)) return true;
        return std::hypot(point.x-anchor.x,point.y-anchor.y,point.z-anchor.z)>radius+eps;
    }),ids.end());
    return ids;
}

struct Samples {
    std::vector<std::vector<double>> values;
    std::vector<std::vector<size_t>> node_ids;
    std::vector<double> tuple_weights;
    size_t anchor_count = 0;
};

Samples collect(const PointSetData& data, const PointSpatialIndex& spatial,
    const std::vector<Direction>& directions, const std::vector<int>& positions,
    size_t max_tuples, size_t max_visits, bool record_ids) {
    Samples result;
    const size_t variables=directions.size()+1;
    std::vector<Vec3> targets(variables, {{0,0,0}});
    std::vector<int> same(variables,-1);
    for (size_t j=1;j<variables;++j) {
        targets[j]=displacement(directions[j-1],positions[j-1]);
        for(size_t k=0;k<j;++k) {
            const double scale=std::max({1.0,std::abs(targets[j][0]),std::abs(targets[j][1]),std::abs(targets[j][2])});
            if(std::hypot(std::hypot(targets[j][0]-targets[k][0],targets[j][1]-targets[k][1]),targets[j][2]-targets[k][2])<=64*std::numeric_limits<double>::epsilon()*scale){same[j]=static_cast<int>(k);break;}
        }
    }
    size_t visits=0;
    for(size_t anchor=0;anchor<data.nodes.size();++anchor) {
        if(!data.nodes[anchor].informed) continue;
        std::vector<std::vector<size_t>> candidates(variables);
        candidates[0].push_back(anchor);
        bool possible=true;
        for(size_t j=1;j<variables;++j) {
            if(positions[j-1]==0) candidates[j].push_back(anchor);
            else {
                const Direction& d=directions[j-1];
                const auto& a=data.nodes[anchor];
                for(size_t id:candidate_query(data,spatial,a,d,positions[j-1],targets[j]))
                    if(id!=anchor && in_direction(a,data.nodes[id],d,positions[j-1],data.dimensions)) candidates[j].push_back(id);
            }
            if(candidates[j].empty()){possible=false;break;}
        }
        if(!possible) continue;
        std::vector<size_t> selected(variables,anchor);
        const size_t before=result.values.size();
        std::function<void(size_t)> visit = [&](size_t j) {
            if(max_visits) {
                if(visits>=max_visits) throw std::runtime_error("Point-set template candidate enumeration exceeds max_candidate_visits; no tuples were silently truncated.");
                ++visits;
            }
            if(j==variables) {
                if(max_tuples && result.values.size()>=max_tuples) throw std::runtime_error("Point-set template exceeds max_template_tuples; no tuples were silently truncated.");
                // Existing Bell/AVL and weighted-core APIs use int counts.
                if(result.values.size()>=static_cast<size_t>(std::numeric_limits<int>::max()))
                    throw std::runtime_error("Point-set template count exceeds the evaluator integer range.");
                std::vector<double> values(variables);
                for(size_t k=0;k<variables;++k) values[k]=data.nodes[selected[k]].value;
                result.values.push_back(std::move(values));
                if(record_ids) result.node_ids.push_back(selected);
                return;
            }
            if(same[j]>=0) {
                const size_t id=selected[static_cast<size_t>(same[j])];
                if(std::binary_search(candidates[j].begin(),candidates[j].end(),id)) {selected[j]=id;visit(j+1);}
                return;
            }
            for(size_t id:candidates[j]) {
                if(std::find(selected.begin(),selected.begin()+j,id)!=selected.begin()+j) continue;
                selected[j]=id;visit(j+1);
            }
        };
        visit(1);
        if(result.values.size()>before) ++result.anchor_count;
    }
    return result;
}

void assign_tuple_weights(Samples& samples, const std::vector<double>& point_weights) {
    samples.tuple_weights.reserve(samples.node_ids.size());
    for (const auto& ids : samples.node_ids) {
        // NOD counts variable occurrences, not sampled locations. Coincident
        // target slots reuse one source record and therefore one point weight.
        std::vector<size_t> distinct(ids);
        std::sort(distinct.begin(), distinct.end());
        distinct.erase(std::unique(distinct.begin(), distinct.end()), distinct.end());
        long double log_sum = 0;
        for (size_t id : distinct) {
            const double weight = point_weights.at(id);
            if (!(weight > 0) || !std::isfinite(weight))
                throw std::runtime_error("A complete informed template contains an invalid declustering weight.");
            log_sum += std::log(static_cast<long double>(weight));
        }
        const double weight = static_cast<double>(std::exp(log_sum / distinct.size()));
        if (!(weight > 0) || !std::isfinite(weight))
            throw std::runtime_error("Template geometric-mean declustering weight is not finite and positive.");
        samples.tuple_weights.push_back(weight);
    }
    if (samples.tuple_weights.size() != samples.values.size())
        throw std::runtime_error("Template declustering weights do not match the complete-tuple collection.");
}

std::vector<wmc::WeightedSample> weighted_views(const Samples& samples) {
    std::vector<wmc::WeightedSample> views(samples.values.size());
    for (size_t i = 0; i < views.size(); ++i) {
        views[i].values = samples.values[i].data();
        views[i].value_count = static_cast<int>(samples.values[i].size());
        views[i].weight = samples.tuple_weights.empty() ? 1.0 : samples.tuple_weights[i];
    }
    return views;
}

// Fill the complete multiplicity down-set with the SAME tuple weights. The
// legacy evaluator keeps its AVL dependency reuse and recurrence unchanged.
void precompute_weighted_moments(const std::vector<wmc::WeightedSample>& samples,
    const std::vector<int>& index, const wmc::WeightedMomentOptions& options) {
    const int variables = static_cast<int>(index.size()) - 1;
    size_t count = 1;
    for (int j = 1; j <= variables; ++j) count *= static_cast<size_t>(index[j] + 1);
    G_SeriesNum = static_cast<int>(count);
    G_Moment = static_cast<double*>(std::calloc(count, sizeof(double)));
    G_MomentIndex = static_cast<char**>(std::calloc(count, sizeof(char*)));
    if (!G_Moment || !G_MomentIndex) throw std::bad_alloc();
    auto store = [&](size_t position, const std::vector<int>& powers) {
        G_Moment[position] = wmc::weighted_moment(samples.data(), static_cast<int>(samples.size()),
            powers.data() + 1, variables, options);
        if (!std::isfinite(G_Moment[position])) throw std::runtime_error("Weighted template moment overflowed.");
        G_MomentIndex[position] = static_cast<char*>(std::calloc(static_cast<size_t>(variables) + 1, sizeof(char)));
        if (!G_MomentIndex[position]) throw std::bad_alloc();
        int last = variables;
        while (last > 0 && powers[last] == 0) --last;
        for (int j = 1; j <= last; ++j) G_MomentIndex[position][j - 1] = static_cast<char>('0' + powers[j]);
    };
    // RootMomentOrDirect in the legacy AVL route expects the target at slot 0.
    store(0, index);
    for (size_t code = 0; code + 1 < count; ++code) {
        size_t remainder = code;
        std::vector<int> powers(static_cast<size_t>(variables) + 1, 0);
        for (int j = variables; j >= 1; --j) {
            powers[j] = static_cast<int>(remainder % (index[j] + 1));
            powers[0] += powers[j];
            remainder /= static_cast<size_t>(index[j] + 1);
        }
        store(code + 1, powers);
    }
    G_CurrMomentPos = G_SeriesNum;
    G_MomentFull = true;
    G_MomentHashReady = false;
    if (BuildSuppliedMomentHashForDLL() != 0)
        throw std::runtime_error("Cannot build weighted template moment lookup index.");
}

class EvaluationWorkspace {
public:
    explicit EvaluationWorkspace(int variables) : variables_(variables) {
        ClearOneEvaluationWorkspaceForDLL(variables_);
    }
    ~EvaluationWorkspace() {
        SetLowOrderCumulantMomentCacheMode(0);
        ResetFactorSwitchesForRun();
        ClearOneEvaluationWorkspaceForDLL(variables_);
        G_Template = nullptr;
        G_TemplateNum = 0;
    }
private:
    int variables_;
};

double evaluate(Samples& samples, const std::vector<Direction>& directions, bool moment, bool bell,
    const wmc::WeightedMomentOptions& options) {
    const int variables=static_cast<int>(directions.size()+1);
    if(samples.values.empty()) return std::numeric_limits<double>::quiet_NaN();
    const bool weighted = options.mode == wmc::WeightMode::External;
    std::vector<int> index(variables+1,1);
    index[0]=1;
    for(int j=1;j<variables;++j){index[j+1]=directions[j-1].occurrences;index[0]+=index[j+1];}
    if (moment && weighted) {
        const auto views = weighted_views(samples);
        return wmc::weighted_moment(views.data(), static_cast<int>(views.size()), index.data() + 1, variables, options);
    }
    if(moment) {
        long double sum=0;
        for(const auto& row:samples.values) {
            long double product=1;
            for(int j=0;j<variables;++j) for(int k=0;k<index[j+1];++k) product*=row[j];
            sum+=product;
        }
        return static_cast<double>(sum/samples.values.size());
    }
    // Center each participating variable over exactly these complete tuples,
    // using their fixed normalized contribution when declustering is enabled.
    // Cumulants of total order >=2 are invariant to these translations.
    for(int j=0;j<variables;++j) {
        long double sum=0, weight_sum=0;
        for(size_t i=0;i<samples.values.size();++i) {
            const long double weight=weighted?samples.tuple_weights[i]:1.0L;
            sum+=weight*samples.values[i][j];weight_sum+=weight;
        }
        const long double mean=sum/weight_sum;
        for(auto& row:samples.values) row[j]=static_cast<double>(row[j]-mean);
    }
    if (weighted && bell) {
        const auto views = weighted_views(samples);
        return wmc::weighted_cumulant_from_multiplicities(views.data(), static_cast<int>(views.size()),
            index.data() + 1, variables, options);
    }
    std::vector<std::vector<double>> values(samples.values.size(),std::vector<double>(variables+1,0));
    std::vector<Template> views(samples.values.size());
    for(size_t i=0;i<views.size();++i) {
        std::copy(samples.values[i].begin(),samples.values[i].end(),values[i].begin()+1);
        views[i].Value=values[i].data();views[i].Lag=nullptr;views[i].Mean=0;views[i].index=index.data();
    }
    EvaluationWorkspace workspace(variables);
    G_Template=views.data();G_TemplateNum=static_cast<int>(views.size());
    double result=myEmpty;
    if(bell) result=BellCumulantFromMultiplicityIndex(index.data(),variables);
    else {
        SetFactorSwitchesForRun(1,1,1,1);
        SetLowOrderCumulantMomentCacheMode(1);
        if (weighted) precompute_weighted_moments(weighted_views(samples), index, options);
        if((weighted || PrecomputeMomentsForIndex(index.data(),variables)!=myEmpty) && InitAVLWorkspaceForDLL(variables)==0)
            result=BetterRecCumulant3(index.data(),variables);
    }
    if(result==myEmpty || !std::isfinite(result)) throw std::runtime_error("Point-set template Bell/AVL evaluation failed or overflowed.");
    return result;
}

std::string sequence(const std::vector<int>& values) {
    std::ostringstream out;
    for(size_t i=0;i<values.size();++i){if(i)out<<';';out<<values[i];}
    return out.str();
}
}

double evaluate_complete_template_samples(std::vector<std::vector<double>>& values,
    const std::vector<double>& tuple_weights, const std::vector<int>& nod,
    bool moment, bool bell, const wmc::WeightedMomentOptions& options) {
    Samples samples;
    samples.values.swap(values);
    samples.tuple_weights = tuple_weights;
    std::vector<Direction> directions(nod.size());
    for (size_t j = 0; j < nod.size(); ++j) directions[j].occurrences = nod[j];
    try {
        const double result = evaluate(samples, directions, moment, bell, options);
        samples.values.swap(values);
        return result;
    } catch (...) {
        samples.values.swap(values);
        throw;
    }
}

int run_point_set_template(const IniConfig& config, const fs::path& config_path) {
    const auto weights=parse_template_weight_options(config);
    const bool decluster=uses_source_weights(config);
    for (const char* key : {"source_origin", "source_spacing"})
        if (config.has("general", key) || config.has("template_shape", key))
            throw std::runtime_error(std::string("Point-set source coordinates are already physical; remove ") + key + ".");
    if(weights.zero_first_order_moment)
        throw std::runtime_error("zero_first_order_moment=true is not supported by point-set templates; cumulants use translation-invariant centering and moments use raw samples.");
    const std::string sampling=to_lower(option(config,"template_shape",{"template_sampling","sampling_mode"},"point_set_tolerance"));
    if(sampling!="point_set_tolerance" && sampling!="pointset_tolerance")
        throw std::runtime_error("Point-set input requires template_sampling=point_set_tolerance.");
    const Shape shape=read_shape(config,config_path);
    const fs::path source=relative_to(config_path,option(config,"general",{"source_grid","source_grid_file"}));
    const std::string property=option(config,"general",{"source_property"});
    const int dimensions=integer(option(config,"general",{"coordinate_dimensions"},std::to_string(shape.dimensions)),"coordinate_dimensions",0,3);
    if(dimensions==1) throw std::runtime_error("Point-set template coordinate_dimensions must be 2 or 3 (or 0 for automatic).");
    if(shape.dimensions && dimensions && shape.dimensions!=dimensions) throw std::runtime_error("INI coordinate_dimensions conflicts with keys_file DIM.");
    const std::string missing=option(config,"general",{"missing_value"});
    const PointSetData data=load_point_set(source,property,dimensions,
        missing.empty()?std::numeric_limits<double>::quiet_NaN():number(missing,"missing_value"), external_weight_property(config));
    if(data.dimensions < 2) throw std::runtime_error("Point-set templates require two or three coordinate dimensions.");
    if(data.dimensions==2) for(const auto& d:shape.directions) if(d.dip!=0) throw std::runtime_error("Two-dimensional point sets require DIP=0.");
    TemplateDeclusteringResult declustering;
    if(decluster) {
        declustering=calculate_template_declustering(config,data);
    }
    const PointSpatialIndex spatial(data);
    const std::string target=option(config,"general",{"target_property","props"},"template_cumulant");
    const std::string calculation=to_lower(option(config,"general",{"calculate_target","target_button"},"cumulant"));
    const bool moment=calculation=="moment" || calculation=="moments";
    if(!moment && calculation!="cumulant" && calculation!="cumulants") throw std::runtime_error("calculate_target must be Moment or Cumulant.");
    const std::string algorithm=to_lower(option(config,"general",{"calculate_algorithm","algorithm_button"},"avl"));
    const bool bell=algorithm.find("bell")!=std::string::npos;
    if(!bell && algorithm.find("avl")==std::string::npos && algorithm.find("hosc")==std::string::npos)
        throw std::runtime_error("Point-set template algorithm must identify Bell or AVL/HOSC.");
    const size_t max_tuples=optional_count_limit(option(config,"template_shape",{"max_template_tuples"},"0"),"max_template_tuples");
    const size_t max_visits=optional_count_limit(option(config,"template_shape",{"max_candidate_visits"},"0"),"max_candidate_visits");
    const size_t minimum=integer(option(config,"template_shape",{"minimum_template_count"},"1"),"minimum_template_count",1,10000000);
    const bool audit=parse_weight_boolean(option(config,"runtime",{"write_template_tuples"},"false"),"write_template_tuples");
    std::vector<double> origin;
    if(config.has("general","target_origin")) origin=numbers(config.get("general","target_origin"),"target_origin");
    else origin={number(option(config,"general",{"ox","xmn_"},"0"),"ox"),
                 number(option(config,"general",{"oy","ymn_"},"0"),"oy"),
                 number(option(config,"general",{"oz","zmn_"},"0"),"oz")};
    if(origin.size()!=3) throw std::runtime_error("target_origin must contain exactly three finite values.");
    const size_t directions=shape.directions.size();
    std::vector<int> displayed;
    const std::string display=option(config,"target_cross_sections",{"displayed_dirs"});
    if(!display.empty()) for(double v:numbers(display,"displayed_dirs")) {
        if(v!=std::floor(v)||v<1||v>directions) throw std::runtime_error("displayed_dirs contains an invalid direction.");
        displayed.push_back(static_cast<int>(v)-1);
    }
    else if(config.has("target_cross_sections","direction_x")||config.has("target_cross_sections","fst_dir")) {
        for(const std::string& text:std::vector<std::string>{
                option(config,"target_cross_sections",{"direction_x","fst_dir"},"0"),
                option(config,"target_cross_sections",{"direction_y","snd_dir"},"0"),
                option(config,"target_cross_sections",{"direction_z","trd_dir"},"0")}) {
            const int direction=integer(text,"display direction",0,static_cast<int>(directions));
            if(direction) displayed.push_back(direction-1);
        }
    }
    else for(size_t j=0;j<std::min<size_t>(3,directions);++j) displayed.push_back(static_cast<int>(j));
    if(displayed.empty() || displayed.size()>3 || std::set<int>(displayed.begin(),displayed.end()).size()!=displayed.size())
        throw std::runtime_error("displayed_dirs requires one to three distinct directions.");
    std::vector<int> fixed(directions,-1);
    const std::string fixed_text=option(config,"target_cross_sections",{"fixed_positions","position","cross_pos"});
    const auto fixed_values=fixed_text.empty()?std::vector<double>():numbers(fixed_text,"fixed_positions");
    if(fixed_values.size()!=directions-displayed.size()) throw std::runtime_error("Provide one one-based fixed_positions entry per hidden direction.");
    size_t f=0,total=1;
    for(size_t j=0;j<directions;++j) {
        if(std::find(displayed.begin(),displayed.end(),j)==displayed.end()) {
            const double v=fixed_values[f++];
            if(v!=std::floor(v)||v<1||v>shape.directions[j].count) throw std::runtime_error("fixed_positions must be one-based valid lag-bin positions.");
            fixed[j]=static_cast<int>(v)-1;
        }
        else {
            total*=static_cast<size_t>(shape.directions[j].count);
            if(total>100000) throw std::runtime_error("Point-set template request exceeds 100000 output lag cells.");
        }
    }
    std::vector<std::vector<int>> positions;
    std::vector<int> current(directions,0);
    std::function<void(size_t)> enumerate=[&](size_t j){
        if(j==directions){positions.push_back(current);return;}
        const int start=fixed[j]<0?0:fixed[j], stop=fixed[j]<0?shape.directions[j].count:fixed[j]+1;
        for(int k=start;k<stop;++k){current[j]=k;enumerate(j+1);}
    };
    enumerate(0);
    const fs::path output=relative_to(config_path,option(config,"runtime",{"output_dir"},"output/template_pointset"));
    fs::create_directories(output);
    std::ostringstream support,tuples;
    support<<"lag_cell,lag_positions,lag_distances,anchor_count,support_count,value,status";
    if(decluster) support<<",weight_sum,ESS";
    support<<'\n'<<std::setprecision(17);
    if(audit) {
        tuples<<"lag_cell,tuple,node_ids,values";
        if(decluster) tuples<<",tuple_weight";
        tuples<<'\n'<<std::setprecision(17);
    }
    std::vector<PointSample> output_nodes;
    for(size_t row=0;row<positions.size();++row) {
        Samples samples=collect(data,spatial,shape.directions,positions[row],max_tuples,max_visits,audit||decluster);
        if(decluster) assign_tuple_weights(samples,declustering.weights);
        wmc::WeightDiagnostics diagnostics;
        if(decluster && !samples.values.empty()) {
            const auto views=weighted_views(samples);
            diagnostics=wmc::compute_weight_diagnostics(views.data(),static_cast<int>(views.size()),weights);
        }
        if(audit) for(size_t n=0;n<samples.values.size();++n) {
            tuples<<row<<','<<n<<",\"";
            for(size_t j=0;j<samples.node_ids[n].size();++j){if(j)tuples<<';';tuples<<data.nodes[samples.node_ids[n][j]].index;}
            tuples<<"\",\"";
            for(size_t j=0;j<samples.values[n].size();++j){if(j)tuples<<';';tuples<<samples.values[n][j];}
            tuples<<'"';
            if(decluster) tuples<<','<<samples.tuple_weights[n];
            tuples<<'\n';
        }
        const double value=samples.values.size()<minimum?std::numeric_limits<double>::quiet_NaN():evaluate(samples,shape.directions,moment,bell,weights);
        if(samples.values.size()>=minimum && !std::isfinite(value)) throw std::runtime_error("Point-set template statistic overflowed.");
        PointSample out;out.index=static_cast<int>(row)+1;out.value=value;out.informed=std::isfinite(value);
        out.x=origin[0]+positions[row][displayed[0]]*shape.directions[displayed[0]].lag;
        out.y=origin[1];out.z=origin[2];
        if(displayed.size()>1)out.y+=positions[row][displayed[1]]*shape.directions[displayed[1]].lag;
        if(displayed.size()>2)out.z+=positions[row][displayed[2]]*shape.directions[displayed[2]].lag;
        output_nodes.push_back(out);
        support<<row<<",\""<<sequence(positions[row])<<"\",\"";
        for(size_t j=0;j<directions;++j){if(j)support<<';';support<<positions[row][j]*shape.directions[j].lag;}
        support<<"\","<<samples.anchor_count<<','<<samples.values.size()<<','<<value<<','<<(out.informed?"ok":"insufficient_support");
        if(decluster) support<<','<<diagnostics.weight_sum<<','<<diagnostics.effective_sample_size;
        support<<'\n';
    }
    write_point_set(output/(target+".txt"),data.title,target,output_nodes);
    std::ofstream support_file(output/(target+".support.csv"));support_file<<support.str();
    if(audit){std::ofstream file(output/(target+".tuples.csv"));file<<tuples.str();if(!file)throw std::runtime_error("Cannot write template tuple audit.");}
    if(decluster) {
        std::ofstream file(output/(target+".point_weights.csv"));
        file<<"node_id,x,y,z,value,informed,weight\n"<<std::setprecision(17);
        for(size_t i=0;i<data.nodes.size();++i) {
            const auto& point=data.nodes[i];
            file<<point.index<<','<<point.x<<','<<point.y<<','<<point.z<<','<<point.value<<','
                <<(point.informed?1:0)<<','<<declustering.weights[i]<<'\n';
        }
        if(!file) throw std::runtime_error("Cannot write point declustering weight audit.");
    }
    std::ofstream meta(output/(target+".meta.txt"));
    meta<<std::setprecision(17)<<"source_type=point_set\nsource_grid="<<source.string()<<"\nsource_property="<<property
        <<"\ncoordinate_dimensions="<<data.dimensions<<"\ntemplate_sampling=point_set_tolerance\ncalculate_target="<<calculation
        <<"\ncalculate_algorithm="<<algorithm<<"\nweight_mode="<<template_weight_mode(config)<<"\nkeys_file="<<shape.keys_file
        <<"\ncoordinate_units=source_physical_coordinates\noutput_coordinates=target_origin_plus_selected_direction_signed_lag_distances"
        <<"\ntarget_origin="<<origin[0]<<','<<origin[1]<<','<<origin[2]
        <<"\nazimuth_convention=degrees_counterclockwise_from_positive_X_toward_positive_Y"
        <<"\ndip_convention=degrees_above_XY_toward_positive_Z\nlag_tolerance=absolute_euclidean_distance_difference"
        <<"\nbandwidth_convention=nonnegative_half_widths_zero_disables_additional_band_limit"
        <<"\ntuple_policy=all_complete_distinct_node_combinations_per_anchor_in_input_order"
        <<"\ncoincident_target_policy=reuse_same_selected_sample\nzero_lag_policy=anchor_sample_only"
        <<"\nnod_policy=repeated_occurrences_of_same_selected_value\nsubmoment_support=same_complete_tuple_collection"
        <<"\ncumulant_centering=per_variable_over_complete_tuples\nmax_template_tuples="<<max_tuples
        <<"\nmax_candidate_visits="<<max_visits<<"\ncount_limit_zero_means=unlimited\nminimum_template_count="<<minimum<<"\n";
    if(decluster) {
        meta<<declustering.metadata
            <<"point_weight_normalization=sum_one_over_informed_selected_property_points"
            <<"\ntuple_weight=geometric_mean_over_distinct_source_record_ids_in_complete_target_tuple"
            <<"\nsubmoment_weights=same_fixed_target_tuple_weights_for_all_dependency_moments"
            <<"\nmoment_normalization=divide_by_sum_of_complete_tuple_weights"
            <<"\ncumulant_weighted_centering=per_variable_using_fixed_complete_tuple_weights"
            <<"\nevaluator="<<(moment?"weighted_raw_product_average":(bell?"weighted_set_partition_reference":"legacy_AVL_recurrence_with_complete_weighted_moment_cache"))
            <<"\nmoment_cache_lookup="<<(!moment&&!bell?"existing_HOSC_hash_index":"not_applicable")<<'\n';
    }
    for(size_t j=0;j<directions;++j){const auto& d=shape.directions[j];meta<<"direction_"<<j+1<<"="<<d.occurrences<<','<<d.count<<','<<d.lag<<','<<d.radial_tolerance<<','<<d.azimuth<<','<<d.azimuth_tolerance<<','<<d.horizontal_band<<','<<d.dip<<','<<d.dip_tolerance<<','<<d.vertical_band<<'\n';}
    std::vector<int> displayed_one_based=displayed;
    for(int& value:displayed_one_based)++value;
    meta<<"displayed_dirs="<<sequence(displayed_one_based)<<"\nfixed_lag_positions_zero_based="<<sequence(fixed)<<'\n';
    if(!support_file||!meta) throw std::runtime_error("Cannot write point-set template output metadata.");
    std::cout<<"Wrote point-set tolerance templates to "<<output.string()<<"\n";
    return 0;
}
}
