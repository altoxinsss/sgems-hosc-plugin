#ifndef HOSC_ENGINE_PLUGIN_H
#define HOSC_ENGINE_PLUGIN_H

#include <GsTLAppli/geostat/geostat_algo.h>
#include <string>
#include <map>
#include <vector>
class Geostat_grid;
class GsTLGridProperty;
class Cartesian_grid;

namespace hosc_sgems {
// The x64 SGeMS adapter and numerical core use the host-compatible VC10 ABI.
// Numerical work runs inside this DLL, with no helper executable.
class EnginePlugin : public Geostat_algo {
public:
    explicit EnginePlugin(bool is_template);
    virtual ~EnginePlugin();
    virtual bool initialize(const Parameters_handler*, Error_messages_handler*);
    virtual int execute(GsTL_project* project=0);
    virtual std::string name() const;
private:
    EnginePlugin(const EnginePlugin&);
    EnginePlugin& operator=(const EnginePlugin&);
    std::string get(const std::string&,const std::string& fallback="") const;
    std::string effective_weight_mode() const;
    void configure();
    void export_source();
    void run_engine(const std::string& config="config.ini", const std::string& log="engine.log");
    void run_standardized();
    void import_result(GsTL_project*);
    void import_slice(GsTL_project*,Geostat_grid*,const std::string&);
    bool template_,point_input_,source_is_points_,point_regular_output_;
    Geostat_grid* source_;
    Cartesian_grid* cartesian_;
    GsTLGridProperty* property_;
    std::map<std::string,std::string> parameters_;
    std::string property_name_,source_name_,run_dir_,ini_;
    int source_dims_[3],output_dims_[3],window_[3],step_[3];
    double source_origin_[3],source_spacing_[3],output_origin_[3],output_spacing_[3],file_origin_[3];
};
}
#endif
