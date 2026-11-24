#define NOMINMAX
#include "EnginePlugin.h"
#include "../native/NativeCore.h"
#include "../native/NativeTypes.h"
#include "../native/WindowGeometry.h"
#include "ValidationFeedback.h"
#include "Qt4ParameterCompat.h"
#include "RunPerformance.h"
#include <GsTLAppli/geostat/parameters_handler.h>
#include <GsTLAppli/grid/grid_model/cartesian_grid.h>
#include <GsTLAppli/grid/grid_model/point_set.h>
#include <GsTLAppli/grid/grid_model/grid_property.h>
#include <GsTLAppli/appli/manager_repository.h>
#include <GsTLAppli/appli/project.h>
#include <GsTLAppli/utils/error_messages_handler.h>
#include <windows.h>
#ifdef interface
#undef interface
#endif
#include <float.h>
#include <climits>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <algorithm>
#include <limits>

namespace hosc_sgems {
namespace {
struct Row {double xyz[3];double value;};
std::string number(double x){std::ostringstream s;s<<std::setprecision(17)<<x;return s.str();}
std::string lower(std::string x){for(size_t i=0;i<x.size();++i)x[i]=(char)tolower((unsigned char)x[i]);return x;}
bool flag(const std::string& s){const std::string v=lower(s);if(v=="true"||v=="1"||v=="enable"||v=="yes")return true;if(v=="false"||v=="0"||v=="disable"||v=="no"||v.empty())return false;throw std::runtime_error("Invalid Boolean: "+s);}
double value(const std::string& s){const std::string v=lower(s);if(v.find("nan")!=std::string::npos||v.find("inf")!=std::string::npos)return std::numeric_limits<double>::quiet_NaN();char* end=0;double x=strtod(s.c_str(),&end);if(end==s.c_str()||*end)throw std::runtime_error("Invalid number: "+s);return x;}
std::vector<double> numbers(std::string s){for(size_t i=0;i<s.size();++i)if(s[i]==','||s[i]==';')s[i]=' ';std::istringstream f(s);std::string t;std::vector<double> v;while(f>>t){double x=value(t);if(!_finite(x))throw std::runtime_error("A finite number is required.");v.push_back(x);}return v;}
int integer(const std::string& s,int minimum,int maximum=INT_MAX){double x=value(s);if(!_finite(x)||x!=floor(x)||x<minimum||x>maximum)throw std::runtime_error("Invalid integer: "+s);return (int)x;}
bool exists(const std::string& s){return GetFileAttributesA(s.c_str())!=INVALID_FILE_ATTRIBUTES;}
std::string parent(const std::string& s){return s.substr(0,s.find_last_of("/\\"));}
std::string escaped(std::string s){std::string out;for(size_t i=0;i<s.size();++i){switch(s[i]){case '&':out+="&amp;";break;case '<':out+="&lt;";break;case '>':out+="&gt;";break;case '"':out+="&quot;";break;default:out+=s[i];}}return out;}
std::string contents(const std::string& p){std::ifstream f(p.c_str());std::ostringstream s;s<<f.rdbuf();return s.str();}
std::string join3(const double* p){return number(p[0])+","+number(p[1])+","+number(p[2]);}
std::string new_run_directory(const std::string& requested){char temp[MAX_PATH],shortpath[MAX_PATH];std::string base=requested;if(base.empty()){if(!GetTempPathA(MAX_PATH,temp))throw std::runtime_error("Cannot locate the temporary directory.");base=temp;}if(!exists(base))throw std::runtime_error("Audit directory must already exist.");if(GetShortPathNameA(base.c_str(),shortpath,MAX_PATH))base=shortpath;for(size_t i=0;i<base.size();++i)if((unsigned char)base[i]>=128)throw std::runtime_error("Use an ASCII audit directory (for example D:\\HOSC_runs) when Windows short paths are unavailable.");for(int attempt=0;attempt<1000;++attempt){std::string path=base+"\\HOSC_"+number(GetCurrentProcessId())+"_"+number(GetTickCount())+"_"+number(attempt);if(CreateDirectoryA(path.c_str(),0))return path;}throw std::runtime_error("Cannot create a unique HOSC working directory.");}
std::vector<Row> read_rows(const std::string& path){std::ifstream in(path.c_str());if(!in)throw std::runtime_error("Missing engine result: "+path);std::string line;for(int j=0;j<6;++j)if(!std::getline(in,line))throw std::runtime_error("Incomplete GSLIB result.");std::vector<Row> out;while(std::getline(in,line)){std::istringstream s(line);std::string t[4];if(!(s>>t[0]))continue;if(!(s>>t[1]>>t[2]>>t[3]))throw std::runtime_error("Invalid engine result row.");Row r;for(int k=0;k<3;++k){r.xyz[k]=value(t[k]);if(!_finite(r.xyz[k]))throw std::runtime_error("Non-finite output location.");}r.value=value(t[3]);out.push_back(r);}return out;}
std::vector<std::string> csv(const std::string& line){std::vector<std::string> out;std::string cell;bool quote=false;for(size_t i=0;i<line.size();++i){if(line[i]=='"')quote=!quote;else if(line[i]==','&&!quote){out.push_back(cell);cell.clear();}else if(line[i]!='\r')cell+=line[i];}out.push_back(cell);return out;}
std::string unique_name(const std::string& base){for(int i=0;i<100000;++i){std::string n=base+(i?"_"+number(i):"");if(!Root::instance()->interface(gridModels_manager+"/"+n).raw_ptr())return n;}throw std::runtime_error("Cannot allocate a unique output name.");}
void store(GsTLGridProperty* p,int i,double v){if(_finite(v)){if(fabs(v)>FLT_MAX)throw std::runtime_error("Result exceeds SGeMS float property range; double results are in the audit directory.");p->set_value((float)v,i);}else p->set_not_informed(i);}
Geostat_grid* create_grid(GsTL_project* project,const std::string& base,const std::vector<Row>& rows,const int* dims,const double* origin,const double* spacing,bool points,std::string& name){
    name=unique_name(base);std::string path=gridModels_manager+"/"+name;
    std::string type=points?"point_set://"+number((double)rows.size()):"cgrid";
    SmartPtr<Named_interface> ni=Root::instance()->new_interface(type,path);Geostat_grid* grid=dynamic_cast<Geostat_grid*>(ni.raw_ptr());
    if(!grid)throw std::runtime_error("SGeMS could not create the output grid.");
    if(points){Point_set* ps=dynamic_cast<Point_set*>(grid);if(!ps)throw std::runtime_error("Output point-set type mismatch.");std::vector<Geostat_grid::location_type> locations;for(size_t i=0;i<rows.size();++i)locations.push_back(Geostat_grid::location_type(rows[i].xyz[0],rows[i].xyz[1],rows[i].xyz[2]));ps->point_locations(locations);}
    else {Cartesian_grid* cg=dynamic_cast<Cartesian_grid*>(grid);if(!cg)throw std::runtime_error("Output Cartesian type mismatch.");cg->set_dimensions(dims[0],dims[1],dims[2],spacing[0],spacing[1],spacing[2]);cg->origin(GsTLPoint(origin[0],origin[1],origin[2]));cg->set_level(1);}
    project->new_object(name);
    return dynamic_cast<Geostat_grid*>(Root::instance()->interface(path).raw_ptr());
}
}

EnginePlugin::EnginePlugin(bool t):template_(t),point_input_(false),source_is_points_(false),point_regular_output_(false),source_(0),cartesian_(0),property_(0){hosc_qt4::install();}
EnginePlugin::~EnginePlugin(){}
std::string EnginePlugin::name()const{return template_?"HOSC_TEMPLATE":"HOSC_WINDOW";}
std::string EnginePlugin::get(const std::string& k,const std::string& fallback)const{std::map<std::string,std::string>::const_iterator i=parameters_.find(k);return i==parameters_.end()||i->second.empty()?fallback:i->second;}
std::string EnginePlugin::effective_weight_mode()const{
    // Older template files have no switch: their selected mode still applies.
    if(template_&&!flag(get("enable_weighted_moment","true")))return "uniform";
    return get("weight_mode","uniform")=="External weights" ? "external" : get("weight_mode","uniform");
}
bool EnginePlugin::initialize(const Parameters_handler* p,Error_messages_handler* errors){
    hosc_validation::FeedbackGuard feedback(errors,template_?"hosc_template_validation_summary":"hosc_window_validation_summary");
    try {
        parameters_.clear();source_=0;cartesian_=0;property_=0;
        const char* fields[]={"Source_Grid","Source_Property","props","Target_Button","Algorithm_Button","template_sampling","NOD","NLAG","LAG","LTOL","AZM","ATOL","BANDWH","DIP","DTOL","BANDWD","fst_dir","snd_dir","trd_dir","cross_pos","integer_direction_vectors","xmn_","ymn_","zmn_","weight_mode","weight_property","enable_weighted_moment","normalize_weights","zero_first_order_moment","coordinate_dimensions","declustering_dimensions","cell_size","declustering_origin","shifts_per_axis","origin_offsets_file","domain_min","domain_max","kriging_model","kriging_ranges","kriging_nugget_fraction","kriging_integration_cells","minimum_template_count","max_template_tuples","max_candidate_visits","moment_order","windowl","windoww","windowh","stepx","stepy","stepz","center_mode","neighborhood_shape","grid_bounds_mode","grid_bounds_min","grid_bounds_max","grid_start","grid_spacing","neighborhood_radii","neighborhood_angles","min_samples","includes_center","center_cumulant_input_for_order_ge_2","invert_window_output","enable_percentile_clipping","percentile_vmin","percentile_vmax","output_standardized_cumulant","variance_epsilon","output_diagnostics","output_point_weights","audit_output_directory","Output_3D_Slice_Enable","Output_3D_Slice_Source_Grid","Output_3D_Slice_Source_Property","Output_3D_Slice_Plane","Output_3D_Slice_Position"};
        for(size_t i=0;i<sizeof(fields)/sizeof(fields[0]);++i){std::string v=p->value(std::string(fields[i])+".value");if(v.find_first_of("\r\n")!=std::string::npos)throw std::runtime_error("Parameter values must be single-line text.");parameters_[fields[i]]=v;}
        const char* removed[]={"enable_coordinate_gaussian","enable_fixed_value_gaussian","enable_value_gaussian","enable_soft_high_porosity","enable_adaptive_upper_tail","enable_joint_geometric_mean"};
        for(size_t i=0;i<sizeof(removed)/sizeof(removed[0]);++i){std::string v=p->value(std::string(removed[i])+".value");if(!v.empty()&&flag(v))throw std::runtime_error(std::string("Removed custom weighting switch: ")+removed[i]);}
        source_name_=get("Source_Grid");property_name_=get("props","HOSC_result");
        if(property_name_.find_first_of("/\\:\r\n")!=std::string::npos)throw std::runtime_error("Output property must not contain path separators.");
        source_=dynamic_cast<Geostat_grid*>(Root::instance()->interface(gridModels_manager+"/"+source_name_).raw_ptr());
        if(!source_||source_->size()<1)throw std::runtime_error("Select a nonempty Cartesian grid or point set.");
        cartesian_=dynamic_cast<Cartesian_grid*>(source_);source_is_points_=dynamic_cast<Point_set*>(source_)!=0;
        if(!cartesian_&&!source_is_points_)throw std::runtime_error("Only Cartesian grids and point sets are supported.");
        property_=source_->property(get("Source_Property"));if(!property_)throw std::runtime_error("Select an existing source property.");
        if(cartesian_&&cartesian_->current_level()!=1)throw std::runtime_error("Set the source grid to multigrid level 1.");
        configure();return true;
    }catch(const std::exception& e){errors->report("Source_Grid",escaped(e.what()));return false;}
}

void EnginePlugin::configure(){
    std::string target=get("Target_Button","Cumulant");if(target!="Cumulant"&&target!="Moment")throw std::runtime_error("Calculate target must be Cumulant or Moment.");
    std::string algorithm=get("Algorithm_Button","AVL Based algorithm");if(algorithm!="Bell"&&algorithm!="AVL Based algorithm"&&algorithm!="AVL")throw std::runtime_error("Select Bell or AVL Based algorithm.");
    std::string method=get("weight_mode","uniform");
    if(method=="External weights")method="external";
    if(method!="uniform"&&method!="cell_declustering"&&method!="polygonal_declustering"&&method!="kriging_declustering"&&method!="external")throw std::runtime_error("Removed custom weighting mode. Select uniform, a declustering method or External weights.");
    method=effective_weight_mode();
    if(!get("normalize_weights").empty()&&!flag(get("normalize_weights")))throw std::runtime_error("Weights are always normalized.");
    if(flag(get("zero_first_order_moment","false")))throw std::runtime_error("Raw moments use original property values.");
    std::string sampling=get("template_sampling","auto");
    if(sampling=="auto")sampling=source_is_points_?"point_set_tolerance":"legacy_search";
    if(sampling!="point_set_tolerance"&&sampling!="legacy_search"&&sampling!="exhaustive_integer_offsets")throw std::runtime_error("Unknown template sampling mode.");
    point_input_=source_is_points_||(template_&&sampling=="point_set_tolerance");
    point_regular_output_=false;
    if(template_&&source_is_points_&&sampling!="point_set_tolerance")throw std::runtime_error("Point sets require point_set_tolerance sampling.");
    int dimensions=3;
    for(int k=0;k<3;++k){source_dims_[k]=1;source_origin_[k]=0;source_spacing_[k]=1;output_dims_[k]=1;output_origin_[k]=file_origin_[k]=0;output_spacing_[k]=1;}
    if(cartesian_){source_dims_[0]=cartesian_->nx();source_dims_[1]=cartesian_->ny();source_dims_[2]=cartesian_->nz();for(int k=0;k<3;++k){source_origin_[k]=cartesian_->origin()[k];source_spacing_[k]=cartesian_->cell_dimensions()[k];}dimensions=source_dims_[2]==1?2:3;}
    if(point_input_)dimensions=integer(get("coordinate_dimensions","3"),2,3);
    std::ostringstream c;c<<std::setprecision(17)<<"[general]\nsource_type="<<(point_input_?"point_set":"cartesian")<<"\nsource_grid=source.gslib\nsource_property=value\ncoordinate_dimensions=3\ncalculate_target="<<target<<"\ncalculate_algorithm="<<algorithm<<"\ntarget_property=result\n";
    // All exported files contain XYZ plus value. Point-set dimensionality is
    // a separate geometry setting; three coordinate columns are always read.
    if(point_input_&&dimensions==2){c.str("");c.clear();c<<"[general]\nsource_type=point_set\nsource_grid=source.gslib\nsource_property=value\ncoordinate_dimensions=2\ncalculate_target="<<target<<"\ncalculate_algorithm="<<algorithm<<"\ntarget_property=result\n";}
    if(!point_input_)c<<"source_origin="<<join3(source_origin_)<<"\nsource_spacing="<<join3(source_spacing_)<<"\n";
    if(template_){
        const std::vector<double> nod=numbers(get("NOD","1,1,1")),nlag=numbers(get("NLAG","2,20,20")),lag=numbers(get("LAG","3,3,3"));
        if(nod.empty()||nod.size()>8||nod.size()!=nlag.size()||nod.size()!=lag.size())throw std::runtime_error("NOD, NLAG and LAG must have matching lengths (1 to 8 directions).");
        int order=1;for(size_t j=0;j<nod.size();++j){order+=integer(number(nod[j]),1,8);integer(number(nlag[j]),1);if(lag[j]==0)throw std::runtime_error("Lag separation cannot be zero.");}
        std::vector<int> display;const char* fields[]={"fst_dir","snd_dir","trd_dir"};const char* defaults[]={"3","2","0"};
        for(int k=0;k<3;++k){int d=integer(get(fields[k],defaults[k]),0,(int)nod.size());if(d){if(std::find(display.begin(),display.end(),d)!=display.end())throw std::runtime_error("Displayed template directions must be distinct.");display.push_back(d);}}
        if(display.empty())throw std::runtime_error("Select at least one display direction.");
        std::string fixed_text=parameters_["cross_pos"];
        if(fixed_text.empty()&&nod.size()>display.size()){for(size_t j=display.size();j<nod.size();++j){if(!fixed_text.empty())fixed_text+=",";fixed_text+="1";}}
        std::vector<double> fixed=numbers(fixed_text);if(fixed.size()!=nod.size()-display.size())throw std::runtime_error("Provide one fixed position for each hidden direction; clear Position if every direction is displayed.");
        const char* originfields[]={"xmn_","ymn_","zmn_"};for(int k=0;k<3;++k){output_origin_[k]=value(get(originfields[k],"0"));if(!_finite(output_origin_[k]))throw std::runtime_error("Output origin must be finite.");if(k<(int)display.size()){output_dims_[k]=(int)nlag[display[k]-1];output_spacing_[k]=fabs(lag[display[k]-1]);if(lag[display[k]-1]<0)output_origin_[k]+=(output_dims_[k]-1)*lag[display[k]-1];}file_origin_[k]=point_input_?output_origin_[k]:0;}
        c<<"target_origin="<<get("xmn_","0")<<','<<get("ymn_","0")<<','<<get("zmn_","0")<<"\n\n[template_shape]\ntemplate_sampling="<<sampling<<"\n";
        if(!point_input_)c<<"evaluate_raw_complete_tuples=true\n";
        const char* names[]={"NOD","NLAG","LAG","LTOL","AZM","ATOL","BANDWH","DIP","DTOL","BANDWD"};const char* defaults2[]={"1,1,1","2,20,20","3,3,3","0,0,0","180,90,0","0,0,0","0,0,0","0,0,0","0,0,0","0,0,0"};
        for(int k=0;k<10;++k){std::string t=get(names[k],defaults2[k]);if(numbers(t).size()!=nod.size())throw std::runtime_error(std::string(names[k])+" must contain one value per direction.");c<<names[k]<<'='<<t<<'\n';}
        if(sampling=="exhaustive_integer_offsets")c<<"integer_direction_vectors="<<get("integer_direction_vectors")<<'\n';
        c<<"minimum_template_count="<<get("minimum_template_count","1")<<"\nmax_template_tuples="<<get("max_template_tuples","0")<<"\nmax_candidate_visits="<<get("max_candidate_visits","0")<<"\n\n[target_cross_sections]\ndisplayed_dirs=";
        for(size_t j=0;j<display.size();++j){if(j)c<<',';c<<display[j];}c<<"\nfixed_positions="<<fixed_text<<"\n";
    }else{
        int order=integer(get("moment_order","4"),1);
        if(!get("center_cumulant_input_for_order_ge_2").empty()&&!flag(get("center_cumulant_input_for_order_ge_2")))throw std::runtime_error("Cumulants of order >=2 use within-window centering.");
        if(flag(get("output_standardized_cumulant","false"))){double eps=value(get("variance_epsilon","1e-12"));if(target!="Cumulant"||order<2||!_finite(eps)||eps<0)throw std::runtime_error("Standardization requires cumulant order >=2 and a finite nonnegative variance threshold.");}
        c<<"\n[sliding_window]\nmoment_cumulant_order="<<order<<"\ncenter_cumulant_input_for_order_ge_2=true\n";
        if(point_input_){
            const std::string mode=get("center_mode","sample_points");
            const std::string shape=get("neighborhood_shape","ellipsoid");
            const std::string bounds=get("grid_bounds_mode","auto");
            hosc_window_geometry::Shape geometry(shape,numbers(get("neighborhood_radii","1,1,1")),
                numbers(get("neighborhood_angles","90,0,0")),dimensions);
            without_sgems::PointSetData locations;locations.dimensions=dimensions;
            if(mode=="regular_grid"){
                locations.nodes.reserve(source_->size());
                for(int i=0;i<source_->size();++i){
                    Geostat_grid::location_type p=source_->location(i);without_sgems::PointSample node;
                    node.x=p[0];node.y=p[1];node.z=p[2];locations.nodes.push_back(node);
                }
            }
            hosc_window_geometry::CenterGrid centers;
            const std::vector<double> empty;
            centers.configure(locations,geometry,mode,bounds,
                mode=="regular_grid"&&bounds=="manual"?numbers(get("grid_bounds_min")):empty,
                mode=="regular_grid"&&bounds=="manual"?numbers(get("grid_bounds_max")):empty,
                mode=="regular_grid"?numbers(get("grid_start")):empty,
                mode=="regular_grid"?numbers(get("grid_spacing")):empty);
            point_regular_output_=centers.regular;
            if(point_regular_output_){
                if(!flag(get("includes_center","true")))throw std::runtime_error("Regular centers have no center record to exclude; set Include center record=true.");
                for(int k=0;k<3;++k){
                    output_dims_[k]=centers.dims[k];output_origin_[k]=file_origin_[k]=centers.origin[k];
                    output_spacing_[k]=centers.spacing[k]>0?centers.spacing[k]:1;
                }
            }
            c<<"neighborhood_radii="<<get("neighborhood_radii","1,1,1")<<"\nneighborhood_angles="<<get("neighborhood_angles","90,0,0")
             <<"\nneighborhood_shape="<<shape<<"\ncenter_mode="<<mode<<"\nmin_samples="<<integer(get("min_samples","2"),1)
             <<"\nincludes_center="<<get("includes_center","true")<<'\n';
            if(point_regular_output_){
                c<<"grid_bounds_mode="<<bounds<<"\ngrid_spacing="<<get("grid_spacing")<<'\n';
                if(bounds=="manual")c<<"grid_bounds_min="<<get("grid_bounds_min")<<"\ngrid_bounds_max="<<get("grid_bounds_max")<<'\n';
                if(!get("grid_start").empty())c<<"grid_start="<<get("grid_start")<<'\n';
            }
        }
        else{const char* wf[]={"windowl","windoww","windowh"};const char* sf[]={"stepx","stepy","stepz"};
            for(int k=0;k<3;++k){window_[k]=integer(get(wf[k],k==2&&dimensions==2?"1":"3"),1);step_[k]=integer(get(sf[k],k==2&&dimensions==2?"0":"1"),k==2&&dimensions==2?0:1);if(k==2&&dimensions==2){if(window_[k]!=1||step_[k]!=0)throw std::runtime_error("For a 2D grid set window Z=1 and step Z=0.");}if(window_[k]>source_dims_[k])throw std::runtime_error("Window exceeds source grid.");output_dims_[k]=step_[k]?1+(source_dims_[k]-window_[k])/step_[k]:1;output_origin_[k]=source_origin_[k]+.5*(window_[k]-1)*source_spacing_[k];output_spacing_[k]=source_spacing_[k]*(step_[k]?step_[k]:1);file_origin_[k]=(window_[k]-1)/2;}
            c<<"window_size="<<window_[0]<<','<<window_[1]<<','<<window_[2]<<"\nmove_step="<<step_[0]<<','<<step_[1]<<','<<step_[2]<<"\nrequire_full_window=true\nrequire_complete_window=true\nrequire_anchor_present=true\nskip_missing_samples=false\noutput_at_window_center=true\n";
        }
        c<<"invert_window_output="<<get("invert_window_output","false")<<"\nenable_percentile_clipping="<<get("enable_percentile_clipping","false")<<"\npercentile_vmin="<<get("percentile_vmin","2")<<"\npercentile_vmax="<<get("percentile_vmax","98")<<'\n';
    }
    c<<"\n[weighted_moment]\nenabled="<<(method=="uniform"?"false":"true")<<"\nweight_mode="<<method<<"\nnormalize_weights=true\nzero_first_order_moment=false\n";
    if(template_)c<<"weighting_switch="<<(flag(get("enable_weighted_moment","true"))?"true":"false")<<"\nselected_weight_mode="<<get("weight_mode","uniform")<<'\n';
    if(method=="external") {
        if(get("weight_property").empty()||!source_->property(get("weight_property")))throw std::runtime_error("Select an existing Weight property on the source object.");
        c<<"weight_property=observation_weight\nsource_weight_property="<<get("weight_property")<<"\nsource_weight_object="<<source_name_<<'\n';
    }
    if(method!="uniform"&&method!="external"){
        c<<"\n[declustering]\n";if(!get("declustering_dimensions").empty()&&get("declustering_dimensions")!="auto")c<<"dimensions="<<get("declustering_dimensions")<<'\n';
        if(method=="cell_declustering"){if(get("cell_size").empty())throw std::runtime_error("Cell declustering requires cell sizes in physical coordinate units.");c<<"cell_size="<<get("cell_size")<<"\nshifts_per_axis="<<get("shifts_per_axis","3")<<'\n';if(!get("declustering_origin").empty())c<<"origin="<<get("declustering_origin")<<'\n';if(!get("origin_offsets_file").empty())c<<"origin_offsets_file="<<get("origin_offsets_file")<<'\n';}
        else {if(get("domain_min").empty()||get("domain_max").empty())throw std::runtime_error("Polygonal and global kriging methods require explicit physical study-domain bounds.");c<<"domain_min="<<get("domain_min")<<"\ndomain_max="<<get("domain_max")<<'\n';if(method=="kriging_declustering"){if(get("kriging_ranges").empty()||get("kriging_integration_cells").empty())throw std::runtime_error("Global kriging requires covariance ranges and integration counts.");c<<"kriging_model="<<get("kriging_model","exponential")<<"\nkriging_ranges="<<get("kriging_ranges")<<"\nkriging_nugget_fraction="<<get("kriging_nugget_fraction","0.1")<<"\nkriging_integration_cells="<<get("kriging_integration_cells")<<'\n';}}
    }
    double output_count=1;for(int k=0;k<3;++k){output_count*=output_dims_[k];if(!_finite(output_origin_[k])||!_finite(output_spacing_[k])||output_spacing_[k]<=0)throw std::runtime_error("Invalid output geometry.");}if(output_count>INT_MAX)throw std::runtime_error("Output grid exceeds SGeMS index capacity.");
    c<<"\n[runtime]\noutput_dir=output\nwrite_template_tuples=false\n";ini_=c.str();
}

void EnginePlugin::export_source(){
    const bool external=effective_weight_mode()=="external";
    GsTLGridProperty* weights=external?source_->property(get("weight_property")):0;
    if(external&&!weights)throw std::runtime_error("External Weight property is unavailable.");
    std::ofstream out((run_dir_+"\\source.gslib").c_str());
    out<<"SGeMS selected property\n"<<(external?5:4)<<"\nx\ny\nz\nvalue\n";
    if(external)out<<"observation_weight\n";
    out<<std::setprecision(17);
    for(int i=0;i<source_->size();++i){if(point_input_){Geostat_grid::location_type p=source_->location(i);for(int k=0;k<3;++k){if(!_finite(p[k]))throw std::runtime_error("Source coordinates must be finite.");out<<p[k]<<' ';}}else{out<<i%source_dims_[0]<<' '<<(i/source_dims_[0])%source_dims_[1]<<' '<<i/(source_dims_[0]*source_dims_[1])<<' ';}if(property_->is_informed(i)&&_finite(property_->get_value(i)))out<<property_->get_value(i);else out<<"nan";
        if(external){
            if(property_->is_informed(i)&&_finite(property_->get_value(i))){
                if(!weights->is_informed(i)||!_finite(weights->get_value(i))||weights->get_value(i)<=0)
                    throw std::runtime_error("External weight must be informed, finite and strictly positive at source node "+number(i)+".");
                out<<' '<<weights->get_value(i);
            }else out<<" nan";
        }
        out<<'\n';}
    if(!out)throw std::runtime_error("Could not export the selected source property.");
    std::ofstream config((run_dir_+"\\config.ini").c_str());config<<ini_;if(!config)throw std::runtime_error("Could not write the resolved calculation parameters.");
}
void EnginePlugin::run_engine(const std::string& config,const std::string& logname){
    std::ofstream log((run_dir_+"\\"+logname).c_str());
    try{hosc_native::run_config(template_,run_dir_,config);log<<"Native in-process calculation completed.\n";}
    catch(const std::exception& e){log<<e.what()<<'\n';throw;}
}

void EnginePlugin::run_standardized(){
    if(template_||!flag(get("output_standardized_cumulant","false")))return;
    // Re-evaluate the unmodified statistic and variance on exactly the same
    // neighborhoods. Display inversion/clipping must not enter this ratio.
    for(int pass=0;pass<2;++pass){std::string raw=ini_;std::string key="moment_cumulant_order="+get("moment_order","4");if(pass)raw.replace(raw.find(key),key.size(),"moment_cumulant_order=2");
        const char* settings[]={"invert_window_output=","enable_percentile_clipping="};for(int k=0;k<2;++k){size_t start=raw.find(settings[k]);size_t end=raw.find('\n',start);raw.replace(start,end-start,std::string(settings[k])+"false");}
        std::string folder=pass?"variance":"raw";size_t at=raw.find("output_dir=output");raw.replace(at,17,"output_dir="+folder);std::string filename=folder+".ini";std::ofstream f((run_dir_+"\\"+filename).c_str());f<<raw;f.close();if(!f)throw std::runtime_error("Cannot write standardization configuration.");run_engine(filename,folder+".log");
    }
}

void EnginePlugin::import_result(GsTL_project* project){
    const std::string prefix=run_dir_+"\\output\\result";std::vector<Row> rows=read_rows(prefix+".txt");
    const bool point_output=!template_&&point_input_&&!point_regular_output_;
    if(point_output&&rows.size()!=(size_t)source_->size())throw std::runtime_error("Point-window output must retain every source record.");
    std::vector<int> nodes;nodes.reserve(rows.size());
    for(size_t i=0;i<rows.size();++i){int node=(int)i;if(!point_output){int xyz[3];for(int k=0;k<3;++k){double spacing=(template_||point_regular_output_)?output_spacing_[k]:(step_[k]?step_[k]:1);double d=(rows[i].xyz[k]-file_origin_[k])/spacing;xyz[k]=(int)floor(d+.5);if(fabs(d-xyz[k])>1e-6||xyz[k]<0||xyz[k]>=output_dims_[k])throw std::runtime_error("Engine result does not match the requested output geometry.");}node=(xyz[2]*output_dims_[1]+xyz[1])*output_dims_[0]+xyz[0];}nodes.push_back(node);}
    std::vector<unsigned char> seen(point_output?rows.size():output_dims_[0]*output_dims_[1]*output_dims_[2],0);
    for(size_t i=0;i<rows.size();++i){if(seen[nodes[i]])throw std::runtime_error("Duplicate output node.");seen[nodes[i]]=1;if(_finite(rows[i].value)&&fabs(rows[i].value)>FLT_MAX)throw std::runtime_error("Result exceeds the SGeMS float range.");}
    std::vector<Row> standardized;
    if(!template_&&flag(get("output_standardized_cumulant","false"))){std::vector<Row> raw=read_rows(run_dir_+"\\raw\\result.txt"),variance=read_rows(run_dir_+"\\variance\\result.txt");if(raw.size()!=rows.size()||variance.size()!=rows.size())throw std::runtime_error("Standardization support mismatch.");double eps=value(get("variance_epsilon","1e-12"));double exponent=.5*integer(get("moment_order","4"),2);standardized=raw;for(size_t i=0;i<raw.size();++i){for(int k=0;k<3;++k)if(raw[i].xyz[k]!=rows[i].xyz[k]||variance[i].xyz[k]!=rows[i].xyz[k])throw std::runtime_error("Standardization location mismatch.");standardized[i].value=_finite(raw[i].value)&&_finite(variance[i].value)&&variance[i].value>eps?raw[i].value/pow(variance[i].value,exponent):std::numeric_limits<double>::quiet_NaN();if(_finite(standardized[i].value)&&fabs(standardized[i].value)>FLT_MAX)throw std::runtime_error("Standardized result exceeds SGeMS float range.");}}
    std::string name;Geostat_grid* output=create_grid(project,this->name()+"_"+property_name_,rows,output_dims_,output_origin_,output_spacing_,point_output,name);
    if(!output)throw std::runtime_error("SGeMS output registration failed.");
    GsTLGridProperty* prop=output->add_property(property_name_);if(!prop)throw std::runtime_error("Cannot add output property.");
    GsTLGridProperty* valid=output->add_property(property_name_+"_valid");for(int i=0;i<output->size();++i)valid->set_value(0,i);
    for(size_t i=0;i<rows.size();++i){store(prop,nodes[i],rows[i].value);valid->set_value(_finite(rows[i].value)?1:0,nodes[i]);}
    if(!standardized.empty()){GsTLGridProperty* p=output->add_property(property_name_+"_standardized");for(size_t i=0;i<standardized.size();++i)store(p,nodes[i],standardized[i].value);}
    output->select_property(property_name_);
    if(flag(get("output_diagnostics","true"))){
        std::string path=prefix+(template_?".support.csv":".window_support.csv");
        if(exists(path)){std::ifstream in(path.c_str());std::string line;std::getline(in,line);std::vector<std::string> header=csv(line);int count=-1,ess=-1;for(size_t j=0;j<header.size();++j){if(header[j]=="support_count"||header[j]=="sample_count")count=(int)j;if(header[j]=="ESS"||header[j]=="effective_sample_size")ess=(int)j;}
            GsTLGridProperty* cp=count<0?0:output->add_property(property_name_+"_support");GsTLGridProperty* ep=ess<0?0:output->add_property(property_name_+"_ESS");size_t row=0;
            while(std::getline(in,line)){if(line.empty())continue;std::vector<std::string> fields=csv(line);int node=(int)row;if(template_||point_output){if(row>=nodes.size())throw std::runtime_error("Diagnostic row mismatch.");node=nodes[row];}if(node>=output->size())throw std::runtime_error("Diagnostic index exceeds output size.");if(cp&&count<(int)fields.size())store(cp,node,value(fields[count]));if(ep&&ess<(int)fields.size())store(ep,node,value(fields[ess]));++row;}
        }
    }
    if(effective_weight_mode()!="uniform"&&flag(get("output_point_weights","true"))&&exists(prefix+".point_weights.csv")){
        std::ifstream in((prefix+".point_weights.csv").c_str());std::string line;std::getline(in,line);std::vector<Row> weights;while(std::getline(in,line)){std::vector<std::string> f=csv(line);if(f.size()!=7&&f.size()!=8)throw std::runtime_error("Invalid point-weight audit row.");Row r;for(int k=0;k<3;++k)r.xyz[k]=value(f[k+1]);r.value=value(f[6]);weights.push_back(r);}std::string weightname;Geostat_grid* wg=create_grid(project,name+"_point_weights",weights,output_dims_,output_origin_,output_spacing_,true,weightname);const std::string weight_property=effective_weight_mode()=="external"?"observation_weight":"declustering_weight";GsTLGridProperty* wp=wg->add_property(weight_property);for(size_t i=0;i<weights.size();++i)store(wp,(int)i,weights[i].value);wg->select_property(weight_property);
    }
    if(template_&&flag(get("Output_3D_Slice_Enable","false")))import_slice(project,output,property_name_);
    GsTLlog<<"HOSC output: "<<name<<"; double results and resolved parameters: "<<run_dir_<<gstlIO::end;
    std::ofstream record((run_dir_+"\\sgems_result.txt").c_str());record<<"output_grid="<<name<<"\noutput_property="<<property_name_<<"\noutput_nodes="<<output->size()<<'\n';
    project->update();
}
void EnginePlugin::import_slice(GsTL_project* project,Geostat_grid* result,const std::string& resultprop){
    Geostat_grid* grid=result;std::string pname=resultprop;if(!get("Output_3D_Slice_Source_Grid").empty()){grid=dynamic_cast<Geostat_grid*>(Root::instance()->interface(gridModels_manager+"/"+get("Output_3D_Slice_Source_Grid")).raw_ptr());pname=get("Output_3D_Slice_Source_Property");}
    Cartesian_grid* cg=dynamic_cast<Cartesian_grid*>(grid);if(!cg||!cg->property(pname))throw std::runtime_error("Slice output requires a Cartesian source and property.");
    std::string plane=get("Output_3D_Slice_Plane","xy");int a=0,b=1,k=2;if(plane=="xz"){b=2;k=1;}else if(plane=="yz"){a=1;b=2;k=0;}else if(plane!="xy")throw std::runtime_error("Slice plane must be xy, xz or yz.");
    int dim[3]={cg->nx(),cg->ny(),cg->nz()};int pos=integer(get("Output_3D_Slice_Position","1"),1,dim[k])-1;std::vector<Row> rows;
    for(int j=0;j<dim[b];++j)for(int i=0;i<dim[a];++i){int ids[3];ids[a]=i;ids[b]=j;ids[k]=pos;int n=(ids[2]*dim[1]+ids[1])*dim[0]+ids[0];Row row;for(int q=0;q<3;++q)row.xyz[q]=cg->location(n)[q];row.value=cg->property(pname)->is_informed(n)?cg->property(pname)->get_value(n):std::numeric_limits<double>::quiet_NaN();rows.push_back(row);}
    int od[3]={dim[0],dim[1],dim[2]};od[k]=1;double origin[3],spacing[3];for(int q=0;q<3;++q){origin[q]=cg->origin()[q];spacing[q]=cg->cell_dimensions()[q];}origin[k]+=pos*spacing[k];std::string n;Geostat_grid* out=create_grid(project,"HOSC_slice_"+pname,rows,od,origin,spacing,false,n);GsTLGridProperty* p=out->add_property(pname);for(size_t i=0;i<rows.size();++i)store(p,(int)i,rows[i].value);out->select_property(pname);
}
int EnginePlugin::execute(GsTL_project* project){
    Error_messages_handler_xml runtime_errors;
    hosc_validation::FeedbackGuard feedback(&runtime_errors,template_?"hosc_template_validation_summary":"hosc_window_validation_summary");
    RunPerformance performance;
    std::map<std::string,std::string> metrics;
    metrics["module"]=name();metrics["algorithm"]=get("Algorithm_Button","AVL Based algorithm");
    metrics["target"]=get("Target_Button","Cumulant");metrics["case_name"]=property_name_;
    metrics["source_grid"]=source_name_;metrics["source_property"]=get("Source_Property");
    metrics["weight_mode"]=effective_weight_mode();metrics["moment_order"]=get("moment_order");
    metrics["NOD"]=get("NOD");metrics["NLAG"]=get("NLAG");metrics["resolved_config"]="config.ini";
    metrics["plugin_build"]=__DATE__ " " __TIME__;
    run_dir_.clear();
    try{
        if(!project||!source_||!property_)throw std::runtime_error("Plugin is not initialized.");
        run_dir_=new_run_directory(get("audit_output_directory"));
        export_source();performance.checkpoint();
        run_engine();performance.checkpoint();
        run_standardized();performance.checkpoint();
        import_result(project);performance.stop();
        if(!performance.write(run_dir_,metrics,"success",""))
            GsTLcerr<<name()<<": calculation completed but performance records could not be saved: "<<run_dir_<<gstlIO::end;
        else GsTLlog<<"HOSC performance: "<<performance.elapsed_seconds()<<" s; sampled peak working set "<<performance.peak_mib()<<" MiB; records: "<<run_dir_<<"\\performance.csv"<<gstlIO::end;
        return 0;
    }
    catch(const std::exception& e){performance.stop();if(!run_dir_.empty())performance.write(run_dir_,metrics,"failed",e.what());runtime_errors.report("audit_output_directory",escaped(std::string(e.what())+"; records: "+run_dir_));GsTLcerr<<name()<<": "<<e.what()<<"; work directory: "<<run_dir_<<gstlIO::end;return 1;}
    catch(...){performance.stop();if(!run_dir_.empty())performance.write(run_dir_,metrics,"failed","Unexpected calculation failure");runtime_errors.report("audit_output_directory","Unexpected calculation failure; inspect the run records.");GsTLcerr<<name()<<": unexpected failure; work directory: "<<run_dir_<<gstlIO::end;return 1;}
}
}
