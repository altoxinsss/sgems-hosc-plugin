#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#endif

namespace fs=std::filesystem;
namespace {
using Row=std::map<std::string,std::string>;
struct Sample { std::vector<double> values; double weight; };
struct Grid { int nx,ny,nz; std::map<int,double> informed; };
void require(bool ok,const std::string& message) { if(!ok) throw std::runtime_error(message); }
bool close_value(double a,double b) {
    return std::isfinite(a)&&std::isfinite(b)&&std::abs(a-b)<=1e-9+5e-9*std::max(std::abs(a),std::abs(b));
}
int run(const fs::path& exe,const fs::path& config) {
#ifdef _WIN32
    std::wstring command=L"\""+exe.wstring()+L"\" \""+config.wstring()+L"\"";
    std::vector<wchar_t> buffer(command.begin(),command.end());buffer.push_back(0);
    STARTUPINFOW startup{};startup.cb=sizeof(startup);PROCESS_INFORMATION process{};
    if(!CreateProcessW(exe.c_str(),buffer.data(),nullptr,nullptr,FALSE,CREATE_NO_WINDOW,nullptr,nullptr,&startup,&process))
        throw std::runtime_error("Cannot start Cartesian template process.");
    const DWORD waited=WaitForSingleObject(process.hProcess,60000);
    if(waited!=WAIT_OBJECT_0) {
        TerminateProcess(process.hProcess,99);CloseHandle(process.hThread);CloseHandle(process.hProcess);
        throw std::runtime_error("Cartesian declustering check timed out.");
    }
    DWORD code=1;GetExitCodeProcess(process.hProcess,&code);
    CloseHandle(process.hThread);CloseHandle(process.hProcess);return static_cast<int>(code);
#else
    return std::system(("\""+exe.string()+"\" \""+config.string()+"\"").c_str());
#endif
}
std::vector<std::string> split(const std::string& line,char delimiter=',') {
    std::vector<std::string> result;std::string value;bool quoted=false;
    for(char c:line) {
        if(c=='"') quoted=!quoted;
        else if(c==delimiter&&!quoted) {result.push_back(value);value.clear();}
        else if(c!='\r') value+=c;
    }
    result.push_back(value);return result;
}
std::vector<Row> table(const fs::path& path) {
    std::ifstream in(path);require(static_cast<bool>(in),"Missing CSV: "+path.string());
    std::string line;std::getline(in,line);const auto names=split(line);std::vector<Row> rows;
    while(std::getline(in,line)) {
        if(line.empty())continue;const auto cells=split(line);require(cells.size()==names.size(),"Malformed CSV: "+path.string());
        Row row;for(size_t i=0;i<names.size();++i)row[names[i]]=cells[i];rows.push_back(row);
    }
    return rows;
}
double number(const Row& row,const std::string& key) {
    const auto it=row.find(key);require(it!=row.end(),"Missing diagnostic column: "+key);return std::stod(it->second);
}
std::string contents(const fs::path& path) {std::ifstream in(path);return {std::istreambuf_iterator<char>(in),{}};}
std::vector<double> grid_values(const fs::path& path) {
    std::ifstream in(path);require(static_cast<bool>(in),"Missing result grid.");std::string line;
    for(int i=0;i<4;++i)std::getline(in,line);
    std::vector<double> result;
    while(std::getline(in,line)) {
        std::istringstream stream(line);std::string token,last;
        while(stream>>token)last=token;
        if(!last.empty())result.push_back(std::stod(last));
    }
    return result;
}
void write_grid(const fs::path& path,const Grid& grid,bool sparse=false,
                std::array<int,3> index_origin={{0,0,0}},const std::string& missing="nan") {
    std::ofstream out(path);out<<"Cartesian independent fixture\n"<<(sparse?grid.informed.size():grid.nx*grid.ny*grid.nz)
        <<" 4\n# INDEX X Y Z VAL\n# Missing nodes are not observations\n"<<std::setprecision(17);
    for(int z=0;z<grid.nz;++z)for(int y=0;y<grid.ny;++y)for(int x=0;x<grid.nx;++x) {
        const int id=1+x+grid.nx*(y+grid.ny*z);const auto it=grid.informed.find(id);
        if(sparse&&it==grid.informed.end())continue;
        out<<id<<' '<<x+index_origin[0]<<' '<<y+index_origin[1]<<' '<<z+index_origin[2]<<' ';
        if(it==grid.informed.end())out<<missing;else out<<it->second;out<<'\n';
    }
}
std::string weights(const std::string& params,const std::string& mode="cell_declustering") {
    return "\n[weighted_moment]\nenabled=true\nnormalize_weights=true\nweight_mode="+mode+"\n[declustering]\n"+params;
}
fs::path config(const fs::path& root,const std::string& name,const fs::path& source,
                const std::string& route,const std::string& algorithm,const std::string& weighting,
                bool dimension3=false,bool audit=true,const std::string& general="",
                const std::string& shape="",int order=4,const std::string& cross_section="") {
    const fs::path path=root/(name+".ini");std::ofstream out(path);
    out<<"[general]\nsource_type=cartesian\nsource_grid="<<source.string()<<"\nsource_property=value\ncoordinate_dimensions="
        <<(dimension3?3:2)<<"\ncalculate_target=Cumulant\ncalculate_algorithm="<<algorithm<<"\ntarget_property=result\n"<<general
        <<"\n[template_shape]\ntemplate_sampling="<<route<<"\nNOD="<<order-2<<",1\nNLAG=2,2\nLAG=1,1\nLTOL=0,0\nAZM=0,0\n"
        <<"ATOL=0,0\nBANDWH=0,0\nDIP=0,0\nDTOL=0,0\nBANDWD=0,0\ninteger_direction_vectors="
        <<(dimension3?"1,0,0;1,0,0":"1,0;1,0")<<'\n'<<shape
        <<"\n[target_cross_sections]\ndisplayed_dirs=1,2\nfixed_positions=\n"<<cross_section<<"\n[runtime]\noutput_dir="<<(root/name).string()
        <<"\nwrite_template_tuples="<<(audit?"true":"false")<<'\n'<<weighting;
    return path;
}
long double moment(const std::vector<Sample>& samples,const std::vector<int>& powers) {
    long double numerator=0,denominator=0;
    for(const auto& sample:samples) {
        long double product=1;
        for(size_t j=0;j<powers.size();++j)product*=std::pow(static_cast<long double>(sample.values[j]),powers[j]);
        numerator+=sample.weight*product;denominator+=sample.weight;
    }
    return numerator/denominator;
}
// Independent raw-moment partition identity. Fixed hand-enumerated tuples,
// with no production search, centering, moment cache, or recurrence calls.
double cumulant(const std::vector<Sample>& samples,const std::vector<int>& powers) {
    std::vector<int> occurrence;for(size_t j=0;j<powers.size();++j)for(int k=0;k<powers[j];++k)occurrence.push_back(static_cast<int>(j));
    std::vector<std::vector<int>> blocks;long double result=0;
    std::function<void(size_t)> visit=[&](size_t next) {
        if(next==occurrence.size()) {
            long double term=1;
            for(const auto& block:blocks) {
                std::vector<int> beta(powers.size(),0);for(int variable:block)++beta[variable];term*=moment(samples,beta);
            }
            for(size_t k=1;k<blocks.size();++k)term*=-static_cast<long double>(k);
            result+=term;return;
        }
        const size_t count=blocks.size();
        for(size_t b=0;b<count;++b) {blocks[b].push_back(occurrence[next]);visit(next+1);blocks[b].pop_back();}
        blocks.push_back({occurrence[next]});visit(next+1);blocks.pop_back();
    };visit(0);return static_cast<double>(result);
}
void check_weights(const fs::path& output,const Grid& grid,const std::map<int,double>& expected) {
    const auto rows=table(output/"result.point_weights.csv");size_t informed=0;double sum=0;
    for(const auto& row:rows) {
        const int id=static_cast<int>(number(row,"node_id"));const auto it=expected.find(id);const double wanted=it==expected.end()?0:it->second;
        require(close_value(number(row,"weight"),wanted),"Cell occupancy/physical coordinates gave wrong weight for dense node "+std::to_string(id));
        if(wanted>0)++informed;
        sum+=number(row,"weight");
    }
    require(informed==grid.informed.size()&&close_value(sum,1),"Missing or omitted grid cells affected occupancy.");
}
void check_audit(const fs::path& output) {
    const auto point_rows=table(output/"result.point_weights.csv");std::map<int,double> point_weights;
    for(const auto& row:point_rows)point_weights[static_cast<int>(number(row,"node_id"))]=number(row,"weight");
    const auto tuples=table(output/"result.tuples.csv");
    for(const auto& row:tuples) {
        std::set<int> ids;for(const auto& id:split(row.at("node_ids"),';'))ids.insert(std::stoi(id));
        long double product=1;for(int id:ids)product*=point_weights.at(id);
        require(close_value(number(row,"tuple_weight"),std::pow(static_cast<double>(product),1.0/ids.size())),"Coincident locations or NOD counted a node weight repeatedly.");
    }
    for(const auto& row:table(output/"result.support.csv")) {
        double sum=0,square=0;size_t count=0;
        for(const auto& tuple:tuples)if(tuple.at("lag_cell")==row.at("lag_cell")) {
            ++count;const double w=number(tuple,"tuple_weight");sum+=w;square+=w*w;
        }
        require(count==number(row,"support_count")&&close_value(sum,number(row,"weight_sum"))&&
            close_value(square>0?sum*sum/square:0,number(row,"ESS")),"Support, tuple-weight sum, or ESS differs from the tuple record.");
    }
}
void check_pair_result(const fs::path& output,const std::vector<Sample>& scalar,const std::vector<Sample>& pairs,int order=4) {
    const auto rows=table(output/"result.support.csv");require(rows.size()==4,"Wrong lag cross-section size.");
    require(number(rows[0],"support_count")==scalar.size()&&close_value(number(rows[0],"value"),cumulant(scalar,{order})),"Zero-lag weighted cumulant differs from independent raw moments.");
    require(number(rows[3],"support_count")==pairs.size()&&close_value(number(rows[3],"value"),cumulant(pairs,{1,order-1})),"Weighted repeated-lag cumulant differs from fixed hand-enumerated support.");
    std::vector<double> expected={cumulant(pairs,{order-1,1}),cumulant(pairs,{2,order-2})};
    require((close_value(number(rows[1],"value"),expected[0])&&close_value(number(rows[2],"value"),expected[1]))||
        (close_value(number(rows[1],"value"),expected[1])&&close_value(number(rows[2],"value"),expected[0])),"Lower-order dependency support/weights changed across the target.");
    require(rows[3].at("status")=="ok","Supported nonuniform result has incorrect status.");
}
}

int main(int argc,char** argv) {
    try {
        require(argc==2,"Usage: cartesian_declustering_checks <template_cli>");const fs::path exe=fs::absolute(argv[1]);
        const fs::path root=fs::current_path()/"test_output"/"cartesian_declustering_checks";fs::create_directories(root);
        // Four occupied cells contain three, two, one, and one observations. The
        // three disjoint +x pairs are explicitly listed, not inferred from output.
        const Grid planar={6,3,1,{{1,1},{2,2},{7,3},{8,6},{13,4},{17,7},{18,10}}};
        const std::map<int,double> planar_weights={{1,1.0/12},{2,.125},{7,1.0/12},{8,.125},{13,1.0/12},{17,.25},{18,.25}};
        const std::vector<Sample> scalar={{{1},1.0/12},{{2},.125},{{3},1.0/12},{{6},.125},{{4},1.0/12},{{7},.25},{{10},.25}};
        const std::vector<Sample> pairs={{{1,2},std::sqrt(1.0/96)},{{3,6},std::sqrt(1.0/96)},{{7,10},.25}};
        const auto weighted=weights("cell_size=1,3\norigin=0,0\nshifts_per_axis=1\n");
        const fs::path source=root/"clustered.grid";write_grid(source,planar);
        for(const std::string route:{"legacy_search","exhaustive_integer_offsets"}) {
            for(int order:{3,4,6})for(const std::string algorithm:{"AVL","Bell"}) {
                const std::string name=route+algorithm+std::to_string(order);
                require(run(exe,config(root,name,source,route,algorithm,weighted,false,true,"","",order))==0,"Cartesian nonuniform run failed: "+name);
                check_weights(root/name,planar,planar_weights);check_pair_result(root/name,scalar,pairs,order);check_audit(root/name);
                const auto meta=contents(root/name/"result.meta.txt");
                require(meta.find("template_sampling="+route)!=std::string::npos&&meta.find("cell_declustering")!=std::string::npos,
                    "Cartesian provenance lost the actual sampler or weight mode.");
            }
            const std::string no_audit=route+"_no_audit";
            require(run(exe,config(root,no_audit,source,route,"AVL",weighted,false,false))==0,"Weighted run with audit disabled failed.");
            check_pair_result(root/no_audit,scalar,pairs);
            require(!fs::exists(root/no_audit/"result.tuples.csv"),"Optional tuple audit was emitted despite false switch.");
            if(route=="legacy_search") {
                const std::string placeholder_name=route+"_fixed_placeholder";
                require(run(exe,config(root,placeholder_name,source,route,"AVL",weighted,false,true,"","",4,"fixed_positions=1\n"))==0,
                    "Legacy Cartesian cell declustering rejected the existing no-hidden-direction fixed_positions=1 placeholder.");
                check_pair_result(root/placeholder_name,scalar,pairs);
                require(contents(root/placeholder_name/"result.support.csv")==contents(root/no_audit/"result.support.csv"),
                    "Legacy fixed_positions=1 placeholder changed full displayed-direction output.");
            }
            for(const std::string algorithm:{"AVL","Bell"}) {
                const std::string raw_name=route+"_raw_"+algorithm;
                require(run(exe,config(root,raw_name,source,route,algorithm,weighted,false,true,"calculate_target=Moment\n"))==0,"Raw weighted Cartesian moment failed.");
                const auto raw_rows=table(root/raw_name/"result.support.csv");
                require(raw_rows.size()==4&&number(raw_rows[0],"support_count")==scalar.size()&&
                    close_value(number(raw_rows[0],"value"),static_cast<double>(moment(scalar,{4})))&&
                    number(raw_rows[3],"support_count")==pairs.size()&&
                    close_value(number(raw_rows[3],"value"),static_cast<double>(moment(pairs,{1,3}))),
                    "Raw Cartesian moment was centered or did not use normalized complete-tuple weights.");
                check_audit(root/raw_name);
            }
            const std::string minimum_name=route+"_minimum";
            require(run(exe,config(root,minimum_name,source,route,"AVL",weighted,false,true,"","minimum_template_count=4\n"))==0,"Minimum complete-template support failed.");
            const auto minimum_rows=table(root/minimum_name/"result.support.csv");
            require(minimum_rows[0].at("status")=="ok"&&minimum_rows[3].at("status")=="insufficient_support"&&
                number(minimum_rows[3],"support_count")==3&&std::isnan(number(minimum_rows[3],"value")),
                "Support threshold changed the sample collection or retained an unsupported estimate.");
            check_audit(root/minimum_name);
            const std::string empty_name=route+"_empty_lag";
            require(run(exe,config(root,empty_name,source,route,"Bell",weighted,false,true,"","LAG=50,50\n"))==0,"Unsupported lag prevented valid zero-lag output.");
            const auto empty_rows=table(root/empty_name/"result.support.csv");
            require(empty_rows[0].at("status")=="ok"&&empty_rows[3].at("status")=="insufficient_support"&&
                number(empty_rows[3],"support_count")==0&&std::isnan(number(empty_rows[3],"value"))&&
                number(empty_rows[3],"weight_sum")==0&&number(empty_rows[3],"ESS")==0,
                "Empty lag did not produce explicit missing output and zero diagnostics.");

            // Omitting unobserved nodes has exactly the same meaning as NaN.
            const fs::path sparse=root/"sparse.grid";write_grid(sparse,planar,true);
            const std::string sparse_name=route+"_sparse";
            require(run(exe,config(root,sparse_name,sparse,route,"AVL",weighted))==0,"Sparse Cartesian declustering failed.");
            check_weights(root/sparse_name,planar,planar_weights);check_pair_result(root/sparse_name,scalar,pairs);
            const fs::path sentinel=root/"sentinel.grid";write_grid(sentinel,planar,false,{{0,0,0}},"-999");
            const std::string sentinel_name=route+"_sentinel";
            require(run(exe,config(root,sentinel_name,sentinel,route,"Bell",weighted,false,true,"missing_value=-999\n"))==0,"Configured missing sentinel failed.");
            check_pair_result(root/sentinel_name,scalar,pairs);

            // Physical coordinates are explicitly different from input indices.
            // Ignoring the transform groups different observed columns together.
            const fs::path translated=root/"translated.grid";write_grid(translated,planar,false,{{10,20,30}});
            const std::string physical_name=route+"_physical";
            require(run(exe,config(root,physical_name,translated,route,"AVL",
                weights("cell_size=2,9\norigin=100,-30\nshifts_per_axis=1\n"),false,true,"",
                "source_origin=100,-30,7\nsource_spacing=2,3,4\n"))==0,"Physical Cartesian declustering failed.");
            check_weights(root/physical_name,planar,planar_weights);check_pair_result(root/physical_name,scalar,pairs);

            const Grid volume={6,3,2,{{1,1},{2,2},{7,3},{8,6},{13,4},{35,7},{36,10}}};
            const fs::path source3=root/"volume.grid";write_grid(source3,volume);
            const std::string volume_name=route+"_3d";
            require(run(exe,config(root,volume_name,source3,route,"Bell",weights("cell_size=1,3,1\norigin=0,0,0\nshifts_per_axis=1\n"),true))==0,"3D Cartesian declustering failed.");
            check_weights(root/volume_name,volume,{{1,1.0/12},{2,.125},{7,1.0/12},{8,.125},{13,1.0/12},{35,.25},{36,.25}});
            check_pair_result(root/volume_name,scalar,pairs);check_audit(root/volume_name);

            // Complete Cartesian observations and equally occupied cells must
            // leave the original sampler's uniform statistic unchanged.
            // Avoid the historical sampler's centered-value -1 sentinel in
            // this compatibility fixture; informed-value correctness above
            // is independently checked against all finite observations.
            const Grid complete={4,2,1,{{1,1},{2,2},{3,4},{4,8},{5,3},{6,5},{7,7},{8,11}}};
            const fs::path complete_source=root/"complete.grid";write_grid(complete_source,complete);
            const std::string equal_name=route+"_equal",uniform_name=route+"_uniform";
            require(run(exe,config(root,equal_name,complete_source,route,"AVL",weights("cell_size=2,2\norigin=0,0\nshifts_per_axis=1\n")))==0,"Equal-cell Cartesian run failed.");
            require(run(exe,config(root,uniform_name,complete_source,route,"AVL",weights("","uniform")))==0,"Original uniform Cartesian run failed.");
            const auto equal_values=grid_values(root/equal_name/"result.txt"),uniform_values=grid_values(root/uniform_name/"result.txt");
            require(equal_values.size()==uniform_values.size()&&equal_values.size()==4,"Equal/uniform output layout differs.");
            for(size_t i=0;i<equal_values.size();++i)require(close_value(equal_values[i],uniform_values[i]),"Equal-cell weighting changed original Cartesian sampler results.");
            check_audit(root/equal_name);
            for(const std::string params:{"cell_size=0,2\n","cell_size=2\n","cell_size=2,2\nshifts_per_axis=0\n","cell_size=2,2\norigin=nan,0\n"}) {
                require(run(exe,config(root,route+"_invalid",source,route,"AVL",weights(params)))!=0,"Invalid Cartesian cell parameters were silently accepted.");
            }
            for(const std::string transform:{"source_spacing=0,1,1\n","source_spacing=nan,1,1\n","source_origin=1,2\n"})
                require(run(exe,config(root,route+"_invalid_transform",source,route,"AVL",weighted,false,true,"",transform))!=0,"Invalid physical grid transform was silently accepted.");
            const fs::path absent=root/"all_missing.grid";write_grid(absent,{2,2,1,{}});
            require(run(exe,config(root,route+"_all_missing",absent,route,"AVL",weighted))!=0,"All-missing Cartesian property generated a weighted statistic.");
        }
        // Direction 3 has lag two and radial tolerance one: the historical
        // rectangular search first selects x+1. Directions 1 and 2 select
        // x+1 exactly. Moving direction 3 into the hidden first slot must move
        // its LTOL with it. The all-one target has the same three explicitly
        // known pairs, hence the independent fourth-cumulant oracle applies.
        const std::string unequal_tolerances=
            "NOD=1,1,1\nNLAG=2,2,2\nLAG=1,1,2\nLTOL=0,0,1\nAZM=0,0,0\nATOL=0,0,0\n"
            "BANDWH=0,0,0\nDIP=0,0,0\nDTOL=0,0,0\nBANDWD=0,0,0\ninteger_direction_vectors=1,0;1,0;1,0\n";
        std::vector<std::string> fixed_target_instances;
        for(const std::string name:{"tolerance_identity","tolerance_reordered"}) {
            const std::string section=name=="tolerance_identity"?"displayed_dirs=2,3\nfixed_positions=2\n":
                "displayed_dirs=1,2\nfixed_positions=2\n";
            require(run(exe,config(root,name,source,"legacy_search","AVL",weighted,false,true,"",unequal_tolerances,4,section))==0,
                "Cartesian hidden/displayed tolerance configuration failed.");
            const auto rows=table(root/name/"result.support.csv");
            const auto target=std::find_if(rows.begin(),rows.end(),[](const Row& row){return row.at("lag_positions")=="1;1;1";});
            require(target!=rows.end()&&number(*target,"support_count")==3&&
                close_value(number(*target,"value"),cumulant(pairs,{1,3})),
                "LTOL stayed attached to a reordered array slot instead of its configured direction.");
            std::vector<std::string> instances;
            for(const auto& tuple:table(root/name/"result.tuples.csv"))if(tuple.at("lag_cell")==target->at("lag_cell"))
                instances.push_back(tuple.at("node_ids")+"|"+tuple.at("values")+"|"+tuple.at("tuple_weight"));
            std::sort(instances.begin(),instances.end());
            if(fixed_target_instances.empty())fixed_target_instances=instances;
            else require(instances==fixed_target_instances,"Hidden/displayed direction permutation changed the physical template instances or weights.");
            check_audit(root/name);
        }
        for(const std::string mode:{"polygonal_declustering","kriging_declustering"})
        for(const std::string route:{"legacy_search","exhaustive_integer_offsets"})
        for(const std::string algorithm:{"Bell","AVL"}) for(const int order:{4,5}) {
            const std::string name=mode+route+algorithm+std::to_string(order);
            std::string parameters="domain_min=-0.5,-0.5\ndomain_max=5.5,2.5\n";
            if(mode=="kriging_declustering") parameters+="kriging_model=exponential\nkriging_ranges=2,2\nkriging_nugget_fraction=0.2\nkriging_integration_cells=24,12\n";
            require(run(exe,config(root,name,source,route,algorithm,weights(parameters,mode),false,true,"","",order))==0,"New Cartesian mode failed: "+name);
            check_audit(root/name);
            const auto records=table(root/name/"result.tuples.csv");
            const auto rows=table(root/name/"result.support.csv");
            for(const auto& row:rows) {
                std::vector<Sample> observed;
                for(const auto& record:records) if(record.at("lag_cell")==row.at("lag_cell")) {
                    Sample sample;sample.weight=number(record,"tuple_weight");
                    for(const auto& value:split(record.at("values"),';'))sample.values.push_back(std::stod(value));
                    observed.push_back(sample);
                }
                require(!observed.empty()&&close_value(number(row,"value"),cumulant(observed,{1,order-2,1})),
                    "New Cartesian statistic differs from independent raw-partition oracle: "+name);
            }
            const auto meta=contents(root/name/"result.meta.txt");
            require(meta.find("weight_mode="+mode)!=std::string::npos&&meta.find("declustering_cell_size=")==std::string::npos,
                "Cartesian metadata reports the wrong method");
        }
        std::cout<<"Cartesian three-method declustering independent checks: PASS\n";return 0;
    }catch(const std::exception& error) {std::cerr<<"Cartesian declustering checks FAILED: "<<error.what()<<'\n';return 1;}
}
