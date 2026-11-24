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

namespace fs = std::filesystem;
namespace {
struct Point { double x, y, z, value; };
using Row = std::map<std::string, std::string>;
struct Sample { std::vector<double> values; double weight; };
void require(bool ok, const std::string& message) { if (!ok) throw std::runtime_error(message); }
bool close_value(double a, double b) {
    return std::isfinite(a) && std::isfinite(b) &&
        std::abs(a-b) <= 1e-10 + 5e-9 * std::max(std::abs(a), std::abs(b));
}
int run(const fs::path& exe, const fs::path& config) {
#ifdef _WIN32
    std::wstring command = L"\"" + exe.wstring() + L"\" \"" + config.wstring() + L"\"";
    std::vector<wchar_t> buffer(command.begin(), command.end()); buffer.push_back(0);
    STARTUPINFOW startup{}; startup.cb = sizeof(startup); PROCESS_INFORMATION process{};
    if (!CreateProcessW(exe.c_str(), buffer.data(), nullptr, nullptr, FALSE,
                        CREATE_NO_WINDOW, nullptr, nullptr, &startup, &process))
        throw std::runtime_error("Cannot start template_cli.");
    const DWORD waited = WaitForSingleObject(process.hProcess, 60000);
    if (waited != WAIT_OBJECT_0) {
        TerminateProcess(process.hProcess, 99); CloseHandle(process.hThread); CloseHandle(process.hProcess);
        throw std::runtime_error("template_cli declustering test timed out.");
    }
    DWORD code=1; GetExitCodeProcess(process.hProcess, &code);
    CloseHandle(process.hThread); CloseHandle(process.hProcess); return static_cast<int>(code);
#else
    return std::system(("\"" + exe.string() + "\" \"" + config.string() + "\"").c_str());
#endif
}
void write_points(const fs::path& path, const std::vector<Point>& points) {
    std::ofstream out(path); out << "declustering fixture\n4\nx\ny\nz\nphi\n" << std::setprecision(17);
    for (const auto& p : points) out << p.x << ' ' << p.y << ' ' << p.z << ' ' << p.value << '\n';
}
std::vector<std::string> split(const std::string& line, char delimiter=',') {
    std::vector<std::string> result; std::string value; bool quoted=false;
    for (char c : line) {
        if (c=='"') quoted=!quoted;
        else if (c==delimiter && !quoted) { result.push_back(value); value.clear(); }
        else if (c!='\r') value += c;
    }
    result.push_back(value); return result;
}
std::vector<Row> table(const fs::path& path) {
    std::ifstream in(path); require(static_cast<bool>(in), "Missing CSV: " + path.string());
    std::string line; std::getline(in, line); const auto names=split(line); std::vector<Row> rows;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        const auto values=split(line); require(names.size()==values.size(), "Malformed CSV row.");
        Row row; for (size_t i=0; i<names.size(); ++i) row[names[i]]=values[i]; rows.push_back(row);
    }
    return rows;
}
double number(const Row& row, const std::string& name) {
    const auto found=row.find(name); require(found!=row.end(), "Missing CSV column: " + name);
    return std::stod(found->second);
}
std::string shape(int nod, int count=2) {
    return "NOD="+std::to_string(nod)+"\nNLAG="+std::to_string(count)+
        "\nLAG=1.5\nLTOL=0\nAZM=0\nATOL=0\nBANDWH=0\nDIP=0\nDTOL=0\nBANDWD=0\n";
}
std::string weighting(const std::string& parameters="cell_size=2,35,2\nshifts_per_axis=1\norigin=-11,-31,-1\n",
                      const std::string& mode="cell_declustering", const std::string& more="") {
    return "\n[weighted_moment]\nenabled=true\nnormalize_weights=true\nweight_mode="+mode+
        "\n"+more+"\n[declustering]\n"+parameters;
}
fs::path config(const fs::path& root, const std::string& name, const fs::path& points,
    const std::string& algorithm, const std::string& geometry, const std::string& weights,
    const std::string& general="", bool audit=true) {
    const fs::path path=root/(name+".ini"); std::ofstream out(path);
    out << "[general]\nsource_type=point_set\nsource_grid=" << points.string()
        << "\nsource_property=phi\ncoordinate_dimensions=3\ncalculate_target=Cumulant\ncalculate_algorithm="
        << algorithm << "\ntarget_property=result\n" << general
        << "\n[template_shape]\ntemplate_sampling=point_set_tolerance\n" << geometry
        << "\n[runtime]\noutput_dir=" << (root/name).string() << "\nwrite_template_tuples="
        << (audit?"true":"false") << '\n' << weights;
    return path;
}
long double moment(const std::vector<Sample>& samples, const std::vector<int>& powers) {
    long double numerator=0, denominator=0;
    for (const auto& s : samples) {
        long double product=1; for (size_t j=0;j<powers.size();++j) product *= std::pow(static_cast<long double>(s.values[j]), powers[j]);
        numerator += s.weight*product; denominator += s.weight;
    }
    return numerator/denominator;
}
// Independent raw-moment set-partition oracle. There is no production header
// dependency, no production cache, and no centering step in this reference.
double cumulant(const std::vector<Sample>& samples, const std::vector<int>& powers) {
    std::vector<int> occurrences;
    for (size_t j=0;j<powers.size();++j) for (int p=0;p<powers[j];++p) occurrences.push_back(static_cast<int>(j));
    std::vector<std::vector<int>> blocks; long double result=0;
    std::function<void(size_t)> visit = [&](size_t occurrence) {
        if (occurrence==occurrences.size()) {
            long double term=1;
            for (const auto& block : blocks) {
                std::vector<int> beta(powers.size(),0); for (int variable : block) ++beta[variable];
                term *= moment(samples,beta);
            }
            for (size_t factor=1;factor<blocks.size();++factor) term *= -static_cast<long double>(factor);
            result += term; return;
        }
        const size_t count=blocks.size();
        for (size_t b=0;b<count;++b) { blocks[b].push_back(occurrences[occurrence]); visit(occurrence+1); blocks[b].pop_back(); }
        blocks.push_back({occurrences[occurrence]}); visit(occurrence+1); blocks.pop_back();
    };
    visit(0); return static_cast<double>(result);
}
double scalar_cumulant(const std::vector<Sample>& samples, int order) {
    const long double m1=moment(samples,{1}), m2=moment(samples,{2});
    if(order==2) return static_cast<double>(m2-m1*m1);
    const long double m3=moment(samples,{3});
    if(order==3) return static_cast<double>(m3-3*m2*m1+2*m1*m1*m1);
    if(order==4) return static_cast<double>(moment(samples,{4})-4*m3*m1-3*m2*m2+12*m2*m1*m1-6*m1*m1*m1*m1);
    std::vector<long double> centered(order+1,0), k(order+1,0); centered[0]=1;
    long double denominator=0; for(const auto& s:samples) denominator+=s.weight;
    for(int r=1;r<=order;++r) {
        for(const auto& s:samples) centered[r]+=s.weight*std::pow(static_cast<long double>(s.values[0])-m1,r);
        centered[r]/=denominator; k[r]=centered[r]; long double choose=1;
        for(int j=1;j<r;++j) { choose*=static_cast<long double>(r-j)/j; k[r]-=choose*k[r-j]*centered[j]; }
    }
    return static_cast<double>(k[order]);
}
void check_point_weights(const fs::path& output, const std::vector<double>& expected) {
    const auto rows=table(output/"result.point_weights.csv"); require(rows.size()==expected.size(), "Point-weight export lost source records.");
    double sum=0;
    for(size_t i=0;i<rows.size();++i) {
        const double got=number(rows[i],"weight");
        require(close_value(got,expected[i]), "Incorrect occupied-cell point weight at source row "+std::to_string(i)); sum+=got;
        if(expected[i]==0) require(number(rows[i],"informed")==0, "Missing property counted as informed.");
    }
    require(close_value(sum,1), "Point weights do not sum to one.");
}
void check_tuple_audit(const fs::path& output) {
    const auto points=table(output/"result.point_weights.csv"); std::map<int,double> weights;
    for(const auto& row:points) weights[static_cast<int>(number(row,"node_id"))]=number(row,"weight");
    const auto tuples=table(output/"result.tuples.csv");
    for(const auto& tuple:tuples) {
        const auto ids=split(tuple.at("node_ids"),';'); std::set<int> distinct;
        for(const auto& id:ids) distinct.insert(std::stoi(id));
        long double product=1; for(int id:distinct) product*=weights.at(id);
        const double expected=std::pow(static_cast<double>(product),1.0/distinct.size());
        require(close_value(number(tuple,"tuple_weight"),expected), "Tuple weight is not the geometric mean over distinct source records.");
    }
    const auto support=table(output/"result.support.csv");
    for(const auto& row:support) {
        double total=0,squared=0;size_t count=0;
        for(const auto& tuple:tuples) if(tuple.at("lag_cell")==row.at("lag_cell")) {
            const double w=number(tuple,"tuple_weight");total+=w;squared+=w*w;++count;
        }
        require(count==static_cast<size_t>(number(row,"support_count")), "Audit/support cardinality mismatch.");
        require(close_value(number(row,"weight_sum"),total), "Template weight sum is inconsistent with tuple audit.");
        require(close_value(number(row,"ESS"),squared>0?total*total/squared:0), "Template ESS is inconsistent with tuple audit.");
    }
}
}

int main(int argc,char** argv) {
    try {
        require(argc==2,"Usage: template_declustering_checks <template_cli>");
        const fs::path exe=fs::absolute(argv[1]);
        const fs::path root=fs::current_path()/"test_output"/"template_declustering_checks"; fs::create_directories(root);
        const double nan=std::numeric_limits<double>::quiet_NaN();
        // Five occupied cells with counts 4,2,1,1,1; the missing row does
        // not contribute to occupancy. Every complete nonzero-lag tuple is
        // specified explicitly, so the oracle cannot share a search error.
        const std::vector<Point> points={
            {-10.25,-30.5,-.75,1},{-8.75,-30.5,-.75,2},
            {-10.25,-10.5,-.75,3},{-8.75,-10.5,-.75,6},
            {-10.25,9.5,-.75,7},{-8.75,9.5,-.75,10},
            {-9.75,-30.0,-.75,9},{-9.45,-30.0,-.75,4},
            {-10.0,-30.25,-.75,nan},{30,100,4,18}};
        const std::vector<double> weights={.05,.1,.05,.1,.2,.2,.05,.05,0,.2};
        const fs::path source=root/"clustered.dat";write_points(source,points);
        std::vector<Sample> scalar;
        for(size_t i=0;i<points.size();++i) if(weights[i]>0) scalar.push_back({{points[i].value},weights[i]});
        const double pair_weight=std::sqrt(.05*.1);
        const std::vector<Sample> pairs={{{1,2},pair_weight},{{3,6},pair_weight},{{7,10},.2}};
        for(int order=2;order<=4;++order) for(const std::string algorithm:{"Bell","AVL"}) {
            const std::string name="order"+std::to_string(order)+algorithm;
            require(run(exe,config(root,name,source,algorithm,shape(order-1,3),weighting()))==0,"Weighted run failed: "+name);
            const auto rows=table(root/name/"result.support.csv");require(rows.size()==3,"Incorrect lag count.");
            require(number(rows[0],"support_count")==9&&close_value(number(rows[0],"value"),scalar_cumulant(scalar,order)),"Weighted zero-lag cumulant does not match nonzero-mean closed form.");
            require(number(rows[1],"support_count")==3&&close_value(number(rows[1],"value"),cumulant(pairs,{1,order-1})),"Weighted joint cumulant does not use one fixed complete-tuple support and weight set.");
            require(number(rows[2],"support_count")==0&&std::isnan(number(rows[2],"value"))&&
                rows[2].at("status")=="insufficient_support"&&number(rows[2],"weight_sum")==0&&number(rows[2],"ESS")==0,
                "Unsupported weighted lag must be missing with zero support/weight/ESS diagnostics.");
            check_point_weights(root/name,weights);check_tuple_audit(root/name);
        }
        // Weight construction must remain active when the optional tuple audit
        // is disabled; collecting row IDs only for diagnostics is insufficient.
        require(run(exe,config(root,"no_audit",source,"AVL",shape(3),weighting(),"",false))==0,"Weighting without tuple export failed.");
        require(close_value(number(table(root/"no_audit"/"result.support.csv")[1],"value"),cumulant(pairs,{1,3})),"Tuple export switch changes weighted calculation.");
        require(run(exe,config(root,"minimum_support",source,"AVL",shape(3)+"minimum_template_count=4\n",weighting()))==0,"Weighted minimum-support run failed.");
        const auto minimum_rows=table(root/"minimum_support"/"result.support.csv");
        require(minimum_rows[0].at("status")=="ok"&&minimum_rows[1].at("status")=="insufficient_support"&&
            number(minimum_rows[1],"support_count")==3&&std::isnan(number(minimum_rows[1],"value")),
            "Weighted minimum support was not applied independently to each lag.");
        check_tuple_audit(root/"minimum_support");
        for(int order=5;order<=9;++order) for(const std::string algorithm:{"Bell","AVL"}) {
            const std::string name="high"+std::to_string(order)+algorithm;
            require(run(exe,config(root,name,source,algorithm,shape(order-1,1),weighting()))==0,"High-order weighted zero-lag run failed.");
            require(close_value(number(table(root/name/"result.support.csv")[0],"value"),scalar_cumulant(scalar,order)),"High-order weighted cumulant differs from independent univariate recurrence.");
        }
        const std::string coincident="NOD=2,1\nNLAG=2,2\nLAG=1.5,1.5\nLTOL=0,0\nAZM=0,0\n";
        for(const std::string algorithm:{"Bell","AVL"}) {
            const std::string name="coincident"+algorithm;
            require(run(exe,config(root,name,source,algorithm,coincident,weighting()))==0,"Coincident target weighted run failed.");
            const auto rows=table(root/name/"result.support.csv");require(rows.size()==4,"Wrong multidirectional output count.");
            require(close_value(number(rows[0],"value"),scalar_cumulant(scalar,4)),"All-zero target slots changed point weights.");
            require(number(rows[3],"support_count")==3&&close_value(number(rows[3],"value"),cumulant(pairs,{1,3})),"Coincident slots or NOD reweighted a sample more than once.");
            std::vector<double> mixed={cumulant(pairs,{3,1}),cumulant(pairs,{2,2})};
            require((close_value(number(rows[1],"value"),mixed[0])&&close_value(number(rows[2],"value"),mixed[1]))||
                    (close_value(number(rows[1],"value"),mixed[1])&&close_value(number(rows[2],"value"),mixed[0])),"Dependent moments or cache leaked across lag supports.");
            check_tuple_audit(root/name);
        }
        // Three distinct participating locations exercise a truly multivariate
        // target: its lower-order dependencies must retain the full triple's
        // geometric weight rather than recomputing pair-specific weights.
        auto triples_points=points;
        triples_points.push_back({-7.25,-30.5,-.75,5});
        triples_points.push_back({-7.25,-10.5,-.75,1});
        triples_points.push_back({-7.25,9.5,-.75,12});
        const fs::path triples_source=root/"triples.dat";write_points(triples_source,triples_points);
        const std::vector<Sample> triples={{{1,2,5},.05},{{3,6,1},.05},{{7,10,12},std::cbrt(.2*.1*.1)}};
        for(const std::array<int,2> powers:std::vector<std::array<int,2>>{{{1,1}},{{1,2}},{{2,3}}}) for(const std::string algorithm:{"Bell","AVL"}) {
            const int first_power=powers[0],last_power=powers[1];
            const std::string name="triple"+std::to_string(first_power)+"_"+std::to_string(last_power)+algorithm;
            const std::string geometry="NOD="+std::to_string(first_power)+","+std::to_string(last_power)+"\nNLAG=2,2\nLAG=1.5,3\nLTOL=0,0\nAZM=0,0\n";
            require(run(exe,config(root,name,triples_source,algorithm,geometry,weighting()))==0,"Three-location weighted run failed.");
            const auto rows=table(root/name/"result.support.csv");
            require(rows.size()==4&&number(rows[3],"support_count")==3&&
                close_value(number(rows[3],"value"),cumulant(triples,{1,first_power,last_power})),"Multivariate dependencies recomputed reduced-support weights or used incomplete tuples.");
            check_tuple_audit(root/name);
        }
        require(run(exe,config(root,"raw",source,"AVL",shape(3),weighting(),"calculate_target=Moment\n"))==0,"Weighted raw-moment run failed.");
        require(close_value(number(table(root/"raw"/"result.support.csv")[1],"value"),static_cast<double>(moment(pairs,{1,3}))),"Weighted raw moment was centered or used an unweighted denominator.");
        // A cell lattice containing the entire sample supplies equal weights,
        // which must reduce to the uniform calculation with identical support.
        const auto equal_weights=weighting("cell_size=1000,1000,1000\nshifts_per_axis=1\norigin=-100,-100,-100\n");
        require(run(exe,config(root,"equal",source,"AVL",shape(3),equal_weights))==0,"Equal-weight declustering failed.");
        require(run(exe,config(root,"uniform",source,"AVL",shape(3),weighting("","uniform")))==0,"Uniform comparison failed.");
        const auto equal=table(root/"equal"/"result.support.csv"),uniform=table(root/"uniform"/"result.support.csv");
        for(size_t i=0;i<equal.size();++i) require(close_value(number(equal[i],"value"),number(uniform[i],"value")),"Equal cell weights do not reduce to uniform moments/cumulants.");
        const auto shifted_weights=weighting("cell_size=2,35,2\nshifts_per_axis=3\norigin=-11,-31,-1\n");
        for(const std::string name:{"shifted1","shifted2"}) require(run(exe,config(root,name,source,"AVL",shape(3),shifted_weights))==0,"Multi-origin declustering run failed.");
        const auto shifted1=table(root/"shifted1"/"result.point_weights.csv"),shifted2=table(root/"shifted2"/"result.point_weights.csv");
        for(size_t i=0;i<shifted1.size();++i) require(close_value(number(shifted1[i],"weight"),number(shifted2[i],"weight")),"Origin averaging is not reproducible.");
        check_tuple_audit(root/"shifted1");
        // A two-dimensional source requires only two physical cell dimensions.
        std::vector<Point> planar=points;for(auto& p:planar)p.z=0;const fs::path two=root/"planar.dat";write_points(two,planar);
        require(run(exe,config(root,"2d",two,"Bell",shape(3),weighting("cell_size=2,35\nshifts_per_axis=1\norigin=-11,-31\n"),"coordinate_dimensions=2\n"))==0,"2D declustering rejected a valid two-component cell size.");
        require(close_value(number(table(root/"2d"/"result.support.csv")[1],"value"),cumulant(pairs,{1,3})),"2D weighted joint cumulant changed unexpectedly.");
        const std::vector<std::string> invalid={"", "cell_size=2,35\n", "cell_size=2,0,2\n", "cell_size=-2,35,2\n", "cell_size=nan,35,2\n", "cell_size=2,35,2\nshifts_per_axis=0\n", "cell_size=2,35,2\nshifts_per_axis=1.5\n", "cell_size=2,35,2\norigin=0,0\n", "cell_size=2,35,2\nshifts_per_axis=17\n", "cell_size=2units,35,2\n"};
        for(size_t i=0;i<invalid.size();++i) require(run(exe,config(root,"invalid"+std::to_string(i),source,"AVL",shape(1),weighting(invalid[i])))!=0,"Malformed declustering parameters silently accepted.");
        const fs::path missing_source=root/"all_missing.dat";write_points(missing_source,{{0,0,0,nan},{1.5,0,0,nan}});
        require(run(exe,config(root,"all_missing",missing_source,"AVL",shape(1),weighting()))!=0,"All-missing property silently generated declustering statistics.");
        const fs::path cartesian_source=root/"regular.grid";
        {std::ofstream out(cartesian_source);out<<"Cartesian control\n6 4\n# INDEX X Y Z VAL\n# Z=0\n1 0 0 0 1\n2 1 0 0 2\n3 2 0 0 4\n4 0 1 0 3\n5 1 1 0 5\n6 2 1 0 8\n";}
        const std::string cartesian_shape="template_sampling=exhaustive_integer_offsets\nNOD=1\nNLAG=2\nLAG=1\nLTOL=0\nAZM=0\nATOL=0\nBANDWH=0\nDIP=0\nDTOL=0\nBANDWD=0\ninteger_direction_vectors=1,0\n";
        const std::string cartesian_general="source_type=cartesian\ncoordinate_dimensions=2\nsource_property=value\n";
        const std::string cartesian_display="\n[target_cross_sections]\ndisplayed_dirs=1\nfixed_positions=\n";
        require(run(exe,config(root,"cartesian_control",cartesian_source,"AVL",cartesian_shape,weighting("","uniform")+cartesian_display,cartesian_general))==0,"Cartesian uniform control is invalid.");
        require(run(exe,config(root,"cartesian_declustering",cartesian_source,"AVL",cartesian_shape,weighting("cell_size=100,100\nshifts_per_axis=1\norigin=0,0\n")+cartesian_display,cartesian_general))==0,"Cartesian cell declustering rejected a valid equal-cell control.");
        check_point_weights(root/"cartesian_declustering",std::vector<double>(6,1.0/6.0));
        const auto cartesian_rows=table(root/"cartesian_declustering"/"result.support.csv");
        require(cartesian_rows.size()==2&&number(cartesian_rows[0],"support_count")==6&&
            close_value(number(cartesian_rows[0],"value"),scalar_cumulant({{{1},1},{{2},1},{{4},1},{{3},1},{{5},1},{{8},1}},2)),
            "Cartesian equal-cell zero-lag covariance differs from the independent variance.");
        for(const std::string mode:{"coordinate_gaussian","soft_high_porosity","fixed_value_gaussian","adaptive_upper_tail","external"})
            require(run(exe,config(root,"removed_"+mode,source,"AVL",shape(1),weighting("",mode)))!=0,"Unsupported weighting mode silently executed.");
        for(const std::string extra:{"normalize_weights=false\n","enabled=false\n","enable_soft_high_porosity=true\n"})
            require(run(exe,config(root,"invalid_switch",source,"AVL",shape(1),weighting("cell_size=2,35,2\n","cell_declustering",extra)))!=0,"Conflicting weighting switches silently accepted.");
        // New methods feed exactly the same complete tuples into the existing
        // evaluators. Recompute every reported moment/cumulant from the audit
        // using the independent raw-moment partition oracle above.
        for (const std::string mode:{"polygonal_declustering","kriging_declustering"})
        for (const int dims:{2,3}) for (const int order:{4,5})
        for (const std::string algorithm:{"Bell","AVL"}) for (const bool raw:{false,true}) {
            const std::string name=mode+std::to_string(dims)+std::to_string(order)+algorithm+(raw?"raw":"k");
            std::string parameters="dimensions="+std::to_string(dims)+"\ndomain_min=-11,-31"+(dims==3?std::string(",-1"):"")+
                "\ndomain_max=31,101"+(dims==3?std::string(",5"):"")+"\n";
            if(mode=="kriging_declustering") parameters+="kriging_model=exponential\nkriging_nugget_fraction=0.5\nkriging_ranges=2,10"+
                (dims==3?std::string(",2"):"")+"\nkriging_integration_cells=20,30"+(dims==3?std::string(",5"):"")+"\n";
            require(run(exe,config(root,name,source,algorithm,shape(order-1,3),weighting(parameters,mode),raw?"calculate_target=Moment\n":""))==0,
                "New declustering route failed: "+name);
            check_tuple_audit(root/name);
            const auto records=table(root/name/"result.tuples.csv");
            const auto rows=table(root/name/"result.support.csv");
            const auto points_audit=table(root/name/"result.point_weights.csv");
            double total=0;for(const auto& row:points_audit)total+=number(row,"weight");
            require(close_value(total,1)&&number(points_audit[8],"weight")==0,"Invalid point weights/missing mask.");
            for(const auto& row:rows) {
                std::vector<Sample> collected;
                for(const auto& record:records) if(record.at("lag_cell")==row.at("lag_cell")) {
                    Sample sample;sample.weight=number(record,"tuple_weight");
                    for(const auto& value:split(record.at("values"),';'))sample.values.push_back(std::stod(value));
                    collected.push_back(sample);
                }
                if(collected.empty()) {require(std::isnan(number(row,"value")),"Empty lag reported a value");continue;}
                const double expected=raw?static_cast<double>(moment(collected,{1,order-1})):cumulant(collected,{1,order-1});
                require(close_value(number(row,"value"),expected),"New method statistic differs from independent raw-partition reference: "+name);
            }
            std::ifstream metadata(root/name/"result.meta.txt");
            const std::string text{std::istreambuf_iterator<char>(metadata),{}};
            require(text.find("weight_mode="+mode)!=std::string::npos&&text.find("declustering_cell_size=")==std::string::npos,
                "Metadata misidentifies the new weighting method.");
        }
        for(const std::string mode:{"polygonal_declustering","kriging_declustering"}) {
            for(const std::string extra:{"enabled=false\n","normalize_weights=false\n"})
                require(run(exe,config(root,"new_invalid_switch",source,"AVL",shape(3),weighting("",mode,extra)))!=0,"New method bypassed switches.");
            require(run(exe,config(root,"missing_domain",source,"AVL",shape(3),weighting("dimensions=2\n",mode)))!=0,"Domain silently inferred from sample bounds.");
        }
        std::cout << "Template three-method declustering process checks: PASS\n"; return 0;
    } catch(const std::exception& e) { std::cerr << "Template cell-declustering checks FAILED: " << e.what() << '\n'; return 1; }
}
