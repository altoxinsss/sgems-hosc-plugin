#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
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
struct Point {double x,y,z,value;};
struct Result {size_t anchors=0,support=0;double value=0;std::string status;};

void require(bool ok,const std::string& message){if(!ok)throw std::runtime_error(message);}
bool close_value(double a,double b){return std::abs(a-b)<=1e-9*std::max({1.0,std::abs(a),std::abs(b)});}

int run(const fs::path& exe,const fs::path& config) {
#ifdef _WIN32
    std::wstring command=L"\""+exe.wstring()+L"\" \""+config.wstring()+L"\"";
    std::vector<wchar_t> buffer(command.begin(),command.end());buffer.push_back(0);
    STARTUPINFOW startup{};startup.cb=sizeof(startup);PROCESS_INFORMATION process{};
    if(!CreateProcessW(exe.c_str(),buffer.data(),nullptr,nullptr,FALSE,CREATE_NO_WINDOW,nullptr,nullptr,&startup,&process))
        throw std::runtime_error("Cannot start template_cli.");
    const DWORD waited=WaitForSingleObject(process.hProcess,60000);
    if(waited!=WAIT_OBJECT_0){TerminateProcess(process.hProcess,99);CloseHandle(process.hThread);CloseHandle(process.hProcess);throw std::runtime_error("template_cli test timed out.");}
    DWORD code=1;GetExitCodeProcess(process.hProcess,&code);CloseHandle(process.hThread);CloseHandle(process.hProcess);
    return static_cast<int>(code);
#else
    return std::system(("\""+exe.string()+"\" \""+config.string()+"\"").c_str());
#endif
}

void write_points(const fs::path& path,const std::vector<Point>& points) {
    std::ofstream out(path);out<<"physical point set\n4\nx\ny\nz\nphi\n"<<std::setprecision(17);
    for(const auto& p:points)out<<p.x<<' '<<p.y<<' '<<p.z<<' '<<p.value<<'\n';
}

std::vector<std::string> csv(const std::string& line) {
    std::vector<std::string> result;std::string value;bool quoted=false;
    for(char c:line){if(c=='\"')quoted=!quoted;else if(c==','&&!quoted){result.push_back(value);value.clear();}else value+=c;}
    result.push_back(value);return result;
}

std::vector<Result> results(const fs::path& path) {
    std::ifstream in(path);require(static_cast<bool>(in),"Missing support CSV.");
    std::string line;std::getline(in,line);std::vector<Result> result;
    while(std::getline(in,line)){const auto values=csv(line);require(values.size()==7,"Unexpected support schema.");result.push_back({static_cast<size_t>(std::stoull(values[3])),static_cast<size_t>(std::stoull(values[4])),std::stod(values[5]),values[6]});}
    return result;
}

std::string direction(int nod,const std::string& extra="") {
    return "NOD="+std::to_string(nod)+"\nNLAG=2\nLAG=1.5\nLTOL=0\nAZM=0\nATOL=0\nBANDWH=0\nDIP=0\nDTOL=0\nBANDWD=0\n"+extra;
}

fs::path config(const fs::path& root,const std::string& name,const fs::path& points,
    const std::string& algorithm,const std::string& shape,const std::string& general="",const std::string& extra="") {
    const fs::path path=root/(name+".ini");std::ofstream out(path);
    out<<"[general]\nsource_type=point_set\nsource_grid="<<points.string()<<"\nsource_property=phi\ncoordinate_dimensions=3\ncalculate_target=Cumulant\ncalculate_algorithm="<<algorithm<<"\ntarget_property=result\n"<<general
       <<"\n[template_shape]\ntemplate_sampling=point_set_tolerance\n"<<shape<<"\n[runtime]\noutput_dir="<<(root/name).string()<<"\nwrite_template_tuples=true\n"<<extra;
    return path;
}

std::string replace(std::string value,const std::string& from,const std::string& to){const size_t where=value.find(from);require(where!=std::string::npos,"Bad test replacement.");value.replace(where,from.size(),to);return value;}

// Independent low-order definitions: multiplicity (1,r-1), so powers act on
// the same target observation, never on successive spatial positions.
double manual(const std::vector<std::array<double,2>>& pairs,int order) {
    long double mx=0,my=0;
    for(const auto& p:pairs){mx+=p[0];my+=p[1];}mx/=pairs.size();my/=pairs.size();
    long double mixed=0,cov=0,var=0;
    for(const auto& p:pairs){const long double x=p[0]-mx,y=p[1]-my;mixed+=x*std::pow(y,order-1);cov+=x*y;var+=y*y;}
    mixed/=pairs.size();cov/=pairs.size();var/=pairs.size();
    return static_cast<double>(order==4?mixed-3*cov*var:mixed);
}

double univariate(const std::vector<Point>& points,int order) {
    long double mean=0;for(const auto& p:points)mean+=p.value;mean/=points.size();
    std::vector<long double> moments(order+1,0),cumulants(order+1,0);moments[0]=1;
    for(int r=1;r<=order;++r){for(const auto& p:points)moments[r]+=std::pow(p.value-mean,r);moments[r]/=points.size();}
    for(int r=1;r<=order;++r){cumulants[r]=moments[r];long double choose=1;for(int j=1;j<r;++j){choose*=static_cast<long double>(r-j)/j;cumulants[r]-=choose*cumulants[r-j]*moments[j];}}
    return static_cast<double>(cumulants[order]);
}

void audit_count(const fs::path& path,size_t expected) {
    std::ifstream in(path);std::string line;std::getline(in,line);size_t count=0;
    while(std::getline(in,line))++count;
    require(count==expected,"Tuple audit count does not match support.");
}

std::string contents(const fs::path& path) {
    std::ifstream input(path,std::ios::binary);require(static_cast<bool>(input),"Missing comparison file.");
    std::ostringstream result;result<<input.rdbuf();return result.str();
}
}

int main(int argc,char** argv) {
    try {
        require(argc==2,"Usage: template_pointset_checks <template_cli>");
        const fs::path exe=fs::absolute(argv[1]);
        const fs::path root=fs::current_path()/"test_output"/"template_pointset_checks";
        fs::create_directories(root);
        std::vector<Point> pairs;
        const double xvalues[]={1,3,7},yvalues[]={2,6,10};
        std::vector<std::array<double,2>> reference;
        for(int i=0;i<3;++i){pairs.push_back({-10.25,20*i-30.5,-0.75,xvalues[i]});pairs.push_back({-8.75,20*i-30.5,-0.75,yvalues[i]});reference.push_back({{xvalues[i],yvalues[i]}});}
        const fs::path points=root/"pairs.dat";write_points(points,pairs);
        // Nonzero means, negative/fractional physical coordinates, zero lag,
        // repeated variables, and both actual production evaluator routes.
        for(int order=2;order<=4;++order)for(const std::string algorithm:{"Bell","AVL"}){
            const std::string name="order"+std::to_string(order)+algorithm;
            require(run(exe,config(root,name,points,algorithm,direction(order-1)))==0,"Point-set run failed: "+name);
            const auto got=results(root/name/"result.support.csv");require(got.size()==2,"Wrong lag cell count.");
            require(got[1].support==3&&got[1].anchors==3&&close_value(got[1].value,manual(reference,order)),"Nonzero-lag cumulant differs from direct definition.");
            std::vector<std::array<double,2>> zero;for(const auto& p:pairs)zero.push_back({{p.value,p.value}});
            require(got[0].support==6&&close_value(got[0].value,manual(zero,order)),"Zero lag must reuse the anchor, including for repeated-variable orders.");
            audit_count(root/name/"result.tuples.csv",9);
        }
        for(int order=5;order<=9;++order)for(const std::string algorithm:{"Bell","AVL"}){
            const std::string name="zero_order"+std::to_string(order)+algorithm;
            require(run(exe,config(root,name,points,algorithm,replace(direction(order-1),"NLAG=2","NLAG=1")))==0,"High-order zero-lag run failed.");
            const auto got=results(root/name/"result.support.csv");
            require(got[0].support==6&&close_value(got[0].value,univariate(pairs,order)),"Bell/AVL differs from independent univariate recurrence at order "+std::to_string(order));
        }
        require(run(exe,config(root,"raw_moment",points,"AVL",direction(3),"calculate_target=Moment\n"))==0,"Raw moment path failed.");
        {const auto got=results(root/"raw_moment"/"result.support.csv");long double expected=0;for(const auto& p:reference)expected+=p[0]*std::pow(p[1],3);expected/=reference.size();require(close_value(got[1].value,static_cast<double>(expected)),"Moment target was incorrectly centered.");}
        require(run(exe,config(root,"minimum",points,"Bell",direction(1)+"minimum_template_count=4\n"))==0,"Minimum support path failed.");
        {const auto got=results(root/"minimum"/"result.support.csv");require(got[0].status=="ok"&&got[1].support==3&&std::isnan(got[1].value),"Minimum support was not applied per lag cell.");}
        std::vector<Point> planar=pairs;for(auto& p:planar)p.z=0;
        const fs::path two=root/"planar.dat";write_points(two,planar);
        require(run(exe,config(root,"planar",two,"AVL",direction(2),"coordinate_dimensions=2\n"))==0,"2-D point-set path failed.");
        {const auto got=results(root/"planar"/"result.support.csv");require(got[1].support==3&&close_value(got[1].value,manual(reference,3)),"2-D and 3-D planar statistics differ.");}
        // Candidate multiplicity is retained, rather than selecting one arbitrary
        // neighbor. One NaN candidate is excluded from every submoment support.
        std::vector<Point> branching;std::vector<std::array<double,2>> branches;
        for(int i=0;i<3;++i){branching.push_back(pairs[2*i]);branching.push_back({-8.77,20*i-30.5,-0.75,yvalues[i]});branching.push_back({-8.73,20*i-30.5,-0.75,yvalues[i]+1});branches.push_back({{xvalues[i],yvalues[i]}});branches.push_back({{xvalues[i],yvalues[i]+1}});}
        branching.push_back({-8.75,-30.5,-0.75,std::numeric_limits<double>::quiet_NaN()});
        const fs::path branched=root/"branches.dat";write_points(branched,branching);
        const std::string tolerance=replace(direction(3),"LTOL=0","LTOL=0.03");
        require(run(exe,config(root,"branches",branched,"AVL",tolerance))==0,"Candidate enumeration failed.");
        auto got=results(root/"branches"/"result.support.csv");
        require(got[1].support==6&&got[1].anchors==3&&close_value(got[1].value,manual(branches,4)),"All complete candidate tuples must share one common support.");
        // A cap failure must be explicit, never a partial statistical result.
        require(run(exe,config(root,"cap",branched,"AVL",tolerance+"max_template_tuples=2\n"))!=0,"Tuple cap silently truncated candidates.");
        require(run(exe,config(root,"visit_cap",branched,"AVL",tolerance+"max_candidate_visits=2\n"))!=0,"Candidate visit cap was ignored.");
        for(const auto& limits:std::vector<std::pair<std::string,std::string>>{
            {"unlimited","max_template_tuples=0\nmax_candidate_visits=0\n"},
            {"large_limits","max_template_tuples=10000001\nmax_candidate_visits=3000000000\n"}}) {
            require(run(exe,config(root,limits.first,branched,"AVL",tolerance+limits.second))==0,"Optional limit parsing failed.");
            require(contents(root/limits.first/"result.support.csv")==contents(root/"branches"/"result.support.csv"),
                "Changing workload limits changed statistical results.");
        }
        for(const std::string bad:{"-1","1.5","18446744073709551616"})
            require(run(exe,config(root,"invalid_limit",branched,"AVL",tolerance+"max_template_tuples="+bad+"\n"))!=0,
                "Invalid count limit was accepted.");
        // Full Cartesian product of two value collections: 1025^2 complete
        // tuples, >the former default cap; independent joint cumulant is zero.
        std::vector<Point> many;
        for(int i=0;i<1025;++i) many.push_back({0,0,0,static_cast<double>(i%5)});
        for(int i=0;i<1025;++i) many.push_back({1.5,0,0,static_cast<double>(i%7+1)});
        const fs::path large=root/"over_million.dat";write_points(large,many);
        require(run(exe,config(root,"over_million",large,"AVL",direction(3),"","write_template_tuples=false\n"))==0,
            "Default template enumeration still has the old million-tuple cap.");
        const auto large_result=results(root/"over_million"/"result.support.csv");
        require(large_result[1].support==1050625 && large_result[1].anchors==1025 && std::abs(large_result[1].value)<1e-9,
            "Large template collection was truncated or its independent-product cumulant is incorrect.");
        require(run(exe,config(root,"excluded",branched,"Bell",direction(1)))==0,"Empty lag support should produce missing status.");
        got=results(root/"excluded"/"result.support.csv");require(got[1].support==0&&std::isnan(got[1].value)&&got[1].status=="insufficient_support","Unsupported lag was reported as zero.");
        // Signed LAG reverses the template orientation without taking abs(x).
        require(run(exe,config(root,"reverse",points,"AVL",replace(direction(1),"LAG=1.5","LAG=-1.5")))==0,"Signed lag failed.");
        got=results(root/"reverse"/"result.support.csv");require(got[1].support==3&&close_value(got[1].value,manual(reference,2)),"Reverse covariance/support changed.");
        // Identical target slots must share one selected sample, even when two
        // candidates occur in the tolerance region. Different slots cannot.
        const std::string equal="NOD=1,1\nNLAG=2,2\nLAG=1.5,1.5\nLTOL=0.03,0.03\nAZM=0,0\n";
        require(run(exe,config(root,"coincident",branched,"AVL",equal))==0,"Coincident target handling failed.");
        got=results(root/"coincident"/"result.support.csv");require(got.size()==4&&got[3].support==6&&close_value(got[3].value,manual(branches,3)),"Coincident target slots produced distinct-node products or duplicate tuples.");
        require(run(exe,config(root,"cross_section",branched,"AVL",equal,"target_origin=100.5,-20.25,7.5\n","\n[target_cross_sections]\ndirection_x=2\ndirection_y=0\nfixed_positions=2\n"))==0,"Legacy cross-section aliases failed.");
        got=results(root/"cross_section"/"result.support.csv");require(got.size()==2&&got[1].support==6&&close_value(got[1].value,manual(branches,3)),"Fixed cross section changed complete tuple support.");
        {std::ifstream in(root/"cross_section"/"result.txt");std::string line;for(int i=0;i<7;++i)std::getline(in,line);std::istringstream row(line);double x,y,z,value;row>>x>>y>>z>>value;require(close_value(x,100.5)&&close_value(y,-20.25)&&close_value(z,7.5),"Output origin was ignored or source coordinates were changed.");}
        const std::string distinct="NOD=1,1\nNLAG=2,2\nLAG=1.49,1.51\nLTOL=0.02,0.02\nAZM=0,0\n";
        require(run(exe,config(root,"distinct",points,"AVL",distinct))==0,"Distinct target handling failed.");
        got=results(root/"distinct"/"result.support.csv");require(got[3].support==0,"One point was incorrectly reused for two distinct nonzero positions.");
        // Truly 3-D inclined lag. Distances are not rounded to integer cells.
        std::vector<Point> inclined;
        for(int i=0;i<3;++i){inclined.push_back(pairs[2*i]);inclined.push_back({-9.35,20*i-29.3,1.25,yvalues[i]});}
        const fs::path three=root/"inclined.dat";write_points(three,inclined);
        const std::string dip="NOD=2\nNLAG=2\nLAG=2.5\nLTOL=0.000000001\nAZM=53.13010235415598\nATOL=0.000000001\nDIP=53.13010235415598\nDTOL=0.000000001\n";
        require(run(exe,config(root,"inclined",three,"AVL",dip))==0,"Inclined 3-D template failed.");
        got=results(root/"inclined"/"result.support.csv");require(got[1].support==3&&close_value(got[1].value,manual(reference,3)),"3-D dip or physical distance is incorrect.");
        // Exact-target KD boxes must preserve the broad-radius sampler's full
        // support and input ordering, including near-boundary directions.
        // 1e-100 LTOL forces the original general query but rounds away when
        // added to its radius/epsilon, providing the same acceptance region.
        struct ExactCase {const char* name;double az,dip,sign,offset;int dimensions;};
        for(const auto& item:std::vector<ExactCase>{{"east",0,0,1,0,2},{"north",90,0,1,0,2},
                {"reverse",37,23,-1,0,3},{"oblique",37,23,1,0,3},{"pole",173,90,1,0,3},
                {"near_pole",173,89.99999999996,1,0,3},{"large_origin",0,0,1,1e8,3}}) {
            constexpr double rad=3.14159265358979323846/180;
            std::vector<Point> fixture;
            for(int i=0;i<3;++i) {
                const Point anchor{item.offset-10.25,20*i-30.5-item.offset,item.dimensions==2?0:item.offset-0.75,xvalues[i]};
                fixture.push_back(anchor);
                for(int k=0;k<4;++k) {
                    const double az=item.az+(k==1?0.5e-10:(k==2?2e-10:0));
                    const double distance=item.sign*(1.5+(k==3?3e-14:0));
                    fixture.push_back({anchor.x+distance*std::cos(item.dip*rad)*std::cos(az*rad),
                        anchor.y+distance*std::cos(item.dip*rad)*std::sin(az*rad),
                        anchor.z+distance*std::sin(item.dip*rad),yvalues[i]+k});
                }
                if(item.dip>89.9999999999) {
                    fixture.push_back({anchor.x+item.sign*1.5*std::cos(item.dip*rad)*std::cos(12*rad),
                        anchor.y+item.sign*1.5*std::cos(item.dip*rad)*std::sin(12*rad),
                        anchor.z+item.sign*1.5*std::sin(item.dip*rad),20.0+i});
                }
            }
            const fs::path fixture_path=root/(std::string("exact_")+item.name+".dat");write_points(fixture_path,fixture);
            std::ostringstream shape;shape<<std::setprecision(17)<<"NOD=2\nNLAG=2\nLAG="<<item.sign*1.5
                <<"\nLTOL=0\nAZM="<<item.az<<"\nATOL=0\nDIP="<<item.dip<<"\nDTOL=0\n";
            const std::string fast=std::string("exact_fast_")+item.name,slow=std::string("exact_slow_")+item.name;
            const std::string general="coordinate_dimensions="+std::to_string(item.dimensions)+"\n";
            require(run(exe,config(root,fast,fixture_path,"AVL",shape.str(),general))==0,"Exact lookup run failed.");
            require(run(exe,config(root,slow,fixture_path,"AVL",replace(shape.str(),"LTOL=0","LTOL=1e-100"),general))==0,"General lookup comparison failed.");
            require(contents(root/fast/"result.support.csv")==contents(root/slow/"result.support.csv"),"Exact lookup changed support or values: "+std::string(item.name));
            require(contents(root/fast/"result.tuples.csv")==contents(root/slow/"result.tuples.csv"),"Exact lookup changed tuple membership/order: "+std::string(item.name));
            require(results(root/fast/"result.support.csv")[1].support>=3,"Exact lookup comparison fixture has no useful support.");
        }
        // Both angle and bandwidth exclusions are active independent controls.
        std::vector<Point> bandpoints;for(int i=0;i<3;++i){bandpoints.push_back(pairs[2*i]);bandpoints.push_back({-8.75,20*i-30.35,-0.75,yvalues[i]});bandpoints.push_back({-8.75,20*i-30.15,-0.75,100.0+i});}
        const fs::path bandfile=root/"bands.dat";write_points(bandfile,bandpoints);
        std::string band=replace(replace(replace(direction(1),"LTOL=0","LTOL=0.1"),"ATOL=0","ATOL=20"),"BANDWH=0","BANDWH=0.2");
        require(run(exe,config(root,"bands",bandfile,"AVL",band))==0,"Bandwidth filtering failed.");
        got=results(root/"bands"/"result.support.csv");require(got[1].support==3&&close_value(got[1].value,manual(reference,2)),"Horizontal bandwidth did not exclude off-axis candidates.");
        require(run(exe,config(root,"angle",bandfile,"AVL",replace(band,"ATOL=20","ATOL=1")))==0,"Angle exclusion failed.");
        got=results(root/"angle"/"result.support.csv");require(got[1].support==0,"Azimuth tolerance was ignored.");
        std::vector<Point> verticalpoints;for(int i=0;i<3;++i){verticalpoints.push_back(pairs[2*i]);verticalpoints.push_back({-8.75,20*i-30.5,-0.60,yvalues[i]});verticalpoints.push_back({-8.75,20*i-30.5,-0.40,100.0+i});}
        const fs::path verticalfile=root/"vertical.dat";write_points(verticalfile,verticalpoints);
        const std::string vertical=replace(replace(replace(direction(1),"LTOL=0","LTOL=0.1"),"DTOL=0","DTOL=20"),"BANDWD=0","BANDWD=0.2");
        require(run(exe,config(root,"vertical",verticalfile,"Bell",vertical))==0,"Vertical bandwidth filter failed.");
        got=results(root/"vertical"/"result.support.csv");require(got[1].support==3&&close_value(got[1].value,manual(reference,2)),"Vertical bandwidth did not exclude off-axis candidates.");
        // Explicit legacy KEYS import and regular flag rejection.
        const fs::path keys=root/"test.KEYS";
        {std::ofstream out(keys);out<<"# dimension\n3\n# irregular\n0\n6 1 1 10 10 10\n2\n1\n2\n1.5\n0\n0 0 0 0 0 0\n";}
        require(run(exe,config(root,"keys",points,"Bell","keys_file="+keys.string()+"\n"))==0,"Legacy KEYS import failed.");
        got=results(root/"keys"/"result.support.csv");require(got[1].support==3&&close_value(got[1].value,manual(reference,2)),"KEYS changed template meaning.");
        {std::ofstream out(keys);out<<"3\n1\n6 1 1 10 10 10\n2\n1\n2\n1.5\n0\n0 0 0 0 0 0\n";}
        require(run(exe,config(root,"regular_keys",points,"Bell","keys_file="+keys.string()+"\n"))!=0,"Regular KEYS was silently interpreted as pointset.");
        for(const auto& bad:std::vector<std::string>{replace(direction(1),"LTOL=0","LTOL=-1"),replace(direction(1),"LAG=1.5","LAG=oops"),replace(direction(1),"DIP=0","DIP=100"),replace(direction(1),"NOD=1","NOD=10")})
            require(run(exe,config(root,"bad",points,"AVL",bad))!=0,"Malformed geometry was silently accepted.");
        require(run(exe,config(root,"removed",points,"AVL",direction(1),"","\n[weighted_moment]\nenabled=false\nweight_mode=soft_high_porosity\n"))!=0,"Removed weight silently fell back to uniform.");
        require(run(exe,config(root,"unsupported_center",points,"AVL",direction(1),"","\n[weighted_moment]\nzero_first_order_moment=true\n"))!=0,"Unsupported legacy centering was silently ignored.");
        for (const std::string transform : {"source_origin=0,0,0\n", "source_spacing=1,1,1\n"})
            require(run(exe,config(root,"source_transform",points,"AVL",direction(1),transform))!=0,"Source transform must not be silently ignored.");
        std::cout<<"Point-set template process checks: PASS\n";return 0;
    }catch(const std::exception& e){std::cerr<<"Point-set template checks FAILED: "<<e.what()<<'\n';return 1;}
}
