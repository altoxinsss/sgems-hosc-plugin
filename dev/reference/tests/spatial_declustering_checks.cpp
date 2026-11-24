#include "spatial_declustering.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>

using namespace without_sgems;
namespace {
void require(bool ok,const std::string& message) { if(!ok) throw std::runtime_error(message); }
void close(double a,double b,double tolerance=1e-10) {
    require(std::isfinite(a) && std::abs(a-b)<=tolerance,"Expected "+std::to_string(b)+", got "+std::to_string(a));
}
template<class F> void rejects(F operation,const std::string& contains) {
    try {operation();} catch(const std::exception& error) {
        require(std::string(error.what()).find(contains)!=std::string::npos,"Unexpected rejection: "+std::string(error.what())); return;
    }
    throw std::runtime_error("Expected rejection containing: "+contains);
}
PointSetData points(std::initializer_list<std::array<double,3>> locations,int dimensions=2) {
    PointSetData data;data.dimensions=dimensions;
    for (const auto& p:locations) data.nodes.push_back({static_cast<int>(data.nodes.size()+1),p[0],p[1],p[2],1.0+data.nodes.size(),true});
    return data;
}
DeclusteringDomain unit(int dimensions=2) { return {dimensions,{{0,0,0}},{{1,1,1}}}; }
KrigingDeclusteringOptions kriging(int dimensions=2) {
    KrigingDeclusteringOptions k;
    k.domain=unit(dimensions);k.covariance_model="exponential";k.ranges={{1,1,1}};
    k.integration_cells=dimensions==2?std::array<int,3>{{30,30,1}}:std::array<int,3>{{10,10,10}};
    return k;
}
void weights(const SpatialDeclusteringResult& r,const std::vector<double>& expected,double tolerance=1e-10) {
    require(r.weights.size()==expected.size(),"Lost source rows");
    for (std::size_t i=0;i<expected.size();++i) close(r.weights[i],expected[i],tolerance);
    close(std::accumulate(r.weights.begin(),r.weights.end(),0.0),1);
}
}
int main() {
    try {
        const auto d2=unit();
        // Bisector x=.35; unbounded hull cells are clipped at the actual domain.
        auto uneven=points({{{.1,.5,0}},{{.6,.5,0}}});
        weights(calculate_polygonal_declustering(uneven,d2),{.35,.65});
        // Oblique boundary x+y=1.25 cuts off a triangle of area .75^2/2.
        auto diagonal=points({{{.5,.5,0}},{{.75,.75,0}}});
        weights(calculate_polygonal_declustering(diagonal,d2),{.71875,.28125});
        auto rows=points({{{.1,.5,1}},{{.1,.5,2}},{{.6,.5,1}},{{.6,.5,2}}},3);
        rows.nodes.push_back({5,100,100,100,std::numeric_limits<double>::quiet_NaN(),false});
        weights(calculate_polygonal_declustering(rows,d2),{.175,.175,.325,.325,0});
        auto cube=points({{{.25,.25,.25}},{{.75,.25,.25}},{{.25,.75,.25}},{{.75,.75,.25}},
                          {{.25,.25,.75}},{{.75,.25,.75}},{{.25,.75,.75}},{{.75,.75,.75}}},3);
        weights(calculate_polygonal_declustering(cube,unit(3)),std::vector<double>(8,.125));
        // In a unit cube the plane x+y+z=1 leaves a tetrahedron of volume 1/6.
        auto tetra=points({{{1.0/6,1.0/6,1.0/6}},{{.5,.5,.5}}},3);
        weights(calculate_polygonal_declustering(tetra,unit(3)),{1.0/6,5.0/6});
        // True Euclidean metric must survive coordinate translation and common
        // scaling; unequal domain side lengths must not silently scale axes.
        auto translated=diagonal;auto bounds=d2;
        for(auto& p:translated.nodes) {p.x=1000+7*p.x;p.y=-2000+7*p.y;}
        bounds.lower={{1000,-2000,0}};bounds.upper={{1007,-1993,1}};
        weights(calculate_polygonal_declustering(translated,bounds),{.71875,.28125});
        auto wide=points({{{1,.5,0}},{{3,.5,0}}});auto rectangle=d2;rectangle.upper[0]=10;
        weights(calculate_polygonal_declustering(wide,rectangle),{.2,.8});
        auto boundary=points({{{0,0,0}},{{1,1,0}}});
        weights(calculate_polygonal_declustering(boundary,d2),{.5,.5});
        auto one=points({{{.2,.8,0}}});weights(calculate_polygonal_declustering(one,d2),{1});
        auto bad=d2;bad.upper[0]=0;
        rejects([&]{calculate_polygonal_declustering(one,bad);},"bounds");
        bad=d2;bad.upper[0]=.1;
        rejects([&]{calculate_polygonal_declustering(one,bad);},"outside");
        auto missing=one;missing.nodes[0].informed=false;
        rejects([&]{calculate_polygonal_declustering(missing,d2);},"informed samples");

        // Two-site exact OK solution: w0=1/2+(b0-b1)/(2*(1-C01)).
        // With a single y-midpoint coinciding with sample y, b is a 1D
        // exponential integral and can be evaluated independently in closed form.
        auto pair=points({{{.1,.5,0}},{{.4,.5,0}}});auto k=kriging();
        k.integration_cells={{10000,1,1}};
        const auto integral=[](double x){return (2-std::exp(-3*x)-std::exp(-3*(1-x)))/3;};
        const double w=.5+(integral(.1)-integral(.4))/(2*(1-std::exp(-.9)));
        auto result=calculate_kriging_declustering(pair,k);
        weights(result,{w,1-w},2e-9);
        require(result.kriging_max_residual<1e-12,"Kriging equations not satisfied");
        require(result.integration_point_count==10000,"Cubature count incorrect");
        // Property values must not affect weights; missing rows remain excluded.
        pair.nodes[0].value=-9e6;pair.nodes[1].value=8e6;
        weights(calculate_kriging_declustering(pair,k),result.weights);
        std::reverse(pair.nodes.begin(),pair.nodes.end());
        weights(calculate_kriging_declustering(pair,k),{1-w,w},2e-9);
        pair=points({{{.1,.5,0}},{{.4,.5,0}}});
        pair.nodes.push_back(pair.nodes[0]);pair.nodes.back().z=9;pair.dimensions=3;
        pair.nodes.push_back({4,100,100,100,0,false});
        weights(calculate_kriging_declustering(pair,k),{w/2,1-w,w/2,0},2e-9);
        k=kriging();k.nugget_fraction=1;
        weights(calculate_kriging_declustering(uneven,k),{.5,.5});
        for (const std::string model:{"exponential","spherical","gaussian"}) {
            k=kriging(3);k.covariance_model=model;k.nugget_fraction=.1;
            weights(calculate_kriging_declustering(cube,k),std::vector<double>(8,.125));
        }
        k=kriging();weights(calculate_kriging_declustering(one,k),{1});
        auto negative=points({{{.4,.5,0}},{{.5,.5,0}},{{.6,.5,0}}});
        k.covariance_model="gaussian";
        rejects([&]{calculate_kriging_declustering(negative,k);},"non-positive site weight");
        k=kriging();k.nugget_fraction=-.1;
        rejects([&]{calculate_kriging_declustering(one,k);},"nugget");
        k=kriging();k.ranges[0]=0;
        rejects([&]{calculate_kriging_declustering(one,k);},"ranges");
        k=kriging();k.integration_cells[0]=0;
        rejects([&]{calculate_kriging_declustering(one,k);},"positive counts");
        k=kriging();k.covariance_model="invalid";
        rejects([&]{calculate_kriging_declustering(one,k);},"kriging_model");
        k=kriging();k.covariance_model="gaussian";k.ranges={{1e12,1e12,1e12}};
        rejects([&]{calculate_kriging_declustering(uneven,k);},"singular or ill-conditioned");

        // Cross the former site cap with an independent equal-area result.
        PointSetData large_polygon;large_polygon.dimensions=2;
        for(int y=0;y<100;++y) for(int x=0;x<101;++x)
            large_polygon.nodes.push_back({static_cast<int>(large_polygon.nodes.size()+1),
                (x+.5)/101,(y+.5)/100,0,1,true});
        weights(calculate_polygonal_declustering(large_polygon,d2),
            std::vector<double>(10100,1.0/10100));

        // Pure nugget gives exact equal weights, while still executing the
        // dense global solver and >100 million covariance integration pairs.
        PointSetData large_kriging;large_kriging.dimensions=2;
        for(int y=0;y<41;++y) for(int x=0;x<50;++x)
            large_kriging.nodes.push_back({static_cast<int>(large_kriging.nodes.size()+1),
                (x+.5)/50,(y+.5)/41,0,1,true});
        k=kriging();k.nugget_fraction=1;k.integration_cells={{250,200,1}};
        weights(calculate_kriging_declustering(large_kriging,k),std::vector<double>(2050,1.0/2050));
        k.integration_cells={{1000001,1,1}};
        result=calculate_kriging_declustering(one,k);
        weights(result,{1});
        require(result.integration_point_count==1000001,"Large integration count was truncated");
        k=kriging(3);k.integration_cells.fill(std::numeric_limits<int>::max());
        rejects([&]{calculate_kriging_declustering(cube,k);},"addressable integer range");
        std::cout<<"Polygonal and global kriging analytical checks: PASS\n";
        return 0;
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n';return 1; }
}
