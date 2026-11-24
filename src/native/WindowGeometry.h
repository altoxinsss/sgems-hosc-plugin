#ifndef HOSC_WINDOW_GEOMETRY_H
#define HOSC_WINDOW_GEOMETRY_H

// Shared C++03-compatible geometry for the CLI and in-process SGeMS core.
#include <algorithm>
#include <climits>
#include <cmath>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace hosc_window_geometry {
inline bool finite(double x) {
    return x == x && x <= std::numeric_limits<double>::max() &&
           x >= -std::numeric_limits<double>::max();
}
inline void triple(const std::vector<double>& v, const char* name) {
    if (v.size() != 3) throw std::runtime_error(std::string(name) + " requires X,Y,Z values.");
    for (int k=0;k<3;++k)
        if (!finite(v[k])) throw std::runtime_error(std::string(name) + " must be finite.");
}
struct Shape {
    bool box;
    double radii[3], angles[3], axes[3][3], extent[3];
    Shape(const std::string& name, const std::vector<double>& r,
          const std::vector<double>& a, int dimensions) : box(name == "box") {
        if (!box && name != "ellipsoid") throw std::runtime_error("neighborhood_shape must be box or ellipsoid.");
        triple(r,"neighborhood_radii"); triple(a,"neighborhood_angles");
        if (dimensions != 2 && dimensions != 3) throw std::runtime_error("Window dimensions must be 2 or 3.");
        if (r[0]<=0 || r[1]<=0 || r[2]<0 || (dimensions==3 && r[2]==0))
            throw std::runtime_error("Half-sizes/radii must be positive; only 2D allows Z=0.");
        if (dimensions==2 && (a[1]!=0 || a[2]!=0)) throw std::runtime_error("2D requires zero dip and rake.");
        for (int k=0;k<3;++k) { radii[k]=r[k]; angles[k]=a[k]; }
        const double rad=acos(-1.0)/180, az=(90-a[0])*rad, dip=-a[1]*rad, rake=a[2]*rad;
        const double ca=cos(az),sa=sin(az),cb=cos(dip),sb=sin(dip),cc=cos(rake),sc=sin(rake);
        axes[0][0]=ca*cb; axes[0][1]=sa*cb; axes[0][2]=-sb;
        axes[1][0]=ca*sb*sc-sa*cc; axes[1][1]=sa*sb*sc+ca*cc; axes[1][2]=cb*sc;
        axes[2][0]=ca*sb*cc+sa*sc; axes[2][1]=sa*sb*cc-ca*sc; axes[2][2]=cb*cc;
        for (int j=0;j<3;++j) {
            double sum=0;
            for (int k=0;k<3;++k) { double p=r[k]*axes[k][j]; sum+=box?fabs(p):p*p; }
            extent[j]=box?sum:sqrt(sum);
            if (!finite(extent[j])) throw std::runtime_error("Window extent overflow.");
        }
    }
    template<class Point> bool contains(const Point& p,const Point& center) const {
        const double delta[3]={p.x-center.x,p.y-center.y,p.z-center.z};
        double distance=0;
        for (int k=0;k<3;++k) {
            double local=0;
            for (int j=0;j<3;++j) local+=axes[k][j]*delta[j];
            if (radii[k]==0) {
                if (fabs(local)>1e-12*std::max(1.0,std::max(radii[0],radii[1]))) return false;
            } else if (box) {
                if (fabs(local)/radii[k]>1+1e-12) return false;
            } else { double scaled=local/radii[k]; distance+=scaled*scaled; }
        }
        return box || distance<=1+1e-12;
    }
};

struct CenterGrid {
    bool regular;
    int dims[3];
    double lower[3],upper[3],start[3],origin[3],spacing[3];
    size_t count;
    CenterGrid() : regular(false),count(0) {
        for(int k=0;k<3;++k){dims[k]=1;lower[k]=upper[k]=start[k]=origin[k]=0;spacing[k]=1;}
    }
    template<class Data> void configure(const Data& data, const Shape& shape,
        const std::string& centers, const std::string& bounds_mode,
        const std::vector<double>& lo,const std::vector<double>& hi,
        const std::vector<double>& first,const std::vector<double>& step) {
        if (centers!="regular_grid" && centers!="sample_points")
            throw std::runtime_error("center_mode must be regular_grid or sample_points.");
        regular=centers=="regular_grid";
        if (!regular) { count=data.nodes.size(); return; }
        if (data.nodes.empty()) throw std::runtime_error("Cannot bound an empty point set.");
        if (bounds_mode!="auto" && bounds_mode!="manual") throw std::runtime_error("grid_bounds_mode must be auto or manual.");
        if (bounds_mode=="manual") {triple(lo,"grid_bounds_min");triple(hi,"grid_bounds_max");}
        triple(step,"grid_spacing"); if(!first.empty()) triple(first,"grid_start");
        lower[0]=upper[0]=data.nodes[0].x;lower[1]=upper[1]=data.nodes[0].y;lower[2]=upper[2]=data.nodes[0].z;
        for(size_t i=0;i<data.nodes.size();++i){
            const double p[3]={data.nodes[i].x,data.nodes[i].y,data.nodes[i].z};
            for(int k=0;k<3;++k){if(!finite(p[k]))throw std::runtime_error("Point coordinates must be finite.");lower[k]=std::min(lower[k],p[k]);upper[k]=std::max(upper[k],p[k]);}
        }
        count=1;
        for(int k=0;k<3;++k){
            if(bounds_mode=="manual"){lower[k]=lo[k];upper[k]=hi[k];}
            if(upper[k]<lower[k])throw std::runtime_error("Grid upper bounds must not be below lower bounds.");
            start[k]=first.empty()?lower[k]:first[k]; spacing[k]=step[k];
            if(k==2 && data.dimensions==2){
                if(lower[k]!=0 || upper[k]!=0 || start[k]!=0 || step[k]!=0)
                    throw std::runtime_error("2D grid bounds/start Z and spacing Z must be zero.");
                dims[k]=1;origin[k]=0;continue;
            }
            if(step[k]<=0)throw std::runtime_error("Grid spacing must be positive on active axes.");
            if(start[k]<lower[k])throw std::runtime_error("First window corner must lie inside the bounding box.");
            const double end=upper[k]-2*shape.extent[k];
            const double tol=1e-12*std::max(1.0,std::max(fabs(upper[k]),fabs(start[k])));
            if(!finite(end) || start[k]>end+tol)throw std::runtime_error("No complete window fits in the chosen bounding box from grid_start.");
            double intervals=(end-start[k])/step[k];
            const double n=1+floor(std::max(0.0,intervals)+1e-10);
            if(!finite(n) || n>INT_MAX)throw std::runtime_error("Output grid dimension exceeds integer capacity.");
            dims[k]=(int)n; origin[k]=start[k]+shape.extent[k];
            if(!finite(origin[k]) || (dims[k]>1 && origin[k]+step[k]==origin[k]))
                throw std::runtime_error("Grid spacing is too small for the coordinate precision.");
            if(count>(size_t)INT_MAX/(size_t)dims[k])throw std::runtime_error("Output grid exceeds integer node capacity.");
            count*=dims[k];
        }
    }
    template<class Point> Point point(size_t i) const {
        Point p;p.index=(int)i;
        p.x=origin[0]+(i%dims[0])*spacing[0];
        p.y=origin[1]+((i/dims[0])%dims[1])*spacing[1];
        p.z=origin[2]+(i/((size_t)dims[0]*dims[1]))*spacing[2];
        p.informed=false;return p;
    }
    void metadata(std::ostream& out) const {
        out<<"center_mode="<<(regular?"regular_grid":"sample_points")<<'\n';
        out<<"output_geometry="<<(regular?"cartesian":"point_set")<<'\n';
        if(!regular)return;
        out<<"boundary_rule=complete_window\nscan_order=x_fastest_then_y_then_z\n";
        const char* names[]={"grid_bounds_min","grid_bounds_max","grid_start","grid_origin","grid_spacing"};
        const double* arrays[]={lower,upper,start,origin,spacing};
        for(int j=0;j<5;++j)out<<names[j]<<'='<<arrays[j][0]<<','<<arrays[j][1]<<','<<arrays[j][2]<<'\n';
        out<<"grid_dimensions="<<dims[0]<<','<<dims[1]<<','<<dims[2]<<'\n';
    }
};
}
#endif
