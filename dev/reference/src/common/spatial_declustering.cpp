#include "spatial_declustering.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>

namespace without_sgems {
namespace {
using Vec = std::array<double, 3>;
using Face = std::vector<Vec>;
using Solid = std::vector<Face>;

Vec sub(const Vec& a, const Vec& b) { return {{a[0]-b[0], a[1]-b[1], a[2]-b[2]}}; }
Vec add(const Vec& a, const Vec& b) { return {{a[0]+b[0], a[1]+b[1], a[2]+b[2]}}; }
Vec mul(const Vec& a, double s) { return {{a[0]*s, a[1]*s, a[2]*s}}; }
double dot(const Vec& a, const Vec& b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
Vec cross(const Vec& a, const Vec& b) {
    return {{a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]}};
}
double norm(const Vec& a) { return std::hypot(a[0], a[1], a[2]); }

struct Sites {
    int dimensions;
    double scale;
    Vec extent;
    std::vector<Vec> positions; // shifted by domain lower and isotropically scaled
    std::vector<std::vector<std::size_t>> rows;
};

Sites prepare_sites(const PointSetData& data, const DeclusteringDomain& domain) {
    Sites sites{};
    sites.dimensions = domain.dimensions ? domain.dimensions : data.dimensions;
    if ((data.dimensions != 2 && data.dimensions != 3) ||
        (sites.dimensions != 2 && sites.dimensions != 3) || sites.dimensions > data.dimensions)
        throw std::runtime_error("Declustering dimensions must be 2 or 3 and no greater than input dimensions.");
    for (int k=0; k<sites.dimensions; ++k) {
        const double extent = domain.upper[k]-domain.lower[k];
        if (!std::isfinite(domain.lower[k]) || !std::isfinite(domain.upper[k]) ||
            !std::isfinite(extent) || !(extent>0))
            throw std::runtime_error("Declustering domain bounds must define a finite positive rectangle/box.");
        sites.extent[k]=extent;
        sites.scale=std::max(sites.scale,extent);
    }
    for (int k=0; k<sites.dimensions; ++k) sites.extent[k]/=sites.scale;
    std::map<Vec, std::vector<std::size_t>> groups;
    for (std::size_t i=0; i<data.nodes.size(); ++i) {
        const auto& p=data.nodes[i];
        if (!p.informed) continue;
        Vec key{{p.x,p.y,p.z}};
        if (!std::isfinite(p.value) || !std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z))
            throw std::runtime_error("Informed declustering records must have finite coordinates and values.");
        for (int k=0; k<sites.dimensions; ++k)
            if (key[k]<domain.lower[k] || key[k]>domain.upper[k])
                throw std::runtime_error("An informed sample lies outside the specified declustering domain.");
        for (int k=sites.dimensions; k<3; ++k) key[k]=0;
        groups[key].push_back(i);
    }
    if (groups.empty()) throw std::runtime_error("Declustering requires informed samples.");
    for (const auto& item:groups) {
        Vec p{{0,0,0}};
        for (int k=0; k<sites.dimensions; ++k) p[k]=(item.first[k]-domain.lower[k])/sites.scale;
        sites.positions.push_back(p);
        sites.rows.push_back(item.second);
    }
    return sites;
}

void distribute(const PointSetData& data, const Sites& sites, const std::vector<double>& masses,
    SpatialDeclusteringResult& result) {
    result.weights.assign(data.nodes.size(),0);
    result.unique_site_count=sites.positions.size();
    result.min_site_weight=*std::min_element(masses.begin(),masses.end());
    result.max_site_weight=*std::max_element(masses.begin(),masses.end());
    double value_scale=0;
    for (const auto& p:data.nodes) if (p.informed) value_scale=std::max(value_scale,std::abs(p.value));
    long double scaled_mean=0;
    for (std::size_t i=0; i<masses.size(); ++i) {
        const double w=masses[i]/static_cast<double>(sites.rows[i].size());
        if (!(w>0) || !std::isfinite(w))
            throw std::runtime_error("Declustering produced a non-positive or non-finite record weight.");
        for (const auto row:sites.rows[i]) {
            result.weights[row]=w;
            ++result.valid_count;
            if (value_scale>0) scaled_mean+=static_cast<long double>(w)*(data.nodes[row].value/value_scale);
        }
    }
    result.weighted_mean=static_cast<double>(std::clamp(scaled_mean,-1.0L,1.0L)*value_scale);
}

// Retain n.(x-midpoint)<=0. The midpoint form avoids cancellation from squared
// absolute coordinates, including for very close sites near a large origin.
Face clip_face(const Face& face, const Vec& n, const Vec& midpoint, Face* cuts=nullptr) {
    Face out;
    if (face.empty()) return out;
    for (std::size_t i=0; i<face.size(); ++i) {
        const Vec& a=face[i];
        const Vec& b=face[(i+1)%face.size()];
        const double da=dot(n,sub(a,midpoint)), db=dot(n,sub(b,midpoint));
        if (da<=0) out.push_back(a);
        if ((da<0 && db>0) || (da>0 && db<0)) {
            const Vec p=add(a,mul(sub(b,a),da/(da-db)));
            out.push_back(p);
            if (cuts) cuts->push_back(p);
        }
        if (cuts && da==0) cuts->push_back(a);
    }
    return out;
}

Solid box(const Vec& e) {
    const Vec a{{0,0,0}}, b{{e[0],0,0}}, c{{e[0],e[1],0}}, d{{0,e[1],0}};
    const Vec f{{0,0,e[2]}}, g{{e[0],0,e[2]}}, h{{e[0],e[1],e[2]}}, j{{0,e[1],e[2]}};
    return {{a,d,c,b},{f,g,h,j},{a,b,g,f},{d,j,h,c},{a,f,j,d},{b,c,h,g}};
}

Solid clip_solid(const Solid& solid, const Vec& n, const Vec& midpoint) {
    bool outside=false, inside=false;
    for (const auto& face:solid) for (const auto& p:face) {
        const double d=dot(n,sub(p,midpoint));
        outside=outside || d>0;
        inside=inside || d<0;
    }
    if (!outside) return solid;
    if (!inside) return {};
    Solid out;
    Face cuts;
    for (const auto& face:solid) {
        Face clipped=clip_face(face,n,midpoint,&cuts);
        if (clipped.size()>=3) out.push_back(std::move(clipped));
    }
    Face unique;
    for (const auto& p:cuts) {
        bool seen=false;
        for (const auto& q:unique) if (norm(sub(p,q))<1e-13) { seen=true; break; }
        if (!seen) unique.push_back(p);
    }
    if (unique.size()>=3) {
        Vec center{{0,0,0}};
        for (const auto& p:unique) center=add(center,p);
        center=mul(center,1.0/unique.size());
        Vec axis=std::abs(n[0])<0.8?Vec{{1,0,0}}:Vec{{0,1,0}};
        Vec u=cross(n,axis); u=mul(u,1.0/norm(u));
        const Vec v=cross(n,u);
        std::sort(unique.begin(),unique.end(),[&](const Vec& a,const Vec& b) {
            const Vec da=sub(a,center), db=sub(b,center);
            return std::atan2(dot(da,v),dot(da,u))<std::atan2(dot(db,v),dot(db,u));
        });
        out.push_back(std::move(unique));
    }
    return out;
}

double area(const Face& face) {
    if (face.size()<3) return 0;
    long double result=0;
    for (std::size_t i=1; i+1<face.size(); ++i)
        result+=norm(cross(sub(face[i],face[0]),sub(face[i+1],face[0])))/2;
    return static_cast<double>(result);
}
double volume(const Solid& solid) {
    Vec center{{0,0,0}}; std::size_t count=0;
    for (const auto& face:solid) for (const auto& p:face) { center=add(center,p); ++count; }
    if (!count) return 0;
    center=mul(center,1.0/count);
    long double result=0;
    for (const auto& face:solid) for (std::size_t i=1; i+1<face.size(); ++i)
        result+=std::abs(dot(sub(face[0],center),cross(sub(face[i],center),sub(face[i+1],center))))/6;
    return static_cast<double>(result);
}

std::vector<double> cholesky_solve(const std::vector<double>& lower, std::vector<double> rhs) {
    const std::size_t n=rhs.size();
    for (std::size_t i=0; i<n; ++i) {
        long double v=rhs[i];
        for (std::size_t j=0; j<i; ++j) v-=static_cast<long double>(lower[i*n+j])*rhs[j];
        rhs[i]=static_cast<double>(v/lower[i*n+i]);
    }
    for (std::size_t i=n; i-->0;) {
        long double v=rhs[i];
        for (std::size_t j=i+1; j<n; ++j) v-=static_cast<long double>(lower[j*n+i])*rhs[j];
        rhs[i]=static_cast<double>(v/lower[i*n+i]);
    }
    return rhs;
}
}

SpatialDeclusteringResult calculate_polygonal_declustering(const PointSetData& data,
    const DeclusteringDomain& domain) {
    const Sites sites=prepare_sites(data,domain);
    double domain_measure=1;
    for (int k=0; k<sites.dimensions; ++k) domain_measure*=sites.extent[k];
    if (!(domain_measure>0) || !std::isfinite(domain_measure))
        throw std::runtime_error("Declustering domain aspect ratio is numerically unsupported.");
    std::vector<double> masses;
    for (std::size_t i=0; i<sites.positions.size(); ++i) {
        std::vector<std::pair<double,std::size_t>> neighbours;
        for (std::size_t j=0; j<sites.positions.size(); ++j) if (i!=j) {
            const double distance=norm(sub(sites.positions[j],sites.positions[i]));
            if (!(distance>0)) throw std::runtime_error("Distinct sites are numerically indistinguishable after domain scaling.");
            neighbours.emplace_back(distance,j);
        }
        std::sort(neighbours.begin(),neighbours.end());
        Face polygon{{{0,0,0}},{{sites.extent[0],0,0}},{{sites.extent[0],sites.extent[1],0}},{{0,sites.extent[1],0}}};
        Solid polyhedron;
        if (sites.dimensions==3) polyhedron=box(sites.extent);
        for (const auto& neighbour:neighbours) {
            // A bisector farther than all current vertices cannot cut this cell.
            double radius=0;
            if (sites.dimensions==2) for (const auto& p:polygon) radius=std::max(radius,norm(sub(p,sites.positions[i])));
            else for (const auto& f:polyhedron) for (const auto& p:f) radius=std::max(radius,norm(sub(p,sites.positions[i])));
            if (neighbour.first>2*radius*(1+1e-13)) break;
            const auto& other=sites.positions[neighbour.second];
            const Vec normal=mul(sub(other,sites.positions[i]),1/neighbour.first);
            const Vec midpoint=add(sites.positions[i],mul(sub(other,sites.positions[i]),0.5));
            if (sites.dimensions==2) polygon=clip_face(polygon,normal,midpoint);
            else polyhedron=clip_solid(polyhedron,normal,midpoint);
        }
        const double measure=sites.dimensions==2?area(polygon):volume(polyhedron);
        if (!(measure>0) || !std::isfinite(measure))
            throw std::runtime_error("A clipped Voronoi cell has invalid area/volume; check domain and coordinate precision.");
        masses.push_back(measure/domain_measure);
    }
    const double total=std::accumulate(masses.begin(),masses.end(),0.0);
    SpatialDeclusteringResult result;
    result.partition_relative_error=std::abs(total-1);
    if (result.partition_relative_error>1e-8)
        throw std::runtime_error("Voronoi cells do not partition the specified domain within tolerance.");
    for (auto& mass:masses) mass/=total;
    distribute(data,sites,masses,result);
    return result;
}

SpatialDeclusteringResult calculate_kriging_declustering(const PointSetData& data,
    const KrigingDeclusteringOptions& options) {
    const Sites sites=prepare_sites(data,options.domain);
    const std::size_t n=sites.positions.size();
    // Check actual storage before forming n*n; no fixed site-count cap or
    // local approximation is applied to the global system.
    if (n>std::vector<double>().max_size()/n)
        throw std::runtime_error("Global kriging matrix exceeds the addressable vector size.");
    if (options.covariance_model!="spherical" && options.covariance_model!="exponential" && options.covariance_model!="gaussian")
        throw std::runtime_error("kriging_model must be spherical, exponential or gaussian.");
    if (!std::isfinite(options.nugget_fraction) || options.nugget_fraction<0 || options.nugget_fraction>1)
        throw std::runtime_error("kriging_nugget_fraction must be in [0,1].");
    Vec ranges{{1,1,1}};
    std::size_t integration_count=1;
    for (int k=0; k<3; ++k) {
        if (k>=sites.dimensions) {
            if (options.integration_cells[k]!=1)
                throw std::runtime_error("Inactive kriging integration axes must have one cell.");
            continue;
        }
        ranges[k]=options.ranges[k]/sites.scale;
        if (!(ranges[k]>0) || !std::isfinite(ranges[k]))
            throw std::runtime_error("Kriging ranges must be finite positive physical distances.");
        const int count=options.integration_cells[k];
        if (count<1)
            throw std::runtime_error("Kriging integration requires positive counts.");
        if (static_cast<std::size_t>(count)>std::numeric_limits<std::size_t>::max()/integration_count)
            throw std::runtime_error("Kriging integration node count exceeds the addressable integer range.");
        integration_count*=count;
    }
    const auto covariance=[&](const Vec& a,const Vec& b) {
        Vec d=sub(a,b);
        for (int k=0;k<sites.dimensions;++k) d[k]/=ranges[k];
        const double h=norm(d);
        double correlation=0;
        if (options.covariance_model=="spherical") correlation=h>=1?0:1-1.5*h+0.5*h*h*h;
        else if (options.covariance_model=="exponential") correlation=std::exp(-3*h);
        else correlation=std::exp(-3*h*h);
        // A point nugget has zero integral over a continuous volume, even if a
        // midpoint node happens to coincide with a sample. It enters C_ii only.
        return (1-options.nugget_fraction)*correlation;
    };
    std::vector<double> matrix(n*n,0), lower(n*n,0), block(n,0);
    for (std::size_t i=0;i<n;++i) for (std::size_t j=0;j<=i;++j)
        matrix[i*n+j]=matrix[j*n+i]=(i==j?1:covariance(sites.positions[i],sites.positions[j]));
    for (std::size_t i=0;i<n;++i) {
        long double sum=0;
        for (std::size_t code=0;code<integration_count;++code) {
            std::size_t remaining=code; Vec p{{0,0,0}};
            for (int k=0;k<sites.dimensions;++k) {
                const int count=options.integration_cells[k];
                p[k]=sites.extent[k]*(static_cast<double>(remaining%count)+0.5)/count;
                remaining/=count;
            }
            sum+=covariance(sites.positions[i],p);
        }
        block[i]=static_cast<double>(sum/integration_count);
    }
    SpatialDeclusteringResult result;
    result.integration_point_count=integration_count;
    result.kriging_min_cholesky_pivot=1;
    for (std::size_t i=0;i<n;++i) for (std::size_t j=0;j<=i;++j) {
        long double v=matrix[i*n+j];
        for (std::size_t k=0;k<j;++k) v-=static_cast<long double>(lower[i*n+k])*lower[j*n+k];
        if (i==j) {
            if (!(v>1e-12L) || !std::isfinite(v))
                throw std::runtime_error("Kriging covariance matrix is singular or ill-conditioned; no diagonal jitter was added.");
            result.kriging_min_cholesky_pivot=std::min(result.kriging_min_cholesky_pivot,static_cast<double>(v));
            lower[i*n+j]=std::sqrt(static_cast<double>(v));
        } else lower[i*n+j]=static_cast<double>(v/lower[j*n+j]);
    }
    const auto inv_block=cholesky_solve(lower,block);
    const auto inv_one=cholesky_solve(lower,std::vector<double>(n,1));
    const double sum_one=std::accumulate(inv_one.begin(),inv_one.end(),0.0);
    if (!(sum_one>0) || !std::isfinite(sum_one)) throw std::runtime_error("Invalid kriging constraint solve.");
    const double lambda=(std::accumulate(inv_block.begin(),inv_block.end(),0.0)-1)/sum_one;
    std::vector<double> masses(n);
    for (std::size_t i=0;i<n;++i) masses[i]=inv_block[i]-lambda*inv_one[i];
    const double total=std::accumulate(masses.begin(),masses.end(),0.0);
    result.kriging_lagrange_multiplier=lambda;
    result.kriging_max_residual=std::abs(total-1);
    for (std::size_t i=0;i<n;++i) {
        long double value=lambda-block[i];
        for (std::size_t j=0;j<n;++j) value+=static_cast<long double>(matrix[i*n+j])*masses[j];
        result.kriging_max_residual=std::max(result.kriging_max_residual,std::abs(static_cast<double>(value)));
    }
    if (!std::isfinite(result.kriging_max_residual) || result.kriging_max_residual>1e-8)
        throw std::runtime_error("Global kriging equations did not meet the residual tolerance.");
    const double minimum=*std::min_element(masses.begin(),masses.end());
    if (!(minimum>0) || !std::isfinite(minimum)) {
        std::ostringstream message;
        message<<std::setprecision(17)<<"Global ordinary block kriging produced a non-positive site weight (minimum="<<minimum
            <<"). The template geometric-mean rule requires positive weights. No clipping, absolute values or parameter adjustment was applied.";
        throw std::runtime_error(message.str());
    }
    for (auto& mass:masses) mass/=total;
    distribute(data,sites,masses,result);
    return result;
}
}
