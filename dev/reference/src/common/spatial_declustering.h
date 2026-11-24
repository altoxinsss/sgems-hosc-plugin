#ifndef WITHOUT_SGEMS_SPATIAL_DECLUSTERING_H
#define WITHOUT_SGEMS_SPATIAL_DECLUSTERING_H

#include "point_set.h"
#include <array>
#include <string>
#include <vector>

namespace without_sgems {

// Explicit physical study-domain bounds, not automatically inferred sample bounds.
// dimensions=2 on a 3-D well set groups records by XY and splits each site's
// mass equally among its informed records; dimensions=3 uses XYZ sites.
struct DeclusteringDomain {
    int dimensions = 0;
    std::array<double, 3> lower{{0, 0, 0}};
    std::array<double, 3> upper{{0, 0, 0}};
};

struct KrigingDeclusteringOptions {
    DeclusteringDomain domain;
    std::string covariance_model; // spherical, exponential or gaussian
    std::array<double, 3> ranges{{0, 0, 0}};
    double nugget_fraction = 0;
    // Midpoint cubature for mean covariance to the whole study domain.
    // Exponential/gaussian use exp(-3h) / exp(-3h^2), respectively.
    std::array<int, 3> integration_cells{{0, 0, 1}};
};

struct SpatialDeclusteringResult {
    std::vector<double> weights;
    std::size_t valid_count = 0;
    std::size_t unique_site_count = 0;
    double weighted_mean = 0;
    double partition_relative_error = 0;
    std::size_t integration_point_count = 0;
    double kriging_lagrange_multiplier = 0;
    double kriging_max_residual = 0;
    double kriging_min_cholesky_pivot = 0;
    double min_site_weight = 0;
    double max_site_weight = 0;
};

// Exact half-space clipping of Euclidean Voronoi cells by a rectangle (2-D)
// or box (3-D). Weight is clipped area/volume divided by domain area/volume.
SpatialDeclusteringResult calculate_polygonal_declustering(
    const PointSetData& data, const DeclusteringDomain& domain);

// Global ordinary BLOCK kriging of the domain mean, using all distinct sites:
// C*w + 1*lambda = c_D, 1^T*w = 1. c_D is midpoint-integrated over the domain.
// No local neighborhood, covariance fitting, diagonal jitter, clipping or
// absolute-value transformation is applied. Non-positive weights are rejected
// explicitly because the template geometric-mean adapter requires positivity.
SpatialDeclusteringResult calculate_kriging_declustering(
    const PointSetData& data, const KrigingDeclusteringOptions& options);

}
#endif
