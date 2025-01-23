#include "shap/coord.hpp"
#include "shap/metric.hpp"
#include <cmath>
#include <array>

namespace shap {

std::array<double,2> Surface2DMetricTensor::christoffel_first(
    int i, int j, int k, double u, double v, double h
) const noexcept {
    /**
     * First kind Christoffel symbols (Γ_ijk) tell us how the surface curves by measuring 
     * how the metric changes as we move in different directions.
     * 
     * Think of them as describing how "steep" the surface is in different directions:
     * - If they're zero, the surface is flat in that direction
     * - If they're non-zero, the surface is curved and vectors will change direction
     *   as they move along the surface
     * 
     * Formula: Γ_ijk = 1/2 (∂_i g_jk + ∂_j g_ik - ∂_k g_ij)
     * 
     * Where:
     * - ∂_i means "how much it changes as you move in direction i"
     * - g_jk are the metric components that measure distances and angles
     */
    
    if (i < 0 || i > 1 || j < 0 || j > 1 || k < 0 || k > 1) {
        return {0.0, 0.0};
    }

    // Get how the metric changes in each direction
    const double dg_i = (i == 0) ? dg_du(j,k, u,v, h) : dg_dv(j,k, u,v, h);
    const double dg_j = (j == 0) ? dg_du(i,k, u,v, h) : dg_dv(i,k, u,v, h);
    const double dg_k = (k == 0) ? dg_du(i,j, u,v, h) : dg_dv(i,j, u,v, h);
    
    // Combine the changes using the fundamental formula
    const double gamma = 0.5 * (dg_i + dg_j - dg_k);
    
    return {gamma, 0.0};
}

std::array<double,2> Surface2DMetricTensor::christoffel_second(
    int i, double u, double v, double epsilon
) const noexcept {
    /**
     * Second kind Christoffel symbols (Γ^i_jk) are used to compute geodesics 
     * (shortest paths) on the surface. They tell us how vectors change direction
     * when moving along the surface.
     * 
     * For example, on a sphere:
     * - Two particles starting parallel will appear to curve toward each other
     * - The Christoffel symbols quantify exactly how much they curve
     * 
     * Formula: Γ^i_jk = g^im Γ_mjk
     * 
     * Where:
     * - g^im are components of the inverse metric tensor
     * - Γ_mjk are first kind Christoffel symbols
     */
    
    if (i < 0 || i > 1) {
        return {0.0, 0.0};
    }
    
    // Check if the metric is well-behaved at this point
    const double det = determinant(u,v);
    if (std::abs(det) < epsilon) {
        return {0.0, 0.0}; // Surface is degenerate here
    }
    
    // Get the metric components at this point
    const double g11 = metric_fns_[0][0](u,v);
    const double g12 = metric_fns_[0][1](u,v);
    const double g21 = metric_fns_[1][0](u,v);
    const double g22 = metric_fns_[1][1](u,v);
    
    // Compute the inverse metric components
    const double inv_det = 1.0 / det;
    const double g11_inv = g22 * inv_det;  // (g22)/det
    const double g12_inv = -g12 * inv_det; // (-g12)/det
    const double g21_inv = -g21 * inv_det; // (-g21)/det
    const double g22_inv = g11 * inv_det;  // (g11)/det
    
    // Get the first kind symbols we need
    const auto gamma_1 = christoffel_first(0, i, i, u, v);
    const auto gamma_2 = christoffel_first(1, i, i, u, v);
    
    // Transform to second kind symbols using the inverse metric
    return {
        g11_inv * gamma_1[0] + g12_inv * gamma_2[0],  // u component
        g21_inv * gamma_1[0] + g22_inv * gamma_2[0]   // v component
    };
}

} // namespace shap
