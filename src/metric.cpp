#include "shap/metric.hpp"
#include <cmath>
#include <array>

namespace shap {

std::array<double,2> Surface2DMetricTensor::christoffel_first(
    int i, int j, int k, double u, double v, double h
) const {
    // First kind Christoffel symbols
    // Γ_ijk = 1/2 (∂_i g_jk + ∂_j g_ik - ∂_k g_ij)
    
    // Get partial derivatives using exact or numerical approximation
    double dg_i = (i == 0) ? dg_du(j,k, u,v, h) : dg_dv(j,k, u,v, h);
    double dg_j = (j == 0) ? dg_du(i,k, u,v, h) : dg_dv(i,k, u,v, h);
    double dg_k = (k == 0) ? dg_du(i,j, u,v, h) : dg_dv(i,j, u,v, h);
    
    return {0.5 * (dg_i + dg_j - dg_k), 0.0};
}

std::array<double,2> Surface2DMetricTensor::christoffel_second(
    int i, double u, double v, double epsilon
) const {
    // Second kind Christoffel symbols
    // Γ^i_jk = g^im Γ_mjk
    
    // Get inverse metric components
    double det = determinant(u,v);
    if (std::abs(det) < epsilon) {
        return {0.0, 0.0}; // Degenerate case
    }
    
    double g11_inv = g22_fn(u,v) / det;
    double g12_inv = -g12_fn(u,v) / det;
    double g21_inv = -g21_fn(u,v) / det;
    double g22_inv = g11_fn(u,v) / det;
    
    // Compute first kind symbols
    auto gamma_1 = christoffel_first(0, i, i, u, v);
    auto gamma_2 = christoffel_first(1, i, i, u, v);
    
    // Contract with inverse metric
    return {
        g11_inv * gamma_1[0] + g12_inv * gamma_2[0],
        g21_inv * gamma_1[0] + g22_inv * gamma_2[0]
    };
}

} // namespace shap