#include "shap/metric.hpp"
#include <cmath>

namespace shap {

std::array<double,2> MetricTensor::christoffel_first(int i, int j, int k) const {
    // First kind Christoffel symbols
    // Γ_ijk = 1/2 (∂_i g_jk + ∂_j g_ik - ∂_k g_ij)
    // For constant metric, all derivatives are zero
    return {0.0, 0.0};
}

std::array<double,2> MetricTensor::christoffel_second(int i) const {
    // Second kind Christoffel symbols
    // Γ^i_jk = g^im Γ_mjk
    // For constant metric, all symbols are zero
    return {0.0, 0.0};
}

} // namespace shap