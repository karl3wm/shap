#pragma once
#include <array>
#include <cmath>

namespace shap {

// 2x2 metric tensor for Riemannian geometry operations
class MetricTensor {
public:
    MetricTensor(double g11, double g12, double g21, double g22)
        : g11(g11), g12(g12), g21(g21), g22(g22) {}
    
    // Convert tangent vector components between coordinate systems
    std::pair<double,double> raise_indices(double v1, double v2) const {
        double det = determinant();
        if (std::abs(det) < 1e-10) {
            return {v1, v2}; // Fallback for degenerate metric
        }
        return {
            (g22 * v1 - g12 * v2) / det,
            (-g21 * v1 + g11 * v2) / det
        };
    }
    
    // Compute Christoffel symbols for geodesic equations
    std::array<double,2> christoffel_first(int i, int j, int k) const;
    std::array<double,2> christoffel_second(int i) const;
    
    double determinant() const {
        return g11 * g22 - g12 * g21;
    }

private:
    double g11, g12, g21, g22;
};

} // namespace shap