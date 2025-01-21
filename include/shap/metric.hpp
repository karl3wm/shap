#pragma once
#include <array>
#include <cmath>
#include <functional>

namespace shap {

// 2x2 metric tensor for surface parameter spaces
class Surface2DMetricTensor {
public:
    // Constructor for constant coefficient metric
    Surface2DMetricTensor(double g11, double g12, double g21, double g22)
        : g11_fn([g11](double,double){ return g11; }),
          g12_fn([g12](double,double){ return g12; }),
          g21_fn([g21](double,double){ return g21; }),
          g22_fn([g22](double,double){ return g22; }),
          has_derivatives(false) {}

    // Constructor for variable coefficient metric with optional derivatives
    Surface2DMetricTensor(
        std::function<double(double,double)> g11,
        std::function<double(double,double)> g12,
        std::function<double(double,double)> g21,
        std::function<double(double,double)> g22,
        std::function<double(double,double)> dg11_du = nullptr,
        std::function<double(double,double)> dg11_dv = nullptr,
        std::function<double(double,double)> dg12_du = nullptr,
        std::function<double(double,double)> dg12_dv = nullptr,
        std::function<double(double,double)> dg22_du = nullptr,
        std::function<double(double,double)> dg22_dv = nullptr
    ) : g11_fn(g11), g12_fn(g12), g21_fn(g21), g22_fn(g22),
        dg11_du_fn(dg11_du), dg11_dv_fn(dg11_dv),
        dg12_du_fn(dg12_du), dg12_dv_fn(dg12_dv),
        dg22_du_fn(dg22_du), dg22_dv_fn(dg22_dv),
        has_derivatives(dg11_du && dg11_dv && dg12_du &&
                       dg12_dv && dg22_du && dg22_dv) {}

    // Get metric coefficient at given parameters
    double g(int i, int j, double u, double v) const {
        switch (i * 2 + j) {
            case 0: return g11_fn(u,v);
            case 1: return g12_fn(u,v);
            case 2: return g21_fn(u,v);
            case 3: return g22_fn(u,v);
            default: return 0.0;
        }
    }

    // Get partial derivative of metric coefficient
    // h: step size for numerical approximation when exact derivatives not available
    double dg_du(int i, int j, double u, double v, double h = 1e-7) const {
        if (!has_derivatives) {
            // Use numerical approximation if exact derivatives not provided
            return (g(i,j, u+h, v) - g(i,j, u-h, v)) / (2*h);
        }

        switch (i * 2 + j) {
            case 0: return dg11_du_fn(u,v);
            case 1: case 2: return dg12_du_fn(u,v);
            case 3: return dg22_du_fn(u,v);
            default: return 0.0;
        }
    }

    double dg_dv(int i, int j, double u, double v, double h = 1e-7) const {
        if (!has_derivatives) {
            // Use numerical approximation if exact derivatives not provided
            return (g(i,j, u, v+h) - g(i,j, u, v-h)) / (2*h);
        }

        switch (i * 2 + j) {
            case 0: return dg11_dv_fn(u,v);
            case 1: case 2: return dg12_dv_fn(u,v);
            case 3: return dg22_dv_fn(u,v);
            default: return 0.0;
        }
    }

    // Convert tangent vector components between coordinate systems
    // epsilon: tolerance for degenerate metric check
    std::pair<double,double> raise_indices(double v1, double v2, double u, double v, 
                                         double epsilon = 1e-10) const {
        double det = determinant(u,v);
        if (std::abs(det) < epsilon) {
            return {v1, v2}; // Fallback for degenerate metric
        }
        return {
            (g22_fn(u,v) * v1 - g12_fn(u,v) * v2) / det,
            (-g21_fn(u,v) * v1 + g11_fn(u,v) * v2) / det
        };
    }

    // Compute first kind Christoffel symbols for geodesic equations
    std::array<double,2> christoffel_first(int i, int j, int k, double u, double v, 
                                         double h = 1e-7) const;

    // Compute second kind Christoffel symbols for geodesic equations
    std::array<double,2> christoffel_second(int i, double u, double v, 
                                          double epsilon = 1e-10) const;

    // Compute determinant at given parameters
    double determinant(double u, double v) const {
        return g11_fn(u,v) * g22_fn(u,v) - g12_fn(u,v) * g21_fn(u,v);
    }

private:
    std::function<double(double,double)> g11_fn;
    std::function<double(double,double)> g12_fn;
    std::function<double(double,double)> g21_fn;
    std::function<double(double,double)> g22_fn;

    // Optional exact derivatives
    std::function<double(double,double)> dg11_du_fn;
    std::function<double(double,double)> dg11_dv_fn;
    std::function<double(double,double)> dg12_du_fn;
    std::function<double(double,double)> dg12_dv_fn;
    std::function<double(double,double)> dg22_du_fn;
    std::function<double(double,double)> dg22_dv_fn;

    bool has_derivatives;
};

} // namespace shap