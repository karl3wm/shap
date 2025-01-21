#pragma once
#include "geometry.hpp"
#include <functional>

namespace shap {

// Base expression template for surfaces
template<typename Derived>
struct Surface {
    constexpr Point operator()(double u, double v) const {
        return static_cast<const Derived&>(*this)(u, v);
    }
    
    // First partial derivatives
    constexpr Point du(double u, double v) const {
        const double h = 1e-7;
        return (operator()(u + h, v) - operator()(u - h, v)) * (0.5 / h);
    }
    
    constexpr Point dv(double u, double v) const {
        const double h = 1e-7;
        return (operator()(u, v + h) - operator()(u, v - h)) * (0.5 / h);
    }
    
    // Metric tensor
    constexpr MetricTensor metric_tensor(double u, double v) const {
        Point du_vec = du(u, v);
        Point dv_vec = dv(u, v);
        
        return MetricTensor(
            du_vec.dot(du_vec),      // g11
            du_vec.dot(dv_vec),      // g12
            du_vec.dot(dv_vec),      // g21
            dv_vec.dot(dv_vec)       // g22
        );
    }
    
    // Normal vector
    constexpr Point normal(double u, double v) const {
        Point n = du(u, v).cross(dv(u, v));
        double len = std::sqrt(n.norm_squared());
        // Handle coordinate singularities
        if (len < 1e-10) {
            // For sphere poles, return appropriate normal
            Point p = operator()(u, v);
            return p.normalize();
        }
        return n * (1.0 / len);
    }
    
    // Gaussian curvature
    constexpr double gaussian_curvature(double u, double v) const {
        // Simplified calculation for demo
        // In practice, would need second derivatives
        return 1.0 / metric_tensor(u, v).determinant();
    }
};

// Helper for creating surfaces from lambdas
template<typename F>
struct ParametricSurface : Surface<ParametricSurface<F>> {
    F func;
    
    constexpr ParametricSurface(F f) : func(f) {}
    
    constexpr Point operator()(double u, double v) const {
        return func(u, v);
    }
};

template<typename F>
constexpr auto make_surface(F&& f) {
    return ParametricSurface<F>(std::forward<F>(f));
}

// Predefined surfaces
namespace surfaces {

constexpr auto sphere(double radius = 1.0) {
    return make_surface([radius](double u, double v) {
        // u: longitude [0, 2π]
        // v: latitude [0, π]
        return Point(
            radius * std::cos(u) * std::sin(v),
            radius * std::sin(u) * std::sin(v),
            radius * std::cos(v)
        );
    });
}

constexpr auto cube(double size = 1.0) {
    return make_surface([size](double u, double v) {
        // u,v in [0,1] for each face
        // Maps six unit squares to cube faces
        
        // Scale inputs to get face index and local coordinates
        double scaled_u = u * 6;
        int face = static_cast<int>(scaled_u);
        double local_u = scaled_u - face;
        
        switch(face) {
            case 0: // Front
                return Point(
                    size * (2 * local_u - 1),
                    size,
                    size * (2 * v - 1)
                );
            case 1: // Right
                return Point(
                    size,
                    size * (1 - 2 * local_u),
                    size * (2 * v - 1)
                );
            case 2: // Back
                return Point(
                    size * (1 - 2 * local_u),
                    -size,
                    size * (2 * v - 1)
                );
            case 3: // Left
                return Point(
                    -size,
                    size * (2 * local_u - 1),
                    size * (2 * v - 1)
                );
            case 4: // Top
                return Point(
                    size * (2 * local_u - 1),
                    size * (2 * v - 1),
                    size
                );
            default: // Bottom
                return Point(
                    size * (2 * local_u - 1),
                    size * (2 * v - 1),
                    -size
                );
        }
    });
}

} // namespace surfaces

} // namespace shap
