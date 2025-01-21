#pragma once
#include "point.hpp"
#include <string>

namespace shap {

// Surface point with parameter coordinates and geometric data
class SurfacePoint {
public:
    // Default constructor
    SurfacePoint()
        : surface_name(""), u(0), v(0),
          position(), normal(0,0,1),
          du(1,0,0), dv(0,1,0) {}
    
    // Full constructor
    SurfacePoint(
        const std::string& surface_name,
        double u, double v,
        const Point& position,
        const Vector& normal,
        const Vector& du,
        const Vector& dv
    ) : surface_name(surface_name),
        u(u), v(v),
        position(position),
        normal(normal.normalize()),
        du(du), dv(dv) {}
    
    // Get tangent vector in given parameter direction
    Vector tangent(double du_component, double dv_component) const {
        return (du * du_component + dv * dv_component).normalize();
    }

    std::string surface_name;  // Name of containing surface
    double u, v;              // Parameter coordinates
    Point position;           // 3D position
    Vector normal;            // Surface normal (normalized)
    Vector du, dv;           // Tangent vectors
};

} // namespace shap