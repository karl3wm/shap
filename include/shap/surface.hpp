#pragma once
#include "surface_point.hpp"
#include "metric.hpp"
#include <memory>
#include <functional>

namespace shap {

class SurfacePath;

// Function types for surface definition
using PositionFunction = std::function<Point(double u, double v)>;
using DerivativeFunction = std::function<Point(double u, double v)>;

// Geometric properties of a surface point
struct GeometricProperties {
    Point position;
    Point normal;
    Point du;
    Point dv;
    
    // Optional second derivatives
    Point duu;
    Point duv;
    Point dvv;
    bool has_second_derivatives = false;
};

// Abstract base class for parametric surfaces
class Surface {
public:
    virtual ~Surface() = default;
    
    // Surface name for identification
    std::string name;
    
    // Core evaluation method
    virtual SurfacePoint evaluate(double u, double v) const {
        auto props = compute_properties(u, v);
        return SurfacePoint(
            name, u, v,
            props.position,
            props.normal,
            props.du,
            props.dv
        );
    }
    
    // Geometric properties computation
    virtual GeometricProperties compute_properties(double u, double v) const = 0;
    
    // Surface classification for geodesic computation
    enum class SurfaceType {
        Smooth,      // Smooth surface (e.g. sphere) - use geodesic equations
        Developable, // Can be flattened (e.g. cylinder, cube face) - geodesics are straight lines
        NonSmooth    // Has sharp edges/corners - geodesics may be undefined at edges
    };
    
    // Get surface type for geodesic computation
    virtual SurfaceType surface_type() const = 0;
    
    // Create paths on surface
    virtual std::unique_ptr<SurfacePath> create_path(
        const SurfacePoint& start,
        const Vector& direction,
        double length
    ) const;
    
    // Factory method for function-based surface creation
    static std::shared_ptr<Surface> create(
        PositionFunction position_func,
        SurfaceType type = SurfaceType::Smooth
    );
    
    // Optional derivative specification
    static std::shared_ptr<Surface> create_with_derivatives(
        PositionFunction position_func,
        DerivativeFunction du_func,
        DerivativeFunction dv_func,
        SurfaceType type = SurfaceType::Smooth
    );
    
    // Optional full derivative specification including second derivatives
    static std::shared_ptr<Surface> create_with_all_derivatives(
        PositionFunction position_func,
        DerivativeFunction du_func,
        DerivativeFunction dv_func,
        DerivativeFunction duu_func,
        DerivativeFunction duv_func,
        DerivativeFunction dvv_func,
        SurfaceType type = SurfaceType::Smooth
    );

protected:
    // Helper to compute normal from derivatives
    static Point compute_normal(const Point& du, const Point& dv) {
        return du.cross(dv).normalize();
    }
};

} // namespace shap
