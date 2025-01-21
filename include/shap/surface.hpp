#pragma once
#include "surface_point.hpp"
#include "metric.hpp"
#include <memory>

namespace shap {

class SurfacePath;

// Abstract base class for parametric surfaces
class Surface {
public:
    virtual ~Surface() = default;
    
    // Surface name for identification
    std::string name;
    
    // Basic evaluation of surface point
    virtual Point evaluate_position(double u, double v) const = 0;
    
    // First partial derivatives
    virtual Point evaluate_du(double u, double v) const = 0;
    virtual Point evaluate_dv(double u, double v) const = 0;
    
    // Get complete geometric data at a point
    virtual SurfacePoint evaluate(double u, double v) const {
        Point pos = evaluate_position(u, v);
        Point du_vec = evaluate_du(u, v);
        Point dv_vec = evaluate_dv(u, v);
        Point n = du_vec.cross(dv_vec).normalize();
        
        return SurfacePoint(name, u, v, pos, n, du_vec, dv_vec);
    }
    
    // Metric tensor and Riemannian connection
    virtual MetricTensor metric_tensor(double u, double v) const {
        Point du_vec = evaluate_du(u, v);
        Point dv_vec = evaluate_dv(u, v);
        
        return MetricTensor(
            du_vec.dot(du_vec),      // g11
            du_vec.dot(dv_vec),      // g12
            du_vec.dot(dv_vec),      // g21
            dv_vec.dot(dv_vec)       // g22
        );
    }
    
    // Surface classification for geodesic computation
    enum class SurfaceType {
        Smooth,      // Smooth surface (e.g. sphere) - use geodesic equations
        Developable, // Can be flattened (e.g. cylinder, cube face) - geodesics are straight lines
        NonSmooth    // Has sharp edges/corners - geodesics may be undefined at edges
    };
    
    // Get surface type for geodesic computation
    virtual SurfaceType surface_type() const = 0;
    
    // For developable surfaces, get the development map (flattening)
    // Returns UV coordinates in the flattened space
    virtual std::pair<double,double> develop_point(double u, double v) const {
        return {u, v}; // Default is identity map
    }
    
    // Create paths on surface
    virtual std::unique_ptr<SurfacePath> create_geodesic(
        const SurfacePoint& start,
        const SurfacePoint& end
    ) const;
    
    virtual std::unique_ptr<SurfacePath> create_directional_path(
        const SurfacePoint& start,
        const Vector& direction,
        double length
    ) const;
    
    // Parallel transport a vector along a path
    virtual Vector parallel_transport(
        const Vector& v,
        const SurfacePath& path,
        double t_start,
        double t_end
    ) const;
};

} // namespace shap
