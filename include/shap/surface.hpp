#pragma once
#include "surface_point.hpp"
#include "metric.hpp"
#include <memory>
#include <functional>

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
    
    // Second partial derivatives (optional - for exact metric derivatives)
    virtual Point evaluate_duu(double u, double v) const { return Point(); }
    virtual Point evaluate_duv(double u, double v) const { return Point(); }
    virtual Point evaluate_dvv(double u, double v) const { return Point(); }
    
    // Check if surface provides exact second derivatives
    virtual bool has_second_derivatives() const { return false; }
    
    // Get complete geometric data at a point
    virtual SurfacePoint evaluate(double u, double v) const {
        Point pos = evaluate_position(u, v);
        Point du_vec = evaluate_du(u, v);
        Point dv_vec = evaluate_dv(u, v);
        Point n = du_vec.cross(dv_vec).normalize();
        
        return SurfacePoint(name, u, v, pos, n, du_vec, dv_vec);
    }
    
    // Metric tensor and Riemannian connection
    virtual Surface2DMetricTensor metric_tensor(double u, double v) const {
        if (!has_second_derivatives()) {
            // Create basic metric coefficient functions without derivatives
            auto g11 = [this](double u, double v) {
                auto du = evaluate_du(u, v);
                return du.dot(du);
            };
            
            auto g12 = [this](double u, double v) {
                auto du = evaluate_du(u, v);
                auto dv = evaluate_dv(u, v);
                return du.dot(dv);
            };
            
            auto g22 = [this](double u, double v) {
                auto dv = evaluate_dv(u, v);
                return dv.dot(dv);
            };
            
            return Surface2DMetricTensor(g11, g12, g12, g22);
        }
        
        // Create metric coefficient functions with exact derivatives
        auto g11 = [this](double u, double v) {
            auto du = evaluate_du(u, v);
            return du.dot(du);
        };
        
        auto g12 = [this](double u, double v) {
            auto du = evaluate_du(u, v);
            auto dv = evaluate_dv(u, v);
            return du.dot(dv);
        };
        
        auto g22 = [this](double u, double v) {
            auto dv = evaluate_dv(u, v);
            return dv.dot(dv);
        };
        
        // Exact derivatives of metric coefficients
        auto dg11_du = [this](double u, double v) {
            auto du = evaluate_du(u, v);
            auto duu = evaluate_duu(u, v);
            return 2.0 * du.dot(duu);
        };
        
        auto dg11_dv = [this](double u, double v) {
            auto du = evaluate_du(u, v);
            auto duv = evaluate_duv(u, v);
            return 2.0 * du.dot(duv);
        };
        
        auto dg12_du = [this](double u, double v) {
            auto du = evaluate_du(u, v);
            auto dv = evaluate_dv(u, v);
            auto duu = evaluate_duu(u, v);
            auto duv = evaluate_duv(u, v);
            return duu.dot(dv) + du.dot(duv);
        };
        
        auto dg12_dv = [this](double u, double v) {
            auto du = evaluate_du(u, v);
            auto dv = evaluate_dv(u, v);
            auto duv = evaluate_duv(u, v);
            auto dvv = evaluate_dvv(u, v);
            return duv.dot(dv) + du.dot(dvv);
        };
        
        auto dg22_du = [this](double u, double v) {
            auto dv = evaluate_dv(u, v);
            auto duv = evaluate_duv(u, v);
            return 2.0 * dv.dot(duv);
        };
        
        auto dg22_dv = [this](double u, double v) {
            auto dv = evaluate_dv(u, v);
            auto dvv = evaluate_dvv(u, v);
            return 2.0 * dv.dot(dvv);
        };
        
        return Surface2DMetricTensor(
            g11, g12, g12, g22,
            dg11_du, dg11_dv,
            dg12_du, dg12_dv,
            dg22_du, dg22_dv
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
