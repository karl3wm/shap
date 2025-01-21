#pragma once
#include "surface.hpp"
#include "surface_point.hpp"
#include "geometry.hpp"
#include <vector>
#include <memory>

namespace shap {

// Base class for paths on surfaces
class SurfacePath {
public:
    virtual ~SurfacePath() = default;
    
    // Core path evaluation methods
    virtual SurfacePoint evaluate(double t) const = 0;
    virtual Vector tangent(double t) const = 0;
    virtual Vector normal(double t) const = 0;
};

// Geodesic curve between two points
class GeodesicCurve : public SurfacePath {
public:
    GeodesicCurve(
        std::shared_ptr<Surface> surface,
        const SurfacePoint& start,
        const SurfacePoint& end
    ) : surface_(surface) {
        if (surface->surface_type() == Surface::SurfaceType::Smooth) {
            compute_smooth_geodesic(start, end);
        } else if (surface->surface_type() == Surface::SurfaceType::Developable) {
            compute_developable_geodesic(start, end);
        } else {
            throw std::runtime_error("Cannot compute geodesic on non-smooth surface");
        }
    }
    
    SurfacePoint evaluate(double t) const override;
    Vector tangent(double t) const override;
    Vector normal(double t) const override;

private:
    void compute_smooth_geodesic(const SurfacePoint& start, const SurfacePoint& end);
    void compute_developable_geodesic(const SurfacePoint& start, const SurfacePoint& end);
    
    std::shared_ptr<Surface> surface_;
    std::vector<Point2D> points_;  // Points in parameter space (u,v)
    double t_start_;
    double t_end_;
};

// Single segment of a path on one surface
class PathSegment : public SurfacePath {
public:
    explicit PathSegment(std::shared_ptr<Surface> surface)
        : surface_(surface) {}
    
    void add_point(double t, double u, double v);
    
    SurfacePoint evaluate(double t) const override;
    Vector tangent(double t) const override;
    Vector normal(double t) const override;
    
    const std::vector<Point>& points() const { return points_; }
    std::shared_ptr<Surface> surface() const { return surface_; }

private:
    std::shared_ptr<Surface> surface_;
    std::vector<Point> points_;  // Points with (t,u,v) coordinates
};

// Path that can transition between surfaces
class TransitionPath : public SurfacePath {
public:
    void add_segment(
        std::shared_ptr<Surface> surface,
        double t_start, double t_end,
        double u_start, double u_end,
        double v_start, double v_end,
        const Vector& direction
    );
    
    SurfacePoint evaluate(double t) const override;
    Vector tangent(double t) const override;
    Vector normal(double t) const override;

private:
    std::vector<std::unique_ptr<PathSegment>> segments_;
};

} // namespace shap