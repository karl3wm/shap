#pragma once
#include "surface_point.hpp"
#include "surface.hpp"
#include <memory>
#include <vector>

namespace shap {

// Type of transition between surface segments
enum class TransitionType {
    Linear,     // Linear interpolation between segments
    Circular,   // Circular arc transition
    Geodesic    // Follow geodesic curve on surface
};

// Abstract base class for parametric paths on surfaces
class SurfacePath {
public:
    virtual ~SurfacePath() = default;
    
    // Evaluate path at parameter t in [0,1]
    virtual SurfacePoint evaluate(double t) const = 0;
    
    // Get tangent vector at parameter t
    virtual Vector tangent(double t) const = 0;
    
    // Get normal vector at parameter t
    virtual Vector normal(double t) const = 0;
    
    // Create offset path at constant distance with specified transition type
    virtual std::unique_ptr<SurfacePath> offset(
        double distance,
        TransitionType transition = TransitionType::Circular
    ) const = 0;
};

// Geodesic curve between two points on a surface
class GeodesicCurve : public SurfacePath {
public:
    GeodesicCurve(
        std::shared_ptr<Surface> surface,
        const SurfacePoint& start,
        const SurfacePoint& end
    );
    
    // SurfacePath interface implementation
    SurfacePoint evaluate(double t) const override;
    Vector tangent(double t) const override;
    Vector normal(double t) const override {
        return evaluate(t).normal;
    }
    
    std::unique_ptr<SurfacePath> offset(
        double distance,
        TransitionType transition = TransitionType::Circular
    ) const override {
        throw std::runtime_error("Offset not supported for geodesic curves");
    }
    
private:
    std::shared_ptr<Surface> surface_;
    std::vector<Point> points_;    // Discretized geodesic curve points
    std::vector<Vector> tangents_; // Tangent vectors at points
    
    // Different computation methods based on surface type
    void compute_smooth_geodesic(
        const SurfacePoint& start,
        const SurfacePoint& end
    );
    
    void compute_developable_geodesic(
        const SurfacePoint& start,
        const SurfacePoint& end
    );
};

// Circular arc transition between surface segments
class CircularArc {
public:
    CircularArc(
        const SurfacePoint& start,
        const Vector& start_tangent,
        const SurfacePoint& end,
        const Vector& end_tangent,
        double radius
    );
    
    SurfacePoint evaluate(double t) const;
    Vector tangent(double t) const;
    Vector normal(double t) const;
    
private:
    Point center_;        // Arc center
    double radius_;       // Arc radius
    Vector normal_;       // Arc plane normal
    double start_angle_; // Start angle in arc plane
    double sweep_angle_; // Angular extent of arc
    
    // Local coordinate system for arc plane
    Vector x_axis_, y_axis_;
};

// Path segment on a single surface
class PathSegment {
public:
    PathSegment(
        std::shared_ptr<Surface> surface,
        double t_start, double t_end,
        double u_start, double u_end,
        double v_start, double v_end,
        const Vector& direction
    );
    
    SurfacePoint evaluate(double t) const;
    Vector tangent(double t) const;
    Vector normal(double t) const;
    
    std::shared_ptr<Surface> surface;
    double t_start, t_end;      // Path parameter range
    double u_start, u_end;      // Surface parameter range in u
    double v_start, v_end;      // Surface parameter range in v
    Vector direction;           // Direction in surface parameters
    
private:
    // Compute smooth normal field along segment
    void compute_normal_field();
    std::vector<Vector> normal_samples_;
};

// Transition between path segments
class PathTransition {
public:
    PathTransition(
        const PathSegment& seg1,
        const PathSegment& seg2,
        TransitionType type,
        double radius = 0.1
    );
    
    SurfacePoint evaluate(double t) const;
    Vector tangent(double t) const;
    Vector normal(double t) const;
    
private:
    TransitionType type_;
    std::unique_ptr<GeodesicCurve> geodesic_;
    std::unique_ptr<CircularArc> arc_;
    
    // Linear transition data
    SurfacePoint start_, end_;
    Vector start_tangent_, end_tangent_;
};

// Path that follows surfaces with transitions
class TransitionPath : public SurfacePath {
public:
    void add_segment(
        std::shared_ptr<Surface> surface,
        double t_start, double t_end,
        double u_start, double u_end,
        double v_start, double v_end,
        const Vector& direction
    );
    
    // SurfacePath interface implementation
    SurfacePoint evaluate(double t) const override;
    Vector tangent(double t) const override;
    Vector normal(double t) const override;
    
    std::unique_ptr<SurfacePath> offset(
        double distance,
        TransitionType transition = TransitionType::Circular
    ) const override;

private:
    std::vector<PathSegment> segments_;
    std::vector<PathTransition> transitions_;
    
    // Helper to create transitions between segments
    void update_transitions(TransitionType type);
    
    // Compute smooth normal field along entire path
    void compute_normal_field();
    
    // Project point back to surface while maintaining offset
    static std::pair<double,double> project_to_surface_with_offset(
        std::shared_ptr<Surface> surface,
        const Point& point,
        double offset_distance
    );
};

} // namespace shap