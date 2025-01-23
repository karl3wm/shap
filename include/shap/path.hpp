#include "coord.hpp"
#pragma once
#include "geometry_point2.hpp"
#include "surface.hpp"
#include <memory>
#include <vector>
#include <stdexcept>

namespace shap {

/**
 * Abstract base class for paths on surfaces.
 *
 * A path represents a curve that lies on one or more surfaces. The curve
 * is parameterized by t ∈ [0,1], where:
 * - t=0 corresponds to the start point
 * - t=1 corresponds to the end point (or last point before surface boundary)
 */
class SurfacePath {
public:
    virtual ~SurfacePath() = default;
    
    /**
     * Evaluate the path at parameter t.
     *
     * @param t Path parameter in [0,1]
     * @throws std::invalid_argument if t is outside [0,1]
     *
     * Post-conditions:
     * - Returns a point that lies on the surface(s)
     * - Position varies continuously with t
     * - For t=0, returns path start point
     * - For t=1, returns path end point
     */
    [[nodiscard]] virtual GeometryPoint2 evaluate(double t) const = 0;
    
    /**
     * Get path tangent vector at parameter t.
     *
     * @param t Path parameter in [0,1]
     * @throws std::invalid_argument if t is outside [0,1]
     *
     * Post-conditions:
     * - Returns normalized tangent vector
     * - Vector lies in surface tangent plane
     * - Direction matches path orientation
     */
    [[nodiscard]] virtual WorldVector3 tangent(double t) const = 0;
    
    /**
     * Get surface normal at parameter t.
     *
     * @param t Path parameter in [0,1]
     * @throws std::invalid_argument if t is outside [0,1]
     *
     * Post-conditions:
     * - Returns normalized surface normal vector
     * - Vector is perpendicular to path tangent
     */
    [[nodiscard]] virtual WorldVector3 normal(double t) const = 0;

protected:
    // Validate parameter t is in [0,1]
    static void validate_parameter(double t) {
        if (t < 0.0 || t > 1.0) {
            throw std::invalid_argument("Path parameter t must be in [0,1]");
        }
    }
};

/**
 * Geodesic curve between two points on a surface.
 *
 * A geodesic is a curve that locally minimizes path length on the surface.
 * For smooth surfaces, it follows the surface curvature.
 * For developable surfaces, it's a straight line in the developed plane.
 */
class GeodesicCurve final : public SurfacePath {
public:
    GeodesicCurve(
        std::shared_ptr<Surface> surface,
        const GeometryPoint2& start,
        const GeometryPoint2& end
    );
    
    [[nodiscard]] GeometryPoint2 evaluate(double t) const override;
    [[nodiscard]] WorldVector3 tangent(double t) const override;
    [[nodiscard]] WorldVector3 normal(double t) const override;

private:
    void compute_smooth_geodesic(
        const GeometryPoint2& start,
        const GeometryPoint2& end
    );
    
    void compute_developable_geodesic(
        const GeometryPoint2& start,
        const GeometryPoint2& end
    );
    
    std::shared_ptr<Surface> surface_;
    std::vector<GeometryPoint2> points_;
};

// Path segment on a single surface
class PathSegment final : public SurfacePath {
public:
    explicit PathSegment(std::shared_ptr<Surface> surface) noexcept
        : surface_(std::move(surface)) {
        // Pre-allocate space for typical path size
        t_values_.reserve(100);
        u_values_.reserve(100);
        v_values_.reserve(100);
    }
    
    // Move operations
    PathSegment(PathSegment&&) noexcept = default;
    PathSegment& operator=(PathSegment&&) noexcept = default;
    
    // Prevent copying
    PathSegment(const PathSegment&) = delete;
    PathSegment& operator=(const PathSegment&) = delete;
    
    void add_point(double t, double u, double v);
    
    [[nodiscard]] GeometryPoint2 evaluate(double t) const override;
    [[nodiscard]] WorldVector3 tangent(double t) const override;
    [[nodiscard]] WorldVector3 normal(double t) const override;
    
    // Accessors for path data
    [[nodiscard]] const std::vector<double>& t_values() const noexcept { return t_values_; }
    [[nodiscard]] const std::vector<double>& u_values() const noexcept { return u_values_; }
    [[nodiscard]] const std::vector<double>& v_values() const noexcept { return v_values_; }
    [[nodiscard]] std::shared_ptr<Surface> surface() const noexcept { return surface_; }

private:
    std::shared_ptr<Surface> surface_;
    std::vector<double> t_values_;
    std::vector<double> u_values_;
    std::vector<double> v_values_;
};

// Path that transitions between multiple surfaces
class TransitionPath final : public SurfacePath {
public:
    void add_segment(
        std::shared_ptr<Surface> surface,
        double t_start, double t_end,
        double u_start, double u_end,
        double v_start, double v_end,
        const WorldVector3& direction
    );
    
    [[nodiscard]] GeometryPoint2 evaluate(double t) const override;
    [[nodiscard]] WorldVector3 tangent(double t) const override;
    [[nodiscard]] WorldVector3 normal(double t) const override;

    // Access segments
    [[nodiscard]] const std::vector<std::unique_ptr<PathSegment>>& segments() const noexcept { 
        return segments_; 
    }

private:
    std::vector<std::unique_ptr<PathSegment>> segments_;
};

} // namespace shap
