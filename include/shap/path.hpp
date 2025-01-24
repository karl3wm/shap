#pragma once
#include "shap/coord.hpp"
#include "shap/geometric_point.hpp"
#include "shap/path3d.hpp"
#include "shap/manifold.hpp"
#include "shap/surface3d.hpp"
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
class SurfacePath : public Path3D {
public:
    using GeometricPoint3D = GeometricPoint<1, 3, WorldSpaceTag>;
    
    SurfacePath() : Path3D(
        [this](double t) { return this->evaluate_position(t); },
        [this](double t) { return this->evaluate_tangent(t); },
        [this](double t) { return this->evaluate_normal(t); }
    ) {}

    virtual ~SurfacePath() = default;
    
    /**
     * Evaluate the path at parameter t.
     */
    [[nodiscard]] GeometricPoint3D evaluate(double t) const {
        validate_parameter(t);
        TargetVector derivs[1] = {evaluate_tangent(t)};
        return GeometricPoint3D(
            this,
            ParamPoint1(t),
            evaluate_position(t),
            derivs
        );
    }

protected:
    [[nodiscard]] virtual WorldPoint3 evaluate_position(double t) const = 0;
    [[nodiscard]] virtual WorldVector3 evaluate_tangent(double t) const = 0;
    [[nodiscard]] virtual WorldVector3 evaluate_normal(double t) const = 0;

    // Validate parameter t is in [0,1]
    static void validate_parameter(double t) {
        if (t < 0.0 || t > 1.0) {
            throw std::invalid_argument("Path parameter t must be in [0,1]");
        }
    }
};

/**
 * Geodesic curve between two points on a surface.
 */
class GeodesicCurve final : public SurfacePath {
public:
    GeodesicCurve(
        std::shared_ptr<Surface3D> surface,
        const GeometricPoint<2, 3, WorldSpaceTag>& start,
        const GeometricPoint<2, 3, WorldSpaceTag>& end
    );
    
protected:
    [[nodiscard]] WorldPoint3 evaluate_position(double t) const override;
    [[nodiscard]] WorldVector3 evaluate_tangent(double t) const override;
    [[nodiscard]] WorldVector3 evaluate_normal(double t) const override;

private:
    void compute_smooth_geodesic(
        const GeometricPoint<2, 3, WorldSpaceTag>& start,
        const GeometricPoint<2, 3, WorldSpaceTag>& end
    );
    
    void compute_developable_geodesic(
        const GeometricPoint<2, 3, WorldSpaceTag>& start,
        const GeometricPoint<2, 3, WorldSpaceTag>& end
    );
    
    std::shared_ptr<Surface3D> surface_;
    std::vector<GeometricPoint<2, 3, WorldSpaceTag>> points_;
};

/**
 * Path segment on a single surface.
 */
class PathSegment final : public SurfacePath {
public:
    explicit PathSegment(std::shared_ptr<Surface3D> surface) noexcept
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
    
public:
    [[nodiscard]] WorldPoint3 evaluate_position(double t) const override;
    [[nodiscard]] WorldVector3 evaluate_tangent(double t) const override;
    [[nodiscard]] WorldVector3 evaluate_normal(double t) const override;
    
    // Accessors for path data
    [[nodiscard]] const std::vector<double>& t_values() const noexcept { return t_values_; }
    [[nodiscard]] const std::vector<double>& u_values() const noexcept { return u_values_; }
    [[nodiscard]] const std::vector<double>& v_values() const noexcept { return v_values_; }
    [[nodiscard]] std::shared_ptr<Surface3D> surface() const noexcept { return surface_; }

private:
    std::shared_ptr<Surface3D> surface_;
    std::vector<double> t_values_;
    std::vector<double> u_values_;
    std::vector<double> v_values_;
};

/**
 * Path that transitions between multiple surfaces.
 */
class TransitionPath final : public SurfacePath {
public:
    void add_segment(
        std::shared_ptr<Surface3D> surface,
        double t_start, double t_end,
        double u_start, double u_end,
        double v_start, double v_end,
        const WorldVector3& direction
    );
    
protected:
    [[nodiscard]] WorldPoint3 evaluate_position(double t) const override;
    [[nodiscard]] WorldVector3 evaluate_tangent(double t) const override;
    [[nodiscard]] WorldVector3 evaluate_normal(double t) const override;

public:
    // Access segments
    [[nodiscard]] const std::vector<std::unique_ptr<PathSegment>>& segments() const noexcept { 
        return segments_; 
    }

private:
    std::vector<std::unique_ptr<PathSegment>> segments_;
};

} // namespace shap
