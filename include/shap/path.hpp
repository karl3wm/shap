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
 * Base class for paths on surfaces.
 *
 * A path represents a curve that lies on one or more surfaces. The curve
 * is parameterized by t ∈ [0,1], where:
 * - t=0 corresponds to the start point
 * - t=1 corresponds to the end point (or last point before surface boundary)
 */
class SurfacePath : public Path3D {
public:
    using GeometricPoint3D = GeometricPoint<1, 3, WorldSpaceTag>;
    using PositionFunc = std::function<WorldPoint3(ParamPoint1)>;
    using TangentFunc = std::function<WorldVector3(ParamPoint1)>;
    using NormalFunc = std::function<WorldVector3(ParamPoint1)>;
    
    SurfacePath(PositionFunc position, TangentFunc tangent, NormalFunc normal) 
        : Path3D(
            std::bind(std::move(position), std::placeholders::_1),
            std::bind(std::move(tangent), std::placeholders::_1),
            std::bind(std::move(normal), std::placeholders::_1)
        ) {}

    /**
     * Evaluate the path at parameter t.
     */
    [[nodiscard]] GeometricPoint3D evaluate(const ParamPoint1& param) const {
        TargetVector derivs[1] = {Path3D::derivatives(param)[0]};
        return GeometricPoint3D(
            this,
            param,
            Path3D::evaluate(param).world_pos(),
            derivs
        );
    }

protected:
};

/**
 * Geodesic curve between two points on a surface.
 */
class GeodesicCurve final : public Path3D {
public:
    GeodesicCurve(
        std::shared_ptr<const Surface3D> surface,
        const GeometricPoint<2, 3, WorldSpaceTag>& start,
        const GeometricPoint<2, 3, WorldSpaceTag>& end
    );

private:
    void compute_smooth_geodesic(
        const GeometricPoint<2, 3, WorldSpaceTag>& start,
        const GeometricPoint<2, 3, WorldSpaceTag>& end
    );
    
    void compute_developable_geodesic(
        const GeometricPoint<2, 3, WorldSpaceTag>& start,
        const GeometricPoint<2, 3, WorldSpaceTag>& end
    );
    
    // Path evaluation functions
    [[nodiscard]] WorldPoint3 evaluate_position(const ParamPoint1& param) const;
    [[nodiscard]] WorldVector3 evaluate_tangent(const ParamPoint1& param) const;
    [[nodiscard]] WorldVector3 evaluate_normal(const ParamPoint1& param) const;

    std::shared_ptr<const Surface3D> surface_;
    std::vector<GeometricPoint<2, 3, WorldSpaceTag>> points_;
};

/**
 * Path segment on a single surface.
 */
class PathSegment final : public Path3D {
public:
    explicit PathSegment(std::shared_ptr<const Surface3D> surface);

    void add_point(const ParamPoint1& t, const ParamPoint2& uv);
    
    // Accessors for path data
    [[nodiscard]] const std::vector<ParamPoint1>& t_values() const noexcept { return t_values_; }
    [[nodiscard]] const std::vector<ParamPoint2>& uv_values() const noexcept { return uv_values_; }
    [[nodiscard]] std::shared_ptr<const Surface3D> surface() const noexcept { return surface_; }

    // Path evaluation functions
    [[nodiscard]] WorldPoint3 evaluate_position(const ParamPoint1& param) const;
    [[nodiscard]] WorldVector3 evaluate_tangent(const ParamPoint1& param) const;
    [[nodiscard]] WorldVector3 evaluate_normal(const ParamPoint1& param) const;

private:
    std::shared_ptr<const Surface3D> surface_;
    std::vector<ParamPoint1> t_values_;
    std::vector<ParamPoint2> uv_values_;
};

/**
 * Path that transitions between multiple surfaces.
 */
class TransitionPath final : public Path3D {
public:
    TransitionPath();

    void add_segment(
        std::shared_ptr<const Surface3D> surface,
        const ParamPoint1& t_start, const ParamPoint1& t_end,
        const ParamPoint2& uv_start, const ParamPoint2& uv_end,
        const WorldVector3& direction
    );

    // Access segments
    [[nodiscard]] const std::vector<std::unique_ptr<PathSegment>>& segments() const noexcept { 
        return segments_; 
    }

private:
    // Path evaluation functions
    [[nodiscard]] WorldPoint3 evaluate_position(const ParamPoint1& param) const;
    [[nodiscard]] WorldVector3 evaluate_tangent(const ParamPoint1& param) const;
    [[nodiscard]] WorldVector3 evaluate_normal(const ParamPoint1& param) const;

    std::vector<std::unique_ptr<PathSegment>> segments_;
};

} // namespace shap
