#include "coord.hpp"
#pragma once
#include "edge_descriptor.hpp"
#include <optional>
#include <utility>

namespace shap {

/**
 * Represents a point on a surface with complete geometric information.
 * Combines local parameter space coordinates with world space geometric properties
 * including position, derivatives, and curvature information.
 */
class Surface;  // Forward declaration

class GeometryPoint2 {
public:
    /**
     * Construct with basic geometric properties.
     */
    GeometryPoint2(
        const Surface* surface,
        ParamPoint2 local,
        WorldPoint3 world_pos,
        WorldVector3 world_normal,
        WorldVector3 world_du,
        WorldVector3 world_dv
    ) : surface_(surface)
      , local_pos_(std::move(local))
      , world_pos_(std::move(world_pos))
      , world_normal_(std::move(world_normal).normalized())
      , world_du_(std::move(world_du))
      , world_dv_(std::move(world_dv)) {}

    /**
     * Construct with complete geometric properties including derivatives and curvature.
     */
    GeometryPoint2(
        const Surface* surface,
        ParamPoint2 local,
        WorldPoint3 world_pos,
        WorldVector3 world_normal,
        WorldVector3 world_du,
        WorldVector3 world_dv,
        WorldVector3 world_duu,
        WorldVector3 world_duv,
        WorldVector3 world_dvv,
        double gaussian_curvature,
        double mean_curvature,
        std::pair<double, double> principal_curvatures
    ) : surface_(surface)
      , local_pos_(std::move(local))
      , world_pos_(std::move(world_pos))
      , world_normal_(std::move(world_normal).normalized())
      , world_du_(std::move(world_du))
      , world_dv_(std::move(world_dv))
      , world_duu_(std::make_optional(std::move(world_duu)))
      , world_duv_(std::make_optional(std::move(world_duv)))
      , world_dvv_(std::make_optional(std::move(world_dvv)))
      , gaussian_curvature_(std::make_optional(gaussian_curvature))
      , mean_curvature_(std::make_optional(mean_curvature))
      , principal_curvatures_(std::make_optional(std::move(principal_curvatures))) {}

    // Param parameter space access
    [[nodiscard]] const ParamPoint2& local_pos() const noexcept { 
        return local_pos_; 
    }
    
    [[nodiscard]] double u() const noexcept { return local_pos_.u(); }
    [[nodiscard]] double v() const noexcept { return local_pos_.v(); }
    
    // World space access - first order properties
    [[nodiscard]] const WorldPoint3& world_pos() const noexcept { return world_pos_; }
    [[nodiscard]] const WorldVector3& world_normal() const noexcept { return world_normal_; }
    [[nodiscard]] const WorldVector3& world_du() const noexcept { return world_du_; }
    [[nodiscard]] const WorldVector3& world_dv() const noexcept { return world_dv_; }

    // World space access - second order properties
    [[nodiscard]] const std::optional<WorldVector3>& world_duu() const noexcept { return world_duu_; }
    [[nodiscard]] const std::optional<WorldVector3>& world_duv() const noexcept { return world_duv_; }
    [[nodiscard]] const std::optional<WorldVector3>& world_dvv() const noexcept { return world_dvv_; }

    // Curvature information
    [[nodiscard]] const std::optional<double>& gaussian_curvature() const noexcept { return gaussian_curvature_; }
    [[nodiscard]] const std::optional<double>& mean_curvature() const noexcept { return mean_curvature_; }
    [[nodiscard]] const std::optional<std::pair<double, double>>& principal_curvatures() const noexcept { 
        return principal_curvatures_; 
    }

    // Scale factors for space conversion
    [[nodiscard]] std::pair<double, double> get_scale_factors() const {
        return {world_du_.length(), world_dv_.length()};
    }

    // Edge information
    [[nodiscard]] bool is_on_edge() const noexcept {
        return u() == 0.0 || u() == 1.0 || v() == 0.0 || v() == 1.0;
    }

    [[nodiscard]] std::optional<EdgeDescriptor> get_edge_descriptor() const noexcept {
        if (!is_on_edge()) return std::nullopt;

        if (u() == 0.0) return EdgeDescriptor{ParamIndex::U, ParamBound::Lower, v()};
        if (u() == 1.0) return EdgeDescriptor{ParamIndex::U, ParamBound::Upper, v()};
        if (v() == 0.0) return EdgeDescriptor{ParamIndex::V, ParamBound::Lower, u()};
        if (v() == 1.0) return EdgeDescriptor{ParamIndex::V, ParamBound::Upper, u()};
        return std::nullopt;
    }

    // Surface access
    [[nodiscard]] const Surface* surface() const noexcept { return surface_; }

private:
    // Parent surface
    const Surface* surface_;

    // First order properties (always present)
    ParamPoint2 local_pos_;       // Parameter space coordinates
    WorldPoint3 world_pos_;       // Position in world space
    WorldVector3 world_normal_;   // Unit surface normal in world space
    WorldVector3 world_du_;       // First derivative in u direction
    WorldVector3 world_dv_;       // First derivative in v direction

    // Second order properties (optional)
    std::optional<WorldVector3> world_duu_;  // Second derivative in u direction
    std::optional<WorldVector3> world_duv_;  // Mixed second derivative
    std::optional<WorldVector3> world_dvv_;  // Second derivative in v direction

    // Curvature properties (optional)
    std::optional<double> gaussian_curvature_;
    std::optional<double> mean_curvature_;
    std::optional<std::pair<double, double>> principal_curvatures_;
};

} // namespace shap
