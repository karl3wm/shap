#pragma once
#include "coord.hpp"
#include "surface_type.hpp"
#include "path_intersection.hpp"
#include "manifold.hpp"
#include "geometric_point.hpp"
#include <functional>
#include <memory>
#include <optional>
#include <utility>
#include <stdexcept>

namespace shap {

// Forward declarations
class Path3D;
class SurfacePath;

// Function types for surface creation and metric calculations
using PositionFunction = std::function<WorldPoint3(const ParamPoint2&)>;
using DerivativeFunction = std::function<WorldVector3(const ParamPoint2&)>;
using CurvatureFunction = std::function<double(const ParamPoint2&)>;
using ParameterSpaceDerivative = std::function<double(const ParamPoint2&)>;

using PathSolver = std::function<std::optional<PathIntersection>(
    const WorldPoint3& world_start,
    const WorldVector3& world_direction,
    double max_world_distance
)>;

// Forward declarations
class RiemannianMetric;

/**
 * Represents a 2D surface embedded in 3D world space.
 * Implements the Manifold interface for 2D->3D mappings.
 */
class Surface3DImpl final : public Manifold<2, 3, WorldSpaceTag>,
                          public std::enable_shared_from_this<Surface3DImpl> {
    friend class RiemannianMetric;
public:
    // Constructor with all required function objects
    Surface3DImpl(
        PositionFunction position_func,
        DerivativeFunction du_func,
        DerivativeFunction dv_func,
        DerivativeFunction duu_func,
        DerivativeFunction duv_func,
        DerivativeFunction dvv_func,
        CurvatureFunction gaussian_func,
        CurvatureFunction mean_func,
        std::optional<PathSolver> path_solver = std::nullopt,
        SurfaceType type = SurfaceType::Smooth,
        ParameterSpaceDerivative du2_du = nullptr,
        ParameterSpaceDerivative du2_dv = nullptr,
        ParameterSpaceDerivative duv_du = nullptr,
        ParameterSpaceDerivative duv_dv = nullptr,
        ParameterSpaceDerivative dv2_du = nullptr,
        ParameterSpaceDerivative dv2_dv = nullptr
    ) : Manifold<2, 3, WorldSpaceTag>(
            [f = position_func](const ParameterPoint& param) {
                return f(param);
            },
            [du = du_func, dv = dv_func](const ParameterPoint& param, int derivative_index) {
                return derivative_index == 0 ? du(param) : dv(param);
            }
        )
      , position_func_(std::move(position_func))
      , du_func_(std::move(du_func))
      , dv_func_(std::move(dv_func))
      , duu_func_(std::move(duu_func))
      , duv_func_(std::move(duv_func))
      , dvv_func_(std::move(dvv_func))
      , gaussian_curv_func_(std::move(gaussian_func))
      , mean_curv_func_(std::move(mean_func))
      , path_solver_(std::move(path_solver))
      , type_(type)
      , du2_du_fn_(std::move(du2_du))
      , du2_dv_fn_(std::move(du2_dv))
      , duv_du_fn_(std::move(duv_du))
      , duv_dv_fn_(std::move(duv_dv))
      , dv2_du_fn_(std::move(dv2_du))
      , dv2_dv_fn_(std::move(dv2_dv))
    {
        if (!position_func_ || !du_func_ || !dv_func_ || !duu_func_ || !duv_func_ || !dvv_func_ ||
            !gaussian_curv_func_ || !mean_curv_func_) {
            throw std::invalid_argument("Required surface functions cannot be null");
        }
    }
    
    // Prevent copying
    Surface3DImpl(const Surface3DImpl&) = delete;
    Surface3DImpl& operator=(const Surface3DImpl&) = delete;
    
    // Allow moving
    Surface3DImpl(Surface3DImpl&&) noexcept = default;
    Surface3DImpl& operator=(Surface3DImpl&&) noexcept = default;

    // Implement Manifold interface
    using Manifold<2, 3, WorldSpaceTag>::evaluate;  // Use base class implementation
    
    [[nodiscard]] std::array<TargetVector, 2>
    derivatives(const ParameterPoint& param) const override {
        return {du_func_(param), dv_func_(param)};
    }

    // Surface-specific functionality
    [[nodiscard]] std::pair<double, double> get_scale_factors(
        const ParameterPoint& local
    ) const;

    [[nodiscard]] SurfaceType surface_type() const noexcept {
        return type_;
    }

    /**
     * Convert a world space position to local coordinates.
     */
    [[nodiscard]] ParamPoint3 world_to_param(const WorldPoint3& pos) const;
    
    /**
     * Create a path on the surface starting from a point in a given direction.
     */
    [[nodiscard]] std::unique_ptr<SurfacePath> create_path(
        const GeometricPoint<2, 3, WorldSpaceTag>& start,
        const WorldVector3& world_direction,
        double world_length
    ) const;

    // Metric component derivative accessors
    [[nodiscard]] double du2_du(const ParameterPoint& param) const noexcept {
        return du2_du_fn_ ? du2_du_fn_(param) : 0.0;
    }
    [[nodiscard]] double du2_dv(const ParameterPoint& param) const noexcept {
        return du2_dv_fn_ ? du2_dv_fn_(param) : 0.0;
    }
    [[nodiscard]] double duv_du(const ParameterPoint& param) const noexcept {
        return duv_du_fn_ ? duv_du_fn_(param) : 0.0;
    }
    [[nodiscard]] double duv_dv(const ParameterPoint& param) const noexcept {
        return duv_dv_fn_ ? duv_dv_fn_(param) : 0.0;
    }
    [[nodiscard]] double dv2_du(const ParameterPoint& param) const noexcept {
        return dv2_du_fn_ ? dv2_du_fn_(param) : 0.0;
    }
    [[nodiscard]] double dv2_dv(const ParameterPoint& param) const noexcept {
        return dv2_dv_fn_ ? dv2_dv_fn_(param) : 0.0;
    }

    /**
     * Get path solver if available.
     */
    [[nodiscard]] std::optional<PathSolver> get_path_solver() const noexcept {
        return path_solver_;
    }

private:
    // Surface functions
    PositionFunction position_func_;
    DerivativeFunction du_func_;
    DerivativeFunction dv_func_;
    DerivativeFunction duu_func_;
    DerivativeFunction duv_func_;
    DerivativeFunction dvv_func_;
    CurvatureFunction gaussian_curv_func_;
    CurvatureFunction mean_curv_func_;
    std::optional<PathSolver> path_solver_;
    SurfaceType type_;

    // Parameter space derivative functions
    ParameterSpaceDerivative du2_du_fn_;  // d(du·du)/du
    ParameterSpaceDerivative du2_dv_fn_;  // d(du·du)/dv
    ParameterSpaceDerivative duv_du_fn_;  // d(du·dv)/du
    ParameterSpaceDerivative duv_dv_fn_;  // d(du·dv)/dv
    ParameterSpaceDerivative dv2_du_fn_;  // d(dv·dv)/du
    ParameterSpaceDerivative dv2_dv_fn_;  // d(dv·dv)/dv
};

// Type alias for backward compatibility and convenience
using Surface3D = Surface3DImpl;

} // namespace shap
