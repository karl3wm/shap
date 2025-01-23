#include "coord.hpp"
#pragma once
#include "edge_connection.hpp"
#include "edge_descriptor.hpp"
#include "geometry_point2.hpp"
#include "param_bound.hpp"
#include "param_index.hpp"
#include "surface_type.hpp"
#include <functional>
#include <memory>
#include <optional>
#include <utility>
#include <stdexcept>

namespace shap {

class SurfacePoint;
class SurfacePath;

// Function types for surface creation
using PositionFunction = std::function<WorldPoint3(const ParamPoint2&)>;
using DerivativeFunction = std::function<WorldVector3(const ParamPoint2&)>;
using CurvatureFunction = std::function<double(const ParamPoint2&)>;

// Path solver returns intersection with surface boundary
struct PathIntersection {
    double t;                // Distance to intersection in world space
    WorldPoint3 position;    // World space intersection point
    ParamIndex param;        // Which parameter (u/v) hit boundary
    ParamBound bound;        // Which bound (0/1) was hit
    double edge_parameter;   // Parameter along the edge [0,1]

    // Constructor with validation
    PathIntersection(
        double t_,
        WorldPoint3 position_,
        ParamIndex param_,
        ParamBound bound_,
        double edge_parameter_
    ) : t(t_)
      , position(std::move(position_))
      , param(param_)
      , bound(bound_)
      , edge_parameter(edge_parameter_) {
        if (t_ < 0) {
            throw std::invalid_argument("Intersection distance must be non-negative");
        }
        if (edge_parameter_ < 0 || edge_parameter_ > 1) {
            throw std::invalid_argument("Edge parameter must be in [0,1]");
        }
    }
};

using PathSolver = std::function<std::optional<PathIntersection>(
    const WorldPoint3& world_start,
    const WorldVector3& world_direction,
    double max_world_distance
)>;

class Surface {
public:
    virtual ~Surface() = default;
    
    // Prevent copying
    Surface(const Surface&) = delete;
    Surface& operator=(const Surface&) = delete;
    
    // Allow moving
    Surface(Surface&&) noexcept = default;
    Surface& operator=(Surface&&) noexcept = default;

protected:
    Surface() = default;

public:
    /**
     * Evaluate surface at parameter space point.
     * 
     * @param local Parameter space coordinates
     * @return GeometryPoint2 containing full geometric information
     * @throws std::invalid_argument if coordinates are invalid
     */
    [[nodiscard]] virtual GeometryPoint2 evaluate(const ParamPoint2& local) const = 0;
    
    /**
     * Convert a world space position to local coordinates.
     * 
     * This function computes three coordinates that fully describe a point's position
     * relative to the surface:
     * - u,v: Param parameter coordinates in [0,1]×[0,1]
     * - normal: Signed distance along surface normal vector
     *
     * For points on the surface, normal will be 0 (within ValidationConfig::vector_length_epsilon).
     * Positive normal indicates the point is on the positive side of the surface
     * (in the direction of the normal vector).
     *
     * @param pos World space position to convert
     * @return ParamPoint3 containing local coordinates
     * @throws std::invalid_argument if coordinate computation fails
     */
    [[nodiscard]] virtual ParamPoint3 world_to_param(const WorldPoint3& pos) const = 0;
    
    /**
     * Convert a world space position to surface parameter coordinates.
     * Projects the point onto the surface along the normal direction.
     *
     * @param pos World space position to convert
     * @return ParamPoint2 containing parameter coordinates
     * @throws std::invalid_argument if coordinate computation fails
     */
    [[nodiscard]] virtual ParamPoint2 world_to_param_r2(const WorldPoint3& pos) const {
        return world_to_param(pos).to_r2();
    }
    
    /**
     * Create a path on the surface starting from a point in a given direction.
     * 
     * @param start Starting point on the surface
     * @param world_direction Desired world-space direction (will be projected onto surface)
     * @param world_length Desired path length in world space units
     * @throws std::invalid_argument if preconditions are not met
     * @return Unique pointer to path object representing the curve
     */
    [[nodiscard]] virtual std::unique_ptr<SurfacePath> create_path(
        const GeometryPoint2& start,
        const WorldVector3& world_direction,
        double world_length
    ) const;
    
    // Get path solver if available
    [[nodiscard]] virtual std::optional<PathSolver> get_path_solver() const noexcept {
        return std::nullopt;
    }
    
    // Get surface type
    [[nodiscard]] virtual SurfaceType surface_type() const noexcept {
        return SurfaceType::Generic;
    }
    
    
    /**
     * Convert world space direction to parameter space velocity.
     * Accounts for surface metric tensor in the conversion.
     *
     * @param world_direction Direction vector in world space
     * @param world_du First derivative in u direction
     * @param world_dv First derivative in v direction
     * @return Velocity vector in parameter space
     */
    [[nodiscard]] WorldVector3 world_to_parameter_velocity(
        const WorldVector3& world_direction,
        const WorldVector3& world_du,
        const WorldVector3& world_dv
    ) const noexcept;
    
    /**
     * Get scale factors for converting between parameter and world space.
     * These represent how much a unit step in parameter space maps to in world space.
     *
     * @param param Parameter space point to compute scale factors at
     * @return Pair of scale factors (du_scale, dv_scale)
     */
    [[nodiscard]] std::pair<double, double> get_scale_factors(
        const ParamPoint2& local
    ) const {
        const auto geom = evaluate(local);
        return {geom.world_du().length(), geom.world_dv().length()};
    }

protected:
    // Validate parameter values are in [0,1]
    static void validate_parameters([[maybe_unused]] const ParamPoint2& local) {
        // ParamPoint2 constructor handles validation
    }

public:
    // Factory methods
    // Factory methods with updated parameter types
    [[nodiscard]] static std::shared_ptr<Surface> create(
        PositionFunction position_func,
        std::optional<PathSolver> path_solver = std::nullopt,
        SurfaceType type = SurfaceType::Generic
    );
};

} // namespace shap
