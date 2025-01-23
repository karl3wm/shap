#include "shap/coord.hpp"
#pragma once
#include "shap/geometry_point2.hpp"
#include "shap/metric.hpp"
#include "shap/surface.hpp"
#include "shap/validation_config.hpp"
#include <cmath>
#include <stdexcept>

namespace shap {
namespace surfaces {

/**
 * A flat parametric patch - the most fundamental parametric surface.
 * Implements a linear mapping from [0,1]×[0,1] to a planar region in 3D space.
 *
 * Parameter Space Mapping:
 * - Domain: (u,v) ∈ [0,1]×[0,1]
 * - Range: 3D rectangle defined by origin and basis vectors
 * - Formula: P(u,v) = origin + u*world_u + v*world_v
 *
 * Coordinate System:
 * - world_u defines the first coordinate direction in world space
 * - world_v defines the second coordinate direction in world space
 * - These vectors form a basis for the patch's tangent space
 * - Their lengths determine the patch's extent in each direction
 *
 * Properties:
 * - Linear mapping from parameters to world space
 * - Face normal is constant (cross product of basis vectors)
 * - All curvatures are zero (planar surface)
 * - Geodesics are straight lines
 */
class FlatPatch final : public Surface {
public:
    /**
     * Construct a flat parametric patch.
     * @param origin Origin point of the patch
     * @param world_u First basis vector
     * @param world_v Second basis vector
     */
    explicit FlatPatch(
        WorldPoint3 origin,
        WorldVector3 world_u,
        WorldVector3 world_v
    ) : origin_(std::move(origin))
      , world_u_(std::move(world_u))
      , world_v_(std::move(world_v))
      , normal_(0, 0, 0) {
        validate_vectors();  // Check for parallel vectors first
        normal_ = world_u_.cross(world_v_).normalize();
        
        // Setup constant coefficient metric tensor
        const double g11 = world_u_.dot(world_u_);
        const double g12 = world_u_.dot(world_v_);
        const double g22 = world_v_.dot(world_v_);
        metric_tensor_ = std::make_unique<Surface2DMetricTensor>(g11, g12, g12, g22);
    }

    // Move operations
    FlatPatch(FlatPatch&&) noexcept = default;
    FlatPatch& operator=(FlatPatch&&) noexcept = default;

    // Prevent copying
    FlatPatch(const FlatPatch&) = delete;
    FlatPatch& operator=(const FlatPatch&) = delete;

    /**
     * Convert a world space position to parameter space coordinates.
     * @param pos World space position to convert
     * @param vector_length_epsilon Used to handle degenerate cases where basis vectors are nearly parallel
     * @return Parameter space coordinates
     * 
     * Called from:
     * - setup_path_solver() in this file
     */
    [[nodiscard]] ParamPoint3 world_to_param(const WorldPoint3& pos) const override {
        // Solve linear system: pos - origin = u*world_u + v*world_v
        const WorldVector3 rel_pos = pos - origin_;
        
        // Project point onto surface normal to get signed distance
        const double normal_dist = rel_pos.dot(normal_);
        
        // Project point onto surface plane
        const WorldVector3 planar_pos = rel_pos - normal_ * normal_dist;
        
        // Use Cramer's rule for 2x2 system
        const double det = world_u_.cross(world_v_).length();
        if (det < ValidationConfig::instance().vector_length_epsilon()) {
            throw std::invalid_argument(
                "Cannot compute local coordinates: basis vectors are nearly parallel"
            );
        }
        
        // Compute parameter coordinates
        const double u = planar_pos.cross(world_v_).dot(normal_) / det;
        const double v = world_u_.cross(planar_pos).dot(normal_) / det;
        
        return ParamPoint3(u, v, normal_dist);
    }

    /**
     * Evaluate surface at parameter space point.
     * 
     * @param local Parameter space coordinates
     * @return GeometryPoint2 containing full geometric information
     * @throws std::invalid_argument if coordinates are invalid
     */
    [[nodiscard]] GeometryPoint2 evaluate(const ParamPoint2& local) const override {
        // Linear mapping from parameter space to world space
        const WorldPoint3 position = origin_ + 
            world_u_ * local.u() + 
            world_v_ * local.v();
        
        return GeometryPoint2(
            this,
            local,
            position,
            normal_,      // Normal is constant
            world_u_,     // First coordinate basis vector
            world_v_      // Second coordinate basis vector
        );
    }

    [[nodiscard]] std::optional<PathSolver> get_path_solver() const noexcept override {
        return path_solver_;
    }

    [[nodiscard]] SurfaceType surface_type() const noexcept override {
        return SurfaceType::Developable;
    }

    // Access geometry
    [[nodiscard]] const WorldPoint3& origin() const noexcept { return origin_; }
    [[nodiscard]] const WorldVector3& world_u() const noexcept { return world_u_; }
    [[nodiscard]] const WorldVector3& world_v() const noexcept { return world_v_; }
    [[nodiscard]] const WorldVector3& normal() const noexcept { return normal_; }
    [[nodiscard]] const Surface2DMetricTensor& metric_tensor() const noexcept { return *metric_tensor_; }

    /**
     * Setup path solver with given epsilon values.
     * @param vector_length_epsilon Used in world_to_parameter_space_with_epsilon() and for direction projection
     * @param parameter_bound_epsilon Used in check_intersection() for parameter bound checks
     * 
     * Called from:
     * - create_flat_patch() in this file
     */
    void setup_path_solver(double vector_length_epsilon, double parameter_bound_epsilon) noexcept {
        path_solver_ = [this, vector_length_epsilon, parameter_bound_epsilon](
            const WorldPoint3& start, const WorldVector3& dir, double max_t)
            -> std::optional<PathIntersection> {
            
            // Project direction onto face plane
            WorldVector3 planar_dir = dir - dir.dot(normal_) * normal_;
            const double planar_length = planar_dir.length();
            if (planar_length < vector_length_epsilon) {
                return std::nullopt;  // Direction perpendicular to face
            }
            planar_dir = planar_dir * (1.0 / planar_length);
            
            // Convert start point to local coordinates
            const auto start_local = world_to_param(start);
            const auto end_local = world_to_param(start + planar_dir);
            const auto param_dir = std::make_pair(
                end_local.u() - start_local.u(),
                end_local.v() - start_local.v()
            );
            
            // Normalize parameter space direction
            const double param_length = std::sqrt(
                param_dir.first * param_dir.first +
                param_dir.second * param_dir.second
            );
            if (param_length < parameter_bound_epsilon) return std::nullopt;
            
            // Find earliest intersection with parameter bounds
            double min_t = max_t;
            ParamIndex hit_param = ParamIndex::U;
            ParamBound hit_bound = ParamBound::Lower;
            double edge_param = 0.0;
            bool found = false;
            
            // Check all bounds using normalized parameter direction
            found |= check_intersection(
                start_local.u(), param_dir.first / param_length, 0,
                ParamIndex::U, ParamBound::Lower, start_local.v(),
                param_dir, param_length, min_t, hit_param, hit_bound, edge_param,
                parameter_bound_epsilon
            );
            found |= check_intersection(
                start_local.u(), param_dir.first / param_length, 1,
                ParamIndex::U, ParamBound::Upper, start_local.v(),
                param_dir, param_length, min_t, hit_param, hit_bound, edge_param,
                parameter_bound_epsilon
            );
            found |= check_intersection(
                start_local.v(), param_dir.second / param_length, 0,
                ParamIndex::V, ParamBound::Lower, start_local.u(),
                param_dir, param_length, min_t, hit_param, hit_bound, edge_param,
                parameter_bound_epsilon
            );
            found |= check_intersection(
                start_local.v(), param_dir.second / param_length, 1,
                ParamIndex::V, ParamBound::Upper, start_local.u(),
                param_dir, param_length, min_t, hit_param, hit_bound, edge_param,
                parameter_bound_epsilon
            );
            
            if (!found) return std::nullopt;
            
            // Convert parameter space distance to world space
            const double world_t = min_t * (hit_param == ParamIndex::U ? world_u_.length() : world_v_.length());
            
            // Compute intersection position using parameter space mapping
            const double u = hit_param == ParamIndex::U ? 
                static_cast<double>(hit_bound) : start_local.u();
            const double v = hit_param == ParamIndex::V ? 
                static_cast<double>(hit_bound) : start_local.v();
            const WorldPoint3 position = origin_ + world_u_ * u + world_v_ * v;
            
            return PathIntersection(
                world_t,
                position,
                hit_param,
                hit_bound,
                edge_param
            );
        };
    }

private:
    /**
     * Validate basis vectors for degenerate configurations using validation epsilons
     * from ValidationConfig.
     * 
     * Called from:
     * - FlatPatch constructor in this file
     */
    void validate_vectors() {
        const auto& config = ValidationConfig::instance();
        if (world_u_.length_squared() < config.vector_length_epsilon()) {
            throw std::invalid_argument("world_u vector cannot be zero");
        }
        if (world_v_.length_squared() < config.vector_length_epsilon()) {
            throw std::invalid_argument("world_v vector cannot be zero");
        }
        if (std::abs(world_u_.normalize().dot(world_v_.normalize())) > 1.0 - config.vector_parallel_epsilon()) {
            throw std::invalid_argument("world_u and world_v vectors cannot be parallel");
        }
    }

    /**
     * Helper to check intersection with parameter bound.
     * @param curr_param Current parameter value
     * @param d_param Parameter direction
     * @param bound_val Bound value to check against
     * @param param Which parameter (u or v)
     * @param bound Which bound (lower or upper)
     * @param other_param Other parameter value
     * @param param_dir Parameter space direction
     * @param param_length Parameter direction length
     * @param min_t Current minimum intersection time
     * @param hit_param Output: parameter that was hit
     * @param hit_bound Output: bound that was hit
     * @param edge_param Output: parameter value along edge
     * @param parameter_bound_epsilon Threshold for parameter bound checks
     * @return True if intersection found
     * 
     * Called from:
     * - setup_path_solver() in this file
     */
    [[nodiscard]] bool check_intersection(
        double curr_param,
        double d_param,
        double bound_val,
        ParamIndex param,
        ParamBound bound,
        double other_param,
        const std::pair<double, double>& param_dir,
        double param_length,
        double& min_t,
        ParamIndex& hit_param,
        ParamBound& hit_bound,
        double& edge_param,
        double parameter_bound_epsilon
    ) const noexcept {
        if (std::abs(d_param) > parameter_bound_epsilon) {
            const double t = (bound_val - curr_param) / d_param;
            if (t > 0 && t < min_t) {
                // Check if intersection point is within other parameter bounds
                const double other_at_t = other_param + param_dir.second * t / param_length;
                if (other_at_t >= -parameter_bound_epsilon && other_at_t <= 1.0 + parameter_bound_epsilon) {
                    min_t = t;
                    hit_param = param;
                    hit_bound = bound;
                    edge_param = std::clamp(other_at_t, 0.0, 1.0);
                    return true;
                }
            }
        }
        return false;
    }

    WorldPoint3 origin_;
    WorldVector3 world_u_;
    WorldVector3 world_v_;
    WorldVector3 normal_;
    std::unique_ptr<Surface2DMetricTensor> metric_tensor_;
    PathSolver path_solver_;
};

/**
 * Create a flat patch with the given origin and basis vectors.
 * @param origin Origin point of the patch
 * @param world_u First basis vector
 * @param world_v Second basis vector
 * @param vector_length_epsilon Used in world_to_parameter_space_with_epsilon() and setup_path_solver()
 * @param parameter_bound_epsilon Used in setup_path_solver() for parameter bound checks
 * @return Shared pointer to created surface
 * 
 * Called from:
 * - path_length_tests.cpp
 * - space_transformation_tests.cpp
 */
[[nodiscard]] inline std::shared_ptr<Surface> create_flat_patch(
    WorldPoint3 origin,
    WorldVector3 world_u,
    WorldVector3 world_v,
    double vector_length_epsilon,
    double parameter_bound_epsilon
) {
    auto patch = std::make_shared<FlatPatch>(
        std::move(origin),
        std::move(world_u),
        std::move(world_v)
    );
    patch->setup_path_solver(vector_length_epsilon, parameter_bound_epsilon);
    return patch;
}

} // namespace surfaces
} // namespace shap
