#pragma once
#include "shap/coord.hpp"
#include "shap/geometric_point.hpp"
#include "shap/surface3d.hpp"
#include "shap/validation_config.hpp"
#include <cmath>
#include <stdexcept>
#include <functional>
#include <memory>

namespace shap {
namespace surfaces {

namespace {
/**
 * Helper to validate basis vectors for degenerate configurations using validation epsilons
 * from ValidationConfig.
 */
void validate_vectors(const WorldVector3& world_u, const WorldVector3& world_v) {
    const auto& config = ValidationConfig::instance();
    if (world_u.length_squared() < config.vector_length_epsilon()) {
        throw std::invalid_argument("world_u vector cannot be zero");
    }
    if (world_v.length_squared() < config.vector_length_epsilon()) {
        throw std::invalid_argument("world_v vector cannot be zero");
    }
    if (std::abs(world_u.normalized().dot(world_v.normalized())) > 1.0 - config.vector_parallel_epsilon()) {
        throw std::invalid_argument("world_u and world_v vectors cannot be parallel");
    }
}
} // anonymous namespace

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
class FlatPatch final : public Surface3D {
public:
    FlatPatch(WorldPoint3 origin_, WorldVector3 world_u_, WorldVector3 world_v_,
              double vector_length_epsilon = ValidationConfig::instance().vector_length_epsilon(),
              double parameter_bound_epsilon = ValidationConfig::instance().parameter_bound_epsilon())
        : Surface3D(
            std::bind(&FlatPatch::position, this, std::placeholders::_1),
            std::bind(&FlatPatch::du, this, std::placeholders::_1),
            std::bind(&FlatPatch::dv, this, std::placeholders::_1),
            std::bind(&FlatPatch::duu, this, std::placeholders::_1),
            std::bind(&FlatPatch::duv, this, std::placeholders::_1),
            std::bind(&FlatPatch::dvv, this, std::placeholders::_1),
            std::bind(&FlatPatch::gaussian, this, std::placeholders::_1),
            std::bind(&FlatPatch::mean, this, std::placeholders::_1),
            std::bind(&FlatPatch::nearest, this, std::placeholders::_1),
            [this, vector_length_epsilon, parameter_bound_epsilon](const WorldPoint3& start, const WorldVector3& dir, double max_t) {
                return solve_path(start, dir, max_t, vector_length_epsilon, parameter_bound_epsilon);
            },
            SurfaceType::Developable,
            std::bind(&FlatPatch::du2_du, this, std::placeholders::_1),
            std::bind(&FlatPatch::du2_dv, this, std::placeholders::_1),
            std::bind(&FlatPatch::duv_du, this, std::placeholders::_1),
            std::bind(&FlatPatch::duv_dv, this, std::placeholders::_1),
            std::bind(&FlatPatch::dv2_du, this, std::placeholders::_1),
            std::bind(&FlatPatch::dv2_dv, this, std::placeholders::_1)
        )
        , origin(std::move(origin_))
        , world_u(std::move(world_u_))
        , world_v(std::move(world_v_))
    {
        validate_vectors(world_u, world_v);
        normal = world_u.crossed(world_v).normalized();
    }

    WorldPoint3 position(const ParamPoint2& local) const {
        return origin + world_u * local.u() + world_v * local.v();
    }

    WorldVector3 du(const ParamPoint2&) const { return world_u; }
    WorldVector3 dv(const ParamPoint2&) const { return world_v; }
    WorldVector3 duu(const ParamPoint2&) const { return WorldVector3(0, 0, 0); }
    WorldVector3 duv(const ParamPoint2&) const { return WorldVector3(0, 0, 0); }
    WorldVector3 dvv(const ParamPoint2&) const { return WorldVector3(0, 0, 0); }
    double gaussian(const ParamPoint2&) const { return 0.0; }
    double mean(const ParamPoint2&) const { return 0.0; }

    std::optional<PathIntersection> solve_path(
        const WorldPoint3& start, const WorldVector3& dir, double max_t,
        double vector_length_epsilon, double parameter_bound_epsilon) const {
        
        // Project direction onto face plane
        WorldVector3 planar_dir = dir - dir.dot(normal) * normal;
        const double planar_length = planar_dir.length();
        if (planar_length < vector_length_epsilon) {
            return std::nullopt;  // Direction perpendicular to face
        }
        planar_dir = planar_dir * (1.0 / planar_length);
        
        // Convert start point to local coordinates
        const auto start_local = nearest(start);
        const auto end_local = nearest(start + planar_dir);
        const auto param_dir = end_local - start_local;
        
        // Normalize parameter space direction
        const double param_length = param_dir.length();
        if (param_length < parameter_bound_epsilon) return std::nullopt;
        
        // Find earliest intersection with parameter bounds
        double min_t = max_t;
        ParamIndex hit_param = ParamIndex::U;
        ParamBound hit_bound = ParamBound::Lower;
        double edge_param = 0.0;
        bool found = false;
        
        // Check U parameter bounds
        const double d_u = param_dir.u() / param_length;
        if (std::abs(d_u) > parameter_bound_epsilon) {
            // Check U=0 bound
            double t = -start_local.u() / d_u;
            if (t > 0 && t < min_t) {
                const double v_at_t = start_local.v() + param_dir.v() * t / param_length;
                if (v_at_t >= -parameter_bound_epsilon && v_at_t <= 1.0 + parameter_bound_epsilon) {
                    min_t = t;
                    hit_param = ParamIndex::U;
                    hit_bound = ParamBound::Lower;
                    edge_param = std::clamp(v_at_t, 0.0, 1.0);
                    found = true;
                }
            }
            // Check U=1 bound
            t = (1.0 - start_local.u()) / d_u;
            if (t > 0 && t < min_t) {
                const double v_at_t = start_local.v() + param_dir.v() * t / param_length;
                if (v_at_t >= -parameter_bound_epsilon && v_at_t <= 1.0 + parameter_bound_epsilon) {
                    min_t = t;
                    hit_param = ParamIndex::U;
                    hit_bound = ParamBound::Upper;
                    edge_param = std::clamp(v_at_t, 0.0, 1.0);
                    found = true;
                }
            }
        }

        // Check V parameter bounds
        const double d_v = param_dir.v() / param_length;
        if (std::abs(d_v) > parameter_bound_epsilon) {
            // Check V=0 bound
            double t = -start_local.v() / d_v;
            if (t > 0 && t < min_t) {
                const double u_at_t = start_local.u() + param_dir.u() * t / param_length;
                if (u_at_t >= -parameter_bound_epsilon && u_at_t <= 1.0 + parameter_bound_epsilon) {
                    min_t = t;
                    hit_param = ParamIndex::V;
                    hit_bound = ParamBound::Lower;
                    edge_param = std::clamp(u_at_t, 0.0, 1.0);
                    found = true;
                }
            }
            // Check V=1 bound
            t = (1.0 - start_local.v()) / d_v;
            if (t > 0 && t < min_t) {
                const double u_at_t = start_local.u() + param_dir.u() * t / param_length;
                if (u_at_t >= -parameter_bound_epsilon && u_at_t <= 1.0 + parameter_bound_epsilon) {
                    min_t = t;
                    hit_param = ParamIndex::V;
                    hit_bound = ParamBound::Upper;
                    edge_param = std::clamp(u_at_t, 0.0, 1.0);
                    found = true;
                }
            }
        }
        
        if (!found) return std::nullopt;
        
        // Convert parameter space distance to world space
        const double world_t = min_t * (hit_param == ParamIndex::U ? world_u.length() : world_v.length());
        
        // Compute intersection position using parameter space mapping
        const double u = hit_param == ParamIndex::U ? 
            static_cast<double>(hit_bound) : start_local.u();
        const double v = hit_param == ParamIndex::V ? 
            static_cast<double>(hit_bound) : start_local.v();
        const WorldPoint3 position = origin + world_u * u + world_v * v;
        
        return PathIntersection(
            world_t,
            position,
            hit_param,
            hit_bound,
            edge_param
        );
    }

    ParamPoint2 nearest(const WorldPoint3& pos) const {
        // Solve linear system: pos - origin = u*world_u + v*world_v
        const WorldVector3 rel_pos = pos - origin;
        
        // Project point onto surface plane
        const WorldVector3 planar_pos = rel_pos - rel_pos.dot(normal) * normal;
        
        // Use Cramer's rule for 2x2 system
        const double det = world_u.crossed(world_v).length();
        if (det < ValidationConfig::instance().vector_length_epsilon()) {
            throw std::invalid_argument(
                "Cannot compute local coordinates: basis vectors are nearly parallel"
            );
        }
        
        // Compute parameter coordinates
        const double u = planar_pos.crossed(world_v).dot(normal) / det;
        const double v = world_u.crossed(planar_pos).dot(normal) / det;
        
        return ParamPoint2(u, v);
    }

    double du2_du(const ParamPoint2&) const { return 0.0; }
    double du2_dv(const ParamPoint2&) const { return 0.0; }
    double duv_du(const ParamPoint2&) const { return 0.0; }
    double duv_dv(const ParamPoint2&) const { return 0.0; }
    double dv2_du(const ParamPoint2&) const { return 0.0; }
    double dv2_dv(const ParamPoint2&) const { return 0.0; }

private:
    WorldPoint3 origin;
    WorldVector3 world_u;
    WorldVector3 world_v;
    WorldVector3 normal;
};

} // namespace surfaces
} // namespace shap
