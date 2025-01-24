#pragma once
#include "shap/coord.hpp"
#include "shap/geometric_point.hpp"
#include "shap/surface3d.hpp"
#include "shap/validation_config.hpp"
#include <cmath>
#include <stdexcept>
#include <functional>

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
 * Create a flat parametric patch.
 * @param origin Origin point of the patch
 * @param world_u First basis vector
 * @param world_v Second basis vector
 * @param vector_length_epsilon Used in world_to_parameter_space_with_epsilon() and setup_path_solver()
 * @param parameter_bound_epsilon Used in setup_path_solver() for parameter bound checks
 * @return Shared pointer to created surface
 */
[[nodiscard]] std::shared_ptr<Surface3D> create_flat_patch(
    WorldPoint3 origin,
    WorldVector3 world_u,
    WorldVector3 world_v,
    double vector_length_epsilon = ValidationConfig::instance().vector_length_epsilon(),
    double parameter_bound_epsilon = ValidationConfig::instance().parameter_bound_epsilon()
) {
    // Create derived class instance to bind member functions from
    class FlatPatchImpl {
    public:
        FlatPatchImpl(WorldPoint3 origin_, WorldVector3 world_u_, WorldVector3 world_v_)
            : origin(std::move(origin_))
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
            
            // Check U parameter bounds
            const double d_u = param_dir.first / param_length;
            if (std::abs(d_u) > parameter_bound_epsilon) {
                // Check U=0 bound
                double t = -start_local.u() / d_u;
                if (t > 0 && t < min_t) {
                    const double v_at_t = start_local.v() + param_dir.second * t / param_length;
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
                    const double v_at_t = start_local.v() + param_dir.second * t / param_length;
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
            const double d_v = param_dir.second / param_length;
            if (std::abs(d_v) > parameter_bound_epsilon) {
                // Check V=0 bound
                double t = -start_local.v() / d_v;
                if (t > 0 && t < min_t) {
                    const double u_at_t = start_local.u() + param_dir.first * t / param_length;
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
                    const double u_at_t = start_local.u() + param_dir.first * t / param_length;
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

        ParamPoint3 world_to_param(const WorldPoint3& pos) const {
            // Solve linear system: pos - origin = u*world_u + v*world_v
            const WorldVector3 rel_pos = pos - origin;
            
            // Project point onto surface normal to get signed distance
            const double normal_dist = rel_pos.dot(normal);
            
            // Project point onto surface plane
            const WorldVector3 planar_pos = rel_pos - normal * normal_dist;
            
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
            
            return ParamPoint3(u, v, normal_dist);
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

    // Create implementation object
    auto impl = std::make_shared<FlatPatchImpl>(origin, world_u, world_v);

    // Create surface using std::bind to member functions
    using namespace std::placeholders;
    return std::make_shared<Surface3D>(
        std::bind(&FlatPatchImpl::position, impl, _1),
        std::bind(&FlatPatchImpl::du, impl, _1),
        std::bind(&FlatPatchImpl::dv, impl, _1),
        std::bind(&FlatPatchImpl::duu, impl, _1),
        std::bind(&FlatPatchImpl::duv, impl, _1),
        std::bind(&FlatPatchImpl::dvv, impl, _1),
        std::bind(&FlatPatchImpl::gaussian, impl, _1),
        std::bind(&FlatPatchImpl::mean, impl, _1),
        [impl, vector_length_epsilon, parameter_bound_epsilon](const WorldPoint3& start, const WorldVector3& dir, double max_t) {
            return impl->solve_path(start, dir, max_t, vector_length_epsilon, parameter_bound_epsilon);
        },
        SurfaceType::Developable,
        std::bind(&FlatPatchImpl::du2_du, impl, _1),
        std::bind(&FlatPatchImpl::du2_dv, impl, _1),
        std::bind(&FlatPatchImpl::duv_du, impl, _1),
        std::bind(&FlatPatchImpl::duv_dv, impl, _1),
        std::bind(&FlatPatchImpl::dv2_du, impl, _1),
        std::bind(&FlatPatchImpl::dv2_dv, impl, _1)
    );
}

} // namespace surfaces
} // namespace shap
