#pragma once
#include "coord.hpp"
#include "edge_connection.hpp"
#include "edge_descriptor.hpp"
#include "geometry_point2.hpp"
#include "path_intersection.hpp"
#include "surface_type.hpp"
#include "surface3d.hpp"
#include <functional>
#include <memory>
#include <optional>
#include <utility>
#include <stdexcept>

namespace shap {

class SurfacePoint;
class SurfacePath;

/**
 * Legacy surface class that adapts to the new Surface3D interface.
 * This class will be deprecated once all code is migrated to use Surface3D directly.
 */
class Surface {
    friend class RiemannianMetric;
public:
    using Ptr = std::shared_ptr<Surface>;

    explicit Surface(std::shared_ptr<Surface3D> impl) : impl_(std::move(impl)) {
        if (!impl_) {
            throw std::invalid_argument("Surface3D implementation cannot be null");
        }
    }
    
    // Prevent copying
    Surface(const Surface&) = delete;
    Surface& operator=(const Surface&) = delete;
    
    // Allow moving
    Surface(Surface&&) noexcept = default;
    Surface& operator=(Surface&&) noexcept = default;

    // Get underlying Surface3D implementation
    [[nodiscard]] std::shared_ptr<Surface3D> impl() const noexcept { return impl_; }

    /**
     * Evaluate surface at parameter space point.
     */
    [[nodiscard]] GeometryPoint2 evaluate(const ParamPoint2& local) const {
        auto geom3d = impl_->evaluate(local);
        const auto& du = geom3d.derivatives()[0];
        const auto& dv = geom3d.derivatives()[1];
        // Calculate normal as cross product of derivatives
        WorldVector3 normal = du.crossed(dv).normalized();
        return GeometryPoint2(
            this,
            local,
            geom3d.world_pos(),
            normal,
            du,
            dv
        );
    }
    
    /**
     * Convert a world space position to local coordinates.
     */
    [[nodiscard]] ParamPoint3 world_to_param(const WorldPoint3& pos) const {
        return impl_->world_to_param(pos);
    }
    
    /**
     * Create a path on the surface starting from a point in a given direction.
     */
    [[nodiscard]] std::unique_ptr<SurfacePath> create_path(
        const GeometryPoint2& start,
        const WorldVector3& world_direction,
        double world_length
    ) const;
    
    // Get surface type
    [[nodiscard]] SurfaceType surface_type() const noexcept {
        return impl_->surface_type();
    }
    
    /**
     * Get scale factors for converting between parameter and world space.
     */
    [[nodiscard]] std::pair<double, double> get_scale_factors(
        const ParamPoint2& local
    ) const {
        return impl_->get_scale_factors(local);
    }

    // Metric component derivative accessors
    [[nodiscard]] double du2_du(const ParamPoint2& param) const noexcept {
        return impl_->du2_du(param);
    }
    [[nodiscard]] double du2_dv(const ParamPoint2& param) const noexcept {
        return impl_->du2_dv(param);
    }
    [[nodiscard]] double duv_du(const ParamPoint2& param) const noexcept {
        return impl_->duv_du(param);
    }
    [[nodiscard]] double duv_dv(const ParamPoint2& param) const noexcept {
        return impl_->duv_dv(param);
    }
    [[nodiscard]] double dv2_du(const ParamPoint2& param) const noexcept {
        return impl_->dv2_du(param);
    }
    [[nodiscard]] double dv2_dv(const ParamPoint2& param) const noexcept {
        return impl_->dv2_dv(param);
    }

    /**
     * Get path solver if available.
     */
    [[nodiscard]] std::optional<PathSolver> get_path_solver() const noexcept {
        return impl_->get_path_solver();
    }

    /**
     * Convert world space velocity to parameter space.
     */
    [[nodiscard]] ParamVector2 world_to_parameter_velocity(
        const WorldVector3& world_vel,
        const WorldVector3& du,
        const WorldVector3& dv
    ) const {
        // Project velocity onto tangent plane
        const auto normal = du.crossed(dv).normalized();
        const auto planar_vel = world_vel - world_vel.dot(normal) * normal;
        
        // Solve for parameter velocities using derivatives
        const double det = du.crossed(dv).length();
        if (det < 1e-10) {
            throw std::invalid_argument("Surface derivatives are nearly parallel");
        }
        
        const double du_dt = planar_vel.crossed(dv).dot(normal) / det;
        const double dv_dt = du.crossed(planar_vel).dot(normal) / det;
        
        return ParamVector2(du_dt, dv_dt);
    }

private:
    std::shared_ptr<Surface3D> impl_;
};

} // namespace shap
