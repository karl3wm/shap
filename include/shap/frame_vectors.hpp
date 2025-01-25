#pragma once
#include "coord.hpp"

namespace shap {

/**
 * Represents the frame vectors that define the local coordinate system at a point on a manifold.
 * Provides access to basis vectors and normal vector without mixing in metric operations.
 */
class FrameVectors {
public:
    /**
     * Construct frame vectors from surface derivatives.
     * 
     * @param du First basis vector (∂x/∂u)
     * @param dv Second basis vector (∂x/∂v)
     */
    FrameVectors(const WorldVector3& du, const WorldVector3& dv) noexcept
        : du_(du)
        , dv_(dv)
        , normal_(du.crossed(dv).normalized())
    {}

    // Frame vector accessors
    [[nodiscard]] const WorldVector3& du() const noexcept { return du_; }
    [[nodiscard]] const WorldVector3& dv() const noexcept { return dv_; }
    [[nodiscard]] const WorldVector3& normal() const noexcept { return normal_; }

    /**
     * Transform a vector from parameter space to target space (pushforward operation).
     * 
     * @param param_vec Vector in parameter space with normal component
     * @return Vector in target (world) space
     */
    [[nodiscard]] WorldVector3 pushforward(const ParamVector3& param_vec) const {
        return du_ * param_vec.u() + dv_ * param_vec.v() + normal_ * param_vec.w();
    }

    /**
     * Transform a vector from target space to parameter space (pullback operation).
     * 
     * @param world_vec Vector in target (world) space to transform
     * @return Vector in parameter space with normal component
     */
    [[nodiscard]] ParamVector3 pullback(const WorldVector3& world_vec) const {
        // Project onto tangent plane basis vectors
        double u = world_vec.dot(du_);
        double v = world_vec.dot(dv_);
        // Get normal component
        double w = world_vec.dot(normal_);
        return ParamVector3(u, v, w);
    }

private:
    WorldVector3 du_;      // First basis vector (∂x/∂u)
    WorldVector3 dv_;      // Second basis vector (∂x/∂v)
    WorldVector3 normal_;  // Normal vector (du × dv normalized)
};

} // namespace shap
