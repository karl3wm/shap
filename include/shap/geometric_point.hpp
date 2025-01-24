#pragma once
#include "manifold.hpp"
#include <array>
#include <memory>

namespace shap {

/**
 * Bundles geometric information about a point on a manifold.
 * Replaces GeometryPoint2 with a more general template that works for both paths and surfaces.
 * 
 * @tparam ParamDim Dimension of parameter space (1 for paths, 2/3 for surfaces)
 * @tparam TargetDim Dimension of target space (2 or 3 for world space)
 * @tparam SpaceTag WorldSpaceTag or ParamSpaceTag for target space
 */
template<int ParamDim, int TargetDim, typename SpaceTag>
class GeometricPoint {
public:
    using ManifoldType = Manifold<ParamDim, TargetDim, SpaceTag>;
    using ParameterPoint = typename ManifoldType::ParameterPoint;
    using TargetPoint = typename ManifoldType::TargetPoint;
    using TargetVector = typename ManifoldType::TargetVector;

    GeometricPoint(
        const ManifoldType* manifold,
        ParameterPoint local_pos,
        TargetPoint world_pos,
        const TargetVector* derivatives
    ) : manifold_(manifold)
      , local_pos_(std::move(local_pos))
      , world_pos_(std::move(world_pos)) {
        if (!manifold) {
            throw std::invalid_argument("Manifold pointer cannot be null");
        }
        if (!derivatives) {
            throw std::invalid_argument("Derivatives pointer cannot be null");
        }
        for (int i = 0; i < ParamDim; ++i) {
            derivatives_[i] = derivatives[i];
        }
    }

    // Accessors
    [[nodiscard]] const ManifoldType* manifold() const noexcept { return manifold_; }
    [[nodiscard]] const ParameterPoint& local_pos() const noexcept { return local_pos_; }
    [[nodiscard]] const TargetPoint& world_pos() const noexcept { return world_pos_; }
    
    // Get first derivatives
    [[nodiscard]] const TargetVector* derivatives() const noexcept {
        return derivatives_;
    }

    // Convenience accessors for specific derivatives
    template<int Dim = ParamDim>
    [[nodiscard]] std::enable_if_t<(Dim >= 1), const TargetVector&> 
    d_u() const noexcept {
        return derivatives_[0];
    }

    template<int Dim = ParamDim>
    [[nodiscard]] std::enable_if_t<(Dim >= 2), const TargetVector&>
    d_v() const noexcept {
        return derivatives_[1];
    }

    template<int Dim = ParamDim>
    [[nodiscard]] std::enable_if_t<(Dim >= 3), const TargetVector&>
    d_w() const noexcept {
        return derivatives_[2];
    }

private:
    const ManifoldType* manifold_;
    ParameterPoint local_pos_;
    TargetPoint world_pos_;
    TargetVector derivatives_[ParamDim];
};

} // namespace shap
