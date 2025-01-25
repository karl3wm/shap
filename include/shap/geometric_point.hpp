#pragma once
#include "manifold.hpp"
#include "edge_descriptor.hpp"
#include "frame_vectors.hpp"
#include "riemannian_metric.hpp"

namespace shap {
class Surface3D;

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
      , world_pos_(std::move(world_pos))
      , frame_(ParamDim >= 2 ? 
            FrameVectors(derivatives[0], derivatives[1]) :  // 2D case
            FrameVectors(derivatives[0], derivatives[0].orthogonal()))  // 1D case
      , metric_(frame_.du().dot(frame_.du()), // g11
                frame_.du().dot(frame_.dv()), // g12
                frame_.dv().dot(frame_.dv()), // g22
                0.0, 0.0, // du2_du, du2_dv (TODO: compute these from second derivatives)
                0.0, 0.0, // duv_du, duv_dv
                0.0, 0.0) // dv2_du, dv2_dv
    {
        if (!manifold) {
            throw std::invalid_argument("Manifold pointer cannot be null");
        }
        if (!derivatives) {
            throw std::invalid_argument("Derivatives pointer cannot be null");
        }
    }

    // Accessors
    [[nodiscard]] const ManifoldType* manifold() const noexcept { return manifold_; }
    [[nodiscard]] const ParameterPoint& local_pos() const noexcept { return local_pos_; }
    [[nodiscard]] const TargetPoint& world_pos() const noexcept { return world_pos_; }
    
    // Frame vector accessors
    [[nodiscard]] const FrameVectors& frame() const noexcept { return frame_; }
    [[nodiscard]] const RiemannianMetric& metric() const noexcept { return metric_; }

    // Convenience accessors for specific derivatives
    template<int Dim = ParamDim>
    [[nodiscard]] std::enable_if_t<(Dim >= 1), const TargetVector&> 
    d_u() const noexcept { return frame_.du(); }

    template<int Dim = ParamDim>
    [[nodiscard]] std::enable_if_t<(Dim >= 2), const TargetVector&>
    d_v() const noexcept { return frame_.dv(); }

    // Transform operations delegated to frame vectors
    [[nodiscard]] TargetVector pushforward(const ParamVector3& param_vec) const {
        return frame_.pushforward(param_vec);
    }

    [[nodiscard]] ParamVector3 pullback(const TargetVector& target_vec) const {
        return frame_.pullback(target_vec);
    }

    // Surface-specific functionality (only available for 2D manifolds)
    template<int Dim = ParamDim>
    [[nodiscard]] std::enable_if_t<Dim == 2, bool>
    is_on_edge() const noexcept {
        const auto& param = local_pos_;
        return param.u() <= 0.0 || param.u() >= 1.0 ||
               param.v() <= 0.0 || param.v() >= 1.0;
    }

    template<int Dim = ParamDim>
    [[nodiscard]] std::enable_if_t<Dim == 2, std::optional<EdgeDescriptor>>
    get_edge_descriptor() const noexcept {
        const auto& param = local_pos_;
        if (param.u() <= 0.0) {
            return EdgeDescriptor{ParamIndex::U, ParamBound::Lower, param.v()};
        }
        if (param.u() >= 1.0) {
            return EdgeDescriptor{ParamIndex::U, ParamBound::Upper, param.v()};
        }
        if (param.v() <= 0.0) {
            return EdgeDescriptor{ParamIndex::V, ParamBound::Lower, param.u()};
        }
        if (param.v() >= 1.0) {
            return EdgeDescriptor{ParamIndex::V, ParamBound::Upper, param.u()};
        }
        return std::nullopt;
    }

    template<int Dim = ParamDim>
    [[nodiscard]] std::enable_if_t<Dim == 2, const Surface3D*>
    surface() const noexcept {
        return dynamic_cast<const Surface3D*>(manifold_);
    }

    template<int Dim = ParamDim>
    [[nodiscard]] std::enable_if_t<Dim == 2, TargetVector>
    world_normal() const noexcept {
        return frame_.normal();
    }

private:
    const ManifoldType* manifold_;
    ParameterPoint local_pos_;
    TargetPoint world_pos_;
    FrameVectors frame_;
    RiemannianMetric metric_;
};

} // namespace shap
