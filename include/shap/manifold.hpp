#pragma once
#include "coord.hpp"
#include <stdexcept>
#include <functional>

namespace shap {

// Forward declarations
template<int ParamDim, int TargetDim, typename SpaceTag>
class GeometricPoint;

class Path3D;
class SurfacePath;

/** 
 * Base class for manifolds that map from parameter space to target space.
 * 
 * @tparam ParamDim Dimension of parameter space (1 for paths, 2/3 for surfaces)
 * @tparam TargetDim Dimension of target space (2 or 3 for world space)
 * @tparam SpaceTag WorldSpaceTag or ParamSpaceTag for target space
 */
template<int ParamDim, int TargetDim, typename SpaceTag>
class Manifold {
public:
    static_assert(ParamDim >= 1 && ParamDim <= 3, "Parameter space must be 1D, 2D, or 3D");
    static_assert(TargetDim >= 2 && TargetDim <= 3, "Target space must be 2D or 3D");
    static_assert(ParamDim <= TargetDim || std::is_same_v<SpaceTag, ParamSpaceTag>,
                 "Parameter dimension cannot exceed target dimension in world space");

    using ParameterPoint = Coord<ParamDim, PointTag, ParamSpaceTag>;
    using TargetPoint = Coord<TargetDim, PointTag, SpaceTag>;
    using TargetVector = Coord<TargetDim, VectorTag, SpaceTag>;

    using WorldPositionFunc = std::function<TargetPoint(const ParameterPoint& param)>;
    using DerivativesFunc = std::function<TargetVector(const ParameterPoint& param, int derivative_index)>;

    Manifold(WorldPositionFunc world_position_func, DerivativesFunc derivatives_func)
        : world_position_func_(world_position_func), derivatives_func_(derivatives_func) {}

    [[nodiscard]] virtual GeometricPoint<ParamDim, TargetDim, SpaceTag> 
    evaluate(const ParameterPoint& param) const {
        TargetPoint world_position = world_position_func_(param);
        TargetVector derivs[ParamDim];
        for (int i = 0; i < ParamDim; ++i) {
            derivs[i] = derivatives_func_(param, i);
        }
        return GeometricPoint<ParamDim, TargetDim, SpaceTag>(this, param, world_position, derivs);
    }

    [[nodiscard]] virtual std::array<TargetVector, ParamDim>
    derivatives(const ParameterPoint& param) const {
        std::array<TargetVector, ParamDim> derivs;
        for (int i = 0; i < ParamDim; ++i) {
            derivs[i] = derivatives_func_(param, i);
        }
        return derivs;
    }

    virtual ~Manifold() = default;

private:
    WorldPositionFunc world_position_func_;
    DerivativesFunc derivatives_func_;
};

// Common manifold type aliases
using WorldPath3D = Manifold<1, 3, WorldSpaceTag>;    // 1D -> 3D world
using WorldPath2D = Manifold<1, 2, WorldSpaceTag>;    // 1D -> 2D world
using Surface2D = Manifold<2, 2, WorldSpaceTag>;      // 2D -> 2D world
using ParamPath2D = Manifold<1, 2, ParamSpaceTag>;    // 1D -> 2D param
using ParamPath3D = Manifold<1, 3, ParamSpaceTag>;    // 1D -> 3D param
using ParamSurface2D = Manifold<2, 2, ParamSpaceTag>; // 2D -> 2D param
using ParamSurface3D = Manifold<3, 3, ParamSpaceTag>; // 3D -> 3D param

} // namespace shap
