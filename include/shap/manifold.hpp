#pragma once
#include "coord.hpp"
#include <stdexcept>
#include <functional>
#include <memory>

namespace shap {

// Forward declarations
template<int ParamDim, int TargetDim, typename SpaceTag>
class GeometricPoint;

class Path3D;
class SurfacePath;

/**
 * Base class for manifolds that map from parameter space to target space.
 * 
 * Note on lifetime management: This prototype uses enable_shared_from_this to handle
 * relationships between manifolds (e.g., when one manifold's geometry depends on another).
 * Manifolds must be allocated as shared_ptr to support these relationships.
 * 
 * @tparam ParamDim Dimension of parameter space
 * @tparam TargetDim Dimension of target space
 * @tparam SpaceTag Tag indicating the coordinate space type of the target space
 */
template<int ParamDim, int TargetDim, typename SpaceTag>
class Manifold : public std::enable_shared_from_this<Manifold<ParamDim, TargetDim, SpaceTag>> {
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
    using NearestFunc = std::function<ParameterPoint(const TargetPoint& target)>;

    Manifold(
        WorldPositionFunc world_position_func,
        DerivativesFunc derivatives_func,
        NearestFunc nearest_func
    )
        : world_position_func_(world_position_func)
        , derivatives_func_(derivatives_func)
        , nearest_func_(nearest_func)
    {
    }

    /**
     * Find the nearest point in parameter space to the given target space point.
     */
    [[nodiscard]] ParameterPoint nearest(const TargetPoint& target) const {
        return nearest_func_(target);
    }

    virtual ~Manifold() = default;

    [[nodiscard]] GeometricPoint<ParamDim, TargetDim, SpaceTag> 
    evaluate(const ParameterPoint& param) const {
        TargetPoint world_position = world_position_func_(param);
        TargetVector derivs[ParamDim];
        for (int i = 0; i < ParamDim; ++i) {
            derivs[i] = derivatives_func_(param, i);
        }
        return GeometricPoint<ParamDim, TargetDim, SpaceTag>(this, param, world_position, derivs);
    }

    [[nodiscard]] std::array<TargetVector, ParamDim>
    derivatives(const ParameterPoint& param) const {
        std::array<TargetVector, ParamDim> derivs;
        for (int i = 0; i < ParamDim; ++i) {
            derivs[i] = derivatives_func_(param, i);
        }
        return derivs;
    }
    WorldPositionFunc world_position_func_;
    DerivativesFunc derivatives_func_;
    NearestFunc nearest_func_;
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
