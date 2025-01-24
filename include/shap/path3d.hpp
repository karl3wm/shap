#pragma once
#include "shap/coord.hpp"
#include "shap/manifold.hpp"
#include "shap/geometric_point.hpp"
#include "shap/validation_config.hpp"
#include <functional>
#include <memory>
#include <vector>

namespace shap {

/**
 * Represents a 1D path in 3D world space.
 * 
 * A path maps from a 1D parameter space [0,1] to 3D world space.
 * This class provides a concrete implementation of Manifold<1,3,WorldSpaceTag>
 * that uses function objects to define the mapping.
 */
using WorldPath3D = Manifold<1, 3, WorldSpaceTag>;  // 1D -> 3D world
class Path3D : public WorldPath3D {
public:
    using PositionFunc = std::function<WorldPoint3(double)>;
    using TangentFunc = std::function<WorldVector3(double)>;
    using NormalFunc = std::function<WorldVector3(double)>;

    Path3D(
        PositionFunc position,
        TangentFunc tangent,
        NormalFunc normal
    ) : WorldPath3D(
            [f = position](const ParameterPoint& param) {
                return f(get_param_value(param));
            },
            [t = tangent](const ParameterPoint& param, int) {
                return t(get_param_value(param));
            }
        )
      , position_(std::move(position))
      , tangent_(std::move(tangent))
      , normal_(std::move(normal))
    {
        if (!position_ || !tangent_ || !normal_) {
            throw std::invalid_argument("Path functions cannot be null");
        }
    }

    // Use base class implementation of evaluate
    using WorldPath3D::evaluate;

    /**
     * Get path derivatives at parameter point.
     * For a path, this returns a single vector representing the tangent direction.
     * 
     * @param param 1D parameter in [0,1]
     * @return Vector containing single tangent vector
     * @throws std::invalid_argument if parameter is invalid
     */
    [[nodiscard]] std::array<TargetVector, 1>
    derivatives(const ParameterPoint& param) const override {
        return {tangent_(get_param_value(param))};
    }

    /**
     * Get the normal vector at a parameter point.
     * This provides additional geometric information useful for path operations.
     * 
     * @param param 1D parameter in [0,1]
     * @return Normal vector at the point
     * @throws std::invalid_argument if parameter is invalid
     */
    [[nodiscard]] TargetVector 
    normal(const ParameterPoint& param) const;

    /**
     * Get the binormal vector at a parameter point.
     * Computed as cross product of tangent and normal vectors.
     * 
     * @param param 1D parameter in [0,1]
     * @return Binormal vector at the point
     * @throws std::invalid_argument if parameter is invalid
     */
    [[nodiscard]] TargetVector
    binormal(const ParameterPoint& param) const;

private:
    PositionFunc position_;
    TangentFunc tangent_;
    NormalFunc normal_;

    // Helper to extract parameter value and validate
    [[nodiscard]] static double get_param_value(const ParameterPoint& param);
};

} // namespace shap
