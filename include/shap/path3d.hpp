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
class Path3D : public Manifold<1, 3, WorldSpaceTag> {
public:
    using Base = Manifold<1, 3, WorldSpaceTag>;
    using PositionFunc = std::function<WorldPoint3(const ParamPoint1&)>;
    using TangentFunc = std::function<WorldVector3(const ParamPoint1&)>;
    using NormalFunc = std::function<WorldVector3(const ParamPoint1&)>;

    static std::shared_ptr<Path3D> create(
        PositionFunc position,
        TangentFunc tangent,
        NormalFunc normal
    ) {
        // Use make_shared to ensure proper enable_shared_from_this setup
        return std::make_shared<Path3D>(std::move(position), std::move(tangent), std::move(normal));
    }

    Path3D(
        PositionFunc position,
        TangentFunc tangent,
        NormalFunc normal
    ) : Base(
            std::bind(&Path3D::world_position, this, std::placeholders::_1),
            std::bind(&Path3D::world_derivative, this, std::placeholders::_1, std::placeholders::_2),
            std::bind(&Path3D::nearest_point_not_implemented, this, std::placeholders::_1)
        )
      , position_(std::move(position))
      , tangent_(std::move(tangent))
      , normal_(std::move(normal))
    {
        if (!position_ || !tangent_ || !normal_) {
            throw std::invalid_argument("Path functions cannot be null");
        }
    }

    /**
     * Get the normal vector at a parameter point.
     * This provides additional geometric information useful for path operations.
     * 
     * @param param 1D parameter in [0,1]
     * @return Normal vector at the point
     * @throws std::invalid_argument if parameter is invalid
     */
    [[nodiscard]] WorldVector3 
    normal(const ParamPoint1& param) const;

    /**
     * Get the binormal vector at a parameter point.
     * Computed as cross product of tangent and normal vectors.
     * 
     * @param param 1D parameter in [0,1]
     * @return Binormal vector at the point
     * @throws std::invalid_argument if parameter is invalid
     */
    [[nodiscard]] WorldVector3
    binormal(const ParamPoint1& param) const;

private:
    // Manifold interface implementation
    [[nodiscard]] WorldPoint3 world_position(const ParamPoint1& param) const {
        return position_(param);
    }

    [[nodiscard]] WorldVector3 world_derivative(const ParamPoint1& param, int) const {
        return tangent_(param);
    }

    PositionFunc position_;
    TangentFunc tangent_;
    NormalFunc normal_;

    // Helper to extract parameter value and validate
    [[nodiscard]] static double get_param_value(const ParamPoint1& param);

    [[nodiscard]] ParamPoint1 nearest_point_not_implemented(const WorldPoint3&) const {
        throw std::runtime_error("Nearest point calculation not implemented for paths");
    }
};

} // namespace shap
