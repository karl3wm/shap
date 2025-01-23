#pragma once
#include "coord.hpp"
#include "validation_config.hpp"
#include <array>
#include <cmath>
#include <functional>
#include <vector>

namespace shap {

/**
 * The metric tensor describes how distances and angles in parameter space (u,v) 
 * relate to distances and angles in world space (x,y,z).
 * 
 * For a surface embedded in 3D space:
 * - Parameter space is a 2D flat grid with coordinates (u,v)
 * - At each point (u,v), we have two fundamental vectors:
 *   ∂x/∂u: How position changes as u increases (keeping v constant)
 *   ∂x/∂v: How position changes as v increases (keeping u constant)
 * 
 * The metric components measure how these vectors interact:
 *   g11 = (∂x/∂u)•(∂x/∂u): Square of how much world distance changes per unit u
 *   g12 = (∂x/∂u)•(∂x/∂v): How u and v directions interact (0 if perpendicular)
 *   g22 = (∂x/∂v)•(∂x/∂v): Square of how much world distance changes per unit v
 */
class Surface2DMetricTensor {
public:
    // Function type for metric components that vary with position
    using MetricFunction = std::function<double(double,double)>;
    using MetricArray = std::array<std::array<MetricFunction, 2>, 2>;
    using DerivativeArray = std::array<std::array<MetricFunction, 2>, 2>;

    /**
     * Constructor for constant coefficient metric.
     * Use this when the relationship between parameter space and world space is uniform.
     * 
     * @param g11 (∂x/∂u)•(∂x/∂u): How u parameter affects world distance
     * @param g12 (∂x/∂u)•(∂x/∂v): Interaction between u and v directions
     * @param g21 Same as g12 (metric tensor is symmetric)
     * @param g22 (∂x/∂v)•(∂x/∂v): How v parameter affects world distance
     */
    Surface2DMetricTensor(double g11, double g12, double g21, double g22) noexcept
        : metric_fns_{{
            {[g11](double,double) noexcept { return g11; },
             [g12](double,double) noexcept { return g12; }},
            {[g21](double,double) noexcept { return g21; },
             [g22](double,double) noexcept { return g22; }}
          }},
          has_derivatives_(false) {}

    /**
     * Constructor for variable coefficient metric with optional derivatives.
     * Use this when the relationship between spaces varies with position.
     * 
     * @param g11_fn Function computing (∂x/∂u)•(∂x/∂u) at each point
     * @param g12_fn Function computing (∂x/∂u)•(∂x/∂v) at each point
     * @param g21_fn Same as g12_fn (metric tensor is symmetric)
     * @param g22_fn Function computing (∂x/∂v)•(∂x/∂v) at each point
     * @param dg11_du Optional: How g11 changes with u
     * @param dg11_dv Optional: How g11 changes with v
     * @param dg12_du Optional: How g12 changes with u
     * @param dg12_dv Optional: How g12 changes with v
     * @param dg22_du Optional: How g22 changes with u
     * @param dg22_dv Optional: How g22 changes with v
     */
    Surface2DMetricTensor(
        MetricFunction g11_fn,
        MetricFunction g12_fn,
        MetricFunction g21_fn,
        MetricFunction g22_fn,
        MetricFunction dg11_du = nullptr,
        MetricFunction dg11_dv = nullptr,
        MetricFunction dg12_du = nullptr,
        MetricFunction dg12_dv = nullptr,
        MetricFunction dg22_du = nullptr,
        MetricFunction dg22_dv = nullptr
    ) noexcept
        : metric_fns_{{
            {std::move(g11_fn), std::move(g12_fn)},
            {std::move(g21_fn), std::move(g22_fn)}
          }},
          du_fns_{{
            {std::move(dg11_du), std::move(dg12_du)},
            {std::move(dg12_du), std::move(dg22_du)}
          }},
          dv_fns_{{
            {std::move(dg11_dv), std::move(dg12_dv)},
            {std::move(dg12_dv), std::move(dg22_dv)}
          }},
          has_derivatives_(dg11_du && dg11_dv && dg12_du &&
                         dg12_dv && dg22_du && dg22_dv) {}

    /**
     * Get metric coefficient at given parameters.
     * 
     * @param i First index (0=u, 1=v)
     * @param j Second index (0=u, 1=v)
     * @param u First parameter coordinate
     * @param v Second parameter coordinate
     * @return The metric component gij at (u,v)
     */
    [[nodiscard]] double g(int i, int j, double u, double v) const noexcept {
        if (i < 0 || i > 1 || j < 0 || j > 1) return 0.0;
        return metric_fns_[i][j](u, v);
    }

    /**
     * Get partial derivative of metric coefficient with respect to u.
     * Uses numerical approximation if exact derivatives not provided.
     */
    [[nodiscard]] double dg_du(int i, int j, double u, double v, double h = 1e-7) const noexcept {
        if (i < 0 || i > 1 || j < 0 || j > 1) return 0.0;
        
        if (!has_derivatives_) {
            return (g(i,j, u+h, v) - g(i,j, u-h, v)) / (2*h);
        }
        
        return du_fns_[i][j](u, v);
    }

    /**
     * Get partial derivative of metric coefficient with respect to v.
     * Uses numerical approximation if exact derivatives not provided.
     */
    [[nodiscard]] double dg_dv(int i, int j, double u, double v, double h = 1e-7) const noexcept {
        if (i < 0 || i > 1 || j < 0 || j > 1) return 0.0;
        
        if (!has_derivatives_) {
            return (g(i,j, u, v+h) - g(i,j, u, v-h)) / (2*h);
        }
        
        return dv_fns_[i][j](u, v);
    }

    /**
     * Convert vector components from covariant to contravariant form.
     * This is the inverse operation of lower_indices.
     * 
     * In parameter space:
     * - Covariant: Components measured along the surface (what you physically measure)
     * - Contravariant: Components in terms of parameter changes (how much u,v to move)
     * 
     * @param covariant_vec Vector with covariant components
     * @param param_point Point where the conversion happens
     * @param epsilon Small number for numerical comparisons
     * @return Vector with contravariant components
     */
    [[nodiscard]] ParamVector2 raise_indices(
        const ParamVector2& covariant_vec,
        const ParamPoint2& param_point
    ) const {
        const double det = determinant(param_point.u(), param_point.v());
        const double validation_epsilon = ValidationConfig::instance().vector_length_epsilon();
        if (std::abs(det) < validation_epsilon) {
            throw std::runtime_error("Degenerate metric tensor");
        }

        const double g11 = metric_fns_[0][0](param_point.u(), param_point.v());
        const double g12 = metric_fns_[0][1](param_point.u(), param_point.v());
        const double g22 = metric_fns_[1][1](param_point.u(), param_point.v());

        return ParamVector2(
            (g22 * covariant_vec.u() - g12 * covariant_vec.v()) / det,
            (-g12 * covariant_vec.u() + g11 * covariant_vec.v()) / det
        );
    }

    /**
     * Convert vector components from contravariant to covariant form.
     * This is the inverse operation of raise_indices.
     * 
     * In parameter space:
     * - Contravariant: Components in terms of parameter changes (how much u,v to move)
     * - Covariant: Components measured along the surface (what you physically measure)
     * 
     * @param contravariant_vec Vector with contravariant components
     * @param param_point Point where the conversion happens
     * @return Vector with covariant components
     */
    [[nodiscard]] ParamVector2 lower_indices(
        const ParamVector2& contravariant_vec,
        const ParamPoint2& param_point
    ) const noexcept {
        const double g11 = metric_fns_[0][0](param_point.u(), param_point.v());
        const double g12 = metric_fns_[0][1](param_point.u(), param_point.v());
        const double g22 = metric_fns_[1][1](param_point.u(), param_point.v());

        return ParamVector2(
            g11 * contravariant_vec.u() + g12 * contravariant_vec.v(),
            g12 * contravariant_vec.u() + g22 * contravariant_vec.v()
        );
    }

    /**
     * Compute first kind Christoffel symbols for geodesic equations.
     * These describe how the metric tensor changes as you move on the surface.
     */
    [[nodiscard]] std::array<double,2> christoffel_first(
        int i, 
        int j, 
        int k, 
        double u, 
        double v,
        double h = 1e-7
    ) const noexcept;

    /**
     * Compute second kind Christoffel symbols for geodesic equations.
     * These determine how tangent vectors change as you move along the surface.
     */
    [[nodiscard]] std::array<double,2> christoffel_second(
        int i, 
        double u, 
        double v,
        double epsilon = 1e-10
    ) const noexcept;

    /**
     * Compute determinant of metric tensor at given parameters.
     * This measures how much area is distorted between parameter and world space.
     */
    [[nodiscard]] double determinant(double u, double v) const noexcept {
        const double g11 = metric_fns_[0][0](u, v);
        const double g12 = metric_fns_[0][1](u, v);
        const double g21 = metric_fns_[1][0](u, v);
        const double g22 = metric_fns_[1][1](u, v);
        return g11 * g22 - g12 * g21;
    }

private:
    // Verify that the metric tensor matches the surface derivatives
    void verify_metric_consistency(
        const WorldVector3& world_du,
        const WorldVector3& world_dv,
        const ParamPoint2& param_point
    ) const {
        const double g11 = metric_fns_[0][0](param_point.u(), param_point.v());
        const double g12 = metric_fns_[0][1](param_point.u(), param_point.v());
        const double g22 = metric_fns_[1][1](param_point.u(), param_point.v());

        const double computed_g11 = world_du.dot(world_du);
        const double computed_g12 = world_du.dot(world_dv);
        const double computed_g22 = world_dv.dot(world_dv);

        const double validation_epsilon = ValidationConfig::instance().vector_length_epsilon();
        if (std::abs(g11 - computed_g11) > validation_epsilon ||
            std::abs(g12 - computed_g12) > validation_epsilon ||
            std::abs(g22 - computed_g22) > validation_epsilon) {
            throw std::runtime_error("Metric tensor inconsistent with surface derivatives");
        }
    }

public:
    /**
     * Transform a world space vector to parameter space (pullback operation).
     * This decomposes a 3D vector into tangential and normal components, then
     * converts the tangential part to parameter space coordinates.
     * 
     * @param world_vec The vector in 3D world space
     * @param world_du Surface derivative ∂x/∂u at the point
     * @param world_dv Surface derivative ∂x/∂v at the point
     * @param world_normal Unit normal vector at the point
     * @param param_point The parameter point where this happens
     * @return ParamVector3 with (u,v) as tangential components and w as normal component
     */
    [[nodiscard]] ParamVector3 pullback_vector(
        const WorldVector3& world_vec,
        const WorldVector3& world_du,
        const WorldVector3& world_dv,
        const WorldVector3& world_normal,
        const ParamPoint2& param_point
    ) const {
        // Get normal component
        const double normal_component = world_vec.dot(world_normal);
        
        // Get tangential component
        const WorldVector3 tangent_vec = world_vec - world_normal * normal_component;

        verify_metric_consistency(world_du, world_dv, param_point);

        // Get contravariant components in parameter space using dot products
        const ParamVector2 tangent_params(
            tangent_vec.dot(world_du),
            tangent_vec.dot(world_dv)
        );

        // Convert to parameter space using raise_indices
        const ParamVector2 param_vec = raise_indices(tangent_params, param_point);

        return ParamVector3(param_vec.u(), param_vec.v(), normal_component);
    }

    /**
     * Transform a parameter space vector to world space (pushforward operation).
     * This combines tangential movement in parameter space with an optional
     * normal component to create a full 3D vector.
     * 
     * @param param_vec Vector in parameter space (u,v components for tangential, w for normal)
     * @param world_du Surface derivative ∂x/∂u at the point
     * @param world_dv Surface derivative ∂x/∂v at the point
     * @param world_normal Unit normal vector at the point
     * @param param_point The parameter point where this happens
     * @return The vector in world space
     */
    [[nodiscard]] WorldVector3 pushforward_vector(
        const ParamVector3& param_vec,
        const WorldVector3& world_du,
        const WorldVector3& world_dv,
        const WorldVector3& world_normal,
        const ParamPoint2& param_point
    ) const {
        verify_metric_consistency(world_du, world_dv, param_point);

        // Convert tangential components to world space
        const WorldVector3 tangent_vec = 
            world_du * param_vec.u() + 
            world_dv * param_vec.v();

        // Add normal component if present
        return tangent_vec + world_normal * param_vec.w();
    }

private:
    MetricArray metric_fns_;      // The metric tensor components
    DerivativeArray du_fns_;      // Derivatives with respect to u
    DerivativeArray dv_fns_;      // Derivatives with respect to v
    bool has_derivatives_;        // Whether exact derivatives are available
};

} // namespace shap
