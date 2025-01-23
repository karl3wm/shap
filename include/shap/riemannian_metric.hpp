#pragma once
#include "coord.hpp"
#include "validation_config.hpp"
#include <array>
#include <cmath>
#include <stdexcept>

namespace shap {

class Surface;

/**
 * RiemannianMetric describes how distances and angles in parameter space (u,v) 
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
class RiemannianMetric {
public:
    /**
     * Constructor for direct component specification.
     * Use this when you know the exact metric components and their derivatives.
     * 
     * @param du2 (∂x/∂u)•(∂x/∂u): Square of how u parameter affects world distance
     * @param duv (∂x/∂u)•(∂x/∂v): Interaction between u and v directions
     * @param dv2 (∂x/∂v)•(∂x/∂v): Square of how v parameter affects world distance
     * @param du2_du Derivative of (∂x/∂u)•(∂x/∂u) with respect to u
     * @param du2_dv Derivative of (∂x/∂u)•(∂x/∂u) with respect to v
     * @param duv_du Derivative of (∂x/∂u)•(∂x/∂v) with respect to u
     * @param duv_dv Derivative of (∂x/∂u)•(∂x/∂v) with respect to v
     * @param dv2_du Derivative of (∂x/∂v)•(∂x/∂v) with respect to u
     * @param dv2_dv Derivative of (∂x/∂v)•(∂x/∂v) with respect to v
     */
    RiemannianMetric(
        double du2, double duv, double dv2,
        double du2_du, double du2_dv,
        double duv_du, double duv_dv,
        double dv2_du, double dv2_dv
    ) noexcept;

    /**
     * Constructor from surface derivatives.
     * Computes metric components and derivatives from the surface geometry.
     * 
     * @param surface Reference to the surface
     * @param param_point Point where metric is evaluated
     */
    RiemannianMetric(const Surface& surface, const ParamPoint2& param_point);

    /**
     * Get metric component at given indices.
     * 
     * @param i First index (0=u, 1=v)
     * @param j Second index (0=u, 1=v)
     * @return The metric component gij
     */
    [[nodiscard]] double g(int i, int j) const noexcept;

    /**
     * Get partial derivative of metric component with respect to u.
     * 
     * @param i First index (0=u, 1=v)
     * @param j Second index (0=u, 1=v)
     * @return The derivative ∂gij/∂u
     */
    [[nodiscard]] double dg_du(int i, int j) const noexcept;

    /**
     * Get partial derivative of metric component with respect to v.
     * 
     * @param i First index (0=u, 1=v)
     * @param j Second index (0=u, 1=v)
     * @return The derivative ∂gij/∂v
     */
    [[nodiscard]] double dg_dv(int i, int j) const noexcept;

    /**
     * Convert vector components from covariant to contravariant form.
     * This is the inverse operation of lower_indices.
     * 
     * In parameter space:
     * - Covariant: Components measured along the surface (what you physically measure)
     * - Contravariant: Components in terms of parameter changes (how much u,v to move)
     * 
     * @param covariant_vec Vector with covariant components
     * @return Vector with contravariant components
     */
    [[nodiscard]] ParamVector2 raise_indices(const ParamVector2& covariant_vec) const;

    /**
     * Convert vector components from contravariant to covariant form.
     * This is the inverse operation of raise_indices.
     * 
     * In parameter space:
     * - Contravariant: Components in terms of parameter changes (how much u,v to move)
     * - Covariant: Components measured along the surface (what you physically measure)
     * 
     * @param contravariant_vec Vector with contravariant components
     * @return Vector with covariant components
     */
    [[nodiscard]] ParamVector2 lower_indices(const ParamVector2& contravariant_vec) const noexcept;

    /**
     * Compute first kind Christoffel symbols for geodesic equations.
     * These describe how the metric tensor changes as you move on the surface.
     * 
     * Formula: Γ_ijk = 1/2 (∂_i g_jk + ∂_j g_ik - ∂_k g_ij)
     * 
     * Where:
     * - ∂_i means "how much it changes as you move in direction i"
     * - g_jk are the metric components that measure distances and angles
     */
    [[nodiscard]] std::array<double,2> christoffel_first(int i, int j, int k) const noexcept;

    /**
     * Compute second kind Christoffel symbols for geodesic equations.
     * These determine how tangent vectors change as you move along the surface.
     * 
     * Formula: Γ^i_jk = g^im Γ_mjk
     * 
     * Where:
     * - g^im are components of the inverse metric tensor
     * - Γ_mjk are first kind Christoffel symbols
     */
    [[nodiscard]] std::array<double,2> christoffel_second(int i) const noexcept;

    /**
     * Compute determinant of metric tensor.
     * This measures how much area is distorted between parameter and world space.
     */
    [[nodiscard]] double determinant() const noexcept;

    /**
     * Transform a world space vector to parameter space (pullback operation).
     * This decomposes a 3D vector into tangential and normal components, then
     * converts the tangential part to parameter space coordinates.
     */
    [[nodiscard]] ParamVector3 pullback_vector(
        const WorldVector3& world_vec,
        const WorldVector3& world_du,
        const WorldVector3& world_dv,
        const WorldVector3& world_normal
    ) const;

    /**
     * Transform a parameter space vector to world space (pushforward operation).
     * This combines tangential movement in parameter space with an optional
     * normal component to create a full 3D vector.
     */
    [[nodiscard]] WorldVector3 pushforward_vector(
        const ParamVector3& param_vec,
        const WorldVector3& world_du,
        const WorldVector3& world_dv,
        const WorldVector3& world_normal
    ) const;

private:
    // Metric components (squared lengths and dot products of partial derivatives)
    double g11_, g12_, g22_;  // Kept for compatibility with existing code
    
    // Metric derivatives (how squared lengths and dot products change)
    double dg11_du_, dg11_dv_;  // d(∂x/∂u • ∂x/∂u)/du, d(∂x/∂u • ∂x/∂u)/dv
    double dg12_du_, dg12_dv_;  // d(∂x/∂u • ∂x/∂v)/du, d(∂x/∂u • ∂x/∂v)/dv
    double dg22_du_, dg22_dv_;  // d(∂x/∂v • ∂x/∂v)/du, d(∂x/∂v • ∂x/∂v)/dv

    // Verify metric consistency with surface derivatives
    void verify_metric_consistency(
        const WorldVector3& world_du,
        const WorldVector3& world_dv
    ) const;
};

} // namespace shap
