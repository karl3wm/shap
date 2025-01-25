#pragma once
#include "frame_vectors.hpp"
#include "coord.hpp"
#include <array>
#include <cmath>

namespace shap {

/**
 * RiemannianMetric provides the fundamental geometric structure for measuring
 * distances, angles, and converting between parameter and target spaces.
 * 
 * Key concepts from differential geometry:
 * - Tangent vectors (contravariant) specify directions in parameter space
 * - Cotangent vectors (covariant) measure change in target space
 * - The metric tensor enables conversion between these dual spaces
 * 
 * For a surface embedded in 3D space:
 * - Parameter space is a 2D flat grid with coordinates (u,v)
 * - The tangent space at each point has basis vectors:
 *   ∂x/∂u: How position changes as u increases (keeping v constant)
 *   ∂x/∂v: How position changes as v increases (keeping u constant)
 * 
 * The metric components measure how these basis vectors interact:
 *   g11 = (∂x/∂u)•(∂x/∂u): Square of how much world distance changes per unit u
 *   g12 = (∂x/∂u)•(∂x/∂v): How u and v directions interact (0 if perpendicular)
 *   g22 = (∂x/∂v)•(∂x/∂v): Square of how much world distance changes per unit v
 * 
 * These components enable:
 * - Converting between tangent and cotangent vectors (raise/lower indices)
 * - Measuring distances and angles in parameter space
 * - Computing geodesics through Christoffel symbols
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
     * Convert a cotangent vector (covariant components) to a tangent vector (contravariant components).
     * 
     * In differential geometry:
     * - Tangent vectors (contravariant) specify directions in parameter space
     * - Cotangent vectors (covariant) measure change in target space
     * - "Raising indices" refers to converting from covariant to contravariant components
     * 
     * @param cotangent_vec Vector with covariant components in the cotangent space
     * @return Vector with contravariant components in the tangent space
     */
    [[nodiscard]] ParamVector2 raise_indices(const ParamVector2& cotangent_vec) const;

    /**
     * Convert a tangent vector (contravariant components) to a cotangent vector (covariant components).
     * 
     * In differential geometry:
     * - Tangent vectors (contravariant) transform like basis vectors under coordinate changes
     * - Cotangent vectors (covariant) transform like differentials under coordinate changes
     * - "Lowering indices" refers to converting from contravariant to covariant components
     * 
     * @param tangent_vec Vector with contravariant components in the tangent space
     * @return Vector with covariant components in the cotangent space
     */
    [[nodiscard]] ParamVector2 lower_indices(const ParamVector2& tangent_vec) const noexcept;

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

private:
    // Metric components (squared lengths and dot products of partial derivatives)
    double g11_, g12_, g22_;
    
    // Metric derivatives (how squared lengths and dot products change)
    double dg11_du_, dg11_dv_;  // d(∂x/∂u • ∂x/∂u)/du, d(∂x/∂u • ∂x/∂u)/dv
    double dg12_du_, dg12_dv_;  // d(∂x/∂u • ∂x/∂v)/du, d(∂x/∂u • ∂x/∂v)/dv
    double dg22_du_, dg22_dv_;  // d(∂x/∂v • ∂x/∂v)/du, d(∂x/∂v • ∂x/∂v)/dv

    // Helper to verify metric consistency with frame vectors
    void verify_metric_consistency(const FrameVectors& frame) const {
        const auto& du = frame.du();
        const auto& dv = frame.dv();
        
        // Verify metric components match frame vectors
        const double computed_g11 = du.dot(du);
        const double computed_g12 = du.dot(dv);
        const double computed_g22 = dv.dot(dv);

        if (std::abs(computed_g11 - g11_) > ValidationConfig::instance().epsilon() ||
            std::abs(computed_g12 - g12_) > ValidationConfig::instance().epsilon() ||
            std::abs(computed_g22 - g22_) > ValidationConfig::instance().epsilon()) {
            throw std::runtime_error("Metric components inconsistent with frame vectors");
        }
    }
};

} // namespace shap
