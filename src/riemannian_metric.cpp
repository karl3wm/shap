#include "shap/riemannian_metric.hpp"
#include "shap/surface.hpp"
#include <cmath>

namespace shap {

RiemannianMetric::RiemannianMetric(
    double du2, double duv, double dv2,
    double du2_du, double du2_dv,
    double duv_du, double duv_dv,
    double dv2_du, double dv2_dv
) noexcept
    : g11_(du2), g12_(duv), g22_(dv2)
    , dg11_du_(du2_du), dg11_dv_(du2_dv)
    , dg12_du_(duv_du), dg12_dv_(duv_dv)
    , dg22_du_(dv2_du), dg22_dv_(dv2_dv)
{}

RiemannianMetric::RiemannianMetric(const Surface& surface, const ParamPoint2& param_point) {
    const auto geom = surface.evaluate(param_point);
    const auto& world_du = geom.world_du();
    const auto& world_dv = geom.world_dv();

    // Compute metric components from surface derivatives
    g11_ = world_du.dot(world_du);
    g12_ = world_du.dot(world_dv);
    g22_ = world_dv.dot(world_dv);

    // Get metric component derivatives from surface
    dg11_du_ = surface.du2_du(param_point);
    dg11_dv_ = surface.du2_dv(param_point);
    dg12_du_ = surface.duv_du(param_point);
    dg12_dv_ = surface.duv_dv(param_point);
    dg22_du_ = surface.dv2_du(param_point);
    dg22_dv_ = surface.dv2_dv(param_point);
}

double RiemannianMetric::g(int i, int j) const noexcept {
    if (i < 0 || i > 1 || j < 0 || j > 1) return 0.0;
    if (i == 0 && j == 0) return g11_;
    if ((i == 0 && j == 1) || (i == 1 && j == 0)) return g12_;
    if (i == 1 && j == 1) return g22_;
    return 0.0;
}

double RiemannianMetric::dg_du(int i, int j) const noexcept {
    if (i < 0 || i > 1 || j < 0 || j > 1) return 0.0;
    if (i == 0 && j == 0) return dg11_du_;
    if ((i == 0 && j == 1) || (i == 1 && j == 0)) return dg12_du_;
    if (i == 1 && j == 1) return dg22_du_;
    return 0.0;
}

double RiemannianMetric::dg_dv(int i, int j) const noexcept {
    if (i < 0 || i > 1 || j < 0 || j > 1) return 0.0;
    if (i == 0 && j == 0) return dg11_dv_;
    if ((i == 0 && j == 1) || (i == 1 && j == 0)) return dg12_dv_;
    if (i == 1 && j == 1) return dg22_dv_;
    return 0.0;
}

ParamVector2 RiemannianMetric::raise_indices(const ParamVector2& covariant_vec) const {
    const double det = determinant();
    const double validation_epsilon = ValidationConfig::instance().vector_length_epsilon();
    if (std::abs(det) < validation_epsilon) {
        throw std::runtime_error("Degenerate metric tensor");
    }

    return ParamVector2(
        (g22_ * covariant_vec.u() - g12_ * covariant_vec.v()) / det,
        (-g12_ * covariant_vec.u() + g11_ * covariant_vec.v()) / det
    );
}

ParamVector2 RiemannianMetric::lower_indices(const ParamVector2& contravariant_vec) const noexcept {
    return ParamVector2(
        g11_ * contravariant_vec.u() + g12_ * contravariant_vec.v(),
        g12_ * contravariant_vec.u() + g22_ * contravariant_vec.v()
    );
}

std::array<double,2> RiemannianMetric::christoffel_first(
    int i, int j, int k
) const noexcept {
    if (i < 0 || i > 1 || j < 0 || j > 1 || k < 0 || k > 1) {
        return {0.0, 0.0};
    }

    // Get how the metric changes in each direction
    const double dg_i = (i == 0) ? dg_du(j,k) : dg_dv(j,k);
    const double dg_j = (j == 0) ? dg_du(i,k) : dg_dv(i,k);
    const double dg_k = (k == 0) ? dg_du(i,j) : dg_dv(i,j);
    
    // Combine the changes using the fundamental formula
    const double gamma = 0.5 * (dg_i + dg_j - dg_k);
    
    return {gamma, 0.0};
}

std::array<double,2> RiemannianMetric::christoffel_second(int i) const noexcept {
    if (i < 0 || i > 1) {
        return {0.0, 0.0};
    }
    
    // Check if the metric is well-behaved
    const double det = determinant();
    const double validation_epsilon = ValidationConfig::instance().vector_length_epsilon();
    if (std::abs(det) < validation_epsilon) {
        return {0.0, 0.0}; // Surface is degenerate here
    }
    
    // Compute the inverse metric components
    const double inv_det = 1.0 / det;
    const double g11_inv = g22_ * inv_det;
    const double g12_inv = -g12_ * inv_det;
    const double g21_inv = -g12_ * inv_det;
    const double g22_inv = g11_ * inv_det;
    
    // Get the first kind symbols we need
    const auto gamma_1 = christoffel_first(0, i, i);
    const auto gamma_2 = christoffel_first(1, i, i);
    
    // Transform to second kind symbols using the inverse metric
    return {
        g11_inv * gamma_1[0] + g12_inv * gamma_2[0],  // u component
        g21_inv * gamma_1[0] + g22_inv * gamma_2[0]   // v component
    };
}

double RiemannianMetric::determinant() const noexcept {
    return g11_ * g22_ - g12_ * g12_;
}

void RiemannianMetric::verify_metric_consistency(
    const WorldVector3& world_du,
    const WorldVector3& world_dv
) const {
    const double computed_g11 = world_du.dot(world_du);
    const double computed_g12 = world_du.dot(world_dv);
    const double computed_g22 = world_dv.dot(world_dv);

    const double validation_epsilon = ValidationConfig::instance().vector_length_epsilon();
    if (std::abs(g11_ - computed_g11) > validation_epsilon ||
        std::abs(g12_ - computed_g12) > validation_epsilon ||
        std::abs(g22_ - computed_g22) > validation_epsilon) {
        throw std::runtime_error("Metric tensor inconsistent with surface derivatives");
    }
}

ParamVector3 RiemannianMetric::pullback_vector(
    const WorldVector3& world_vec,
    const WorldVector3& world_du,
    const WorldVector3& world_dv,
    const WorldVector3& world_normal
) const {
    // Get normal component
    const double normal_component = world_vec.dot(world_normal);
    
    // Get tangential component
    const WorldVector3 tangent_vec = world_vec - world_normal * normal_component;

    verify_metric_consistency(world_du, world_dv);

    // Get contravariant components in parameter space using dot products
    const ParamVector2 tangent_params(
        tangent_vec.dot(world_du),
        tangent_vec.dot(world_dv)
    );

    // Convert to parameter space using raise_indices
    const ParamVector2 param_vec = raise_indices(tangent_params);

    return ParamVector3(param_vec.u(), param_vec.v(), normal_component);
}

WorldVector3 RiemannianMetric::pushforward_vector(
    const ParamVector3& param_vec,
    const WorldVector3& world_du,
    const WorldVector3& world_dv,
    const WorldVector3& world_normal
) const {
    verify_metric_consistency(world_du, world_dv);

    // Convert tangential components to world space
    const WorldVector3 tangent_vec = 
        world_du * param_vec.u() + 
        world_dv * param_vec.v();

    // Add normal component if present
    return tangent_vec + world_normal * param_vec.w();
}

} // namespace shap
