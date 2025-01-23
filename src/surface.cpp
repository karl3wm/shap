#include "shap/coord.hpp"
#include "shap/surface.hpp"
#include "shap/geometry_point2.hpp"
#include "shap/path.hpp"
#include <limits>
#include <cmath>
#include <array>
#include <iostream>

namespace shap {

namespace {
    constexpr double EPSILON = 1e-10;
    constexpr int GRID_SIZE = 10;
    constexpr int MAX_ITERATIONS = 20;
    constexpr double GRADIENT_STEP = 0.01;
    constexpr int PATH_POINTS = 20;  // Increased for better accuracy

    // Adaptive step size for numerical derivatives
    [[nodiscard]] constexpr double compute_step_size(double x) noexcept {
        const double eps = std::numeric_limits<double>::epsilon();
        return std::cbrt(eps) * (1.0 + std::abs(x));
    }

    // Helper to compute normal from derivatives
    [[nodiscard]] WorldVector3 compute_normal(
        const WorldVector3& du, 
        const WorldVector3& dv
    ) noexcept {
        return du.cross(dv).normalize();
    }

    // Helper to compute curvature coefficients
    struct CurvatureCoefficients {
        double E, F, G;  // First fundamental form
        double L, M, N;  // Second fundamental form
        double det;      // EG - F²

        [[nodiscard]] static CurvatureCoefficients compute(
            const WorldVector3& du,
            const WorldVector3& dv,
            const WorldVector3& duu,
            const WorldVector3& duv,
            const WorldVector3& dvv,
            const WorldVector3& normal
        ) noexcept {
            CurvatureCoefficients coeff;
            coeff.E = du.dot(du);
            coeff.F = du.dot(dv);
            coeff.G = dv.dot(dv);
            coeff.L = duu.dot(normal);
            coeff.M = duv.dot(normal);
            coeff.N = dvv.dot(normal);
            coeff.det = coeff.E * coeff.G - coeff.F * coeff.F;
            return coeff;
        }

        [[nodiscard]] std::optional<double> gaussian_curvature() const noexcept {
            if (std::abs(det) <= EPSILON) return std::nullopt;
            return (L * N - M * M) / det;
        }

        [[nodiscard]] std::optional<double> mean_curvature() const noexcept {
            if (std::abs(det) <= EPSILON) return std::nullopt;
            return (E * N - 2.0 * F * M + G * L) / (2.0 * det);
        }
    };
} // anonymous namespace

class FunctionSurface final : public Surface {
public:
    FunctionSurface(
        PositionFunction pos,
        std::optional<DerivativeFunction> du = std::nullopt,
        std::optional<DerivativeFunction> dv = std::nullopt,
        std::optional<DerivativeFunction> duu = std::nullopt,
        std::optional<DerivativeFunction> duv = std::nullopt,
        std::optional<DerivativeFunction> dvv = std::nullopt,
        std::optional<CurvatureFunction> gaussian = std::nullopt,
        std::optional<CurvatureFunction> mean = std::nullopt,
        std::optional<PathSolver> path_solver = std::nullopt,
        SurfaceType type = SurfaceType::Smooth,
        std::optional<MetricDerivativeFunction> du2_du = std::nullopt,
        std::optional<MetricDerivativeFunction> du2_dv = std::nullopt,
        std::optional<MetricDerivativeFunction> duv_du = std::nullopt,
        std::optional<MetricDerivativeFunction> duv_dv = std::nullopt,
        std::optional<MetricDerivativeFunction> dv2_du = std::nullopt,
        std::optional<MetricDerivativeFunction> dv2_dv = std::nullopt
    ) noexcept
        : position_func_(std::move(pos))
        , du_func_(std::move(du))
        , dv_func_(std::move(dv))
        , duu_func_(std::move(duu))
        , duv_func_(std::move(duv))
        , dvv_func_(std::move(dvv))
        , gaussian_curv_func_(std::move(gaussian))
        , mean_curv_func_(std::move(mean))
        , path_solver_(std::move(path_solver))
        , type_(type) {
        // Initialize metric component derivative functions
        du2_du_fn_ = std::move(du2_du).value_or(nullptr);
        du2_dv_fn_ = std::move(du2_dv).value_or(nullptr);
        duv_du_fn_ = std::move(duv_du).value_or(nullptr);
        duv_dv_fn_ = std::move(duv_dv).value_or(nullptr);
        dv2_du_fn_ = std::move(dv2_du).value_or(nullptr);
        dv2_dv_fn_ = std::move(dv2_dv).value_or(nullptr);
    }

    [[nodiscard]] GeometryPoint2 evaluate(const ParamPoint2& local) const override {
        WorldVector3 du(0.0, 0.0, 0.0), dv(0.0, 0.0, 0.0);
        
        // Compute first derivatives
        if (du_func_ && dv_func_) {
            du = (*du_func_)(local);
            dv = (*dv_func_)(local);
        } else {
            // Adaptive step size numerical derivatives
            const double hu = compute_step_size(local.u());
            const double hv = compute_step_size(local.v());
            
            const auto u_plus = ParamPoint2(local.u() + hu, local.v());
            const auto u_minus = ParamPoint2(local.u() - hu, local.v());
            const auto v_plus = ParamPoint2(local.u(), local.v() + hv);
            const auto v_minus = ParamPoint2(local.u(), local.v() - hv);
            
            du = (position_func_(u_plus) - position_func_(u_minus)) * (0.5 / hu);  // Point subtraction returns vector
            dv = (position_func_(v_plus) - position_func_(v_minus)) * (0.5 / hv);  // Point subtraction returns vector
        }

        const WorldVector3 normal = compute_normal(du, dv);
        const WorldPoint3 position = position_func_(local);
        
        // For smooth surfaces, compute second derivatives and curvature
        if (type_ == SurfaceType::Smooth) {
            WorldVector3 duu(0.0, 0.0, 0.0), duv(0.0, 0.0, 0.0), dvv(0.0, 0.0, 0.0);
            
            if (duu_func_) {
                duu = (*duu_func_)(local);
            } else {
                const double hu = compute_step_size(local.u());
                const auto u_plus = ParamPoint2(local.u() + hu, local.v());
                const auto u_minus = ParamPoint2(local.u() - hu, local.v());
                duu = ((position_func_(u_plus) - position) - (position - position_func_(u_minus))) * (1.0 / (hu * hu));
            }
            
            if (duv_func_) {
                duv = (*duv_func_)(local);
            } else {
                const double hu = compute_step_size(local.u());
                const double hv = compute_step_size(local.v());
                const auto uv_plus = ParamPoint2(local.u() + hu, local.v() + hv);
                const auto uv_minus_u = ParamPoint2(local.u() + hu, local.v() - hv);
                const auto uv_minus_v = ParamPoint2(local.u() - hu, local.v() + hv);
                const auto uv_minus = ParamPoint2(local.u() - hu, local.v() - hv);
                duv = ((position_func_(uv_plus) - position_func_(uv_minus_u)) -
                      (position_func_(uv_minus_v) - position_func_(uv_minus))) * 
                     (0.25 / (hu * hv));  // Point subtraction returns vector
            }
            
            if (dvv_func_) {
                dvv = (*dvv_func_)(local);
            } else {
                const double hv = compute_step_size(local.v());
                const auto v_plus = ParamPoint2(local.u(), local.v() + hv);
                const auto v_minus = ParamPoint2(local.u(), local.v() - hv);
                dvv = ((position_func_(v_plus) - position) - (position - position_func_(v_minus))) * (1.0 / (hv * hv));
            }

            // Compute curvature
            const auto coeffs = CurvatureCoefficients::compute(
                du, dv, duu, duv, dvv, normal);

            double gaussian = 0.0;
            double mean = 0.0;
            std::pair<double, double> principal{0.0, 0.0};

            if (gaussian_curv_func_) {
                gaussian = (*gaussian_curv_func_)(local);
            } else if (auto k = coeffs.gaussian_curvature()) {
                gaussian = *k;
            }

            if (mean_curv_func_) {
                mean = (*mean_curv_func_)(local);
            } else if (auto h = coeffs.mean_curvature()) {
                mean = *h;
            }

            // Compute principal curvatures
            const double disc = mean*mean - gaussian;
            if (disc >= 0) {
                const double sqrt_disc = std::sqrt(disc);
                if (mean >= 0) {
                    const double k1 = mean + sqrt_disc;
                    const double k2 = gaussian / k1;  // More stable than mean - sqrt_disc
                    principal = std::make_pair(k1, k2);
                } else {
                    const double k2 = mean - sqrt_disc;
                    const double k1 = gaussian / k2;  // More stable than mean + sqrt_disc
                    principal = std::make_pair(k1, k2);
                }
            }

            return GeometryPoint2(
                this,
                local,
                position,
                normal,
                du,
                dv,
                duu,
                duv,
                dvv,
                gaussian,
                mean,
                principal
            );
        }
        
        // For non-smooth surfaces, return just first derivatives
        return GeometryPoint2(
            this,
            local,
            position,
            normal,
            du,
            dv
        );
    }

    [[nodiscard]] ParamPoint3 world_to_param(const WorldPoint3& pos) const override {
        // Grid search for initial guess
        double best_u = 0, best_v = 0;
        double min_dist = std::numeric_limits<double>::max();
        
        for (int i = 0; i <= GRID_SIZE; ++i) {
            const double u = static_cast<double>(i) / GRID_SIZE;
            for (int j = 0; j <= GRID_SIZE; ++j) {
                const double v = static_cast<double>(j) / GRID_SIZE;
                const auto local = ParamPoint2(u, v);
                const WorldPoint3 surface_pt = position_func_(local);
                const double dist = (surface_pt - pos).length_squared();
                if (dist < min_dist) {
                    min_dist = dist;
                    best_u = u;
                    best_v = v;
                }
            }
        }
        
        // Gradient descent refinement
        for (int iter = 0; iter < MAX_ITERATIONS; ++iter) {
            const auto local = ParamPoint2(best_u, best_v);
            const WorldPoint3 curr_pos = position_func_(local);
            const WorldVector3 diff = pos - curr_pos;
            if (diff.length_squared() < EPSILON) break;
            
            // Compute numerical derivatives
            const double hu = compute_step_size(best_u);
            const double hv = compute_step_size(best_v);
            
            const auto u_plus = ParamPoint2(best_u + hu, best_v);
            const auto v_plus = ParamPoint2(best_u, best_v + hv);
            
            const WorldVector3 du = (position_func_(u_plus) - curr_pos) * (1.0/hu);
            const WorldVector3 dv = (position_func_(v_plus) - curr_pos) * (1.0/hv);
            
            // Update parameters
            best_u = std::clamp(best_u + GRADIENT_STEP * diff.dot(du), 0.0, 1.0);
            best_v = std::clamp(best_v + GRADIENT_STEP * diff.dot(dv), 0.0, 1.0);
        }
        
        // Get final position and compute normal distance
        const auto local = ParamPoint2(best_u, best_v);
        const auto geom = evaluate(local);
        const WorldVector3 diff = pos - geom.world_pos();
        const double normal_dist = diff.dot(geom.world_normal());
        
        return ParamPoint3(best_u, best_v, normal_dist);
    }

    [[nodiscard]] std::optional<PathSolver> get_path_solver() const noexcept override {
        return path_solver_;
    }

    [[nodiscard]] SurfaceType surface_type() const noexcept override {
        return type_;
    }

private:
    PositionFunction position_func_;
    std::optional<DerivativeFunction> du_func_;
    std::optional<DerivativeFunction> dv_func_;
    std::optional<DerivativeFunction> duu_func_;
    std::optional<DerivativeFunction> duv_func_;
    std::optional<DerivativeFunction> dvv_func_;
    std::optional<CurvatureFunction> gaussian_curv_func_;
    std::optional<CurvatureFunction> mean_curv_func_;
    std::optional<PathSolver> path_solver_;
    SurfaceType type_;
};

std::unique_ptr<SurfacePath> Surface::create_path(
    const GeometryPoint2& start,
    const WorldVector3& world_direction,
    double world_length
) const {
    if (world_length <= 0) {
        throw std::invalid_argument("Path length must be positive");
    }
    if (world_direction.length_squared() < EPSILON) {
        throw std::invalid_argument("Direction vector cannot be zero");
    }

    auto path = std::make_unique<PathSegment>(
        std::shared_ptr<Surface>(const_cast<Surface*>(this), [](Surface*){})
    );
    
    // Project direction onto surface tangent plane
    const auto start_geom = evaluate(start.local_pos());
    
    std::cout << "\nCreate Path Analysis:\n"
              << "World direction: " << world_direction.x() << ", "
              << world_direction.y() << ", " << world_direction.z()
              << " (length=" << world_direction.length() << ")\n"
              << "Surface derivatives at start:\n"
              << "  du = (" << start_geom.world_du().x() << ", "
              << start_geom.world_du().y() << ", " << start_geom.world_du().z()
              << ") length=" << start_geom.world_du().length() << "\n"
              << "  dv = (" << start_geom.world_dv().x() << ", "
              << start_geom.world_dv().y() << ", " << start_geom.world_dv().z()
              << ") length=" << start_geom.world_dv().length() << "\n"
              << "  normal = (" << start_geom.world_normal().x() << ", "
              << start_geom.world_normal().y() << ", " << start_geom.world_normal().z()
              << ")\n";
    
    WorldVector3 tangent_dir = world_direction - 
        world_direction.dot(start_geom.world_normal()) * start_geom.world_normal();
    if (tangent_dir.length_squared() < EPSILON) {
        throw std::runtime_error("Direction is perpendicular to surface");
    }
    
    std::cout << "Projected direction: " << tangent_dir.x() << ", "
              << tangent_dir.y() << ", " << tangent_dir.z()
              << " (length=" << tangent_dir.length() << ")\n";
    
    tangent_dir = tangent_dir.normalize();
    std::cout << "Normalized direction: " << tangent_dir.x() << ", "
              << tangent_dir.y() << ", " << tangent_dir.z() << "\n";

    // Convert direction to parameter space using metric tensor
    std::cout << "\nComputing parameter velocity for direction...\n";
    const auto param_vel = world_to_parameter_velocity(
        tangent_dir,  // Convert direction first, then scale by length
        start_geom.world_du(),
        start_geom.world_dv()
    );
    
    // Scale parameter velocity by world length
    const auto scaled_vel = param_vel * world_length;
    
    std::cout << "Parameter velocity (includes world length):\n"
              << "  du/dt = " << param_vel.x() << "\n"
              << "  dv/dt = " << param_vel.y() << "\n";
    
    // Compute end parameters using scaled velocity
    const auto& start_local = start.local_pos();
    const auto end_local = ParamPoint2(
        start_local.u() + scaled_vel.x(),
        start_local.v() + scaled_vel.y()
    );
    
    std::cout << "\nParameter space coordinates:\n"
              << "Start: u=" << start_local.u() << " v=" << start_local.v() << "\n"
              << "End: u=" << end_local.u() << " v=" << end_local.v() << "\n"
              << "Delta: du=" << (end_local.u() - start_local.u())
              << " dv=" << (end_local.v() - start_local.v()) << "\n";

    // Check for surface transitions
    double transition_t = 1.0;  // Normalized t value
    ParamPoint2 transition_local = end_local;

    if (auto solver = get_path_solver()) {
        if (auto intersection = (*solver)(start.world_pos(), tangent_dir, 1.0)) {
            transition_t = intersection->t;
            transition_local = world_to_param_r2(intersection->position);
        }
    }

    // Add start point
    path->add_point(0.0, start_local.u(), start_local.v());
    
    // Get metric tensor at start point for proper scaling
    const double du_scale = start_geom.world_du().length();
    const double dv_scale = start_geom.world_dv().length();
    
    std::cout << "\nPath sampling analysis:\n"
              << "Surface scale factors:\n"
              << "  |du| = " << du_scale << "\n"
              << "  |dv| = " << dv_scale << "\n";
    
    // Sample points with metric-aware interpolation
    WorldPoint3 prev_pos = start.world_pos();
    double accumulated_length = 0.0;
    
    for (int i = 1; i <= PATH_POINTS; ++i) {
        // Use normalized parameter
        const double alpha = static_cast<double>(i) / PATH_POINTS;
        if (alpha > transition_t) break;
        
        // Scale parameter interpolation by metric
        const double u = start_local.u() + param_vel.x() * alpha;
        const double v = start_local.v() + param_vel.y() * alpha;
        
        // Compute actual world position and length
        const auto curr_geom = evaluate(ParamPoint2(u, v));
        const auto curr_pos = curr_geom.world_pos();
        accumulated_length += (curr_pos - prev_pos).length();
        const double t = accumulated_length / world_length;
        
        std::cout << "Sample point " << i << ":\n"
                  << "  alpha = " << alpha << "\n"
                  << "  t = " << t << "\n"
                  << "  u = " << u << "\n"
                  << "  v = " << v << "\n"
                  << "  pos = (" << curr_pos.x() << ", " 
                  << curr_pos.y() << ", " << curr_pos.z() << ")\n"
                  << "  accumulated_length = " << accumulated_length << "\n";
        
        path->add_point(t, u, v);
        prev_pos = curr_pos;
    }
    
    return path;
}

std::shared_ptr<Surface> Surface::create(
    PositionFunction position_func,
    std::optional<PathSolver> path_solver,
    SurfaceType type,
    std::optional<MetricDerivativeFunction> du2_du,
    std::optional<MetricDerivativeFunction> du2_dv,
    std::optional<MetricDerivativeFunction> duv_du,
    std::optional<MetricDerivativeFunction> duv_dv,
    std::optional<MetricDerivativeFunction> dv2_du,
    std::optional<MetricDerivativeFunction> dv2_dv
) {
    if (!position_func) {
        throw std::invalid_argument("Position function cannot be null");
    }
    return std::make_shared<FunctionSurface>(
        std::move(position_func),
        std::nullopt, std::nullopt,
        std::nullopt, std::nullopt, std::nullopt,
        std::nullopt, std::nullopt,
        std::move(path_solver),
        type,
        std::move(du2_du),
        std::move(du2_dv),
        std::move(duv_du),
        std::move(duv_dv),
        std::move(dv2_du),
        std::move(dv2_dv)
    );
}

WorldVector3 Surface::world_to_parameter_velocity(
    const WorldVector3& world_direction,
    const WorldVector3& world_du,
    const WorldVector3& world_dv
) const noexcept {
    // Solve linear system to convert world direction to parameter velocity
    const double det = world_du.cross(world_dv).length();
    if (det < EPSILON) {
        return WorldVector3(0, 0, 0);  // Degenerate case
    }
    
    // Use Cramer's rule to solve the system:
    // world_direction = du_dt * world_du + dv_dt * world_dv
    const double du_dt = world_direction.cross(world_dv).dot(world_du.cross(world_dv).normalize()) / det;
    const double dv_dt = world_du.cross(world_direction).dot(world_du.cross(world_dv).normalize()) / det;
    
    return WorldVector3(du_dt, dv_dt, 0);
}

} // namespace shap
