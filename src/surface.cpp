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
    constexpr int PATH_POINTS = 20;  // Number of points to sample along path

    // Helper to compute normal from derivatives
    [[nodiscard]] WorldVector3 compute_normal(
        const WorldVector3& du, 
        const WorldVector3& dv
    ) noexcept {
        return du.cross(dv).normalize();
    }
} // anonymous namespace

Surface::Surface(
    PositionFunction position_func,
    DerivativeFunction du_func,
    DerivativeFunction dv_func,
    DerivativeFunction duu_func,
    DerivativeFunction duv_func,
    DerivativeFunction dvv_func,
    CurvatureFunction gaussian_func,
    CurvatureFunction mean_func,
    WorldToParamFunction world_to_param_func,
    std::optional<PathSolver> path_solver,
    SurfaceType type,
    ParameterSpaceDerivative du2_du,
    ParameterSpaceDerivative du2_dv,
    ParameterSpaceDerivative duv_du,
    ParameterSpaceDerivative duv_dv,
    ParameterSpaceDerivative dv2_du,
    ParameterSpaceDerivative dv2_dv
) : position_func_(std::move(position_func))
  , du_func_(std::move(du_func))
  , dv_func_(std::move(dv_func))
  , duu_func_(std::move(duu_func))
  , duv_func_(std::move(duv_func))
  , dvv_func_(std::move(dvv_func))
  , gaussian_curv_func_(std::move(gaussian_func))
  , mean_curv_func_(std::move(mean_func))
  , world_to_param_func_(std::move(world_to_param_func))
  , path_solver_(std::move(path_solver))
  , type_(type)
  , du2_du_fn_(std::move(du2_du))
  , du2_dv_fn_(std::move(du2_dv))
  , duv_du_fn_(std::move(duv_du))
  , duv_dv_fn_(std::move(duv_dv))
  , dv2_du_fn_(std::move(dv2_du))
  , dv2_dv_fn_(std::move(dv2_dv)) {
    if (!position_func_ || !du_func_ || !dv_func_ || 
        !duu_func_ || !duv_func_ || !dvv_func_ ||
        !gaussian_curv_func_ || !mean_curv_func_ ||
        !world_to_param_func_) {
        throw std::invalid_argument("Required functions cannot be null");
    }
}

[[nodiscard]] GeometryPoint2 Surface::evaluate(const ParamPoint2& local) const {
        const WorldPoint3 position = position_func_(local);
        const WorldVector3 du = du_func_(local);
        const WorldVector3 dv = dv_func_(local);
        const WorldVector3 normal = compute_normal(du, dv);
        
        if (type_ == SurfaceType::Smooth) {
            const WorldVector3 duu = duu_func_(local);
            const WorldVector3 duv = duv_func_(local);
            const WorldVector3 dvv = dvv_func_(local);
            
            const double gaussian = gaussian_curv_func_(local);
            const double mean = mean_curv_func_(local);
            
            // Compute principal curvatures
            const double disc = mean*mean - gaussian;
            std::pair<double, double> principal{0.0, 0.0};
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
        
        return GeometryPoint2(
            this,
            local,
            position,
            normal,
            du,
            dv
        );
    }

std::unique_ptr<SurfacePath> Surface::create_path(
    const GeometryPoint2& start,
    const WorldVector3& world_direction,
    double world_length
) const {
    if (world_length <= 0) {
        throw std::invalid_argument("Path length must be positive");
    }
    if (world_direction.length_squared() < ValidationConfig::instance().vector_length_epsilon()) {
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
    if (tangent_dir.length_squared() < ValidationConfig::instance().vector_length_epsilon()) {
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

WorldVector3 Surface::world_to_parameter_velocity(
    const WorldVector3& world_direction,
    const WorldVector3& world_du,
    const WorldVector3& world_dv
) const noexcept {
    // Solve linear system to convert world direction to parameter velocity
    const WorldVector3 normal = world_du.cross(world_dv);
    const double det = normal.length();
    if (det < ValidationConfig::instance().vector_length_epsilon()) {
        return WorldVector3(0, 0, 0);  // Degenerate case
    }
    
    // Use Cramer's rule to solve the system:
    // world_direction = du_dt * world_du + dv_dt * world_dv
    const WorldVector3 normalized_normal = normal * (1.0 / det);
    const double du_dt = world_direction.cross(world_dv).dot(normalized_normal);
    const double dv_dt = world_du.cross(world_direction).dot(normalized_normal);
    
    return WorldVector3(du_dt, dv_dt, 0);
}

} // namespace shap
