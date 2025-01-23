#include "coord.hpp"
#pragma once
#include "../geometry_point2.hpp"
#include "../surface.hpp"
#include <cmath>
#include <numbers>

namespace shap {
namespace surfaces {

/**
 * A sphere surface centered at the origin with given radius.
 *
 * The sphere is parameterized using spherical coordinates:
 * u ∈ [0,1] maps to longitude [0,2π]
 * v ∈ [0,1] maps to colatitude [0,π]
 *
 * Properties:
 * - Constant Gaussian curvature K = 1/r²
 * - Constant mean curvature H = 1/r
 * - Geodesics are great circles
 * - Singularities at poles (v=0 and v=1)
 */
class SphereSurface final : public Surface {
public:
    /**
     * Create a sphere with given radius and tolerances.
     *
     * @param r Sphere radius (must be positive)
     * @param tangent_epsilon Tolerance for tangent vector length (default: 1e-10)
     * @param surface_distance_epsilon Tolerance for point-to-surface distance (default: 1e-6)
     * @throws std::invalid_argument if r <= 0 or if any epsilon <= 0
     */
    explicit SphereSurface(
        double r,
        double tangent_epsilon = 1e-10,
        double surface_distance_epsilon = 1e-6
    ) {
        if (r <= 0) {
            throw std::invalid_argument("Sphere radius must be positive");
        }
        if (tangent_epsilon <= 0 || surface_distance_epsilon <= 0) {
            throw std::invalid_argument("Epsilon values must be positive");
        }
        radius_ = r;
        tangent_epsilon_ = tangent_epsilon;
        surface_distance_epsilon_ = surface_distance_epsilon;
        setup_path_solver();
    }

    // Move operations
    SphereSurface(SphereSurface&&) noexcept = default;
    SphereSurface& operator=(SphereSurface&&) noexcept = default;

    // Prevent copying
    SphereSurface(const SphereSurface&) = delete;
    SphereSurface& operator=(const SphereSurface&) = delete;

    [[nodiscard]] GeometryPoint2 evaluate(const ParamPoint2& local) const override {
        const double u = local.u();
        const double v = local.v();

        // Map parameters to angles
        const double phi = u * TWO_PI;    // longitude [0,2π]
        const double theta = v * PI;       // colatitude [0,π]

        // Cache trigonometric values
        const auto [sin_phi, cos_phi] = sincos(phi);
        const auto [sin_theta, cos_theta] = sincos(theta);
        
        // Position (r * sin(θ)cos(φ), r * sin(θ)sin(φ), r * cos(θ))
        const double r_sin_theta = radius_ * sin_theta;
        const WorldPoint3 world_pos(
            r_sin_theta * cos_phi,
            r_sin_theta * sin_phi,
            radius_ * cos_theta
        );
        
        // Normal points outward from origin (unit vector in radial direction)
        const WorldVector3 world_normal(
            sin_theta * cos_phi,
            sin_theta * sin_phi,
            cos_theta
        );
        
        // First derivatives
        // ∂/∂φ = r * sin(θ) * (-sin(φ), cos(φ), 0)
        const WorldVector3 world_du(
            -r_sin_theta * sin_phi,
            r_sin_theta * cos_phi,
            0
        );
        
        // ∂/∂θ = r * (cos(θ)cos(φ), cos(θ)sin(φ), -sin(θ))
        const double r_cos_theta = radius_ * cos_theta;
        const WorldVector3 world_dv(
            r_cos_theta * cos_phi,
            r_cos_theta * sin_phi,
            -r_sin_theta
        );
        
        // Second derivatives
        // ∂²/∂φ² = -r * sin(θ) * (cos(φ), sin(φ), 0)
        const WorldVector3 world_duu(
            -r_sin_theta * cos_phi,
            -r_sin_theta * sin_phi,
            0
        );
        
        // ∂²/∂φ∂θ = r * cos(θ) * (-sin(φ), cos(φ), 0)
        const WorldVector3 world_duv(
            -r_cos_theta * sin_phi,
            r_cos_theta * cos_phi,
            0
        );
        
        // ∂²/∂θ² = -r * (sin(θ)cos(φ), sin(θ)sin(φ), cos(θ))
        const WorldVector3 world_dvv(
            -r_sin_theta * cos_phi,
            -r_sin_theta * sin_phi,
            -r_cos_theta
        );
        
        // Constant curvature values
        const double inv_r = 1.0 / radius_;
        const double inv_r2 = inv_r * inv_r;
        
        return GeometryPoint2(
            local,
            world_pos,
            world_normal,
            world_du,
            world_dv,
            world_duu,
            world_duv,
            world_dvv,
            inv_r2,                        // Gaussian curvature
            inv_r,                         // Mean curvature
            std::make_pair(inv_r, inv_r)   // Principal curvatures
        );
    }

    [[nodiscard]] ParamPoint3 world_to_param(const WorldPoint3& pos) const override {
        // Get distance from origin
        const double r = pos.length();
        
        // Get signed distance from sphere surface
        const double normal_dist = r - radius_;
        
        // Normalize position to unit sphere for parameter computation
        const double inv_r = 1.0 / r;
        const double x = pos.x() * inv_r;
        const double y = pos.y() * inv_r;
        const double z = pos.z() * inv_r;
        
        // Convert to spherical coordinates
        double v = std::acos(std::clamp(z, -1.0, 1.0));  // colatitude [0,π]
        double u = std::atan2(y, x);                      // longitude [-π,π]
        
        // Normalize u to [0,2π]
        if (u < 0) u += TWO_PI;
        
        // Convert to parameter space [0,1]×[0,1]
        return ParamPoint3(u / TWO_PI, v / PI, normal_dist);
    }

    [[nodiscard]] std::optional<PathSolver> get_path_solver() const noexcept override {
        return path_solver_;
    }

    [[nodiscard]] SurfaceType surface_type() const noexcept override {
        return SurfaceType::Smooth;
    }

    // Access radius
    [[nodiscard]] double radius() const noexcept { return radius_; }

private:
    static constexpr double PI = std::numbers::pi;
    static constexpr double TWO_PI = 2 * PI;
    static constexpr double HALF_PI = PI / 2;

    // Helper to compute both sin and cos
    [[nodiscard]] static std::pair<double, double> sincos(double x) noexcept {
        return {std::sin(x), std::cos(x)};
    }

    void setup_path_solver() noexcept {
        path_solver_ = [this](const WorldPoint3& start, const WorldVector3& dir, double max_t)
            -> std::optional<PathIntersection> {
            // Project direction onto tangent plane at start point
            const WorldVector3 surface_normal = start * (1.0 / start.length());
            WorldVector3 tangent = dir - dir.dot(surface_normal) * surface_normal;
            const double tangent_length = tangent.length();
            if (tangent_length < tangent_epsilon_) {
                return std::nullopt;  // Direction perpendicular to surface
            }
            tangent = tangent * (1.0 / tangent_length);
            
            // Great circle radius = sphere radius
            // Distance = radius * angle
            const double angle = max_t / radius_;
            
            // No intersection if we don't complete half circle
            if (angle <= PI) {
                return std::nullopt;
            }
            
            // Convert start point to spherical coordinates
            const double v = std::acos(std::clamp(start.z() / radius_, -1.0, 1.0));
            double u = std::atan2(start.y(), start.x());
            if (u < 0) u += TWO_PI;
            
            // Find intersection parameters
            const ParamBound bound = (v < HALF_PI) ? ParamBound::Upper : ParamBound::Lower;
            const double pole_z = (v < HALF_PI) ? radius_ : -radius_;
            
            return PathIntersection(
                HALF_PI * radius_,           // Time to reach pole
                WorldPoint3(0, 0, pole_z),   // Pole position
                ParamIndex::V,               // Vertical parameter
                bound,                       // Upper/lower bound
                u / TWO_PI                   // Normalized longitude
            );
        };
    }

    double radius_;
    double tangent_epsilon_;
    double surface_distance_epsilon_;
    PathSolver path_solver_;
};

/**
 * Create a sphere surface with the given radius.
 * 
 * @param radius Sphere radius (must be positive)
 * @param tangent_epsilon Tolerance for tangent vector length (default: 1e-10)
 * @param surface_distance_epsilon Tolerance for point-to-surface distance (default: 1e-6)
 * @return Shared pointer to sphere surface
 * @throws std::invalid_argument if radius <= 0 or if any epsilon <= 0
 */
[[nodiscard]] inline std::shared_ptr<Surface> create_sphere(
    double radius = 1.0,
    double tangent_epsilon = 1e-10,
    double surface_distance_epsilon = 1e-6
) {
    return std::make_shared<SphereSurface>(radius, tangent_epsilon, surface_distance_epsilon);
}

} // namespace surfaces
} // namespace shap
