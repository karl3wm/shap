#pragma once
#include "shap/coord.hpp"
#include "shap/geometry_point2.hpp"
#include "shap/surface.hpp"
#include <cmath>
#include <numbers>
#include <functional>

namespace shap {
namespace surfaces {

namespace {

static constexpr double PI = std::numbers::pi;
static constexpr double TWO_PI = 2 * PI;
static constexpr double HALF_PI = PI / 2;

// Helper to compute both sin and cos
[[nodiscard]] static std::pair<double, double> sincos(double x) noexcept {
    return {std::sin(x), std::cos(x)};
}

} // anonymous namespace

/**
 * Create a sphere surface centered at the origin with given radius.
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
 *
 * @param radius Sphere radius (must be positive)
 * @param tangent_epsilon Tolerance for tangent vector length (default: 1e-10)
 * @param surface_distance_epsilon Tolerance for point-to-surface distance (default: 1e-6)
 * @return Shared pointer to sphere surface
 * @throws std::invalid_argument if radius <= 0 or if any epsilon <= 0
 */
[[nodiscard]] std::shared_ptr<Surface> create_sphere(
    double radius = 1.0,
    double tangent_epsilon = 1e-10,
    double surface_distance_epsilon = 1e-6
) {
    if (radius <= 0) {
        throw std::invalid_argument("Sphere radius must be positive");
    }
    if (tangent_epsilon <= 0 || surface_distance_epsilon <= 0) {
        throw std::invalid_argument("Epsilon values must be positive");
    }

    // Create derived class instance to bind member functions from
    class SphereSurfaceImpl {
    public:
        explicit SphereSurfaceImpl(double r) : radius(r) {}

        WorldPoint3 position(const ParamPoint2& local) const {
            const double u = local.u();
            const double v = local.v();

            // Map parameters to angles
            const double phi = u * TWO_PI;    // longitude [0,2π]
            const double theta = v * PI;       // colatitude [0,π]

            // Cache trigonometric values
            const auto [sin_phi, cos_phi] = sincos(phi);
            const auto [sin_theta, cos_theta] = sincos(theta);
            
            // Position (r * sin(θ)cos(φ), r * sin(θ)sin(φ), r * cos(θ))
            const double r_sin_theta = radius * sin_theta;
            return WorldPoint3(
                r_sin_theta * cos_phi,
                r_sin_theta * sin_phi,
                radius * cos_theta
            );
        }

        WorldVector3 du(const ParamPoint2& local) const {
            const double phi = local.u() * TWO_PI;
            const double theta = local.v() * PI;
            const auto [sin_phi, cos_phi] = sincos(phi);
            const double r_sin_theta = radius * std::sin(theta);
            return WorldVector3(
                -r_sin_theta * sin_phi,
                r_sin_theta * cos_phi,
                0
            );
        }

        WorldVector3 dv(const ParamPoint2& local) const {
            const double phi = local.u() * TWO_PI;
            const double theta = local.v() * PI;
            const auto [sin_phi, cos_phi] = sincos(phi);
            const auto [sin_theta, cos_theta] = sincos(theta);
            const double r_cos_theta = radius * cos_theta;
            return WorldVector3(
                r_cos_theta * cos_phi,
                r_cos_theta * sin_phi,
                -radius * sin_theta
            );
        }

        WorldVector3 duu(const ParamPoint2& local) const {
            const double phi = local.u() * TWO_PI;
            const double theta = local.v() * PI;
            const auto [sin_phi, cos_phi] = sincos(phi);
            const double r_sin_theta = radius * std::sin(theta);
            return WorldVector3(
                -r_sin_theta * cos_phi,
                -r_sin_theta * sin_phi,
                0
            );
        }

        WorldVector3 duv(const ParamPoint2& local) const {
            const double phi = local.u() * TWO_PI;
            const double theta = local.v() * PI;
            const auto [sin_phi, cos_phi] = sincos(phi);
            const double r_cos_theta = radius * std::cos(theta);
            return WorldVector3(
                -r_cos_theta * sin_phi,
                r_cos_theta * cos_phi,
                0
            );
        }

        WorldVector3 dvv(const ParamPoint2& local) const {
            const double phi = local.u() * TWO_PI;
            const double theta = local.v() * PI;
            const auto [sin_phi, cos_phi] = sincos(phi);
            const auto [sin_theta, cos_theta] = sincos(theta);
            const double r_sin_theta = radius * sin_theta;
            return WorldVector3(
                -r_sin_theta * cos_phi,
                -r_sin_theta * sin_phi,
                -radius * cos_theta
            );
        }

        double gaussian(const ParamPoint2&) const {
            return 1.0 / (radius * radius);
        }

        double mean(const ParamPoint2&) const {
            return 1.0 / radius;
        }

        ParamPoint3 world_to_param(const WorldPoint3& pos) const {
            // Get distance from origin
            const double r = pos.length();
            if (r < ValidationConfig::instance().vector_length_epsilon()) {
                throw std::invalid_argument("Cannot compute parameters for zero position vector");
            }
            
            // Get signed distance from sphere surface
            const double normal_dist = r - radius;
            
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

        std::optional<PathIntersection> solve_path(
            const WorldPoint3& start, const WorldVector3& dir, double max_t, double tangent_epsilon) const {
            // Project direction onto tangent plane at start point
            const double start_length = start.length();
            if (start_length < tangent_epsilon) {
                return std::nullopt;  // Start point too close to origin
            }
            const WorldVector3 surface_normal = start * (1.0 / start_length);
            WorldVector3 tangent = dir - dir.dot(surface_normal) * surface_normal;
            const double tangent_length = tangent.length();
            if (tangent_length < tangent_epsilon) {
                return std::nullopt;  // Direction perpendicular to surface
            }
            tangent = tangent * (1.0 / tangent_length);
            
            // Great circle radius = sphere radius
            // Distance = radius * angle
            const double angle = max_t / radius;
            
            // No intersection if we don't complete half circle
            if (angle <= PI) {
                return std::nullopt;
            }
            
            // Convert start point to spherical coordinates
            const double v = std::acos(std::clamp(start.z() / radius, -1.0, 1.0));
            double u = std::atan2(start.y(), start.x());
            if (u < 0) u += TWO_PI;
            
            // Find intersection parameters
            const ParamBound bound = (v < HALF_PI) ? ParamBound::Upper : ParamBound::Lower;
            const double pole_z = (v < HALF_PI) ? radius : -radius;
            
            return PathIntersection(
                HALF_PI * radius,           // Time to reach pole
                WorldPoint3(0, 0, pole_z),   // Pole position
                ParamIndex::V,               // Vertical parameter
                bound,                       // Upper/lower bound
                u / TWO_PI                   // Normalized longitude
            );
        }

        double du2_du(const ParamPoint2&) const { return 0.0; }
        double du2_dv(const ParamPoint2& p) const {
            const double theta = p.v() * PI;
            return radius * radius * std::sin(2 * theta);
        }
        double duv_du(const ParamPoint2&) const { return 0.0; }
        double duv_dv(const ParamPoint2&) const { return 0.0; }
        double dv2_du(const ParamPoint2&) const { return 0.0; }
        double dv2_dv(const ParamPoint2&) const { return 0.0; }

    private:
        double radius;
    };

    // Create implementation object
    auto impl = std::make_shared<SphereSurfaceImpl>(radius);

    // Create surface using std::bind to member functions
    using namespace std::placeholders;
    return std::make_shared<Surface>(
        std::bind(&SphereSurfaceImpl::position, impl, _1),
        std::bind(&SphereSurfaceImpl::du, impl, _1),
        std::bind(&SphereSurfaceImpl::dv, impl, _1),
        std::bind(&SphereSurfaceImpl::duu, impl, _1),
        std::bind(&SphereSurfaceImpl::duv, impl, _1),
        std::bind(&SphereSurfaceImpl::dvv, impl, _1),
        std::bind(&SphereSurfaceImpl::gaussian, impl, _1),
        std::bind(&SphereSurfaceImpl::mean, impl, _1),
        std::bind(&SphereSurfaceImpl::world_to_param, impl, _1),
        [impl, tangent_epsilon](const WorldPoint3& start, const WorldVector3& dir, double max_t) {
            return impl->solve_path(start, dir, max_t, tangent_epsilon);
        },
        SurfaceType::Smooth,
        std::bind(&SphereSurfaceImpl::du2_du, impl, _1),
        std::bind(&SphereSurfaceImpl::du2_dv, impl, _1),
        std::bind(&SphereSurfaceImpl::duv_du, impl, _1),
        std::bind(&SphereSurfaceImpl::duv_dv, impl, _1),
        std::bind(&SphereSurfaceImpl::dv2_du, impl, _1),
        std::bind(&SphereSurfaceImpl::dv2_dv, impl, _1)
    );
}

} // namespace surfaces
} // namespace shap
