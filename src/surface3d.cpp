#include "shap/surface3d.hpp"
#include "shap/geometric_point.hpp"
#include "shap/surface.hpp"
#include "shap/path.hpp"

namespace shap {

namespace {
    // Helper function to validate parameter bounds
    void validate_parameters(const ParamPoint2& param) {
        if (param.u() < 0.0 || param.u() > 1.0 || param.v() < 0.0 || param.v() > 1.0) {
            throw std::invalid_argument("Parameter values must be in [0,1]");
        }
    }
}

std::pair<double, double>
Surface3DImpl::get_scale_factors(const ParameterPoint& local) const {
    validate_parameters(local);
    
    // Scale factors are lengths of first derivatives
    WorldVector3 du = du_func_(local);
    WorldVector3 dv = dv_func_(local);
    
    return {du.length(), dv.length()};
}

ParamPoint3 Surface3DImpl::world_to_param(const WorldPoint3& pos) const {
    // Get basis vectors from derivatives at origin
    const auto du = du_func_(ParamPoint2(0, 0));
    const auto dv = dv_func_(ParamPoint2(0, 0));
    const auto normal = du.crossed(dv).normalized();
    const auto origin = position_func_(ParamPoint2(0, 0));
    
    // Project point onto surface plane
    const WorldVector3 rel_pos = pos - origin;
    const double normal_dist = rel_pos.dot(normal);
    const WorldVector3 planar_pos = rel_pos - normal * normal_dist;
    
    // Use Cramer's rule to solve for parameters
    const double det = du.crossed(dv).length();
    if (det < 1e-10) {
        throw std::invalid_argument("Cannot compute local coordinates: basis vectors are nearly parallel");
    }
    
    const double u = planar_pos.crossed(dv).dot(normal) / det;
    const double v = du.crossed(planar_pos).dot(normal) / det;
    
    return ParamPoint3(u, v, normal_dist);
}

std::unique_ptr<SurfacePath> Surface3DImpl::create_path(
    const GeometricPoint<2, 3, WorldSpaceTag>& start,
    const WorldVector3& world_direction,
    double world_length
) const {
    // Project direction onto tangent plane
    const auto& derivs = start.derivatives();
    const auto normal = derivs[0].crossed(derivs[1]).normalized();
    WorldVector3 planar_dir = world_direction - world_direction.dot(normal) * normal;
    
    // Normalize projected direction
    if (planar_dir.length_squared() < 1e-10) {
        throw std::invalid_argument("Direction is nearly perpendicular to surface");
    }
    planar_dir = planar_dir.normalized() * world_length;
    
    // Create path in parameter space
    const auto end_pos = start.world_pos() + planar_dir;
    const auto end_params = world_to_param(end_pos);
    
    // Create a path segment for the linear path
    auto path = std::make_unique<PathSegment>(std::const_pointer_cast<Surface3D>(shared_from_this()));
    
    // Add start and end points
    path->add_point(0.0, start.local_pos().u(), start.local_pos().v());
    path->add_point(1.0, end_params.u(), end_params.v());
    
    return path;
}

} // namespace shap
