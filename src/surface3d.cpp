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
Surface3D::get_scale_factors(const ParameterPoint& local) const {
    validate_parameters(local);
    
    // Scale factors are lengths of first derivatives
    WorldVector3 du = du_func_(local);
    WorldVector3 dv = dv_func_(local);
    
    return {du.length(), dv.length()};
}

std::unique_ptr<SurfacePath> Surface3D::create_path(
    const GeometricPoint<2, 3, WorldSpaceTag>& start,
    const WorldVector3& world_direction,
    double world_length
) const {
    // Project direction onto tangent plane and normalize
    const auto& derivs = start.derivatives();
    const auto normal = derivs[0].crossed(derivs[1]).normalized();
    WorldVector3 planar_dir = world_direction - world_direction.dot(normal) * normal;
    
    // Validate direction is not perpendicular
    if (planar_dir.length_squared() < 1e-10) {
        throw std::invalid_argument("Direction is nearly perpendicular to surface");
    }
    
    // Scale direction to preserve world length after projection
    const double scale = world_length / planar_dir.length();
    planar_dir *= scale;
    
    // Create path in parameter space
    const auto end_pos = start.world_pos() + planar_dir;
    const auto end_params = nearest(end_pos);
    
    // Create a path segment for the linear path
    auto segment = std::make_shared<PathSegment>(
        std::const_pointer_cast<Surface3D>(
            std::dynamic_pointer_cast<const Surface3D>(shared_from_this())));
    
    // Add start and end points
    segment->add_point(ParamPoint1(0.0), ParamPoint2(start.local_pos().u(), start.local_pos().v()));
    segment->add_point(ParamPoint1(1.0), end_params);
    
    // Create path using segment's evaluation functions
    return std::make_unique<SurfacePath>(
        std::bind(&PathSegment::evaluate_position, segment, std::placeholders::_1),
        std::bind(&PathSegment::evaluate_tangent, segment, std::placeholders::_1),
        std::bind(&PathSegment::evaluate_normal, segment, std::placeholders::_1)
    );
}

} // namespace shap
