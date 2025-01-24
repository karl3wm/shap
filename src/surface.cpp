#include "shap/surface.hpp"
#include "shap/path.hpp"
#include <memory>

namespace shap {

std::unique_ptr<SurfacePath> Surface::create_path(
    const GeometryPoint2& start,
    const WorldVector3& world_direction,
    double world_length
) const {
    // Convert GeometryPoint2 to GeometricPoint<2,3,WorldSpaceTag>
    auto start3d = impl_->evaluate(start.local_pos());

    // Create path using Surface3D implementation
    auto path = std::make_unique<PathSegment>(impl_);
    
    // Forward to Surface3D implementation
    auto surface3d_path = impl_->create_path(start3d, world_direction, world_length);
    
    // Cast Surface3D path to PathSegment to access point data
    auto* segment = dynamic_cast<PathSegment*>(surface3d_path.get());
    if (!segment) {
        throw std::runtime_error("Surface3D path must be a PathSegment");
    }
    
    // Copy points from Surface3D path to legacy path
    const auto& t_values = segment->t_values();
    const auto& u_values = segment->u_values();
    const auto& v_values = segment->v_values();
    
    for (size_t i = 0; i < t_values.size(); ++i) {
        path->add_point(t_values[i], u_values[i], v_values[i]);
    }
    
    return path;
}

} // namespace shap
