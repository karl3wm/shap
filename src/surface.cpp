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
    
    // Delegate to Surface3D implementation
    return impl_->create_path(start3d, world_direction, world_length);
}

} // namespace shap
