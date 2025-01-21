#include "shap/surface.hpp"
#include "shap/path.hpp"

namespace shap {

std::unique_ptr<SurfacePath> Surface::create_geodesic(
    const SurfacePoint& start,
    const SurfacePoint& end
) const {
    // Default implementation creates straight line in parameter space
    auto path = std::make_unique<TransitionPath>();
    
    path->add_segment(
        std::shared_ptr<Surface>(const_cast<Surface*>(this), [](Surface*) {}),
        0.0, 1.0,
        start.u, end.u,
        start.v, end.v,
        Vector(end.u - start.u, end.v - start.v, 0).normalize()
    );
    
    return path;
}

std::unique_ptr<SurfacePath> Surface::create_directional_path(
    const SurfacePoint& start,
    const Vector& direction,
    double length
) const {
    auto path = std::make_unique<TransitionPath>();
    
    // Convert direction to parameter space using metric tensor
    auto metric = metric_tensor(start.u, start.v);
    auto [du, dv] = metric.raise_indices(direction.x, direction.y);
    Vector param_dir(du, dv, 0);
    
    path->add_segment(
        std::shared_ptr<Surface>(const_cast<Surface*>(this), [](Surface*) {}),
        0.0, length,
        start.u, start.u + param_dir.x * length,
        start.v, start.v + param_dir.y * length,
        param_dir.normalize()
    );
    
    return path;
}

Vector Surface::parallel_transport(
    const Vector& v,
    const SurfacePath& path,
    double t_start,
    double t_end
) const {
    // Default implementation - no parallel transport
    return v;
}

} // namespace shap