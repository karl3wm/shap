#include "shap/coord.hpp"
#include "shap/surface_collection.hpp"
#include "shap/manifold.hpp"
#include "shap/riemannian_metric.hpp"
#include <algorithm>
#include <iostream>

namespace shap {

GeometricPoint<2, 3, WorldSpaceTag> SurfaceConnection::map_point(const GeometricPoint<2, 3, WorldSpaceTag>& point) const {
    // Get edge descriptor for source point
    const auto edge_desc = point.get_edge_descriptor();
    if (!edge_desc) {
        throw std::invalid_argument("Point is not on an edge");
    }
    
    // Map parameter along edge using orientation
    const double edge_param = edge_desc->edge_param;
    const double mapped_param = orientation_ > 0 ? edge_param : 1.0 - edge_param;
    
    // Convert to target surface parameters
    double u, v;
    if (target_edge_.param == ParamIndex::U) {
        u = target_edge_.bound == ParamBound::Lower ? 0.0 : 1.0;
        v = mapped_param;
    } else {
        u = mapped_param;
        v = target_edge_.bound == ParamBound::Lower ? 0.0 : 1.0;
    }
    
    // Create point on target surface
    const auto local = ParamPoint2(u, v);
    auto target_point = target_->evaluate(local);
    
    // Check normal orientation
    if (target_point.world_normal().dot(point.world_normal()) < 0) {
        throw std::runtime_error("Surface normals have opposite orientation");
    }
    
    return target_point;
}

std::unique_ptr<SurfacePath> SurfaceCollection::create_path(
    const GeometricPoint<2, 3, WorldSpaceTag>& start,
    const WorldVector3& world_direction,
    double world_length
) const {
    if (world_length <= 0) {
        throw std::invalid_argument("Path length must be positive");
    }
    if (world_direction.length_squared() < 1e-10) {
        throw std::invalid_argument("Direction vector cannot be zero");
    }

    // Create transition path to handle multiple surfaces
    auto path = std::make_shared<TransitionPath>();
    double t = 0.0;
    auto current = start;
    auto current_dir = world_direction;
    
    while (t < world_length) {
        // Get current surface
        auto current_surface = current.surface();
        if (!current_surface) {
            throw std::runtime_error("Invalid surface pointer");
        }
        
        // Try path solver first for surface transitions
        if (auto solver = current_surface->get_path_solver()) {
            auto intersection = (*solver)(current.world_pos(), current_dir, world_length - t);
            if (intersection) {
                // Convert end point to parameter space
                const auto end_local = current_surface->nearest(intersection->position);
                const auto& start_local = current.local_pos();
                
                // Add segment up to intersection
                path->add_segment(
                    std::const_pointer_cast<Surface3D>(
                        std::dynamic_pointer_cast<const Surface3D>(current_surface->shared_from_this())
                    ),
                    ParamPoint1(t), ParamPoint1(t + intersection->t),
                    start_local, end_local,
                    current_dir
                );
                
                // Find connection at intersection point
                const auto trans_point = current_surface->evaluate(end_local);
                auto connection = find_connection(trans_point);
                if (!connection) {
                    // End of path at surface boundary
                    break;
                }
                
                // Map point to next surface
                current = connection->map_point(trans_point);
                t += intersection->t;
                
                // Update direction for next surface
                const auto new_geom = current_surface->evaluate(current.local_pos());
                current_dir = current_dir - current_dir.dot(new_geom.world_normal()) * new_geom.world_normal();
                if (current_dir.length_squared() < 1e-10) {
                    throw std::runtime_error("Direction became perpendicular to surface");
                }
                current_dir = current_dir.normalized();
                continue;
            }
        }
        
        // No intersection found, continue to end of path
        const double remaining = world_length - t;
        const auto& start_local = current.local_pos();
        
        // Convert direction to parameter space using metric tensor
        const auto geom = current_surface->evaluate(start_local);
        const auto& derivs = geom.derivatives();
        const auto normal = derivs[0].crossed(derivs[1]).normalized();
        
        // Create metric tensor at current point
        RiemannianMetric metric(
            derivs[0].dot(derivs[0]),  // g11
            derivs[0].dot(derivs[1]),  // g12
            derivs[1].dot(derivs[1]),  // g22
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0  // Derivatives not needed for pullback
        );
        
        // Use pullback to convert world velocity to parameter space
        const auto param_vec = metric.pullback(current_dir, derivs[0], derivs[1], normal);
        const auto param_vel = ParamVector2(param_vec.u(), param_vec.v());
        
        // Scale parameter derivatives by inverse of surface scale factors
        const auto [du_scale, dv_scale] = current_surface->get_scale_factors(start_local);
        const double scaled_du = param_vel.u() / (du_scale > 1e-10 ? du_scale : 1.0);
        const double scaled_dv = param_vel.v() / (dv_scale > 1e-10 ? dv_scale : 1.0);
        
        // Compute end parameters
        const auto end_local = ParamPoint2(
            start_local.u() + scaled_du * remaining,
            start_local.v() + scaled_dv * remaining
        );
        
        // Add final segment
        path->add_segment(
            std::const_pointer_cast<Surface3D>(
                std::dynamic_pointer_cast<const Surface3D>(current_surface->shared_from_this())
            ),
            ParamPoint1(t), ParamPoint1(world_length),
            start_local, end_local,
            current_dir
        );
        break;
    }
    
    // Create path with evaluation functions from the transition path
    return std::make_unique<SurfacePath>(
        [p = path](const ParamPoint1& param) { return p->evaluate(param).world_pos(); },
        [p = path](const ParamPoint1& param) { return p->derivatives(param)[0]; },
        [p = path](const ParamPoint1& param) { 
            const auto derivs = p->derivatives(param);
            return derivs[0].crossed(derivs[1]).normalized();
        }
    );
}

} // namespace shap
