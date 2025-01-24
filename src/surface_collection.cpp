#include "shap/coord.hpp"
#include "shap/surface_collection.hpp"
#include <algorithm>
#include <iostream>

namespace shap {

GeometryPoint2 SurfaceConnection::map_point(const GeometryPoint2& point) const {
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
    const GeometryPoint2& start,
    const WorldVector3& world_direction,
    double world_length
) const {
    if (world_length <= 0) {
        throw std::invalid_argument("Path length must be positive");
    }
    if (world_direction.length_squared() < 1e-10) {
        throw std::invalid_argument("Direction vector cannot be zero");
    }

    auto path = std::make_unique<TransitionPath>();
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
                const auto end_local = current_surface->world_to_param(intersection->position).uv();
                const auto& start_local = current.local_pos();
                
                // Add segment up to intersection
                path->add_segment(
                    current_surface->impl(),
                    t, t + intersection->t,
                    start_local.u(), end_local.u(),
                    start_local.v(), end_local.v(),
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
        
        // Convert direction to parameter space and scale by surface metric
        const auto geom = current_surface->evaluate(start_local);
        const auto param_vel = current_surface->world_to_parameter_velocity(
            current_dir, geom.world_du(), geom.world_dv());
        
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
            current_surface->impl(),
            t, world_length,
            start_local.u(), end_local.u(),
            start_local.v(), end_local.v(),
            current_dir
        );
        break;
    }
    
    return path;
}

} // namespace shap
