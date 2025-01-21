#pragma once
#include "surface.hpp"
#include "path.hpp"
#include <functional>
#include <unordered_map>
#include <vector>

namespace shap {

// Collection of connected surfaces with transition handling
class SurfaceCollection {
public:
    // Add surface with optional name
    void add_surface(std::shared_ptr<Surface> surface, const std::string& name = "") {
        if (!name.empty()) {
            surface->name = name;
            surface_map[name] = surface;
        }
        surfaces.push_back(surface);
    }
    
    // Get surface by name
    std::shared_ptr<Surface> get_surface(const std::string& name) const {
        auto it = surface_map.find(name);
        if (it == surface_map.end()) {
            throw std::runtime_error("Surface not found: " + name);
        }
        return it->second;
    }
    
    // Define connection between surfaces
    struct Connection {
        std::shared_ptr<Surface> surface1;
        std::shared_ptr<Surface> surface2;
        std::function<bool(const SurfacePoint&, const Vector&)> transition_test;
        std::function<SurfacePoint(const SurfacePoint&)> map_point;
    };
    
    void add_connection(
        std::shared_ptr<Surface> s1,
        std::shared_ptr<Surface> s2,
        std::function<bool(const SurfacePoint&, const Vector&)> test,
        std::function<SurfacePoint(const SurfacePoint&)> map
    ) {
        connections.push_back({s1, s2, test, map});
    }
    
    // Find surface and next transition info for a point
    struct SurfaceInfo {
        std::shared_ptr<Surface> surface;
        const Connection* next_connection;
        double distance_to_transition;
        
        SurfaceInfo(
            std::shared_ptr<Surface> s,
            const Connection* conn = nullptr,
            double dist = std::numeric_limits<double>::max()
        ) : surface(s), next_connection(conn), distance_to_transition(dist) {}
    };
    
    SurfaceInfo find_surface_info(const SurfacePoint& pt, const Vector& dir) const {
        // Find current surface by name
        auto current_surface = get_surface(pt.surface_name);
        
        // Find next transition
        for (const auto& conn : connections) {
            if (conn.surface1 == current_surface && conn.transition_test(pt, dir)) {
                // Calculate approximate distance to transition
                double dist_to_edge = 1.0 - pt.u; // Simple approximation
                return SurfaceInfo(current_surface, &conn, dist_to_edge);
            }
        }
        
        return SurfaceInfo(current_surface);
    }
    
    // Create path that can transition between surfaces
    std::unique_ptr<SurfacePath> create_path(
        const SurfacePoint& start,
        const Vector& direction,
        double length
    ) {
        auto path = std::make_unique<TransitionPath>();
        
        double t = 0.0;
        SurfacePoint current = start;
        Vector current_dir = direction;
        
        while (t < length) {
            // Find current surface and check for transition
            auto info = find_surface_info(current, current_dir);
            
            if (info.next_connection) {
                // Add segment up to transition
                double segment_length = std::min(0.25, info.distance_to_transition);
                path->add_segment(
                    info.surface, t, t + segment_length,
                    current.u, std::min(1.0, current.u + current_dir.x * segment_length),
                    current.v, current.v + current_dir.y * segment_length,
                    current_dir
                );
                
                // Transition to next surface
                current = info.next_connection->map_point(current);
                t += segment_length;
            } else {
                // Add segment on current surface
                double segment_length = std::min(0.25, length - t);
                path->add_segment(
                    info.surface, t, t + segment_length,
                    current.u, std::min(1.0, current.u + current_dir.x * segment_length),
                    current.v, current.v + current_dir.y * segment_length,
                    current_dir
                );
                
                // Update current point
                current.u = std::min(1.0, current.u + current_dir.x * segment_length);
                current.v += current_dir.y * segment_length;
                current = info.surface->evaluate(current.u, current.v);
                t += segment_length;
            }
        }
        
        return path;
    }

private:
    std::vector<std::shared_ptr<Surface>> surfaces;
    std::vector<Connection> connections;
    std::unordered_map<std::string, std::shared_ptr<Surface>> surface_map;
};

} // namespace shap