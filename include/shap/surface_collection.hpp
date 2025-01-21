#pragma once
#include "surface.hpp"
#include "path.hpp"
#include <functional>
#include <unordered_map>
#include <vector>

namespace shap {

// Connection type between surfaces
enum class ConnectionType {
    Linear,    // Simple linear transition
    Circular,  // Smooth circular arc transition
    Geodesic   // Follow geodesic curves (for smooth surfaces)
};

// Edge specification for connections
enum class Edge {
    Left,
    Right,
    Top,
    Bottom
};

// Collection of connected surfaces with simplified connection handling
class SurfaceCollection {
public:
    // Add surface to collection with fluent interface
    SurfaceCollection& add(std::shared_ptr<Surface> surface, const std::string& name) {
        if (!name.empty()) {
            surface->name = name;
            surface_map[name] = surface;
        }
        surfaces.push_back(surface);
        return *this;
    }
    
    // Get surface by name
    std::shared_ptr<Surface> get(const std::string& name) const {
        auto it = surface_map.find(name);
        if (it == surface_map.end()) {
            throw std::runtime_error("Surface not found: " + name);
        }
        return it->second;
    }
    
    // Connection builder for fluent interface
    class ConnectionBuilder {
    public:
        ConnectionBuilder& along(Edge edge) {
            edge_ = edge;
            return *this;
        }
        
        ConnectionBuilder& with_type(ConnectionType type) {
            type_ = type;
            return *this;
        }
        
        void build() {
            if (!collection_ || !surface1_ || !surface2_) {
                throw std::runtime_error("Invalid connection specification");
            }
            
            // Create connection based on edge and type
            collection_->create_connection(
                surface1_, surface2_,
                create_transition_test(edge_),
                create_point_mapping(edge_, type_)
            );
        }
        
    private:
        friend class SurfaceCollection;
        ConnectionBuilder(
            SurfaceCollection* collection,
            std::shared_ptr<Surface> s1,
            std::shared_ptr<Surface> s2
        ) : collection_(collection), surface1_(s1), surface2_(s2) {}
        
        SurfaceCollection* collection_;
        std::shared_ptr<Surface> surface1_;
        std::shared_ptr<Surface> surface2_;
        Edge edge_ = Edge::Right;
        ConnectionType type_ = ConnectionType::Linear;
        
        // Helper to create transition test based on edge
        static std::function<bool(const SurfacePoint&, const Vector&)>
        create_transition_test(Edge edge) {
            switch (edge) {
                case Edge::Right:
                    return [](const SurfacePoint& pt, const Vector& dir) {
                        return pt.u >= 0.95 && dir.x > 0;
                    };
                case Edge::Left:
                    return [](const SurfacePoint& pt, const Vector& dir) {
                        return pt.u <= 0.05 && dir.x < 0;
                    };
                // Add cases for other edges
                default:
                    throw std::runtime_error("Unsupported edge type");
            }
        }
        
        // Helper to create point mapping based on edge and transition type
        static std::function<SurfacePoint(const SurfacePoint&)>
        create_point_mapping(Edge edge, ConnectionType type) {
            // Basic mapping for linear transition
            switch (edge) {
                case Edge::Right:
                    return [](const SurfacePoint& pt) {
                        return SurfacePoint(
                            "next",             // Updated by collection
                            0.0, pt.v,          // Map to left edge
                            pt.position,        // Keep position
                            Vector(1, 0, 0),    // Update normal
                            Vector(0, -1, 0),   // Update du
                            Vector(0, 0, 1)     // Keep dv
                        );
                    };
                // Add cases for other edges
                default:
                    throw std::runtime_error("Unsupported edge type");
            }
        }
    };
    
    // Start connection specification with fluent interface
    ConnectionBuilder connect(const std::string& surface1, const std::string& surface2) {
        auto s1 = get(surface1);
        auto s2 = get(surface2);
        return ConnectionBuilder(this, s1, s2);
    }
    
    // Create path that can transition between surfaces
    std::unique_ptr<SurfacePath> create_path(
        const SurfacePoint& start,
        const Vector& direction,
        double length
    ) const {
        auto path = std::make_unique<TransitionPath>();
        
        double t = 0.0;
        SurfacePoint current = start;
        Vector current_dir = direction;
        
        while (t < length) {
            // Find current surface and next transition
            auto current_surface = get(current.surface_name);
            const Connection* next_connection = nullptr;
            double dist_to_transition = std::numeric_limits<double>::max();
            
            // Check for transitions
            for (const auto& conn : connections) {
                if (conn.surface1 == current_surface && 
                    conn.transition_test(current, current_dir)) {
                    next_connection = &conn;
                    // Simple approximation of distance to edge
                    dist_to_transition = 1.0 - current.u;
                    break;
                }
            }
            
            if (next_connection) {
                // Add segment up to transition
                double segment_length = std::min(0.25, dist_to_transition);
                path->add_segment(
                    current_surface, t, t + segment_length,
                    current.u, std::min(1.0, current.u + current_dir.x * segment_length),
                    current.v, current.v + current_dir.y * segment_length,
                    current_dir
                );
                
                // Transition to next surface
                current = next_connection->map_point(current);
                current.surface_name = next_connection->surface2->name;
                t += segment_length;
            } else {
                // Add segment on current surface
                double segment_length = std::min(0.25, length - t);
                path->add_segment(
                    current_surface, t, t + segment_length,
                    current.u, std::min(1.0, current.u + current_dir.x * segment_length),
                    current.v, current.v + current_dir.y * segment_length,
                    current_dir
                );
                
                // Update current point
                current.u = std::min(1.0, current.u + current_dir.x * segment_length);
                current.v += current_dir.y * segment_length;
                current = current_surface->evaluate(current.u, current.v);
                t += segment_length;
            }
        }
        
        return path;
    }

private:
    struct Connection {
        std::shared_ptr<Surface> surface1;
        std::shared_ptr<Surface> surface2;
        std::function<bool(const SurfacePoint&, const Vector&)> transition_test;
        std::function<SurfacePoint(const SurfacePoint&)> map_point;
    };
    
    void create_connection(
        std::shared_ptr<Surface> s1,
        std::shared_ptr<Surface> s2,
        std::function<bool(const SurfacePoint&, const Vector&)> test,
        std::function<SurfacePoint(const SurfacePoint&)> map
    ) {
        connections.push_back({s1, s2, test, map});
    }
    
    std::vector<std::shared_ptr<Surface>> surfaces;
    std::vector<Connection> connections;
    std::unordered_map<std::string, std::shared_ptr<Surface>> surface_map;
};

} // namespace shap