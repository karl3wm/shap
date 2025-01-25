#pragma once
#include "coord.hpp"
#include "edge_descriptor.hpp"
#include "geometric_point.hpp"
#include "path.hpp"
#include "surface3d.hpp"
#include <vector>
#include <memory>

namespace shap {

/**
 * Represents a connection between two surfaces along their edges.
 * Handles mapping points from one surface to another across the connection.
 */
class SurfaceConnection {
public:
    SurfaceConnection(
        Surface3D* target_surface,
        EdgeDescriptor target_edge_desc,
        int orientation_sign
    ) : target_(target_surface)
      , target_edge_(target_edge_desc)
      , orientation_(orientation_sign) {}

    // Map a point from source surface to target surface
    [[nodiscard]] GeometricPoint<2, 3, WorldSpaceTag> map_point(const GeometricPoint<2, 3, WorldSpaceTag>& point) const;

private:
    Surface3D* target_;              // Target surface for connection
    EdgeDescriptor target_edge_;   // Edge descriptor on target surface
    int orientation_;             // +1 if parameters map directly, -1 if reversed
};

/**
 * Represents a collection of connected surfaces.
 * Manages transitions between surfaces and path creation across multiple surfaces.
 */
class SurfaceCollection {
public:
    // Add a surface to the collection
    void add_surface(std::shared_ptr<Surface3D> surface) {
        surfaces_.push_back(std::move(surface));
    }

    // Get a surface by index
    [[nodiscard]] Surface3D* get_surface(size_t index) const {
        if (index >= surfaces_.size()) return nullptr;
        return surfaces_[index].get();
    }

    // Add a connection between surfaces
    void add_connection(
        Surface3D* source,
        EdgeDescriptor source_edge,
        Surface3D* target,
        EdgeDescriptor target_edge,
        int orientation
    ) {
        connections_.emplace_back(
            source,
            source_edge,
            std::make_unique<SurfaceConnection>(target, target_edge, orientation)
        );
    }

    // Create a path across multiple surfaces
    [[nodiscard]] std::unique_ptr<SurfacePath> create_path(
        const GeometricPoint<2, 3, WorldSpaceTag>& start,
        const WorldVector3& world_direction,
        double world_length
    ) const;

protected:
    // Find connection for a point on a surface edge
    [[nodiscard]] const SurfaceConnection* find_connection(const GeometricPoint<2, 3, WorldSpaceTag>& point) const {
        if (!point.is_on_edge()) return nullptr;
        
        const auto edge_desc = point.get_edge_descriptor();
        if (!edge_desc) return nullptr;
        
        // Find matching connection
        for (const auto& conn : connections_) {
            if (conn.source == point.surface() && 
                conn.source_edge.param == edge_desc->param &&
                conn.source_edge.bound == edge_desc->bound) {
                return conn.connection.get();
            }
        }
        return nullptr;
    }

private:
    // Connection between two surfaces
    struct Connection {
        Surface3D* source;
        EdgeDescriptor source_edge;
        std::unique_ptr<SurfaceConnection> connection;

        Connection(
            Surface3D* src,
            EdgeDescriptor src_edge,
            std::unique_ptr<SurfaceConnection> conn
        ) : source(src)
          , source_edge(src_edge)
          , connection(std::move(conn)) {}
    };

    std::vector<std::shared_ptr<Surface3D>> surfaces_;
    std::vector<Connection> connections_;
};

} // namespace shap
