#include "coord.hpp"
#pragma once
#include "edge_descriptor.hpp"

namespace shap {

/**
 * Describes a connection between edges of two surfaces.
 * Provides information about how the edges align and how to map
 * parameters between them.
 */
struct EdgeConnection {
    EdgeDescriptor edge1;  // Edge descriptor for first surface
    EdgeDescriptor edge2;  // Edge descriptor for second surface
    int orientation;       // +1 if parameters flow in same direction, -1 if opposite

    /**
     * Map a parameter value from edge1 to edge2.
     * Takes into account the orientation of the connection.
     * 
     * @param param Parameter value on edge1 [0,1]
     * @return Corresponding parameter value on edge2 [0,1]
     */
    [[nodiscard]] double map_parameter(double param) const noexcept {
        return orientation > 0 ? param : 1.0 - param;
    }

    /**
     * Check if two edge connections are equal.
     * @param other Edge connection to compare with
     * @return true if all fields match
     */
    bool operator==(const EdgeConnection& other) const noexcept {
        return edge1 == other.edge1 &&
               edge2 == other.edge2 &&
               orientation == other.orientation;
    }
};

} // namespace shap
