#include "coord.hpp"
#pragma once
#include <cstdint>

namespace shap {

/**
 * Classification of surface types based on their geometric properties.
 * Used to optimize algorithms and handle special cases appropriately.
 */
enum class SurfaceType : uint8_t {
    Generic,     // Base type for surfaces with no special properties
    Smooth,      // Surface with no singularities or edges
    Developable, // Surface with zero Gaussian curvature (can be flattened)
    Singular     // Surface containing singularities or edges
};

} // namespace shap
