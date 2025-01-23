#include "coord.hpp"
#pragma once
#include "param_index.hpp"
#include "param_bound.hpp"
#include <utility>

namespace shap {

/**
 * Describes a point on a surface edge in parameter space.
 * Provides information about which parameter is on a boundary,
 * which boundary it's on, and the position along that edge.
 */
struct EdgeDescriptor {
    ParamIndex param;     // Which parameter (u or v) is on boundary
    ParamBound bound;     // Which boundary (lower=0 or upper=1)
    double edge_param;    // Position along the edge [0,1]

    /**
     * Check if two edge descriptors are equal.
     * @param other Edge descriptor to compare with
     * @return true if all fields match
     */
    bool operator==(const EdgeDescriptor& other) const noexcept {
        return param == other.param &&
               bound == other.bound &&
               edge_param == other.edge_param;
    }

    /**
     * Get parameter space coordinates for a point on this edge.
     * @param t Parameter value along edge [0,1]
     * @return Pair of (u,v) coordinates
     */
    [[nodiscard]] std::pair<double, double> get_params(double t) const noexcept {
        double u = param == ParamIndex::U ? 
            (bound == ParamBound::Lower ? 0.0 : 1.0) : t;
        double v = param == ParamIndex::V ? 
            (bound == ParamBound::Lower ? 0.0 : 1.0) : t;
        return {u, v};
    }

    /**
     * Get the parameter that varies along the edge.
     * @return ParamIndex of the free parameter
     */
    [[nodiscard]] ParamIndex free_param() const noexcept {
        return param == ParamIndex::U ? ParamIndex::V : ParamIndex::U;
    }
};

} // namespace shap
