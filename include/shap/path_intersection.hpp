#pragma once
#include "coord.hpp"
#include "param_bound.hpp"
#include "param_index.hpp"
#include <utility>
#include <stdexcept>

namespace shap {

// Path solver returns intersection with surface boundary
struct PathIntersection {
    double t;                // Distance to intersection in world space
    WorldPoint3 position;    // World space intersection point
    ParamIndex param;        // Which parameter (u/v) hit boundary
    ParamBound bound;        // Which bound (0/1) was hit
    double edge_parameter;   // Parameter along the edge [0,1]

    // Constructor with validation
    PathIntersection(
        double t_,
        WorldPoint3 position_,
        ParamIndex param_,
        ParamBound bound_,
        double edge_parameter_
    ) : t(t_)
      , position(std::move(position_))
      , param(param_)
      , bound(bound_)
      , edge_parameter(edge_parameter_) {
        if (t_ < 0) {
            throw std::invalid_argument("Intersection distance must be non-negative");
        }
        if (edge_parameter_ < 0 || edge_parameter_ > 1) {
            throw std::invalid_argument("Edge parameter must be in [0,1]");
        }
    }
};

} // namespace shap
