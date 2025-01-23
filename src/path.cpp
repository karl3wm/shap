#include "shap/coord.hpp"
#include "shap/geometry_point2.hpp"
#include "shap/path.hpp"
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <array>
#include <iostream>

namespace shap {

namespace {
    // Constants for numerical integration
    constexpr int GEODESIC_STEPS = 50;
    constexpr double GEODESIC_DT = 1.0 / GEODESIC_STEPS;
    constexpr double CURVATURE_EPSILON = 1e-10;
    constexpr int BASE_TRANSITION_POINTS = 10;
    
    // Helper for RK4 integration
    struct RK4State {
        double u, v;      // Position
        double up, vp;    // Velocity
    };
}

void PathSegment::add_point(double t, double u, double v) {
    t_values_.push_back(t);
    u_values_.push_back(u);
    v_values_.push_back(v);
}

GeometryPoint2 PathSegment::evaluate(double t) const {
    validate_parameter(t);
    
    if (t_values_.empty()) {
        throw std::runtime_error("Path segment has no points");
    }

    // Handle exact endpoints to avoid interpolation issues
    if (t <= t_values_.front()) {
        const auto local = ParamPoint2(u_values_.front(), v_values_.front());
        return surface_->evaluate(local);
    }
    if (t >= t_values_.back()) {
        const auto local = ParamPoint2(u_values_.back(), v_values_.back());
        return surface_->evaluate(local);
    }
    
    // Find segment containing t
    auto it = std::upper_bound(t_values_.begin(), t_values_.end(), t);
    if (it == t_values_.begin() || it == t_values_.end()) {
        throw std::runtime_error("Path parameter t outside stored range");
    }
    
    const size_t segment_idx = std::distance(t_values_.begin(), it) - 1;
    
    // Linear interpolation
    const double dt = t_values_[segment_idx+1] - t_values_[segment_idx];
    const double alpha = (t - t_values_[segment_idx]) / dt;
    
    const double u = u_values_[segment_idx] + (u_values_[segment_idx+1] - u_values_[segment_idx]) * alpha;
    const double v = v_values_[segment_idx] + (v_values_[segment_idx+1] - v_values_[segment_idx]) * alpha;
    
    const auto local = ParamPoint2(u, v);
    auto geom = surface_->evaluate(local);
    
    // Diagnostic: Log evaluation details
    std::cout << "\nPath Evaluation Diagnostics:\n"
              << "t = " << t << "\n"
              << "Segment: " << segment_idx << " of " << (t_values_.size() - 1) << "\n"
              << "t range: [" << t_values_[segment_idx] << ", " << t_values_[segment_idx+1] << "]\n"
              << "alpha = " << alpha << "\n"
              << "Parameters: u=" << u << " v=" << v << "\n"
              << "Position: " << geom.world_pos().x() << ", "
              << geom.world_pos().y() << ", " << geom.world_pos().z() << "\n"
              << "Distance from start: " 
              << (geom.world_pos() - surface_->evaluate(
                     ParamPoint2(u_values_.front(), v_values_.front())
                 ).world_pos()).length() << "\n";
              
    return geom;
}

void GeodesicCurve::compute_developable_geodesic(
    const GeometryPoint2& start,
    const GeometryPoint2& end
) {
    constexpr int steps = 20;
    points_.clear();
    points_.reserve(steps + 1);
    
    const auto& start_local = start.local_pos();
    const auto& end_local = end.local_pos();
    const double du = end_local.u() - start_local.u();
    const double dv = end_local.v() - start_local.v();
    
    for (int i = 0; i <= steps; ++i) {
        const double t = static_cast<double>(i) / steps;
        const double u = start_local.u() + t * du;
        const double v = start_local.v() + t * dv;
        const auto local = ParamPoint2(u, v);
        points_.push_back(surface_->evaluate(local));
    }
}

GeometryPoint2 GeodesicCurve::evaluate(double t) const {
    validate_parameter(t);
    
    if (points_.empty()) {
        throw std::runtime_error("Geodesic curve has no points");
    }
    
    // Find segment containing t
    const auto num_segments = static_cast<double>(points_.size() - 1);
    const double scaled_t = t * num_segments;
    const size_t idx = std::min(static_cast<size_t>(scaled_t), points_.size() - 2);
    const double alpha = scaled_t - static_cast<double>(idx);
    
    // Linear interpolation
    const auto& p0 = points_[idx];
    const auto& p1 = points_[idx + 1];
    
    const auto& p0_local = p0.local_pos();
    const auto& p1_local = p1.local_pos();
    
    const double u = p0_local.u() + (p1_local.u() - p0_local.u()) * alpha;
    const double v = p0_local.v() + (p1_local.v() - p0_local.v()) * alpha;
    
    const auto local = ParamPoint2(u, v);
    return surface_->evaluate(local);
}

WorldVector3 GeodesicCurve::tangent(double t) const {
    validate_parameter(t);
    
    if (points_.size() < 2) {
        throw std::runtime_error("Geodesic curve has insufficient points for tangent computation");
    }
    
    // Find segment containing t
    const auto num_segments = static_cast<double>(points_.size() - 1);
    const double scaled_t = t * num_segments;
    const size_t idx = std::min(static_cast<size_t>(scaled_t), points_.size() - 2);
    
    // Use central difference for interior points
    if (idx > 0 && idx < points_.size() - 2) {
        const WorldVector3 diff = points_[idx+1].world_pos() - points_[idx-1].world_pos();
        return diff.normalized();
    }
    
    // Use forward/backward difference at endpoints
    const WorldVector3 diff = points_[idx+1].world_pos() - points_[idx].world_pos();
    return diff.normalized();
}

WorldVector3 GeodesicCurve::normal(double t) const {
    validate_parameter(t);
    return evaluate(t).world_normal();
}

WorldVector3 PathSegment::tangent(double t) const {
    validate_parameter(t);
    
    if (t_values_.size() < 2) {
        throw std::runtime_error("Path segment has insufficient points for tangent computation");
    }
    
    // Find segment containing t
    auto it = std::lower_bound(t_values_.begin(), t_values_.end(), t);
    const size_t idx = std::min(
        static_cast<size_t>(std::distance(t_values_.begin(), it)),
        t_values_.size() - 2
    );
    
    // Use central difference for interior points
    if (idx > 0 && idx < t_values_.size() - 2) {
        const auto p1_local = ParamPoint2(u_values_[idx+1], v_values_[idx+1]);
        const auto p0_local = ParamPoint2(u_values_[idx-1], v_values_[idx-1]);
        const auto p1 = surface_->evaluate(p1_local);
        const auto p0 = surface_->evaluate(p0_local);
        const WorldVector3 diff = p1.world_pos() - p0.world_pos();
        return diff.normalized();
    }
    
    // Use forward/backward difference at endpoints
    const auto p1_local = ParamPoint2(u_values_[idx+1], v_values_[idx+1]);
    const auto p0_local = ParamPoint2(u_values_[idx], v_values_[idx]);
    const auto p1 = surface_->evaluate(p1_local);
    const auto p0 = surface_->evaluate(p0_local);
    const WorldVector3 diff = p1.world_pos() - p0.world_pos();
    return diff.normalized();
}

WorldVector3 PathSegment::normal(double t) const {
    validate_parameter(t);
    return evaluate(t).world_normal();
}

void TransitionPath::add_segment(
    std::shared_ptr<Surface> surface,
    double t_start, double t_end,
    double u_start, double u_end,
    double v_start, double v_end,
    const WorldVector3& /*direction*/  // Used by derived classes
) {
    if (!surface) {
        throw std::invalid_argument("Surface pointer cannot be null");
    }

    auto segment = std::make_unique<PathSegment>(
        std::shared_ptr<Surface>(const_cast<Surface*>(surface.get()), [](Surface*){})
    );
    
    // Adaptive sampling based on surface curvature
    int num_points = BASE_TRANSITION_POINTS;
    
    // Get surface properties at start
    const auto start_local = ParamPoint2(u_start, v_start);
    const auto geom = surface->evaluate(start_local);
    if (geom.gaussian_curvature()) {
        const double curvature = std::abs(*geom.gaussian_curvature());
        num_points += static_cast<int>(5.0 * std::sqrt(curvature));
    }
    
    // Pre-compute parameter deltas
    const double dt = t_end - t_start;
    const double du = u_end - u_start;
    const double dv = v_end - v_start;
    
    // Linear interpolation for transition paths
    for (int i = 0; i < num_points; ++i) {
        const double alpha = static_cast<double>(i) / (num_points - 1);
        segment->add_point(
            t_start + dt * alpha,
            u_start + du * alpha,
            v_start + dv * alpha
        );
    }
    
    segments_.push_back(std::move(segment));
}

GeometryPoint2 TransitionPath::evaluate(double t) const {
    validate_parameter(t);
    
    if (segments_.empty()) {
        throw std::runtime_error("Transition path has no segments");
    }
    
    // Find segment containing t
    for (const auto& segment : segments_) {
        if (t <= segment->t_values().back()) {
            return segment->evaluate(t);
        }
    }
    
    // If t is beyond last segment, evaluate at end of last segment
    return segments_.back()->evaluate(segments_.back()->t_values().back());
}

WorldVector3 TransitionPath::tangent(double t) const {
    validate_parameter(t);
    
    if (segments_.empty()) {
        throw std::runtime_error("Transition path has no segments");
    }
    
    // Find segment containing t
    for (const auto& segment : segments_) {
        if (t <= segment->t_values().back()) {
            return segment->tangent(t);
        }
    }
    
    // If t is beyond last segment, use tangent at end of last segment
    return segments_.back()->tangent(segments_.back()->t_values().back());
}

WorldVector3 TransitionPath::normal(double t) const {
    validate_parameter(t);
    
    if (segments_.empty()) {
        throw std::runtime_error("Transition path has no segments");
    }
    
    // Find segment containing t
    for (const auto& segment : segments_) {
        if (t <= segment->t_values().back()) {
            return segment->normal(t);
        }
    }
    
    // If t is beyond last segment, use normal at end of last segment
    return segments_.back()->normal(segments_.back()->t_values().back());
}

} // namespace shap
