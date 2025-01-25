#include "shap/path.hpp"
#include "shap/coord.hpp"
#include "shap/geometric_point.hpp"
#include "shap/surface3d.hpp"
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <array>
#include <iostream>

namespace shap {

namespace {
    // Constants for SLERP
    constexpr double SLERP_EPSILON = 1e-10;
}

GeodesicCurve::GeodesicCurve(
    std::shared_ptr<const Surface3D> surface,
    const GeometricPoint<2, 3, WorldSpaceTag>& start,
    const GeometricPoint<2, 3, WorldSpaceTag>& end
) : Path3D(
        std::bind(&GeodesicCurve::evaluate_position, this, std::placeholders::_1),
        std::bind(&GeodesicCurve::evaluate_tangent, this, std::placeholders::_1),
        std::bind(&GeodesicCurve::evaluate_normal, this, std::placeholders::_1)
    )
    , surface_(std::move(surface))
{
    compute_developable_geodesic(start, end);
}

WorldPoint3 GeodesicCurve::evaluate_position(const ParamPoint1& param) const {
    if (points_.empty()) {
        throw std::runtime_error("Geodesic curve has no points");
    }
    const auto num_segments = static_cast<double>(points_.size() - 1);
    const double scaled_t = static_cast<double>(param) * num_segments;
    const size_t idx = std::min(static_cast<size_t>(scaled_t), points_.size() - 2);
    const double alpha = scaled_t - static_cast<double>(idx);
    const auto& p0 = points_[idx];
    const auto& p1 = points_[idx + 1];
    const auto& p0_local = p0.local_pos();
    const auto& p1_local = p1.local_pos();
    const double u = p0_local[0] + (p1_local[0] - p0_local[0]) * alpha;
    const double v = p0_local[1] + (p1_local[1] - p0_local[1]) * alpha;
    const auto local = ParamPoint2(u, v);
    return surface_->evaluate(local).world_pos();
}

WorldVector3 GeodesicCurve::evaluate_tangent(const ParamPoint1& param) const {
    if (points_.size() < 2) {
        throw std::runtime_error("Geodesic curve has insufficient points for tangent computation");
    }
    const auto num_segments = static_cast<double>(points_.size() - 1);
    const double scaled_t = static_cast<double>(param) * num_segments;
    const size_t idx = std::min(static_cast<size_t>(scaled_t), points_.size() - 2);
    if (idx > 0 && idx < points_.size() - 2) {
        const WorldVector3 diff = points_[idx+1].world_pos() - points_[idx-1].world_pos();
        return diff.normalized();
    }
    const WorldVector3 diff = points_[idx+1].world_pos() - points_[idx].world_pos();
    return diff.normalized();
}

WorldVector3 GeodesicCurve::evaluate_normal(const ParamPoint1& param) const {
    // Use param to get appropriate point
    const auto num_segments = static_cast<double>(points_.size() - 1);
    const double scaled_t = static_cast<double>(param) * num_segments;
    const size_t idx = std::min(static_cast<size_t>(scaled_t), points_.size() - 1);
    const auto local = points_[idx].local_pos();
    const auto geom = surface_->evaluate(local);
    const auto& derivs = geom.derivatives();
    return derivs[0].crossed(derivs[1]).normalized();
}

void GeodesicCurve::compute_developable_geodesic(
    const GeometricPoint<2, 3, WorldSpaceTag>& start,
    const GeometricPoint<2, 3, WorldSpaceTag>& end
) {
    switch (surface_->surface_type()) {
        case SurfaceType::Smooth: {
            // For spheres (only smooth surface currently supported)
            compute_smooth_geodesic(start, end);
            break;
        }
        case SurfaceType::Developable: {
            // For flat patches and cube faces
            points_.clear();
            points_.reserve(2);  // Only need start and end for exact line
            points_.push_back(start);
            points_.push_back(end);
            break;
        }
        default:
            throw std::runtime_error("Surface type does not support exact geodesics");
    }
}

void GeodesicCurve::compute_smooth_geodesic(
    const GeometricPoint<2, 3, WorldSpaceTag>& start,
    const GeometricPoint<2, 3, WorldSpaceTag>& end
) {
    // For spheres, compute great circle arc
    const WorldPoint3 origin(0, 0, 0);
    const WorldVector3 start_vec = start.world_pos() - origin;
    const WorldVector3 end_vec = end.world_pos() - origin;
    const double radius = start_vec.length();
    
    // Get normalized vectors
    const WorldVector3 v1 = start_vec.normalized();
    const WorldVector3 v2 = end_vec.normalized();
    
    // Compute angle between vectors
    const double cos_angle = std::clamp(v1.dot(v2), -1.0, 1.0);
    const double angle = std::acos(cos_angle);
    
    // Sample points along great circle
    constexpr int steps = 20;  // Can be fewer points since this is exact
    points_.clear();
    points_.reserve(steps + 1);
    
    for (int i = 0; i <= steps; ++i) {
        const double t = static_cast<double>(i) / steps;
        const double current_angle = angle * t;
        const double sin_angle = std::sin(angle);
        
        // Spherical linear interpolation (SLERP)
        WorldPoint3 interp_point;
        if (sin_angle > SLERP_EPSILON) {
            const double sin_current = std::sin(current_angle);
            const double sin_remaining = std::sin(angle - current_angle);
            const WorldVector3 interp_vec = 
                (v1 * sin_remaining + v2 * sin_current) * (radius / sin_angle);
            interp_point = origin + interp_vec;
        } else {
            // Points are very close, use linear interpolation
            interp_point = WorldPoint3(
                start.world_pos().x() * (1.0 - t) + end.world_pos().x() * t,
                start.world_pos().y() * (1.0 - t) + end.world_pos().y() * t,
                start.world_pos().z() * (1.0 - t) + end.world_pos().z() * t
            );
        }
        const auto param_point = surface_->nearest(interp_point);
        points_.push_back(surface_->evaluate(param_point));
    }
}

void PathSegment::add_point(const ParamPoint1& t, const ParamPoint2& uv) {
    t_values_.push_back(t);
    uv_values_.push_back(uv);
}

PathSegment::PathSegment(std::shared_ptr<const Surface3D> surface)
    : Path3D(
        std::bind(&PathSegment::evaluate_position, this, std::placeholders::_1),
        std::bind(&PathSegment::evaluate_tangent, this, std::placeholders::_1),
        std::bind(&PathSegment::evaluate_normal, this, std::placeholders::_1)
    )
    , surface_(std::move(surface))
{
    // Pre-allocate space for typical path size
    t_values_.reserve(100);
    uv_values_.reserve(100);
}

WorldPoint3 PathSegment::evaluate_position(const ParamPoint1& param) const {
    if (t_values_.empty()) {
        throw std::runtime_error("Path segment has no points");
    }
    if (param <= t_values_.front()) {
        const auto local = uv_values_.front();
        return surface_->evaluate(local).world_pos();
    }
    if (param >= t_values_.back()) {
        const auto local = uv_values_.back();
        return surface_->evaluate(local).world_pos();
    }
    auto it = std::upper_bound(t_values_.begin(), t_values_.end(), param);
    if (it == t_values_.begin() || it == t_values_.end()) {
        throw std::runtime_error("Path parameter t outside stored range");
    }
    const size_t segment_idx = std::distance(t_values_.begin(), it) - 1;
    const double dt = static_cast<double>(t_values_[segment_idx+1]) - static_cast<double>(t_values_[segment_idx]);
    const double alpha = (static_cast<double>(param) - static_cast<double>(t_values_[segment_idx])) / dt;
    const auto& uv0 = uv_values_[segment_idx];
    const auto& uv1 = uv_values_[segment_idx+1];
    const auto local = ParamPoint2(
        uv0.u() + (uv1.u() - uv0.u()) * alpha,
        uv0.v() + (uv1.v() - uv0.v()) * alpha
    );
    auto geom = surface_->evaluate(local);
    return geom.world_pos();
}

WorldVector3 PathSegment::evaluate_tangent(const ParamPoint1& param) const {
    if (t_values_.size() < 2) {
        throw std::runtime_error("Path segment has insufficient points for tangent computation");
    }
    auto it = std::lower_bound(t_values_.begin(), t_values_.end(), param);
    const size_t idx = std::min(
        static_cast<size_t>(std::distance(t_values_.begin(), it)),
        t_values_.size() - 2
    );
    if (idx > 0 && idx < t_values_.size() - 2) {
        const auto& p1_local = uv_values_[idx+1];
        const auto& p0_local = uv_values_[idx-1];
        const auto p1 = surface_->evaluate(p1_local);
        const auto p0 = surface_->evaluate(p0_local);
        const WorldVector3 diff = p1.world_pos() - p0.world_pos();
        return diff.normalized();
    }
    const auto& p1_local = uv_values_[idx+1];
    const auto& p0_local = uv_values_[idx];
    const auto p1 = surface_->evaluate(p1_local);
    const auto p0 = surface_->evaluate(p0_local);
    const WorldVector3 diff = p1.world_pos() - p0.world_pos();
    return diff.normalized();
}

WorldVector3 PathSegment::evaluate_normal(const ParamPoint1& param) const {
    if (t_values_.empty()) {
        throw std::runtime_error("Path segment has no points");
    }
    // Use param to get appropriate point
    auto it = std::lower_bound(t_values_.begin(), t_values_.end(), param);
    const size_t idx = std::min(
        static_cast<size_t>(std::distance(t_values_.begin(), it)),
        t_values_.size() - 1
    );
    const auto& local = uv_values_[idx];
    const auto geom = surface_->evaluate(local);
    const auto& derivs = geom.derivatives();
    return derivs[0].crossed(derivs[1]).normalized();
}

void TransitionPath::add_segment(
    std::shared_ptr<const Surface3D> surface,
    const ParamPoint1& t_start, const ParamPoint1& t_end,
    const ParamPoint2& uv_start, const ParamPoint2& uv_end,
    const WorldVector3& direction
) {
    if (!surface) {
        throw std::invalid_argument("Surface pointer cannot be null");
    }

    auto segment = std::make_unique<PathSegment>(surface);
    
    switch (surface->surface_type()) {
        case SurfaceType::Developable: {
            // For flat patches and cube faces
            segment->add_point(t_start, uv_start);
            segment->add_point(t_end, uv_end);
            break;
        }
        case SurfaceType::Smooth: {
            // For spheres, must have path solver
            if (auto solver = surface->get_path_solver()) {
                const auto start_pos = surface->evaluate(uv_start).world_pos();
                if (auto intersection = (*solver)(start_pos, direction, 1000.0)) {
                    segment->add_point(t_start, uv_start);
                    const auto hit_params = surface->nearest(intersection->position);
                    segment->add_point(t_end, hit_params);
                    break;
                }
            }
            throw std::runtime_error("Smooth surface missing required path solver");
        }
        default:
            throw std::runtime_error("Surface type does not support exact paths");
    }
    
    segments_.push_back(std::move(segment));
}

TransitionPath::TransitionPath()
    : Path3D(
        std::bind(&TransitionPath::evaluate_position, this, std::placeholders::_1),
        std::bind(&TransitionPath::evaluate_tangent, this, std::placeholders::_1),
        std::bind(&TransitionPath::evaluate_normal, this, std::placeholders::_1)
    )
{}

WorldPoint3 TransitionPath::evaluate_position(const ParamPoint1& param) const {
    if (segments_.empty()) {
        throw std::runtime_error("Transition path has no segments");
    }
    for (const auto& segment : segments_) {
        if (param <= segment->t_values().back()) {
            const auto& local = segment->uv_values().back();
            return segment->surface()->evaluate(local).world_pos();
        }
    }
    const auto& last = segments_.back();
    const auto& local = last->uv_values().back();
    return last->surface()->evaluate(local).world_pos();
}

WorldVector3 TransitionPath::evaluate_tangent(const ParamPoint1& param) const {
    if (segments_.empty()) {
        throw std::runtime_error("Transition path has no segments");
    }
    for (const auto& segment : segments_) {
        if (param <= segment->t_values().back()) {
            const auto& t_vals = segment->t_values();
            const auto& uv_vals = segment->uv_values();
            auto it = std::lower_bound(t_vals.begin(), t_vals.end(), param);
            const size_t idx = std::min(
                static_cast<size_t>(std::distance(t_vals.begin(), it)),
                t_vals.size() - 2
            );
            const auto& p1_local = uv_vals[idx+1];
            const auto& p0_local = uv_vals[idx];
            const auto p1 = segment->surface()->evaluate(p1_local);
            const auto p0 = segment->surface()->evaluate(p0_local);
            return (p1.world_pos() - p0.world_pos()).normalized();
        }
    }
    const auto& last = segments_.back();
    const auto& t_vals = last->t_values();
    const auto& uv_vals = last->uv_values();
    const size_t idx = t_vals.size() - 2;
    const auto& p1_local = uv_vals[idx+1];
    const auto& p0_local = uv_vals[idx];
    const auto p1 = last->surface()->evaluate(p1_local);
    const auto p0 = last->surface()->evaluate(p0_local);
    return (p1.world_pos() - p0.world_pos()).normalized();
}

WorldVector3 TransitionPath::evaluate_normal(const ParamPoint1& param) const {
    if (segments_.empty()) {
        throw std::runtime_error("Transition path has no segments");
    }
    for (const auto& segment : segments_) {
        if (param <= segment->t_values().back()) {
            const auto& local = segment->uv_values().front();
            const auto geom = segment->surface()->evaluate(local);
            const auto& derivs = geom.derivatives();
            return derivs[0].crossed(derivs[1]).normalized();
        }
    }
    const auto& last = segments_.back();
    const auto& local = last->uv_values().front();
    const auto geom = last->surface()->evaluate(local);
    const auto& derivs = geom.derivatives();
    return derivs[0].crossed(derivs[1]).normalized();
}

} // namespace shap
