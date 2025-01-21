#include "../include/shap/path.hpp"
#include <stdexcept>
#include <cmath>

namespace shap {

void GeodesicCurve::compute_smooth_geodesic(
    const SurfacePoint& start,
    const SurfacePoint& end
) {
    // Initialize points array with start point
    const int num_points = 100;
    points_.resize(num_points);
    points_[0] = Point2D(start.u, start.v);
    
    // Simple straight line interpolation in parameter space for now
    // TODO: Implement proper geodesic computation using metric
    double du = (end.u - start.u) / (num_points - 1);
    double dv = (end.v - start.v) / (num_points - 1);
    
    for (int i = 1; i < num_points; ++i) {
        points_[i] = Point2D(
            start.u + i * du,
            start.v + i * dv
        );
    }
    
    // Store surface and parameter range
    t_start_ = 0.0;
    t_end_ = 1.0;
}

void GeodesicCurve::compute_developable_geodesic(
    const SurfacePoint& start,
    const SurfacePoint& end
) {
    // For developable surfaces, geodesics are straight lines in the developed space
    // For now, just use parameter space straight line like smooth case
    compute_smooth_geodesic(start, end);
}

SurfacePoint GeodesicCurve::evaluate(double t) const {
    if (t < t_start_ || t > t_end_) {
        throw std::out_of_range("Path parameter t out of range");
    }
    
    // Interpolate between points
    double normalized_t = (t - t_start_) / (t_end_ - t_start_);
    double index = normalized_t * (points_.size() - 1);
    int i = static_cast<int>(index);
    double frac = index - i;
    
    // Handle endpoint cases
    if (i >= points_.size() - 1) {
        auto& p = points_.back();
        return surface_->evaluate(p.x, p.y);
    }
    
    // Interpolate between points
    auto& p1 = points_[i];
    auto& p2 = points_[i + 1];
    double u = p1.x + frac * (p2.x - p1.x);
    double v = p1.y + frac * (p2.y - p1.y);
    
    return surface_->evaluate(u, v);
}

Vector GeodesicCurve::tangent(double t) const {
    // Compute tangent using central difference
    const double h = 1e-7;
    auto pt1 = evaluate(t - h);
    auto pt2 = evaluate(t + h);
    return (pt2.position - pt1.position).normalize();
}

Vector GeodesicCurve::normal(double t) const {
    return evaluate(t).normal;
}

// Path segment implementation
void PathSegment::add_point(double t, double u, double v) {
    points_.push_back(Point(t, u, v));
}

SurfacePoint PathSegment::evaluate(double t) const {
    if (points_.empty()) {
        throw std::runtime_error("No points in path segment");
    }
    
    // Find surrounding points
    auto it = std::lower_bound(
        points_.begin(), points_.end(), t,
        [](const Point& p, double val) { return p.x < val; }
    );
    
    if (it == points_.begin()) {
        return surface_->evaluate(it->y, it->z);
    }
    if (it == points_.end()) {
        auto& last = points_.back();
        return surface_->evaluate(last.y, last.z);
    }
    
    // Interpolate between points
    auto& p1 = *(it - 1);
    auto& p2 = *it;
    double frac = (t - p1.x) / (p2.x - p1.x);
    
    double u = p1.y + frac * (p2.y - p1.y);
    double v = p1.z + frac * (p2.z - p1.z);
    
    return surface_->evaluate(u, v);
}

Vector PathSegment::tangent(double t) const {
    // Use surface derivatives for tangent
    auto pt = evaluate(t);
    auto props = surface_->compute_properties(pt.u, pt.v);
    
    // Find velocity in parameter space
    auto it = std::lower_bound(
        points_.begin(), points_.end(), t,
        [](const Point& p, double val) { return p.x < val; }
    );
    
    double du_dt, dv_dt;
    if (it == points_.begin() || it == points_.end()) {
        // Use one-sided difference at endpoints
        if (points_.size() < 2) {
            throw std::runtime_error("Need at least 2 points for tangent");
        }
        if (it == points_.begin()) {
            auto& p1 = points_[0];
            auto& p2 = points_[1];
            double dt = p2.x - p1.x;
            du_dt = (p2.y - p1.y) / dt;
            dv_dt = (p2.z - p1.z) / dt;
        } else {
            auto& p1 = points_[points_.size() - 2];
            auto& p2 = points_[points_.size() - 1];
            double dt = p2.x - p1.x;
            du_dt = (p2.y - p1.y) / dt;
            dv_dt = (p2.z - p1.z) / dt;
        }
    } else {
        // Use central difference
        auto& prev = *(it - 1);
        auto& next = *it;
        double dt = next.x - prev.x;
        du_dt = (next.y - prev.y) / dt;
        dv_dt = (next.z - prev.z) / dt;
    }
    
    // Compute tangent vector
    return (props.du * du_dt + props.dv * dv_dt).normalize();
}

Vector PathSegment::normal(double t) const {
    return evaluate(t).normal;
}

// Transition path implementation
void TransitionPath::add_segment(
    std::shared_ptr<Surface> surface,
    double t_start, double t_end,
    double u_start, double u_end,
    double v_start, double v_end,
    const Vector& direction
) {
    auto segment = std::make_unique<PathSegment>(surface);
    
    // Add points along segment
    const int num_points = 10;
    for (int i = 0; i < num_points; ++i) {
        double t = t_start + (t_end - t_start) * i / (num_points - 1);
        double u = u_start + (u_end - u_start) * i / (num_points - 1);
        double v = v_start + (v_end - v_start) * i / (num_points - 1);
        segment->add_point(t, u, v);
    }
    
    segments_.push_back(std::move(segment));
}

SurfacePoint TransitionPath::evaluate(double t) const {
    // Find segment containing t
    for (const auto& segment : segments_) {
        auto& points = segment->points();
        if (!points.empty() && t <= points.back().x) {
            return segment->evaluate(t);
        }
    }
    
    // If t is past end, return last point
    if (!segments_.empty()) {
        auto& last_segment = segments_.back();
        auto& points = last_segment->points();
        if (!points.empty()) {
            auto& last = points.back();
            return last_segment->surface()->evaluate(last.y, last.z);
        }
    }
    
    throw std::runtime_error("Invalid path parameter t");
}

Vector TransitionPath::tangent(double t) const {
    // Find segment containing t
    for (const auto& segment : segments_) {
        auto& points = segment->points();
        if (!points.empty() && t <= points.back().x) {
            return segment->tangent(t);
        }
    }
    
    // If t is past end, use last segment
    if (!segments_.empty()) {
        return segments_.back()->tangent(t);
    }
    
    throw std::runtime_error("Invalid path parameter t");
}

Vector TransitionPath::normal(double t) const {
    // Find segment containing t
    for (const auto& segment : segments_) {
        auto& points = segment->points();
        if (!points.empty() && t <= points.back().x) {
            return segment->normal(t);
        }
    }
    
    // If t is past end, use last segment
    if (!segments_.empty()) {
        return segments_.back()->normal(t);
    }
    
    throw std::runtime_error("Invalid path parameter t");
}

} // namespace shap