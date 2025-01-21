#include "shap/path.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace shap {

// GeodesicCurve implementation
GeodesicCurve::GeodesicCurve(
    std::shared_ptr<Surface> surface,
    const SurfacePoint& start,
    const SurfacePoint& end
) : surface_(surface) {
    switch (surface_->surface_type()) {
        case Surface::SurfaceType::Smooth:
            compute_smooth_geodesic(start, end);
            break;
            
        case Surface::SurfaceType::Developable:
            compute_developable_geodesic(start, end);
            break;
            
        case Surface::SurfaceType::NonSmooth:
            throw std::runtime_error(
                "Cannot compute geodesics on non-smooth surface: " + 
                surface_->name
            );
    }
}

void GeodesicCurve::compute_smooth_geodesic(
    const SurfacePoint& start,
    const SurfacePoint& end
) {
    // Solve geodesic equation using RK4 method
    const int steps = 100;
    points_.resize(steps + 1);
    tangents_.resize(steps + 1);
    
    // Initial conditions
    points_[0] = Point(start.u, start.v, 0);
    Vector dir = end.position - start.position;
    tangents_[0] = Vector(dir.x, dir.y, 0).normalize();
    
    for (int i = 1; i <= steps; ++i) {
        double t = static_cast<double>(i) / steps;
        
        // Get metric tensor and Christoffel symbols
        auto metric = surface_->metric_tensor(points_[i-1].x, points_[i-1].y);
        
        // RK4 step for geodesic equation
        auto rk4_step = [&](double u, double v, double du, double dv) {
            std::array<double,2> accel = metric.christoffel_second(0, u, v);
            double d2u = -accel[0] * du * du - 2 * accel[1] * du * dv;
            accel = metric.christoffel_second(1, u, v);
            double d2v = -accel[0] * du * du - 2 * accel[1] * du * dv;
            return std::make_tuple(du, dv, d2u, d2v);
        };
        
        // RK4 integration
        double h = 1.0 / steps;
        auto [k1_du, k1_dv, k1_d2u, k1_d2v] = rk4_step(
            points_[i-1].x, points_[i-1].y,
            tangents_[i-1].x, tangents_[i-1].y
        );
        
        auto [k2_du, k2_dv, k2_d2u, k2_d2v] = rk4_step(
            points_[i-1].x + 0.5*h*k1_du,
            points_[i-1].y + 0.5*h*k1_dv,
            tangents_[i-1].x + 0.5*h*k1_d2u,
            tangents_[i-1].y + 0.5*h*k1_d2v
        );
        
        auto [k3_du, k3_dv, k3_d2u, k3_d2v] = rk4_step(
            points_[i-1].x + 0.5*h*k2_du,
            points_[i-1].y + 0.5*h*k2_dv,
            tangents_[i-1].x + 0.5*h*k2_d2u,
            tangents_[i-1].y + 0.5*h*k2_d2v
        );
        
        auto [k4_du, k4_dv, k4_d2u, k4_d2v] = rk4_step(
            points_[i-1].x + h*k3_du,
            points_[i-1].y + h*k3_dv,
            tangents_[i-1].x + h*k3_d2u,
            tangents_[i-1].y + h*k3_d2v
        );
        
        // Update position and velocity
        points_[i].x = points_[i-1].x + (h/6) * (k1_du + 2*k2_du + 2*k3_du + k4_du);
        points_[i].y = points_[i-1].y + (h/6) * (k1_dv + 2*k2_dv + 2*k3_dv + k4_dv);
        tangents_[i].x = tangents_[i-1].x + (h/6) * (k1_d2u + 2*k2_d2u + 2*k3_d2u + k4_d2u);
        tangents_[i].y = tangents_[i-1].y + (h/6) * (k1_d2v + 2*k2_d2v + 2*k3_d2v + k4_d2v);
        tangents_[i] = tangents_[i].normalize();
    }
}

void GeodesicCurve::compute_developable_geodesic(
    const SurfacePoint& start,
    const SurfacePoint& end
) {
    // For developable surfaces, geodesics are straight lines in the developed space
    auto [u1, v1] = surface_->develop_point(start.u, start.v);
    auto [u2, v2] = surface_->develop_point(end.u, end.v);
    
    // Create straight line in parameter space
    const int steps = 100;
    points_.resize(steps + 1);
    tangents_.resize(steps + 1);
    
    Vector dir(u2 - u1, v2 - v1, 0);
    Vector tangent = dir.normalize();
    
    for (int i = 0; i <= steps; ++i) {
        double t = static_cast<double>(i) / steps;
        points_[i] = Point(u1, v1, 0) + dir * t;
        tangents_[i] = tangent;
    }
}

SurfacePoint GeodesicCurve::evaluate(double t) const {
    // Interpolate between stored points
    double idx = t * (points_.size() - 1);
    int i = static_cast<int>(idx);
    double frac = idx - i;
    
    if (i >= points_.size() - 1) {
        return surface_->evaluate(points_.back().x, points_.back().y);
    }
    
    // Linear interpolation between points
    Point p = points_[i] * (1-frac) + points_[i+1] * frac;
    return surface_->evaluate(p.x, p.y);
}

Vector GeodesicCurve::tangent(double t) const {
    // Interpolate between stored tangents
    double idx = t * (tangents_.size() - 1);
    int i = static_cast<int>(idx);
    double frac = idx - i;
    
    if (i >= tangents_.size() - 1) {
        return tangents_.back();
    }
    
    // Spherical linear interpolation between tangents
    Vector v1 = tangents_[i];
    Vector v2 = tangents_[i+1];
    double dot = v1.dot(v2);
    
    if (dot > 0.9999) {
        return (v1 * (1-frac) + v2 * frac).normalize();
    }
    
    double angle = std::acos(dot);
    double sin_angle = std::sin(angle);
    return (v1 * std::sin((1-frac)*angle) + v2 * std::sin(frac*angle)) * (1.0/sin_angle);
}

// CircularArc implementation
CircularArc::CircularArc(
    const SurfacePoint& start,
    const Vector& start_tangent,
    const SurfacePoint& end,
    const Vector& end_tangent,
    double radius
) : radius_(radius) {
    // Compute arc center and orientation
    Vector chord = end.position - start.position;
    normal_ = start_tangent.cross(end_tangent).normalize();
    
    // Set up local coordinate system
    x_axis_ = start_tangent.normalize();
    y_axis_ = normal_.cross(x_axis_).normalize();
    
    // Find arc parameters
    double chord_length = chord.norm();
    double angle = std::acos(start_tangent.dot(end_tangent));
    
    // Compute center point
    double dist_to_center = radius_ / std::tan(angle/2);
    center_ = start.position + x_axis_ * dist_to_center + y_axis_ * radius_;
    
    // Store arc angles
    start_angle_ = 0;
    sweep_angle_ = angle;
}

SurfacePoint CircularArc::evaluate(double t) const {
    double angle = start_angle_ + t * sweep_angle_;
    Point pos = center_ + 
        x_axis_ * (radius_ * std::cos(angle)) +
        y_axis_ * (radius_ * std::sin(angle));
    
    // Create surface point with approximate parameters
    return SurfacePoint(
        "transition",
        t, 0,  // Approximate parameters
        pos,
        normal_,
        tangent(t),
        normal_.cross(tangent(t))
    );
}

Vector CircularArc::tangent(double t) const {
    double angle = start_angle_ + t * sweep_angle_;
    return (
        x_axis_ * (-std::sin(angle)) +
        y_axis_ * (std::cos(angle))
    ).normalize();
}

Vector CircularArc::normal(double t) const {
    return normal_;
}

// PathSegment implementation
PathSegment::PathSegment(
    std::shared_ptr<Surface> surface,
    double t_start, double t_end,
    double u_start, double u_end,
    double v_start, double v_end,
    const Vector& direction
) : surface(surface),
    t_start(t_start), t_end(t_end),
    u_start(u_start), u_end(u_end),
    v_start(v_start), v_end(v_end),
    direction(direction) {
    compute_normal_field();
}

void PathSegment::compute_normal_field() {
    // Sample normals along segment
    const int samples = 50;
    normal_samples_.resize(samples);
    
    for (int i = 0; i < samples; ++i) {
        double t = static_cast<double>(i) / (samples - 1);
        double u = u_start + t * (u_end - u_start);
        double v = v_start + t * (v_end - v_start);
        normal_samples_[i] = surface->evaluate(u, v).normal;
    }
}

SurfacePoint PathSegment::evaluate(double t) const {
    double local_t = (t - t_start) / (t_end - t_start);
    double u = u_start + local_t * (u_end - u_start);
    double v = v_start + local_t * (v_end - v_start);
    return surface->evaluate(u, v);
}

Vector PathSegment::tangent(double t) const {
    double local_t = (t - t_start) / (t_end - t_start);
    double u = u_start + local_t * (u_end - u_start);
    double v = v_start + local_t * (v_end - v_start);
    auto metric = surface->metric_tensor(u, v);
    auto [du, dv] = metric.raise_indices(direction.x, direction.y, u, v);
    return Vector(du, dv, 0).normalize();
}

Vector PathSegment::normal(double t) const {
    double local_t = (t - t_start) / (t_end - t_start);
    int idx = static_cast<int>(local_t * (normal_samples_.size() - 1));
    idx = std::clamp(idx, 0, static_cast<int>(normal_samples_.size()) - 1);
    return normal_samples_[idx];
}

// PathTransition implementation
PathTransition::PathTransition(
    const PathSegment& seg1,
    const PathSegment& seg2,
    TransitionType type,
    double radius
) : type_(type),
    start_(seg1.evaluate(seg1.t_end)),
    end_(seg2.evaluate(seg2.t_start)),
    start_tangent_(seg1.tangent(seg1.t_end)),
    end_tangent_(seg2.tangent(seg2.t_start)) {
    
    // Check if geodesic transition is possible
    if (type == TransitionType::Geodesic) {
        auto surface_type = seg1.surface->surface_type();
        if (surface_type == Surface::SurfaceType::NonSmooth ||
            (seg1.surface != seg2.surface && 
             surface_type != Surface::SurfaceType::Developable)) {
            throw std::runtime_error(
                "Geodesic transition not possible between surfaces: " +
                seg1.surface->name + " and " + seg2.surface->name
            );
        }
    }
    
    switch (type) {
        case TransitionType::Circular:
            arc_ = std::make_unique<CircularArc>(
                start_, start_tangent_,
                end_, end_tangent_,
                radius
            );
            break;
            
        case TransitionType::Geodesic:
            if (seg1.surface == seg2.surface) {
                geodesic_ = std::make_unique<GeodesicCurve>(
                    seg1.surface,
                    start_,
                    end_
                );
            }
            break;
            
        default: // Linear
            break;
    }
}

SurfacePoint PathTransition::evaluate(double t) const {
    switch (type_) {
        case TransitionType::Circular:
            return arc_->evaluate(t);
        case TransitionType::Geodesic:
            return geodesic_->evaluate(t);
        default:
            // Linear interpolation
            return SurfacePoint(
                "transition",
                t, 0,  // Approximate parameters
                start_.position * (1-t) + end_.position * t,
                start_.normal * (1-t) + end_.normal * t,
                start_tangent_ * (1-t) + end_tangent_ * t,
                Vector()  // Dummy dv vector
            );
    }
}

Vector PathTransition::tangent(double t) const {
    switch (type_) {
        case TransitionType::Circular:
            return arc_->tangent(t);
        case TransitionType::Geodesic:
            return geodesic_->tangent(t);
        default:
            return (start_tangent_ * (1-t) + end_tangent_ * t).normalize();
    }
}

Vector PathTransition::normal(double t) const {
    switch (type_) {
        case TransitionType::Circular:
            return arc_->normal(t);
        case TransitionType::Geodesic:
            return start_.normal * (1-t) + end_.normal * t;
        default:
            return start_.normal * (1-t) + end_.normal * t;
    }
}

// TransitionPath implementation
void TransitionPath::add_segment(
    std::shared_ptr<Surface> surface,
    double t_start, double t_end,
    double u_start, double u_end,
    double v_start, double v_end,
    const Vector& direction
) {
    segments_.emplace_back(
        surface, t_start, t_end,
        u_start, u_end, v_start, v_end,
        direction
    );
}

SurfacePoint TransitionPath::evaluate(double t) const {
    // Find relevant segment or transition
    for (size_t i = 0; i < segments_.size(); ++i) {
        const auto& seg = segments_[i];
        if (t >= seg.t_start && t <= seg.t_end) {
            return seg.evaluate(t);
        }
        if (i < transitions_.size()) {
            const auto& trans = transitions_[i];
            if (t > seg.t_end && t < segments_[i+1].t_start) {
                double local_t = (t - seg.t_end) / 
                    (segments_[i+1].t_start - seg.t_end);
                return trans.evaluate(local_t);
            }
        }
    }
    throw std::runtime_error("Invalid path parameter");
}

Vector TransitionPath::tangent(double t) const {
    for (size_t i = 0; i < segments_.size(); ++i) {
        const auto& seg = segments_[i];
        if (t >= seg.t_start && t <= seg.t_end) {
            return seg.tangent(t);
        }
        if (i < transitions_.size()) {
            const auto& trans = transitions_[i];
            if (t > seg.t_end && t < segments_[i+1].t_start) {
                double local_t = (t - seg.t_end) / 
                    (segments_[i+1].t_start - seg.t_end);
                return trans.tangent(local_t);
            }
        }
    }
    throw std::runtime_error("Invalid path parameter");
}

Vector TransitionPath::normal(double t) const {
    for (size_t i = 0; i < segments_.size(); ++i) {
        const auto& seg = segments_[i];
        if (t >= seg.t_start && t <= seg.t_end) {
            return seg.normal(t);
        }
        if (i < transitions_.size()) {
            const auto& trans = transitions_[i];
            if (t > seg.t_end && t < segments_[i+1].t_start) {
                double local_t = (t - seg.t_end) / 
                    (segments_[i+1].t_start - seg.t_end);
                return trans.normal(local_t);
            }
        }
    }
    throw std::runtime_error("Invalid path parameter");
}

std::unique_ptr<SurfacePath> TransitionPath::offset(
    double distance,
    TransitionType transition
) const {
    auto result = std::make_unique<TransitionPath>();
    
    // Offset each segment
    for (const auto& seg : segments_) {
        // Sample points along segment
        const int samples = 50;
        std::vector<Point> offset_points;
        offset_points.reserve(samples);
        
        for (int i = 0; i < samples; ++i) {
            double t = static_cast<double>(i) / (samples - 1);
            double local_t = seg.t_start + t * (seg.t_end - seg.t_start);
            
            // Get point and normal
            auto pt = seg.evaluate(local_t);
            auto n = seg.normal(local_t);
            
            // Offset point
            Point offset_pos = pt.position + n * distance;
            
            // Project back to surface maintaining offset
            auto [u, v] = project_to_surface_with_offset(
                seg.surface, offset_pos, distance
            );
            
            offset_points.push_back(Point(u, v, 0));
        }
        
        // Create new segment from offset points
        result->add_segment(
            seg.surface,
            seg.t_start, seg.t_end,
            offset_points.front().x, offset_points.back().x,
            offset_points.front().y, offset_points.back().y,
            seg.direction
        );
    }
    
    // Update transitions
    result->update_transitions(transition);
    
    return result;
}

void TransitionPath::update_transitions(TransitionType type) {
    transitions_.clear();
    
    for (size_t i = 0; i < segments_.size() - 1; ++i) {
        transitions_.emplace_back(
            segments_[i],
            segments_[i+1],
            type,
            0.1  // Default radius
        );
    }
}

std::pair<double,double> TransitionPath::project_to_surface_with_offset(
    std::shared_ptr<Surface> surface,
    const Point& point,
    double offset_distance
) {
    // Find closest point on surface that maintains offset distance
    double min_u = 0, min_v = 0;
    double min_error = std::numeric_limits<double>::max();
    
    const int samples = 10;
    for (int i = 0; i <= samples; ++i) {
        double u = static_cast<double>(i) / samples;
        for (int j = 0; j <= samples; ++j) {
            double v = static_cast<double>(j) / samples;
            
            auto pt = surface->evaluate(u, v);
            double actual_offset = (point - pt.position).norm();
            double error = std::abs(actual_offset - offset_distance);
            
            if (error < min_error) {
                min_error = error;
                min_u = u;
                min_v = v;
            }
        }
    }
    
    return {min_u, min_v};
}

} // namespace shap