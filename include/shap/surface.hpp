#pragma once
#include "geometry.hpp"
#include <functional>
#include <memory>
#include <vector>
#include <optional>
#include <cmath>

namespace shap {

// Forward declarations
class SurfaceCollection;
class SurfacePath;

// Base type for storing any surface type
struct SurfaceBase {
    virtual ~SurfaceBase() = default;
    virtual Point operator()(double u, double v) const = 0;
    virtual SurfacePoint evaluate(double u, double v) const = 0;
    virtual Point du(double u, double v) const = 0;
    virtual Point dv(double u, double v) const = 0;
    virtual MetricTensor metric_tensor(double u, double v) const = 0;
};

// Base class for parametric surfaces with Riemannian geometry
template<typename Derived>
class Surface : public SurfaceBase {
public:
    // Basic evaluation
    Point operator()(double u, double v) const override {
        return static_cast<const Derived&>(*this)(u, v);
    }
    
    // Surface name for identification
    std::string name;
    
    // Get complete geometric data at a point
    SurfacePoint evaluate(double u, double v) const override {
        Point pos = operator()(u, v);
        Point du_vec = du(u, v);
        Point dv_vec = dv(u, v);
        Point n = du_vec.cross(dv_vec).normalize();
        
        return SurfacePoint(name, u, v, pos, n, du_vec, dv_vec);
    }
    
    // First partial derivatives
    Point du(double u, double v) const override {
        const double h = 1e-7;
        return (operator()(u + h, v) - operator()(u - h, v)) * (0.5 / h);
    }
    
    Point dv(double u, double v) const override {
        const double h = 1e-7;
        return (operator()(u, v + h) - operator()(u, v - h)) * (0.5 / h);
    }
    
    // Metric tensor and Riemannian connection
    MetricTensor metric_tensor(double u, double v) const override {
        Point du_vec = du(u, v);
        Point dv_vec = dv(u, v);
        
        return MetricTensor(
            du_vec.dot(du_vec),      // g11
            du_vec.dot(dv_vec),      // g12
            du_vec.dot(dv_vec),      // g21
            dv_vec.dot(dv_vec)       // g22
        );
    }
    
    // Create paths on surface
    std::unique_ptr<SurfacePath> create_geodesic(
        const SurfacePoint& start,
        const SurfacePoint& end
    ) const {
        return create_geodesic_path(start, end);
    }
    
    std::unique_ptr<SurfacePath> create_directional_path(
        const SurfacePoint& start,
        const Vector& direction,
        double length
    ) const {
        return create_directional_path(start, direction, length);
    }
    
    // Parallel transport a vector along a path
    Vector parallel_transport(
        const Vector& v,
        const SurfacePath& path,
        double t_start,
        double t_end
    ) const;
    
    // Find transition to adjacent surface (if any)
    virtual std::optional<SurfacePoint> find_transition(
        const SurfacePoint& point,
        const Vector& direction
    ) const {
        return std::nullopt;
    }
};

// Helper for creating surfaces from lambdas
template<typename F>
struct ParametricSurface : Surface<ParametricSurface<F>> {
    F func;
    
    ParametricSurface(F f) : func(std::move(f)) {}
    
    Point operator()(double u, double v) const {
        return func(u, v);
    }
};

template<typename F>
auto make_surface(F&& f) {
    return ParametricSurface<F>(std::forward<F>(f));
}

// Path that follows a surface with transitions
class TransitionPath : public SurfacePath {
    struct Segment {
        std::shared_ptr<SurfaceBase> surface;
        double t_start, t_end;  // Path parameter range
        double u_start, u_end;  // Surface parameter range in u
        double v_start, v_end;  // Surface parameter range in v
        Vector direction;       // Direction in surface parameters
    };
    
    std::vector<Segment> segments;
    
public:
    SurfacePoint evaluate(double t) const override {
        // Find segment containing t
        for (const auto& seg : segments) {
            if (t >= seg.t_start && t <= seg.t_end) {
                double local_t = (t - seg.t_start) / (seg.t_end - seg.t_start);
                
                // Get base point on surface
                double u = seg.u_start + local_t * (seg.u_end - seg.u_start);
                double v = seg.v_start + local_t * (seg.v_end - seg.v_start);
                return seg.surface->evaluate(u, v);
            }
        }
        throw std::runtime_error("Invalid path parameter");
    }
    
    Vector tangent(double t) const override {
        for (const auto& seg : segments) {
            if (t >= seg.t_start && t <= seg.t_end) {
                // Get tangent from surface metric
                auto metric = seg.surface->metric_tensor(seg.u_start, seg.v_start);
                auto [du, dv] = metric.raise_indices(seg.direction.x, seg.direction.y);
                return Vector(du, dv, 0).normalize();
            }
        }
        throw std::runtime_error("Invalid path parameter");
    }
    
    std::unique_ptr<SurfacePath> offset(double distance) const override {
        auto result = std::make_unique<TransitionPath>();
        
        // Offset each segment along surface normal
        for (const auto& seg : segments) {
            auto new_seg = seg;
            auto normal = seg.surface->evaluate(seg.u_start, seg.v_start).normal;
            new_seg.u_start += normal.x * distance;
            new_seg.v_start += normal.y * distance;
            new_seg.u_end += normal.x * distance;
            new_seg.v_end += normal.y * distance;
            result->segments.push_back(new_seg);
        }
        
        return result;
    }
    
    std::unique_ptr<SurfacePath> smooth(double radius) const override {
        auto result = std::make_unique<TransitionPath>();
        
        // Add circular arcs at segment transitions
        for (size_t i = 0; i < segments.size(); ++i) {
            const auto& seg = segments[i];
            result->segments.push_back(seg);
            
            if (i < segments.size() - 1) {
                const auto& next = segments[i + 1];
                
                // Create circular arc between segments
                double arc_t_start = seg.t_end - radius;
                double arc_t_end = next.t_start + radius;
                
                auto arc_seg = seg;
                arc_seg.t_start = arc_t_start;
                arc_seg.t_end = arc_t_end;
                arc_seg.direction = (seg.direction + next.direction).normalize();
                result->segments.push_back(arc_seg);
            }
        }
        
        return result;
    }
    
    void add_segment(
        std::shared_ptr<SurfaceBase> surface,
        double t_start, double t_end,
        double u_start, double u_end,
        double v_start, double v_end,
        const Vector& direction
    ) {
        segments.push_back({
            surface,
            t_start, t_end,
            u_start, u_end,
            v_start, v_end,
            direction
        });
    }
};

// Collection of connected surfaces
class SurfaceCollection {
    struct Connection {
        std::shared_ptr<SurfaceBase> surface1;
        std::shared_ptr<SurfaceBase> surface2;
        std::function<bool(const SurfacePoint&, const Vector&)> transition_test;
        std::function<SurfacePoint(const SurfacePoint&)> map_point;
    };
    
    std::vector<std::shared_ptr<SurfaceBase>> surfaces;
    std::vector<Connection> connections;
    
public:
    template<typename S>
    void add_surface(S&& surface) {
        surfaces.push_back(std::make_shared<S>(std::forward<S>(surface)));
    }
    
    template<typename S1, typename S2>
    void add_connection(
        S1&& s1, S2&& s2,
        std::function<bool(const SurfacePoint&, const Vector&)> test,
        std::function<SurfacePoint(const SurfacePoint&)> map
    ) {
        auto wrapped1 = std::make_shared<S1>(std::forward<S1>(s1));
        auto wrapped2 = std::make_shared<S2>(std::forward<S2>(s2));
        connections.push_back({wrapped1, wrapped2, test, map});
    }
    
    // Store surfaces by name for lookup
    std::unordered_map<std::string, std::shared_ptr<SurfaceBase>> surface_map;
    
    // Add named surface
    template<typename S>
    void add_surface(std::string name, S&& surface) {
        auto wrapped = std::make_shared<S>(std::forward<S>(surface));
        surfaces.push_back(wrapped);
        surface_map[name] = wrapped;
    }
    
    // Get surface by name
    std::shared_ptr<SurfaceBase> get_surface(const std::string& name) const {
        auto it = surface_map.find(name);
        if (it == surface_map.end()) {
            throw std::runtime_error("Surface not found: " + name);
        }
        return it->second;
    }
    
    // Find surface and connection for a point
    struct SurfaceInfo {
        std::shared_ptr<SurfaceBase> surface;
        const Connection* next_connection;
        double distance_to_transition;
        
        SurfaceInfo(
            std::shared_ptr<SurfaceBase> s,
            const Connection* conn,
            double dist
        ) : surface(s), next_connection(conn), distance_to_transition(dist) {}
    };
    
    SurfaceInfo find_surface_info(const SurfacePoint& pt, const Vector& dir) const {
        // Find current surface
        std::shared_ptr<SurfaceBase> current_surface;
        double min_dist = std::numeric_limits<double>::max();
        
        for (const auto& surface : surfaces) {
            Point sp = surface->operator()(pt.u, pt.v);
            double dist = (sp.x - pt.position.x) * (sp.x - pt.position.x) +
                         (sp.y - pt.position.y) * (sp.y - pt.position.y) +
                         (sp.z - pt.position.z) * (sp.z - pt.position.z);
            if (dist < min_dist) {
                min_dist = dist;
                current_surface = surface;
            }
        }
        
        // Find next transition
        for (const auto& conn : connections) {
            if (conn.surface1 == current_surface && conn.transition_test(pt, dir)) {
                // Calculate distance to transition
                double dist_to_edge = 1.0 - pt.u;  // Assuming transition at u=1
                return SurfaceInfo(current_surface, &conn, dist_to_edge);
            }
        }
        
        return SurfaceInfo(current_surface, nullptr, std::numeric_limits<double>::max());
    }
    
    // Create path that can transition between surfaces
    std::unique_ptr<SurfacePath> create_path(
        const SurfacePoint& start,
        const Vector& direction,
        double length
    ) {
        auto path = std::make_unique<TransitionPath>();
        
        // Start with first segment
        double t = 0.0;
        SurfacePoint current = start;
        Vector current_dir = direction;
        
        while (t < length) {
            // Find current surface and check for transition
            auto info = find_surface_info(current, current_dir);
            
            if (info.next_connection) {
                // Add segment up to transition
                double segment_length = std::min(0.25, info.distance_to_transition);
                path->add_segment(
                    info.surface, t, t + segment_length,
                    current.u, std::min(1.0, current.u + current_dir.x * segment_length),
                    current.v, current.v + current_dir.y * segment_length,
                    current_dir
                );
                
                // Transition to next surface
                current = info.next_connection->map_point(current);
                t += segment_length;
            } else {
                // Add segment on current surface
                double segment_length = 0.25;
                path->add_segment(
                    info.surface, t, t + segment_length,
                    current.u, std::min(1.0, current.u + current_dir.x * segment_length),
                    current.v, current.v + current_dir.y * segment_length,
                    current_dir
                );
                
                // Update current point
                current.u = std::min(1.0, current.u + current_dir.x * segment_length);
                current.v += current_dir.y * segment_length;
                current = info.surface->evaluate(current.u, current.v);
                t += segment_length;
            }
        }
        
        return path;
    }
};

namespace surfaces {

// Helper to create a square face
inline auto square_face(auto transform) {
    return make_surface([transform](double u, double v) {
        // u,v in [0,1]
        return transform(u, v);
    });
}

inline auto sphere(double radius = 1.0) {
    return make_surface([radius](double u, double v) {
        // u: longitude [0, 2π]
        // v: latitude [0, π]
        return Point(
            radius * std::cos(u) * std::sin(v),
            radius * std::sin(u) * std::sin(v),
            radius * std::cos(v)
        );
    });
}

// Create a cube as a collection of connected faces
inline SurfaceCollection cube(double size = 1.0) {
    SurfaceCollection cube;
    
    using SurfaceFunc = std::function<Point(double,double)>;
    
    // Create faces
    auto make_face = [](SurfaceFunc f) {
        return make_surface(std::move(f));
    };
    
    auto front = make_face(
        [size](double u, double v) {
            return Point(size * (2*u - 1), size, size * (2*v - 1));
        }
    );
    
    auto right = make_face(
        [size](double u, double v) {
            return Point(size, size * (1 - 2*u), size * (2*v - 1));
        }
    );
    
    auto back = make_face(
        [size](double u, double v) {
            return Point(size * (1 - 2*u), -size, size * (2*v - 1));
        }
    );
    
    auto left = make_face(
        [size](double u, double v) {
            return Point(-size, size * (2*u - 1), size * (2*v - 1));
        }
    );
    
    // Add faces to collection with names
    cube.add_surface("front", std::move(front));
    cube.add_surface("right", std::move(right));
    cube.add_surface("back", std::move(back));
    cube.add_surface("left", std::move(left));
    
    // Add all face connections
    // Front -> Right
    cube.add_connection(
        cube.get_surface("front"), cube.get_surface("right"),
        [](const SurfacePoint& pt, const Vector& dir) {
            return pt.u >= 0.95 && dir.x > 0;
        },
        [](const SurfacePoint& pt) {
            return SurfacePoint(
                "right",            // Next surface
                0.0, pt.v,         // Map to left edge
                pt.position,       // Keep position
                Vector(1, 0, 0),   // Normal points right
                Vector(0, -1, 0),  // du points back
                Vector(0, 0, 1)    // dv points up
            );
        }
    );
    
    // Right -> Back
    cube.add_connection(
        cube.get_surface("right"), cube.get_surface("back"),
        [](const SurfacePoint& pt, const Vector& dir) {
            return pt.u >= 0.95 && dir.x < 0;
        },
        [](const SurfacePoint& pt) {
            return SurfacePoint(
                "back",            // Next surface
                0.0, pt.v,        // Map to left edge
                pt.position,      // Keep position
                Vector(0, -1, 0), // Normal points back
                Vector(-1, 0, 0), // du points left
                Vector(0, 0, 1)   // dv points up
            );
        }
    );
    
    // Back -> Left
    cube.add_connection(
        cube.get_surface("back"), cube.get_surface("left"),
        [](const SurfacePoint& pt, const Vector& dir) {
            return pt.u >= 0.95 && dir.x < 0;
        },
        [](const SurfacePoint& pt) {
            return SurfacePoint(
                "left",           // Next surface
                0.0, pt.v,       // Map to left edge
                pt.position,     // Keep position
                Vector(-1, 0, 0), // Normal points left
                Vector(0, 1, 0),  // du points front
                Vector(0, 0, 1)   // dv points up
            );
        }
    );
    
    // Left -> Front
    cube.add_connection(
        cube.get_surface("left"), cube.get_surface("front"),
        [](const SurfacePoint& pt, const Vector& dir) {
            return pt.u >= 0.95 && dir.x > 0;
        },
        [](const SurfacePoint& pt) {
            return SurfacePoint(
                "front",          // Next surface
                0.0, pt.v,       // Map to left edge
                pt.position,     // Keep position
                Vector(0, 1, 0),  // Normal points front
                Vector(1, 0, 0),  // du points right
                Vector(0, 0, 1)   // dv points up
            );
        }
    );
    
    return cube;
}

} // namespace surfaces

} // namespace shap
