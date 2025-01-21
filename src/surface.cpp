#include "../include/shap/surface.hpp"
#include "../include/shap/path.hpp"
#include <cmath>

namespace shap {

namespace {

// Implementation class for function-based surfaces
class FunctionSurface : public Surface {
public:
    FunctionSurface(
        PositionFunction pos,
        DerivativeFunction du = nullptr,
        DerivativeFunction dv = nullptr,
        DerivativeFunction duu = nullptr,
        DerivativeFunction duv = nullptr,
        DerivativeFunction dvv = nullptr,
        SurfaceType type = SurfaceType::Smooth
    ) : pos_func_(std::move(pos)),
        du_func_(std::move(du)),
        dv_func_(std::move(dv)),
        duu_func_(std::move(duu)),
        duv_func_(std::move(duv)),
        dvv_func_(std::move(dvv)),
        type_(type) {}

    GeometricProperties compute_properties(double u, double v) const override {
        GeometricProperties props;
        props.position = pos_func_(u, v);
        
        // Compute first derivatives
        if (du_func_ && dv_func_) {
            props.du = du_func_(u, v);
            props.dv = dv_func_(u, v);
        } else {
            // Numerical derivatives if not provided
            const double h = 1e-7;
            props.du = (pos_func_(u + h, v) - pos_func_(u - h, v)) * (0.5 / h);
            props.dv = (pos_func_(u, v + h) - pos_func_(u, v - h)) * (0.5 / h);
        }
        
        props.normal = compute_normal(props.du, props.dv);
        
        // Compute second derivatives if available
        if (duu_func_ && duv_func_ && dvv_func_) {
            props.duu = duu_func_(u, v);
            props.duv = duv_func_(u, v);
            props.dvv = dvv_func_(u, v);
            props.has_second_derivatives = true;
        }
        
        return props;
    }
    
    SurfaceType surface_type() const override {
        return type_;
    }

private:
    PositionFunction pos_func_;
    DerivativeFunction du_func_;
    DerivativeFunction dv_func_;
    DerivativeFunction duu_func_;
    DerivativeFunction duv_func_;
    DerivativeFunction dvv_func_;
    SurfaceType type_;
};

} // anonymous namespace

// Factory method implementations
std::shared_ptr<Surface> Surface::create(
    PositionFunction position_func,
    SurfaceType type
) {
    return std::make_shared<FunctionSurface>(
        std::move(position_func),
        nullptr, nullptr,  // No derivative functions
        nullptr, nullptr, nullptr,  // No second derivatives
        type
    );
}

std::shared_ptr<Surface> Surface::create_with_derivatives(
    PositionFunction position_func,
    DerivativeFunction du_func,
    DerivativeFunction dv_func,
    SurfaceType type
) {
    return std::make_shared<FunctionSurface>(
        std::move(position_func),
        std::move(du_func),
        std::move(dv_func),
        nullptr, nullptr, nullptr,  // No second derivatives
        type
    );
}

std::shared_ptr<Surface> Surface::create_with_all_derivatives(
    PositionFunction position_func,
    DerivativeFunction du_func,
    DerivativeFunction dv_func,
    DerivativeFunction duu_func,
    DerivativeFunction duv_func,
    DerivativeFunction dvv_func,
    SurfaceType type
) {
    return std::make_shared<FunctionSurface>(
        std::move(position_func),
        std::move(du_func),
        std::move(dv_func),
        std::move(duu_func),
        std::move(duv_func),
        std::move(dvv_func),
        type
    );
}

// Default path creation implementation
std::unique_ptr<SurfacePath> Surface::create_path(
    const SurfacePoint& start,
    const Vector& direction,
    double length
) const {
    auto path = std::make_unique<PathSegment>(
        std::const_pointer_cast<Surface>(
            std::shared_ptr<const Surface>(this, [](const Surface*){})
        )
    );
    
    // Create simple straight line path in parameter space
    const int num_points = 10;
    for (int i = 0; i < num_points; ++i) {
        double t = length * i / (num_points - 1);
        double u = start.u + direction.x * t;
        double v = start.v + direction.y * t;
        path->add_point(t, u, v);
    }
    
    return path;
}

} // namespace shap