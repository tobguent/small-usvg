#include <usvg/diffusion_curve.hpp>

namespace usvg
{
    bool DiffusionCurve::isValid() const
    {
        // positions valid?
        if (!position.isValid())
            return false;

        // colors valid?
        if (!colorLeft.isValid() || !colorRight.isValid())
            return false;

        // all good
        return true;
    }
}
