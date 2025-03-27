#include <usvg/poisson_curve.hpp>

namespace usvg
{
    bool PoissonCurve::isValid() const
    {
        // positions valid?
        if (!position.isValid())
            return false;

        // weights valid?
        if (!weights.isValid())
            return false;

        // all good
        return true;
    }
}
