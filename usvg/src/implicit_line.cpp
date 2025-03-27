#include <usvg/implicit_line.hpp>

namespace usvg
{
    ImplicitLine2d::ImplicitLine2d(const Eigen::Vector2d& pnt1, const Eigen::Vector2d& pnt2) noexcept
    {
        Eigen::Vector2d delta = pnt2 - pnt1;

        a = -delta.y();
        b = delta.x();
        c = delta.y() * pnt1.x() - delta.x() * pnt1.y();

        // normalize such that we a^2 + b^2 = 1
        // Eq. (1) from T.W. Sederberg, T. Nishita, "Curve intersection using Bézier clipping" Computer-Aided Design 22(9), 538-549, 1990.
        double s = std::sqrt(a * a + b * b);
        if (s > 0)
        {
            a /= s;
            b /= s;
            c /= s;
        }
    }

    double ImplicitLine2d::distance(const Eigen::Vector2d& pnt) const noexcept
    {
        // Eq. (2) from T.W. Sederberg, T. Nishita, "Curve intersection using Bézier clipping" Computer-Aided Design 22(9), 538-549, 1990.
        return a * pnt.x() + b * pnt.y() + c;
    }

    Eigen::Vector4d ImplicitLine2d::distance(const std::array<Eigen::Vector2d, 4>& controlPoints) const noexcept
    {
        // We get the control points of the distance polynomial simply by computing the distance of each curve control point.
        // This is Eq. (18) from T.W. Sederberg, T. Nishita, "Curve intersection using Bézier clipping" Computer-Aided Design 22(9), 538-549, 1990.
        return Eigen::Vector4d(
            distance(controlPoints[0]),
            distance(controlPoints[1]),
            distance(controlPoints[2]),
            distance(controlPoints[3]));
    }
}
