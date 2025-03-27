#pragma once

#include "bezier_curve.hpp"

namespace usvg
{
    /**
     * @brief Represents a normalized implicit line in 2D. The line is infinitely long.
     */
    class ImplicitLine2d
    {
    public:
        /**
         * @brief Constructor from two points. Note that the line is infinitely long.
         * @param pnt1 First point on the line.
         * @param pnt2 Second point on the line.
         */
        ImplicitLine2d(const Eigen::Vector2d& pnt1, const Eigen::Vector2d& pnt2) noexcept;

        /**
         * @brief Evaluates the implicit line equation. This gives the Euclidean distance to the implicit line.
         * @param pnt Point at which to evaluate the implicit equation.
         * @return Signed distance.
         */
        [[nodiscard]] double distance(const Eigen::Vector2d& pnt) const noexcept;

        /**
         * @brief Computes the distance function from P(t) to the line Q and converts this distance into Bezier Bernstein basis.
         * @param controlPoints Control points of a cubic Bezier curve P(t) to compute the distance to.
         * @return Control points of the cubic distance Bezier curve.
         */
        [[nodiscard]] Eigen::Vector4d distance(const std::array<Eigen::Vector2d, 4>& controlPoints) const noexcept;

        /**
         * @brief Coefficient of the implicit line: a*x + b*y + c = 0,  s.t. a^2+b^2=1
         */
        double a, b, c;
    };
}
