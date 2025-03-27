#pragma once

#include "types.hpp"

namespace usvg
{
    /**
     * @brief Base class for parametric curves.
     * @tparam TValue Type of values stored on the curve.
     */
    template <typename TValue>
    class IParametricCurve
    {
    public:
        /**
         * @brief Coordinate in space.
         */
        using DomainCoord = double;

        /**
         * @brief Type of values stored on the curve.
         */
        using Value = TValue;

        /**
         * @brief Default constructor for a parameteric curve with the parameter range t \in [0,1].
         */
        IParametricCurve() noexcept
            : domain(0, 1)
        {
        }

        /**
         * @brief Evaluates the parametric curve at a given t coordinate.
         * @param t Domain location where to sample the curve.
         * @return Position of the curve.
         */
        [[nodiscard]] virtual Value sample(const DomainCoord& t) const noexcept = 0;

        /**
         * @brief Evaluates the first-order derivative at a given t coordinate.
         * @param t Domain location where to sample the curve.
         * @return Tangent of the curve.
         */
        [[nodiscard]] virtual Value sample_dt(const DomainCoord& t) const noexcept = 0;

        /**
         * @brief Evaluates the second-order derivative at a given t coordinate.
         * @param t Domain location where to sample the curve.
         * @return Acceleration of the curve.
         */
        [[nodiscard]] virtual Value sample_dtt(const DomainCoord& t) const noexcept = 0;

        /**
         * @brief Checks if the curve is properly initialized.
         * @return True if valid.
         */
        [[nodiscard]] bool isValid() const
        {
            // is the parameter range not collapsing to a point?
            if (domain.min().x() == domain.max().x() || domain.isEmpty())
                return false;

            // all good
            return true;
        }

        /**
         * @brief Parameter domain over which the curve is defined.
         */
        Eigen::AlignedBox1d domain;
    };

    /**
     * @brief Base class for 1-variate parametric curves.
     */
    using IParametricCurve1d = IParametricCurve<Eigen::Vector1d>;

    /**
     * @brief Base class for 2-variate parametric curves.
     */
    using IParametricCurve2d = IParametricCurve<Eigen::Vector2d>;

    /**
     * @brief Base class for 3-variate parametric curves.
     */
    using IParametricCurve3d = IParametricCurve<Eigen::Vector3d>;

    /**
     * @brief Base class for 4-variate parametric curves.
     */
    using IParametricCurve4d = IParametricCurve<Eigen::Vector4d>;
}
