#pragma once

#include "array.hpp"
#include "iparametric_curve.hpp"

namespace usvg
{
    /**
     * @brief Base class for parameterizable curves.
     * @tparam TValue Type of values stored along the curve.
     */
    template <typename TValue>
    class IParameterizableCurve : public IParametricCurve<TValue>
    {
    public:
        /**
         * @brief Coordinate in space.
         */
        using DomainCoord = double;

        /**
         * @brief Constructor.
         */
        IParameterizableCurve()
            : uniform(true)
        {
        }

        /**
         * @brief Locates a given domain coordinate and returns information for interpolation. If the sample point is outside the domain, then the last segment is returned and the interpolation factor is outside of [0,1] leading to an extrapolation of the last segment.
         * @param domainCoord Domain coordinate to look for.
         * @return Lower index in array, upper index in array, interpolation factor in [0,1] between lower and upper, and the parameterization distance between lower and upper.
         */
        [[nodiscard]] std::tuple<Eigen::Index, Eigen::Index, double, double> locate(const DomainCoord& domainCoord) const noexcept
        {
            if (uniform)
                return locateUniform(domainCoord);
            else
                return locateExplicit(domainCoord);
        }

        /**
         * @brief Checks if the parameterization is properly initialized.
         * @return True if valid.
         */
        [[nodiscard]] bool isValid() const
        {
            if (!IParametricCurve<TValue>::isValid())
                return false;

            // if not uniform
            if (!uniform)
            {
                // parameterization monotonically increasing?
                for (Eigen::Index i = 0; i < parameters.getSize() - 1; ++i)
                    if (parameters.getValue(i).x() > parameters.getValue(i + 1).x())
                        return false;

                // does the parameterization start and end at the domain boundary?
                if (parameters.getSize() == 0 || parameters.first().x() != this->domain.min().x() || parameters.last().x() != this->domain.max().x())
                    return false;
            }
            else
            {
                // parameterization actually uniform?
                for (Eigen::Index i = 0; i < parameters.getSize(); ++i)
                    if (std::abs(parameters.getValue(i).x() - i - this->domain.min().x()) >= 1E-5)
                        return false;
            }

            // all good
            return true;
        }

        /**
         * @brief Sets a flag to use explicit parameterization from now on. Sets the domain according to the first and last parameter value. Thus, make sure to set the parameters before calling this function.
         */
        void setExplicitParameterization() noexcept
        {
            uniform = false;
            if (parameters.getSize() != 0)
            {
                this->domain = Eigen::AlignedBox1d(parameters.first(), parameters.last());
            }
        }

        /**
         * @brief Is the parameterization uniform?
         * @return True, if the parameterization is uniform.
         */
        [[nodiscard]] bool isUniform() const { return uniform; }

        /**
         * @brief Array containing an explicit parameterization.
         */
        Array1d parameters;

    protected:
        /**
         * @brief Sets a flag to use uniform parameterization from now on.
         * @tparam TArray Type of the array, which stores the positions.
         * @param positions Position array from which the parameterization is derived.
         */
        template <typename TArray>
        void setUniformParameterization(const TArray& positions, Eigen::Index stride) noexcept
        {
            uniform      = true;
            double param = this->domain.min().x();
            parameters.setSize((positions.getSize() - 1) / stride + 1);
            for (Eigen::Index i = 0; i < parameters.getSize(); ++i)
                parameters.setValue(i, Eigen::Vector1d(param + i));
            this->domain = Eigen::AlignedBox1d(this->domain.min().x(), parameters.last());
        }

    private:
        /**
         * @brief Is uniform parameterization
         */
        bool uniform;

        /**
         * @brief Locates a given domain coordinate assuming uniform parameterization and returns information for interpolation. If the sample point is outside the domain, then the last segment is returned and the interpolation factor is outside of [0,1] leading to an extrapolation of the last segment.
         * @param domainCoord Domain coordinate to look for.
         * @return Lower index in array, upper index in array, interpolation factor in [0,1] between lower and upper, and the parameterization distance between lower and upper.
         */
        [[nodiscard]] std::tuple<Eigen::Index, Eigen::Index, double, double> locateUniform(const DomainCoord& domainCoord) const noexcept
        {
            Eigen::Index lower = (Eigen::Index)(domainCoord - this->domain.min().x());
            // clamp lower to [0, N-2]
            lower = std::min(std::max(Eigen::Index(0), lower), parameters.getSize() - 2);
            // get upper in [1, N-1]
            Eigen::Index upper = lower + 1;
            // t value might be outside of [0,1] range, which would lead to extrapolation
            double t = (domainCoord - this->domain.min().x()) - lower;
            return std::make_tuple(lower, upper, t, 1.);
        }

        /**
         * @brief Locates a given domain coordinate assuming explicit parameterization and returns information for interpolation. If the sample point is outside the domain, then the last segment is returned and the interpolation factor is outside of [0,1] leading to an extrapolation of the last segment.
         * @param domainCoord Domain coordinate to look for.
         * @return Lower index in array, upper index in array, interpolation factor in [0,1] between lower and upper if in range (otherwise it extrapolates), and the parameterization distance between lower and upper.
         */
        [[nodiscard]] std::tuple<Eigen::Index, Eigen::Index, double, double> locateExplicit(const DomainCoord& domainCoord) const noexcept
        {
            Eigen::Index L = 0;
            Eigen::Index R = parameters.getSize() - 1;
            while (L < R)
            {
                Eigen::Index m = (L + R) / 2;
                double h       = parameters.getValue(m).x();
                if (h < domainCoord)
                    L = m + 1;
                else
                    R = m;
            }
            Eigen::Index lower = std::max(Eigen::Index(0), L - 1);
            Eigen::Index upper = std::min(lower + 1, parameters.getSize() - 1);
            double h0          = parameters.getValue(lower).x();
            double h1          = parameters.getValue(upper).x();
            if (h0 > domainCoord || domainCoord > h1)
            {
                int test = 0;
                ++test;
            }
            assert(h0 <= domainCoord && domainCoord <= h1);
            double rat = h1 == h0 ? 0 : ((domainCoord - h0) / (h1 - h0));
            // rat        = std::min(std::max(0., rat), 1.);
            return std::make_tuple(lower, upper, rat, h1 - h0);
        }
    };

    /**
     * @brief Base class for 1-variate parameterizable curves.
     */
    using IParameterizableCurve1d = IParameterizableCurve<Eigen::Vector1d>;

    /**
     * @brief Base class for 2-variate parameterizable curves.
     */
    using IParameterizableCurve2d = IParameterizableCurve<Eigen::Vector2d>;

    /**
     * @brief Base class for 3-variate parameterizable curves.
     */
    using IParameterizableCurve3d = IParameterizableCurve<Eigen::Vector3d>;

    /**
     * @brief Base class for 4-variate parameterizable curves.
     */
    using IParameterizableCurve4d = IParameterizableCurve<Eigen::Vector4d>;
}
