#pragma once

#include "array.hpp"
#include "iparameterizable_curve.hpp"

namespace usvg
{
    /**
     * @brief Class for piecewise linear curves.
     * @tparam TArray Type of array that stores the values.
     */
    template <typename TArray>
    class PiecewiseLinearCurve : public IParameterizableCurve<typename TArray::Element>
    {
    public:
        /**
         * @brief Coordinate in parameter space.
         */
        using DomainCoord = double;

        /**
         * @brief Type of array containing the values.
         */
        using Array = TArray;

        /**
         * @brief Type of value a single value.
         */
        using Value = typename TArray::Element;

        /**
         * @brief Evaluates the value of a curve using linear interpolation at a given domain location.
         * @param domainCoord Domain location where to sample the curve.
         * @return Value of the curve.
         */
        [[nodiscard]] Value sample(const DomainCoord& domainCoord) const noexcept override
        {
            auto [lower, upper, t, h] = this->locate(domainCoord);
            return (1 - t) * values.getValue(lower) + t * values.getValue(upper);
        }

        /**
         * @brief Evaluates the first-order derivative of the values at a given domain coordinate.
         * @param t Domain location where to sample the curve.
         * @return First derivative of the value.
         */
        [[nodiscard]] Value sample_dt(const DomainCoord& domainCoord) const noexcept override
        {
            auto [lower, upper, t, h] = this->locate(domainCoord);
            return (values.getValue(upper) - values.getValue(lower)) / h;
        }

        /**
         * @brief Evaluates the second-order derivative of the values at a given domain coordinate. Note that this is always zero, since the curve is piecewise linear.
         * @param t Domain location where to sample the curve.
         * @return Second derivative of the value.
         */
        [[nodiscard]] Value sample_dtt(const DomainCoord& domainCoord) const noexcept override
        {
            return Value::Zero();
        }

        /**
         * @brief Subdivides the piecewise linear curve into two pieces.
         * @param curve1 First part of the curve.
         * @param curve2 Second part of the curve.
         * @param domainCoord Domain coordinate where to split the input curve.
         */
        void subdivide(PiecewiseLinearCurve<TArray>& curve1, PiecewiseLinearCurve<TArray>& curve2, DomainCoord domainCoord) const
        {
            auto [startIndex, endIndex, t, h] = this->locate(domainCoord);

            // copy the first half into curve1
            curve1.values.setSize(startIndex + 2);
            curve1.parameters.setSize(startIndex + 2);
            for (Eigen::Index i = 0; i <= startIndex; ++i)
            {
                curve1.values.setValue(i, this->values.getValue(i));
                curve1.parameters.setValue(i, this->parameters.getValue(i));
            }
            Value interpVal             = this->values.getValue(startIndex) * (1 - t) + this->values.getValue(endIndex) * t;
            Eigen::Vector1d interpParam = this->parameters.getValue(startIndex) * (1 - t) + this->parameters.getValue(endIndex) * t;
            curve1.values.setValue(startIndex + 1, interpVal);
            curve1.parameters.setValue(startIndex + 1, interpParam);
            curve1.domain = Eigen::AlignedBox1d(this->domain.min()[0], domainCoord);
            curve1.setExplicitParameterization();
            curve1.recomputeBoundingBox();

            // copy the second half into curve2
            int N = this->values.getSize();
            curve2.values.setSize(N - endIndex + 1);
            curve2.parameters.setSize(N - endIndex + 1);
            for (Eigen::Index i = 0; i < N - endIndex; ++i)
            {
                curve2.values.setValue(i + 1, this->values.getValue(endIndex + i));
                curve2.parameters.setValue(i + 1, this->parameters.getValue(endIndex + i));
            }
            curve2.values.setValue(0, interpVal);
            curve2.parameters.setValue(0, interpParam);
            curve2.setExplicitParameterization();
            curve2.domain = Eigen::AlignedBox1d(domainCoord, this->domain.max()[0]);
            curve2.recomputeBoundingBox();
        }

        /**
         * @brief Checks if the curve is properly initialized.
         * @return True if valid.
         */
        [[nodiscard]] bool isValid() const
        {
            if (!IParameterizableCurve<typename TArray::Element>::isValid())
                return false;

            // is parameterization valid?
            if (values.getSize() != this->parameters.getSize())
                return false;

            // bounding box computed?
            if (this->mBoundingBox.isEmpty())
                return false;

            // all good
            return true;
        }

        /**
         * @brief Recomputes the bounding box of the values.
         */
        void recomputeBoundingBox() noexcept
        {
            mBoundingBox.setEmpty();
            for (int i = 0; i < values.getSize(); ++i)
                mBoundingBox.extend(values.getValue(i));
        }

        /**
         * @brief Gets the bounding box of the values. Note that recomputeBoundingBox has to be called first!
         * @return Bounding box of the values.
         */
        [[nodiscard]] inline const Eigen::AlignedBox<double, TArray::Element::RowsAtCompileTime>& getBoundingBox() const noexcept { return mBoundingBox; }

        /**
         * @brief Computes a uniform parameterization. The start of the parameter domain is retained. The end of the parameter domain is set accordingly.
         */
        void setUniformParameterization() noexcept
        {
            IParameterizableCurve<typename TArray::Element>::setUniformParameterization(values, 1);
        }

        /**
         * @brief Computes a Chordal parameterization from the given positions. The start of the parameter domain is retained. The end of the parameter domain is set according to the total arc length of the curve.
         */
        void setChordalParameterization() noexcept
        {
            IParameterizableCurve<typename TArray::Element>::setChordalParameterization(values, 1);
        }

        /**
         * @brief Array containing the values.
         */
        Array values;

    private:
        /**
         * @brief Bounding box of the curve.
         */
        Eigen::AlignedBox<double, Value::RowsAtCompileTime> mBoundingBox;
    };

    /**
     * @brief Uni-variate piecewise linear curve.
     */
    using PiecewiseLinearCurve1d = PiecewiseLinearCurve<Array1d>;

    /**
     * @brief Bi-variate piecewise linear curve.
     */
    using PiecewiseLinearCurve2d = PiecewiseLinearCurve<Array2d>;

    /**
     * @brief Tri-variate piecewise linear curve.
     */
    using PiecewiseLinearCurve3d = PiecewiseLinearCurve<Array3d>;

    /**
     * @brief Tetra-variate piecewise linear curve.
     */
    using PiecewiseLinearCurve4d = PiecewiseLinearCurve<Array4d>;
}
