#pragma once

#include "array.hpp"
#include "bernstein_basis.hpp"
#include "bezier_curve.hpp"
#include "iparameterizable_curve.hpp"

namespace usvg
{
    /**
     * @brief Class that stores and manipulates bezier splines of general degree, having at least C0 continuity.
     * @tparam TDegree Degree of the individual Bezier curve.
     * @tparam TValue Type of values stored at the control points.
     */
    template <int TDegree, typename TValue>
    class BezierSpline : public IParameterizableCurve<TValue>
    {
    public:
        /**
         * @brief Degree of the curve.
         */
        static constexpr int Degree = TDegree;

        /**
         * @brief Type of values stored on the curve.
         */
        using Value = TValue;

        /**
         * @brief Coordinate in space.
         */
        using DomainCoord = typename IParametricCurve<TValue>::DomainCoord;

        /**
         * @brief Evaluates the parametric curve algebraically at a given t coordinate.
         * @param t Domain location where to sample the curve.
         * @return Position of the curve.
         */
        [[nodiscard]] TValue sample(const DomainCoord& t) const noexcept override
        {
            // find the line segment
            auto [lower, upper, tvalue, h] = this->locate(t);
            // get entry index into control points
            int offset = lower * TDegree;

            // algebraically evaluate the Bezier curve
            TValue result = TValue::Zero();
            for (int i = 0; i <= Degree; ++i)
            {
                int cpIndex = std::min(std::max(0, offset + i), (int)controlPoints.getSize() - 1); // clamp to not accidentally run out
                result += controlPoints.getValue(cpIndex) * BernsteinBasis::sample(tvalue, i, Degree);
            }
            return result;
        }

        /**
         * @brief Evaluates the first-order derivative algebraically at a given t coordinate.
         * @param t Domain location where to sample the curve.
         * @return Tangent of the curve.
         */
        [[nodiscard]] TValue sample_dt(const DomainCoord& t) const noexcept override
        {
            // find the line segment
            auto [lower, upper, tvalue, h] = this->locate(t);
            // get entry index into control points
            int offset = lower * TDegree;
            // rate of change to unit parameterization
            double dt_ds = (this->parameters.getValue(upper).x() - this->parameters.getValue(lower).x());

            TValue result = TValue::Zero();
            for (int i = 0; i <= Degree; ++i)
                result += controlPoints.getValue(offset + i) * BernsteinBasis::sample_dt(tvalue, i, Degree) / dt_ds;
            return result;
        }

        /**
         * @brief Evaluates the second-order derivative algebraically at a given t coordinate.
         * @param t Domain location where to sample the curve.
         * @return Acceleration of the curve.
         */
        [[nodiscard]] TValue sample_dtt(const DomainCoord& t) const noexcept override
        {
            // find the line segment
            auto [lower, upper, tvalue, h] = this->locate(t);
            // get entry index into control points
            int offset = lower * TDegree;
            // rate of change to unit parameterization
            double dt_ds = (this->parameters.getValue(upper).x() - this->parameters.getValue(lower).x());

            TValue result = TValue::Zero();
            for (int i = 0; i <= Degree; ++i)
                result += controlPoints.getValue(offset + i) * BernsteinBasis::sample_dtt(tvalue, i, Degree) / (dt_ds * dt_ds);
            return result;
        }

        /**
         * @brief Subdivides a Bezier spline into two splines at a specified parameter value t.
         * @param spline1 First part of the split spline.
         * @param spline2 Second part of the split spline.
         * @param t Parameter location where to split.
         */
        void subdivide(
            BezierSpline<TDegree, TValue>& spline1,
            BezierSpline<TDegree, TValue>& spline2,
            DomainCoord domainCoord = DomainCoord(0.5)) const noexcept
        {
            //        start      end
            // |        |         |        |
            // 0  1  2  3  4 | 5  6  7  8  9
            // a0       a1        b0       b1
            // --------------------            spline1
            //          --------------------   spline2
            auto [startIndex, endIndex, t, h] = this->locate(domainCoord);

            int a0 = 0;
            int a1 = startIndex * TDegree;
            int b0 = endIndex * TDegree;
            int b1 = this->controlPoints.getSize() - 1;

            // get the curve that needs to split
            auto midCurve = this->getCurve(startIndex);
            BezierCurve<TDegree, TValue> part1, part2;
            midCurve.subdivide(part1, part2, midCurve.domain.min()[0] + (midCurve.domain.max()[0] - midCurve.domain.min()[0]) * t);

            // copy the first half into spline1
            spline1.controlPoints.setSize(b0 - a0 + 1);
            spline1.parameters.setSize(endIndex + 1);
            for (Eigen::Index i = a0; i <= a1; ++i)
                spline1.controlPoints.setValue(i, this->controlPoints.getValue(i));
            for (Eigen::Index i = 0; i <= startIndex; ++i)
                spline1.parameters.setValue(i, this->parameters.getValue(i));
            for (int d = 0; d <= TDegree; ++d)
                spline1.controlPoints.setValue(a1 + d, part1.controlPoints[d]);
            spline1.parameters.setValue(startIndex + 1, Eigen::Vector1d(domainCoord));
            spline1.domain = Eigen::AlignedBox1d(this->domain.min()[0], domainCoord);
            // remove the last points if it is a duplicate (i.e., we split exactly on the control point)
            if (spline1.controlPoints.getSize() > TDegree)
            {
                bool isequal = std::abs(spline1.parameters.getValue(spline1.parameters.getSize() - 2)[0] - spline1.parameters.getValue(spline1.parameters.getSize() - 1)[0]) < 1E-10;
                for (int d = 0; d < TDegree; ++d)
                {
                    if ((spline1.controlPoints.getValue(spline1.controlPoints.getSize() - 1 - d) - spline1.controlPoints.getValue(spline1.controlPoints.getSize() - 1)).stableNorm() > 1E-10)
                    {
                        isequal = false;
                        break;
                    }
                }
                if (isequal)
                {
                    spline1.controlPoints.removeLast(TDegree);
                    spline1.parameters.removeLast();
                }
            }
            spline1.setExplicitParameterization();
            spline1.recomputeBoundingBox();

            // copy the second half into spline2
            int N = this->getNumCurves();
            spline2.controlPoints.setSize(b1 - a1 + 1);
            spline2.parameters.setSize(N - startIndex + 1);
            for (Eigen::Index i = b0; i <= b1; ++i)
                spline2.controlPoints.setValue(i - b0 + TDegree, this->controlPoints.getValue(i));
            for (Eigen::Index i = endIndex; i <= N; ++i)
                spline2.parameters.setValue(i - endIndex + 1, this->parameters.getValue(i));
            for (int d = 0; d <= TDegree; ++d)
                spline2.controlPoints.setValue(d, part2.controlPoints[d]);
            spline2.parameters.setValue(0, Eigen::Vector1d(domainCoord));
            // remove the first points if it is a duplicate (i.e., we split exactly on the control point)
            if (spline2.controlPoints.getSize() > TDegree)
            {
                bool isequal = std::abs(spline2.parameters.getValue(0)[0] - spline2.parameters.getValue(1)[0]) < 1E-10;
                for (int d = 0; d < TDegree; ++d)
                {
                    if ((spline2.controlPoints.getValue(0) - spline2.controlPoints.getValue(d + 1)).stableNorm() > 1E-10)
                    {
                        isequal = false;
                        break;
                    }
                }
                if (isequal)
                {
                    spline2.controlPoints.removeFirst(TDegree);
                    spline2.parameters.removeFirst();
                }
            }
            spline2.setExplicitParameterization();
            spline2.domain = Eigen::AlignedBox1d(domainCoord, this->domain.max()[0]);
            spline2.recomputeBoundingBox();
        }

        /**
         * @brief Checks if the curve is properly initialized.
         * @return True if the bounding box is valid.
         */
        bool isValid() const
        {
            if (!IParameterizableCurve<TValue>::isValid())
                return false;

            // is parameterization valid?
            if (this->parameters.getSize() != getNumCurves() + 1)
                return false;

            // bounding box computed?
            if (this->mBoundingBox.isEmpty())
                return false;

            // no dangling control points?
            if (((controlPoints.getSize() - 1) % TDegree) != 0)
                return false;

            // all good
            return true;
        }

        /**
         * @brief Recomputes the bounding box of the Bezier curve.
         */
        void recomputeBoundingBox() noexcept
        {
            mBoundingBox.setEmpty();
            for (int i = 0; i < controlPoints.getSize(); ++i)
                mBoundingBox.extend(controlPoints.getValue(i));
        }

        /**
         * @brief Gets the bounding box of the vertices. Note that recomputeBoundingBox has to be called first!
         * @return Bounding box of this geometry.
         */
        [[nodiscard]] inline const Eigen::AlignedBox<double, TValue::SizeAtCompileTime>& getBoundingBox() const noexcept { return mBoundingBox; }

        /**
         * @brief Returns the number Bezier curves segments in the spline.
         * @return Number of Bezier curves.
         */
        [[nodiscard]] int getNumCurves() const noexcept
        {
            return (controlPoints.getSize() - 1) / TDegree;
        }

        /**
         * @brief Gets the i'th Bezier curve of the spline.
         * @param i Index of Bezier curve to get.
         * @return Bezier curve.
         */
        [[nodiscard]] BezierCurve<TDegree, TValue> getCurve(int index) const
        {
            BezierCurve<TDegree, TValue> curve;
            for (int cp = 0; cp <= TDegree; ++cp)
                curve.controlPoints[cp] = controlPoints.getValue(index * TDegree + cp);
            curve.domain = Eigen::AlignedBox1d(this->parameters.getValue(index), this->parameters.getValue(index + 1));
            curve.recomputeBoundingBox();
            return curve;
        }

        /**
         * @brief Computes a uniform parameterization. The start of the parameter domain is retained. The end of the parameter domain is set accordingly.
         */
        void setUniformParameterization() noexcept
        {
            IParameterizableCurve<TValue>::setUniformParameterization(controlPoints, TDegree);
        }

        /**
         * @brief Array of Bezier control points.
         */
        Array<TValue> controlPoints;

    private:
        /**
         * @brief Bounding box of the curve.
         */
        Eigen::AlignedBox<double, TValue::SizeAtCompileTime> mBoundingBox;
    };

    /**
     * @brief Cubic bezier spline.
     * @tparam TValue Coordinate type of the curve.
     */
    template <typename TValue>
    using CubicBezierSpline = BezierSpline<3, TValue>;

    /**
     * @brief Cubic uni-variate spline with Bernstein basis functions.
     */
    using CubicBezierSpline1d = CubicBezierSpline<Eigen::Vector<double, 1>>;

    /**
     * @brief Cubic bi-variate spline with Bernstein basis functions.
     */
    using CubicBezierSpline2d = CubicBezierSpline<Eigen::Vector2d>;

    /**
     * @brief Cubic tri-variate spline with Bernstein basis functions.
     */
    using CubicBezierSpline3d = CubicBezierSpline<Eigen::Vector3d>;

    /**
     * @brief Cubic tetra-variate spline with Bernstein basis functions.
     */
    using CubicBezierSpline4d = CubicBezierSpline<Eigen::Vector4d>;
}
