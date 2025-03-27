#pragma once

#include "bernstein_basis.hpp"
#include "iparametric_curve.hpp"

namespace usvg
{
    /**
     * @brief Class that stores and manipulates bezier curves of general degree.
     * @tparam TDegree Degree of the curve.
     * @tparam TValue Type of values stored at the control points.
     */
    template <int TDegree, typename TValue>
    class BezierCurve : public IParametricCurve<TValue>
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
         * @brief Constructor with zero initialization.
         */
        BezierCurve() noexcept
        {
            for (int i = 0; i <= Degree; ++i)
                controlPoints[i] = TValue::Zero();
        }

        /**
         * @brief Constructor with initial weights.
         * @param initial Initial set of control points.
         */
        explicit BezierCurve(const std::array<TValue, Degree + 1>& initial) noexcept
            : controlPoints(initial)
        {
            recomputeBoundingBox();
        }

        /**
         * @brief Evaluates the parametric curve algebraically at a given t coordinate.
         * @param t Domain location where to sample the curve.
         * @return Position of the curve.
         */
        [[nodiscard]] TValue sample(const DomainCoord& t) const noexcept override
        {
            double tmin = this->domain.min().x(), tmax = this->domain.max().x();
            double trel   = (t - tmin) / (tmax - tmin);
            TValue result = TValue::Zero();
            for (int i = 0; i <= Degree; ++i)
                result += controlPoints[i] * BernsteinBasis::sample(trel, i, Degree);
            return result;
        }

        /**
         * @brief Evaluates the first-order derivative algebraically at a given t coordinate.
         * @param t Domain location where to sample the curve.
         * @return Tangent of the curve.
         */
        [[nodiscard]] TValue sample_dt(const DomainCoord& t) const noexcept override
        {
            double tmin = this->domain.min().x(), tmax = this->domain.max().x();
            double trel   = (t - tmin) / (tmax - tmin);
            TValue result = TValue::Zero();
            for (int i = 0; i <= Degree; ++i)
                result += controlPoints[i] * BernsteinBasis::sample_dt(trel, i, Degree);
            return result / (tmax - tmin);
        }

        /**
         * @brief Evaluates the second-order derivative algebraically at a given t coordinate.
         * @param t Domain location where to sample the curve.
         * @return Acceleration of the curve.
         */
        [[nodiscard]] TValue sample_dtt(const DomainCoord& t) const noexcept override
        {
            double tmin = this->domain.min().x(), tmax = this->domain.max().x();
            double trel   = (t - tmin) / (tmax - tmin);
            TValue result = TValue::Zero();
            for (int i = 0; i <= Degree; ++i)
                result += controlPoints[i] * BernsteinBasis::sample_dtt(trel, i, Degree);
            return result / std::pow(tmax - tmin, 2);
        }

        /**
         * @brief Subdivides a Bezier curve into two curves at a specified parameter value t.
         * @param output1 First part of the split curve.
         * @param output2 Second part of the split curve.
         * @param t Parameter location where to split.
         */
        void subdivide(
            BezierCurve<TDegree, TValue>& output1,
            BezierCurve<TDegree, TValue>& output2,
            typename BezierCurve<TDegree, TValue>::DomainCoord t = typename BezierCurve<TDegree, TValue>::DomainCoord(0.5)) const noexcept
        {
            double tmin = this->domain.min().x(), tmax = this->domain.max().x();
            double trel = (t - tmin) / (tmax - tmin);
            for (int n = 0; n <= TDegree; ++n)
            {
                output1.controlPoints[n] = intermediateControlPoint(0, n, trel);
                output2.controlPoints[n] = intermediateControlPoint(n, TDegree - n, trel);
            }
            output1.recomputeBoundingBox();
            output2.recomputeBoundingBox();
            output1.domain = Eigen::AlignedBox1d(Eigen::Vector1d(tmin), Eigen::Vector1d(t));
            output2.domain = Eigen::AlignedBox1d(Eigen::Vector1d(t), Eigen::Vector1d(tmax));
        }

        /**
         * @brief Recomputes the bounding box of the Bezier curve.
         */
        void recomputeBoundingBox() noexcept
        {
            mBoundingBox.setEmpty();
            for (int i = 0; i <= TDegree; ++i)
                mBoundingBox.extend(controlPoints[i]);
        }

        /**
         * @brief Gets the bounding box of the vertices. Note that recomputeBoundingBox has to be called first!
         * @return Bounding box of this geometry.
         */
        [[nodiscard]] inline const Eigen::AlignedBox<double, TValue::SizeAtCompileTime>& getBoundingBox() const noexcept { return mBoundingBox; }

        /**
         * @brief Checks if the curve is properly initialized.
         * @return True if valid.
         */
        [[nodiscard]] bool isValid() const
        {
            // check base class
            if (!IParametricCurve<TValue>::isValid())
                return false;

            // bounding box computed?
            if (this->mBoundingBox.isEmpty())
                return false;

            // all good
            return true;
        }

        /**
         * @brief Array of Degree+1 control points.
         */
        std::array<TValue, Degree + 1> controlPoints;

    private:
        /**
         * @brief Calculates the intermediate control point of the de Casteljau algorithm.
         * @param i Index of the control point.
         * @param r Degree of the intermediate point.
         * @param t Parameter location where to evaluate.
         * @return Position of the intermediate control point.
         */
        [[nodiscard]] TValue intermediateControlPoint(int i, int r, typename BezierCurve<TDegree, TValue>::DomainCoord t) const noexcept
        {
            TValue value = TValue::Zero();
            for (int j = 0; j <= r; ++j)
                value += BernsteinBasis::sample(t, j, r) * controlPoints[i + j];
            return value;
        }

        /**
         * @brief Bounding box of the curve.
         */
        Eigen::AlignedBox<double, TValue::SizeAtCompileTime> mBoundingBox;
    };

    /**
     * @brief Cubic bezier curve.
     * @tparam TValue Coordinate type of the curve.
     */
    template <typename TValue>
    using CubicBezierCurve = BezierCurve<3, TValue>;

    /**
     * @brief Cubic uni-variate polynomial with Bernstein basis functions.
     */
    using CubicBezierCurve1d = CubicBezierCurve<Eigen::Vector<double, 1>>;

    /**
     * @brief Cubic bi-variate polynomial with Bernstein basis functions.
     */
    using CubicBezierCurve2d = CubicBezierCurve<Eigen::Vector2d>;

    /**
     * @brief Cubic tri-variate polynomial with Bernstein basis functions.
     */
    using CubicBezierCurve3d = CubicBezierCurve<Eigen::Vector3d>;

    /**
     * @brief Cubic tetra-variate polynomial with Bernstein basis functions.
     */
    using CubicBezierCurve4d = CubicBezierCurve<Eigen::Vector4d>;
}
