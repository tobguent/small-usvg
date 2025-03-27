#pragma once

#include "bezier_curve.hpp"

namespace usvg
{
    /**
     * @brief Class for tensor product surfaces with Bernstein basis functions.
     * @tparam TValue Type of values stored at the control points.
     * @tparam TN Degree of the surface in u direction.
     * @tparam TM Degree of the surface in v direction.
     */
    template <int TN, int TM, typename TValue>
    class BezierSurface
    {
    public:
        /**
         * @brief Degree of the curve in u direction.
         */
        static constexpr int N = TN;

        /**
         * @brief Degree of the curve in v direction.
         */
        static constexpr int M = TM;

        /**
         * @brief Type of values stored on the surface.
         */
        using Value = TValue;

        /**
         * @brief Coordinate in uv space.
         */
        using DomainCoord = Eigen::Vector2d;

        /**
         * @brief Constructor with zero initialization.
         */
        BezierSurface() noexcept
            : domain(Eigen::Vector2d::Zero(), Eigen::Vector2d::Ones())
        {
            for (int j = 0; j <= M; ++j)
                for (int i = 0; i <= N; ++i)
                    controlPoints[i][j] = TValue::Zero();
        }

        /**
         * @brief Constructor with initial weights.
         * @param initial Initial set of control points.
         */
        explicit BezierSurface(const std::array<std::array<TValue, M + 1>, N + 1>& initial) noexcept
            : controlPoints(initial)
            , domain(Eigen::Vector2d::Zero(), Eigen::Vector2d::Ones())
        {
        }

        /**
         * @brief Evaluates the tensor-product surface algebraically at a given uv coordinate.
         * @param uv Domain location where to evaluate the surface.
         * @return Position on the surface.
         */
        [[nodiscard]] TValue sample(const DomainCoord& uv) const noexcept
        {
            DomainCoord trel = (uv - this->domain.min()).cwiseQuotient(this->domain.max() - this->domain.min());
            TValue result    = TValue::Zero();
            for (int j = 0; j <= M; ++j)
                for (int i = 0; i <= N; ++i)
                    result += controlPoints[i][j] * BernsteinBasis::sample(trel.x(), i, N) * BernsteinBasis::sample(trel.y(), j, M);
            return result;
        }

        /**
         * @brief Evaluates the u-partial derivative of a tensor-product surface algebraically at a given uv coordinate.
         * @param uv Domain location where to evaluate the surface.
         * @return First u-partial derivative.
         */
        [[nodiscard]] TValue sample_du(const DomainCoord& uv) const noexcept
        {
            DomainCoord trel = (uv - this->domain.min()).cwiseQuotient(this->domain.max() - this->domain.min());
            TValue result    = TValue::Zero();
            for (int j = 0; j <= M; ++j)
                for (int i = 0; i <= N; ++i)
                    result += controlPoints[i][j] * BernsteinBasis::sample_dt(trel.x(), i, N) * BernsteinBasis::sample(trel.y(), j, M);
            return result;
        }

        /**
         * @brief Evaluates the u-partial derivative of a tensor-product surface algebraically at a given uv coordinate.
         * @param uv Domain location where to evaluate the surface.
         * @return First v-partial derivative.
         */
        [[nodiscard]] TValue sample_dv(const DomainCoord& uv) const noexcept
        {
            DomainCoord trel = (uv - this->domain.min()).cwiseQuotient(this->domain.max() - this->domain.min());
            TValue result    = TValue::Zero();
            for (int j = 0; j <= M; ++j)
                for (int i = 0; i <= N; ++i)
                    result += controlPoints[i][j] * BernsteinBasis::sample(trel.x(), i, N) * BernsteinBasis::sample_dt(trel.y(), j, M);
            return result;
        }

        /**
         * @brief Evaluates the second-order u-partial derivative of a tensor-product surface algebraically at a given uv coordinate.
         * @param uv Domain location where to evaluate the surface.
         * @return Second u-partial derivative.
         */
        [[nodiscard]] TValue sample_duu(const DomainCoord& uv) const noexcept
        {
            DomainCoord trel = (uv - this->domain.min()).cwiseQuotient(this->domain.max() - this->domain.min());
            TValue result    = TValue::Zero();
            for (int j = 0; j <= M; ++j)
                for (int i = 0; i <= N; ++i)
                    result += controlPoints[i][j] * BernsteinBasis::sample_dtt(trel.x(), i, N) * BernsteinBasis::sample(trel.y(), j, M);
            return result;
        }

        /**
         * @brief Evaluates the second-order mixed-partial derivative of a tensor-product surface algebraically at a given uv coordinate.
         * @param uv Domain location where to evaluate the surface.
         * @return Mixed second u- and v-partial derivative.
         */
        [[nodiscard]] TValue sample_duv(const DomainCoord& uv) const noexcept
        {
            DomainCoord trel = (uv - this->domain.min()).cwiseQuotient(this->domain.max() - this->domain.min());
            TValue result    = TValue::Zero();
            for (int j = 0; j <= M; ++j)
                for (int i = 0; i <= N; ++i)
                    result += controlPoints[i][j] * BernsteinBasis::sample_dt(trel.x(), i, N) * BernsteinBasis::sample_dt(trel.y(), j, M);
            return result;
        }

        /**
         * @brief Evaluates the second-order v-partial derivative of a tensor-product surface algebraically at a given uv coordinate.
         * @param uv Domain location where to evaluate the surface.
         * @return Second v-partial derivative.
         */
        [[nodiscard]] TValue sample_dvv(const DomainCoord& uv) const noexcept
        {
            DomainCoord trel = (uv - this->domain.min()).cwiseQuotient(this->domain.max() - this->domain.min());
            TValue result    = TValue::Zero();
            for (int j = 0; j <= M; ++j)
                for (int i = 0; i <= N; ++i)
                    result += controlPoints[i][j] * BernsteinBasis::sample(trel.x(), i, N) * BernsteinBasis::sample_dtt(trel.y(), j, M);
            return result;
        }

        /**
         * @brief Evaluates the Laplacian with respect to the patch parameterization at a given uv coordinated.
         * @param uv Domain location where to evaluate the surface.
         * @return Laplacian of the surface.
         */
        [[nodiscard]] Value sampleLaplacian(const DomainCoord& uv) const noexcept
        {
            return sample_duu(uv) + sample_dvv(uv);
        }

        /**
         * @brief Subdivides a Bezier tensor product surface at a specified domain coordinate into four patches.
         * @param output00 First part of the split surface.
         * @param output01 Second part of the split surface.
         * @param output10 Third part of the split surface.
         * @param output11 Fourth part of the split surface.
         * @param uv Parameter location where to split.
         */
        void subdivide(
            BezierSurface<TN, TM, TValue>& output00,
            BezierSurface<TN, TM, TValue>& output01,
            BezierSurface<TN, TM, TValue>& output10,
            BezierSurface<TN, TM, TValue>& output11,
            typename BezierSurface<TN, TM, TValue>::DomainCoord uv = typename BezierSurface<TN, TM, TValue>::DomainCoord(0.5, 0.5)) const noexcept
        {
            double umin = this->domain.min().x(), umax = this->domain.max().x();
            double vmin = this->domain.min().y(), vmax = this->domain.max().y();
            Eigen::Vector2d uvrel(
                (uv.x() - umin) / (umax - umin),
                (uv.y() - vmin) / (vmax - vmin));
            for (int i = 0; i <= TN; ++i)
                for (int j = 0; j <= TM; ++j)
                {
                    output00.controlPoints[i][j] = intermediateControlPoint(0, 0, i, j, uvrel);
                    output01.controlPoints[i][j] = intermediateControlPoint(0, j, i, TM - j, uvrel);
                    output10.controlPoints[i][j] = intermediateControlPoint(i, 0, TN - i, j, uvrel);
                    output11.controlPoints[i][j] = intermediateControlPoint(i, j, TN - i, TM - j, uvrel);
                }
            Eigen::Vector2d c00 = this->domain.min();
            Eigen::Vector2d c11 = this->domain.max();

            output00.domain = Eigen::AlignedBox2d(
                Eigen::Vector2d(this->domain.min().x(), this->domain.min().y()),
                Eigen::Vector2d(uv.x(), uv.y()));
            output01.domain = Eigen::AlignedBox2d(
                Eigen::Vector2d(this->domain.min().x(), uv.y()),
                Eigen::Vector2d(uv.x(), this->domain.max().y()));
            output10.domain = Eigen::AlignedBox2d(
                Eigen::Vector2d(uv.x(), this->domain.min().y()),
                Eigen::Vector2d(this->domain.max().x(), uv.y()));
            output11.domain = Eigen::AlignedBox2d(
                Eigen::Vector2d(uv.x(), uv.y()),
                Eigen::Vector2d(this->domain.max().x(), this->domain.max().y()));

            output00.recomputeBoundingBox();
            output01.recomputeBoundingBox();
            output10.recomputeBoundingBox();
            output11.recomputeBoundingBox();
        }

        /**
         * @brief Subdivides a Bezier tensor product surface at a specified u domain coordinate into two patches.
         * @param output0 First part of the split surface.
         * @param output1 Second part of the split surface.
         * @param u Parameter location where to split.
         */
        void subdivide_u(
            BezierSurface<TN, TM, TValue>& output0,
            BezierSurface<TN, TM, TValue>& output1,
            typename BezierCurve<TM, TValue>::DomainCoord u = typename BezierCurve<TM, TValue>::DomainCoord(0.5)) const noexcept
        {
            double umin = this->domain.min().x(), umax = this->domain.max().x();
            double urel = (u - umin) / (umax - umin);
            for (int i = 0; i <= TN; ++i)
                for (int j = 0; j <= TM; ++j)
                {
                    output0.controlPoints[i][j] = intermediateControlPoint(0, 0, i, j, Eigen::Vector2d(urel, 1.));
                    output1.controlPoints[i][j] = intermediateControlPoint(i, 0, TN - i, j, Eigen::Vector2d(urel, 1.));
                }
            output0.domain = Eigen::AlignedBox2d(
                Eigen::Vector2d(this->domain.min().x(), this->domain.min().y()),
                Eigen::Vector2d(u, this->domain.max().y()));
            output1.domain = Eigen::AlignedBox2d(
                Eigen::Vector2d(u, this->domain.min().y()),
                Eigen::Vector2d(this->domain.max().x(), this->domain.max().y()));
            output0.recomputeBoundingBox();
            output1.recomputeBoundingBox();
        }

        /**
         * @brief Subdivides a Bezier tensor product surface at a specified v domain coordinate into two patches.
         * @param output0 First part of the split surface.
         * @param output1 Second part of the split surface.
         * @param v Parameter location where to split.
         */
        void subdivide_v(
            BezierSurface<TN, TM, TValue>& output0,
            BezierSurface<TN, TM, TValue>& output1,
            typename BezierCurve<TN, TValue>::DomainCoord v = typename BezierCurve<TN, TValue>::DomainCoord(0.5)) const noexcept
        {
            double vmin = this->domain.min().y(), vmax = this->domain.max().y();
            double vrel = (v - vmin) / (vmax - vmin);
            for (int i = 0; i <= TN; ++i)
                for (int j = 0; j <= TM; ++j)
                {
                    output0.controlPoints[i][j] = intermediateControlPoint(0, 0, i, j, Eigen::Vector2d(1., vrel));
                    output1.controlPoints[i][j] = intermediateControlPoint(0, j, i, TM - j, Eigen::Vector2d(1., vrel));
                }
            output0.domain = Eigen::AlignedBox2d(
                Eigen::Vector2d(this->domain.min().x(), this->domain.min().y()),
                Eigen::Vector2d(this->domain.max().x(), v));
            output1.domain = Eigen::AlignedBox2d(
                Eigen::Vector2d(this->domain.min().x(), v),
                Eigen::Vector2d(this->domain.max().x(), this->domain.max().y()));
            output0.recomputeBoundingBox();
            output1.recomputeBoundingBox();
        }

        /**
         * @brief Gets the boundary curve along u (from 1->0) for v=0.
         * @param curve Resulting Bezier curve.
         */
        void getBoundaryCurveSouth(BezierCurve<TN, TValue>& curve) const noexcept
        {
            for (int i = 0; i <= TN; ++i)
                curve.controlPoints[TN - i] = this->controlPoints[i][0];
            curve.domain = Eigen::AlignedBox1d(this->domain.min().x(), this->domain.max().x());
            curve.recomputeBoundingBox();
        }

        /**
         * @brief Gets the boundary curve along u (from 0->1) for v=1.
         * @param curve Resulting Bezier curve.
         */
        void getBoundaryCurveNorth(BezierCurve<TN, TValue>& curve) const noexcept
        {
            for (int i = 0; i <= TN; ++i)
                curve.controlPoints[i] = this->controlPoints[i][TM];
            curve.domain = Eigen::AlignedBox1d(this->domain.min().x(), this->domain.max().x());
            curve.recomputeBoundingBox();
        }

        /**
         * @brief Gets the boundary curve along v (from 0->1) for u=0.
         * @param curve Resulting Bezier curve.
         */
        void getBoundaryCurveWest(BezierCurve<TM, TValue>& curve) const noexcept
        {
            for (int j = 0; j <= TM; ++j)
                curve.controlPoints[j] = this->controlPoints[0][j];
            curve.domain = Eigen::AlignedBox1d(this->domain.min().y(), this->domain.max().y());
            curve.recomputeBoundingBox();
        }

        /**
         * @brief Gets the boundary curve along v (from 1->0) for u=1.
         * @param curve Resulting Bezier curve.
         */
        void getBoundaryCurveEast(BezierCurve<TM, TValue>& curve) const noexcept
        {
            for (int j = 0; j <= TM; ++j)
                curve.controlPoints[TM - j] = this->controlPoints[TM][j];
            curve.domain = Eigen::AlignedBox1d(this->domain.min().y(), this->domain.max().y());
            curve.recomputeBoundingBox();
        }

        /**
         * @brief Recomputes the bounding box of the Bezier surface.
         */
        void recomputeBoundingBox() noexcept
        {
            mBoundingBox.setEmpty();
            for (int i = 0; i <= TN; ++i)
                for (int j = 0; j <= TM; ++j)
                    mBoundingBox.extend(this->controlPoints[i][j]);
        }

        /**
         * @brief Gets the bounding box of the Bezier surface.
         * @return Bounding box of this geometry.
         */
        [[nodiscard]] const Eigen::AlignedBox<double, TValue::SizeAtCompileTime>& getBoundingBox() const noexcept { return mBoundingBox; }

        /**
         * @brief Checks if the curve is properly initialized.
         * @return True if the bounding box is valid.
         */
        [[nodiscard]] bool isValid() const
        {
            // bounding box computed?
            if (mBoundingBox.isEmpty())
                return false;

            // is the parameter range not collapsing to a point?
            if (domain.min().x() == domain.max().x() || domain.min().y() == domain.max().y() || domain.isEmpty())
                return false;

            // all good
            return true;
        }

        /**
         * @brief Two-dimensional array of control points.
         */
        std::array<std::array<TValue, M + 1>, N + 1> controlPoints;

        /**
         * @brief Parameter domain [u,v] over which the surface is defined.
         */
        Eigen::AlignedBox2d domain;

    private:
        /**
         * @brief Calculates the intermediate control point of the de Casteljau algorithm.
         * @param i Index of the control point in u-direction.
         * @param j Index of the control point in v-direction.
         * @param r Degree of the intermediate point in u-direction.
         * @param s Degree of the intermediate point in v-direction.
         * @param uv Parameter location where to evaluate.
         * @return Position of the intermediate control point.
         */
        [[nodiscard]] TValue intermediateControlPoint(int i, int j, int r, int s, const typename BezierSurface<TN, TM, TValue>::DomainCoord& uv) const noexcept
        {
            TValue value = TValue::Zero();
            for (int ii = 0; ii <= r; ++ii)
                for (int jj = 0; jj <= s; ++jj)
                    value += BernsteinBasis::sample(uv.x(), ii, r) * BernsteinBasis::sample(uv.y(), jj, s) * this->controlPoints[i + ii][j + jj];
            return value;
        }

        /**
         * @brief Bounding box of the surface.
         */
        Eigen::AlignedBox<double, TValue::SizeAtCompileTime> mBoundingBox;
    };

    /**
     * @brief Bicubic tensor product surface with Bernstein basis functions.
     * @tparam TValue Coordinate type of the curve.
     */
    template <typename TValue>
    using BicubicBezierSurface = BezierSurface<3, 3, TValue>;

    /**
     * @brief Bicubic uni-variate tensor product surface with Bernstein basis functions.
     */
    using BicubicBezierSurface1d = BicubicBezierSurface<Eigen::Vector<double, 1>>;

    /**
     * @brief Bicubic bi-variate tensor product surface with Bernstein basis functions.
     */
    using BicubicBezierSurface2d = BicubicBezierSurface<Eigen::Vector2d>;

    /**
     * @brief Bicubic tri-variate tensor product surface with Bernstein basis functions.
     */
    using BicubicBezierSurface3d = BicubicBezierSurface<Eigen::Vector3d>;

    /**
     * @brief Bicubic tetra-variate tensor product surface with Bernstein basis functions.
     */
    using BicubicBezierSurface4d = BicubicBezierSurface<Eigen::Vector4d>;
}
