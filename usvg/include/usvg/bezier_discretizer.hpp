#pragma once

#include "bezier_curve.hpp"
#include "bezier_spline.hpp"
#include "implicit_line.hpp"
#include "piecewise_linear_curve.hpp"

namespace usvg
{
    /**
     * @brief Class that discretizes a cubic Bezier curve into a polyline.
     */
    class BezierDiscretizer2d
    {
    public:
        /**
         * @brief Computes the largest distance from any point along the curve to the straight line that connects its start and end point.
         * @param curve Curve to compute distance for.
         * @return Largest distance and the curve parameter where the largest distance occured.
         */
        [[nodiscard]] static std::pair<double, double> maxDist(const CubicBezierCurve2d& _curve);

        /**
         * @brief Discretizes the given curve into a polyline. The discretization happens recursively and proceeds until the largest point of the curve is not further away than a given threshold from the straight line connecting the end points.
         * @param _curve Cubic Bezier curve to discretize.
         * @param _distanceThreshold Distance threshold that determines the discretization accuracy.
         * @param _polyline Resulting polyline that discretely approximates the given curve.
         */
        [[nodiscard]] static void topDown(const CubicBezierCurve2d& _curve, double _distanceThreshold, PiecewiseLinearCurve2d& _polyline);

        /**
         * @brief Discretizes the given spline into a polyline. The discretization happens recursively and proceeds until the largest point of the curve is not further away than a given threshold from the straight line connecting the end points.
         * @param _spline Cubic Bezier spline to discretize.
         * @param _distanceThreshold Distance threshold that determines the discretization accuracy.
         * @param _polyline Resulting polyline that discretely approximates the given curve.
         */
        [[nodiscard]] static void topDown(const CubicBezierSpline2d& _spline, double _distanceThreshold, PiecewiseLinearCurve2d& _polyline);

    private:
        /**
         * @brief Computes the largest distance from any point along the curve to the straight line that connects its start and end point.
         * @param _nodes Set of control points defining a cubic curve.
         * @param _implicit Implicit line to test against.
         * @return Largest distance and the curve parameter where the largest distance occured.
         */
        [[nodiscard]] static std::pair<double, double> maxDist_signed(const std::array<Eigen::Vector2d, 4>& _nodes, const ImplicitLine2d& _implicit);

        /**
         * @brief Computes the largest distance from any point along the curve to the straight line that connects its start and end point - in the positive half space only!
         * @param _nodes Set of control points defining a cubic curve.
         * @param _implicit Implicit line to test against.
         * @return Largest distance and the curve parameter where the largest distance occured.
         */
        [[nodiscard]] static std::pair<double, double> maxDist_positive(const std::array<Eigen::Vector2d, 4>& _nodes, const ImplicitLine2d& _implicit);

        /**
         * @brief Recursively subdivides a given curve into two pieces if the largest distance is above a given threshold.
         * @param _curve Cubic Bezier curve to discretize.
         * @param _distanceThreshold Distance threshold that determines the discretization accuracy.
         * @param _polyline Resulting polyline that discretely approximates the given curve.
         */
        static void topDownRecursive(const CubicBezierCurve2d& _curve, Eigen::Vector2d _range, double _distanceThreshold, PiecewiseLinearCurve2d& _polyline);
    };
}
