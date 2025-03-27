#pragma once

#include "bezier_spline.hpp"
#include "piecewise_linear_curve.hpp"

namespace usvg
{
    /**
     * @brief Class that represents a Poisson curve.
     *
     * Poisson Vector Graphics (PVG)
     * Fei Hou, Qian Sun, Zheng Fang, Yong-Jin Liu, Shi-Min Hu, Hong Qin, Aimin Hao, Ying He
     * IEEE Transactions on Visualization and Computer Graphics 26 (2), 1361-1371, 2020, doi:10.1109/TVCG.2018.2867478
     */
    class PoissonCurve
    {
    public:
        /**
         * @brief Checks if the curve is properly initialized.
         * @return True if valid.
         */
        [[nodiscard]] bool isValid() const;

        /**
         * @brief Position curve.
         */
        CubicBezierSpline2d position;

        /**
         * @brief Laplacian along the curve.
         */
        PiecewiseLinearCurve3d weights;
    };
}
