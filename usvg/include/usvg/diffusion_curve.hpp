#pragma once

#include "bezier_spline.hpp"
#include "boundary_condition.hpp"
#include "piecewise_linear_curve.hpp"

namespace usvg
{
    /**
     * @brief Class that represents a single diffusion curve.
     *
     * The concept or diffusion curves goes back to Orzan et al.
     *   Diffusion curves: a vector representation for smooth-shaded images.
     *   Alexandrina Orzan, Adrien Bousseau, Pascal Barla, Holger Winnemöller, Joëlle Thollot, and David Salesin.
     *   Communications of the ACM 56 (7), 101–108, 2013. https://doi.org/10.1145/2483852.2483873
     * Follow-up papers discarded the smoothing post-processes and modeled diffusion curves as Dirichlet boundary conditions. We use this formulation.
     *   A GPU Laplacian solver for diffusion curves and Poisson image editing.
     *   Stefan Jeschke, David Cline, and Peter Wonka.
     *   In ACM SIGGRAPH Asia 2009 papers. ACM, New York, NY, USA, Article 116, 1–8, 2009. https://doi.org/10.1145/1661412.1618462
     */
    class DiffusionCurve
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
         * @brief Colors on the left side of the curve.
         */
        PiecewiseLinearCurve3d colorLeft;

        /**
         * @brief Colors on the right side of the curve.
         */
        PiecewiseLinearCurve3d colorRight;

        /**
         * @brief Type of boundary condition on the left side of the curve.
         */
        EBoundaryCondition boundaryConditionLeft;

        /**
         * @brief Type of boundary condition on the right side of the curve.
         */
        EBoundaryCondition boundaryConditionRight;
    };
}
