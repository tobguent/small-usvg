#pragma once

#include "bezier_curve.hpp"
#include "bezier_spline.hpp"
#include "boundary_condition.hpp"
#include "edge_graph.hpp"
#include "piecewise_linear_curve.hpp"

namespace usvg
{
    class GradientMesh;
    class PoissonCurve;

    /**
     * @brief Represents a patch, consisting of curves for the boundary conditions and source terms for the interior, originating from gradient meshes and Poisson curves.
     */
    class Patch
    {
    public:
        /**
         * @brief Patch boundary curve segment.
         */
        class BoundaryCurve
        {
        public:
            /**
             * @brief Constructor.
             */
            BoundaryCurve() = default;

            /**
             * @brief Constructs a patch boundary curve.
             * @param _position Position of the curve.
             * @param _color Color of the curve on the right side.
             * @param _boundaryCondition Boundary condition of the curve on the right side.
             * @param _gradientMesh Optional parent gradient mesh.
             */
            BoundaryCurve(const PiecewiseLinearCurve2d& _position, const CubicBezierSpline3d& _color, const EBoundaryCondition _boundaryCondition, std::shared_ptr<GradientMesh> _gradientMesh);

            /**
             * @brief Checks if the curve is properly initialized.
             * @return True if valid.
             */
            [[nodiscard]] bool isValid() const;

            /**
             * @brief Position of the curve. The orientation aligns with the direction of the boundary curve.
             */
            PiecewiseLinearCurve2d position;

            /**
             * @brief Color of the curve on the right side. The orientation aligns with the direction of the boundary curve.
             */
            CubicBezierSpline3d color;

            /**
             * @brief Boundary condition of the curve on the right side.
             */
            EBoundaryCondition boundaryCondition;

            /**
             * @brief Optional reference to a parent gradient mesh in case this curve bounds a gradient mesh.
             */
            std::shared_ptr<GradientMesh> gradientMesh;
        };

        /**
         * @brief Discrete representation of a Poisson curve.
         */
        class PoissonCurveDiscrete
        {
        public:
            /**
             * @brief Constructor.
             */
            PoissonCurveDiscrete() = default;

            /**
             * @brief Constructs a patch boundary curve.
             * @param _position Position of the curve.
             * @param _weights Laplacian on both sides of the curve.
             */
            PoissonCurveDiscrete(const PiecewiseLinearCurve2d& _position, const PiecewiseLinearCurve3d& _weights);

            /**
             * @brief Checks if the curve is properly initialized.
             * @return True if valid.
             */
            [[nodiscard]] bool isValid() const;

            /**
             * @brief Positions
             */
            PiecewiseLinearCurve2d position;

            /**
             * @brief Laplacian on both sides of the curve.
             */
            PiecewiseLinearCurve3d weights;
        };

        /**
         * @brief Container that stores a loop of boundary curves.
         */
        class Loop
        {
        public:
            /**
             * @brief Constructor.
             */
            Loop();

            /**
             * @brief Constructor that computes the bounding box on its own.
             * @param _boundaryCurves Closed sequence of boundary curves.
             * @param _turningNumber Turning number of this closed loop.
             */
            Loop(const std::vector<BoundaryCurve>& _boundaryCurves, double _turningNumber);

            /**
             * @brief Checks if the curve is properly initialized.
             * @return True if valid.
             */
            [[nodiscard]] bool isValid() const;

            /**
             * @brief Tests if this boundary curve is on the interior of a closed loop with non-empty interior.
             * @return True if this boundary curve is on the interior of a closed loop with non-empty interior.
             */
            [[nodiscard]] bool isInterior() const noexcept;

            /**
             * @brief Tests if this boundary curve is on the exterior of a closed loop with non-empty interior.
             * @return True if this boundary curve is on the exterior of a closed loop with non-empty interior.
             */
            [[nodiscard]] bool isExterior() const noexcept;

            /**
             * @brief Tests if this boundary curve has an empty interior.
             * @return True if this boundary curve has an empty interior.
             */
            [[nodiscard]] bool isFlat() const noexcept;

            /**
             * @brief Tests if a point is contained in the loop using the winding number theorem. For the computation, the discrete representation is used!
             * @param _point Point to test.
             * @param _epsilon Allowed numerical epsilon on the winding number.
             * @return True if the point is contained.
             */
            [[nodiscard]] bool contains(const Eigen::Vector2d& _point, double _epsilon) const;

            /**
             * @brief Tests if an entire loop is contained in this loop using the winding number theorem.
             * @param _other Loop to test for containment.
             * @param _epsilon Allowed numerical epsilon on the winding number.
             * @return True if the other loop is fully contained in this loop.
             */
            [[nodiscard]] bool contains(const Loop& _other, double _epsilon) const;

            /**
             * @brief Gets the bounding box of the boundary curves. Note that recomputeBoundingBox has to be called first!
             * @return Bounding box of this geometry.
             */
            [[nodiscard]] const Eigen::AlignedBox2d& getBoundingBox() const noexcept;

            /**
             * @brief Recomputes the bounding box of the boundary curves.
             */
            void recomputeBoundingBox() noexcept;

            /**
             * @brief Turning number of this closed loop.
             */
            double turningNumber;

            /**
             * @brief Closed sequence of boundary curves.
             */
            std::vector<BoundaryCurve> boundaryCurves;

            /**
             * @brief Calculates the winding number for a point and a given loop.
             * @param _point Point to compute the winding number for.
             * @return Signed winding number.
             */
            [[nodiscard]] double windingNumber(const Eigen::Vector2d& _point) const;

            /**
             * @brief Computes the signed turning number between the two vectors (point3-point2) and (point2-point1)
             * @param _point1 First point of the triplet.
             * @param _point2 Second point of the triplet.
             * @param _point3 Third point of the triplet.
             * @return Signed integer-valued turning number.
             */
            [[nodiscard]] static double computeTurningNumber(const Eigen::Vector2d& _point1, const Eigen::Vector2d& _point2, const Eigen::Vector2d& _point3) noexcept;

            /**
             * @brief Computes the turning number for a sequence of edges in the discretized patch boundary curve.
             * @param _boundaryCurves Closed loop of patch boundary curves. Note that subsequent curves in reverse direction are not allowed, since the turning angle is not unique!
             * @return If the line is closed, we get a signed integer number (+/- numerical noise).
             */
            [[nodiscard]] static double computeTurningNumber(const std::vector<BoundaryCurve>& _boundaryCurves);

        private:
            /**
             * @brief Bounding box of this loop.
             */
            Eigen::AlignedBox2d mBoundingBox;
        };

        /**
         * @brief Checks if the curve is properly initialized.
         * @return True if valid.
         */
        [[nodiscard]] bool isValid() const;

        /**
         * @brief Tests if a point is contained inside the patch using the winding number theorem. This test uses the discrete representation.
         * @param _point Point to test for containment.
         * @param _epsilon Allowed numerical epsilon on the winding number.
         * @return True if the point is contained in the patch.
         */
        [[nodiscard]] bool contains(const Eigen::Vector2d& _point, double _epsilon) const;

        /**
         * @brief Creates the left patch boundary curve for a given edge.
         * @param _edge Edge to create boundary curve for.
         * @return Left patch boundary curve.
         */
        static Patch::BoundaryCurve makeLeftPatchBoundaryCurve(const EdgeGraph::Edge& _edge);

        /**
         * @brief Creates the right patch boundary curve for a given edge.
         * @param _edge Edge to create boundary curve for.
         * @return Right patch boundary curve.
         */
        static Patch::BoundaryCurve makeRightPatchBoundaryCurve(const EdgeGraph::Edge& _edge);

        /**
         * @brief Traverses an edge graph to find all closed loops.
         * @param _graph Edge graph.
         * @return Vector of closed loops.
         */
        static std::vector<std::shared_ptr<Loop>> computeLoops(const EdgeGraph& _graph);

        /**
         * @brief Assembles a collection of loops into patches.
         * @param _loops Loops to form the boundaries of patches from.
         * @param _poissonCurves Set of poisson curves to insert into the patches.
         * @param _discretizationResidual Discretization residual when discreting the curves.
         * @return Set of patches.
         */
        static std::vector<std::shared_ptr<Patch>> computePatches(const std::vector<std::shared_ptr<Loop>>& _loops, const std::vector<std::shared_ptr<PoissonCurve>>& _poissonCurves, double _discretizationResidual);

        /**
         * @brief Vector of loops surrounding the patch.
         */
        std::vector<std::shared_ptr<Loop>> loops;

        /**
         * @brief Optional gradient meshes that determine the fill color. If none is set, the source term is zero.
         */
        std::vector<std::shared_ptr<GradientMesh>> gradientMeshes;

        /**
         * @brief Poisson curves that exist inside the domain.
         */
        std::vector<std::shared_ptr<PoissonCurveDiscrete>> poissonCurves;
    };
}
