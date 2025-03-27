#pragma once

#include "bezier_spline.hpp"
#include "boundary_condition.hpp"
#include "bvh.hpp"
#include "piecewise_linear_curve.hpp"

#include <stack>

namespace usvg
{
    class DiffusionCurve;
    class GradientMesh;

    /**
     * @brief Data structure of an undirected graph that is used in the construction of patches.
     */
    class EdgeGraph
    {
    public:
        /**
         * @brief Representation of an input boundary curve.
         */
        class InputBoundaryCurve
        {
        public:
            /**
             * @brief Intersects the input boundary curve with another one and returns all intersections.
             * @param _other Other input boundary curve to intersect with.
             * @param _intersections List of intersections. The first component is the t-value of this curve, and the second component is the t-value of the other curve.
             * @return True if there was at least one inetersection.
             */
            [[nodiscard]] bool intersect(const InputBoundaryCurve& _other, std::vector<Eigen::Vector2d>& _intersections) const;

            /**
             * @brief Subdivides this curve into two pieces.
             * @param _t Parameter location where to split.
             * @return Pair of the two curve pieces.
             */
            [[nodiscard]] std::pair<std::shared_ptr<InputBoundaryCurve>, std::shared_ptr<InputBoundaryCurve>> subdivide(double _t) const;

            /**
             * @brief Finds the closest point to a given point.
             * @param _position Point to find the closest point for.
             * @return Position of the closest point, distance to the curve, curve parameter.
             */
            std::tuple<Eigen::Vector2d, double, double> closestPoint(const Eigen::Vector2d& _point);

            /**
             * @brief Rebuilds the bounding volume hierarchy.
             */
            void update();

            /**
             * @brief Position curve.
             */
            PiecewiseLinearCurve2d position;

            /**
             * @brief Colors on the left side of the curve.
             */
            CubicBezierSpline3d colorLeft;

            /**
             * @brief Colors on the right side of the curve.
             */
            CubicBezierSpline3d colorRight;

            /**
             * @brief Type of boundary condition on the left side of the curve.
             */
            EBoundaryCondition boundaryConditionLeft;

            /**
             * @brief Type of boundary condition on the right side of the curve.
             */
            EBoundaryCondition boundaryConditionRight;

            /**
             * @brief Optional parent gradient mesh in case this was constructed from a gradient mesh.
             */
            std::shared_ptr<GradientMesh> gradientMesh;

            /**
             * @brief Bounding volume hierarchy of the position.
             */
            Bvh2d bvh;
        };

        class Vertex;
        class Edge;

        /**
         * @brief Vertex of the edge graph data structure.
         */
        class Vertex
        {
        public:
            /**
             * @brief Constructor.
             * @param _position Position of the vertex.
             */
            Vertex(const Eigen::Vector2d& _position);

            /**
             * @brief Gets the next edge, turning right.
             * @param _incoming Edge from which we are coming.
             * @return Next edge, turning right from the incoming direction.
             */
            std::shared_ptr<Edge> nextEdge(std::shared_ptr<const Edge> _incoming);

            /**
             * @brief Position of the vertex.
             */
            Eigen::Vector2d position;

            /**
             * @brief List of edges connected to this vertex.
             */
            std::list<std::weak_ptr<Edge>> edges;
        };

        /**
         * @brief Edge of the edge graph data structure.
         */
        class Edge : public std::enable_shared_from_this<Edge>
        {
        public:
            /**
             * @brief Constructor.
             * @param _boundaryCurve Input boundary curve underneath the edge.
             */
            Edge(std::shared_ptr<InputBoundaryCurve> _boundaryCurve);

            /**
             * @brief Intersects the edge with another edge and reports all intersections.
             * @param _other Other edge to intersect with.
             * @param _intersections List of intersections, containing the parameter values for both curves. The first component is the t-value of this edge, and the second component is the t-value of the other edge.
             * @return True if there is an intersection.
             */
            [[nodiscard]] bool intersect(const Edge& _other, std::vector<Eigen::Vector2d>& _intersections) const;

            /**
             * @brief Subdivides the edge into two pieces. This edge becomes the first piece and the second piece as well as the split vertex are returned.
             * @param _t Parameter location where to split.
             * @return The second edge segment and the split vertex.
             */
            [[nodiscard]] std::pair<std::shared_ptr<Edge>, std::shared_ptr<Vertex>> subdivide(double _t);

            /**
             * @brief Vertex at which the edge begins.
             */
            std::weak_ptr<Vertex> vertexA;

            /**
             * @brief Vertex at which the edge ends.
             */
            std::weak_ptr<Vertex> vertexB;

            /**
             * @brief Boundary curve that is underneath the edge.
             */
            std::shared_ptr<InputBoundaryCurve> boundaryCurve;
        };

        /**
         * @brief Constructor.
         * @param _diffusionCurves Vector of diffusion curves.
         * @param _gradientMeshes Vector of gradient meshes.
         * @param _vertexMergeThreshold Edges are merged if their end points have a distance below this threshold.
         * @param _epsilon Residual epsilon during Bezier clipping. Should be smaller than vertexMergeThreshold. Otherwise intersection points might not snap.
         * @param _distanceThreshold Threshold for discretization error.
         */
        EdgeGraph(const std::vector<std::shared_ptr<DiffusionCurve>>& _diffusionCurves,
                  const std::vector<std::shared_ptr<GradientMesh>>& _gradientMeshes,
                  double _vertexMergeThreshold, double _epsilon, double _distanceThreshold);

        /**
         * @brief Vector of vertices in the graph.
         */
        std::vector<std::shared_ptr<Vertex>> vertices;

        /**
         * @brief Vector of edges in the graph.
         */
        std::vector<std::shared_ptr<Edge>> edges;

        /**
         * @brief Vector of input boundary curves that the edge graph is constructed from.
         */
        std::vector<std::shared_ptr<InputBoundaryCurve>> inputBoundaryCurves;

    private:
        /**
         * @brief Constructs an input boundary curve from a diffusion curve.
         * @param _curve Diffusion curve to convert.
         * @param _distanceThreshold Threshold for discretization error.
         * @return Input boundary curve.
         */
        [[nodiscard]] static std::shared_ptr<InputBoundaryCurve> makeInputBoundaryCurveFromDiffusionCurve(std::shared_ptr<const DiffusionCurve> _curve, double _distanceThreshold);

        /**
         * @brief Constructs four input boundary curves from a gradient mesh.
         * @param _gradientMesh Gradient mesh to convert.
         * @param _distanceThreshold Threshold for discretization error.
         * @return Four input boundary curves, ordered North, East, South, West.
         */
        [[nodiscard]] static std::tuple<std::shared_ptr<InputBoundaryCurve>, std::shared_ptr<InputBoundaryCurve>, std::shared_ptr<InputBoundaryCurve>, std::shared_ptr<InputBoundaryCurve>> makeInputBoundaryCurveFromGradientMesh(std::shared_ptr<GradientMesh> _gradientMesh, double _distanceThreshold);

        /**
         * @brief Assembles the set of input boundary curves from the given diffusion curves and gradient meshes.
         * @param _diffusionCurves Vector of diffusion curves.
         * @param _gradientMeshes Vector of gradient meshes.
         * @param _distanceThreshold Threshold for discretization error.
         */
        void assembleInputBoundaryCurves(const std::vector<std::shared_ptr<DiffusionCurve>>& _diffusionCurves, const std::vector<std::shared_ptr<GradientMesh>>& _gradientMeshes, double _distanceThreshold);

        /**
         * @brief Builds the edge graph incrementally from the input boundary curves.
         * @param _vertexMergeThreshold Edges are merged if their end points have a distance below this threshold.
         * @param _epsilon Residual epsilon during Bezier clipping.
         */
        void constructEdgeGraph(double _vertexMergeThreshold, double _epsilon);

        /**
         * @brief Deletes all invalid vertices and edges from the given vectors.
         * @param _vertices Vertex array to delete invalid vertices from.
         * @param _edges Edge array to delete invalid edges from.
         * @param _epsilon Minimal distance for edge lengths.
         */
        static void deleteInvalid(std::vector<std::shared_ptr<EdgeGraph::Vertex>>& _vertices, std::vector<std::shared_ptr<EdgeGraph::Edge>>& _edges, double _epsilon);

        /**
         * @brief Redirects all the incoming edge of "from" to vertex "to".
         * @param _from Vertex to be disconnected from all its edges.
         * @param _to Vertex which receives the connections.
         */
        static void redirectEdges(std::shared_ptr<EdgeGraph::Vertex> _from, std::shared_ptr<EdgeGraph::Vertex> _to);
    };
}
