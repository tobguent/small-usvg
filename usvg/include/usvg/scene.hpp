#pragma once

#include <usvg/array.hpp>
#include <usvg/jacobi_iteration.hpp>
#include <usvg/piecewise_linear_curve.hpp>

#include "tinyxml2.h"
#include <Eigen/Eigen>

#include <fstream>
#include <iostream>

namespace usvg
{
    class DiffusionCurve;
    class EdgeGraph;
    class GradientMesh;
    class Image;
    class Patch;
    class PoissonCurve;

    /**
     * @brief Unified smooth vector graphics scene.
     *
     * Unified Smooth Vector Graphics: Modeling Gradient Meshes and Curve-based Approaches Jointly as Poisson Problem
     * Xingze Tian, Tobias G�nther
     * ArXiv, 17 Aug 2024, doi:10.48550/arXiv.2408.09211
     */
    class Scene
    {
    public:
        /**
         * @brief In case multiple gradient meshes overlap, the Laplacians are blended. This enum says which blend operation to apply.
         */
        enum class EBlendOperation
        {
            /**
             * @brief If multiple gradient meshes overlap, none of them is used. The Laplacian is set to zero.
             */
            Zero,
            /**
             * @brief The Laplacians of the overlapping gradient meshes are added up.
             */
            Sum,
            /**
             * @brief The Laplacians of the overlapping gradient meshes are averages.
             */
            Average,
            /**
             * @brief The Laplacian of the first gradient mesh is used.
             */
            First
        };

        /**
         * @brief Represents a 2D point (x, y).
         */
        using point_type = Eigen::Vector2d;

        /**
         * @brief Represents a 3D color (r, g, b).
         */
        using color_type = Eigen::Vector3d;

        /**
         * @brief Represents a 3D color at a certain parameter location (r, g, b, t).
         */
        using color_point_type = Eigen::Vector4d;

        /**
         * @brief Constructor.
         */
        Scene();

        /**
         * @brief Reads diffusion curves from a given XML element.
         * @param _parent_element Parent element to read from.
         * @param _swap Swaps the x,y coordinates and the red and blue color channel. This is for compatibility with Orzan's file format.
         * @param _scene Scene to add the diffusion curves to.
         * @return True if successful.
         */
        [[nodiscard]] static bool readDiffusionCurves(const tinyxml2::XMLElement* _parent_element, bool _swap, Scene* _scene);

        /**
         * @brief Reads Poisson curves from a given XML element.
         * @param _parent_element Parent element to read from.
         * @param _scene Scene to add the Poisson curves to.
         * @return True if successful.
         */
        [[nodiscard]] static bool readPoissonCurves(const tinyxml2::XMLElement* _parent_element, Scene* _scene);

        /**
         * @brief Reads gradient meshes from a given XML element.
         * @param _parent_element Parent element to read from.
         * @param _scene Scene to add the gradient mesh to.
         * @return True if successful.
         */
        [[nodiscard]] static bool readGradientMeshes(const tinyxml2::XMLElement* _parent_element, Scene* _scene);

        /**
         * @brief Reads a scene from an XML file.
         * @param path Path to the file to read.
         * @return True if reading was successful.
         */
        bool readXML(const std::string& _path);

        /**
         * @brief Read 2D positions (x,y) from the XML file.
         * @param _points Output array that receives the points.
         * @param _parent Parent XML element to read from.
         * @param _child_name Name of the child to read from.
         * @param _num_points Number of points to read.
         * @param _is_normalized Flag that determines whether the data comes in normalized.
         * @param _domain Domain over which the image is defined.
         * @param _swap Flag that enables swapping of the x,y coordinates.
         * @return True if successful.
         */
        [[nodiscard]] static bool readPositions(Array2d& _points, const tinyxml2::XMLElement* _parent, const std::string& _child_name, int _num_points, bool _is_normalized, const Eigen::AlignedBox2d& _domain, bool _swap);

        /**
         * @brief Read 3D colors (R,G,B) from the XML file
         * @param _colors Output vector that receives the colors.
         * @param _parent Parent XML element to read from.
         * @param _child_name Name of the child to read from.
         * @param _num_colors Number of colors to read.
         * @param _swap Flag that enables swapping of the R and B channel.
         * @return True if successful.
         */
        [[nodiscard]] static bool readColors(Array3d& _colors, const tinyxml2::XMLElement* _parent, const std::string& _child_name, int _num_colors, bool _swap);

        /**
         * @brief Reads piecewise linear tri-variate curve, storing for example 3D colors (r,g,b).
         * @param _polyline Output polyline into which the values are written.
         * @param _parent Parent XML element to read from.
         * @param _child_name Name of the child to read from.
         * @param _num_points Number of points to read.
         * @param _swap Flag that enables swapping of the R and B channel.
         * @return True if successful.
         */
        [[nodiscard]] static bool readPiecewiseLinearCurve3d(PiecewiseLinearCurve3d& _polyline, const tinyxml2::XMLElement* _parent, const std::string& _child_name, int _num_points, bool _swap);

        /**
         * @brief Checks if the scene is properly initialized.
         * @return True if valid.
         */
        [[nodiscard]] bool isValid() const;

        /**
         * @brief Performs the patch construction.
         * @param vertexMergeThreshold Edges are merged if their end points have a distance below this threshold.
         * @param clippingEpsilon Residual epsilon during Bezier clipping.
         * @param discretizationResidual Discretization residual when discreting the curves.
         */
        void buildPatches(double _vertexMergeThreshold, double _clippingEpsilon, double _discretizationResidual);

        /**
         * @brief Numerically estimates the Laplacian of the gradient mesh.
         * @param _mesh Gradient mesh to compute the Laplacian for.
         * @param _coord Location where to compute the Laplacian.
         * @param _maxDepth Maximum recursion depth for the Bezier clipping.
         * @param _epsilon Numerical epsilon that determines the accuracy of the Bezier clipping.
         * @param _spacing Grid spacing between pixels.
         * @param _open_x0 Flag that determines the left neighbor is not separated through a curve.
         * @param _open_x1 Flag that determines the right neighbor is not separated through a curve.
         * @param _open_y0 Flag that determines the bottom neighbor is not separated through a curve.
         * @param _open_y1 Flag that determines the top neighbor is not separated through a curve.
         * @return
         */
        static std::tuple<bool, Eigen::Vector3d, Eigen::Vector3d> sampleGradientMeshLaplacian(
            const GradientMesh& _mesh,
            Eigen::Vector2d& _coord,
            int _maxDepth,
            double _epsilon,
            const Eigen::Vector2d& _spacing,
            bool _open_x0, bool _open_x1, bool _open_y0, bool _open_y1);

        /**
         * @brief computes the intersection of two line segments. t contains the relative positions on the two segments.
         * @param _a Start of the first segment.
         * @param _b End of the first segment.
         * @param _c Start of the second segment.
         * @param _d End of the second segment.
         * @param _t Relative parameter between ab and cd where the intersection occurs.
         * @return True if there is an intersection.
         */
        static bool intersect_segment_segment(const Eigen::Vector2d& _a, const Eigen::Vector2d& _b, const Eigen::Vector2d& _c, const Eigen::Vector2d& _d, Eigen::Vector2d& _t);

        /**
         * @brief Returns the closest point to p on a segment with endpoints b0 and b1
         * @param _b0 Start of the segment.
         * @param _b1 End of the segment.
         * @param _p Query point to find closest point on the segment for.
         * @param _min_dist Output variable that stores the distance.
         * @param _min_t Output variable that stores the relative location on the segment where the closest point was found.
         */
        static void closest_point_on_segment(const Eigen::Vector2d& _b0, const Eigen::Vector2d& _b1, const Eigen::Vector2d& _p, double& _min_dist, double& _min_t);

        /**
         * @brief Initializes the Jacobi solver from a given scene.
         * @param scene Scene that contains the patches.
         * @param jacobi Jacobi solver to initialize.
         */
        static void initializeFromPatches(const Scene& scene, JacobiIteration3d& jacobi);

        /**
         * @brief Computes an image by solving the scene PDE with Jacobi relaxation.
         * @param resolution Resolution of the output image.
         * @param numIterations Number of Jacobi iterations.
         * @param useMultigrid Flag that enables the multi-grid solver.
         * @return Image containing the resulting raster image.
         */
        std::shared_ptr<Image> solvePDE(const Eigen::Vector2i& resolution, int numIterations, bool useMultigrid) const;

        /**
         * @brief Vector of input diffusion curves.
         */
        std::vector<std::shared_ptr<DiffusionCurve>> diffusionCurves;

        /**
         * @brief Vector of input Poisson curves.
         */
        std::vector<std::shared_ptr<PoissonCurve>> poissonCurves;

        /**
         * @brief Vector of input gradient meshes.
         */
        std::vector<std::shared_ptr<GradientMesh>> gradientMeshes;

        /**
         * @brief Domain over which the image is defined.
         */
        Eigen::AlignedBox2d domain;

        /**
         * @brief Blending operation to apply when multiple gradient meshes overlap.
         */
        EBlendOperation blendOperation;

        /**
         * @brief Patches that were build from the input diffusion curves, gradient meshes, and poisson curves.
         */
        std::vector<std::shared_ptr<Patch>> patches;
    };
}
