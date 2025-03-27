#pragma once

#include "bezier_curve.hpp"
#include "bezier_spline.hpp"
#include "bezier_surface.hpp"

namespace usvg
{
    /**
     * @brief Implementation of the Bezier clipping algorithm of Sederberg and Nishita.
     * T.W. Sederberg, T. Nishita, "Curve intersection using Bézier clipping" Computer-Aided Design 22(9), 538-549, 1990.
     */
    class BezierClipping2d
    {
    public:
        /**
         * @brief Locates a two-dimensional point inside a bicubic Bezier surface and returns the uv coordinates. If the surface has self-overlaps, multiple results are reported. If an empty list is returned, then the point was not inside the Bezier surface.
         * @param _point Query point to find in the Bezier surface.
         * @param _surface Bicubic Bezier surface to search the point in.
         * @param _maxDepth Maximum recursion depth in the Bezier clipping search.
         * @param _epsilon Numerical epsilon that determines the precision of the uv coordinate.
         * @param _uvs List of found uv coordinates.
         * @return True if at least one uv coordinate was found.
         */
        [[nodiscard]] static bool locate(const Eigen::Vector2d& _point, const BicubicBezierSurface2d& _surface, int _maxDepth, double _epsilon, std::vector<Eigen::Vector2d>& _uvs);

    private:
        /**
         * @brief Returns the clipping range for the 16 distance values, given by a bicubic tensor product surface. The distance values are lifted, where each columns is assumed to have the same t-parameter value.
         * @param _nodesD Matrix of 16 distance values, one for each control point.
         * @return The clipping range [smin,smax] to discard parts of the surface that are certainly not cutting the zero distance line.
         */
        [[nodiscard]] static std::pair<double, double> clip(const Eigen::Matrix4d& _nodesD);

        /**
         * @brief Recursive Bezier clipping search for finding the uv coordinate in a bicubic Bezier surface that maps to a given position.
         * @param _point Query point to find in the Bezier surface.
         * @param _surface Bicubic Bezier surface to search the point in.
         * @param _rangeU Currently tested range of u-coordinates.
         * @param _rangeV Currently tested range of v-coordinates.
         * @param _depth Current recursion depth. Once this reaches zero, the recursion terminates.
         * @param _epsilon Numerical epsilon that determines the precision of the uv coordinate.
         * @param _initial_surface The initial bicubic Bezier surface, which is used to find the cutting directions, in case the numerical precision difficulties.
         * @param _uvs List of found uv coordinates.
         * @return True if at least one uv coordinate was found.
         */
        [[nodiscard]] static bool locateRecursive(const Eigen::Vector2d& _point, const BicubicBezierSurface2d& _surface, Eigen::Vector2d _rangeU, Eigen::Vector2d _rangeV, int _depth, double _epsilon, const BicubicBezierSurface2d& _initial_surface, std::vector<Eigen::Vector2d>& _uvs);
    };
}
