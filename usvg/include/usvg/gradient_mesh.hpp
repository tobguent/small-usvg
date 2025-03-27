#pragma once

#include "array.hpp"
#include "ferguson_patch.hpp"

namespace usvg
{
    /**
     * @brief Class that stores a gradient mesh. A gradient mesh is formed from a 2D array of Ferguson patches. The patches are in a regular layout and each vertex has positions/colors and tangent handles for both.
     */
    class GradientMesh
    {
    public:
        /*
           The indexing is not following the mathematical convention (row first, column second),
           but instead first indexes the u-dimension (columns), then the v-dimension (rows).
           Here, the resolution is (3,2).

           03 - - - 33  03 - - - 33  03 - - - 33
            | Tile3 |    | Tile4 |    | Tile5 |
            |  0,1  |    |  1,1  |    |  2,1  |
            |       |    |       |    |       |
           00 - - - 30  00 - - - 30  00 - - - 30
           03 - - - 33  03 - - - 33  03 - - - 33
            | Tile0 |    | Tile1 |    | Tile2 |
         v  |  0,0  |    |  1,0  |    |  2,0  |
            |       |    |       |    |       |
           00 - - - 30  00 - - - 30  00 - - - 30
                u
        */

        /**
         * @brief Coordinate in uv space.
         */
        using DomainCoord = Eigen::Vector2d;

        /**
         * @brief Representation of a single tile in the gradient mesh, containing the colors and positions represented via Ferguson patches.
         */
        class Tile
        {
        public:
            /**
             * @brief Type of color values stored at the control points.
             */
            using Value = Eigen::Vector3d;

            /**
             * @brief Computes the inverse of the coordinate partials.
             * @param _uv Domain location where to evaluate the surface.
             * @param _uv_x x-partial derivative of the uv coordinate.
             * @param _uv_y y-partial derivative of the uv coordinate.
             * @param _uv_xx xx-partial derivative of the uv coordinate.
             * @param _uv_yy yy-partial derivative of the uv coordinate.
             * @param _uv_xy xy-partial derivative of the uv coordinate.
             */
            void sample_InverseCoordinatePartials(const DomainCoord& _uv, DomainCoord& _uv_x, DomainCoord& _uv_y, DomainCoord& _uv_xx, DomainCoord& _uv_yy, DomainCoord& _uv_xy) const noexcept;

            /**
             * @brief Samples the Laplacian of color with respect to spatial coordinates at a given uv coordinate.
             * @param _uv Domain location where to evaluate the surface.
             * @return Laplacian of color with respect to xy coordinates.
             */
            [[nodiscard]] Value sample_ColorLaplacian_xy(const DomainCoord& _uv) const noexcept;

            /**
             * @brief Checks if the Ferguson patches of the gradient mesh are valid.
             * @return True if the Ferguson patches are valid.
             */
            [[nodiscard]] bool isValid() const;

            /**
             * @brief Positions of the patch, represented by a Ferguson patch.
             */
            FergusonPatch2d position;

            /**
             * @brief Color of the patch, represented by a Ferguson patch.
             */
            FergusonPatch3d color;
        };

        GradientMesh() noexcept;

        /**
         * @brief Gets a specific tile.
         * @param _i column index of the tile (u-direction).
         * @param _j row index of the tile (v-direction).
         * @return Selected tile.
         */
        [[nodiscard]] Tile& getTile(Eigen::Index _i, Eigen::Index _j);

        /**
         * @brief Gets a specific tile.
         * @param _i column index of the tile (u-direction).
         * @param _j row index of the tile (v-direction).
         * @return Selected tile.
         */
        [[nodiscard]] const Tile& getTile(Eigen::Index _i, Eigen::Index _j) const;

        /**
         * @brief Checks if the Ferguson patches of the gradient mesh are valid.
         * @return True if the Ferguson patches are valid.
         */
        [[nodiscard]] bool isValid() const;

        /**
         * @brief Gets the position of a specific vertex.
         * @param _i column index of the vertex (u-direction).
         * @param _j row index of the vertex (v-direction).
         * @return Position value.
         */
        const Eigen::Vector2d& getPosition(Eigen::Index _i, Eigen::Index _j) const;

        /**
         * @brief Sets the position of a specific vertex. The value is changed in all adjacent Ferguson patches.
         * @param _i column index of the vertex (u-direction).
         * @param _j row index of the vertex (v-direction).
         * @param _pos Position value to set.
         */
        void setPosition(Eigen::Index _i, Eigen::Index _j, const Eigen::Vector2d& _pos);

        /**
         * @brief Gets the u-tangent position of a specific vertex.
         * @param _i column index of the vertex (u-direction).
         * @param _j row index of the vertex (v-direction).
         * @return u-tangent of the position value.
         */
        const Eigen::Vector2d& getPosition_du(Eigen::Index _i, Eigen::Index _j) const;

        /**
         * @brief Sets the u-tangent position of a specific vertex. The value is changed in all adjacent Ferguson patches.
         * @param _i column index of the vertex (u-direction).
         * @param _j row index of the vertex (v-direction).
         * @param _tangentU u-tangent of the position value to set.
         */
        void setPosition_du(Eigen::Index _i, Eigen::Index _j, const Eigen::Vector2d& _tangentU);

        /**
         * @brief Gets the v-tangent position of a specific vertex.
         * @param _i column index of the vertex (u-direction).
         * @param _j row index of the vertex (v-direction).
         * @return v-tangent of the position value.
         */
        const Eigen::Vector2d& getPosition_dv(Eigen::Index _i, Eigen::Index _j) const;

        /**
         * @brief Sets the v-tangent position of a specific vertex. The value is changed in all adjacent Ferguson patches.
         * @param _i column index of the vertex (u-direction).
         * @param _j row index of the vertex (v-direction).
         * @param _tangentV v-tangent of the position value to set.
         */
        void setPosition_dv(Eigen::Index _i, Eigen::Index _j, const Eigen::Vector2d& _tangentV);

        /**
         * @brief Updates the Bezier representation of the position for all Ferguson patches that are adjacent to a specific vertex.
         * @param _i column index of the vertex (u-direction).
         * @param _j row index of the vertex (v-direction).
         */
        void updateBezierPosition(Eigen::Index _i, Eigen::Index _j);

        /**
         * @brief Updates the Bezier representation of the position for all Ferguson patches.
         */
        void updateBezierPosition();

        /**
         * @brief Gets the color of a specific vertex.
         * @param _i column index of the vertex (u-direction).
         * @param _j row index of the vertex (v-direction).
         * @return Color value.
         */
        const Eigen::Vector3d& getColor(Eigen::Index _i, Eigen::Index _j) const;

        /**
         * @brief Sets the color of a specific vertex. The value is changed in all adjacent Ferguson patches.
         * @param _i column index of the vertex (u-direction).
         * @param _j row index of the vertex (v-direction).
         * @param _color Color value to set.
         */
        void setColor(Eigen::Index _i, Eigen::Index _j, const Eigen::Vector3d& _color);

        /**
         * @brief Gets the u-tangent color of a specific vertex.
         * @param _i column index of the vertex (u-direction).
         * @param _j row index of the vertex (v-direction).
         * @return u-tangent of the color value.
         */
        const Eigen::Vector3d& getColor_du(Eigen::Index _i, Eigen::Index _j) const;

        /**
         * @brief Sets the u-tangent color of a specific vertex. The value is changed in all adjacent Ferguson patches.
         * @param _i column index of the vertex (u-direction).
         * @param _j row index of the vertex (v-direction).
         * @param _tangentU u-tangent of the color value to set.
         */
        void setColor_du(Eigen::Index _i, Eigen::Index _j, const Eigen::Vector3d& _tangentU);

        /**
         * @brief Gets the v-tangent color of a specific vertex.
         * @param _i column index of the vertex (u-direction).
         * @param _j row index of the vertex (v-direction).
         * @return v-tangent of the color value.
         */
        const Eigen::Vector3d& getColor_dv(Eigen::Index _i, Eigen::Index _j) const;

        /**
         * @brief Sets the v-tangent color of a specific vertex. The value is changed in all adjacent Ferguson patches.
         * @param _i column index of the vertex (u-direction).
         * @param _j row index of the vertex (v-direction).
         * @param _tangentV v-tangent of the color value to set.
         */
        void setColor_dv(Eigen::Index _i, Eigen::Index _j, const Eigen::Vector3d& _tangentV);

        /**
         * @brief Updates the Bezier representation of the color for all Ferguson patches that are adjacent to a specific vertex.
         * @param _i column index of the vertex (u-direction).
         * @param _j row index of the vertex (v-direction).
         */
        void updateBezierColor(Eigen::Index _i, Eigen::Index _j);

        /**
         * @brief Updates the Bezier representation of the color for all Ferguson patches.
         */
        void updateBezierColor();

        /**
         * @brief Given a position coordinate, this function first locates the uv coordinates in all tiles using Bezier clipping and then interpolates for each uv coordinate the colors. If the position is outside of the gradient mesh, an empty list is returned. Multiple colors occur if the gradient mesh has self-overlaps. Note that updateBezierPosition must have been called in advance.
         * @param _position Position at which to locate and interpolate colors.
         * @param _maxDepth Maximum recursion depth for the Bezier clipping.
         * @param _epsilon Numerical epsilon that determines the accuracy of the Bezier clipping.
         * @param _colors List of retrieved colors.
         */
        void sampleColors(const Eigen::Vector2d& _position, int _maxDepth, double _epsilon, std::vector<Eigen::Vector3d>& _colors) const;

        /**
         * @brief Given a position coordinate, this function first locates the uv coordinates in all tiles using Bezier clipping and then interpolates for each uv coordinate the color Laplacians. If the position is outside of the gradient mesh, an empty list is returned. Multiple colorLaplacians occur if the gradient mesh has self-overlaps. Note that updateBezierPosition must have been called in advance.
         * @param _position Position at which to locate and interpolate colors.
         * @param _maxDepth Maximum recursion depth for the Bezier clipping.
         * @param _epsilon Numerical epsilon that determines the accuracy of the Bezier clipping.
         * @param _colorLaplacians List of retrieved color Laplacians.
         */
        void sampleColorLaplacian(const Eigen::Vector2d& _position, int _maxDepth, double _epsilon, std::vector<Eigen::Vector3d>& _colorLaplacians) const;

        /**
         * @brief Number of tiles in column (u-axis) and row (v-axis) direction.
         */
        Eigen::Vector2i resolution;

        /**
         * @brief Tiles of the gradient mesh (row-wise).
         */
        std::vector<Tile> tiles;
    };
}
