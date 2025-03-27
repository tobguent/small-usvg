#pragma once

#include "regular_grid.hpp"

namespace usvg
{
    /**
     * @brief Class that performs a Jacobi iteration in two dimensions.
     * @tparam TValue Value that is stored at each grid point. Must be an Eigen::Vector type.
     */
    template <typename TValue>
    class JacobiIteration
    {
    public:
        /**
         * @brief Type of values being stored at each grid point.
         */
        using Value = TValue;

        /**
         * @brief Types of pixels.
         */
        enum EMask : uint8_t
        {
            /**
             * @brief Dirichlet boundary condition is set at this grid point.
             */
            Dirichlet,
            /**
             * @brief Neumann boundary condition is set at this grid point.
             */
            Neumann,
            /**
             * @brief The value at this grid point is unknown.
             */
            Unknown
        };

        /**
         * @brief Constructor that allocates the interal fields for the given grid.
         * @param _grid Grid that provides resolution and spacing.
         */
        JacobiIteration(std::shared_ptr<RegularGrid2d> _grid)
            : ping(0)
            , grid(_grid)
        {
            Eigen::Vector2i resolution = _grid->getResolution();
            openEast.resize((resolution.x() - 1) * resolution.y(), true);
            openNorth.resize(resolution.x() * (resolution.y() - 1), true);
            source.resize(resolution.prod(), TValue::Zero());
            mask.resize(resolution.prod(), EMask::Unknown);
            field[0].resize(resolution.prod(), TValue::Zero());
            field[1].resize(resolution.prod(), TValue::Zero());
        }

        /**
         * @brief Resets the output field to zero.
         */
        void reset()
        {
            std::memset(field[ping].data(), 0, field[ping].size() * sizeof(TValue));
            std::memset(field[1 - ping].data(), 0, field[1 - ping].size() * sizeof(TValue));
        }

        /**
         * @brief Performs a certain number of Jacobi iterations
         * @param numIterations Number of Jacobi iterations to perform.
         */
        void iterate(int numIterations)
        {
            Eigen::Vector2i resolution = grid->getResolution();
            Eigen::Index numPixels     = resolution.prod();

            // iterate the Jacobi relaxations
            for (Eigen::Index it = 0; it < numIterations; ++it)
            {
#ifdef NDEBUG
#pragma omp parallel for
#endif
                for (Eigen::Index linearIndex = 0; linearIndex < numPixels; ++linearIndex)
                {
                    // compute pixel index
                    Eigen::Vector2i gridCoord   = grid->getGridCoord(linearIndex);
                    Eigen::Index linearIndex_x0 = linearIndex - 1;
                    Eigen::Index linearIndex_x1 = linearIndex + 1;
                    Eigen::Index linearIndex_y0 = linearIndex - resolution.x();
                    Eigen::Index linearIndex_y1 = linearIndex + resolution.x();

                    // get source and mask
                    TValue s = source[linearIndex];
                    EMask m  = mask[linearIndex];

                    // if Dirichlet condition
                    if (m == EMask::Dirichlet)
                    {
                        field[1 - ping][linearIndex] = s;
                    }
                    else
                    {
                        bool open_x0 = getOpenWest(gridCoord);
                        bool open_x1 = getOpenEast(gridCoord);
                        bool open_y0 = getOpenSouth(gridCoord);
                        bool open_y1 = getOpenNorth(gridCoord);

                        // Poisson update
                        int numOpen = (open_x0 ? 1 : 0) + (open_x1 ? 1 : 0) + (open_y0 ? 1 : 0) + (open_y1 ? 1 : 0);
                        if (numOpen > 0)
                        {
                            TValue clr_x0 = open_x0 ? ((mask[linearIndex_x0] == EMask::Dirichlet) ? source[linearIndex_x0] : field[ping][linearIndex_x0]) : TValue::Zero();
                            TValue clr_x1 = open_x1 ? ((mask[linearIndex_x1] == EMask::Dirichlet) ? source[linearIndex_x1] : field[ping][linearIndex_x1]) : TValue::Zero();
                            TValue clr_y0 = open_y0 ? ((mask[linearIndex_y0] == EMask::Dirichlet) ? source[linearIndex_y0] : field[ping][linearIndex_y0]) : TValue::Zero();
                            TValue clr_y1 = open_y1 ? ((mask[linearIndex_y1] == EMask::Dirichlet) ? source[linearIndex_y1] : field[ping][linearIndex_y1]) : TValue::Zero();

                            Eigen::Vector3d new_clr      = (clr_x0 + clr_x1 + clr_y0 + clr_y1 - s * grid->getSpacing().prod()) / numOpen;
                            field[1 - ping][linearIndex] = new_clr;
                        }
                    }
                }
                ping = 1 - ping;
            }
        }

        /**
         * @brief Sets the is open flag for the east neighbor. If true, information may diffuse this way.
         * @param _gridCoord Grid coord at which the openness is viewed.
         * @param _open True if information can diffuse this way.
         */
        void setOpenEast(const Eigen::Vector2i& _gridCoord, bool _open)
        {
            if (0 <= _gridCoord.x() && _gridCoord.x() < grid->getResolution().x() - 1 &&
                0 <= _gridCoord.y() && _gridCoord.y() < grid->getResolution().y())
            {
                openEast[_gridCoord.y() * (grid->getResolution().x() - 1) + _gridCoord.x()] = _open;
            }
        }

        /**
         * @brief Gets the is open flag for the east neighbor. If true, information may diffuse this way.
         * @param _gridCoord Grid coord at which the openness is viewed.
         * @return True if information can diffuse this way.
         */
        [[nodiscard]] bool getOpenEast(const Eigen::Vector2i& _gridCoord) const
        {
            if (0 <= _gridCoord.x() && _gridCoord.x() < grid->getResolution().x() - 1 &&
                0 <= _gridCoord.y() && _gridCoord.y() < grid->getResolution().y())
            {
                return openEast[_gridCoord.y() * (grid->getResolution().x() - 1) + _gridCoord.x()];
            }
            else
                return false;
        }

        /**
         * @brief Gets the is open flag for the west neighbor. If true, information may diffuse this way.
         * @param _gridCoord Grid coord at which the openness is viewed.
         * @return True if information can diffuse this way.
         */
        [[nodiscard]] bool getOpenWest(const Eigen::Vector2i& _gridCoord) const
        {
            if (1 <= _gridCoord.x() && _gridCoord.x() < grid->getResolution().x() &&
                0 <= _gridCoord.y() && _gridCoord.y() < grid->getResolution().y())
            {
                return openEast[_gridCoord.y() * (grid->getResolution().x() - 1) + _gridCoord.x() - 1];
            }
            else
                return false;
        }

        /**
         * @brief Sets the is open flag for the north neighbor. If true, information may diffuse this way.
         * @param _gridCoord Grid coord at which the openness is viewed.
         * @param _open True if information can diffuse this way.
         */
        void setOpenNorth(const Eigen::Vector2i& _gridCoord, bool _open)
        {
            if (0 <= _gridCoord.x() && _gridCoord.x() < grid->getResolution().x() &&
                0 <= _gridCoord.y() && _gridCoord.y() < grid->getResolution().y() - 1)
            {
                openNorth[_gridCoord.y() * grid->getResolution().x() + _gridCoord.x()] = _open;
            }
        }

        /**
         * @brief Gets the is open flag for the north neighbor. If true, information may diffuse this way.
         * @param _gridCoord Grid coord at which the openness is viewed.
         * @return True if information can diffuse this way.
         */
        [[nodiscard]] bool getOpenNorth(const Eigen::Vector2i& _gridCoord) const
        {
            if (0 <= _gridCoord.x() && _gridCoord.x() < grid->getResolution().x() &&
                0 <= _gridCoord.y() && _gridCoord.y() < grid->getResolution().y() - 1)
            {
                return openNorth[_gridCoord.y() * grid->getResolution().x() + _gridCoord.x()];
            }
            else
                return false;
        }

        /**
         * @brief Gets the is open flag for the south neighbor. If true, information may diffuse this way.
         * @param _gridCoord Grid coord at which the openness is viewed.
         * @return True if information can diffuse this way.
         */
        [[nodiscard]] bool getOpenSouth(const Eigen::Vector2i& _gridCoord) const
        {
            if (0 <= _gridCoord.x() && _gridCoord.x() < grid->getResolution().x() &&
                1 <= _gridCoord.y() && _gridCoord.y() < grid->getResolution().y())
            {
                return openNorth[(_gridCoord.y() - 1) * grid->getResolution().x() + _gridCoord.x()];
            }
            else
                return false;
        }

        /**
         * @brief Sets the type of grid point at a given grid coordinate.
         * @param _gridCoord Grid coordinate at which to set the grid point type.
         * @param _mask Grid point type to set.
         */
        void setMask(const Eigen::Vector2i& _gridCoord, EMask _mask)
        {
            mask[grid->getLinearIndex(_gridCoord)] = _mask;
        }

        /**
         * @brief Gets the type of grid point at a given grid coordinate.
         * @param _gridCoord Grid coordinate at which to get the grid point type.
         * @return Grid point type to get.
         */
        [[nodiscard]] EMask getMask(const Eigen::Vector2i& _gridCoord) const
        {
            return mask[grid->getLinearIndex(_gridCoord)];
        }

        /**
         * @brief Sets the source at a given grid coordinate. For Dirichlet grid points, this is the Dirichlet value. For all other pixels, this is the source value.
         * @param _gridCoord Grid coordinate at which to set the Dirichlet or source value.
         * @param _value Dirichlet or source value to set.
         */
        void setSource(const Eigen::Vector2i& _gridCoord, TValue _value)
        {
            source[grid->getLinearIndex(_gridCoord)] = _value;
        }

        /**
         * @brief Gets the source at a given grid coordinate. For Dirichlet grid points, this is the Dirichlet value. For all other pixels, this is the source value.
         * @param _gridCoord Grid coordinate at which to get the Dirichlet or source value.
         * @return Dirichlet or source value to get.
         */
        [[nodiscard]] TValue getSource(const Eigen::Vector2i& _gridCoord) const
        {
            return source[grid->getLinearIndex(_gridCoord)];
        }

        /**
         * @brief Gets the field buffer containing the final field value for each grid point.
         * @return Vector with solved field values.
         */
        [[nodiscard]] const std::vector<TValue>& getField() const
        {
            return field[ping];
        }

        /**
         * @brief Gets the field value at a given grid coordinate.
         * @param _gridCoord Grid coordinate at which to read the field value.
         * @return Final field value at requested grid point.
         */
        [[nodiscard]] const TValue& getField(const Eigen::Vector2i& _gridCoord) const
        {
            return field[ping][grid->getLinearIndex(_gridCoord)];
        }

        /**
         * @brief Sets the field value at a given grid coordinate.
         * @param _gridCoord Grid coordinate at which to read the field value.
         * @param _value Field value to set.
         */
        void setField(const Eigen::Vector2i& _gridCoord, const TValue& _value)
        {
            field[ping][grid->getLinearIndex(_gridCoord)] = _value;
        }

        /**
         * @brief Gets the grid over which the PDE is solved.
         * @return Grid over which the PDE is solved.
         */
        [[nodiscard]] std::shared_ptr<const RegularGrid2d> getGrid() const noexcept
        {
            return grid;
        }

    private:
        /**
         * @brief Ping pong state in {0,1}
         */
        int ping;

        /**
         * @brief Ping-pong buffer storing the field that is solved for.
         */
        std::vector<TValue> field[2];

        /**
         * @brief Grid over which the PDE is solved.
         */
        std::shared_ptr<RegularGrid2d> grid;

        /**
         * @brief Mask value per grid point, identifying what kind of grid point it is.
         */
        std::vector<EMask> mask;

        /**
         * @brief Contains Dirichlet conditions for Dirichlet grid points and source values for all other grid points.
         */
        std::vector<TValue> source;

        /**
         * @brief Staggered grid [X-1,Y] storing a flag that indicates whether information can flow to the east neighbor.
         */
        std::vector<bool> openEast;

        /**
         * @brief Staggered grid [X,Y-1] storing a flag that indicates whether information can flow to the north neighbor.
         */
        std::vector<bool> openNorth;
    };

    /**
     * @brief Two-dimensional Jacobi iteration solver for tri-variate fields.
     */
    using JacobiIteration3d = JacobiIteration<Eigen::Vector3d>;
}
