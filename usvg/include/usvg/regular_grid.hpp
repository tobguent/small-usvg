#pragma once

#include <Eigen/Eigen>

namespace usvg
{
    /**
     * @brief Basic interface for regular grids in a certain dimension.
     * @tparam TScalar Scalar value used to store data on the grid.
     * @tparam TDimensions Number of dimensions.
     */
    template <typename TScalar, int64_t TDimensions>
    class RegularGrid
    {
    public:
        /**
         * @brief Number of dimensions.
         */
        static const int64_t Dimensions = TDimensions;

        /**
         * @brief Scalar type of domain coordinates.
         */
        using Scalar = TScalar;

        /**
         * @brief Grid index coordinate.
         */
        using GridCoord = typename Eigen::Vector<int, TDimensions>;

        /**
         * @brief Axis-aligned bounding box type.
         */
        using BoundingBox = typename Eigen::AlignedBox<TScalar, TDimensions>;

        /**
         * @brief Physical domain coordinate.
         */
        using DomainCoord = typename Eigen::Vector<TScalar, TDimensions>;

        /**
         * @brief Stride type.
         */
        using Strides = typename Eigen::Vector<int, TDimensions>;

        /**
         * @brief Gets the grid resolution.
         * @return Resolution of the grid.
         */
        [[nodiscard]] inline const GridCoord& getResolution() const { return this->mResolution; }

        /**
         * @brief Sets the grid resolution.
         * @param resolution Resolution to set.
         */
        inline void setResolution(const GridCoord& resolution)
        {
            this->mResolution = resolution;

            // compute the strides
            computeStrides();
        }

        /**
         * @brief Gets the domain bounding box.
         * @return Bounding box of domain.
         */
        [[nodiscard]] inline const BoundingBox& getDomain() const { return this->mDomain; }

        /**
         * @brief Sets the domain bounding box.
         * @param domain Bounding box to set.
         */
        inline void setDomain(const BoundingBox& domain) { this->mDomain = domain; }

        /**
         * @brief Gets the total number of grid points.
         * @return Total number of grid points.
         */
        [[nodiscard]] inline Eigen::Index getNumGridPoints() const
        {
            return mResolution.prod();
        }

        /**
         * @brief Gets the physical coordinate of a grid point, identified by its linear index.
         * @param linearIndex Linear index of grid point.
         * @return Physical domain coordinate of grid point.
         */
        [[nodiscard]] inline DomainCoord getCoordAt(Eigen::Index linearIndex) const
        {
            return getCoordAt(getGridCoord(linearIndex));
        }

        /**
         * @brief Gets the linear array index based on a grid coordinate index.
         * @param gridCoord Grid coordinate.
         * @return Corresponding linear index.
         */
        [[nodiscard]] Eigen::Index getLinearIndex(const GridCoord& gridCoord) const
        {
            return gridCoord.dot(this->mPointStrides);
        }

        /**
         * @brief Gets the spatial location of a grid vertex.
         * @param gridCoord Grid coordinate.
         * @return Physical domain coordiante of grid point.
         */
        [[nodiscard]] DomainCoord getCoordAt(const GridCoord& gridCoord) const
        {
            DomainCoord s;
            for (int i = 0; i < Dimensions; ++i)
            {
                s[i] = mResolution[i] < 2 ? 0.5 : gridCoord[i] / (this->mResolution[i] - Scalar(1.));
            }
            return this->mDomain.min() + (this->mDomain.max() - this->mDomain.min()).cwiseProduct(s);
        }

        /**
         * @brief Gets the grid coordinate based on the linear array index.
         * @param linearIndex Linear array index.
         * @return Corresponding grid coordinate.
         */
        [[nodiscard]] GridCoord getGridCoord(Eigen::Index linearIndex) const
        {
            return getGridCoord(linearIndex, mResolution, mPointStrides);
        }

        /**
         * @brief Gets the number of cells.
         * @return Number of cells.
         */
        [[nodiscard]] Eigen::Index getNumCells() const
        {
            return (this->mResolution - GridCoord::Ones()).prod();
        }

        /**
         * @brief Gets the linear grid point indices of a given cell.
         * @param cellIndex Index of the cell to get.
         * @return Grid point coordinates of the cell.
         */
        [[nodiscard]] Eigen::VectorX<Eigen::Index> getCell(Eigen::Index cellIndex) const
        {
            assert(cellIndex < getNumCells());

            // find the grid index of the "lower-left" corner in cell coordinates.
            GridCoord cellOrigin = getGridCoord(cellIndex, mResolution - GridCoord::Ones(), mCellStrides);

            // collect the corners
            const int64_t numCorners = 1 << Dimensions;
            Eigen::VectorX<Eigen::Index> result;
            result.resize(numCorners);
            for (int64_t i = 0; i < numCorners; ++i)
            {
                GridCoord currCoord = cellOrigin;
                for (int64_t dim = 0; dim < Dimensions; ++dim)
                {
                    currCoord[dim] += (i & (1ll << dim)) ? 1 : 0;
                }
                result[i] = getLinearIndex(currCoord);
            }
            return result;
        }

        /**
         * @brief Gets the strides per dimension for linear indexing of grid points.
         * @return Strides for iterating grid points.
         */
        [[nodiscard]] const Strides& getPointStrides() const
        {
            return mPointStrides;
        }

        /**
         * @brief Gets the strides per dimension for linear indexing of cells.
         * @return Strides for iterating cells.
         */
        [[nodiscard]] const Strides& getCellStrides() const
        {
            return mCellStrides;
        }

        /**
         * @brief Gets the origin of the domain, i.e., the "lower-left" corner of the domain.
         * @return Origin of the domain.
         */
        [[nodiscard]] inline DomainCoord getOrigin() const
        {
            return mDomain.min();
        }

        /**
         * @brief Gets the spacing between adjacent grid points.
         * @return Spacing of the domain.
         */
        [[nodiscard]] inline DomainCoord getSpacing() const
        {
            return (mDomain.max() - mDomain.min()).cwiseQuotient((mResolution - GridCoord::Ones()).template cast<TScalar>());
        }

        /**
         * @brief Tests if the resolution is set properly (no dimension is zero) and that the domain is not invalid (zero domain). If this function returns false, the object has probably not been initialized correctly with setResolution and setDomain.
         * @return True if the grid is valid.
         */
        [[nodiscard]] bool isValid() const
        {
            for (int64_t d = 0; d < Dimensions; ++d)
                if (mResolution[d] <= 0)
                    return false;
            return !mDomain.isNull();
        }

        /**
         * @brief Tests if the other grid has the same domain and the same resolution.
         * @param other Other grid to compare to.
         * @return True if the domain and the resolution are the same.
         */
        [[nodiscard]] bool isEqual(const RegularGrid& other) const
        {
            return mDomain.min() == other.mDomain.min() && mDomain.max() == other.mDomain.max() && mResolution == other.mResolution;
        }

    private:
        /**
         * @brief Computes the strides for a given resolution.
         * @param resolution Resolution to compute strides for.
         * @return Stride for computing linear indices.
         */
        [[nodiscard]] static Strides computeStrides(const GridCoord& resolution)
        {
            Strides strides;
            strides[0] = 1;
            for (int64_t d = 0; d < Dimensions - 1; ++d)
            {
                strides[d + 1] = strides[d] * resolution[d];
            }
            return strides;
        }

        /**
         * @brief Computes the strides for the current resolution.
         */
        void computeStrides()
        {
            this->mPointStrides = computeStrides(this->mResolution);
            this->mCellStrides  = computeStrides(this->mResolution - GridCoord::Ones());
        }

        /**
         * @brief Gets the grid coordinate based on the linear array index for a given resolution and precomputed strides.
         * @param linearIndex Linear array index.
         * @param resolution Resolution of the grid.
         * @param strides Precomputed strides.
         * @return Corresponding grid coordinate.
         */
        [[nodiscard]] static GridCoord getGridCoord(Eigen::Index linearIndex, const GridCoord& resolution, const Strides& strides)
        {
            if constexpr (Dimensions == 2)
            {
                return GridCoord(linearIndex % resolution.x(), linearIndex / resolution.x());
            }
            else
            {
                GridCoord result;
                Eigen::Index t = linearIndex;
                for (int64_t d = Dimensions - 1; d >= 0; --d)
                {
                    result[d] = (int)(t / strides[d]);
                    t         = t % strides[d];
                }
                return result;
            }
        }

        /**
         * @brief Domain over which the regular grid is defined.
         */
        BoundingBox mDomain;

        /**
         * @brief Grid point resolution of this regular grid.
         */
        GridCoord mResolution;

        /**
         * @brief Strides per dimension for linear indexing of grid points.
         */
        Strides mPointStrides;

        /**
         * @brief Strides per dimension for linear indexing of cells.
         */
        Strides mCellStrides;
    };

    /**
     * @brief One-dimensional regular grid.
     */
    using RegularGrid1d = RegularGrid<double, 1>;

    /**
     * @brief Two-dimensional regular grid.
     */
    using RegularGrid2d = RegularGrid<double, 2>;

    /**
     * @brief Three-dimensional regular grid.
     */
    using RegularGrid3d = RegularGrid<double, 3>;

    /**
     * @brief Four-dimensional regular grid.
     */
    using RegularGrid4d = RegularGrid<double, 4>;
}
