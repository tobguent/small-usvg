#pragma once

#include "bezier_surface.hpp"

namespace usvg
{
    /**
     * @brief Implementation of a Ferguson patch. Internally, the Ferguson patch stores a Bezier representation of its data.
     * @tparam TValue Type of values stored at the control points.
     */
    template <typename TValue>
    class FergusonPatch
    {
    public:
        /*

           03 - - - 33
            |       |
         v  |       |
            |       |
           00 - - - 30
                u
        */

        /**
         * @brief Enumeration of patch corners for which data is given.
         */
        enum EIndex
        {
            I00, // lower left corner (u=0, v=0)
            I30, // lower right corner (u=1, v=0)
            I03, // upper left corner (u=0, v=1)
            I33  // upper right corner (u=1, v=1)
        };

        /**
         * @brief Coordinate in uv space.
         */
        using DomainCoord = Eigen::Vector2d;

        /**
         * @brief Type of values stored at the control points.
         */
        using Value = TValue;

        /**
         * @brief Evaluates the Ferguson patch algebraically.
         * @param uv Domain location where to sample the surface.
         * @return Value of the patch.
         */
        [[nodiscard]] Value sample(const DomainCoord& uv) const noexcept
        {
            // Compute monomial basis vectors
            Eigen::RowVector4d uu(1, uv.x(), uv.x() * uv.x(), uv.x() * uv.x() * uv.x());
            Eigen::Vector4d vv(1, uv.y(), uv.y() * uv.y(), uv.y() * uv.y() * uv.y());
            // Coefficient matrix
            Eigen::Matrix4d C;
            C << 1, 0, -3, 2,
                0, 0, 3, -2,
                0, 1, -2, 1,
                0, 0, -1, 1;
            // Compute per component.
            Value result = Value::Zero();
            for (int dim = 0; dim < Value::RowsAtCompileTime; ++dim)
            {
                Eigen::Matrix4d Q;
                Q << value[EIndex::I00][dim], value[EIndex::I03][dim], tangentV[EIndex::I00][dim], tangentV[EIndex::I03][dim],
                    value[EIndex::I30][dim], value[EIndex::I33][dim], tangentV[EIndex::I30][dim], tangentV[EIndex::I33][dim],
                    tangentU[EIndex::I00][dim], tangentU[EIndex::I03][dim], 0, 0,
                    tangentU[EIndex::I30][dim], tangentU[EIndex::I33][dim], 0, 0;
                result[dim] = uu * C.transpose() * Q * C * vv;
            }
            return result;
        }

        /**
         * @brief Converts the Ferguson patch to a bicubic tensor product surface in Bezier Bernstein representation.
         */
        void updateBezierSurface() noexcept
        {
            bezierSurface.controlPoints[0][0] = value[EIndex::I00];                                                             // b[0,0]=p[0,0]
            bezierSurface.controlPoints[0][1] = value[EIndex::I00] + tangentV[EIndex::I00] / 3.0;                               // b[0,1]=(3*p[0,0]+dv[0,0])/3
            bezierSurface.controlPoints[0][2] = value[EIndex::I03] - tangentV[EIndex::I03] / 3.0;                               // b[0,2]=(3*p[0,1]-dv[0,1])/3
            bezierSurface.controlPoints[0][3] = value[EIndex::I03];                                                             // b[0,3]=p[0,1]
            bezierSurface.controlPoints[1][0] = value[EIndex::I00] + tangentU[EIndex::I00] / 3.0;                               // b[1,0]=(3*p[0,0]+du[0,0])/3
            bezierSurface.controlPoints[1][1] = value[EIndex::I00] + tangentV[EIndex::I00] / 3.0 + tangentU[EIndex::I00] / 3.0; // b[1,1]=(3*p[0,0]+dv[0,0]+du[0,0])/3
            bezierSurface.controlPoints[1][2] = value[EIndex::I03] - tangentV[EIndex::I03] / 3.0 + tangentU[EIndex::I03] / 3.0; // b[1,2]=(3*p[0,1]-dv[0,1]+du[0,1])/3
            bezierSurface.controlPoints[1][3] = value[EIndex::I03] + tangentU[EIndex::I03] / 3.0;                               // b[1,3]=(3*p[0,1]+du[0,1])/3
            bezierSurface.controlPoints[2][0] = value[EIndex::I30] - tangentU[EIndex::I30] / 3.0;                               // b[2,0]=(3*p[1,0]-du[1,0])/3
            bezierSurface.controlPoints[2][1] = value[EIndex::I30] + tangentV[EIndex::I30] / 3.0 - tangentU[EIndex::I30] / 3.0; // b[2,1]=(3*p[1,0]+dv[1,0]-du[1,0])/3
            bezierSurface.controlPoints[2][2] = value[EIndex::I33] - tangentV[EIndex::I33] / 3.0 - tangentU[EIndex::I33] / 3.0; // b[2,2]=(3*p[1,1]-dv[1,1]-du[1,1])/3
            bezierSurface.controlPoints[2][3] = value[EIndex::I33] - tangentU[EIndex::I33] / 3.0;                               // b[2,3]=(3*p[1,1]-du[1,1])/3
            bezierSurface.controlPoints[3][0] = value[EIndex::I30];                                                             // b[3,0]=p[1,0]
            bezierSurface.controlPoints[3][1] = value[EIndex::I30] + tangentV[EIndex::I30] / 3.0;                               // b[3,1]=(3*p[1,0]+dv[1,0])/3
            bezierSurface.controlPoints[3][2] = value[EIndex::I33] - tangentV[EIndex::I33] / 3.0;                               // b[3,2]=(3*p[1,1]-dv[1,1])/3
            bezierSurface.controlPoints[3][3] = value[EIndex::I33];                                                             // b[3,3]=p[1,1]
            bezierSurface.recomputeBoundingBox();
        }

        /**
         * @brief Gets the Bezier representation of the Ferguson patch. Note that updateBezierSurface() has to be called first.
         * @return Bicubic Bezier surface.
         */
        const BicubicBezierSurface<TValue>& getBezierSurface() const { return bezierSurface; }

        /**
         * @brief Checks if the Bezier representation is valid. Note that this function does not test if the Bezier corresponds to the values and tangents!
         * @return True if the Bezier representation is valid.
         */
        [[nodiscard]] bool isValid() const
        {
            // is Bezier representation valid?
            if (!bezierSurface.isValid())
                return false;

            // all good
            return true;
        }

        /**
         * @brief Values at the corners.
         */
        std::array<Value, 4> value;

        /**
         * @brief u-partial at the corners.
         */
        std::array<Value, 4> tangentU;

        /**
         * @brief v-partial at the corners.
         */
        std::array<Value, 4> tangentV;

    private:
        /**
         * @brief Bezier representation of the Ferguson patch.
         */
        BicubicBezierSurface<TValue> bezierSurface;
    };

    /**
     * @brief Uni-variate Ferguson patch.
     */
    using FergusonPatch1d = FergusonPatch<Eigen::Vector<double, 1>>;

    /**
     * @brief Bi-variate Ferguson patch.
     */
    using FergusonPatch2d = FergusonPatch<Eigen::Vector2d>;

    /**
     * @brief Tri-variate Ferguson patch.
     */
    using FergusonPatch3d = FergusonPatch<Eigen::Vector3d>;

    /**
     * @brief Tetra-variate Ferguson patch.
     */
    using FergusonPatch4d = FergusonPatch<Eigen::Vector4d>;
}
