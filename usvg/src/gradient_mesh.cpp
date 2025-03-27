#include <usvg/gradient_mesh.hpp>

#include <usvg/bezier_clipping.hpp>

namespace usvg
{
    void GradientMesh::Tile::sample_InverseCoordinatePartials(const DomainCoord& _uv, DomainCoord& _uv_x, DomainCoord& _uv_y, DomainCoord& _uv_xx, DomainCoord& _uv_yy, DomainCoord& _uv_xy) const noexcept
    {
        Eigen::Vector2d xxu  = position.getBezierSurface().sample_du(_uv);
        Eigen::Vector2d xxv  = position.getBezierSurface().sample_dv(_uv);
        Eigen::Vector2d xxuu = position.getBezierSurface().sample_duu(_uv);
        Eigen::Vector2d xxuv = position.getBezierSurface().sample_duv(_uv);
        Eigen::Vector2d xxvv = position.getBezierSurface().sample_dvv(_uv);

        double xu = xxu.x(), yu = xxu.y();
        double xv = xxv.x(), yv = xxv.y();
        double xuu = xxuu.x(), yuu = xxuu.y();
        double xvv = xxvv.x(), yvv = xxvv.y();
        double xuv = xxuv.x(), yuv = xxuv.y();

        Eigen::Matrix<double, 5, 5> M;
        M << xu, yu, 0, 0, 0,
            xv, yv, 0, 0, 0,
            xuu, yuu, xu * xu, yu * yu, 2 * xu * yu,
            xvv, yvv, xv * xv, yv * yv, 2 * xv * yv,
            xuv, yuv, xu * xv, yu * yv, xu * yv + xv * yu;
        auto Minv = M.inverse();

        _uv_x  = DomainCoord(Minv(0, 0), Minv(0, 1));
        _uv_y  = DomainCoord(Minv(1, 0), Minv(1, 1));
        _uv_xx = DomainCoord(Minv(2, 0), Minv(2, 1));
        _uv_yy = DomainCoord(Minv(3, 0), Minv(3, 1));
        _uv_xy = DomainCoord(Minv(4, 0), Minv(4, 1));
    }

    GradientMesh::Tile::Value GradientMesh::Tile::sample_ColorLaplacian_xy(const DomainCoord& _uv) const noexcept
    {
        DomainCoord uv_x, uv_y, uv_xx, uv_yy, uv_xy;
        sample_InverseCoordinatePartials(_uv, uv_x, uv_y, uv_xx, uv_yy, uv_xy);

        Eigen::Vector3d clr_u  = color.getBezierSurface().sample_du(_uv);
        Eigen::Vector3d clr_v  = color.getBezierSurface().sample_dv(_uv);
        Eigen::Vector3d clr_uu = color.getBezierSurface().sample_duu(_uv);
        Eigen::Vector3d clr_uv = color.getBezierSurface().sample_duv(_uv);
        Eigen::Vector3d clr_vv = color.getBezierSurface().sample_dvv(_uv);

        Eigen::Vector3d laplace = Eigen::Vector3d::Zero();
        for (int dim = 0; dim < 3; ++dim)
        {
            Eigen::Matrix2d H;
            H << clr_uu[dim], clr_uv[dim],
                clr_uv[dim], clr_vv[dim];
            Eigen::Vector2d grad(clr_u[dim], clr_v[dim]);
            laplace[dim] = uv_x.dot(H * uv_x) + grad.dot(uv_xx) + uv_y.dot(H * uv_y) + grad.dot(uv_yy);
        }
        return laplace;
    }

    bool GradientMesh::Tile::isValid() const
    {
        // check if Ferguson patches are valid.
        if (!position.isValid())
            return false;
        if (!color.isValid())
            return false;

        // all good
        return true;
    }

    GradientMesh::GradientMesh() noexcept
        : resolution(0, 0)
    {
    }

    GradientMesh::Tile& GradientMesh::getTile(Eigen::Index _i, Eigen::Index _j)
    {
        return tiles[_j * resolution.x() + _i];
    }

    const GradientMesh::Tile& GradientMesh::getTile(Eigen::Index _i, Eigen::Index _j) const
    {
        return tiles[_j * resolution.x() + _i];
    }

    bool GradientMesh::isValid() const
    {
        // is the number of position and color patches consistent with the number of rows and columns?
        if (tiles.size() != resolution.x() * resolution.y())
            return false;

        // check if all Ferguson patches are valid.
        for (auto& tile : tiles)
            if (!tile.isValid())
                return false;

        // all good
        return true;
    }

    const Eigen::Vector2d& GradientMesh::getPosition(Eigen::Index _i, Eigen::Index _j) const
    {
        if (_i == resolution.x() && _j == resolution.y())
            return getTile(_i - 1, _j - 1).position.value[FergusonPatch2d::I33];
        else if (_i != resolution.x() && _j == resolution.y())
            return getTile(_i, _j - 1).position.value[FergusonPatch2d::I03];
        else if (_i == resolution.x() && _j != resolution.y())
            return getTile(_i - 1, _j).position.value[FergusonPatch2d::I30];
        else
            return getTile(_i, _j).position.value[FergusonPatch2d::I00];
    }

    void GradientMesh::setPosition(Eigen::Index _i, Eigen::Index _j, const Eigen::Vector2d& _pos)
    {
        if (_i > 0 && _j > 0)
            getTile(_i - 1, _j - 1).position.value[FergusonPatch2d::I33] = _pos;
        if (_i > 0 && _j < resolution.y())
            getTile(_i - 1, _j).position.value[FergusonPatch2d::I30] = _pos;
        if (_i < resolution.x() && _j > 0)
            getTile(_i, _j - 1).position.value[FergusonPatch2d::I03] = _pos;
        if (_i < resolution.x() && _j < resolution.y())
            getTile(_i, _j).position.value[FergusonPatch2d::I00] = _pos;
    }

    const Eigen::Vector2d& GradientMesh::getPosition_du(Eigen::Index _i, Eigen::Index _j) const
    {
        if (_i == resolution.x() && _j == resolution.y())
            return getTile(_i - 1, _j - 1).position.tangentU[FergusonPatch2d::I33];
        else if (_i != resolution.x() && _j == resolution.y())
            return getTile(_i, _j - 1).position.tangentU[FergusonPatch2d::I03];
        else if (_i == resolution.x() && _j != resolution.y())
            return getTile(_i - 1, _j).position.tangentU[FergusonPatch2d::I30];
        else
            return getTile(_i, _j).position.tangentU[FergusonPatch2d::I00];
    }

    void GradientMesh::setPosition_du(Eigen::Index _i, Eigen::Index _j, const Eigen::Vector2d& _tangentU)
    {
        if (_i > 0 && _j > 0)
            getTile(_i - 1, _j - 1).position.tangentU[FergusonPatch2d::I33] = _tangentU;
        if (_i > 0 && _j < resolution.y())
            getTile(_i - 1, _j).position.tangentU[FergusonPatch2d::I30] = _tangentU;
        if (_i < resolution.x() && _j > 0)
            getTile(_i, _j - 1).position.tangentU[FergusonPatch2d::I03] = _tangentU;
        if (_i < resolution.x() && _j < resolution.y())
            getTile(_i, _j).position.tangentU[FergusonPatch2d::I00] = _tangentU;
    }

    const Eigen::Vector2d& GradientMesh::getPosition_dv(Eigen::Index _i, Eigen::Index _j) const
    {
        if (_i == resolution.x() && _j == resolution.y())
            return getTile(_i - 1, _j - 1).position.tangentV[FergusonPatch2d::I33];
        else if (_i != resolution.x() && _j == resolution.y())
            return getTile(_i, _j - 1).position.tangentV[FergusonPatch2d::I03];
        else if (_i == resolution.x() && _j != resolution.y())
            return getTile(_i - 1, _j).position.tangentV[FergusonPatch2d::I30];
        else
            return getTile(_i, _j).position.tangentV[FergusonPatch2d::I00];
    }

    void GradientMesh::setPosition_dv(Eigen::Index _i, Eigen::Index _j, const Eigen::Vector2d& _tangentV)
    {
        if (_i > 0 && _j > 0)
            getTile(_i - 1, _j - 1).position.tangentV[FergusonPatch2d::I33] = _tangentV;
        if (_i > 0 && _j < resolution.y())
            getTile(_i - 1, _j).position.tangentV[FergusonPatch2d::I30] = _tangentV;
        if (_i < resolution.x() && _j > 0)
            getTile(_i, _j - 1).position.tangentV[FergusonPatch2d::I03] = _tangentV;
        if (_i < resolution.x() && _j < resolution.y())
            getTile(_i, _j).position.tangentV[FergusonPatch2d::I00] = _tangentV;
    }

    void GradientMesh::updateBezierPosition(Eigen::Index _i, Eigen::Index _j)
    {
        if (_i > 0 && _j > 0)
            getTile(_i - 1, _j - 1).position.updateBezierSurface();
        if (_i > 0 && _j < resolution.y() - 1)
            getTile(_i - 1, _j).position.updateBezierSurface();
        if (_i < resolution.x() - 1 && _j > 0)
            getTile(_i, _j - 1).position.updateBezierSurface();
        if (_i < resolution.x() - 1 && _j < resolution.y() - 1)
            getTile(_i, _j).position.updateBezierSurface();
    }

    void GradientMesh::updateBezierPosition()
    {
        for (auto& tile : tiles)
            tile.position.updateBezierSurface();
    }

    const Eigen::Vector3d& GradientMesh::getColor(Eigen::Index _i, Eigen::Index _j) const
    {
        if (_i == resolution.x() && _j == resolution.y())
            return getTile(_i - 1, _j - 1).color.value[FergusonPatch2d::I33];
        else if (_i != resolution.x() && _j == resolution.y())
            return getTile(_i, _j - 1).color.value[FergusonPatch2d::I03];
        else if (_i == resolution.x() && _j != resolution.y())
            return getTile(_i - 1, _j).color.value[FergusonPatch2d::I30];
        else
            return getTile(_i, _j).color.value[FergusonPatch2d::I00];
    }

    void GradientMesh::setColor(Eigen::Index _i, Eigen::Index _j, const Eigen::Vector3d& _color)
    {
        if (_i > 0 && _j > 0)
            getTile(_i - 1, _j - 1).color.value[FergusonPatch2d::I33] = _color;
        if (_i > 0 && _j < resolution.y())
            getTile(_i - 1, _j).color.value[FergusonPatch2d::I30] = _color;
        if (_i < resolution.x() && _j > 0)
            getTile(_i, _j - 1).color.value[FergusonPatch2d::I03] = _color;
        if (_i < resolution.x() && _j < resolution.y())
            getTile(_i, _j).color.value[FergusonPatch2d::I00] = _color;
    }

    const Eigen::Vector3d& GradientMesh::getColor_du(Eigen::Index _i, Eigen::Index _j) const
    {
        if (_i == resolution.x() && _j == resolution.y())
            return getTile(_i - 1, _j - 1).color.tangentU[FergusonPatch2d::I33];
        else if (_i != resolution.x() && _j == resolution.y())
            return getTile(_i, _j - 1).color.tangentU[FergusonPatch2d::I03];
        else if (_i == resolution.x() && _j != resolution.y())
            return getTile(_i - 1, _j).color.tangentU[FergusonPatch2d::I30];
        else
            return getTile(_i, _j).color.tangentU[FergusonPatch2d::I00];
    }

    void GradientMesh::setColor_du(Eigen::Index _i, Eigen::Index _j, const Eigen::Vector3d& _tangentU)
    {
        if (_i > 0 && _j > 0)
            getTile(_i - 1, _j - 1).color.tangentU[FergusonPatch2d::I33] = _tangentU;
        if (_i > 0 && _j < resolution.y())
            getTile(_i - 1, _j).color.tangentU[FergusonPatch2d::I30] = _tangentU;
        if (_i < resolution.x() && _j > 0)
            getTile(_i, _j - 1).color.tangentU[FergusonPatch2d::I03] = _tangentU;
        if (_i < resolution.x() && _j < resolution.y())
            getTile(_i, _j).color.tangentU[FergusonPatch2d::I00] = _tangentU;
    }

    const Eigen::Vector3d& GradientMesh::getColor_dv(Eigen::Index _i, Eigen::Index _j) const
    {
        if (_i == resolution.x() && _j == resolution.y())
            return getTile(_i - 1, _j - 1).color.tangentV[FergusonPatch2d::I33];
        else if (_i != resolution.x() && _j == resolution.y())
            return getTile(_i, _j - 1).color.tangentV[FergusonPatch2d::I03];
        else if (_i == resolution.x() && _j != resolution.y())
            return getTile(_i - 1, _j).color.tangentV[FergusonPatch2d::I30];
        else
            return getTile(_i, _j).color.tangentV[FergusonPatch2d::I00];
    }

    void GradientMesh::setColor_dv(Eigen::Index _i, Eigen::Index _j, const Eigen::Vector3d& _tangentV)
    {
        if (_i > 0 && _j > 0)
            getTile(_i - 1, _j - 1).color.tangentV[FergusonPatch2d::I33] = _tangentV;
        if (_i > 0 && _j < resolution.y())
            getTile(_i - 1, _j).color.tangentV[FergusonPatch2d::I30] = _tangentV;
        if (_i < resolution.x() && _j > 0)
            getTile(_i, _j - 1).color.tangentV[FergusonPatch2d::I03] = _tangentV;
        if (_i < resolution.x() && _j < resolution.y())
            getTile(_i, _j).color.tangentV[FergusonPatch2d::I00] = _tangentV;
    }

    void GradientMesh::updateBezierColor(Eigen::Index _i, Eigen::Index _j)
    {
        if (_i > 0 && _j > 0)
            getTile(_i - 1, _j - 1).color.updateBezierSurface();
        if (_i > 0 && _j < resolution.y() - 1)
            getTile(_i - 1, _j).color.updateBezierSurface();
        if (_i < resolution.x() - 1 && _j > 0)
            getTile(_i, _j - 1).color.updateBezierSurface();
        if (_i < resolution.x() - 1 && _j < resolution.y() - 1)
            getTile(_i, _j).color.updateBezierSurface();
    }

    void GradientMesh::updateBezierColor()
    {
        for (auto& tile : tiles)
            tile.color.updateBezierSurface();
    }

    void GradientMesh::sampleColors(const Eigen::Vector2d& _position, int _maxDepth, double _epsilon, std::vector<Eigen::Vector3d>& _colors) const
    {
        // clear output vector
        _colors.clear();

        // for each tile
        for (auto& tile : tiles)
        {
            // locate the uv coordinate of the given position
            std::vector<Eigen::Vector2d> uvs;
            if (BezierClipping2d::locate(_position, tile.position.getBezierSurface(), _maxDepth, _epsilon, uvs))
            {
                // for each uv
                for (auto& uv : uvs)
                {
                    Eigen::Vector2d found_pos = tile.position.sample(uv);

                    // interpolate the color and add to list
                    _colors.push_back(tile.color.sample(uv));
                }
            }
        }
    }

    void GradientMesh::sampleColorLaplacian(const Eigen::Vector2d& _position, int _maxDepth, double _epsilon, std::vector<Eigen::Vector3d>& _colorLaplacians) const
    {
        // clear output vector
        _colorLaplacians.clear();

        // for each tile
        for (auto& tile : tiles)
        {
            // locate the uv coordinate of the given position
            std::vector<Eigen::Vector2d> uvs;
            if (BezierClipping2d::locate(_position, tile.position.getBezierSurface(), _maxDepth, _epsilon, uvs))
            {
                // for each uv
                for (auto& uv : uvs)
                {
                    Eigen::Vector2d found_pos = tile.position.sample(uv);

                    // interpolate the color and add to list
                    _colorLaplacians.push_back(tile.sample_ColorLaplacian_xy(uv));
                }
            }
        }
    }
}
