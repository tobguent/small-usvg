#include <usvg/bezier_clipping.hpp>

#include <usvg/implicit_line.hpp>

namespace usvg
{
    bool BezierClipping2d::locate(const Eigen::Vector2d& _point, const BicubicBezierSurface2d& _surface, int _maxDepth, double _epsilon, std::vector<Eigen::Vector2d>& _uvs)
    {
        // clear content of output vector
        _uvs.clear();

        // do the bounding boxes overlap?
        if (!_surface.getBoundingBox().contains(_point))
            return false;

        // get ranges and call recursion
        Eigen::Vector2d rangeU(_surface.domain.min()[0], _surface.domain.max()[0]);
        Eigen::Vector2d rangeV(_surface.domain.min()[1], _surface.domain.max()[1]);
        return locateRecursive(_point, _surface, rangeU, rangeV, _maxDepth, _epsilon, _surface, _uvs);
    }

    std::pair<double, double> BezierClipping2d::clip(const Eigen::Matrix4d& _nodesD)
    {
        // reset clip range
        double smin = 1;
        double smax = 0;

        // traverse the rows if the control polygon and collect the smallest and largest distance for each column.
        Eigen::Vector4d minD(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
        Eigen::Vector4d maxD(-std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(), -std::numeric_limits<double>::max());
        for (int irow = 0; irow < 4; ++irow)
            for (int icol = 0; icol < 4; ++icol)
            {
                minD[icol] = std::min(minD[icol], _nodesD(irow, icol));
                maxD[icol] = std::max(maxD[icol], _nodesD(irow, icol));
            }

        // if the first column has a positive and negative distance in it, then this left boundary of the convex hull intersects the x-axis at 0.
        if (minD[0] < 0 && 0 < maxD[0])
            smin = 0;
        // same for the last column
        if (minD[3] < 0 && 0 < maxD[3])
            smax = 1;

        // Let's call the entries of nodesD = {e_ij}
        // The columns of nodesD are mapped to t-values[t1, t2, t3, t4] = [0, 1/3, 2/3, 1]
        // This gives a coordinate for each point in nodesD: (t1,e_11), (t1,e21), (t1,e31), (t1,e41), (t2,e12), ..., (t4,e44)
        // We need to intersect the convex hull of these 16 points with the x-axis at y=0.
        // Instead of finding the convex hull, we get the lowest and highest values per column, and test each connecting edge among the loweest values, and the connecting edge among the highest values per column.
        // This will include the edges of the convex hull, but it is also contains a few unnecessary checks.
        for (int i = 0; i < 3; ++i)
            for (int j = i + 1; j < 4; ++j)
            {
                // the two points of the segment are(t1, e1) and (t2, e2)
                double t1 = i / 3.0;
                double t2 = j / 3.0;
                double e1 = minD[i];
                double e2 = minD[j];
                // horizontal segment cannot intersect with the line
                if (e1 != e2)
                {
                    // does the segment intersect with 0?
                    double u = (0 - e1) / (e2 - e1);
                    if (0 <= u && u <= 1)
                    {
                        double t = u * (t2 - t1) + t1;
                        smin     = std::min(t, smin);
                        smax     = std::max(t, smax);
                    }
                }
                e1 = maxD[i];
                e2 = maxD[j];
                // horizontal segment cannot intersect with the line dmin or the line dmax
                if (e1 != e2)
                {
                    // does the segment intersect with 0?
                    double u = (0 - e1) / (e2 - e1);
                    if (0 <= u && u <= 1)
                    {
                        double t = u * (t2 - t1) + t1;
                        smin     = std::min(t, smin);
                        smax     = std::max(t, smax);
                    }
                }
            }
        return std::make_pair(smin, smax);
    }

    bool BezierClipping2d::locateRecursive(const Eigen::Vector2d& _point, const BicubicBezierSurface2d& _surface, Eigen::Vector2d _rangeU, Eigen::Vector2d _rangeV, int _depth, double _epsilon, const BicubicBezierSurface2d& _initial_surface, std::vector<Eigen::Vector2d>& _uvs)
    {
        // by default, we toggle between subdividing in u and v direction
        bool subdivide_v_next = (_depth % 2) == 0;

        // have we finished yet in the u or v directions?
        bool done_with_u = std::abs(_rangeU.y() - _rangeU.x()) < _epsilon;
        bool done_with_v = std::abs(_rangeV.y() - _rangeV.x()) < _epsilon;

        // narrowed down far enough?
        if (done_with_u && done_with_v)
        {
            _uvs.push_back(Eigen::Vector2d(_rangeU.mean(), _rangeV.mean()));
            return true;
        }

        // not successful; recursion went too deep
        if (_depth == 0)
            return false;

        // if finished in one direction, continue on the other.
        if (done_with_u)
            subdivide_v_next = true;
        if (done_with_v)
            subdivide_v_next = false;

        // make local copy of the curves and set the range to [0,1]
        BicubicBezierSurface2d surface = _surface;
        surface.domain                 = Eigen::AlignedBox2d(Eigen::Vector2d(0, 0), Eigen::Vector2d(1, 1));

        // get offset to define the line orientation
        Eigen::Vector2d offset;
        // align the implicit line with the control polygon. See Eq.(12) in Alexander Efremov, Vlastimil Havran, and Hans-Peter Seidel. 2005. Robust and numerically stable Bézier clipping method for ray tracing NURBS surfaces. In Proceedings of the 21st Spring Conference on Computer Graphics (SCCG '05). Association for Computing Machinery, New York, NY, USA, 127–135. https://doi.org/10.1145/1090122.1090144
        if (subdivide_v_next)
        {
            // horizontal line: subdivide v
            Eigen::Vector2d upperEdge = surface.controlPoints[3][3] - surface.controlPoints[0][3];
            Eigen::Vector2d lowerEdge = surface.controlPoints[3][0] - surface.controlPoints[0][0];
            offset                    = (upperEdge + lowerEdge) / 2;
            if (offset.stableNorm() < 1E-10)
            {
                // fall back to initial directions in case the surface already collapsed too far upon convergence
                upperEdge = _initial_surface.controlPoints[3][3] - _initial_surface.controlPoints[0][3];
                lowerEdge = _initial_surface.controlPoints[3][0] - _initial_surface.controlPoints[0][0];
                offset    = (upperEdge + lowerEdge) / 2;
            }
        }
        else
        {
            // vertical line: subdivide u
            Eigen::Vector2d leftEdge  = surface.controlPoints[0][3] - surface.controlPoints[0][0];
            Eigen::Vector2d rightEdge = surface.controlPoints[3][3] - surface.controlPoints[3][0];
            offset                    = (leftEdge + rightEdge) / 2;
            if (offset.stableNorm() < 1E-10)
            {
                // fall back to initial directions in case the surface already collapsed too far upon convergence
                leftEdge  = _initial_surface.controlPoints[0][3] - _initial_surface.controlPoints[0][0];
                rightEdge = _initial_surface.controlPoints[3][3] - _initial_surface.controlPoints[3][0];
                offset    = (leftEdge + rightEdge) / 2;
            }
        }
        // define an implicit line.
        ImplicitLine2d line(_point, _point + offset);

        // compute the distance to all control points.
        Eigen::Vector4d dist0 = line.distance(surface.controlPoints[0]);
        Eigen::Vector4d dist1 = line.distance(surface.controlPoints[1]);
        Eigen::Vector4d dist2 = line.distance(surface.controlPoints[2]);
        Eigen::Vector4d dist3 = line.distance(surface.controlPoints[3]);
        Eigen::Matrix4d dists;
        dists << dist0, dist1, dist2, dist3;
        if (subdivide_v_next)
        {
            dists.transposeInPlace();
        }

        // all control points of nodesD are above / below the threshold
        bool any_above = false, any_below = false;
        for (int i = 0; i < 4 && !(any_above && any_below); ++i)
            for (int j = 0; j < 4 && !(any_above && any_below); ++j)
            {
                if (dists(i, j) > 0)
                    any_above = true;
                if (dists(i, j) < 0)
                    any_below = true;
            }
        if (!any_above || !any_below)
        {
            // done, certainly not an intersection here
            return false;
        }
        else
        {
            // find range for clipping of surface
            auto [smin, smax] = BezierClipping2d::clip(dists);

            // no progress?
            // if (smin == 0 && smax == 1)
            // if (smin <= 1E-10 && smax >= 1 - 1E-10)
            if (smax - smin > 0.9) // reduction is less than 10%
            {
                if (done_with_u && !done_with_v)
                {
                    // subdivide the surface and try recursively on the four pieces
                    BicubicBezierSurface2d clip0, clip1;
                    surface.subdivide_v(clip0, clip1, 0.5);
                    bool side1 = locateRecursive(_point, clip0, Eigen::Vector2d(_rangeU.x(), _rangeU.y()), Eigen::Vector2d(_rangeV.x(), _rangeV.mean()), _depth - 1, _epsilon, _initial_surface, _uvs);
                    bool side2 = locateRecursive(_point, clip1, Eigen::Vector2d(_rangeU.x(), _rangeU.y()), Eigen::Vector2d(_rangeV.mean(), _rangeV.y()), _depth - 1, _epsilon, _initial_surface, _uvs);
                    return side1 || side2;
                }

                if (!done_with_u && done_with_v)
                {
                    // subdivide the surface and try recursively on the four pieces
                    BicubicBezierSurface2d clip0, clip1;
                    surface.subdivide_u(clip0, clip1, 0.5);
                    bool side1 = locateRecursive(_point, clip0, Eigen::Vector2d(_rangeU.x(), _rangeU.mean()), Eigen::Vector2d(_rangeV.x(), _rangeV.y()), _depth - 1, _epsilon, _initial_surface, _uvs);
                    bool side2 = locateRecursive(_point, clip1, Eigen::Vector2d(_rangeU.mean(), _rangeU.y()), Eigen::Vector2d(_rangeV.x(), _rangeV.y()), _depth - 1, _epsilon, _initial_surface, _uvs);
                    return side1 || side2;
                }

                // subdivide the surface and try recursively on the four pieces
                BicubicBezierSurface2d clip00, clip01, clip10, clip11;
                surface.subdivide(clip00, clip01, clip10, clip11, Eigen::Vector2d(0.5, 0.5));
                bool side1 = locateRecursive(_point, clip00, Eigen::Vector2d(_rangeU.x(), _rangeU.mean()), Eigen::Vector2d(_rangeV.x(), _rangeV.mean()), _depth - 1, _epsilon, _initial_surface, _uvs);
                bool side2 = locateRecursive(_point, clip01, Eigen::Vector2d(_rangeU.x(), _rangeU.mean()), Eigen::Vector2d(_rangeV.mean(), _rangeV.y()), _depth - 1, _epsilon, _initial_surface, _uvs);
                bool side3 = locateRecursive(_point, clip10, Eigen::Vector2d(_rangeU.mean(), _rangeU.y()), Eigen::Vector2d(_rangeV.x(), _rangeV.mean()), _depth - 1, _epsilon, _initial_surface, _uvs);
                bool side4 = locateRecursive(_point, clip11, Eigen::Vector2d(_rangeU.mean(), _rangeU.y()), Eigen::Vector2d(_rangeV.mean(), _rangeV.y()), _depth - 1, _epsilon, _initial_surface, _uvs);
                return side1 || side2 || side3 || side4;
            }
            else
            {
                // adjust global range, since [smin,smax] are in [0,1] range
                if (subdivide_v_next)
                {
                    _rangeV = Eigen::Vector2d(
                        _rangeV[0] + (_rangeV[1] - _rangeV[0]) * smin,
                        _rangeV[0] + (_rangeV[1] - _rangeV[0]) * smax);
                }
                else
                {
                    _rangeU = Eigen::Vector2d(
                        _rangeU[0] + (_rangeU[1] - _rangeU[0]) * smin,
                        _rangeU[0] + (_rangeU[1] - _rangeU[0]) * smax);
                }

                // subdivide if necessary
                if (smin > 0)
                {
                    BicubicBezierSurface2d clip1, clip2;
                    if (subdivide_v_next)
                        surface.subdivide_v(clip1, clip2, smin);
                    else
                        surface.subdivide_u(clip1, clip2, smin);
                    surface = clip2;
                }
                if (smax < 1)
                {
                    BicubicBezierSurface2d clip1, clip2;
                    if (subdivide_v_next)
                        surface.subdivide_v(clip1, clip2, smax);
                    else
                        surface.subdivide_u(clip1, clip2, smax);
                    surface = clip1;
                }
                surface.domain = Eigen::AlignedBox2d(Eigen::Vector2d(0, 0), Eigen::Vector2d(1, 1));
            }
        }

        // call recursively
        return locateRecursive(_point, surface, _rangeU, _rangeV, _depth - 1, _epsilon, _initial_surface, _uvs);
    }
}
