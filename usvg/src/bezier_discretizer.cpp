#include <usvg/bezier_discretizer.hpp>

namespace usvg
{
    std::pair<double, double> BezierDiscretizer2d::maxDist(const CubicBezierCurve2d& _curve)
    {
        // get implicit line through first and last point
        ImplicitLine2d implicit_signed(_curve.controlPoints.front(), _curve.controlPoints.back());
        auto [maxDist_sig, maxT_sig] = maxDist_signed(_curve.controlPoints, implicit_signed);

        // get left-rotated line through first point
        Eigen::Vector2d dir = _curve.controlPoints.back() - _curve.controlPoints.front();
        dir                 = Eigen::Vector2d(-dir.y(), dir.x());
        ImplicitLine2d implicit_left(_curve.controlPoints.front(), _curve.controlPoints.front() + dir);
        auto [maxDist_pos, maxT_pos] = maxDist_positive(_curve.controlPoints, implicit_left);

        // get right-rotated line through last point
        dir = _curve.controlPoints.back() - _curve.controlPoints.front();
        dir = Eigen::Vector2d(dir.y(), -dir.x());
        ImplicitLine2d implicit_right(_curve.controlPoints.back(), _curve.controlPoints.back() + dir);
        auto [maxDist_neg, maxT_neg] = maxDist_positive(_curve.controlPoints, implicit_right);

        // return the largest distance
        if (maxDist_sig >= maxDist_pos && maxDist_sig >= maxDist_neg)
            return std::make_pair(maxDist_sig, maxT_sig);
        else if (maxDist_pos >= maxDist_sig && maxDist_pos >= maxDist_neg)
            return std::make_pair(maxDist_pos, maxT_pos);
        else
            return std::make_pair(maxDist_neg, maxT_neg);
    }

    void BezierDiscretizer2d::topDown(const CubicBezierCurve2d& _curve, double _distanceThreshold, PiecewiseLinearCurve2d& _polyline)
    {
        // clear output
        _polyline.values.clear();
        _polyline.domain = _curve.domain;

        // does line connect onto itself?
        if ((_curve.controlPoints[0] - _curve.controlPoints[3]).stableNorm() < 1E-5)
        {
            // subdivide and discretize subparts
            CubicBezierCurve2d curve1, curve2;
            _curve.subdivide(curve1, curve2, 0.5);
            PiecewiseLinearCurve2d polyline1, polyline2;
            topDownRecursive(curve1, Eigen::Vector2d(_curve.domain.min()[0], _curve.domain.max()[0]), _distanceThreshold, polyline1);
            topDownRecursive(curve2, Eigen::Vector2d(_curve.domain.min()[0], _curve.domain.max()[0]), _distanceThreshold, polyline2);

            polyline1.values.removeLast();
            polyline1.parameters.removeLast();
            _polyline.values.append(polyline1.values);
            _polyline.values.append(polyline2.values);
            _polyline.parameters.append(polyline1.parameters);
            _polyline.parameters.append(polyline2.parameters);
        }
        else
        {
            topDownRecursive(_curve, Eigen::Vector2d(_curve.domain.min()[0], _curve.domain.max()[0]), _distanceThreshold, _polyline);
        }

        // append last point to finish the line strip
        _polyline.values.append(_curve.controlPoints.back());
        _polyline.parameters.append(_curve.domain.max());

        // recompute parameterization
        for (Eigen::Index ip = 0; ip < _polyline.parameters.getSize(); ++ip)
        {
            double t              = ip / (_polyline.parameters.getSize() - 1.0);
            Eigen::Vector1d param = _curve.domain.min() + t * (_curve.domain.max() - _curve.domain.min());
            _polyline.parameters.setValue(ip, param);
        }
        _polyline.setExplicitParameterization();

        _polyline.recomputeBoundingBox();
    }

    void BezierDiscretizer2d::topDown(const CubicBezierSpline2d& _spline, double _distanceThreshold, PiecewiseLinearCurve2d& _polyline)
    {
        // clear output
        _polyline.values.clear();
        _polyline.domain = _spline.domain;

        // loop over all curves
        int numCurves = _spline.getNumCurves();
        for (int icurve = 0; icurve < numCurves; ++icurve)
        {
            // discretize the curve
            auto curve = _spline.getCurve(icurve);
            PiecewiseLinearCurve2d discrete;
            topDown(curve, _distanceThreshold, discrete);
            // if not the last curve, remove the last vertex, since the next curve will start with the same
            if (icurve != numCurves - 1)
            {
                discrete.values.removeLast();
                discrete.parameters.removeLast();
            }
            // append at the end of the polyline.
            _polyline.values.append(discrete.values);
            _polyline.parameters.append(discrete.parameters);
        }
        _polyline.setExplicitParameterization();
        _polyline.recomputeBoundingBox();
    }

    std::pair<double, double> BezierDiscretizer2d::maxDist_signed(const std::array<Eigen::Vector2d, 4>& _nodes, const ImplicitLine2d& _implicit)
    {
        // compute the distance function to the implicit line in Bezier form
        Eigen::Vector4d nodesD = _implicit.distance(_nodes);

        // grab the distance values for node 1 and 2. for nodes 0 and 3, we know the distance is 0.
        double d1 = nodesD[1];
        double d2 = nodesD[2];
        if (std::abs(d1 - d2) < 1E-10)
        {
            double t1  = 0.5;
            double dt1 = 0.375 * (d1 + d2); // Eq. (7) simplifies quite a bit when t1=0.5 and when d1==d2
            return std::make_pair(std::abs(dt1), t1);
        }
        else
        {
            double desc  = d1 * d1 + d2 * d2 - d1 * d2;
            double t1    = (2 * d1 - d2 + std::sqrt(desc)) / (3 * (d1 - d2)); // Eq. (10)
            double tmax  = 0;
            double dtmax = 0;
            if (0 <= t1 && t1 <= 1)
            {
                double dt1 = 3 * t1 * (1 - t1) * ((1 - t1) * d1 + t1 * d2); // insert into Eq. (7)
                if (std::abs(dt1) > dtmax)
                {
                    dtmax = std::abs(dt1);
                    tmax  = t1;
                }
            }
            double t2 = (2 * d1 - d2 - std::sqrt(desc)) / (3 * (d1 - d2));
            if (0 <= t2 && t2 <= 1)
            {
                double dt2 = 3 * t2 * (1 - t2) * ((1 - t2) * d1 + t2 * d2);
                if (std::abs(dt2) > dtmax)
                {
                    dtmax = std::abs(dt2);
                    tmax  = t2;
                }
            }
            return std::make_pair(dtmax, tmax);
        }
    }

    std::pair<double, double> BezierDiscretizer2d::maxDist_positive(const std::array<Eigen::Vector2d, 4>& _nodes, const ImplicitLine2d& _implicit)
    {
        // compute the distance function to the implicit line in Bezier form
        Eigen::Vector4d nodesD = _implicit.distance(_nodes);

        // we need at least a positive number, which is why dtmax=0 initially. tmax is initially infinite to know if no point is moving into the positive sector of the implicit line
        double tmax  = std::numeric_limits<double>::max();
        double dtmax = 0;

        // grab the distance values for node 1 and 2. for nodes 0 and 3, we know one distance is 0 and the other is on the wrong side.
        double d1 = nodesD[1];
        double d2 = nodesD[2];
        if (std::abs(d1 - d2) < 1E-10)
        {
            double t1  = 0.5;
            double dt1 = 0.375 * (d1 + d2); // Eq. (7) simplifies quite a bit when t1=0.5 and when d1==d2
            if (dt1 > dtmax)
            {
                dtmax = dt1;
                tmax  = t1;
            }
        }
        else
        {
            double desc = d1 * d1 + d2 * d2 - d1 * d2;
            double t1   = (2 * d1 - d2 + std::sqrt(desc)) / (3 * (d1 - d2)); // Eq. (10)
            if (0 <= t1 && t1 <= 1)
            {
                double dt1 = 3 * t1 * (1 - t1) * ((1 - t1) * d1 + t1 * d2); // insert into Eq. (7)
                if (dt1 > dtmax)
                {
                    dtmax = dt1;
                    tmax  = t1;
                }
            }
            double t2 = (2 * d1 - d2 - std::sqrt(desc)) / (3 * (d1 - d2));
            if (0 <= t2 && t2 <= 1)
            {
                double dt2 = 3 * t2 * (1 - t2) * ((1 - t2) * d1 + t2 * d2);
                if (dt2 > dtmax)
                {
                    dtmax = dt2;
                    tmax  = t2;
                }
            }
        }
        return std::make_pair(dtmax, tmax);
    }

    void BezierDiscretizer2d::topDownRecursive(const CubicBezierCurve2d& _curve, Eigen::Vector2d _range, double _distanceThreshold, PiecewiseLinearCurve2d& _polyline)
    {
        CubicBezierCurve2d curve = _curve;
        curve.domain             = Eigen::AlignedBox1d(0, 1);

        // compute largest distance from implicit line going through start and end, and get the relative position in [0,1] where this occurs
        auto [dist, t] = maxDist(curve);

        // subdivide the line if dist is too large
        if (dist > _distanceThreshold)
        {
            // determine the domain location where to subdivide
            // double tt = curve.domain.min()[0] + t * (curve.domain.max()[0] - curve.domain.min()[0]);

            // subdivide with de Casteljau
            CubicBezierCurve2d clip1, clip2;
            curve.subdivide(clip1, clip2, t);

            // continue recursively
            topDownRecursive(clip1, Eigen::Vector2d(_range.x(), _range.x() + t * (_range.y() - _range.x())), _distanceThreshold, _polyline);
            topDownRecursive(clip2, Eigen::Vector2d(_range.x() + t * (_range.y() - _range.x()), _range.y()), _distanceThreshold, _polyline);
        }
        else
        {
            // add line to output (only the start point, since we get by construction a line strip)
            _polyline.values.append(curve.controlPoints.front());
            _polyline.parameters.append(curve.domain.min());
        }
    }
}
