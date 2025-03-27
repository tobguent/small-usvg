#pragma once

#include <usvg/bvh.hpp>

#include <algorithm>
#include <assert.h>

#include <Eigen/Eigen>

#include <stack>

namespace usvg
{
    void Bvh2d::build(const PiecewiseLinearCurve2d& _segments)
    {
        segments = _segments;

        int num_segments = (int)segments.values.getSize() - 1;
        for (int k = 0; k < 2; ++k)
            sorted_segments[k].resize(num_segments);
        is_left.resize(num_segments);

        // init sorted index lists
        for (int i = 0; i < num_segments; i++)
        {
            for (int k = 0; k < 2; ++k)
                sorted_segments[k][i] = i;
        }

        for (int k = 0; k < 2; ++k)
            sort_segments(sorted_segments[k], 0, num_segments - 1, k);

        nodes.clear();
        nodes.push_back(node()); // insert root

        for (int i = 0; i < num_segments; i++)
        {
            nodes[0].bounds.extend(segments.values.getValue(i));
            nodes[0].bounds.extend(segments.values.getValue(i + 1));
        }

        create_hierarchy(0, nodes[0].bounds, sorted_segments, num_segments);
        recompute_bounds(nodes[0]);

        is_left.clear();
    }

    Bvh2d::node::node()
        : left(invalid)
        , right(invalid)
    {
        bounds.setEmpty();
    }

    bool Bvh2d::intersect_segment_segment(const Eigen::Vector2d& _a, const Eigen::Vector2d& _b, const Eigen::Vector2d& _c, const Eigen::Vector2d& _d, Eigen::Vector2d& _t, bool testLeft, bool testRight)
    {
        // if one side is not being tested
        if (!testLeft || !testRight)
        {
            // determine the normal
            Eigen::Vector2d normal(_b.y() - _a.y(), _a.x() - _b.x());
            double dot = normal.dot(_d - _c);
            if ((!testLeft && dot > 0) || (!testRight && dot < 0))
                return false;
        }

        Eigen::Vector2d ba = _b - _a, dc = _d - _c;
        double disc = _a.x() * (_d.y() - _c.y()) + _b.x() * (_c.y() - _d.y()) + (_b.y() - _a.y()) * _d.x() + (_a.y() - _b.y()) * _c.x();
        if (abs(disc) < 1E-10)
            return false;

        _t.x() = (_a.x() * (_d.y() - _c.y()) + _c.x() * (_a.y() - _d.y()) + (_c.y() - _a.y()) * _d.x()) / (disc != 0 ? disc : 1);
        if (_t.x() < 0 || 1 < _t.x())
            return false;

        _t.y() = -(_a.x() * (_c.y() - _b.y()) + _b.x() * (_a.y() - _c.y()) + (_b.y() - _a.y()) * _c.x()) / (disc != 0 ? disc : 1);
        return 0 <= _t.y() && _t.y() <= 1;
    }

    bool Bvh2d::any_segment_segment(const Eigen::Vector2d& _b0, const Eigen::Vector2d& _b1, const Eigen::Vector2d& _a0, const Eigen::Vector2d& _a1, bool testLeft, bool testRight)
    {
        Eigen::Vector2d t;
        return intersect_segment_segment(_b0, _b1, _a0, _a1, t, testLeft, testRight);
    }

    bool Bvh2d::any_segment_AABB(const Eigen::Vector2d& _b0, const Eigen::Vector2d& _b1, const Eigen::AlignedBox2d& bounds)
    {
        Eigen::Vector2d min_bounds(std::min(_b0.x(), _b1.x()), std::min(_b0.y(), _b1.y()));
        Eigen::Vector2d max_bounds(std::max(_b0.x(), _b1.x()), std::max(_b0.y(), _b1.y()));
        if ((min_bounds.x() > bounds.max().x()) || (min_bounds.y() > bounds.max().y()) || (max_bounds.x() < bounds.min().x()) || (max_bounds.y() < bounds.min().y()))
            return false;
        if ((bounds.min().x() <= _b0.x()) && (bounds.min().y() <= _b0.y()) &&
            (_b0.x() <= bounds.max().x()) && (_b0.y() <= bounds.max().y()) &&
            (bounds.min().x() <= _b1.x()) && (bounds.min().y() <= _b1.y()) &&
            (_b1.x() <= bounds.max().x()) && (_b1.y() <= bounds.max().y()))
            return true;
        if (any_segment_segment(_b0, _b1, Eigen::Vector2d(bounds.min().x(), bounds.min().y()), Eigen::Vector2d(bounds.max().x(), bounds.min().y()), true, true))
            return true;
        if (any_segment_segment(_b0, _b1, Eigen::Vector2d(bounds.max().x(), bounds.min().y()), Eigen::Vector2d(bounds.max().x(), bounds.max().y()), true, true))
            return true;
        if (any_segment_segment(_b0, _b1, Eigen::Vector2d(bounds.min().x(), bounds.max().y()), Eigen::Vector2d(bounds.max().x(), bounds.max().y()), true, true))
            return true;
        if (any_segment_segment(_b0, _b1, Eigen::Vector2d(bounds.min().x(), bounds.min().y()), Eigen::Vector2d(bounds.min().x(), bounds.max().y()), true, true))
            return true;
        return false;
    }

    void Bvh2d::closest_point_on_segment(const Eigen::Vector2d& _b0, const Eigen::Vector2d& _b1, const Eigen::Vector2d& _p, int _curve, double& _min_dist, double& _min_t, int& _min_curve)
    {
        Eigen::Vector2d u = _b1 - _b0;
        double t          = (_p - _b0).dot(u) / u.dot(u);
        t                 = std::min(std::max(0., t), 1.);
        Eigen::Vector2d q = (1 - t) * _b0 + t * _b1;
        double dist       = (_p - q).stableNorm();
        if (dist < _min_dist)
        {
            _min_dist  = dist;
            _min_t     = t;
            _min_curve = _curve;
        }
    }

    bool Bvh2d::anyHit(const Eigen::Vector2d& _a0, const Eigen::Vector2d& _a1, bool testLeft, bool testRight) const
    {
        int _node = 0;
        std::stack<int> stack;
        int num_stack = 0;
        stack.push(_node);

        while (!stack.empty())
        {
            const node& n = nodes[stack.top()];
            stack.pop();
            if (!any_segment_AABB(_a0, _a1, n.bounds))
                continue;

            if (is_segment(n.left))
            {
                if (any_segment_segment(segments.values.getValue(-n.left), segments.values.getValue(-n.left + 1), _a0, _a1, testLeft, testRight))
                {
                    return true;
                }
            }

            if (is_segment(n.right))
            {
                if (any_segment_segment(segments.values.getValue(-n.right), segments.values.getValue(-n.right + 1), _a0, _a1, testLeft, testRight))
                {
                    return true;
                }
            }

            if (is_node(n.left))
                stack.push(n.left);

            if (is_node(n.right))
                stack.push(n.right);
        }
        return false;
    }

    bool Bvh2d::closestHit(const Eigen::Vector2d& _a0, const Eigen::Vector2d& _a1, int& _line_segment, bool testLeft, bool testRight, Eigen::Vector2d& t) const
    {
        t.y() = std::numeric_limits<double>::max();
        closestHitRecursive(0, _a0, _a1, _line_segment, testLeft, testRight, t);
        return t.y() != std::numeric_limits<double>::max();
    }

    void Bvh2d::closestHitRecursive(int _node, const Eigen::Vector2d& _a0, const Eigen::Vector2d& _a1, int& _line_segment, bool testLeft, bool testRight, Eigen::Vector2d& t) const
    {
        const node& n = nodes[_node];
        if (!any_segment_AABB(_a0, _a1, n.bounds))
            return;

        if (is_segment(n.left))
        {
            Eigen::Vector2d test_t;
            if (intersect_segment_segment(segments.values.getValue(-n.left), segments.values.getValue(-n.left + 1), _a0, _a1, test_t, testLeft, testRight))
            {
                if (test_t.y() < t.y())
                {
                    _line_segment = -n.left;
                    t             = test_t;
                }
            }
        }

        if (is_segment(n.right))
        {
            Eigen::Vector2d test_t;
            if (intersect_segment_segment(segments.values.getValue(-n.right), segments.values.getValue(-n.right + 1), _a0, _a1, test_t, testLeft, testRight))
            {
                if (test_t.y() < t.y())
                {
                    _line_segment = -n.right;
                    t             = test_t;
                }
            }
        }

        if (is_node(n.left))
            closestHitRecursive(n.left, _a0, _a1, _line_segment, testLeft, testRight, t);

        if (is_node(n.right))
            closestHitRecursive(n.right, _a0, _a1, _line_segment, testLeft, testRight, t);
    }

    const Bvh2d::node& Bvh2d::get_node(int _index) const
    {
        assert(is_node(_index));
        return nodes[_index];
    }

    bool Bvh2d::is_segment(int _index) const
    {
        return _index <= 0;
    }

    bool Bvh2d::is_node(int _index) const
    {
        return _index > 0 && _index != node::invalid;
    }

    bool Bvh2d::is_invalid(int _index) const
    {
        return _index == node::invalid;
    }

    void Bvh2d::recompute_bounds(node& _node)
    {
        _node.bounds.setEmpty();
        if (is_segment(_node.left))
        {
            _node.bounds.extend(segments.values.getValue(-_node.left));
            _node.bounds.extend(segments.values.getValue(-_node.left + 1));
        }
        if (is_segment(_node.right))
        {
            _node.bounds.extend(segments.values.getValue(-_node.right));
            _node.bounds.extend(segments.values.getValue(-_node.right + 1));
        }
        if (is_node(_node.left))
        {
            _node.bounds.extend(get_node(_node.left).bounds);
        }
        if (is_node(_node.right))
        {
            _node.bounds.extend(get_node(_node.right).bounds);
        }
    }

    Eigen::Vector2d Bvh2d::get_segment_center(int index) const
    {
        return (segments.values.getValue(index) + segments.values.getValue(index + 1)) * 0.5;
    }

    void Bvh2d::sort_segments(std::vector<int>& _sorted_segments, int _left, int _right, int _axis)
    {
        if (_right <= _left)
            return;

        double comparison_center = get_segment_center(_sorted_segments[_left])[_axis];
        int i                    = _left + 1;
        int j                    = _right;

        // special case of all coordinates being identical -> pick the middle
        if (get_segment_center(_sorted_segments[i])[_axis] == get_segment_center(_sorted_segments[j])[_axis])
        {
            i = (_left + _right) / 2;
            j = (_left + _right) / 2;
        }
        else
        {
            do
            {
                while ((i < _right) && (get_segment_center(_sorted_segments[i])[_axis] <= comparison_center))
                    i++;

                while (get_segment_center(_sorted_segments[j])[_axis] > comparison_center)
                    j--;

                if (i < j)
                    std::swap(_sorted_segments[i], _sorted_segments[j]);
            } while (i < j);
        }

        if (j != _left)
        {
            std::swap(_sorted_segments[j], _sorted_segments[_left]);
        }

        sort_segments(_sorted_segments, _left, j - 1, _axis);
        sort_segments(_sorted_segments, j + 1, _right, _axis);
    }

    void Bvh2d::create_hierarchy(int _node, const Eigen::AlignedBox2d& _bounds, std::vector<int> _sorted_segments_list[2], int _num_segments)
    {
        int split_axis = 0, split_index = 0, i, left_index = 0, right_index = 0;
        Eigen::AlignedBox2d left_box, right_box;
        std::vector<int> left_sorted_segments_list[2];
        std::vector<int> right_sorted_segments_list[2];

        if (_num_segments < 1)
        {
            return;
        }
        if (_num_segments == 1)
        {
            nodes[_node].left = -_sorted_segments_list[0][0];
            recompute_bounds(nodes[_node]);
            return;
        }
        if (_num_segments == 2)
        {
            nodes[_node].left  = -_sorted_segments_list[0][0];
            nodes[_node].right = -_sorted_segments_list[0][1];
            recompute_bounds(nodes[_node]);
            return;
        }

        std::vector<Eigen::AlignedBox2d> left_boxes(_num_segments - 1);
        std::vector<Eigen::AlignedBox2d> right_boxes(_num_segments - 1);

        double total_area = compute_area(_bounds);
        double min_cost   = std::numeric_limits<double>::max();

        for (int axis = 0; axis < 2; axis++)
        {
            create_left_boxes(left_boxes, _num_segments, _sorted_segments_list[axis]);
            create_right_boxes(right_boxes, _num_segments, _sorted_segments_list[axis]);

            for (i = 0; i <= _num_segments - 2; i++)
            {
                double left_area  = compute_area(left_boxes[i]);
                double right_area = compute_area(right_boxes[i]);
                double split_cost = (left_area * (i + 1) + right_area * (_num_segments - i - 1)) / total_area;
                if (split_cost < min_cost)
                {
                    min_cost    = split_cost;
                    split_index = i;
                    split_axis  = axis;
                    left_box    = left_boxes[i];
                    right_box   = right_boxes[i];
                }
            }
        }

        left_sorted_segments_list[split_axis].resize(_num_segments);
        right_sorted_segments_list[split_axis].resize(_num_segments);

        for (i = 0; i <= split_index; i++)
        {
            is_left[_sorted_segments_list[split_axis][i]] = 1;
            left_sorted_segments_list[split_axis][i]      = _sorted_segments_list[split_axis][i];
        }
        for (i = split_index + 1; i <= _num_segments - 1; i++)
        {
            is_left[_sorted_segments_list[split_axis][i]]               = 0;
            right_sorted_segments_list[split_axis][i - split_index - 1] = _sorted_segments_list[split_axis][i];
        }

        for (int axis = 0; axis < 2; axis++)
        {
            if (axis != split_axis)
            {
                left_sorted_segments_list[axis].resize(_num_segments);
                right_sorted_segments_list[axis].resize(_num_segments);

                left_index  = 0;
                right_index = 0;
                for (i = 0; i < _num_segments; i++)
                {
                    if (is_left[_sorted_segments_list[axis][i]] == 1)
                    {
                        left_sorted_segments_list[axis][left_index] = _sorted_segments_list[axis][i];
                        left_index++;
                    }
                    else
                    {
                        right_sorted_segments_list[axis][right_index] = _sorted_segments_list[axis][i];
                        right_index++;
                    }
                }
            }
        }
        assert(left_index != 0);
        assert(right_index != 0);

        if (left_index > 1)
        {
            // allocate a new left node
            nodes.push_back(node());
            int left_node = (int)nodes.size() - 1;
            // set link from parent
            nodes[_node].left = left_node;
            // recursively compute for children
            create_hierarchy(left_node, left_box, left_sorted_segments_list, split_index + 1);
            // compute bounds for left node
            recompute_bounds(nodes[_node]);
        }
        else
        {
            nodes[_node].left = -left_sorted_segments_list[split_axis][0];
            recompute_bounds(nodes[_node]);
        }

        if (right_index > 1)
        {
            // allocate a new right node
            nodes.push_back(node());
            int right_node = (int)nodes.size() - 1;
            // set link from parent
            nodes[_node].right = right_node;
            create_hierarchy(right_node, right_box, right_sorted_segments_list, _num_segments - split_index - 1);
            // compute bounds for right node
            recompute_bounds(nodes[_node]);
        }
        else
        {
            nodes[_node].right = -right_sorted_segments_list[split_axis][0];
            recompute_bounds(nodes[_node]);
        }
    }

    void Bvh2d::create_left_boxes(std::vector<Eigen::AlignedBox2d>& _left_boxes, int _num_segments, const std::vector<int>& _sorted_segments)
    {
        Eigen::AlignedBox2d box;
        box.setEmpty();
        for (int i = 0; i < _num_segments - 1; i++)
        {
            box.extend(segments.values.getValue(_sorted_segments[i]));
            box.extend(segments.values.getValue(_sorted_segments[i] + 1));
            _left_boxes[i] = box;
        }
    }

    void Bvh2d::create_right_boxes(std::vector<Eigen::AlignedBox2d>& _right_boxes, int _num_segments, const std::vector<int>& _sorted_segments)
    {
        Eigen::AlignedBox2d box;
        box.setEmpty();
        for (int i = _num_segments - 1; i > 0; i--)
        {
            box.extend(segments.values.getValue(_sorted_segments[i]));
            box.extend(segments.values.getValue(_sorted_segments[i] + 1));
            _right_boxes[i - 1] = box;
        }
    }

    double Bvh2d::compute_area(const Eigen::AlignedBox2d& _bounds) const
    {
        if (_bounds.isEmpty())
            return std::numeric_limits<double>::max() / 10;
        Eigen::Vector2d delta = _bounds.max() - _bounds.min();
        return 2.0f * (delta.x() + delta.y());
    }
}
