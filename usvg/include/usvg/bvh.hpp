#pragma once

#include "piecewise_linear_curve.hpp"

#include <algorithm>
#include <assert.h>

#include <Eigen/Eigen>

#include <stack>

namespace usvg
{
    /**
     * @brief Class that forms a bounding volume hierarchy over line segments.
     */
    class Bvh2d
    {
    public:
        /**
         * @brief Builds the tree from the given line segments.
         * @param _segments Line segments to build tree over.
         */
        void build(const PiecewiseLinearCurve2d& _segments);

        /**
         * @brief Node in the bvh tree.
         */
        struct node
        {
            static const int invalid = 0xFFFFFFF; // 7*F on purpose!

            /**
             * @brief Constructor.
             */
            node();

            /**
             * @brief Bounding box of this node.
             */
            Eigen::AlignedBox2d bounds;

            /**
             * @brief Index on the left. If positive, this is the index of the child. If negtive, this is the negated index of the line segment.
             */
            int left;

            /**
             * @brief Index on the right. If positive, this is the index of the child. If negtive, this is the negated index of the line segment.
             */
            int right;
        };

        // computes the intersection of two line segments. t contains the relative positions on the two segments.
        static bool intersect_segment_segment(const Eigen::Vector2d& _a, const Eigen::Vector2d& _b, const Eigen::Vector2d& _c, const Eigen::Vector2d& _d, Eigen::Vector2d& _t, bool testLeft, bool testRight);

        // tests if two line segments have an intersection
        static bool any_segment_segment(const Eigen::Vector2d& _b0, const Eigen::Vector2d& _b1, const Eigen::Vector2d& _a0, const Eigen::Vector2d& _a1, bool testLeft, bool testRight);

        // tests if a segment intersects with an axis-aligned bounding box
        static bool any_segment_AABB(const Eigen::Vector2d& _b0, const Eigen::Vector2d& _b1, const Eigen::AlignedBox2d& bounds);

        // returns the closest point to x on a segment with endpoints b0 and b1
        static void closest_point_on_segment(const Eigen::Vector2d& _b0, const Eigen::Vector2d& _b1, const Eigen::Vector2d& _p, int _curve, double& _min_dist, double& _min_t, int& _min_curve);

        // test if the line segment a0-a1 intersects any edge
        [[nodiscard]] bool anyHit(const Eigen::Vector2d& _a0, const Eigen::Vector2d& _a1, bool testLeft, bool testRight) const;

        // Finds the edge that has the closest hit with the line segment a0-a1.
        [[nodiscard]] bool closestHit(const Eigen::Vector2d& _a0, const Eigen::Vector2d& _a1, int& _line_segment, bool testLeft, bool testRight, Eigen::Vector2d& t) const;

        /**
         * @brief Nodes of the bvh. the first node is the root.
         */
        std::vector<node> nodes;

    private:
        // Finds the edge that has the closest hit with the line segment a0-a1.
        [[nodiscard]] void closestHitRecursive(int _node, const Eigen::Vector2d& _a0, const Eigen::Vector2d& _a1, int& _line_segment, bool testLeft, bool testRight, Eigen::Vector2d& t) const;

        /**
         * @brief Segments for which the tree is built.
         */
        PiecewiseLinearCurve2d segments;

        /**
         * @brief Gets a node by linear index.
         * @param _index Index of node to get. This is a number > 0.
         * @return BVH tree node.
         */
        const node& get_node(int _index) const;

        /**
         * @brief Tests if a given index belongs to a line segments, i.e., _index <= 0.
         * @param _index Index to check.
         * @return True if this is a line segment.
         */
        bool is_segment(int _index) const;

        /**
         * @brief Tests if a given index belongs to a BVH tree node, i.e., index > 0.
         * @param _index Index to check.
         * @return True if this is a node.
         */
        bool is_node(int _index) const;

        /**
         * @brief Checks if this is an invalid node.
         * @param _index Index to check.
         * @return True if invalid.
         */
        bool is_invalid(int _index) const;

        /**
         * @brief Recomputes the bounding box of a given node.
         * @param _node Node to compute the bounding box for.
         */
        void recompute_bounds(node& _node);

    private:
        /**
         * @brief Computes the center of a given line segment.
         * @param index Index of line segment.
         * @return Position of the center.
         */
        Eigen::Vector2d get_segment_center(int index) const;

        /**
         * @brief Sorts the segments.
         * @param _sorted_segments Index buffer that holds the sorting.
         * @param _left Starting index.
         * @param _right Ending index.
         * @param _axis Axis over which to sort.
         */
        void sort_segments(std::vector<int>& _sorted_segments, int _left, int _right, int _axis);

        /**
         * @brief Creates the bounding volume hierarchy for the nodes under _node.
         * @param _node Index of node to compute the hierarchy of.
         * @param _bounds Bounding box of the data inside this node.
         * @param _sorted_segments_list List of all the sorted segments.
         * @param _num_segments Number of segments.
         */
        void create_hierarchy(int _node, const Eigen::AlignedBox2d& _bounds, std::vector<int> _sorted_segments_list[2], int _num_segments);

        /**
         * @brief Creates the bounding boxes on the left.
         * @param _left_boxes Pointer where to store the bounding boxes.
         * @param _num_segments Number of segments.
         * @param _sorted_segments Reference to the sorted segments.
         */
        void create_left_boxes(std::vector<Eigen::AlignedBox2d>& _left_boxes, int _num_segments, const std::vector<int>& _sorted_segments);

        /**
         * @brief Creates the bounding boxes on the right.
         * @param _right_boxes Pointer where to store the bounding boxes.
         * @param _num_segments Number of segments.
         * @param _sorted_segments Reference to the sorted segments.
         */
        void create_right_boxes(std::vector<Eigen::AlignedBox2d>& _right_boxes, int _num_segments, const std::vector<int>& _sorted_segments);

        /**
         * @brief Computes the circumference of a 2D bounding box.
         * @param _bounds Bounding box to compute circumference for.
         * @return Circumference of bounding box.
         */
        [[nodiscard]] double compute_area(const Eigen::AlignedBox2d& _bounds) const;

        /**
         * @brief Index array of sorted objects.
         */
        std::vector<int> sorted_segments[2];

        /**
         * @brief Flag that determines whether an element is left or right of the center.
         */
        std::vector<int> is_left;
    };
}
