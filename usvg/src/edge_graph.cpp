#include <usvg/edge_graph.hpp>

#include <usvg/bezier_discretizer.hpp>
#include <usvg/diffusion_curve.hpp>
#include <usvg/gradient_mesh.hpp>

namespace usvg
{
    bool EdgeGraph::InputBoundaryCurve::intersect(const InputBoundaryCurve& _other, std::vector<Eigen::Vector2d>& _intersections) const
    {
        for (int ivertex = 0; ivertex < _other.position.values.getSize() - 1; ++ivertex)
        {
            Eigen::Vector2d q0 = _other.position.values.getValue(ivertex);
            Eigen::Vector2d q1 = _other.position.values.getValue(ivertex + 1);
            int isegment       = 0;
            Eigen::Vector2d t;
            if (bvh.closestHit(q0, q1, isegment, true, true, t))
            {
                Eigen::Vector2d p0 = position.values.getValue(isegment);
                Eigen::Vector2d p1 = position.values.getValue(isegment + 1);
                double param_p     = (1 - t.x()) * position.parameters.getValue(isegment).x() + t.x() * position.parameters.getValue(isegment + 1).x();
                double param_q     = (1 - t.y()) * _other.position.parameters.getValue(ivertex).x() + t.y() * _other.position.parameters.getValue(ivertex + 1).x();
                _intersections.push_back(Eigen::Vector2d(param_p, param_q));
            }
        }
        return !_intersections.empty();
    }

    std::pair<std::shared_ptr<EdgeGraph::InputBoundaryCurve>, std::shared_ptr<EdgeGraph::InputBoundaryCurve>> EdgeGraph::InputBoundaryCurve::subdivide(double _t) const
    {
        auto partA = std::make_shared<EdgeGraph::InputBoundaryCurve>();
        auto partB = std::make_shared<EdgeGraph::InputBoundaryCurve>();
        this->position.subdivide(partA->position, partB->position, _t);
        this->colorLeft.subdivide(partA->colorLeft, partB->colorLeft, _t);
        this->colorRight.subdivide(partA->colorRight, partB->colorRight, _t);
        partA->boundaryConditionLeft  = this->boundaryConditionLeft;
        partB->boundaryConditionLeft  = this->boundaryConditionLeft;
        partA->boundaryConditionRight = this->boundaryConditionRight;
        partB->boundaryConditionRight = this->boundaryConditionRight;
        partA->gradientMesh           = this->gradientMesh;
        partB->gradientMesh           = this->gradientMesh;
        partA->update();
        partB->update();
        return std::make_pair(partA, partB);
    }

    std::tuple<Eigen::Vector2d, double, double> EdgeGraph::InputBoundaryCurve::closestPoint(const Eigen::Vector2d& _point)
    {
        double min_dist = std::numeric_limits<double>::max();
        double min_t    = 0;
        int min_segment = 0;
        for (int ivertex = 0; ivertex < position.values.getSize() - 1; ++ivertex)
        {
            Bvh2d::closest_point_on_segment(position.values.getValue(ivertex), position.values.getValue(ivertex + 1), _point, ivertex, min_dist, min_t, min_segment);
        }
        double param              = (1 - min_t) * position.parameters.getValue(min_segment).x() + min_t * position.parameters.getValue(min_segment + 1).x();
        Eigen::Vector2d min_point = (1 - min_t) * position.values.getValue(min_segment) + min_t * position.values.getValue(min_segment + 1);
        return std::make_tuple(min_point, min_dist, param);
    }

    void EdgeGraph::InputBoundaryCurve::update()
    {
        bvh.build(position);
    }

    EdgeGraph::Vertex::Vertex(const Eigen::Vector2d& _position)
        : position(_position)
    {
    }

    std::shared_ptr<EdgeGraph::Edge> EdgeGraph::Vertex::nextEdge(std::shared_ptr<const Edge> _incoming)
    {
        // if there is only one edge, return it regardless of whether we are coming from it
        if (edges.size() == 1)
            return edges.front().lock();

        // travers all edges
        size_t numEdges = edges.size();
        for (auto itedge = edges.begin(); itedge != edges.end(); ++itedge)
        {
            // find the incoming edge
            auto edge = *itedge;
            if (edge.lock() == _incoming)
            {
                // the next edge is turning right
                itedge++;
                if (itedge == edges.end())
                    // if we had the last edge, then the next one is the first edge of the list
                    return edges.begin()->lock();
                else
                    return itedge->lock();
            }
        }
        return nullptr;
    }

    EdgeGraph::Edge::Edge(std::shared_ptr<InputBoundaryCurve> _boundaryCurve)
        : boundaryCurve(_boundaryCurve)
    {
    }

    bool EdgeGraph::Edge::intersect(const Edge& _other, std::vector<Eigen::Vector2d>& _intersections) const
    {
        return boundaryCurve->intersect(*_other.boundaryCurve, _intersections);
    }

    std::pair<std::shared_ptr<EdgeGraph::Edge>, std::shared_ptr<EdgeGraph::Vertex>> EdgeGraph::Edge::subdivide(double _t)
    {
        // subdivide curve into two pieces
        auto [ibcA, ibcB] = this->boundaryCurve->subdivide(_t);

        // allocate a new vertex
        auto vert = std::make_shared<EdgeGraph::Vertex>(ibcA->position.values.last());

        // make a new curve for the second half
        auto partB           = std::make_shared<EdgeGraph::Edge>(ibcB);
        partB->boundaryCurve = ibcB;
        partB->vertexA       = vert;
        partB->vertexB       = this->vertexB.lock();

        // update the incident edge on old vertex B: remove old edge
        std::shared_ptr<Edge> partA = shared_from_this();

        // if cyclic, don't take out the connection
        if (this->vertexA.lock() != this->vertexB.lock())
        {
            partB->vertexB.lock()->edges.remove_if([partA](std::weak_ptr<EdgeGraph::Edge> e)
                                                   { return !(partA.owner_before(e) || e.owner_before(partA)); });
        }
        partB->vertexB.lock()->edges.push_back(partB); // add the new edge

        // update this curve to the first half
        this->boundaryCurve = ibcA;
        this->vertexB       = vert;

        // connect the new vertex to its two edges
        vert->edges.push_back(partA);
        vert->edges.push_back(partB);

        // return the new two objects for insertion
        return std::make_pair(partB, vert);
    }

    EdgeGraph::EdgeGraph(const std::vector<std::shared_ptr<DiffusionCurve>>& _diffusionCurves,
                         const std::vector<std::shared_ptr<GradientMesh>>& _gradientMeshes,
                         double _vertexMergeThreshold, double _epsilon, double _distanceThreshold)
    {
        assembleInputBoundaryCurves(_diffusionCurves, _gradientMeshes, _distanceThreshold);
        constructEdgeGraph(_vertexMergeThreshold, _epsilon);
    }

    std::shared_ptr<EdgeGraph::InputBoundaryCurve> EdgeGraph::makeInputBoundaryCurveFromDiffusionCurve(std::shared_ptr<const DiffusionCurve> _curve, double _distanceThreshold)
    {
        auto result          = std::make_shared<EdgeGraph::InputBoundaryCurve>();
        result->gradientMesh = nullptr;
        BezierDiscretizer2d::topDown(_curve->position, _distanceThreshold, result->position);
        result->boundaryConditionLeft  = _curve->boundaryConditionLeft;
        result->boundaryConditionRight = _curve->boundaryConditionRight;
        result->colorLeft.domain       = _curve->colorLeft.domain;
        result->colorRight.domain      = _curve->colorRight.domain;
        result->colorLeft.controlPoints.setSize((_curve->colorLeft.values.getSize() - 1) * 3 + 1);
        result->colorRight.controlPoints.setSize((_curve->colorRight.values.getSize() - 1) * 3 + 1);
        for (Eigen::Index ivert = 0; ivert < _curve->colorLeft.values.getSize() - 1; ++ivert)
        {
            const Eigen::Vector3d& c0 = _curve->colorLeft.values.getValue(ivert);
            const Eigen::Vector3d& c1 = _curve->colorLeft.values.getValue(ivert + 1);
            result->colorLeft.controlPoints.setValue(ivert * 3 + 0, c0);
            result->colorLeft.controlPoints.setValue(ivert * 3 + 1, c0 + (c1 - c0) * 1. / 3.);
            result->colorLeft.controlPoints.setValue(ivert * 3 + 2, c0 + (c1 - c0) * 2. / 3.);
            result->colorLeft.controlPoints.setValue(ivert * 3 + 3, c1);
        }
        for (Eigen::Index ivert = 0; ivert < _curve->colorRight.values.getSize() - 1; ++ivert)
        {
            const Eigen::Vector3d& c0 = _curve->colorRight.values.getValue(ivert);
            const Eigen::Vector3d& c1 = _curve->colorRight.values.getValue(ivert + 1);
            result->colorRight.controlPoints.setValue(ivert * 3 + 0, c0);
            result->colorRight.controlPoints.setValue(ivert * 3 + 1, c0 + (c1 - c0) * 1. / 3.);
            result->colorRight.controlPoints.setValue(ivert * 3 + 2, c0 + (c1 - c0) * 2. / 3.);
            result->colorRight.controlPoints.setValue(ivert * 3 + 3, c1);
        }
        result->position.recomputeBoundingBox();
        result->colorLeft.recomputeBoundingBox();
        result->colorRight.recomputeBoundingBox();
        result->colorLeft.parameters  = _curve->colorLeft.parameters;
        result->colorRight.parameters = _curve->colorRight.parameters;
        result->position.setExplicitParameterization();
        result->colorLeft.setExplicitParameterization();
        result->colorRight.setExplicitParameterization();

        bool isValid = result->position.isValid();
        isValid      = result->colorLeft.isValid();
        isValid      = result->colorRight.isValid();

        result->update();

        return result;
    }

    std::tuple<std::shared_ptr<EdgeGraph::InputBoundaryCurve>, std::shared_ptr<EdgeGraph::InputBoundaryCurve>, std::shared_ptr<EdgeGraph::InputBoundaryCurve>, std::shared_ptr<EdgeGraph::InputBoundaryCurve>> EdgeGraph::makeInputBoundaryCurveFromGradientMesh(std::shared_ptr<GradientMesh> _gradientMesh, double _distanceThreshold)
    {
        auto c1                    = std::make_shared<EdgeGraph::InputBoundaryCurve>();
        auto c2                    = std::make_shared<EdgeGraph::InputBoundaryCurve>();
        auto c3                    = std::make_shared<EdgeGraph::InputBoundaryCurve>();
        auto c4                    = std::make_shared<EdgeGraph::InputBoundaryCurve>();
        c1->gradientMesh           = _gradientMesh;
        c2->gradientMesh           = _gradientMesh;
        c3->gradientMesh           = _gradientMesh;
        c4->gradientMesh           = _gradientMesh;
        c1->boundaryConditionLeft  = EBoundaryCondition::Neumann;
        c2->boundaryConditionLeft  = EBoundaryCondition::Neumann;
        c3->boundaryConditionLeft  = EBoundaryCondition::Neumann;
        c4->boundaryConditionLeft  = EBoundaryCondition::Neumann;
        c1->boundaryConditionRight = EBoundaryCondition::Dirichlet;
        c2->boundaryConditionRight = EBoundaryCondition::Dirichlet;
        c3->boundaryConditionRight = EBoundaryCondition::Dirichlet;
        c4->boundaryConditionRight = EBoundaryCondition::Dirichlet;
        CubicBezierSpline2d c1_position, c2_position, c3_position, c4_position;
        c1_position.controlPoints.setSize((_gradientMesh->resolution.x()) * 3 + 1);    // N
        c2_position.controlPoints.setSize((_gradientMesh->resolution.y()) * 3 + 1);    // E
        c3_position.controlPoints.setSize((_gradientMesh->resolution.x()) * 3 + 1);    // S
        c4_position.controlPoints.setSize((_gradientMesh->resolution.y()) * 3 + 1);    // W
        c1->colorRight.controlPoints.setSize((_gradientMesh->resolution.x()) * 3 + 1); // N
        c2->colorRight.controlPoints.setSize((_gradientMesh->resolution.y()) * 3 + 1); // E
        c3->colorRight.controlPoints.setSize((_gradientMesh->resolution.x()) * 3 + 1); // S
        c4->colorRight.controlPoints.setSize((_gradientMesh->resolution.y()) * 3 + 1); // W
        c1->colorLeft.controlPoints.setSize(4);
        c2->colorLeft.controlPoints.setSize(4);
        c3->colorLeft.controlPoints.setSize(4);
        c4->colorLeft.controlPoints.setSize(4);
        for (int i = 0; i < 4; ++i)
        {
            c1->colorLeft.controlPoints.setValue(i, Eigen::Vector3d::Zero());
            c2->colorLeft.controlPoints.setValue(i, Eigen::Vector3d::Zero());
            c3->colorLeft.controlPoints.setValue(i, Eigen::Vector3d::Zero());
            c4->colorLeft.controlPoints.setValue(i, Eigen::Vector3d::Zero());
        }
        c1->colorLeft.setUniformParameterization();
        c2->colorLeft.setUniformParameterization();
        c3->colorLeft.setUniformParameterization();
        c4->colorLeft.setUniformParameterization();
        for (Eigen::Index ix = 0; ix < _gradientMesh->resolution.x(); ++ix)
        {
            CubicBezierCurve2d posCurve;
            CubicBezierCurve3d colorCurve;
            // N
            _gradientMesh->getTile(ix, _gradientMesh->resolution.y() - 1).position.getBezierSurface().getBoundaryCurveNorth(posCurve);
            _gradientMesh->getTile(ix, _gradientMesh->resolution.y() - 1).color.getBezierSurface().getBoundaryCurveNorth(colorCurve);
            c1_position.controlPoints.setValue(ix * 3 + 0, posCurve.controlPoints[0]);
            c1_position.controlPoints.setValue(ix * 3 + 1, posCurve.controlPoints[1]);
            c1_position.controlPoints.setValue(ix * 3 + 2, posCurve.controlPoints[2]);
            c1_position.controlPoints.setValue(ix * 3 + 3, posCurve.controlPoints[3]);
            c1->colorRight.controlPoints.setValue(ix * 3 + 0, colorCurve.controlPoints[0]);
            c1->colorRight.controlPoints.setValue(ix * 3 + 1, colorCurve.controlPoints[1]);
            c1->colorRight.controlPoints.setValue(ix * 3 + 2, colorCurve.controlPoints[2]);
            c1->colorRight.controlPoints.setValue(ix * 3 + 3, colorCurve.controlPoints[3]);

            // S
            _gradientMesh->getTile(_gradientMesh->resolution.x() - 1 - ix, 0).position.getBezierSurface().getBoundaryCurveSouth(posCurve);
            _gradientMesh->getTile(_gradientMesh->resolution.x() - 1 - ix, 0).color.getBezierSurface().getBoundaryCurveSouth(colorCurve);
            c3_position.controlPoints.setValue(ix * 3 + 0, posCurve.controlPoints[0]);
            c3_position.controlPoints.setValue(ix * 3 + 1, posCurve.controlPoints[1]);
            c3_position.controlPoints.setValue(ix * 3 + 2, posCurve.controlPoints[2]);
            c3_position.controlPoints.setValue(ix * 3 + 3, posCurve.controlPoints[3]);
            c3->colorRight.controlPoints.setValue(ix * 3 + 0, colorCurve.controlPoints[0]);
            c3->colorRight.controlPoints.setValue(ix * 3 + 1, colorCurve.controlPoints[1]);
            c3->colorRight.controlPoints.setValue(ix * 3 + 2, colorCurve.controlPoints[2]);
            c3->colorRight.controlPoints.setValue(ix * 3 + 3, colorCurve.controlPoints[3]);
        }
        for (Eigen::Index iy = 0; iy < _gradientMesh->resolution.y(); ++iy)
        {
            CubicBezierCurve2d posCurve;
            CubicBezierCurve3d colorCurve;
            // E
            _gradientMesh->getTile(_gradientMesh->resolution.x() - 1, _gradientMesh->resolution.y() - 1 - iy).position.getBezierSurface().getBoundaryCurveEast(posCurve);
            _gradientMesh->getTile(_gradientMesh->resolution.x() - 1, _gradientMesh->resolution.y() - 1 - iy).color.getBezierSurface().getBoundaryCurveEast(colorCurve);
            c2_position.controlPoints.setValue(iy * 3 + 0, posCurve.controlPoints[0]);
            c2_position.controlPoints.setValue(iy * 3 + 1, posCurve.controlPoints[1]);
            c2_position.controlPoints.setValue(iy * 3 + 2, posCurve.controlPoints[2]);
            c2_position.controlPoints.setValue(iy * 3 + 3, posCurve.controlPoints[3]);
            c2->colorRight.controlPoints.setValue(iy * 3 + 0, colorCurve.controlPoints[0]);
            c2->colorRight.controlPoints.setValue(iy * 3 + 1, colorCurve.controlPoints[1]);
            c2->colorRight.controlPoints.setValue(iy * 3 + 2, colorCurve.controlPoints[2]);
            c2->colorRight.controlPoints.setValue(iy * 3 + 3, colorCurve.controlPoints[3]);

            // W
            _gradientMesh->getTile(0, iy).position.getBezierSurface().getBoundaryCurveWest(posCurve);
            _gradientMesh->getTile(0, iy).color.getBezierSurface().getBoundaryCurveWest(colorCurve);
            c4_position.controlPoints.setValue(iy * 3 + 0, posCurve.controlPoints[0]);
            c4_position.controlPoints.setValue(iy * 3 + 1, posCurve.controlPoints[1]);
            c4_position.controlPoints.setValue(iy * 3 + 2, posCurve.controlPoints[2]);
            c4_position.controlPoints.setValue(iy * 3 + 3, posCurve.controlPoints[3]);
            c4->colorRight.controlPoints.setValue(iy * 3 + 0, colorCurve.controlPoints[0]);
            c4->colorRight.controlPoints.setValue(iy * 3 + 1, colorCurve.controlPoints[1]);
            c4->colorRight.controlPoints.setValue(iy * 3 + 2, colorCurve.controlPoints[2]);
            c4->colorRight.controlPoints.setValue(iy * 3 + 3, colorCurve.controlPoints[3]);
        }
        for (Eigen::Index ix = 0; ix <= _gradientMesh->resolution.x(); ++ix)
        {
            c1_position.parameters.append(Eigen::Vector1d(ix / double(_gradientMesh->resolution.x())));
            c1->colorRight.parameters.append(Eigen::Vector1d(ix / double(_gradientMesh->resolution.x())));
            c3_position.parameters.append(Eigen::Vector1d(ix / double(_gradientMesh->resolution.x())));
            c3->colorRight.parameters.append(Eigen::Vector1d(ix / double(_gradientMesh->resolution.x())));
        }
        for (Eigen::Index iy = 0; iy <= _gradientMesh->resolution.y(); ++iy)
        {
            c2_position.parameters.append(Eigen::Vector1d(iy / double(_gradientMesh->resolution.y())));
            c2->colorRight.parameters.append(Eigen::Vector1d(iy / double(_gradientMesh->resolution.y())));
            c4_position.parameters.append(Eigen::Vector1d(iy / double(_gradientMesh->resolution.y())));
            c4->colorRight.parameters.append(Eigen::Vector1d(iy / double(_gradientMesh->resolution.y())));
        }
        c1_position.setExplicitParameterization();
        c2_position.setExplicitParameterization();
        c3_position.setExplicitParameterization();
        c4_position.setExplicitParameterization();
        c1->colorRight.setExplicitParameterization();
        c2->colorRight.setExplicitParameterization();
        c3->colorRight.setExplicitParameterization();
        c4->colorRight.setExplicitParameterization();
        c1_position.recomputeBoundingBox();
        c2_position.recomputeBoundingBox();
        c3_position.recomputeBoundingBox();
        c4_position.recomputeBoundingBox();
        c1->colorRight.recomputeBoundingBox();
        c2->colorRight.recomputeBoundingBox();
        c3->colorRight.recomputeBoundingBox();
        c4->colorRight.recomputeBoundingBox();
        c1->colorLeft.recomputeBoundingBox();
        c2->colorLeft.recomputeBoundingBox();
        c3->colorLeft.recomputeBoundingBox();
        c4->colorLeft.recomputeBoundingBox();

        BezierDiscretizer2d::topDown(c1_position, _distanceThreshold, c1->position);
        BezierDiscretizer2d::topDown(c2_position, _distanceThreshold, c2->position);
        BezierDiscretizer2d::topDown(c3_position, _distanceThreshold, c3->position);
        BezierDiscretizer2d::topDown(c4_position, _distanceThreshold, c4->position);

        c1->update();
        c2->update();
        c3->update();
        c4->update();

        return std::make_tuple(c1, c2, c3, c4);
    }

    void EdgeGraph::assembleInputBoundaryCurves(const std::vector<std::shared_ptr<DiffusionCurve>>& _diffusionCurves,
                                                const std::vector<std::shared_ptr<GradientMesh>>& _gradientMeshes, double _distanceThreshold)
    {
        // delete old ones
        inputBoundaryCurves.clear();

        // make 4 curves per gradient mesh
        for (auto& gradientMesh : _gradientMeshes)
        {
            auto [c1, c2, c3, c4] = makeInputBoundaryCurveFromGradientMesh(gradientMesh, _distanceThreshold);
            inputBoundaryCurves.push_back(c1);
            inputBoundaryCurves.push_back(c2);
            inputBoundaryCurves.push_back(c3);
            inputBoundaryCurves.push_back(c4);
        }

        // make 1 curve per diffusion curve, unless it is a closed curve
        for (auto& diffusionCurve : _diffusionCurves)
        {
            if ((diffusionCurve->position.controlPoints.first() - diffusionCurve->position.controlPoints.last()).stableNorm() < 1E-10)
            {
                CubicBezierSpline2d position1, position2;
                diffusionCurve->position.subdivide(position1, position2, 0.5);
                PiecewiseLinearCurve3d colorLeft1, colorLeft2, colorRight1, colorRight2;
                diffusionCurve->colorLeft.subdivide(colorLeft1, colorLeft2, 0.5);
                diffusionCurve->colorRight.subdivide(colorRight1, colorRight2, 0.5);

                auto diffCurve1                    = std::make_shared<DiffusionCurve>();
                diffCurve1->position               = position1;
                diffCurve1->colorLeft              = colorLeft1;
                diffCurve1->colorRight             = colorRight1;
                diffCurve1->boundaryConditionLeft  = diffusionCurve->boundaryConditionLeft;
                diffCurve1->boundaryConditionRight = diffusionCurve->boundaryConditionRight;

                auto diffCurve2                    = std::make_shared<DiffusionCurve>();
                diffCurve2->position               = position2;
                diffCurve2->colorLeft              = colorLeft2;
                diffCurve2->colorRight             = colorRight2;
                diffCurve2->boundaryConditionLeft  = diffusionCurve->boundaryConditionLeft;
                diffCurve2->boundaryConditionRight = diffusionCurve->boundaryConditionRight;

                auto incurve1 = makeInputBoundaryCurveFromDiffusionCurve(diffCurve1, _distanceThreshold);
                auto incurve2 = makeInputBoundaryCurveFromDiffusionCurve(diffCurve2, _distanceThreshold);
                inputBoundaryCurves.push_back(incurve1);
                inputBoundaryCurves.push_back(incurve2);
            }
            else
            {
                auto incurve = makeInputBoundaryCurveFromDiffusionCurve(diffusionCurve, _distanceThreshold);
                inputBoundaryCurves.push_back(incurve);
            }
        }
    }

    void EdgeGraph::constructEdgeGraph(double _vertexMergeThreshold, double _epsilon)
    {
        // begin with empty graph
        vertices.clear();
        edges.clear();

        // insert the input boundary curves one at a time
        for (size_t ibc = 0; ibc < inputBoundaryCurves.size(); ++ibc)
        {
            // if the edge is too short, ignore it
            if (inputBoundaryCurves[ibc]->position.getBoundingBox().diagonal().stableNorm() < _vertexMergeThreshold)
            {
                continue;
            }

            // we have a stack of edges that are currently under testing.
            // those are edges that have recently been added.
            std::stack<std::shared_ptr<Edge>> edgesToTest;

            // allocate an edge for the next input boundary curve
            {
                auto edge = std::make_shared<EdgeGraph::Edge>(inputBoundaryCurves[ibc]);
                if ((inputBoundaryCurves[ibc]->position.values.first() - inputBoundaryCurves[ibc]->position.values.last()).stableNorm() < _vertexMergeThreshold)
                {
                    auto vertexA = std::make_shared<EdgeGraph::Vertex>(inputBoundaryCurves[ibc]->position.values.first());
                    auto vertexB = vertexA;
                    vertexA->edges.push_back(edge);
                    edge->vertexA = vertexA;
                    edge->vertexB = vertexB;
                    // insert vertices and insert edge
                    vertices.push_back(vertexA);
                    edges.push_back(edge);
                    edgesToTest.push(edge);
                }
                else
                {
                    auto vertexA = std::make_shared<EdgeGraph::Vertex>(inputBoundaryCurves[ibc]->position.values.first());
                    auto vertexB = std::make_shared<EdgeGraph::Vertex>(inputBoundaryCurves[ibc]->position.values.last());
                    vertexA->edges.push_back(edge);
                    vertexB->edges.push_back(edge);
                    edge->vertexA = vertexA;
                    edge->vertexB = vertexB;
                    // insert vertices and insert edge
                    vertices.push_back(vertexA);
                    vertices.push_back(vertexB);
                    edges.push_back(edge);
                    edgesToTest.push(edge);
                }
            }

            // while there is another edge to test against all others
            while (!edgesToTest.empty())
            {
                // delete invalid
                deleteInvalid(vertices, edges, _epsilon);

                // pop the next edge
                auto nextEdge = edgesToTest.top();
                edgesToTest.pop();

                // invalid edges and vertices are removed from the collections, but the edgesToTest may contain invalid edges. Identify and skip those.
                if (!nextEdge->vertexA.lock() || !nextEdge->vertexB.lock())
                    continue;

                // iterate all existing vertices and see if we can connect the next edge to an existing vertex
                bool performedVertexOperation = false;
                for (auto otherVertex : vertices)
                {
                    auto vertA = nextEdge->vertexA.lock();
                    if (vertA != otherVertex && (vertA->position - otherVertex->position).stableNorm() < _vertexMergeThreshold)
                    {
                        // redirect edges
                        redirectEdges(nextEdge->vertexA.lock(), otherVertex);
                        // connect next edge to the vertex of the other edge
                        nextEdge->vertexA = otherVertex;
                        // reinsert this edge to start over
                        edgesToTest.push(nextEdge);
                        // mark that something happened with this edge
                        performedVertexOperation = true;
                        break;
                    }
                    auto vertB = nextEdge->vertexB.lock();
                    if (vertB != otherVertex && (vertB->position - otherVertex->position).stableNorm() < _vertexMergeThreshold)
                    {
                        // redirect edges
                        redirectEdges(nextEdge->vertexB.lock(), otherVertex);
                        // connect next edge to the vertex of the other edge
                        nextEdge->vertexB = otherVertex;
                        // reinsert this edge to start over
                        edgesToTest.push(nextEdge);
                        // mark that something happened with this edge
                        performedVertexOperation = true;
                        break;
                    }
                }
                // if something happened to nextEdge, start over.
                if (performedVertexOperation)
                    continue;

                // see if the edge intersects any existing edge
                for (size_t ie = 0; ie < edges.size(); ++ie)
                {
                    // get the other edge to test against
                    auto otherEdge = edges[ie];

                    // avoid comparison against self
                    if (nextEdge == otherEdge)
                        continue;

                    // see if the edge end point (vertexA) is coming close enough to another edge to snap onto it
                    {
                        auto vertA = nextEdge->vertexA.lock();
                        if (vertA != otherEdge->vertexA.lock() && vertA != otherEdge->vertexB.lock())
                        {
                            auto [closest_pnt, closest_dist, closest_t] = otherEdge->boundaryCurve->closestPoint(vertA->position);
                            // if the point is close enough to the curve and not close to an end point
                            if (closest_dist < _vertexMergeThreshold)
                            {
                                auto [otherEdge2, newVert] = otherEdge->subdivide(closest_t);
                                // add the new edge and new vertex to the collection
                                edges.push_back(otherEdge2);
                                vertices.push_back(newVert);
                                // redirect edges
                                auto from = nextEdge->vertexA.lock();
                                redirectEdges(from, newVert);
                                // connect this edge to the new vertex
                                nextEdge->vertexA = newVert;
                                // enqueue everyone who was touched for another round of testing
                                for (auto e : newVert->edges)
                                    edgesToTest.push(e.lock());
                                break;
                            }
                        }
                    }

                    // see if the edge end point (vertexB) is coming close enough to another edge to snap onto it
                    {
                        auto vertB = nextEdge->vertexB.lock();
                        if (vertB != otherEdge->vertexA.lock() && vertB != otherEdge->vertexB.lock())
                        {
                            auto [closest_pnt, closest_dist, closest_t] = otherEdge->boundaryCurve->closestPoint(vertB->position);
                            // if the point is close enough to the curve and not close to an end point
                            if (closest_dist < _vertexMergeThreshold)
                            {
                                auto [otherEdge2, newVert] = otherEdge->subdivide(closest_t);
                                // add the new edge to the collection
                                edges.push_back(otherEdge2);
                                vertices.push_back(newVert);
                                // redirect edges
                                auto from = nextEdge->vertexB.lock();
                                redirectEdges(from, newVert);
                                // connect this edge to the new vertex
                                nextEdge->vertexB = newVert;
                                // enqueue everyone who was touched for another round of testing
                                for (auto e : newVert->edges)
                                    edgesToTest.push(e.lock());
                                break;
                            }
                        }
                    }

                    // see if other line (vertex A) ends on the next edge to test
                    {
                        auto vertA = otherEdge->vertexA.lock();
                        if (vertA != nextEdge->vertexA.lock() && vertA != nextEdge->vertexB.lock())
                        {
                            auto [closest_pnt, closest_dist, closest_t] = nextEdge->boundaryCurve->closestPoint(vertA->position);
                            // if the point is close enough to the curve and not close to an end point
                            if (closest_dist < _vertexMergeThreshold)
                            {
                                auto [nextEdge2, newVert] = nextEdge->subdivide(closest_t);
                                // add the new edge to the collection
                                edges.push_back(nextEdge2);
                                // redirect edges
                                redirectEdges(newVert, vertA);
                                // connect the edges to the other vertex
                                nextEdge->vertexB  = vertA;
                                nextEdge2->vertexA = vertA;
                                // enqueue everyone who was touched for another round of testing
                                for (auto e : vertA->edges)
                                    edgesToTest.push(e.lock());
                                break;
                            }
                        }
                    }

                    // see if other line (vertex B) ends on the next edge to test
                    {
                        auto vertB = otherEdge->vertexB.lock();
                        if (vertB != nextEdge->vertexA.lock() && vertB != nextEdge->vertexB.lock())
                        {
                            auto [closest_pnt, closest_dist, closest_t] = nextEdge->boundaryCurve->closestPoint(vertB->position);
                            // if the point is close enough to the curve and not close to an end point
                            if (closest_dist < _vertexMergeThreshold)
                            {
                                auto [nextEdge2, newVert] = nextEdge->subdivide(closest_t);
                                // add the new edge to the collection
                                edges.push_back(nextEdge2);
                                // redirect edges
                                redirectEdges(newVert, vertB);
                                // connect the edges to the other vertex
                                nextEdge->vertexB  = vertB;
                                nextEdge2->vertexA = vertB;
                                // enqueue everyone who was touched for another round of testing
                                for (auto e : vertB->edges)
                                    edgesToTest.push(e.lock());
                                break;
                            }
                        }
                    }

                    // we find all intersections, but only ever process the first intersection
                    std::vector<Eigen::Vector2d> intersections;
                    if (nextEdge->intersect(*otherEdge, intersections))
                    {
                        for (auto intersection : intersections)
                        {
                            // skip endpoint intersections, since we already handled those
                            if (std::abs(intersection.x() - nextEdge->boundaryCurve->position.domain.min()[0]) < _epsilon ||
                                std::abs(intersection.x() - nextEdge->boundaryCurve->position.domain.max()[0]) < _epsilon ||
                                std::abs(intersection.y() - otherEdge->boundaryCurve->position.domain.min()[0]) < _epsilon ||
                                std::abs(intersection.y() - otherEdge->boundaryCurve->position.domain.max()[0]) < _epsilon)
                                continue;

                            // skip all invalid intersections
                            Eigen::Vector2d t0 = nextEdge->boundaryCurve->position.sample(intersection.x());
                            Eigen::Vector2d t1 = otherEdge->boundaryCurve->position.sample(intersection.y());
                            if ((t0 - t1).stableNorm() > _epsilon)
                                continue;

                            // subdivide the edges, each giving us a new vertex and a new edge
                            auto [nextEdge1, nextVertex]   = nextEdge->subdivide(intersection.x());
                            auto [otherEdge1, otherVertex] = otherEdge->subdivide(intersection.y());

                            // redirect edges
                            redirectEdges(otherVertex, nextVertex);
                            // connect the edges to the other vertex
                            otherEdge->vertexB  = nextVertex;
                            otherEdge1->vertexA = nextVertex;

                            // add the newly created vertices and edges
                            vertices.push_back(nextVertex);
                            edges.push_back(nextEdge1);
                            edges.push_back(otherEdge1);

                            // enqueue everyone who was touched for another round of testing
                            for (auto e : nextVertex->edges)
                                edgesToTest.push(e.lock());

                            // we only process one valid intersection per iteration
                            break;
                        }
                    }
                }
            }

            // delete invalid
            deleteInvalid(vertices, edges, _epsilon);
        }

        // clean-up collapsed bezier curves
        for (auto edge : edges)
        {
            auto& cps = edge->boundaryCurve->position.values;
            if (cps.getSize() > 4 &&
                (cps.getValue(3) - cps.getValue(0)).stableNorm() < 1E-10 &&
                (cps.getValue(2) - cps.getValue(0)).stableNorm() < 1E-10 &&
                (cps.getValue(1) - cps.getValue(0)).stableNorm() < 1E-10)
            {
                cps.removeFirst(3);
                edge->boundaryCurve->position.parameters.removeFirst(1);
            }
            if (cps.getSize() > 4 &&
                (cps.getValue(cps.getSize() - 4) - cps.getValue(cps.getSize() - 1)).stableNorm() < 1E-10 &&
                (cps.getValue(cps.getSize() - 3) - cps.getValue(cps.getSize() - 1)).stableNorm() < 1E-10 &&
                (cps.getValue(cps.getSize() - 2) - cps.getValue(cps.getSize() - 1)).stableNorm() < 1E-10)
            {
                cps.removeLast(3);
                edge->boundaryCurve->position.parameters.removeLast(1);
            }
        }

        // for each vertex, sort the exitant edges in ascending order by their angle to the positive x-axis (counter-clockwise).
        for (auto& vertex : vertices)
        {
            // make tuples <angle, edge>
            std::vector<std::pair<double, std::shared_ptr<Edge>>> edges;
            for (auto a : vertex->edges)
            {
                auto sa = a.lock();
                Eigen::Vector2d tangent_a;
                if (sa->vertexA.lock() == vertex)
                    tangent_a = sa->boundaryCurve->position.sample_dt(sa->boundaryCurve->position.domain.min()[0]);
                else if (sa->vertexB.lock() == vertex)
                    tangent_a = -sa->boundaryCurve->position.sample_dt(sa->boundaryCurve->position.domain.max()[0]);
                // angle is in [-pi, pi] with 0 being (1,0)
                double angle_a = std::atan2(tangent_a.y(), tangent_a.x());
                edges.push_back(std::make_pair(angle_a, sa));
            }
            // sort them
            std::sort(edges.begin(), edges.end());
            // remove duplicate edges at vertex if necessary
            vertex->edges.clear();
            double last_angle = std::numeric_limits<double>::max();
            for (size_t ie = 0; ie < edges.size(); ++ie)
            {
                if (edges[ie].first != last_angle)
                {
                    vertex->edges.push_back(edges[ie].second);
                    last_angle = edges[ie].first;
                }
            }
        }
    }

    void EdgeGraph::deleteInvalid(std::vector<std::shared_ptr<EdgeGraph::Vertex>>& _vertices,
                                  std::vector<std::shared_ptr<EdgeGraph::Edge>>& _edges,
                                  double _epsilon)
    {
        // delete invalidated edges
        for (auto itEdge = _edges.begin(); itEdge != _edges.end();)
        {
            auto vertA = (*itEdge)->vertexA.lock();
            auto vertB = (*itEdge)->vertexB.lock();
            if ((vertA->position - vertB->position).stableNorm() < _epsilon && (*itEdge)->boundaryCurve->position.getBoundingBox().diagonal().stableNorm() < _epsilon)
            {
                auto edge = *itEdge;
                vertA->edges.remove_if([edge](std::weak_ptr<EdgeGraph::Edge> e)
                                       { return !(edge.owner_before(e) || e.owner_before(edge)); });
                vertB->edges.remove_if([edge](std::weak_ptr<EdgeGraph::Edge> e)
                                       { return !(edge.owner_before(e) || e.owner_before(edge)); });
                itEdge = _edges.erase(itEdge);
            }
            else
                ++itEdge;
        }

        // delete invalidated vertices
        for (auto itVertex = _vertices.begin(); itVertex != _vertices.end();)
        {
            if ((*itVertex)->edges.empty())
                itVertex = _vertices.erase(itVertex);
            else
                ++itVertex;
        }
    }

    void EdgeGraph::redirectEdges(std::shared_ptr<EdgeGraph::Vertex> _from, std::shared_ptr<EdgeGraph::Vertex> _to)
    {
        // redirect the edges that "from" is connected to
        for (auto from_edge : _from->edges)
        {
            auto from_edge_s = from_edge.lock();
            if (!from_edge_s)
                continue;
            if (from_edge_s->vertexA.lock() == _from)
            {
                from_edge_s->vertexA = _to;
                from_edge_s->boundaryCurve->position.values.setValue(0, _to->position);
                from_edge_s->boundaryCurve->position.recomputeBoundingBox();
            }
            if (from_edge_s->vertexB.lock() == _from)
            {
                from_edge_s->vertexB = _to;
                from_edge_s->boundaryCurve->position.values.setValue(from_edge_s->boundaryCurve->position.values.getSize() - 1, _to->position);
                from_edge_s->boundaryCurve->position.recomputeBoundingBox();
            }
            if (std::find_if(_to->edges.begin(), _to->edges.end(), [from_edge_s](std::weak_ptr<EdgeGraph::Edge> e)
                             { return !(from_edge_s.owner_before(e) || e.owner_before(from_edge_s)); }) == _to->edges.end())
                _to->edges.push_back(from_edge_s);
        }
        _from->edges.clear();
    }
}
