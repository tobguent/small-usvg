#include <usvg/patch.hpp>

#include <usvg/bezier_discretizer.hpp>
#include <usvg/gradient_mesh.hpp>
#include <usvg/poisson_curve.hpp>

#include <set>

namespace usvg
{

    Patch::BoundaryCurve::BoundaryCurve(const PiecewiseLinearCurve2d& _position, const CubicBezierSpline3d& _color, const EBoundaryCondition _boundaryCondition, std::shared_ptr<GradientMesh> _gradientMesh)
        : position(_position)
        , color(_color)
        , boundaryCondition(_boundaryCondition)
        , gradientMesh(_gradientMesh)
    {
    }

    bool Patch::BoundaryCurve::isValid() const
    {
        // positions valid?
        if (!position.isValid())
            return false;

        // colors valid?
        if (!color.isValid())
            return false;

        // gradient mesh valid?
        if (gradientMesh && !gradientMesh->isValid())
            return false;

        // all good
        return true;
    }

    Patch::PoissonCurveDiscrete::PoissonCurveDiscrete(const PiecewiseLinearCurve2d& _position, const PiecewiseLinearCurve3d& _weights)
        : position(_position)
        , weights(_weights)
    {
    }

    bool Patch::PoissonCurveDiscrete::isValid() const
    {
        // positions valid?
        if (!position.isValid())
            return false;

        // weight valid?
        if (!weights.isValid())
            return false;

        // all good
        return true;
    }

    Patch::Loop::Loop()
        : turningNumber(std::numeric_limits<double>::infinity())
    {
    }

    Patch::Loop::Loop(const std::vector<BoundaryCurve>& _boundaryCurves, double _turningNumber)
        : boundaryCurves(_boundaryCurves)
        , turningNumber(_turningNumber)
    {
        recomputeBoundingBox();
    }

    bool Patch::Loop::isValid() const
    {
        // boundary curves valid?
        for (auto& curves : boundaryCurves)
            if (!curves.isValid())
                return false;

        // is a closed loop?
        if (!boundaryCurves.empty())
        {
            if ((boundaryCurves.front().position.values.first() - boundaryCurves.back().position.values.last()).stableNorm() > 1E-10)
                return false;
        }

        // bounding box valid?
        if (mBoundingBox.isEmpty())
            return false;

        // turning number initialized?
        if (turningNumber == std::numeric_limits<double>::infinity())
            return false;

        // all good
        return true;
    }

    bool Patch::Loop::isInterior() const noexcept
    {
        return turningNumber <= -0.5;
    }

    bool Patch::Loop::isExterior() const noexcept
    {
        return 0.5 <= turningNumber;
    }

    bool Patch::Loop::isFlat() const noexcept
    {
        return -0.5 < turningNumber && turningNumber < 0.5;
    }

    bool Patch::Loop::contains(const Eigen::Vector2d& _point, double _epsilon) const
    {
        // take absolute since we don't mind the orientation of the traversal
        double theta = std::abs(windingNumber(_point));
        return 1 - _epsilon < theta && theta < 1 + _epsilon;
    }

    bool Patch::Loop::contains(const Loop& _other, double _epsilon) const
    {
        // is bounding box fully contained?
        if (!mBoundingBox.contains(_other.mBoundingBox))
            return false;
        // loops are by construction crossing-free. If one point is contained, then all points are contained
        if (_other.boundaryCurves.empty() || _other.boundaryCurves[0].position.values.getSize() == 0)
            return false;
        return contains(_other.boundaryCurves[0].position.values.getValue(0), _epsilon);
    }

    const Eigen::AlignedBox2d& Patch::Loop::getBoundingBox() const noexcept
    {
        return mBoundingBox;
    }

    /**
     * @brief Recomputes the bounding box of the boundary curves.
     */
    void Patch::Loop::recomputeBoundingBox() noexcept
    {
        mBoundingBox.setEmpty();
        for (auto& curve : boundaryCurves)
            mBoundingBox.extend(curve.position.getBoundingBox());
    }

    /**
     * @brief Calculates the winding number for a point and a given loop.
     * @param _point Point to compute the winding number for.
     * @return Signed winding number.
     */
    [[nodiscard]] double Patch::Loop::windingNumber(const Eigen::Vector2d& _point) const
    {
        // early out
        if (!mBoundingBox.contains(_point))
            return 0;

        // loop the boundary
        double theta = 0;
        for (auto& line : boundaryCurves)
        {
            auto vertices = line.position.values;
            for (Eigen::Index i = 0; i < vertices.getSize() - 1; ++i)
            {
                Eigen::Vector2d a = (vertices.getValue(i) - _point).stableNormalized();
                Eigen::Vector2d b = (vertices.getValue(i + 1) - _point).stableNormalized();
                theta += std::atan2(a.x() * b.y() - a.y() * b.x(), a.x() * b.x() + a.y() * b.y()) / (2 * EIGEN_PI);
            }
        }
        return theta;
    }

    double Patch::Loop::computeTurningNumber(const Eigen::Vector2d& _point1, const Eigen::Vector2d& _point2, const Eigen::Vector2d& _point3) noexcept
    {
        // in the event that we travel from point1 via point2 back to point1, we have a half loop:
        if ((_point1 - _point3).stableNorm() < 1E-10)
            return 0.5;

        Eigen::Vector2d dir1 = (_point2 - _point1).stableNormalized();
        Eigen::Vector2d dir2 = (_point3 - _point2).stableNormalized();
        double angle         = std::acos(dir1.dot(dir2));
        if (dir1.x() * dir2.y() - dir1.y() * dir2.x() < 0)
            angle *= -1;
        return angle / (2. * EIGEN_PI);
    }

    double Patch::Loop::computeTurningNumber(const std::vector<BoundaryCurve>& _boundaryCurves)
    {
        size_t numCurves     = _boundaryCurves.size();
        double turningNumber = 0;
        std::vector<Eigen::Vector2d> points;
        for (auto& boundaryCurve : _boundaryCurves)
        {
            for (auto& vertex : boundaryCurve.position.values.getData())
            {
                points.push_back(vertex);
            }
            points.pop_back();
        }
        points.push_back(points[0]);
        points.push_back(points[1]);
        for (size_t i = 0; i < points.size() - 2; ++i)
        {
            turningNumber += computeTurningNumber(points[i], points[i + 1], points[i + 2]);
        }
        return turningNumber;
    }

    bool Patch::isValid() const
    {
        // closed loops valid?
        for (auto& loop : loops)
            if (!loop->isValid())
                return false;

        // gradient meshes valid?
        for (auto& gradientMesh : gradientMeshes)
            if (!gradientMesh->isValid())
                return false;

        // Poisson curves valid?
        for (auto& poissonCurve : poissonCurves)
            if (!poissonCurve->isValid())
                return false;

        return true;
    }

    bool Patch::contains(const Eigen::Vector2d& _point, double _epsilon) const
    {
        for (auto& loop : loops)
        {
            if (loop->isInterior() && !loop->contains(_point, _epsilon))
                return false;
            if (loop->isExterior() && loop->contains(_point, _epsilon))
                return false;
        }
        return true;
    }

    Patch::BoundaryCurve Patch::makeLeftPatchBoundaryCurve(const EdgeGraph::Edge& _edge)
    {
        // get position and left color
        auto position = _edge.boundaryCurve->position;
        auto color    = _edge.boundaryCurve->colorLeft;

        // reverse the order, since patch boundary curves are oriented such that the color is on the right
        position.values.reverse();
        position.parameters.reverse();
        color.controlPoints.reverse();
        color.parameters.reverse();
        for (Eigen::Index ip = 0; ip < position.parameters.getSize(); ++ip)
        {
            Eigen::Vector1d parameter = position.parameters.getValue(ip);
            double t                  = (parameter[0] - position.domain.min()[0]) / (position.domain.max()[0] - position.domain.min()[0]);
            parameter                 = position.domain.min() + (1 - t) * (position.domain.max() - position.domain.min());
            position.parameters.setValue(ip, parameter);
        }
        for (Eigen::Index ip = 0; ip < color.parameters.getSize(); ++ip)
        {
            Eigen::Vector1d parameter = color.parameters.getValue(ip);
            double t                  = (parameter[0] - color.domain.min()[0]) / (color.domain.max()[0] - color.domain.min()[0]);
            parameter                 = color.domain.min() + (1 - t) * (color.domain.max() - color.domain.min());
            color.parameters.setValue(ip, parameter);
        }

        // return the patch boundary curve
        return Patch::BoundaryCurve(
            position,
            color,
            _edge.boundaryCurve->boundaryConditionLeft,
            nullptr);
    }

    Patch::BoundaryCurve Patch::makeRightPatchBoundaryCurve(const EdgeGraph::Edge& _edge)
    {
        // get position and right color
        auto position = _edge.boundaryCurve->position;
        auto color    = _edge.boundaryCurve->colorRight;

        // return the patch boundary curve
        return Patch::BoundaryCurve(
            position,
            color,
            _edge.boundaryCurve->boundaryConditionRight,
            _edge.boundaryCurve->gradientMesh);
    }

    std::vector<std::shared_ptr<Patch::Loop>> Patch::computeLoops(const EdgeGraph& _graph)
    {
        // ------------------------------------------
        // prune edges for turning number computation
        // ------------------------------------------

        // We make a map to have an indirection from the edge to its index in the array. This is because we later build sets of those indices.
        // Those sets are automatically sorted by the index. If we directly used the edge pointers, their order would be non-deterministic since the order would depend on their address in memory.
        std::map<std::shared_ptr<EdgeGraph::Edge>, int> edgeToIndex;
        int index = 0;
        for (auto& e : _graph.edges)
            edgeToIndex.insert(std::make_pair(e, index++));

        // first traverse from all vertices with valence one until we reach a point with valence !=2 and mark the edges along the way as prunedEdges
        std::set<int> prunedEdges;
        for (auto vertex : _graph.vertices)
        {
            if (vertex->edges.size() != 1)
                continue;
            std::shared_ptr<EdgeGraph::Vertex> currVertex = vertex;
            std::shared_ptr<EdgeGraph::Edge> currEdge     = vertex->edges.front().lock();

            // if the edge directly connects back onto itself, then it should not be pruned
            if (currEdge->vertexA.lock() == currEdge->vertexB.lock())
            {
                continue;
            }

            while (true)
            {
                // prune the edge
                prunedEdges.insert(edgeToIndex[currEdge]);

                // get next vertex
                std::shared_ptr<EdgeGraph::Vertex> nextVertex;
                if (currEdge->vertexA.lock() == currVertex)
                    nextVertex = currEdge->vertexB.lock();
                else
                    nextVertex = currEdge->vertexA.lock();

                // found a vertex with valence != 2
                if (nextVertex->edges.size() != 2)
                    break;

                // get next edge
                std::shared_ptr<EdgeGraph::Edge> nextEdge = nullptr;
                if (nextVertex->edges.front().lock() == currEdge)
                    nextEdge = nextVertex->edges.back().lock();
                else
                    nextEdge = nextVertex->edges.front().lock();

                // advance state
                currVertex = nextVertex;
                currEdge   = nextEdge;
            };
        }

        // during traversal all edges are visited twice. at the edges to pending lists
        std::set<int> pendingLeftVisit;
        std::set<int> pendingRightVisit;
        for (auto e : _graph.edges)
        {
            pendingLeftVisit.insert(edgeToIndex[e]);
            pendingRightVisit.insert(edgeToIndex[e]);
        }

        // --------------------------
        // detect loops in edge graph
        // --------------------------

        // traverse the graph and detected all loops
        std::vector<std::shared_ptr<Patch::Loop>> loops;
        std::set<int> pendingLeftVisitTurning  = pendingLeftVisit;
        std::set<int> pendingRightVisitTurning = pendingRightVisit;
        while (!pendingLeftVisit.empty() || !pendingRightVisit.empty())
        {
            // get the next edge to start traversal from
            std::shared_ptr<EdgeGraph::Edge> currEdge = nullptr;
            std::shared_ptr<EdgeGraph::Vertex> currVertex1, currVertex2; // start and end point of oriented edge
            std::vector<Patch::BoundaryCurve> patchBoundaryCurves;

            if (!pendingLeftVisit.empty())
            {
                // pop next edge
                currEdge = _graph.edges[*pendingLeftVisit.begin()];
                pendingLeftVisit.erase(edgeToIndex[currEdge]);
                pendingLeftVisitTurning.erase(edgeToIndex[currEdge]);

                // make patch boundary curve
                patchBoundaryCurves.push_back(
                    makeLeftPatchBoundaryCurve(*currEdge));

                // get edge info and orient consistently
                currVertex1 = currEdge->vertexB.lock();
                currVertex2 = currEdge->vertexA.lock();
            }
            else
            {
                // pop next edge
                currEdge = _graph.edges[*pendingRightVisit.begin()];
                pendingRightVisit.erase(edgeToIndex[currEdge]);
                pendingRightVisitTurning.erase(edgeToIndex[currEdge]);

                // make patch boundary curve
                patchBoundaryCurves.push_back(
                    makeRightPatchBoundaryCurve(*currEdge));

                // get edge info and orient consistently
                currVertex1 = currEdge->vertexA.lock();
                currVertex2 = currEdge->vertexB.lock();
            }

            // loop around and record all non-pruned boundaries
            std::vector<Patch::BoundaryCurve> nonPrunedBoundaries;
            do
            {
                if (prunedEdges.find(edgeToIndex[currEdge]) == prunedEdges.end())
                {
                    nonPrunedBoundaries.push_back(patchBoundaryCurves.back());
                }

                // get the next edge, turning right
                auto nextEdge = currVertex2->nextEdge(currEdge);

                // is the edge oriented the right way? (if not, the interior we are interested in "isLeft")
                bool isLeft = nextEdge->vertexB.lock() == currVertex2;

                // make patch boundary curve
                if (isLeft)
                {
                    if (pendingLeftVisit.find(edgeToIndex[nextEdge]) == pendingLeftVisit.end())
                        break;
                    patchBoundaryCurves.push_back(
                        makeLeftPatchBoundaryCurve(*nextEdge));
                    currVertex1 = nextEdge->vertexB.lock();
                    currVertex2 = nextEdge->vertexA.lock();

                    // pop next edge
                    currEdge = nextEdge;
                    pendingLeftVisit.erase(edgeToIndex[currEdge]);
                }
                else
                {
                    if (pendingRightVisit.find(edgeToIndex[nextEdge]) == pendingRightVisit.end())
                        break;
                    patchBoundaryCurves.push_back(
                        makeRightPatchBoundaryCurve(*nextEdge));
                    currVertex1 = nextEdge->vertexA.lock();
                    currVertex2 = nextEdge->vertexB.lock();

                    // pop next edge
                    currEdge = nextEdge;
                    pendingRightVisit.erase(edgeToIndex[currEdge]);
                }

            } while (true);

            // compute the turning number from the non-pruned boundaries
            double turningNumber = nonPrunedBoundaries.empty() ? 0 : Patch::Loop::computeTurningNumber(nonPrunedBoundaries);
            // add loop
            loops.push_back(std::make_shared<Patch::Loop>(patchBoundaryCurves, turningNumber));
        }
        return loops;
    }

    std::vector<std::shared_ptr<Patch>> Patch::computePatches(const std::vector<std::shared_ptr<Loop>>& _loops, const std::vector<std::shared_ptr<PoissonCurve>>& _poissonCurves, double _discretizationResidual)
    {
        // allocate a patch for the domain
        auto domainPatch = std::make_shared<Patch>();
        // allocate set for gradient meshes for the domain
        std::set<std::shared_ptr<GradientMesh>> domainPatch_gradientMeshes;

        // make a vector of patch pointers
        std::vector<std::shared_ptr<Patch>> patch_candiates;
        patch_candiates.resize(_loops.size());
        for (int i = 0; i < _loops.size(); ++i)
        {
            // if a loop is an inner loop, then it will generate a patch
            bool inner_i = _loops[i]->isInterior();
            if (inner_i)
                patch_candiates[i] = std::make_shared<Patch>();
        }
        // also allocate sets of gradient meshes for each loop
        std::vector<std::set<std::shared_ptr<GradientMesh>>> patch_gradientMeshes;
        patch_gradientMeshes.resize(_loops.size());

        // sort all the loops by bounding box volume in ascending order.
        // if a loop is contained in another one, then so is the bounding box
        std::vector<std::pair<double, int>> loop_sorting;
        loop_sorting.resize(_loops.size());
        for (int l = 0; l < _loops.size(); ++l)
        {
            double volume = _loops[l]->getBoundingBox().volume();
            if (_loops[l]->isExterior())
                volume += 1E-10;
            loop_sorting[l] = std::make_pair(volume, l);
        }
        std::sort(loop_sorting.begin(), loop_sorting.end());

        // loop over all patches in ascending order by their volume
        for (int _i = 0; _i < loop_sorting.size(); ++_i)
        {
            int i = loop_sorting[_i].second;

            // if this is an inner loop, then insert into self
            if (patch_candiates[i])
            {
                patch_candiates[i]->loops.push_back(_loops[i]);
                for (auto curve : _loops[i]->boundaryCurves)
                {
                    if (curve.gradientMesh)
                        patch_gradientMeshes[i].insert(curve.gradientMesh);
                }
                continue;
            }

            bool inserted = false;
            // loop over all other patches
            for (int _j = 0; _j < loop_sorting.size(); ++_j)
            {
                // only consider those that are inner loops
                int j = loop_sorting[_j].second;
                if (patch_candiates[j] == nullptr)
                    continue;

                // if the inner loop contains patch
                double windingEpsilon = 1E-5;
                if (_loops[j]->contains(*_loops[i], windingEpsilon))
                {
                    // insert all the domain boundary curves and gradient meshes into the loop j (inner loop)
                    patch_candiates[j]->loops.push_back(_loops[i]);
                    for (auto curve : _loops[i]->boundaryCurves)
                    {
                        if (curve.gradientMesh)
                            patch_gradientMeshes[j].insert(curve.gradientMesh);
                    }
                    inserted = true;
                    break;
                }
            }
            // if i could not be inserted into any other patch
            if (!inserted)
            {
                // insert all the domain boundary curves and gradient meshes into the domain patch
                domainPatch->loops.push_back(_loops[i]);
                for (auto curve : _loops[i]->boundaryCurves)
                {
                    if (curve.gradientMesh)
                        domainPatch_gradientMeshes.insert(curve.gradientMesh);
                }
            }
        }

        // assign the unique set of gradient meshes into the patches (we did this to avoid duplicates)
        for (int i = 0; i < patch_candiates.size(); ++i)
        {
            if (patch_candiates[i])
                patch_candiates[i]->gradientMeshes.insert(
                    patch_candiates[i]->gradientMeshes.end(),
                    patch_gradientMeshes[i].begin(),
                    patch_gradientMeshes[i].end());
        }
        domainPatch->gradientMeshes.insert(
            domainPatch->gradientMeshes.end(),
            domainPatch_gradientMeshes.begin(),
            domainPatch_gradientMeshes.end());

        // lastly, add the patches to the scene
        std::vector<std::shared_ptr<Patch>> patches;
        for (auto& pc : patch_candiates)
            if (pc != nullptr)
                patches.push_back(pc);
        patches.push_back(domainPatch);

        // for each input poisson curve
        for (auto& pc : _poissonCurves)
        {
            // discreteize the poisson curve
            PiecewiseLinearCurve2d discrete;
            BezierDiscretizer2d::topDown(pc->position, _discretizationResidual, discrete);
            auto pcd = std::make_shared<PoissonCurveDiscrete>(
                discrete,
                pc->weights);

            // add to all patches for now
            for (auto& patch : patches)
            {
                patch->poissonCurves.push_back(pcd);
            }
        }

        return patches;
    }
}
