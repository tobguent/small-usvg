#include <usvg/scene.hpp>

#include <usvg/diffusion_curve.hpp>
#include <usvg/edge_graph.hpp>
#include <usvg/gradient_mesh.hpp>
#include <usvg/image.hpp>
#include <usvg/jacobi_iteration.hpp>
#include <usvg/patch.hpp>
#include <usvg/poisson_curve.hpp>

#include "tinyxml2.h"
#include <fstream>
#include <iostream>

#include <Eigen/Eigen>

namespace usvg
{
    Scene::Scene()
        : blendOperation(EBlendOperation::Zero)
    {
    }

    bool Scene::readDiffusionCurves(const tinyxml2::XMLElement* _parent_element, bool _swap, Scene* _scene)
    {
        // read the individual curves from file
        int num_curves = 0;
        _parent_element->QueryIntAttribute("nb_curves", &num_curves);
        const tinyxml2::XMLElement* curve_element = _parent_element->FirstChildElement("curve");
        for (int i = 0; i < num_curves; i++)
        {
            DiffusionCurve curve;
            if (curve_element == nullptr)
            {
                std::cout << "Cannot read curve " + std::to_string(i) << std::endl;
                return false;
            }

            // -------------------- Read control points
            int num_control_points = 0;
            curve_element->QueryIntAttribute("nb_control_points", &num_control_points);
            const tinyxml2::XMLElement* point_set_spec_element = curve_element->FirstChildElement("control_points_set");
            if (point_set_spec_element == nullptr)
            {
                std::cout << "Cannot read control points of diffusion curve " + std::to_string(i) << std::endl;
                return false;
            }
            if (!readPositions(curve.position.controlPoints, point_set_spec_element, "control_point", num_control_points, false, _scene->domain, _swap))
                return false;
            curve.position.setUniformParameterization();
            // normalize such that the parameter is in [0,1]
            for (int i = 0; i < curve.position.parameters.getSize(); ++i)
                curve.position.parameters.setValue(i, curve.position.parameters.getValue(i) / curve.position.parameters.last());
            curve.position.setExplicitParameterization();
            curve.position.recomputeBoundingBox();

            // -------------------- Read left colors
            int num_colors_left = 0;
            curve_element->QueryIntAttribute("nb_left_colors", &num_colors_left);
            const tinyxml2::XMLElement* left_color_set_spec_element = curve_element->FirstChildElement("left_colors_set");
            if (left_color_set_spec_element == nullptr)
            {
                std::cout << "Cannot read left colors of diffusion curve " + std::to_string(i) << std::endl;
                return false;
            }

            // by default, the boundary type is Dirichlet
            EBoundaryCondition& boundary_left = _swap ? curve.boundaryConditionRight : curve.boundaryConditionLeft;
            boundary_left                     = EBoundaryCondition::Dirichlet;
            if (left_color_set_spec_element->Attribute("boundary") != nullptr)
            {
                std::string boundary = left_color_set_spec_element->Attribute("boundary");
                if (boundary == "Neumann")
                    boundary_left = EBoundaryCondition::Neumann;
            }

            // read color control points
            if (!readPiecewiseLinearCurve3d(_swap ? curve.colorRight : curve.colorLeft, left_color_set_spec_element, "left_color", num_colors_left, _swap))
                return false;

            // -------------------- Read right colors
            int num_colors_right = 0;
            curve_element->QueryIntAttribute("nb_right_colors", &num_colors_right);
            const tinyxml2::XMLElement* right_color_spec_element = curve_element->FirstChildElement("right_colors_set");
            if (right_color_spec_element == nullptr)
            {
                std::cout << "Cannot read right colors of diffusion curve " + std::to_string(i) << std::endl;
                return false;
            }

            // by default, the boundary type is Dirichlet
            EBoundaryCondition& boundary_right = _swap ? curve.boundaryConditionLeft : curve.boundaryConditionRight;
            boundary_right                     = EBoundaryCondition::Dirichlet;
            if (right_color_spec_element->Attribute("boundary") != nullptr)
            {
                std::string boundary = right_color_spec_element->Attribute("boundary");
                if (boundary == "Neumann")
                    boundary_right = EBoundaryCondition::Neumann;
            }

            // read color control points
            if (!readPiecewiseLinearCurve3d(_swap ? curve.colorLeft : curve.colorRight, right_color_spec_element, "right_color", num_colors_right, _swap))
                return false;

            // -------------------- Set diffusion curve
            if (!curve.isValid())
            {
                std::cout << "Diffusion curve is not valid." << std::endl;
                return false;
            }

            _scene->diffusionCurves.push_back(std::make_shared<DiffusionCurve>(curve));
            curve_element = curve_element->NextSiblingElement("curve");
        }
        return true;
    }

    bool Scene::readPoissonCurves(const tinyxml2::XMLElement* _parent_element, Scene* _scene)
    {
        int num_curves = 0;
        _parent_element->QueryIntAttribute("nb_curves", &num_curves);
        const tinyxml2::XMLElement* curve_element = _parent_element->FirstChildElement("poisson_curve");
        for (int i = 0; i < num_curves; i++)
        {
            PoissonCurve curve;
            if (curve_element == nullptr)
            {
                std::cout << "Cannot read curve " + std::to_string(i) << std::endl;
                return false;
            }

            int num_control_points = 0;
            curve_element->QueryIntAttribute("nb_control_points", &num_control_points);
            const tinyxml2::XMLElement* point_set_spec_element = curve_element->FirstChildElement("control_points_set");
            if (point_set_spec_element == nullptr)
            {
                std::cout << "Cannot read control points of diffusion curve " + std::to_string(i) << std::endl;
                return false;
            }

            // -------------------- Read control points
            const tinyxml2::XMLElement* point = point_set_spec_element->FirstChildElement("control_point");
            if (!readPositions(curve.position.controlPoints, point_set_spec_element, "control_point", num_control_points, false, _scene->domain, false))
                return false;

            curve.position.setUniformParameterization();
            curve.position.recomputeBoundingBox();

            // -------------------- Read Poisson weights
            int num_weights = 0;
            curve_element->QueryIntAttribute("nb_weights", &num_weights);
            const tinyxml2::XMLElement* weight_set_spec_element = curve_element->FirstChildElement("weights_set");
            if (weight_set_spec_element == nullptr)
            {
                std::cout << "Cannot read weights of Poisson curve " + std::to_string(i) << std::endl;
                return false;
            }

            if (!readPiecewiseLinearCurve3d(curve.weights, weight_set_spec_element, "weight", num_weights, false))
                return false;

            // -------------------- Set Poisson curve
            _scene->poissonCurves.push_back(std::make_shared<PoissonCurve>(curve));
            curve_element = curve_element->NextSiblingElement("poisson_curve");
        }
        return true;
    }

    bool Scene::readGradientMeshes(const tinyxml2::XMLElement* _parent_element, Scene* _scene)
    {
        int num_meshes = 0;
        _parent_element->QueryIntAttribute("nb_meshes", &num_meshes);
        const tinyxml2::XMLElement* mesh_spec_element = _parent_element->FirstChildElement("mesh");
        for (int i = 0; i < num_meshes; i++)
        {
            GradientMesh mesh;
            if (mesh_spec_element == nullptr)
            {
                std::cout << "Cannot read mesh " + std::to_string(i) << std::endl;
                return false;
            }
            int num_tiles_rows = 0, num_tiles_cols = 0;
            mesh_spec_element->QueryIntAttribute("nb_rows", &num_tiles_rows);
            mesh_spec_element->QueryIntAttribute("nb_cols", &num_tiles_cols);
            int num_rows = num_tiles_rows + 1;
            int num_cols = num_tiles_cols + 1;

            // if the mesh is normalized, its positions are scaled by the image size
            bool is_normalized = false;
            mesh_spec_element->QueryBoolAttribute("normalized", &is_normalized);

            // -------------------- Read positions of the mesh vertices
            int num_positions = 0;
            mesh_spec_element->QueryIntAttribute("nb_positions", &num_positions);
            if (num_positions != num_rows * num_cols)
            {
                std::cout << "Number of positions does not match the mesh size in mesh " + std::to_string(i) << std::endl;
                return false;
            }

            const tinyxml2::XMLElement* vertex_set_spec_element = mesh_spec_element->FirstChildElement("position_set");
            if (vertex_set_spec_element == nullptr)
            {
                std::cout << "Cannot read positions of mesh " + std::to_string(i) << std::endl;
                return false;
            }

            const tinyxml2::XMLElement* position = vertex_set_spec_element->FirstChildElement("position");
            Array2d positions;
            if (!readPositions(positions, vertex_set_spec_element, "position", num_positions, is_normalized, _scene->domain, false))
                return false;

            // -------------------- Read colors of the mesh vertices
            int num_colors = 0;
            mesh_spec_element->QueryIntAttribute("nb_colors", &num_colors);
            if (num_colors != num_rows * num_cols)
            {
                std::cout << "Number of colors does not match the mesh size in mesh " + std::to_string(i) << std::endl;
                return false;
            }

            const tinyxml2::XMLElement* color_set_spec_element = mesh_spec_element->FirstChildElement("color_set");
            if (color_set_spec_element == nullptr)
            {
                std::cout << "Cannot read colors of mesh " + std::to_string(i) << std::endl;
                return false;
            }

            Array3d colors;
            if (!readColors(colors, color_set_spec_element, "color", num_colors, false))
                return false;

            // Read optional position tangents
            Array2d tangents_u;
            Array2d tangents_v;
            const tinyxml2::XMLElement* tangent_set_spec_element = mesh_spec_element->FirstChildElement("pos_tangent_set");
            if (tangent_set_spec_element != nullptr)
            {
                if (!readPositions(tangents_u, tangent_set_spec_element, "positionU", num_positions, is_normalized, _scene->domain, false))
                    return false;
                if (!readPositions(tangents_v, tangent_set_spec_element, "positionV", num_positions, is_normalized, _scene->domain, false))
                    return false;
            }

            // -------------------- Set gradient mesh
            mesh.resolution = Eigen::Vector2i(num_tiles_cols, num_tiles_rows);
            mesh.tiles.resize(num_tiles_rows * num_tiles_cols, GradientMesh::Tile());

            // copy positions and colors into patch
            for (int iy = 0; iy < num_rows; ++iy)
                for (int ix = 0; ix < num_cols; ++ix)
                {
                    Eigen::Index linearIndex = iy * num_cols + ix;
                    mesh.setPosition(ix, iy, positions.getValue(linearIndex));
                    mesh.setColor(ix, iy, colors.getValue(linearIndex));
                }
            // set derivatives or numerically estimate them
            for (int iy = 0; iy < num_rows; ++iy)
                for (int ix = 0; ix < num_cols; ++ix)
                {
                    Eigen::Index linearIndex = iy * num_cols + ix;
                    if (tangents_u.getSize() != 0)
                    {
                        // use read tangent
                        mesh.setPosition_du(ix, iy, tangents_u.getValue(linearIndex));
                    }
                    else
                    {
                        // numerically estimate
                        int ix0                = std::max(0, ix - 1);
                        int ix1                = std::min(ix + 1, num_cols - 1);
                        Eigen::Vector2d pos0   = mesh.getPosition(ix0, iy);
                        Eigen::Vector2d pos1   = mesh.getPosition(ix1, iy);
                        Eigen::Vector2d pos_du = (pos1 - pos0) / (ix1 - ix0);
                        mesh.setPosition_du(ix, iy, pos_du);
                    }

                    if (tangents_v.getSize() != 0)
                    {
                        // use read tangent
                        mesh.setPosition_dv(ix, iy, tangents_v.getValue(linearIndex));
                    }
                    else
                    {
                        // numerically estimate
                        int iy0                = std::max(0, iy - 1);
                        int iy1                = std::min(iy + 1, num_rows - 1);
                        Eigen::Vector2d pos0   = mesh.getPosition(ix, iy0);
                        Eigen::Vector2d pos1   = mesh.getPosition(ix, iy1);
                        Eigen::Vector2d pos_dv = (pos1 - pos0) / (iy1 - iy0);
                        mesh.setPosition_dv(ix, iy, pos_dv);
                    }

                    // numerically estimate color gradients
                    {
                        // numerically estimate
                        int ix0                  = std::max(0, ix - 1);
                        int ix1                  = std::min(ix + 1, num_cols - 1);
                        Eigen::Vector3d color0   = mesh.getColor(ix0, iy);
                        Eigen::Vector3d color1   = mesh.getColor(ix1, iy);
                        Eigen::Vector3d color_du = (color1 - color0) / (ix1 - ix0);
                        mesh.setColor_du(ix, iy, color_du);
                    }
                    {
                        // numerically estimate
                        int iy0                  = std::max(0, iy - 1);
                        int iy1                  = std::min(iy + 1, num_rows - 1);
                        Eigen::Vector3d color0   = mesh.getColor(ix, iy0);
                        Eigen::Vector3d color1   = mesh.getColor(ix, iy1);
                        Eigen::Vector3d color_dv = (color1 - color0) / (iy1 - iy0);
                        mesh.setColor_dv(ix, iy, color_dv);
                    }
                }

            mesh.updateBezierPosition();
            mesh.updateBezierColor();
            if (!mesh.isValid())
            {
                std::cout << "Gradient mesh is not valid." << std::endl;
                return false;
            }
            _scene->gradientMeshes.push_back(std::make_shared<GradientMesh>(mesh));

            mesh_spec_element = mesh_spec_element->NextSiblingElement("mesh");
        }
        return true;
    }

    bool Scene::readXML(const std::string& _path)
    {
        // get output
        this->diffusionCurves.clear();
        this->poissonCurves.clear();
        this->gradientMeshes.clear();
        this->domain.setEmpty();

        // read the doctype
        std::ifstream infile(_path);
        if (!infile.good())
        {
            std::cout << "Cannot load XML file: " + _path << std::endl;
            return false;
        }
        std::string doc_type;
        std::getline(infile, doc_type);
        infile.close();

        // read xml file
        tinyxml2::XMLDocument doc;
        if (doc.LoadFile(_path.c_str()) != tinyxml2::XML_SUCCESS)
        {
            std::cout << "Cannot load XML file: " + _path << std::endl;
            return false;
        }

        // unified scene reader
        if (doc_type == "<!DOCTYPE SceneXML>")
        {
            const tinyxml2::XMLElement* root_element = doc.FirstChildElement("scene");
            if (root_element == nullptr)
            {
                std::cout << "Cannot find scene in XML file" << std::endl;
                return false;
            }

            int image_width, image_height;
            root_element->QueryIntAttribute("image_width", &image_width);
            root_element->QueryIntAttribute("image_height", &image_height);
            this->domain = Eigen::AlignedBox2d(Eigen::Vector2d(0, 0), Eigen::Vector2d(image_width, image_height));

            const tinyxml2::XMLElement* diffusion_curves_element = root_element->FirstChildElement("curve_set");
            if (diffusion_curves_element != nullptr)
            {
                if (!readDiffusionCurves(diffusion_curves_element, false, this))
                    return false;
            }

            const tinyxml2::XMLElement* poisson_curves_element = root_element->FirstChildElement("poisson_curve_set");
            if (poisson_curves_element != nullptr)
            {
                if (!readPoissonCurves(poisson_curves_element, this))
                    return false;
            }

            const tinyxml2::XMLElement* gradient_meshes_element = root_element->FirstChildElement("mesh_set");
            if (gradient_meshes_element != nullptr)
            {
                if (!readGradientMeshes(gradient_meshes_element, this))
                    return false;
            }
        }
        // Orzan reader
        else if (doc_type == "<!DOCTYPE CurveSetXML>")
        {
            const tinyxml2::XMLElement* root_element = doc.RootElement();
            if (root_element == nullptr)
            {
                std::cout << "Cannot find scene in XML file" << std::endl;
                return false;
            }

            int image_width, image_height;
            root_element->QueryIntAttribute("image_width", &image_width);
            root_element->QueryIntAttribute("image_height", &image_height);
            this->domain = Eigen::AlignedBox2d(Eigen::Vector2d(0, 0), Eigen::Vector2d(image_width, image_height));

            if (!readDiffusionCurves(root_element, true, this))
                return false;
        }
        else
        {
            std::cout << "Unrecognized DOCTYPE in XML" << std::endl;
            return false;
        }
        return true;
    }

    bool Scene::readPositions(Array2d& _points, const tinyxml2::XMLElement* _parent, const std::string& _child_name, int _num_points, bool _is_normalized, const Eigen::AlignedBox2d& _domain, bool _swap)
    {
        _points.setSize(_num_points);
        const tinyxml2::XMLElement* child_element = _parent->FirstChildElement(_child_name.c_str());
        if (child_element == nullptr)
        {
            std::cout << "Cannot read " + _child_name << std::endl;
            return false;
        }
        for (int i = 0; i < _num_points; i++)
        {
            if (child_element == nullptr)
            {
                std::cout << "Cannot read " + _child_name + " " + std::to_string(i) << std::endl;
                return false;
            }

            point_type p;
            child_element->QueryDoubleAttribute("x", &p[0]);
            child_element->QueryDoubleAttribute("y", &p[1]);

            if (_is_normalized)
            {
                p[1] = _domain.min().x() + (_domain.max().x() - _domain.min().x()) * p[1];
                p[0] = _domain.min().y() + (_domain.max().y() - _domain.min().y()) * p[0];
            }

            if (_swap)
            {
                // swap x and y components
                std::swap(p[0], p[1]);
            }

            _points.setValue(i, p);
            child_element = child_element->NextSiblingElement(_child_name.c_str());
        }
        return true;
    }

    bool Scene::readColors(Array3d& _colors, const tinyxml2::XMLElement* _parent, const std::string& _child_name, int _num_colors, bool _swap)
    {
        _colors.setSize(_num_colors);
        const tinyxml2::XMLElement* child_element = _parent->FirstChildElement(_child_name.c_str());
        if (child_element == nullptr)
        {
            std::cout << "Cannot read " + _child_name << std::endl;
            return false;
        }

        for (int i = 0; i < _num_colors; i++)
        {
            if (child_element == nullptr)
            {
                std::cout << "Cannot read " + _child_name + " " + std::to_string(i) << std::endl;
                return false;
            }

            color_type c;
            if (child_element->Attribute("R") != nullptr)
            {
                child_element->QueryDoubleAttribute("R", &c[0]);
                c[0] /= 255.0;
            }
            else
            {
                child_element->QueryDoubleAttribute("r", &c[0]);
            }

            if (child_element->Attribute("G") != nullptr)
            {
                child_element->QueryDoubleAttribute("G", &c[1]);
                c[1] /= 255.0;
            }
            else
            {
                child_element->QueryDoubleAttribute("g", &c[1]);
            }

            if (child_element->Attribute("B") != nullptr)
            {
                child_element->QueryDoubleAttribute("B", &c[2]);
                c[2] /= 255.0;
            }
            else
            {
                child_element->QueryDoubleAttribute("b", &c[2]);
            }

            if (_swap)
            {
                // swap R and B components
                std::swap(c[0], c[2]);
            }

            _colors.setValue(i, c);
            child_element = child_element->NextSiblingElement(_child_name.c_str());
        }
        return true;
    }

    bool Scene::readPiecewiseLinearCurve3d(PiecewiseLinearCurve3d& _polyline, const tinyxml2::XMLElement* _parent, const std::string& _child_name, int _num_points, bool _swap)
    {
        const tinyxml2::XMLElement* child_element = _parent->FirstChildElement(_child_name.c_str());
        if (child_element == nullptr)
        {
            std::cout << "Cannot read " + _child_name << std::endl;
            return false;
        }

        double maxT = -std::numeric_limits<double>::max();
        for (int i = 0; i < _num_points; i++)
        {
            if (child_element == nullptr)
            {
                std::cout << "Cannot read " + _child_name + " " + std::to_string(i) << std::endl;
                return false;
            }

            color_point_type cp;
            if (child_element->Attribute("R") != nullptr)
            {
                child_element->QueryDoubleAttribute("R", &cp[0]);
                cp[0] /= 255.0;
            }
            else
            {
                child_element->QueryDoubleAttribute("r", &cp[0]);
            }

            if (child_element->Attribute("G") != nullptr)
            {
                child_element->QueryDoubleAttribute("G", &cp[1]);
                cp[1] /= 255.0;
            }
            else
            {
                child_element->QueryDoubleAttribute("g", &cp[1]);
            }

            if (child_element->Attribute("B") != nullptr)
            {
                child_element->QueryDoubleAttribute("B", &cp[2]);
                cp[2] /= 255.0;
            }
            else
            {
                child_element->QueryDoubleAttribute("b", &cp[2]);
            }

            if (_swap)
            {
                // swap R and B components
                std::swap(cp[0], cp[2]);
            }

            child_element->QueryDoubleAttribute("globalID", &cp[3]);
            if (cp[3] > 1)
            {
                maxT = std::max(maxT, cp[3]);
            }

            _polyline.values.append(cp.xyz());
            _polyline.parameters.append(Eigen::Vector1d(cp.w()));
            child_element = child_element->NextSiblingElement(_child_name.c_str());
        }

        // Normalize the globalID
        if (maxT > 1)
        {
            for (Eigen::Vector1d& cp : _polyline.parameters.getData())
            {
                cp /= maxT;
            }
        }
        _polyline.setExplicitParameterization();
        _polyline.recomputeBoundingBox();
        if (!_polyline.isValid())
        {
            std::cout << "Piecewise linear curve is not valid." << std::endl;
            return false;
        }
        return true;
    }

    bool Scene::isValid() const
    {
        // diffusion curves valid?
        for (auto& diffusionCurve : diffusionCurves)
            if (!diffusionCurve->isValid())
                return false;

        // Poisson curves valid?
        for (auto& poissonCurve : poissonCurves)
            if (!poissonCurve->isValid())
                return false;

        // gradient meshes valid?
        for (auto& gradientMesh : gradientMeshes)
            if (!gradientMesh->isValid())
                return false;

        // domain defined? (this check requires at least some geometry)
        if (!diffusionCurves.empty() || !poissonCurves.empty() || !gradientMeshes.empty())
            if (this->domain.isEmpty())
                return false;

        // patches valid?
        for (auto& patch : patches)
            if (!patch->isValid())
                return false;

        // all good
        return true;
    }

    void Scene::buildPatches(double _vertexMergeThreshold, double _clippingEpsilon, double _discretizationResidual)
    {
        // build the discrete edge graph
        EdgeGraph graph(diffusionCurves, gradientMeshes, _vertexMergeThreshold, _clippingEpsilon, _discretizationResidual);

        // compute the loops
        const std::vector<std::shared_ptr<Patch::Loop>>& loops = Patch::computeLoops(graph);

        // compute the patches
        patches = Patch::computePatches(loops, poissonCurves, _discretizationResidual);
    }

    std::tuple<bool, Eigen::Vector3d, Eigen::Vector3d> Scene::sampleGradientMeshLaplacian(
        const GradientMesh& _mesh,
        Eigen::Vector2d& _coord,
        int _maxDepth,
        double _epsilon,
        const Eigen::Vector2d& _spacing,
        bool _open_x0, bool _open_x1, bool _open_y0, bool _open_y1)
    {
        // use a finite difference estimation of the Laplacian instead
        std::vector<Eigen::Vector3d> colors;
        _mesh.sampleColors(_coord, _maxDepth, _epsilon, colors);
        if (!colors.empty())
        {
            Eigen::Vector3d c_x0 = Eigen::Vector3d::Zero();
            int num_open         = 0;
            if (_open_x0)
            {
                std::vector<Eigen::Vector3d> colors;
                _mesh.sampleColors(_coord - Eigen::Vector2d(_spacing.x(), 0), _maxDepth, _epsilon, colors);
                if (!colors.empty())
                {
                    c_x0 = colors[0];
                    num_open++;
                }
            }
            Eigen::Vector3d c_x1 = Eigen::Vector3d::Zero();
            if (_open_x1)
            {
                std::vector<Eigen::Vector3d> colors;
                _mesh.sampleColors(_coord + Eigen::Vector2d(_spacing.x(), 0), _maxDepth, _epsilon, colors);
                if (!colors.empty())
                {
                    c_x1 = colors[0];
                    num_open++;
                }
            }
            Eigen::Vector3d c_y0 = Eigen::Vector3d::Zero();
            if (_open_y0)
            {
                std::vector<Eigen::Vector3d> colors;
                _mesh.sampleColors(_coord - Eigen::Vector2d(0, _spacing.y()), _maxDepth, _epsilon, colors);
                if (!colors.empty())
                {
                    c_y0 = colors[0];
                    num_open++;
                }
            }
            Eigen::Vector3d c_y1 = Eigen::Vector3d::Zero();
            if (_open_y1)
            {
                std::vector<Eigen::Vector3d> colors;
                _mesh.sampleColors(_coord + Eigen::Vector2d(0, _spacing.y()), _maxDepth, _epsilon, colors);
                if (!colors.empty())
                {
                    c_y1 = colors[0];
                    num_open++;
                }
            }
            // if we are at a valid pixel, also try to sample the four points around
            {
                // if this was successful, estimate the Laplacian numerically.
                // This discretization matches with the discretization later used by the Jacobi solver
                Eigen::Vector3d estimatedLaplacian =
                    (c_x0 + c_x1 + c_y0 + c_y1 - num_open * colors[0]) / (_spacing.prod());
                return std::make_tuple(true, colors[0], estimatedLaplacian);
            }
        }
        else
            return std::make_tuple(false, Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero());
    }

    bool Scene::intersect_segment_segment(const Eigen::Vector2d& _a, const Eigen::Vector2d& _b, const Eigen::Vector2d& _c, const Eigen::Vector2d& _d, Eigen::Vector2d& _t)
    {
        Eigen::Vector2d ba = _b - _a, dc = _d - _c;
        float disc = _a.x() * (_d.y() - _c.y()) + _b.x() * (_c.y() - _d.y()) + (_b.y() - _a.y()) * _d.x() + (_a.y() - _b.y()) * _c.x();
        if (abs(disc) < 1E-10)
            return false;

        _t.x() = (_a.x() * (_d.y() - _c.y()) + _c.x() * (_a.y() - _d.y()) + (_c.y() - _a.y()) * _d.x()) / (disc != 0 ? disc : 1);
        if (_t.x() < 0 || 1 < _t.x())
            return false;

        _t.y() = -(_a.x() * (_c.y() - _b.y()) + _b.x() * (_a.y() - _c.y()) + (_b.y() - _a.y()) * _c.x()) / (disc != 0 ? disc : 1);
        return 0 <= _t.y() && _t.y() <= 1;
    }

    void Scene::closest_point_on_segment(const Eigen::Vector2d& _b0, const Eigen::Vector2d& _b1, const Eigen::Vector2d& _p, double& _min_dist, double& _min_t)
    {
        Eigen::Vector2d u = _b1 - _b0;
        double t          = (_p - _b0).dot(u) / u.dot(u);
        t                 = std::min(std::max(0., t), 1.);
        Eigen::Vector2d q = (1 - t) * _b0 + t * _b1;
        double dist       = (_p - q).norm();
        if (dist < _min_dist)
        {
            _min_dist = dist;
            _min_t    = t;
        }
    }

    void Scene::initializeFromPatches(const Scene& scene, JacobiIteration3d& jacobi)
    {
        std::shared_ptr<const RegularGrid2d> grid = jacobi.getGrid();
        Eigen::Index numPixels                    = grid->getResolution().prod();
        Eigen::Vector2d spacing                   = grid->getSpacing();

        std::vector<int> patchId(numPixels, -1);
        Image gradientMeshColor(grid->getResolution());
        gradientMeshColor.setZero();

#ifdef NDEBUG
#pragma omp parallel for
#endif
        for (Eigen::Index linearIndex = 0; linearIndex < numPixels; ++linearIndex)
        {
            Eigen::Vector2i gridCoord = grid->getGridCoord(linearIndex);
            Eigen::Vector2d coord     = grid->getCoordAt(gridCoord);

            double closest_t = std::numeric_limits<double>::max();
            Eigen::Vector3d closest_color(1, 0, 0);
            EBoundaryCondition closest_condition = EBoundaryCondition::Neumann;

            for (int ipatch = 0; ipatch < scene.patches.size(); ++ipatch)
            {
                auto& patch = scene.patches[ipatch];
                for (auto& loop : patch->loops)
                {
                    for (auto& curve : loop->boundaryCurves)
                    {
                        int numVert = curve.position.values.getSize();
                        for (int ivert = 0; ivert < numVert - 1; ++ivert)
                        {
                            const Eigen::Vector2d& p0 = curve.position.values.getValue(ivert);
                            const Eigen::Vector2d& p1 = curve.position.values.getValue(ivert + 1);

                            Eigen::Vector2d tangent = p1 - p0;
                            if (tangent.stableNorm() < 1E-10)
                                continue;
                            Eigen::Vector2d normal(tangent.y(), -tangent.x());
                            bool isLeft = normal.dot(p0 - coord) > 0;
                            if (!isLeft)
                            {
                                // east
                                Eigen::Vector2d t;
                                if (intersect_segment_segment(coord, coord + Eigen::Vector2d(spacing.x(), 0), p0, p1, t))
                                {
                                    jacobi.setOpenEast(gridCoord, false);
                                    if (t.x() < closest_t)
                                    {
                                        closest_t = t.x();
                                        if (curve.boundaryCondition == EBoundaryCondition::Dirichlet)
                                        {
                                            double param  = curve.position.parameters.getValue(ivert).x() * (1 - t.y()) + curve.position.parameters.getValue(ivert + 1).x() * t.y();
                                            closest_color = curve.color.sample(param);
                                        }
                                        closest_condition = curve.boundaryCondition;
                                    }
                                }
                                // north
                                if (intersect_segment_segment(coord, coord + Eigen::Vector2d(0, spacing.y()), p0, p1, t))
                                {
                                    jacobi.setOpenNorth(gridCoord, false);
                                    if (t.x() < closest_t)
                                    {
                                        closest_t = t.x();
                                        if (curve.boundaryCondition == EBoundaryCondition::Dirichlet)
                                        {
                                            double param  = curve.position.parameters.getValue(ivert).x() * (1 - t.y()) + curve.position.parameters.getValue(ivert + 1).x() * t.y();
                                            closest_color = curve.color.sample(param);
                                        }
                                        closest_condition = curve.boundaryCondition;
                                    }
                                }
                                // west
                                if (intersect_segment_segment(coord, coord - Eigen::Vector2d(spacing.x(), 0), p0, p1, t))
                                {
                                    if (t.x() < closest_t)
                                    {
                                        closest_t = t.x();
                                        if (curve.boundaryCondition == EBoundaryCondition::Dirichlet)
                                        {
                                            double param  = curve.position.parameters.getValue(ivert).x() * (1 - t.y()) + curve.position.parameters.getValue(ivert + 1).x() * t.y();
                                            closest_color = curve.color.sample(param);
                                        }
                                        closest_condition = curve.boundaryCondition;
                                    }
                                }
                                // south
                                if (intersect_segment_segment(coord, coord - Eigen::Vector2d(0, spacing.y()), p0, p1, t))
                                {
                                    if (t.x() < closest_t)
                                    {
                                        closest_t = t.x();
                                        if (curve.boundaryCondition == EBoundaryCondition::Dirichlet)
                                        {
                                            double param  = curve.position.parameters.getValue(ivert).x() * (1 - t.y()) + curve.position.parameters.getValue(ivert + 1).x() * t.y();
                                            closest_color = curve.color.sample(param);
                                        }
                                        closest_condition = curve.boundaryCondition;
                                    }
                                }
                            }
                        }
                    }
                }

                if (patch->contains(coord, 1E-8))
                {
                    // if there is a gradient mesh, then sample its Laplacian.
                    if (!patch->gradientMeshes.empty())
                    {
                        bool open_x0 = jacobi.getOpenWest(gridCoord);
                        bool open_x1 = jacobi.getOpenEast(gridCoord);
                        bool open_y0 = jacobi.getOpenSouth(gridCoord);
                        bool open_y1 = jacobi.getOpenNorth(gridCoord);

                        const int maxDepth   = 64;
                        const double epsilon = 1E-10;

                        switch (scene.blendOperation)
                        {
                        case Scene::EBlendOperation::Zero:
                        {
                            // only sample gradient mesh when there is exactly one
                            if (patch->gradientMeshes.size() == 1)
                            {
                                auto [valid, initialColor, estimatedLaplacian] = sampleGradientMeshLaplacian(*patch->gradientMeshes[0], coord, maxDepth, epsilon, spacing, open_x0, open_x1, open_y0, open_y1);
                                if (valid)
                                {
                                    jacobi.setField(gridCoord, initialColor);
                                    jacobi.setSource(gridCoord, estimatedLaplacian);
                                }
                            }
                            break;
                        }
                        case Scene::EBlendOperation::Sum:
                        {
                            std::vector<Eigen::Vector3d> colorLaplacians;
                            Eigen::Vector3d source(0, 0, 0);
                            Eigen::Vector3d initial(0, 0, 0);
                            int counter = 0;
                            for (auto gm : patch->gradientMeshes)
                            {
                                auto [valid, initialColor, estimatedLaplacian] = sampleGradientMeshLaplacian(*gm, coord, maxDepth, epsilon, spacing, open_x0, open_x1, open_y0, open_y1);
                                if (valid)
                                {
                                    source += estimatedLaplacian;
                                    initial += initialColor;
                                    counter += 1;
                                }
                            }
                            if (counter > 0)
                            {
                                jacobi.setField(gridCoord, initial / counter);
                                jacobi.setSource(gridCoord, source);
                            }
                            break;
                        }
                        case Scene::EBlendOperation::Average:
                        {
                            std::vector<Eigen::Vector3d> colorLaplacians;
                            Eigen::Vector3d source(0, 0, 0);
                            Eigen::Vector3d initial(0, 0, 0);
                            int counter = 0;
                            for (auto gm : patch->gradientMeshes)
                            {
                                auto [valid, initialColor, estimatedLaplacian] = sampleGradientMeshLaplacian(*gm, coord, maxDepth, epsilon, spacing, open_x0, open_x1, open_y0, open_y1);
                                if (valid)
                                {
                                    source += estimatedLaplacian;
                                    initial += initialColor;
                                    counter++;
                                }
                            }
                            if (counter > 0)
                            {
                                jacobi.setField(gridCoord, initial / counter);
                                jacobi.setSource(gridCoord, source / counter);
                            }
                            break;
                        }
                        case Scene::EBlendOperation::First:
                        {
                            // look for the gradient mesh in the order in which they were specified in the input
                            for (auto& sceneGM : scene.gradientMeshes)
                            {
                                // if the gradient mesh is in the patch, then sample it
                                if (std::find(patch->gradientMeshes.begin(), patch->gradientMeshes.end(), sceneGM) != patch->gradientMeshes.end())
                                {
                                    auto [valid, initialColor, estimatedLaplacian] = sampleGradientMeshLaplacian(*sceneGM, coord, maxDepth, epsilon, spacing, open_x0, open_x1, open_y0, open_y1);
                                    if (valid)
                                    {
                                        jacobi.setField(gridCoord, initialColor);
                                        jacobi.setSource(gridCoord, estimatedLaplacian);
                                        break;
                                    }
                                }
                            }
                            break;
                        }
                        }
                    }
                }
            }

            if (closest_t != std::numeric_limits<double>::max())
            {
                if (closest_condition == EBoundaryCondition::Dirichlet)
                {
                    jacobi.setMask(gridCoord, JacobiIteration3d::EMask::Dirichlet);
                    jacobi.setSource(gridCoord, closest_color);
                }
                else if (closest_condition == EBoundaryCondition::Neumann)
                    jacobi.setMask(gridCoord, JacobiIteration3d::EMask::Neumann);
                else
                    jacobi.setMask(gridCoord, JacobiIteration3d::EMask::Unknown);
            }
            else
                jacobi.setMask(gridCoord, JacobiIteration3d::EMask::Unknown);
        }

        // a sparse matrix mask to store if a pixel is left or right
        Eigen::SparseMatrix<double> poissonMask(grid->getResolution().y(), grid->getResolution().x());
        std::vector<Eigen::Vector2i> left_curve_ids, right_curve_ids;
        std::vector<Eigen::Vector3d> left_curve_colors, right_curve_colors;

        // Iterate each pixel
        for (Eigen::Index linearIndex = 0; linearIndex < numPixels; ++linearIndex)
        {
            Eigen::Vector2i gridCoord = grid->getGridCoord(linearIndex);
            Eigen::Vector2d pos       = grid->getCoordAt(gridCoord);

            JacobiIteration3d::EMask maskID = jacobi.getMask(gridCoord);
            if (maskID == JacobiIteration3d::EMask::Dirichlet)
                continue;

            // try to find boundary
            Eigen::AlignedBox2d pixelBounds(pos - spacing, pos + spacing);
            double pixelRadius         = pixelBounds.diagonal().norm() * 0.5;
            double prevFoundDist       = std::numeric_limits<double>::max();
            Eigen::Vector3d foundColor = Eigen::Vector3d(1, 0, 0);
            int isLeft                 = -1;

            for (auto& poissonCurve : scene.patches[0]->poissonCurves)
            {
                int numVert = poissonCurve->position.values.getSize();
                for (int ivert = 0; ivert < numVert - 1; ++ivert)
                {
                    const Eigen::Vector2d& p0 = poissonCurve->position.values.getValue(ivert);
                    const Eigen::Vector2d& p1 = poissonCurve->position.values.getValue(ivert + 1);

                    // east
                    {
                        double minDist = std::numeric_limits<double>::max();
                        double minT    = 0;
                        closest_point_on_segment(p0, p1, pos, minDist, minT);
                        if (minDist < prevFoundDist && minDist < pixelRadius)
                        {
                            double param            = poissonCurve->position.parameters.getValue(ivert).x() * (1 - minT) + poissonCurve->position.parameters.getValue(ivert + 1).x() * minT;
                            Eigen::Vector2d tangent = poissonCurve->position.sample_dt(param);
                            if (tangent.stableNorm() < 1E-10)
                                continue;
                            Eigen::Vector2d normal(tangent.y(), -tangent.x());
                            isLeft        = normal.dot(p0 - pos) > 0;
                            foundColor    = poissonCurve->weights.sample(param);
                            prevFoundDist = minDist;
                        }
                    }
                    // west
                    {
                        double minDist = std::numeric_limits<double>::max();
                        double minT    = 0;
                        closest_point_on_segment(p0, p1, pos, minDist, minT);
                        if (minDist < prevFoundDist && minDist < pixelRadius)
                        {
                            double param            = poissonCurve->position.parameters.getValue(ivert).x() * (1 - minT) + poissonCurve->position.parameters.getValue(ivert + 1).x() * minT;
                            Eigen::Vector2d tangent = poissonCurve->position.sample_dt(param);
                            if (tangent.stableNorm() < 1E-10)
                                continue;
                            Eigen::Vector2d normal(tangent.y(), -tangent.x());
                            isLeft        = normal.dot(p0 - pos) > 0;
                            foundColor    = poissonCurve->weights.sample(param);
                            prevFoundDist = minDist;
                        }
                    }
                    // north
                    {
                        double minDist = std::numeric_limits<double>::max();
                        double minT    = 0;
                        closest_point_on_segment(p0, p1, pos, minDist, minT);
                        if (minDist < prevFoundDist && minDist < pixelRadius)
                        {
                            double param            = poissonCurve->position.parameters.getValue(ivert).x() * (1 - minT) + poissonCurve->position.parameters.getValue(ivert + 1).x() * minT;
                            Eigen::Vector2d tangent = poissonCurve->position.sample_dt(param);
                            if (tangent.stableNorm() < 1E-10)
                                continue;
                            Eigen::Vector2d normal(tangent.y(), -tangent.x());
                            isLeft        = normal.dot(p0 - pos) > 0;
                            foundColor    = poissonCurve->weights.sample(param);
                            prevFoundDist = minDist;
                        }
                    }
                    // south
                    {
                        double minDist = std::numeric_limits<double>::max();
                        double minT    = 0;
                        closest_point_on_segment(p0, p1, pos, minDist, minT);
                        if (minDist < prevFoundDist && minDist < pixelRadius)
                        {
                            double param            = poissonCurve->position.parameters.getValue(ivert).x() * (1 - minT) + poissonCurve->position.parameters.getValue(ivert + 1).x() * minT;
                            Eigen::Vector2d tangent = poissonCurve->position.sample_dt(param);
                            if (tangent.stableNorm() < 1E-10)
                                continue;
                            Eigen::Vector2d normal(tangent.y(), -tangent.x());
                            isLeft        = normal.dot(p0 - pos) > 0;
                            foundColor    = poissonCurve->weights.sample(param);
                            prevFoundDist = minDist;
                        }
                    }
                }
            }

            if (prevFoundDist < std::numeric_limits<double>::max())
            {
                // the weights on the Poisson curves are with respect to a certain resolution. we need to rescale them to our target resolution.
                double scaleToTargetResolutionScale =
                    (grid->getResolution().x()) * (grid->getResolution().y()) / double(1024 * 1024);
                foundColor *= scaleToTargetResolutionScale;

                if (isLeft == 1)
                {
                    left_curve_ids.push_back(gridCoord);
                    left_curve_colors.push_back(foundColor);
                    poissonMask.insert(gridCoord.x(), gridCoord.y()) = 1;
                }
                else if (isLeft == 0)
                {
                    right_curve_ids.push_back(gridCoord);
                    right_curve_colors.push_back(foundColor);
                    poissonMask.insert(gridCoord.x(), gridCoord.y()) = -1;
                }
            }
        }

        Eigen::Vector3d sumLeftColors(0, 0, 0);
        Eigen::Vector3d sumRightColors(0, 0, 0);

        poissonMask.makeCompressed();
        // Loop each found position
        for (int i = 0; i < left_curve_ids.size(); i++)
        {
            // get pixel position
            Eigen::Vector2i pixelId = left_curve_ids[i];

            // check the four neighbours
            int weight = 0;
            int px     = pixelId.x();
            int py     = pixelId.y();

            // left
            if (px > 0)
            {
                if (poissonMask.coeff(px - 1, py) == -1)
                {
                    weight++;
                }
            }

            // right
            if (px < grid->getResolution().x() - 1)
            {
                if (poissonMask.coeff(px + 1, py) == -1)
                {
                    weight++;
                }
            }

            // top
            if (py > 0)
            {
                if (poissonMask.coeff(px, py - 1) == -1)
                {
                    weight++;
                }
            }

            // bottom
            if (py < grid->getResolution().y() - 1)
            {
                if (poissonMask.coeff(px, py + 1) == -1)
                {
                    weight++;
                }
            }

            if (weight > 0)
            {
                Eigen::Vector3d weightedColor = left_curve_colors[i] * weight;
                jacobi.setSource(pixelId, jacobi.getSource(pixelId) + weightedColor);
                sumLeftColors += weightedColor;
            }
        }

        std::vector<Eigen::Vector2i> weightedRightPositions;
        std::vector<Eigen::Vector3d> weightedRightColors;
        std::vector<int> weightsRight;
        // check for each right pixel
        for (int i = 0; i < right_curve_ids.size(); i++)
        {
            // get pixel position
            Eigen::Vector2i pixelId = right_curve_ids[i];

            // check the four neighbours
            int weight = 0;
            int px     = pixelId.x();
            int py     = pixelId.y();

            // left
            if (px > 0)
            {
                if (poissonMask.coeff(px - 1, py) == 1)
                {
                    weight++;
                }
            }

            // right
            if (px < grid->getResolution().x() - 1)
            {
                if (poissonMask.coeff(px + 1, py) == 1)
                {
                    weight++;
                }
            }

            // top
            if (py > 0)
            {
                if (poissonMask.coeff(px, py - 1) == 1)
                {
                    weight++;
                }
            }

            // bottom
            if (py < grid->getResolution().y() - 1)
            {
                if (poissonMask.coeff(px, py + 1) == 1)
                {
                    weight++;
                }
            }

            if (weight > 0)
            {
                weightedRightPositions.push_back(pixelId);
                Eigen::Vector3d weightedColor = right_curve_colors[i] * weight;
                weightedRightColors.push_back(weightedColor);
                sumRightColors += weightedColor;
            }
        }

        // Normalize right curve
        Eigen::Vector3d normalizationFactor = -sumLeftColors.array() / sumRightColors.array();
        for (int i = 0; i < weightedRightPositions.size(); i++)
        {
            Eigen::Vector2i pixelId         = weightedRightPositions[i];
            Eigen::Vector3d foundColor      = weightedRightColors[i];
            Eigen::Vector3d normalizedColor = foundColor.array() * normalizationFactor.array();
            jacobi.setSource(pixelId, jacobi.getSource(pixelId) + normalizedColor);
        }
    }

    std::shared_ptr<Image> Scene::solvePDE(const Eigen::Vector2i& resolution, int numIterations, bool useMultigrid) const
    {
        // get the output image
        auto image = std::make_shared<Image>();
        image->setResolution(resolution);

        if (useMultigrid)
        {
            // compute the number of mipmap levels
            int num_mip_levels = 1;
            int max_res        = std::max(image->getResolution().x(), image->getResolution().y());
            while (max_res > 1)
            {
                num_mip_levels++;
                max_res /= 2;
            }

            // allocate solver for each level
            std::vector<std::shared_ptr<JacobiIteration3d>> solvers;
            solvers.resize(num_mip_levels);
            for (int ilevel = 0; ilevel < num_mip_levels; ++ilevel)
            {
                // define the co-located grid over which we solve
                Eigen::Vector2i resolution(
                    image->getResolution().x() / (1u << ilevel),
                    image->getResolution().y() / (1u << ilevel));
                auto grid = std::make_shared<RegularGrid2d>();
                Eigen::Vector2d hrspacing( // real spacing to make the staggered grids align
                    (this->domain.max().x() - this->domain.min().x()) / image->getResolution().x(),
                    (this->domain.max().y() - this->domain.min().y()) / image->getResolution().y());
                grid->setResolution(resolution);
                grid->setDomain(Eigen::AlignedBox2d( // offset into pixel center
                    this->domain.min() + hrspacing * 0.5,
                    this->domain.max() - hrspacing * 0.5));

                // create and initialize each level
                solvers[ilevel] = std::make_shared<JacobiIteration3d>(grid);
                initializeFromPatches(*this, *solvers[ilevel]);
            }

            // run solver per level and propagate result upwards
            for (int ilevel = num_mip_levels - 2; ilevel >= 0; ilevel--)
            {
                // run solver iterations on current level
                int depth = num_mip_levels - ilevel;
                solvers[ilevel]->iterate(numIterations * depth);

                // pull the solution one level up to the next level
                if (ilevel != 0)
                {
                    auto sourceGrid                  = solvers[ilevel]->getGrid();
                    auto targetGrid                  = solvers[ilevel - 1]->getGrid();
                    Eigen::Vector2i targetResolution = targetGrid->getResolution();
                    for (int ix = 0; ix < targetResolution.x(); ++ix)
                    {
                        for (int iy = 0; iy < targetResolution.y(); ++iy)
                        {
                            // if this pixel is not pinned down by a Dirichlet condition
                            if (solvers[ilevel - 1]->getMask(Eigen::Vector2i(ix, iy)) == JacobiIteration3d::EMask::Unknown)
                            {
                                // get the value from the lower level
                                Eigen::Index srcLinearIndex = sourceGrid->getLinearIndex(Eigen::Vector2i(ix / 2, iy / 2));
                                Eigen::Vector3d srcColor    = solvers[ilevel]->getField()[srcLinearIndex];

                                // set it at the higher level
                                Eigen::Index targetLinearIndex = targetGrid->getLinearIndex(Eigen::Vector2i(ix, iy));
                                solvers[ilevel - 1]->setField(Eigen::Vector2i(ix, iy), srcColor);
                            }
                        }
                    }
                }
            }

            // copy result to output
            Eigen::Index numPixels = image->getResolution().prod();
#ifdef NDEBUG
#pragma omp parallel for
#endif
            for (Eigen::Index linearIndex = 0; linearIndex < numPixels; ++linearIndex)
            {
                image->setValue(linearIndex, solvers.front()->getField()[linearIndex]);
            }
        }
        else // simple Jacobi iteration
        {
            // define the co-located grid over which we solve
            auto grid = std::make_shared<RegularGrid2d>();
            Eigen::Vector2d spacing(
                (this->domain.max().x() - this->domain.min().x()) / image->getResolution().x(),
                (this->domain.max().y() - this->domain.min().y()) / image->getResolution().y());
            grid->setResolution(image->getResolution());
            grid->setDomain(Eigen::AlignedBox2d( // offset into pixel center
                this->domain.min() + spacing * 0.5,
                this->domain.max() - spacing * 0.5));
            Eigen::Index numPixels = grid->getResolution().prod();

            // create Jacobi solver
            JacobiIteration3d jacobi(grid);

            initializeFromPatches(*this, jacobi);

            // solve
            jacobi.iterate(numIterations);

// copy result to output
#ifdef NDEBUG
#pragma omp parallel for
#endif
            for (Eigen::Index linearIndex = 0; linearIndex < numPixels; ++linearIndex)
            {
                image->setValue(linearIndex, jacobi.getField()[linearIndex]);
            }
        }
        return image;
    }
}
