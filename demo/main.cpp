#include <usvg/image.hpp>
#include <usvg/scene.hpp>

/**
 * @brief Reads and renders a unified smooth vector graphics scene.
 * @param pathIn Path to the XML file that contains the scene.
 * @param pathOut Output where the resulting PDF is written to.
 * @param resolution Resolution of the image to render.
 * @param iterations Number of Jacobi iterations.
 * @param vertexMergeThreshold Threshold for snapping of nearby vertices. The parameter is called 'tau' in the paper.
 * @param discretizationResidual Accepted approximation error for the discretization of the curves into polylines. The parameter is called 'epsilon' in the paper.
 * @param useMultigrid Flag that determines whether a multi-grid Jacobi solver is used. If not, Jacobi relaxation on the target image is used.
 */
void process(const std::string& pathIn, const std::string& pathOut, const Eigen::Vector2i& resolution, int iterations, double vertexMergeThreshold, double discretizationResidual, bool useMultigrid)
{
    using namespace usvg;

    // print current scene path
    std::cout << "Computing " << pathOut << std::endl;

    // read scene
    Scene scene;
    if (!scene.readXML(pathIn))
        return;

    // run the patch builder
    scene.buildPatches(vertexMergeThreshold, 1E-8, discretizationResidual);

    // solve scene
    std::shared_ptr<Image> image = scene.solvePDE(resolution, iterations, useMultigrid);

    // write image
    image->writeBMP(pathOut);
}

int main()
{
    process(std::string(USVG_SCENES_PATH) + "crane.xml", "crane.bmp", Eigen::Vector2i(512, 512), 100, 1, 0.1, true);
    process(std::string(USVG_SCENES_PATH) + "crane_ribbons.xml", "crane_ribbons.bmp", Eigen::Vector2i(512, 512), 100, 1, 0.1, true);
    process(std::string(USVG_SCENES_PATH) + "bubble.xml", "bubble.bmp", Eigen::Vector2i(512, 512), 100, 1, 0.1, true);
    process(std::string(USVG_SCENES_PATH) + "drape.xml", "drape.bmp", Eigen::Vector2i(512, 512), 100, 0.01, 0.1, true);
    process(std::string(USVG_SCENES_PATH) + "ladybug.xml", "ladybug.bmp", Eigen::Vector2i(512, 512), 100, 1, 0.1, true);
    process(std::string(USVG_SCENES_PATH) + "pepper.xml", "pepper.bmp", Eigen::Vector2i(512, 512), 100, 1, 0.1, true);
    process(std::string(USVG_SCENES_PATH) + "portal.xml", "portal.bmp", Eigen::Vector2i(512, 512), 100, 0.01, 0.1, true);
    process(std::string(USVG_SCENES_PATH) + "sunset.xml", "sunset.bmp", Eigen::Vector2i(256, 256), 50000, 1, 0.1, false);

    return 0;
}
