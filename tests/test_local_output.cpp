#include "test_local_output.hpp"

/**
 * localOutput writes, beside every trajectory frame, what the model knows
 * locally: each face's bending and mesh-term energy, mean curvature and
 * area, and each vertex's force by term, the bead's included. The files
 * are only worth anything if their rows are the mesh's own numbers in the
 * mesh's own order, so that is what is checked: shape, and agreement with
 * the Face and Vertex objects and with the recorded totals.
 *
 * @see docs/bead_push_benchmark.md section 8
 */

namespace
{

void configure_clamped_sheet(Param &param)
{
    param.VERBOSE_MODE = false;
    param.sideX = 60.0;
    param.sideY = 60.0;
    param.lFace = 5.0;
    param.kCurv = 83.4;
    param.KBT = 4.17;
    param.uSurf = 250.0;
    param.uVol = 0.0;
    param.isGlobalConstraint = true;
    param.usingRpi = false;
    param.isEnergyHarmonicBondIncluded = false;
    param.randomSeed = 20260923u;
    param.timeStep = 1.0e-3;
    param.diffConst = 1.0;
    param.surfaceSolver = "iterative";
    param.boundaryCondition = BoundaryType::Mixed;
    param.boundaryConditionX = BoundaryType::Fixed;
    param.boundaryConditionY = BoundaryType::Fixed;
    param.fixedBoundaryRings = 2;
    param.localOutput = true;
    param.beadEnabled = true;
    param.beadRadius = 8.0;
    param.beadStiffness = 150.0;
    param.beadCenterZ = 10.0 - 8.0 + 3.0; // the sheet cuts the sphere
}

std::vector<double> read_last_line(const std::string &path, int &nLines)
{
    std::ifstream in(path);
    std::vector<double> last;
    std::string line;
    nLines = 0;
    while (std::getline(in, line))
    {
        if (line.empty())
        {
            continue;
        }
        nLines++;
        last.clear();
        std::stringstream stream(line);
        std::string item;
        while (std::getline(stream, item, ','))
        {
            if (!item.empty())
            {
                last.push_back(std::stod(item));
            }
        }
    }
    return last;
}

double close_to(double value)
{
    return 1e-6 * (1.0 + std::fabs(value));
}

} // namespace

TEST(LocalOutputTest, TheParameterKeyParses)
{
    Param param;
    EXPECT_FALSE(param.localOutput);
    EXPECT_TRUE(import_kv_string("localOutput", "true", param));
    EXPECT_TRUE(param.localOutput);
    EXPECT_TRUE(import_kv_string("localOutput", "false", param));
    EXPECT_FALSE(param.localOutput);
}

TEST(LocalOutputTest, TheFramesAreTheMeshesOwnNumbersInOrder)
{
    Param param;
    configure_clamped_sheet(param);
    DynamicMesh mesh(param);
    mesh.setup_flat();
    for (Vertex &vertex : mesh.vertices)
    {
        // A relief, so the bending energies and forces are not all zero.
        const double x = vertex.coord(0, 0);
        const double y = vertex.coord(1, 0);
        vertex.coord.set(2, 0, 10.0 + 1.5 * std::sin(0.2 * x) * std::cos(0.15 * y));
    }
    mesh.calculate_element_area_volume();
    mesh.sum_membrane_area_and_volume(mesh.param.area0, mesh.param.vol0);
    mesh.param.vol0 = 0.0;
    mesh.update_previous_coord_for_vertex();
    mesh.update_reference_coord_from_previous_coord();
    mesh.Compute_Energy_And_Force();
    mesh.update_previous_force_for_vertex();
    mesh.update_vertices_mat_with_vector();
    mesh.apply_mesh_to_surface();
    Record record(16);
    record.add(mesh.param.area, mesh.param.energy, mesh.calculate_mean_force());
    DynamicModel model(mesh, record);
    model.update_bead(0);
    ASSERT_GT(model.beadContacts, 0);

    const std::string prefix = "test_local_output_tmp";
    const std::string facePath = "faceenergy" + prefix + ".csv";
    const std::string vertexPath = "vertexforce" + prefix + ".csv";
    std::remove(facePath.c_str());
    std::remove(vertexPath.c_str());
    dynamics_create_local_files(mesh, model, prefix);
    dynamics_output_local_frame(mesh, model, prefix);

    int nFaceLines = 0;
    int nVertexLines = 0;
    const std::vector<double> faceRow = read_last_line(facePath, nFaceLines);
    const std::vector<double> vertexRow = read_last_line(vertexPath, nVertexLines);
    std::remove(facePath.c_str());
    std::remove(vertexPath.c_str());

    // Two frames: the one the files were created with, and the one appended.
    EXPECT_EQ(nFaceLines, 2);
    EXPECT_EQ(nVertexLines, 2);
    ASSERT_EQ(faceRow.size(), 4u * mesh.faces.size());
    ASSERT_EQ(vertexRow.size(), 12u * mesh.vertices.size());

    double sumCurvature = 0.0;
    double sumMesh = 0.0;
    for (std::size_t f = 0; f < mesh.faces.size(); f++)
    {
        const Face &face = mesh.faces[f];
        EXPECT_NEAR(faceRow[4 * f + 0], face.energy.energyCurvature, close_to(face.energy.energyCurvature)) << f;
        EXPECT_NEAR(faceRow[4 * f + 1], face.energy.energyRegularization,
                    close_to(face.energy.energyRegularization)) << f;
        EXPECT_NEAR(faceRow[4 * f + 2], face.meanCurvature, close_to(face.meanCurvature)) << f;
        EXPECT_NEAR(faceRow[4 * f + 3], face.elementArea, close_to(face.elementArea)) << f;
        sumCurvature += faceRow[4 * f + 0];
        sumMesh += faceRow[4 * f + 1];
    }
    EXPECT_GT(sumCurvature, 0.0);
    EXPECT_NEAR(sumCurvature, mesh.param.energy.energyCurvature, 1e-5 * (1.0 + mesh.param.energy.energyCurvature));
    EXPECT_NEAR(sumMesh, mesh.param.energy.energyRegularization,
                1e-5 * (1.0 + std::fabs(mesh.param.energy.energyRegularization)));

    int beadRows = 0;
    for (std::size_t v = 0; v < mesh.vertices.size(); v++)
    {
        const Force &force = mesh.vertices[v].force;
        for (int axis = 0; axis < 3; axis++)
        {
            EXPECT_NEAR(vertexRow[12 * v + axis], force.forceCurvature(axis, 0), close_to(force.forceCurvature(axis, 0))) << v;
            EXPECT_NEAR(vertexRow[12 * v + 3 + axis], force.forceArea(axis, 0), close_to(force.forceArea(axis, 0))) << v;
            EXPECT_NEAR(vertexRow[12 * v + 6 + axis], force.forceRegularization(axis, 0),
                        close_to(force.forceRegularization(axis, 0))) << v;
            EXPECT_NEAR(vertexRow[12 * v + 9 + axis], model.beadForce(static_cast<int>(v), axis),
                        close_to(model.beadForce(static_cast<int>(v), axis))) << v;
        }
        if (model.beadForce(static_cast<int>(v), 2) != 0.0)
        {
            beadRows++;
        }
    }
    EXPECT_EQ(beadRows, model.beadContacts);
}
