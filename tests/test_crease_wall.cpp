#include "test_crease_wall.hpp"

/**
 * A face folded flat onto its neighbour is the failure the altitude floor
 * could not stop, and this is the wall against it. The checks: zero on the
 * lattice and on any thermal crease; a hand-made fold charged what the
 * formula says; the force minus the gradient with the wall active on every
 * edge; the flip trial and the force pass reaching one expression; and the
 * configuration guard.
 *
 * @see docs/edge_flip_plan.md work package 7
 */

namespace
{

void configure_sheet(Param &param, double side = 60.0)
{
    param.VERBOSE_MODE = false;
    param.boundaryCondition = BoundaryType::Periodic;
    param.sideX = side;
    param.sideY = side;
    param.lFace = 5.0;
    param.kCurv = 83.4;
    param.KBT = 4.17;
    param.uSurf = 250.0;
    param.uVol = 0.0;
    param.isGlobalConstraint = true;
    param.usingRpi = false;
    param.isEnergyHarmonicBondIncluded = false;
    param.randomSeed = 20260911u;
    param.timeStep = 1.0e-3;
    param.diffConst = 1.0;
    param.surfaceSolver = "iterative";
    param.edgeSpringEnabled = true;
    param.edgeSpringConstant = 83.4;
    param.edgeTetherMinRatio = 0.01; // the tether's walls out of reach: only the
    param.edgeTetherMaxRatio = 100.0; // crease wall is in the regularization slot
    param.creaseWallEnabled = true;
}

void build_sheet(DynamicMesh &mesh)
{
    mesh.setup_flat();
    for (Vertex &vertex : mesh.vertices)
    {
        vertex.coord.set(2, 0, 10.0);
    }
    mesh.calculate_element_area_volume();
    mesh.sum_membrane_area_and_volume(mesh.param.area0, mesh.param.vol0);
    mesh.update_previous_coord_for_vertex();
    mesh.update_reference_coord_from_previous_coord();
    mesh.Compute_Energy_And_Force();
    mesh.update_vertices_mat_with_vector();
}

double total_regularization(const Mesh &mesh)
{
    double sum = 0.0;
    for (const Face &face : mesh.faces)
    {
        sum += face.energy.energyRegularization;
    }
    return sum;
}

double total_energy(Mesh &mesh)
{
    mesh.Compute_Energy_And_Force();
    return mesh.param.energy.energyTotal;
}

/// Cosine of the angle between the normals of an edge's two faces.
double edge_cosine(const Mesh &mesh, const MeshEdge &edge)
{
    double n[2][3];
    for (int side = 0; side < 2; side++)
    {
        const std::vector<int> &c = mesh.faces[edge.face[side]].adjacentVertices;
        double d1[3], d2[3];
        for (int axis = 0; axis < 3; axis++)
        {
            d1[axis] = mesh.vertices[c[1]].coord.get(axis, 0) - mesh.vertices[c[0]].coord.get(axis, 0);
            d2[axis] = mesh.vertices[c[2]].coord.get(axis, 0) - mesh.vertices[c[0]].coord.get(axis, 0);
        }
        n[side][0] = d1[1] * d2[2] - d1[2] * d2[1];
        n[side][1] = d1[2] * d2[0] - d1[0] * d2[2];
        n[side][2] = d1[0] * d2[1] - d1[1] * d2[0];
        const double len = std::sqrt(n[side][0] * n[side][0] + n[side][1] * n[side][1] + n[side][2] * n[side][2]);
        for (int axis = 0; axis < 3; axis++)
        {
            n[side][axis] /= len;
        }
    }
    return n[0][0] * n[1][0] + n[0][1] * n[1][1] + n[0][2] * n[1][2];
}

int a_deep_interior_edge(const Mesh &mesh)
{
    for (const MeshEdge &edge : mesh.edges)
    {
        if (edge.face[0] < 0 || edge.face[1] < 0)
        {
            continue;
        }
        bool deep = true;
        for (int side = 0; side < 2; side++)
        {
            for (int v : mesh.faces[edge.face[side]].adjacentVertices)
            {
                deep = deep && !mesh.vertices[v].isGhost &&
                       (mesh.flipFrozenVertex.empty() || !mesh.flipFrozenVertex[v]) &&
                       std::abs(mesh.vertices[v].coord.get(0, 0)) < 8.0 &&
                       std::abs(mesh.vertices[v].coord.get(1, 0)) < 8.0;
            }
        }
        if (deep)
        {
            return edge.index;
        }
    }
    return -1;
}

int first_flippable_edge(const Mesh &mesh)
{
    for (const MeshEdge &edge : mesh.edges)
    {
        if (edge.flippable)
        {
            return edge.index;
        }
    }
    return -1;
}

} // namespace

TEST(CreaseWallTest, TheLatticeAndAThermalCreaseCarryNoEnergy)
{
    Param param;
    configure_sheet(param);
    DynamicMesh mesh(param);
    ASSERT_NO_FATAL_FAILURE(build_sheet(mesh));
    EXPECT_EQ(total_regularization(mesh), 0.0) << "the flat lattice carries crease energy";

    // A gentle undulation: adjacent normals a few degrees apart, far from 60.
    for (Vertex &vertex : mesh.vertices)
    {
        vertex.coord.set(2, 0, 10.0 + 1.0 * std::sin(0.2 * vertex.coord.get(0, 0)));
    }
    mesh.energy_force_edge_spring();
    mesh.energy_force_crease_wall();
    EXPECT_EQ(total_regularization(mesh), 0.0) << "a thermal-scale crease was charged";
}

TEST(CreaseWallTest, AFoldIsChargedWhatTheFormulaSays)
{
    Param param;
    configure_sheet(param);
    DynamicMesh mesh(param);
    ASSERT_NO_FATAL_FAILURE(build_sheet(mesh));

    const int iEdge = a_deep_interior_edge(mesh);
    ASSERT_GE(iEdge, 0);
    const MeshEdge &edge = mesh.edges[iEdge];
    // Fold face 1 onto face 0: reflect face 1's opposite corner through the
    // edge's line, in the sheet's plane, so the two faces overlap exactly.
    const int a = edge.v[0], b = edge.v[1], t = edge.opposite[1];
    double pa[3], pb[3], pt[3];
    for (int axis = 0; axis < 3; axis++)
    {
        pa[axis] = mesh.vertices[a].coord.get(axis, 0);
        pb[axis] = mesh.vertices[b].coord.get(axis, 0);
        pt[axis] = mesh.vertices[t].coord.get(axis, 0);
    }
    double e[3], w[3];
    double ee = 0.0, we = 0.0;
    for (int axis = 0; axis < 3; axis++)
    {
        e[axis] = pb[axis] - pa[axis];
        w[axis] = pt[axis] - pa[axis];
        ee += e[axis] * e[axis];
        we += w[axis] * e[axis];
    }
    for (int axis = 0; axis < 3; axis++)
    {
        const double along = we / ee * e[axis];
        const double across = w[axis] - along;
        mesh.vertices[t].coord.set(axis, 0, pa[axis] + along - across); // reflected
    }
    // Lift it a hair so the two faces are not coplanar-and-inverted but a fold
    // of very nearly 180 degrees, which is what the run produced.
    mesh.vertices[t].coord.set(2, 0, mesh.vertices[t].coord.get(2, 0) + 0.05);

    const double c = edge_cosine(mesh, edge);
    ASSERT_LT(c, -0.9) << "the fold was not built: cosine " << c;
    const double c0 = std::cos(param.creaseWallAngle * M_PI / 180.0);
    const double expectedEdge = 0.5 * param.creaseWallConstant * (c0 - c) * (c0 - c);

    // Both faces of the folded edge get half; the other edges of the moved
    // corner fold too, so compare the whole slot against the hand formula.
    mesh.energy_force_edge_spring();
    mesh.energy_force_crease_wall();
    double expectedTotal = 0.0;
    for (const MeshEdge &other : mesh.edges)
    {
        if (other.face[0] < 0 || other.face[1] < 0 || !mesh.edge_carries_tether(other))
        {
            continue;
        }
        const double ci = edge_cosine(mesh, other);
        if (ci < c0)
        {
            expectedTotal += 0.5 * param.creaseWallConstant * (c0 - ci) * (c0 - ci);
        }
    }
    EXPECT_GT(expectedEdge, 100.0) << "a full fold should cost of order the wall constant";
    EXPECT_NEAR(total_regularization(mesh), expectedTotal, 1e-9 * expectedTotal);
    // The two faces of the folded edge carry that edge's energy between them,
    // plus half of whatever their other edges owe.
    EXPECT_GE(mesh.face_crease_energy(edge.face[0]) + mesh.face_crease_energy(edge.face[1]),
              expectedEdge - 1e-9 * expectedEdge);
    double sumOverFaces = 0.0;
    for (int iFace = 0; iFace < static_cast<int>(mesh.faces.size()); iFace++)
    {
        sumOverFaces += mesh.face_crease_energy(iFace);
    }
    EXPECT_NEAR(sumOverFaces, expectedTotal, 1e-9 * expectedTotal)
        << "the per-face split does not add back up to the sum over edges";
}

TEST(CreaseWallTest, TheForceIsMinusTheGradientOfTheEnergy)
{
    Param param;
    configure_sheet(param);
    param.kCurv = 0.0;
    param.uSurf = 0.0;
    param.isGlobalConstraint = false;
    param.creaseWallAngle = 2.0; // active on nearly every edge of a rough sheet
    DynamicMesh mesh(param);
    ASSERT_NO_FATAL_FAILURE(build_sheet(mesh));

    std::uint64_t state = 0x9E3779B97F4A7C15ULL;
    const auto next = [&state]() {
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        return static_cast<double>((state >> 11) * (1.0 / 9007199254740992.0)) - 0.5;
    };
    for (Vertex &vertex : mesh.vertices)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            vertex.coord.set(axis, 0, vertex.coord.get(axis, 0) + 1.2 * next());
        }
    }
    mesh.energy_force_edge_spring();
    mesh.energy_force_crease_wall();
    ASSERT_GT(total_regularization(mesh), 0.0) << "no edge is on the wall";

    const double h = 1e-6;
    int checked = 0;
    for (int iVertex = 30; iVertex < static_cast<int>(mesh.vertices.size()) && checked < 15;
         iVertex += 7)
    {
        if (mesh.vertices[iVertex].isGhost)
        {
            continue;
        }
        for (int axis = 0; axis < 3; axis++)
        {
            const double original = mesh.vertices[iVertex].coord.get(axis, 0);
            mesh.vertices[iVertex].coord.set(axis, 0, original + h);
            mesh.energy_force_edge_spring();
            mesh.energy_force_crease_wall();
            const double plus = total_regularization(mesh);
            mesh.vertices[iVertex].coord.set(axis, 0, original - h);
            mesh.energy_force_edge_spring();
            mesh.energy_force_crease_wall();
            const double minus = total_regularization(mesh);
            mesh.vertices[iVertex].coord.set(axis, 0, original);
            mesh.energy_force_edge_spring();
            mesh.energy_force_crease_wall();
            const double analytic = mesh.vertices[iVertex].force.forceRegularization.get(axis, 0);
            const double numeric = -(plus - minus) / (2.0 * h);
            EXPECT_NEAR(analytic, numeric, 1e-4 * std::max(1.0, std::abs(numeric)))
                << "vertex " << iVertex << " axis " << axis;
            checked++;
        }
    }
    EXPECT_GE(checked, 12);
}

TEST(CreaseWallTest, TheTrialAndTheDynamicsShareOneCreaseTerm)
{
    Param param;
    configure_sheet(param);
    param.creaseWallAngle = 2.0;
    DynamicMesh mesh(param);
    ASSERT_NO_FATAL_FAILURE(build_sheet(mesh));
    for (Vertex &vertex : mesh.vertices)
    {
        vertex.coord.set(2, 0, 10.0 + 0.6 * std::sin(0.4 * vertex.coord.get(0, 0)) *
                                          std::cos(0.3 * vertex.coord.get(1, 0)));
    }
    const double fromForcePass = total_energy(mesh) - 0.0;
    double slot = total_regularization(mesh);
    ASSERT_GT(slot, 0.0);

    std::vector<int> everyFace(mesh.faces.size());
    for (int iFace = 0; iFace < static_cast<int>(mesh.faces.size()); iFace++)
    {
        everyFace[iFace] = iFace;
    }
    const FaceSubsetEnergy subset = mesh.evaluate_face_subset(everyFace);
    EXPECT_NEAR(subset.regularization, slot, 1e-9 * slot)
        << "the trial's evaluator returned " << subset.regularization << " where the force pass computes " << slot;

    const int iEdge = first_flippable_edge(mesh);
    ASSERT_GE(iEdge, 0);
    const double before = total_energy(mesh);
    EdgeFlipDelta delta;
    ASSERT_TRUE(mesh.evaluate_edge_flip(iEdge, delta));
    EXPECT_NEAR(total_energy(mesh), before, 1e-9 * std::max(1.0, std::abs(before)));
    mesh.flip_edge(iEdge);
    mesh.param.area += delta.area;
    mesh.param.vol += delta.volume;
    const double after = total_energy(mesh);
    EXPECT_NEAR(delta.energy, after - before, 1e-6 * std::max(1.0, std::abs(delta.energy)))
        << "the sweep's dE was " << delta.energy << " but the mesh's energy moved by " << (after - before);
    (void)fromForcePass;
}

TEST(CreaseWallTest, TheWallWithoutTheTetherIsRefusedAtSetup)
{
    Param param;
    configure_sheet(param);
    param.edgeSpringEnabled = false;
    DynamicMesh mesh(param);
    EXPECT_THROW(mesh.setup_flat(), std::runtime_error);
}

TEST(CreaseWallTest, TheParametersRoundTrip)
{
    Param param;
    EXPECT_FALSE(param.creaseWallEnabled);
    EXPECT_DOUBLE_EQ(param.creaseWallAngle, 60.0);
    EXPECT_TRUE(import_kv_string("creaseWallEnabled", "true", param));
    EXPECT_TRUE(param.creaseWallEnabled);
    EXPECT_TRUE(import_kv_string("creaseWallAngle", "45", param));
    EXPECT_DOUBLE_EQ(param.creaseWallAngle, 45.0);
    EXPECT_TRUE(import_kv_string("creaseWallConstant", "250", param));
    EXPECT_DOUBLE_EQ(param.creaseWallConstant, 250.0);
}
