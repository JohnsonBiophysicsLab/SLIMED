#include "test_irregular_patch_quadrature.hpp"

namespace
{
const double kPi = std::acos(-1.0);

double subdominant(int valence) { return 3.0 / 8.0 + std::cos(2.0 * kPi / valence) / 4.0; }
double saddle_eigenvalue(int valence) { return 3.0 / 8.0 + std::cos(4.0 * kPi / valence) / 4.0; }

/**
 * @brief A two-ring disk around one valence-N vertex: unit spokes, the outer
 * ring at 2 u_k and u_k + u_{k+1}.
 *
 * All N fan faces have complete one-rings, so summing them sums the whole
 * neighbourhood of the extraordinary vertex, where the cross terms between
 * angular modes cancel. `bowl` and `cylinder` set the height.
 */
void build_fan(int valence, double bowl, double cylinder, std::vector<std::vector<double>> &V,
               std::vector<std::vector<int>> &F)
{
    auto height = [&](double x, double y) { return bowl * (x * x + y * y) + cylinder * x * x; };
    V.clear();
    F.clear();
    V.push_back({0.0, 0.0, 0.0});
    std::vector<double> ux(valence), uy(valence);
    for (int k = 0; k < valence; k++)
    {
        ux[k] = std::cos(2.0 * kPi * k / valence);
        uy[k] = std::sin(2.0 * kPi * k / valence);
        V.push_back({ux[k], uy[k], height(ux[k], uy[k])});
    }
    for (int k = 0; k < valence; k++)
    {
        V.push_back({2.0 * ux[k], 2.0 * uy[k], height(2.0 * ux[k], 2.0 * uy[k])});
    }
    for (int k = 0; k < valence; k++)
    {
        const int next = (k + 1) % valence;
        const double x = ux[k] + ux[next];
        const double y = uy[k] + uy[next];
        V.push_back({x, y, height(x, y)});
    }
    for (int k = 0; k < valence; k++)
    {
        const int next = (k + 1) % valence;
        const int a = 1 + k, b = 1 + next;
        const int outerA = 1 + valence + k, outerB = 1 + valence + next;
        const int middle = 1 + 2 * valence + k;
        F.push_back({0, a, b});
        F.push_back({a, outerA, middle});
        F.push_back({a, middle, b});
        F.push_back({b, middle, outerB});
    }
}

Param fan_param()
{
    Param param;
    param.VERBOSE_MODE = false;
    param.boundaryCondition = BoundaryType::Fixed;
    param.uVol = 0.0;
    return param;
}

/// The fan's faces, set up and evaluated once.
struct Fan
{
    explicit Fan(const Param &p, int valence, double bowl, double cylinder) : param(p), mesh(param)
    {
        std::vector<std::vector<double>> V;
        std::vector<std::vector<int>> F;
        build_fan(valence, bowl, cylinder, V, F);
        mesh.setup_from_vertices_faces(V, F);
        for (int i = 0; i < static_cast<int>(mesh.faces.size()); i++)
        {
            if (mesh.faces[i].patchKind == PatchKind::SingleExtraordinary)
            {
                fanFaces.push_back(i);
            }
        }
        mesh.calculate_element_area_volume();
        double area = 0.0;
        double volume = 0.0;
        mesh.sum_membrane_area_and_volume(area, volume);
        mesh.param.area0 = area;
        mesh.param.vol0 = 0.0;
        mesh.update_previous_coord_for_vertex();
        mesh.update_reference_coord_from_previous_coord();
        mesh.Compute_Energy_And_Force();
    }

    double bending() const
    {
        double sum = 0.0;
        for (int i : fanFaces)
        {
            sum += mesh.faces[i].energy.energyCurvature;
        }
        return sum;
    }

    double area() const
    {
        double sum = 0.0;
        for (int i : fanFaces)
        {
            sum += mesh.faces[i].elementArea;
        }
        return sum;
    }

    Param param;
    Mesh mesh;
    std::vector<int> fanFaces;
};

Param with_quadrature(const IrregularPatchDepths &depths, bool tailCorrection)
{
    Param param = fan_param();
    for (int i = 0; i < kExtraordinaryValenceCount; i++)
    {
        param.irregularPatchDepths[i] = depths[i];
    }
    param.irregularPatchTailCorrection = tailCorrection;
    return param;
}

constexpr int kValences[kExtraordinaryValenceCount] = {4, 5, 7, 8};
} // namespace

/// The weights section 9.3 of docs/irregular_patch_results.md tabulates.
TEST(IrregularPatchQuadratureTest, TailWeightsMatchTheDocumentedTable)
{
    // N: bending w_{D-2}, bending w_{D-1}, area w_{D-1}
    const double table[kExtraordinaryValenceCount][3] = {{0.90000, 2.02500, 1.16364},
                                                         {0.92483, 1.68602, 1.25712},
                                                         {0.83837, 2.01378, 1.39242},
                                                         {0.77828, 2.33838, 1.43773}};
    for (int i = 0; i < kExtraordinaryValenceCount; i++)
    {
        const int valence = kValences[i];
        const int depth = 5;
        for (int level = 0; level < depth; level++)
        {
            double bend = 0.0, area = 0.0;
            irregular_tail_weights(valence, depth, level, true, bend, area);
            const double expectedBend =
                level == depth - 2 ? table[i][0] : (level == depth - 1 ? table[i][1] : 1.0);
            const double expectedArea = level == depth - 1 ? table[i][2] : 1.0;
            EXPECT_NEAR(bend, expectedBend, 1e-5) << "N=" << valence << " level " << level;
            EXPECT_NEAR(area, expectedArea, 1e-5) << "N=" << valence << " level " << level;

            irregular_tail_weights(valence, depth, level, false, bend, area);
            EXPECT_EQ(bend, 1.0);
            EXPECT_EQ(area, 1.0);
        }
    }
}

/**
 * @brief The weights sum both geometric bending tails, and the area tail,
 * exactly.
 *
 * That is their definition: whatever the depth, sum_{d<D} w_d r^d must equal
 * sum_{d>=0} r^d for the bowl and the saddle rate alike. They must also stay
 * positive, or the corrected energy could go negative.
 */
TEST(IrregularPatchQuadratureTest, TailWeightsSumTheGeometricTailsExactly)
{
    for (const int valence : kValences)
    {
        const double lambda = subdominant(valence);
        const double mu2 = saddle_eigenvalue(valence);
        const double rates[2] = {1.0 / (16.0 * lambda * lambda), mu2 * mu2 / (lambda * lambda)};
        const double areaRate = lambda * lambda;
        for (int depth = 2; depth <= kMaxIrregularPatchDepth; depth++)
        {
            double bendSums[2] = {0.0, 0.0};
            double areaSum = 0.0;
            for (int level = 0; level < depth; level++)
            {
                double bend = 0.0, area = 0.0;
                irregular_tail_weights(valence, depth, level, true, bend, area);
                EXPECT_GT(bend, 0.0);
                EXPECT_GT(area, 0.0);
                for (int r = 0; r < 2; r++)
                {
                    bendSums[r] += bend * std::pow(rates[r], level);
                }
                areaSum += area * std::pow(areaRate, level);
            }
            for (int r = 0; r < 2; r++)
            {
                EXPECT_NEAR(bendSums[r], 1.0 / (1.0 - rates[r]), 1e-12)
                    << "N=" << valence << " D=" << depth << " rate " << r;
            }
            EXPECT_NEAR(areaSum, 1.0 / (1.0 - areaRate), 1e-12)
                << "N=" << valence << " D=" << depth;
        }
    }
}

/// A Param left alone gets the tail-corrected 6/5/5/5, and a Mesh uses it.
TEST(IrregularPatchQuadratureTest, TheDefaultIsTheTailCorrectedDepths)
{
    Param param = fan_param();
    for (int i = 0; i < kExtraordinaryValenceCount; i++)
    {
        EXPECT_EQ(param.irregularPatchDepths[i], kDefaultIrregularPatchDepths[i]);
    }
    EXPECT_TRUE(param.irregularPatchTailCorrection);

    Mesh mesh(param);
    EXPECT_TRUE(mesh.irregularRows.tail_correction());
    for (int i = 0; i < kExtraordinaryValenceCount; i++)
    {
        EXPECT_EQ(mesh.irregularRows.depth_for(kValences[i]), kDefaultIrregularPatchDepths[i]);
        EXPECT_EQ(recommended_irregular_depth(kValences[i]), kDefaultIrregularPatchDepths[i]);
        EXPECT_EQ(legacy_irregular_depth(kValences[i]), kLegacyIrregularPatchDepths[i]);
    }
}

TEST(IrregularPatchQuadratureTest, TheParametersRoundTripThroughTheParamsFile)
{
    Param param;
    // The form a params file writes after pop_space(), with and without brackets.
    EXPECT_TRUE(import_kv_string("irregularPatchDepths", "12,8,7,7", param));
    EXPECT_EQ(param.irregularPatchDepths[0], 12);
    EXPECT_EQ(param.irregularPatchDepths[1], 8);
    EXPECT_EQ(param.irregularPatchDepths[2], 7);
    EXPECT_EQ(param.irregularPatchDepths[3], 7);
    EXPECT_TRUE(import_kv_string("irregularPatchDepths", "[9,6,5,4]", param));
    EXPECT_EQ(param.irregularPatchDepths[0], 9);
    EXPECT_EQ(param.irregularPatchDepths[3], 4);

    EXPECT_THROW(import_kv_string("irregularPatchDepths", "12,8,7", param), std::runtime_error);
    EXPECT_THROW(import_kv_string("irregularPatchDepths", "12,8,7,7,7", param),
                 std::runtime_error);
    EXPECT_THROW(import_kv_string("irregularPatchDepths", "12,0,7,7", param), std::runtime_error);
    EXPECT_THROW(import_kv_string("irregularPatchDepths", "12,8,7,25", param),
                 std::runtime_error);
    EXPECT_THROW(import_kv_string("irregularPatchDepths", "12,8.5,7,7", param),
                 std::runtime_error);
    // A rejected list leaves the last good one in place.
    EXPECT_EQ(param.irregularPatchDepths[0], 9);

    EXPECT_TRUE(import_kv_string("irregularPatchTailCorrection", "false", param));
    EXPECT_FALSE(param.irregularPatchTailCorrection);
    EXPECT_TRUE(import_kv_string("irregularPatchTailCorrection", "true", param));
    EXPECT_TRUE(param.irregularPatchTailCorrection);
    EXPECT_THROW(import_kv_string("irregularPatchTailCorrection", "True", param),
                 std::runtime_error);

    const IrregularPatchQuadrature quadrature = irregular_patch_quadrature(param);
    EXPECT_EQ(quadrature.depths[0], 9);
    EXPECT_TRUE(quadrature.tailCorrection);
}

/// Depth lists the table cannot build are refused at setup, not at a face.
TEST(IrregularPatchQuadratureTest, AnOutOfRangeDepthIsRejected)
{
    Param param = fan_param();
    IrregularPatchRowTable table;
    IrregularPatchQuadrature quadrature;
    quadrature.depths = {{0, 5, 5, 5}};
    EXPECT_THROW(table.build(param.shapeFunctions, quadrature), std::invalid_argument);
    quadrature.depths = {{6, 5, 5, kMaxIrregularPatchDepth + 1}};
    EXPECT_THROW(table.build(param.shapeFunctions, quadrature), std::invalid_argument);
    quadrature.depths = {{6, 5, 5, 5}};
    quadrature.depthScale = 0.0;
    EXPECT_THROW(table.build(param.shapeFunctions, quadrature), std::invalid_argument);
}

/**
 * @brief `irregularPatchDepths = 12, 8, 7, 7` without the tail correction is
 * the old evaluation, bit for bit.
 *
 * Summed here with the plain kernel over the same rows, the way every run
 * before the correction did it, and compared exactly per face.
 */
TEST(IrregularPatchQuadratureTest, TheLegacyParametersReproduceThePlainTruncatedSum)
{
    for (int i = 0; i < kExtraordinaryValenceCount; i++)
    {
        const int valence = kValences[i];
        Fan fan(with_quadrature(kLegacyIrregularPatchDepths, false), valence, 0.12, 0.05);
        ASSERT_EQ(static_cast<int>(fan.fanFaces.size()), valence);
        const Mesh &mesh = fan.mesh;
        ASSERT_EQ(mesh.irregularRows.depth_for(valence), kLegacyIrregularPatchDepths[i]);

        for (int face : fan.fanFaces)
        {
            const Face &f = mesh.faces[face];
            const int nCtrl = static_cast<int>(f.oneRingVertices.size());
            double coords[slimed::kMaxControlPoints * 3];
            for (int j = 0; j < nCtrl; j++)
            {
                for (int axis = 0; axis < 3; axis++)
                {
                    coords[j * 3 + axis] = mesh.vertices[f.oneRingVertices[j]].coord(axis, 0);
                }
            }
            slimed::PatchParams p;
            p.kCurv = mesh.param.kCurv;
            p.spontCurv = f.spontCurvature;
            p.uSurfPerArea = mesh.param.uSurf / mesh.param.area0;
            p.area = mesh.param.area;
            p.area0 = mesh.param.area0;
            double bending = 0.0, area = 0.0, volume = 0.0;
            double fBend[slimed::kMaxControlPoints * 3] = {0.0};
            double fArea[slimed::kMaxControlPoints * 3] = {0.0};
            double fVolume[slimed::kMaxControlPoints * 3] = {0.0};
            for (int d = 0; d < kLegacyIrregularPatchDepths[i]; d++)
            {
                for (int c = 0; c < kRegularChildrenPerStep; c++)
                {
                    const double *rows = mesh.patchRowsFlat.child(valence, d, c);
                    double childBending = 0.0, meanCurv = 0.0, normal[3] = {0.0, 0.0, 0.0};
                    slimed::element_energy_force_patch_pod(
                        rows, mesh.patchRowsFlat.gaussCoeff(), mesh.patchRowsFlat.nSamples(),
                        coords, nCtrl, p, childBending, meanCurv, normal, fBend, fArea, fVolume);
                    bending += childBending;
                    slimed::element_area_volume_pod(rows, mesh.patchRowsFlat.gaussCoeff(),
                                                    mesh.patchRowsFlat.nSamples(), coords, nCtrl,
                                                    area, volume);
                }
            }
            EXPECT_EQ(f.energy.energyCurvature, bending) << "N=" << valence << " face " << face;
            EXPECT_EQ(f.elementArea, area) << "N=" << valence << " face " << face;
        }
    }
}

/**
 * @brief The default is as accurate as the legacy depths, and the tail
 * correction is what makes it so.
 *
 * Against a deep reference (depth 16 everywhere, corrected), on a bowl and on
 * a bowl plus cylinder, which carries the saddle mode that converges slowest
 * at valences 7 and 8. Bounds from the fan table of
 * docs/irregular_patch_results.md section 9.6.
 */
TEST(IrregularPatchQuadratureTest, TheDefaultMatchesADeepReferenceWithFewerLevels)
{
    for (const double cylinder : {0.0, 0.12})
    {
        for (int i = 0; i < kExtraordinaryValenceCount; i++)
        {
            const int valence = kValences[i];
            SCOPED_TRACE("N=" + std::to_string(valence) + " cylinder " + std::to_string(cylinder));
            const double bowl = cylinder > 0.0 ? 0.0 : 0.12;

            Fan reference(with_quadrature({{16, 16, 16, 16}}, true), valence, bowl, cylinder);
            Fan corrected(fan_param(), valence, bowl, cylinder);
            Fan plain(with_quadrature(kDefaultIrregularPatchDepths, false), valence, bowl,
                      cylinder);
            Fan legacy(with_quadrature(kLegacyIrregularPatchDepths, false), valence, bowl,
                       cylinder);

            const double ref = reference.bending();
            ASSERT_GT(ref, 0.0);
            const double correctedError = std::abs(corrected.bending() - ref) / ref;
            const double plainError = std::abs(plain.bending() - ref) / ref;
            const double legacyError = std::abs(legacy.bending() - ref) / ref;

            EXPECT_LT(correctedError, 1.5e-4);
            EXPECT_GT(plainError, 10.0 * correctedError)
                << "the depths alone, without the correction, should be far worse";
            EXPECT_LT(correctedError, std::max(2.0 * legacyError, 1e-4));

            // Valence 4's area converges fast enough that legacy depth 12 is
            // tighter still; on a fluid sheet as a whole the default is 1e-7
            // against the legacy 8e-5.
            const double refArea = reference.area();
            EXPECT_LT(std::abs(corrected.area() - refArea) / refArea, 1e-5);
        }
    }
}

/**
 * @brief With the correction on, every force is still the exact gradient of
 * its energy.
 *
 * The area and volume weights differ from the bending weights on the last
 * level, so the kernel's area and volume forces are rescaled after the call;
 * this is what catches a rescale that drifts from the area pass. Central
 * differences over the control points of one fan face, per term.
 */
TEST(IrregularPatchQuadratureTest, TailCorrectedForcesAreTheGradientOfTheEnergy)
{
    for (const int valence : kValences)
    {
        SCOPED_TRACE("N=" + std::to_string(valence));
        Fan fan(fan_param(), valence, 0.12, 0.05);
        Mesh &mesh = fan.mesh;
        // Pull both constraints off their minima so they carry force.
        mesh.param.area0 = 0.97 * mesh.param.area;
        mesh.param.uVol = 50.0;
        mesh.param.vol0 = 0.9 * mesh.param.vol;
        ASSERT_NE(mesh.param.vol0, 0.0);
        mesh.Compute_Energy_And_Force();

        const std::vector<int> ring = mesh.faces[fan.fanFaces[0]].oneRingVertices;
        std::vector<double> analytic[3];
        double scale[3] = {0.0, 0.0, 0.0};
        for (int node : ring)
        {
            for (int k = 0; k < 3; k++)
            {
                const double values[3] = {mesh.vertices[node].force.forceCurvature(k, 0),
                                          mesh.vertices[node].force.forceArea(k, 0),
                                          mesh.vertices[node].force.forceVolume(k, 0)};
                for (int term = 0; term < 3; term++)
                {
                    analytic[term].push_back(values[term]);
                    scale[term] = std::max(scale[term], std::abs(values[term]));
                }
            }
        }
        for (int term = 0; term < 3; term++)
        {
            ASSERT_GT(scale[term], 0.0) << "term " << term << " carries no force";
        }

        const double h = 1e-6;
        double worst[3] = {0.0, 0.0, 0.0};
        int component = 0;
        for (int node : ring)
        {
            for (int k = 0; k < 3; k++, component++)
            {
                const double original = mesh.vertices[node].coord(k, 0);
                double plus[3], minus[3];
                mesh.vertices[node].coord.set(k, 0, original + h);
                mesh.Compute_Energy_And_Force();
                plus[0] = mesh.param.energy.energyCurvature;
                plus[1] = mesh.param.energy.energyArea;
                plus[2] = mesh.param.energy.energyVolume;
                mesh.vertices[node].coord.set(k, 0, original - h);
                mesh.Compute_Energy_And_Force();
                minus[0] = mesh.param.energy.energyCurvature;
                minus[1] = mesh.param.energy.energyArea;
                minus[2] = mesh.param.energy.energyVolume;
                mesh.vertices[node].coord.set(k, 0, original);
                for (int term = 0; term < 3; term++)
                {
                    const double gradient = (plus[term] - minus[term]) / (2.0 * h);
                    worst[term] = std::max(
                        worst[term], std::abs(gradient + analytic[term][component]) / scale[term]);
                }
            }
        }
        EXPECT_LT(worst[0], 1e-6) << "bending";
        EXPECT_LT(worst[1], 1e-6) << "area";
        EXPECT_LT(worst[2], 1e-6) << "volume";
    }
}
