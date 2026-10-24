#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

#include "Parameters.hpp"
#include "energy_force/Regularization_kernel.hpp"
#include "mesh/Mesh.hpp"

/**
 * @file test_regularization_blend.cpp
 * @brief The regularization energy is continuous and its force is its gradient.
 *
 * The R-pi cases used to switch hard at gamaShape and gamaArea. The targets
 * of the springs differ on the two sides of a threshold, so a face crossing
 * one changed the energy by a constant, and run_flat stopped wherever a face
 * sat on a threshold with the force pushing it across: no step, however
 * short, lowered the energy. The shape-deformed force also left out that its
 * target follows the face's own area. These tests pin the blended energy:
 *
 * - dE/dl from regularization_energy_and_gradient() matches central
 *   differences of the energy in every case and in both blend bands;
 * - the corner forces are minus the gradient in the coordinates;
 * - across a threshold the blended energy is continuous, where the hard
 *   switch (blendWidth 0) jumps;
 * - with blendWidth 0 the energy is the old case's, number for number;
 * - on a mesh, forceRegularization is minus the gradient of
 *   energyRegularization.
 */

namespace
{

using slimed::RegularizationParams;

RegularizationParams default_params(double blendWidth)
{
    RegularizationParams params;
    params.kCurv = 83.4;
    params.gamaShape = 0.2;
    params.gamaArea = 0.2;
    params.usingRpi = true;
    params.blendWidth = blendWidth;
    return params;
}

double energy_of(const double length[3], const double refLength[3], const RegularizationParams &params)
{
    double dEdl[3];
    unsigned char deformCase = 0;
    return slimed::regularization_energy_and_gradient(length, refLength, params, dEdl, deformCase);
}

/// Side-length scatter, as the regularizer measures it.
double gama_of(const double length[3])
{
    const double mean = (length[0] + length[1] + length[2]) / 3.0;
    double sum = 0.0;
    for (int j = 0; j < 3; j++)
    {
        sum += (length[j] - mean) * (length[j] - mean);
    }
    return sum / (mean * mean);
}

double heron(const double l[3])
{
    const double s = (l[0] + l[1] + l[2]) / 2.0;
    return std::sqrt(s * (s - l[0]) * (s - l[1]) * (s - l[2]));
}

/// The energy as a function of three corner positions.
double energy_of_corners(const double x[9], const double refLength[3], const RegularizationParams &params)
{
    double length[3];
    const int pairs[3][2] = {{0, 1}, {1, 2}, {2, 0}};
    for (int j = 0; j < 3; j++)
    {
        double sum = 0.0;
        for (int axis = 0; axis < 3; axis++)
        {
            const double d = x[pairs[j][0] * 3 + axis] - x[pairs[j][1] * 3 + axis];
            sum += d * d;
        }
        length[j] = std::sqrt(sum);
    }
    return energy_of(length, refLength, params);
}

} // namespace

/**
 * @brief dE/dl matches central differences, case by case and in the bands.
 *
 * Random triangles around a 5 nm equilateral one, with random reference
 * triangles, so that the sample covers the reference-length case, the
 * area-deformed and shape-deformed cases and both blend bands. Every
 * category is asserted to occur.
 */
TEST(RegularizationBlendTest, TheLengthGradientMatchesFiniteDifferences)
{
    std::mt19937 rng(20261019u);
    std::uniform_real_distribution<double> scatter(-2.4, 2.4);
    std::uniform_real_distribution<double> refScale(0.75, 1.3);
    const RegularizationParams params = default_params(0.1);

    int inShapeBand = 0, inAreaBand = 0, deepShape = 0, deepArea = 0, reference = 0;
    for (int sample = 0; sample < 20000; sample++)
    {
        double length[3], refLength[3];
        const double scale = refScale(rng);
        for (int j = 0; j < 3; j++)
        {
            length[j] = 5.0 + scatter(rng);
            refLength[j] = scale * (5.0 + 0.3 * scatter(rng));
        }
        // Keep clear of degenerate triangles, where dA/dl blows up.
        const double area = heron(length);
        const double refArea = heron(refLength);
        if (!(area > 2.0) || !(refArea > 2.0))
        {
            continue;
        }
        const double gama = gama_of(length);
        const double excess = std::abs(area - refArea) / refArea;
        // Stay off the band edges, where the smoothstep's second derivative
        // jumps and central differences lose an order.
        const double tShape = (gama - params.gamaShape) / params.blendWidth;
        const double tArea = (excess - params.gamaArea) / params.blendWidth;
        if (std::abs(tShape) < 1e-3 || std::abs(tShape - 1.0) < 1e-3 ||
            std::abs(tArea) < 1e-3 || std::abs(tArea - 1.0) < 1e-3)
        {
            continue;
        }
        if (tShape > 0.0 && tShape < 1.0) inShapeBand++;
        else if (tShape >= 1.0) deepShape++;
        else if (tArea > 0.0 && tArea < 1.0) inAreaBand++;
        else if (tArea >= 1.0) deepArea++;
        else reference++;

        double dEdl[3];
        unsigned char deformCase = 0;
        slimed::regularization_energy_and_gradient(length, refLength, params, dEdl, deformCase);
        const double h = 1e-6;
        for (int j = 0; j < 3; j++)
        {
            double plus[3] = {length[0], length[1], length[2]};
            double minus[3] = {length[0], length[1], length[2]};
            plus[j] += h;
            minus[j] -= h;
            const double numeric = (energy_of(plus, refLength, params) - energy_of(minus, refLength, params)) / (2.0 * h);
            ASSERT_NEAR(dEdl[j], numeric, 1e-6 * std::max(1.0, std::abs(numeric)))
                << "sample " << sample << " side " << j << " gama " << gama << " excess " << excess;
        }
    }
    EXPECT_GT(inShapeBand, 20);
    EXPECT_GT(inAreaBand, 20);
    EXPECT_GT(deepShape, 20);
    EXPECT_GT(deepArea, 20);
    EXPECT_GT(reference, 20);
}

/// The corner forces are minus the gradient of the energy in the coordinates.
TEST(RegularizationBlendTest, TheCornerForcesAreMinusTheCoordinateGradient)
{
    std::mt19937 rng(7u);
    std::uniform_real_distribution<double> jitter(-1.2, 1.2);
    const double base[9] = {0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 2.5, 4.33, 0.0};
    const double refLength[3] = {5.0, 5.0, 5.0};
    for (const double width : {0.0, 0.05, 0.1})
    {
        const RegularizationParams params = default_params(width);
        for (int sample = 0; sample < 200; sample++)
        {
            double x[9];
            for (int i = 0; i < 9; i++)
            {
                x[i] = base[i] + jitter(rng);
            }
            double length[3], u10[3], u21[3], u02[3];
            const int pairs[3][2] = {{0, 1}, {1, 2}, {2, 0}};
            double *units[3] = {u10, u21, u02};
            for (int j = 0; j < 3; j++)
            {
                double sum = 0.0;
                for (int axis = 0; axis < 3; axis++)
                {
                    units[j][axis] = x[pairs[j][0] * 3 + axis] - x[pairs[j][1] * 3 + axis];
                    sum += units[j][axis] * units[j][axis];
                }
                length[j] = std::sqrt(sum);
                for (int axis = 0; axis < 3; axis++)
                {
                    units[j][axis] /= length[j];
                }
            }
            if (!(heron(length) > 2.0))
            {
                continue;
            }
            // The hard switch has a gradient only away from its thresholds.
            const double tShape = gama_of(length) - params.gamaShape;
            const double tArea = std::abs(heron(length) - heron(refLength)) / heron(refLength) - params.gamaArea;
            if (width == 0.0 && (std::abs(tShape) < 1e-4 || std::abs(tArea) < 1e-4))
            {
                continue;
            }

            double dEdl[3];
            unsigned char deformCase = 0;
            slimed::regularization_energy_and_gradient(length, refLength, params, dEdl, deformCase);
            double force[3][3];
            slimed::regularization_corner_forces(dEdl, u10, u21, u02, force[0], force[1], force[2]);

            const double h = 1e-6;
            for (int i = 0; i < 9; i++)
            {
                double plus[9], minus[9];
                std::copy(x, x + 9, plus);
                std::copy(x, x + 9, minus);
                plus[i] += h;
                minus[i] -= h;
                const double numeric = -(energy_of_corners(plus, refLength, params) -
                                         energy_of_corners(minus, refLength, params)) /
                                       (2.0 * h);
                ASSERT_NEAR(force[i / 3][i % 3], numeric, 1e-6 * std::max(1.0, std::abs(numeric)))
                    << "width " << width << " sample " << sample << " coordinate " << i;
            }
        }
    }
}

/**
 * @brief Across gamaShape the blended energy is continuous; the hard switch jumps.
 *
 * A triangle is stretched along its base at constant area, from equilateral
 * past the shape threshold. Its reference is the stretched triangle a little
 * short of the threshold, as the minimizer leaves it after refreshing: the
 * reference-length energy is small there, the shape-deformed one is not.
 */
TEST(RegularizationBlendTest, TheBlendedEnergyIsContinuousAcrossTheShapeThreshold)
{
    // Corners (0,0), (5s,0), (2.5s, 4.33/s): area 10.825 for every stretch s.
    auto sides = [](double s, double length[3]) {
        const double base = 5.0 * s;
        const double height = 4.330127 / s;
        const double leg = std::sqrt(base * base / 4.0 + height * height);
        length[0] = base;
        length[1] = leg;
        length[2] = leg;
    };
    double refLength[3];
    sides(1.5, refLength);
    for (const double width : {0.0, 0.05})
    {
        const RegularizationParams params = default_params(width);
        double largestJump = 0.0;
        bool crossed = false;
        double previous = 0.0;
        const int nSteps = 100000;
        for (int step = 0; step <= nSteps; step++)
        {
            const double s = 1.0 + 1.0 * step / nSteps;
            double length[3];
            sides(s, length);
            const double energy = energy_of(length, refLength, params);
            if (step > 0)
            {
                largestJump = std::max(largestJump, std::abs(energy - previous));
            }
            crossed = crossed || gama_of(length) > params.gamaShape + width;
            previous = energy;
        }
        ASSERT_TRUE(crossed);
        if (width == 0.0)
        {
            EXPECT_GT(largestJump, 100.0) << "the hard switch should jump";
        }
        else
        {
            EXPECT_LT(largestJump, 0.1) << "the blend should not jump";
        }
    }
}

/// With blendWidth 0 every face is in exactly one of the old cases.
TEST(RegularizationBlendTest, TheHardSwitchReproducesTheOldCases)
{
    const double k = 83.4;
    const RegularizationParams params = default_params(0.0);
    const double refLength[3] = {5.0, 5.2, 4.9};
    const double refArea = heron(refLength);

    // Reference-length case.
    {
        const double length[3] = {5.1, 5.15, 4.95};
        double dEdl[3];
        unsigned char deformCase = 9;
        const double e = slimed::regularization_energy_and_gradient(length, refLength, params, dEdl, deformCase);
        EXPECT_EQ(deformCase, 0);
        const double expected = k / 2.0 * (std::pow(5.1 - 5.0, 2.0) + std::pow(5.15 - 5.2, 2.0) + std::pow(4.95 - 4.9, 2.0));
        EXPECT_DOUBLE_EQ(e, expected);
    }
    // Area-deformed: uniformly scaled up by 1.2 in length, 1.44 in area.
    {
        const double length[3] = {6.0, 6.24, 5.88};
        double dEdl[3];
        unsigned char deformCase = 9;
        const double e = slimed::regularization_energy_and_gradient(length, refLength, params, dEdl, deformCase);
        EXPECT_EQ(deformCase, 1);
        const double target = std::sqrt(4.0 * refArea / std::sqrt(3.0));
        const double expected = k / 2.0 * (std::pow(6.0 - target, 2.0) + std::pow(6.24 - target, 2.0) + std::pow(5.88 - target, 2.0));
        EXPECT_DOUBLE_EQ(e, expected);
    }
    // Shape-deformed.
    {
        const double length[3] = {7.5, 4.0, 4.5};
        ASSERT_GT(gama_of(length), 0.2);
        double dEdl[3];
        unsigned char deformCase = 9;
        const double e = slimed::regularization_energy_and_gradient(length, refLength, params, dEdl, deformCase);
        EXPECT_EQ(deformCase, 2);
        const double target = std::sqrt(4.0 * heron(length) / std::sqrt(3.0));
        const double expected = k / 2.0 * (std::pow(7.5 - target, 2.0) + std::pow(4.0 - target, 2.0) + std::pow(4.5 - target, 2.0));
        EXPECT_DOUBLE_EQ(e, expected);
    }
    // usingRpi = false: the reference-length case whatever the shape.
    {
        RegularizationParams off = params;
        off.usingRpi = false;
        off.blendWidth = 0.05;
        const double length[3] = {7.5, 4.0, 4.5};
        double dEdl[3];
        unsigned char deformCase = 9;
        const double e = slimed::regularization_energy_and_gradient(length, refLength, off, dEdl, deformCase);
        EXPECT_EQ(deformCase, 0);
        const double expected = k / 2.0 * (std::pow(7.5 - 5.0, 2.0) + std::pow(4.0 - 5.2, 2.0) + std::pow(4.5 - 4.9, 2.0));
        EXPECT_DOUBLE_EQ(e, expected);
    }
}

/**
 * @brief On a mesh, forceRegularization is minus the gradient of the energy.
 *
 * A clamped sheet with its interior vertices jittered in and out of plane, so
 * that faces sit in every case and in the bands. Each interior vertex is
 * moved along each axis and the regularization energy differenced. The
 * clamped boundary keeps the ghost handling out of this: it is the term that
 * is checked, not the boundary condition.
 */
TEST(RegularizationBlendTest, TheMeshForceIsMinusTheGradientOfTheEnergy)
{
    Param param;
    param.VERBOSE_MODE = false;
    param.boundaryCondition = BoundaryType::Fixed;
    param.sideX = 40.0;
    param.sideY = 40.0;
    param.lFace = 5.0;
    param.kCurv = 83.4;
    param.uSurf = 250.0;
    param.uVol = 0.0;
    param.isGlobalConstraint = true;
    param.isEnergyHarmonicBondIncluded = false;
    param.usingRpi = true;
    param.regularizationBlendWidth = 0.1;

    Mesh mesh(param);
    mesh.setup_flat();
    mesh.calculate_element_area_volume();
    mesh.sum_membrane_area_and_volume(param.area0, param.vol0);
    mesh.update_previous_coord_for_vertex();
    mesh.update_reference_coord_from_previous_coord();

    std::mt19937 rng(11u);
    std::uniform_real_distribution<double> jitter(-1.3, 1.3);
    std::vector<int> interior;
    for (int i = 0; i < static_cast<int>(mesh.vertices.size()); i++)
    {
        if (!mesh.vertices[i].isBoundary && !mesh.vertices[i].isGhost)
        {
            interior.push_back(i);
            for (int axis = 0; axis < 3; axis++)
            {
                mesh.vertices[i].coord.set(axis, 0, mesh.vertices[i].coord(axis, 0) + jitter(rng));
            }
        }
    }
    ASSERT_GT(interior.size(), 10u);

    mesh.Compute_Energy_And_Force();
    EXPECT_GT(param.deformationCount.shapeDeformCount, 0);
    EXPECT_LT(param.deformationCount.shapeDeformCount, static_cast<int>(mesh.faces.size()));
    std::vector<double> force(mesh.vertices.size() * 3);
    for (int i : interior)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            force[i * 3 + axis] = mesh.vertices[i].force.forceRegularization(axis, 0);
        }
    }

    const double h = 1e-6;
    double largestForce = 0.0;
    for (int i : interior)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            const double x = mesh.vertices[i].coord(axis, 0);
            mesh.vertices[i].coord.set(axis, 0, x + h);
            mesh.Compute_Energy_And_Force();
            const double ePlus = param.energy.energyRegularization;
            mesh.vertices[i].coord.set(axis, 0, x - h);
            mesh.Compute_Energy_And_Force();
            const double eMinus = param.energy.energyRegularization;
            mesh.vertices[i].coord.set(axis, 0, x);
            const double numeric = -(ePlus - eMinus) / (2.0 * h);
            largestForce = std::max(largestForce, std::abs(numeric));
            // The energy is a sum over faces: a reduction under OpenMP, so
            // allow for its rounding as well as the difference's.
            EXPECT_NEAR(force[i * 3 + axis], numeric, 1e-5 * std::max(1.0, std::abs(numeric)))
                << "vertex " << i << " axis " << axis;
        }
    }
    EXPECT_GT(largestForce, 1.0);
}
