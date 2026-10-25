/**
 * @file Regularization_kernel.hpp
 * @brief The per-face regularization energy and its gradient, shared by the
 * CPU loop (Mesh::energy_force_regularization()) and the device backends
 * (regularization_for_face() in Force_kernels.inl).
 *
 * The energy is a spring on each of the three sides of a face, pulling it
 * towards a target length. Three cases pick the target:
 *
 *   0. reference lengths: each side towards its own length in coordRef;
 *   1. area-deformed: every side towards the side of the equilateral triangle
 *      of the reference area, once the area is off the reference area by the
 *      fraction gamaArea;
 *   2. shape-deformed: every side towards the side of the equilateral triangle
 *      of the face's current area, once the side lengths scatter by more than
 *      gamaShape. Wins over case 1.
 *
 * Each case is smooth, but the switch between them was not: the three targets
 * differ at the thresholds, so a face crossing one changed the energy by a
 * constant (a few hundred pN nm on the shipped bud), and a line search that
 * straddled the step found no step, however short, that lowered the energy.
 * With blendWidth > 0 the cases are blended instead. The shape weight rises
 * from 0 at gama = gamaShape to 1 at gamaShape + blendWidth (a cubic
 * smoothstep, so the energy is C1), and the area weight does the same in the
 * fractional area excess. Below the thresholds nothing changes. Past them plus
 * blendWidth a face is in the old case. In between the energy climbs
 * smoothly, which is the soft version of the wall the jump used to be.
 *
 * Everything is a function of the three side lengths, so the force on the
 * corners follows from dE/dl through the unit side vectors. The gradient is
 * exact, including the terms the old shape-deformed force left out: its
 * target follows the current area, so moving a corner also moves the target.
 */
#pragma once

#include <cmath>

#include "energy_force/Patch_kernel.hpp" // SLIMED_HD

namespace slimed
{

struct RegularizationParams
{
    double kCurv = 0.0;      ///< Spring constant of every side
    double gamaShape = 0.2;  ///< Param::gamaShape
    double gamaArea = 0.2;   ///< Param::gamaArea
    bool usingRpi = true;    ///< Param::usingRpi: false leaves every face in case 0
    double blendWidth = 0.0; ///< Param::regularizationBlendWidth; 0 switches hard
};

SLIMED_HD inline double regularization_square(double x)
{
    return x * x;
}

/**
 * @brief The cubic smoothstep 3t^2 - 2t^3 on [0, 1], clamped outside, and its
 * derivative.
 */
SLIMED_HD inline double regularization_smoothstep(double t, double &slope)
{
    if (t <= 0.0)
    {
        slope = 0.0;
        return 0.0;
    }
    if (t >= 1.0)
    {
        slope = 0.0;
        return 1.0;
    }
    slope = 6.0 * t * (1.0 - t);
    return t * t * (3.0 - 2.0 * t);
}

/**
 * @brief Regularization energy of one face and its derivative with respect to
 * the three side lengths.
 *
 * @param length Current side lengths, in the order 10, 21, 02.
 * @param refLength Reference side lengths, same order.
 * @param params Spring constant, thresholds and blend width.
 * @param[out] dEdl dE / d length[j].
 * @param[out] deformCase 2 when the face carries any shape-deformed weight,
 * else 1 when it carries any area-deformed weight, else 0. With blendWidth 0
 * this is the case the face is in, as it has always been reported.
 * @return The energy.
 */
SLIMED_HD inline double regularization_energy_and_gradient(const double length[3],
                                                           const double refLength[3],
                                                           const RegularizationParams &params,
                                                           double dEdl[3],
                                                           unsigned char &deformCase)
{
    const double k = params.kCurv;
    const double l10 = length[0];
    const double l21 = length[1];
    const double l02 = length[2];

    // Heron's formula on the current and the reference triangle.
    double semiperi = (l10 + l21 + l02) / 2.0;
    const double area = std::sqrt(semiperi * (semiperi - l10) * (semiperi - l21) * (semiperi - l02));
    semiperi = (refLength[0] + refLength[1] + refLength[2]) / 2.0;
    const double refArea = std::sqrt(semiperi * (semiperi - refLength[0]) *
                                     (semiperi - refLength[1]) * (semiperi - refLength[2]));

    // Shape deformation: the scatter of the side lengths about their mean,
    // zero for an equilateral triangle.
    const double meanSideLength = (l10 + l21 + l02) / 3.0;
    const double gama = (regularization_square(l10 - meanSideLength) +
                         regularization_square(l21 - meanSideLength) +
                         regularization_square(l02 - meanSideLength)) /
                        regularization_square(meanSideLength);
    const double areaExcess = (area - refArea) / refArea;
    const double absAreaExcess = areaExcess < 0.0 ? -areaExcess : areaExcess;

    // The weights of the shape- and area-deformed cases.
    double shapeWeight = 0.0, shapeSlope = 0.0; // slope: d weight / d gama
    double areaWeight = 0.0, areaSlope = 0.0;   // slope: d weight / d |excess|
    if (params.usingRpi)
    {
        if (params.blendWidth > 0.0)
        {
            shapeWeight = regularization_smoothstep((gama - params.gamaShape) / params.blendWidth,
                                                    shapeSlope);
            shapeSlope /= params.blendWidth;
            areaWeight = regularization_smoothstep(
                (absAreaExcess - params.gamaArea) / params.blendWidth, areaSlope);
            areaSlope /= params.blendWidth;
        }
        else
        {
            shapeWeight = gama > params.gamaShape ? 1.0 : 0.0;
            areaWeight = absAreaExcess >= params.gamaArea ? 1.0 : 0.0;
        }
    }
    deformCase = shapeWeight > 0.0 ? 2 : (areaWeight > 0.0 ? 1 : 0);

    // dA / dl_j = l_j (l_k^2 + l_m^2 - l_j^2) / (8 A), from 16 A^2 =
    // 2(a^2 b^2 + b^2 c^2 + c^2 a^2) - (a^4 + b^4 + c^4). Only needed where a
    // weight moves or a target follows the area.
    double dAdl[3] = {0.0, 0.0, 0.0};
    const bool needAreaGradient =
        (shapeWeight > 0.0 || areaSlope != 0.0) && area > 0.0;
    if (needAreaGradient)
    {
        const double s10 = l10 * l10, s21 = l21 * l21, s02 = l02 * l02;
        dAdl[0] = l10 * (s21 + s02 - s10) / (8.0 * area);
        dAdl[1] = l21 * (s02 + s10 - s21) / (8.0 * area);
        dAdl[2] = l02 * (s10 + s21 - s02) / (8.0 * area);
    }

    // Case 0, towards the reference lengths.
    double energyRef = 0.0;
    double gradRef[3] = {0.0, 0.0, 0.0};
    if (shapeWeight < 1.0 && areaWeight < 1.0)
    {
        energyRef = k / 2.0 * (regularization_square(l10 - refLength[0]) + regularization_square(l21 - refLength[1]) +
                               regularization_square(l02 - refLength[2]));
        for (int j = 0; j < 3; j++)
        {
            gradRef[j] = k * (length[j] - refLength[j]);
        }
    }

    // Case 1, towards the equilateral triangle of the reference area.
    double energyArea = 0.0;
    double gradArea[3] = {0.0, 0.0, 0.0};
    if (shapeWeight < 1.0 && areaWeight > 0.0)
    {
        const double target = std::sqrt(4.0 * refArea / std::sqrt(3.0));
        energyArea = k / 2.0 * (regularization_square(l10 - target) + regularization_square(l21 - target) +
                                regularization_square(l02 - target));
        for (int j = 0; j < 3; j++)
        {
            gradArea[j] = k * (length[j] - target);
        }
    }

    // Case 2, towards the equilateral triangle of the current area. Its
    // target L = sqrt(4 A / sqrt 3) moves with the corners, dL/dA = L / (2A).
    double energyShape = 0.0;
    double gradShape[3] = {0.0, 0.0, 0.0};
    if (shapeWeight > 0.0)
    {
        const double target = std::sqrt(4.0 * area / std::sqrt(3.0));
        energyShape = k / 2.0 * (regularization_square(l10 - target) + regularization_square(l21 - target) +
                                 regularization_square(l02 - target));
        const double sumDeviation = (l10 - target) + (l21 - target) + (l02 - target);
        const double dTargetdA = area > 0.0 ? target / (2.0 * area) : 0.0;
        for (int j = 0; j < 3; j++)
        {
            gradShape[j] = k * (length[j] - target) - k * sumDeviation * dTargetdA * dAdl[j];
        }
    }

    // Blend: below the shape-deformed case, case 1 over case 0 by the area
    // weight; over both, case 2 by the shape weight.
    const double energyLow = areaWeight * energyArea + (1.0 - areaWeight) * energyRef;
    const double energy = shapeWeight * energyShape + (1.0 - shapeWeight) * energyLow;

    // d gama / dl_j = 2 (l_j - m) / m^2 - (2/3) gama / m, and
    // d |excess| / dl_j = sign(excess) dA/dl_j / A_ref.
    const double areaExcessSign = areaExcess < 0.0 ? -1.0 : 1.0;
    for (int j = 0; j < 3; j++)
    {
        double dShapeWeight = 0.0;
        if (shapeSlope != 0.0)
        {
            const double dGama = 2.0 * (length[j] - meanSideLength) / (meanSideLength * meanSideLength) -
                                 2.0 / 3.0 * gama / meanSideLength;
            dShapeWeight = shapeSlope * dGama;
        }
        double dAreaWeight = 0.0;
        if (areaSlope != 0.0)
        {
            dAreaWeight = areaSlope * areaExcessSign * dAdl[j] / refArea;
        }
        const double dLow = areaWeight * gradArea[j] + (1.0 - areaWeight) * gradRef[j] +
                            dAreaWeight * (energyArea - energyRef);
        dEdl[j] = shapeWeight * gradShape[j] + (1.0 - shapeWeight) * dLow +
                  dShapeWeight * (energyShape - energyLow);
    }
    return energy;
}

/**
 * @brief The corner forces of one face from dE/dl and its unit side vectors
 * u10 = (x0 - x1)/l10, u21 = (x1 - x2)/l21, u02 = (x2 - x0)/l02.
 */
SLIMED_HD inline void regularization_corner_forces(const double dEdl[3], const double u10[3],
                                                   const double u21[3], const double u02[3],
                                                   double force0[3], double force1[3],
                                                   double force2[3])
{
    for (int axis = 0; axis < 3; axis++)
    {
        force0[axis] = dEdl[2] * u02[axis] - dEdl[0] * u10[axis];
        force1[axis] = dEdl[0] * u10[axis] - dEdl[1] * u21[axis];
        force2[axis] = dEdl[1] * u21[axis] - dEdl[2] * u02[axis];
    }
}

} // namespace slimed
