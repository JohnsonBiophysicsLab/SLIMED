/**
 * @file Patch_kernel.inl
 * @brief Definition of element_energy_force_patch_pod().
 *
 * Kept out of the header's own text, but included from its foot, so that a
 * CUDA translation unit can compile this same source for the device with
 * SLIMED_HD expanded to __host__ __device__. The functions are inline, so
 * every caller carries a definition and the device needs no relocatable device
 * code. There is no device-side variant of this code to keep in sync -- there
 * is one body.
 */
#pragma once

#include "energy_force/Patch_kernel.hpp"

namespace slimed
{

SLIMED_HD inline void element_energy_force_patch_pod(const double *sampleRows,
                                              const double *gaussCoeff,
                                              int nSamples,
                                              const double *ctrlPts,
                                              int nCtrl,
                                              const PatchParams &p,
                                              double &eBend,
                                              double &meanCurv,
                                              double normVector[3],
                                              double *fBend,
                                              double *fArea,
                                              double *fVolume)
{
    eBend = 0.0;
    meanCurv = 0.0;

    const int rowStride = kShapeRows * nCtrl;

    for (int s = 0; s < nSamples; s++)
    {
        const double halfGaussQuadratureCoeff = 0.5 * gaussCoeff[s];
        const double *rows = sampleRows + s * rowStride;

        // The limit-surface point and its derivatives at this quadrature point.
        //   sfDot[0] : x,    the point itself
        //   sfDot[1] : a_1 = dx/du      sfDot[2] : a_2  = dx/dv
        //   sfDot[3] : a_11             sfDot[4] : a_22
        //   sfDot[5] : a_12             sfDot[6] : a_21
        double sfDot[kShapeRows][3];
        shape_times_control_points(rows, ctrlPts, nCtrl, sfDot);

        const double *x = sfDot[0];
        const double *a_1 = sfDot[1];
        const double *a_2 = sfDot[2];
        const double *a_11 = sfDot[3];
        const double *a_22 = sfDot[4];
        const double *a_12 = sfDot[5];
        const double *a_21 = sfDot[6];

        // Surface metric: xa is normal to the tangent plane and its length is
        // the area element.
        double xa[3];
        v3_cross(a_1, a_2, xa);
        const double sqa = std::sqrt(v3_dot(xa, xa));
        const double sqaSqrInv = 1.0 / sqa / sqa;

        double xa_1[3];
        double xa_2[3];
        v3_cross_plus_cross(a_11, a_2, a_1, a_21, xa_1);
        v3_cross_plus_cross(a_12, a_2, a_1, a_22, xa_2);

        const double sqa_1 = v3_dot(xa, xa_1) / sqa;
        const double sqa_2 = v3_dot(xa, xa_2) / sqa;

        // Unit normal a_3 and its two derivatives, by the quotient rule.
        double a_3[3];
        double a_31[3];
        double a_32[3];
        v3_scale_recip(xa, sqa, a_3);
        for (int i = 0; i < 3; i++)
        {
            a_31[i] = (xa_1[i] * sqa - xa[i] * sqa_1) * sqaSqrInv;
            a_32[i] = (xa_2[i] * sqa - xa[i] * sqa_2) * sqaSqrInv;
        }

        // Contravariant basis a1, a2 and its derivatives.
        double crossA2A3[3];
        double crossA3A1[3];
        v3_cross(a_2, a_3, crossA2A3);
        v3_cross(a_3, a_1, crossA3A1);

        double a1[3];
        double a2[3];
        v3_scale_recip(crossA2A3, sqa, a1);
        v3_scale_recip(crossA3A1, sqa, a2);

        double a11[3];
        double a12[3];
        double a21[3];
        double a22[3];
        v3_quotient_rule(a_21, a_3, a_2, a_31, sqa, crossA2A3, sqa_1, sqaSqrInv, a11);
        v3_quotient_rule(a_22, a_3, a_2, a_32, sqa, crossA2A3, sqa_2, sqaSqrInv, a12);
        v3_quotient_rule(a_31, a_1, a_3, a_11, sqa, crossA3A1, sqa_1, sqaSqrInv, a21);
        v3_quotient_rule(a_32, a_1, a_3, a_12, sqa, crossA3A1, sqa_2, 1.0 / sqa / sqa, a22);

        // Mean curvature, and the bending stress and couple resultants.
        const double H_curv = 0.5 * (v3_dot(a1, a_31) + v3_dot(a2, a_32));

        const double curvatureExcess = 2.0 * H_curv - p.spontCurv;
        const double bendScale = -p.kCurv * curvatureExcess;
        const double bendSquare = p.kCurv * 0.5 * (curvatureExcess * curvatureExcess);

        double n1_be[3];
        double n2_be[3];
        const double a1DotA1 = v3_dot(a1, a1);
        const double a1DotA2 = v3_dot(a1, a2);
        const double a2DotA1 = v3_dot(a2, a1);
        const double a2DotA2 = v3_dot(a2, a2);
        for (int i = 0; i < 3; i++)
        {
            n1_be[i] = (a_31[i] * a1DotA1 + a_32[i] * a1DotA2) * bendScale + a1[i] * bendSquare;
            n2_be[i] = (a_31[i] * a2DotA1 + a_32[i] * a2DotA2) * bendScale + a2[i] * bendSquare;
        }

        double m1_be[3];
        double m2_be[3];
        v3_scale(a1, -bendScale, m1_be);
        v3_scale(a2, -bendScale, m2_be);

        // Area constraint. The original computes a local-mode alternative here
        // and then unconditionally overwrites it with the global-mode value,
        // so local mode has never actually reached the force. Reproduced as
        // is: this kernel is a performance change, not a physics change.
        const double areaScale = p.uSurfPerArea * (p.area - p.area0);
        double n1_cons[3];
        double n2_cons[3];
        v3_scale(a1, areaScale, n1_cons);
        v3_scale(a2, areaScale, n2_cons);

        // Volume constraint.
        const double volScale = p.uVol * (p.vol - p.vol0) / 3.0;
        const double xDotA_3 = v3_dot(x, a_3);
        const double xDotA1 = v3_dot(x, a1);
        const double xDotA2 = v3_dot(x, a2);
        double n1_conv[3];
        double n2_conv[3];
        for (int i = 0; i < 3; i++)
        {
            n1_conv[i] = (a1[i] * xDotA_3 - a_3[i] * xDotA1) * volScale;
            n2_conv[i] = (a2[i] * xDotA_3 - a_3[i] * xDotA2) * volScale;
        }

        const double eBendSample = 0.5 * p.kCurv * sqa * (curvatureExcess * curvatureExcess);

        // Per control point: differentiate the basis with respect to that
        // point's coordinate and contract with the resultants above. The force
        // is the negative derivative of the energy, hence the -sqa.
        for (int j = 0; j < nCtrl; j++)
        {
            const double sf0 = rows[0 * nCtrl + j];
            const double sf1 = rows[1 * nCtrl + j];
            const double sf2 = rows[2 * nCtrl + j];
            const double sf3 = rows[3 * nCtrl + j];
            const double sf4 = rows[4 * nCtrl + j];
            const double sf5 = rows[5 * nCtrl + j];
            const double sf6 = rows[6 * nCtrl + j];

            // da1 and da2 are the outer-product sums the original built with
            // six kron() calls apiece -- each of which allocated a 3x3 matrix,
            // copied it out by value, and threw the copy away.
            double da1[3][3];
            double da2[3][3];
            for (int r = 0; r < 3; r++)
            {
                for (int c = 0; c < 3; c++)
                {
                    da1[r][c] = a1[r] * a_3[c] * -sf3
                              + a11[r] * a_3[c] * -sf1
                              + a1[r] * a_31[c] * -sf1
                              + a2[r] * a_3[c] * -sf6
                              + a21[r] * a_3[c] * -sf2
                              + a2[r] * a_31[c] * -sf2;
                    da2[r][c] = a1[r] * a_3[c] * -sf5
                              + a12[r] * a_3[c] * -sf1
                              + a1[r] * a_32[c] * -sf1
                              + a2[r] * a_3[c] * -sf4
                              + a22[r] * a_3[c] * -sf2
                              + a2[r] * a_32[c] * -sf2;
                }
            }

            double *outBend = fBend + j * 3;
            double *outArea = fArea + j * 3;
            double *outVol = fVolume + j * 3;
            for (int i = 0; i < 3; i++)
            {
                // m_be contracted with da over the first index, the
                // colvec_matrix_multiplication() helper.
                const double couple1 = m1_be[0] * da1[0][i] + m1_be[1] * da1[1][i] + m1_be[2] * da1[2][i];
                const double couple2 = m2_be[0] * da2[0][i] + m2_be[1] * da2[1][i] + m2_be[2] * da2[2][i];
                const double couple = couple1 + couple2;

                const double fBendSample = (couple + n1_be[i] * sf1 + n2_be[i] * sf2) * -sqa;
                const double fAreaSample = (n1_cons[i] * sf1 + n2_cons[i] * sf2) * -sqa;
                const double fVolSample = (n1_conv[i] * sf1 + n2_conv[i] * sf2 + a_3[i] * (volScale * sf0)) * -sqa;

                outBend[i] += halfGaussQuadratureCoeff * fBendSample;
                outArea[i] += halfGaussQuadratureCoeff * fAreaSample;
                outVol[i] += halfGaussQuadratureCoeff * fVolSample;
            }
        }

        meanCurv += halfGaussQuadratureCoeff * H_curv;
        eBend += halfGaussQuadratureCoeff * eBendSample;
        for (int i = 0; i < 3; i++)
        {
            normVector[i] += halfGaussQuadratureCoeff * a_3[i];
        }
    }

    v3_normalize(normVector);
}

SLIMED_HD inline void element_area_volume_pod(const double *sampleRows,
                                       const double *gaussCoeff,
                                       int nSamples,
                                       const double *ctrlPts,
                                       int nCtrl,
                                       double &area,
                                       double &volume)
{
    // Folded straight into the caller's accumulators, one sample at a time.
    // The code this replaces carried an OpenMP reduction over the samples, but
    // that pragma is compiled out of the serial build and runs single-threaded
    // in the parallel one, so a per-call partial sum is not what either build
    // actually computed. Summing locally first would associate the additions
    // differently and drift from the energy the area is constrained against.
    const int rowStride = kShapeRows * nCtrl;
    for (int s = 0; s < nSamples; s++)
    {
        const double *rows = sampleRows + s * rowStride;

        // Only x, a_1 and a_2 -- rows 0, 1 and 2.
        double x[3] = {0.0, 0.0, 0.0};
        double a_1[3] = {0.0, 0.0, 0.0};
        double a_2[3] = {0.0, 0.0, 0.0};
        for (int k = 0; k < nCtrl; k++)
        {
            const double wx = rows[0 * nCtrl + k];
            const double w1 = rows[1 * nCtrl + k];
            const double w2 = rows[2 * nCtrl + k];
            const double *point = ctrlPts + k * 3;
            for (int axis = 0; axis < 3; axis++)
            {
                x[axis] += wx * point[axis];
                a_1[axis] += w1 * point[axis];
                a_2[axis] += w2 * point[axis];
            }
        }

        double a_3[3];
        v3_cross(a_1, a_2, a_3);
        const double sqa = std::sqrt(v3_dot(a_3, a_3));

        const double coeff = gaussCoeff[s];
        area += 0.5 * coeff * sqa;
        volume += kSignedVolumeQuadratureFactor * coeff * v3_dot(x, a_3);
    }
}

} // namespace slimed
