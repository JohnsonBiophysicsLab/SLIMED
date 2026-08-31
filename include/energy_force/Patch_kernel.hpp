/**
 * @file Patch_kernel.hpp
 * @brief The limit-surface energy/force kernel, written over plain arrays.
 *
 * This is the hot loop of the whole program: profiling the default 8400-face
 * flat membrane put 92% of wall time inside Mesh::element_energy_force_patch,
 * and roughly 70% of that in malloc/free and cross-library gsl_matrix_get /
 * gsl_matrix_set calls rather than in arithmetic. The cause was the Matrix
 * class: every 3-vector was a separate gsl_matrix heap allocation, and the
 * kernel declared 43 of them per call, then built and discarded a dozen more
 * per control point through kron()'s return-by-value.
 *
 * Everything here is fixed-size and lives on the stack, so the same source
 * compiles for the host and, with SLIMED_HD expanding to __host__ __device__,
 * for a CUDA device. There are no allocations, no STL containers, no virtual
 * calls, and no library boundaries to cross.
 *
 * The arithmetic deliberately mirrors the Matrix-based original operation for
 * operation, including its association order, so results match to rounding.
 * tests/test_patch_kernel.cpp pins that against the preserved reference
 * implementation. Where the original does something odd -- reciprocal-then-
 * multiply in const_division(), a dead local-area branch -- this reproduces
 * the oddity rather than quietly fixing it; changing the numbers is a separate
 * question from making them fast.
 */
#pragma once

#include <cmath>

#if defined(__CUDACC__)
#define SLIMED_HD __host__ __device__
#else
#define SLIMED_HD
#endif

namespace slimed
{

/// Widest control-point net one patch can carry. A regular patch has 12; an
/// irregular one has valence + 6, and valence is capped at
/// kMaxIrregularValence = 8. Sizing the stack buffers by the maximum keeps the
/// kernel allocation-free for both patch kinds.
constexpr int kMaxControlPoints = 14;

/// Rows in a shape-function block: the value and its six derivatives
/// (du, dv, duu, dvv, duv, dvu), in that order.
constexpr int kShapeRows = 7;

/// (1/3) from the divergence theorem times (1/2) for the reference triangle's
/// own area. See element_area_volume_pod().
constexpr double kSignedVolumeQuadratureFactor = 1.0 / 6.0;

// --------------------------------------------------------------------------
// 3-vector primitives. Each mirrors one free function from Linear_algebra.cpp.
// --------------------------------------------------------------------------

SLIMED_HD inline void v3_zero(double a[3])
{
    a[0] = 0.0;
    a[1] = 0.0;
    a[2] = 0.0;
}

SLIMED_HD inline double v3_dot(const double a[3], const double b[3])
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

SLIMED_HD inline double v3_norm(const double a[3])
{
    return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

SLIMED_HD inline void v3_cross(const double a[3], const double b[3], double out[3])
{
    const double x = a[1] * b[2] - a[2] * b[1];
    const double y = a[2] * b[0] - a[0] * b[2];
    const double z = a[0] * b[1] - a[1] * b[0];
    out[0] = x;
    out[1] = y;
    out[2] = z;
}

SLIMED_HD inline void v3_scale(const double a[3], double s, double out[3])
{
    out[0] = a[0] * s;
    out[1] = a[1] * s;
    out[2] = a[2] * s;
}

/**
 * @brief out = a * (1/s).
 *
 * Not a * / s. The original const_division() forms the reciprocal once and
 * multiplies, which is a different rounding from a true division, so the
 * reciprocal is kept here to preserve the numbers.
 */
SLIMED_HD inline void v3_scale_recip(const double a[3], double s, double out[3])
{
    v3_scale(a, 1.0 / s, out);
}

/// out = a x b + c x d, the a_cross_b_plus_c_cross_d() helper.
SLIMED_HD inline void v3_cross_plus_cross(const double a[3], const double b[3],
                                          const double c[3], const double d[3],
                                          double out[3])
{
    double first[3];
    double second[3];
    v3_cross(a, b, first);
    v3_cross(c, d, second);
    out[0] = first[0] + second[0];
    out[1] = first[1] + second[1];
    out[2] = first[2] + second[2];
}

/**
 * @brief out = ((a x b + c x d) * s - e * t) * u.
 *
 * The second-derivative vectors a_31, a_32, a11, a12, a21 and a22 are all this
 * shape; spelling it once keeps the six call sites readable and keeps their
 * operation order identical to each other and to the original.
 */
SLIMED_HD inline void v3_quotient_rule(const double a[3], const double b[3],
                                       const double c[3], const double d[3],
                                       double s,
                                       const double e[3], double t,
                                       double u,
                                       double out[3])
{
    double sum[3];
    v3_cross_plus_cross(a, b, c, d, sum);
    for (int i = 0; i < 3; i++)
    {
        out[i] = (sum[i] * s - e[i] * t) * u;
    }
}

/// The unit vector of m, in place. A zero vector stays zero rather than
/// becoming NaN, matching get_unit_vector().
SLIMED_HD inline void v3_normalize(double m[3])
{
    const double magnitude = v3_norm(m);
    if (magnitude == 0.0)
    {
        v3_zero(m);
        return;
    }
    v3_scale(m, 1.0 / magnitude, m);
}

/**
 * @brief The physical and constraint constants a patch evaluation needs.
 *
 * Passed by value so the struct can sit in a CUDA kernel's parameter space.
 * uSurfPerArea and uVol arrive already divided by their reference quantity,
 * with the caller responsible for the area0 == 0 and vol0 == 0 guards -- a
 * flat sheet encloses no volume, and dividing by that zero poisons every
 * vertex force with NaN.
 */
struct PatchParams
{
    double kCurv = 0.0;         ///< Bending modulus.
    double spontCurv = 0.0;     ///< Spontaneous curvature of this patch.
    double uSurfPerArea = 0.0;  ///< uSurf / area0, or 0 when area0 == 0.
    double area = 0.0;          ///< Current total membrane area.
    double area0 = 0.0;         ///< Reference (relaxed) area.
    double uVol = 0.0;          ///< uVol / vol0, or 0 when vol0 == 0.
    double vol = 0.0;           ///< Current enclosed volume.
    double vol0 = 0.0;          ///< Reference enclosed volume.
};

/**
 * @brief C = A * B for the shape-function contraction, beta = 0.
 *
 * Replaces gsl_blas_dgemm on a (kShapeRows, nCtrl) x (nCtrl, 3) product. The
 * k-outer loop with the zero skip is the reference CBLAS order GSL itself
 * uses, so the sequence of floating-point additions into each C entry is the
 * same one the original produced.
 *
 * @param[in]  rows    Shape-function block, kShapeRows x nCtrl, row-major.
 * @param[in]  ctrlPts Control points, nCtrl x 3, row-major.
 * @param[in]  nCtrl   Patch width.
 * @param[out] out     kShapeRows x 3: the limit-surface point and derivatives.
 */
SLIMED_HD inline void shape_times_control_points(const double *rows,
                                                 const double *ctrlPts,
                                                 int nCtrl,
                                                 double out[kShapeRows][3])
{
    for (int i = 0; i < kShapeRows; i++)
    {
        out[i][0] = 0.0;
        out[i][1] = 0.0;
        out[i][2] = 0.0;
    }
    for (int k = 0; k < nCtrl; k++)
    {
        for (int i = 0; i < kShapeRows; i++)
        {
            const double temp = rows[i * nCtrl + k];
            if (temp != 0.0)
            {
                out[i][0] += temp * ctrlPts[k * 3 + 0];
                out[i][1] += temp * ctrlPts[k * 3 + 1];
                out[i][2] += temp * ctrlPts[k * 3 + 2];
            }
        }
    }
}

/**
 * @brief Bending, area and volume energy and force on one triangular patch.
 *
 * Integrates over the patch's Gauss quadrature points on the Loop subdivision
 * limit surface. A regular face is evaluated once with a 12-wide net; an
 * irregular face is tiled by regular children and this is called once per
 * child, which is why the force outputs accumulate rather than overwrite.
 *
 * @param[in]     sampleRows  nSamples blocks of kShapeRows x nCtrl, row-major
 *                            and contiguous: block s starts at
 *                            sampleRows + s * kShapeRows * nCtrl.
 * @param[in]     gaussCoeff  nSamples quadrature weights.
 * @param[in]     nSamples    Number of quadrature points.
 * @param[in]     ctrlPts     nCtrl x 3 control-point coordinates, row-major.
 * @param[in]     nCtrl       Patch width; must not exceed kMaxControlPoints.
 * @param[in]     p           Physical and constraint constants.
 * @param[out]    eBend       Bending energy of this patch. Overwritten.
 * @param[out]    meanCurv    Mean curvature of this patch. Overwritten.
 * @param[in,out] normVector  Unit normal. Accumulated from zero, then
 *                            normalized in place before returning.
 * @param[in,out] fBend       nCtrl x 3 bending force. Accumulated; the caller
 *                            zeroes it before the first child.
 * @param[in,out] fArea       nCtrl x 3 area-constraint force. Accumulated.
 * @param[in,out] fVolume     nCtrl x 3 volume-constraint force. Accumulated.
 */
SLIMED_HD void element_energy_force_patch_pod(const double *sampleRows,
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
                                              double *fVolume);

/**
 * @brief Surface area and enclosed signed volume of one patch, by the same
 * quadrature the energy uses.
 *
 * The divergence theorem with F = x gives V = (1/3) * closed_integral(x . n dA),
 * and the Gauss rule sums over the reference triangle with weights summing to
 * one, so each sample carries the reference triangle's own area of 1/2. The
 * product (1/3) * (1/2) = 1/6 is kSignedVolumeQuadratureFactor. It pairs with
 * the full three-component integrand dot(x, a_1 x a_2); see
 * docs/volume_functional_split.md.
 *
 * Only the first three shape-function rows are read -- the value and the two
 * first derivatives -- so this does under half the contraction work
 * element_energy_force_patch_pod() does.
 *
 * @param[in]     sampleRows  As for element_energy_force_patch_pod().
 * @param[in]     gaussCoeff  nSamples quadrature weights.
 * @param[in]     nSamples    Number of quadrature points.
 * @param[in]     ctrlPts     nCtrl x 3 control-point coordinates, row-major.
 * @param[in]     nCtrl       Patch width.
 * @param[in,out] area        Accumulated, one quadrature sample at a time.
 *                            Irregular patches call this once per regular
 *                            child, accumulating across children too.
 * @param[in,out] volume      Accumulated the same way, and signed.
 */
SLIMED_HD void element_area_volume_pod(const double *sampleRows,
                                       const double *gaussCoeff,
                                       int nSamples,
                                       const double *ctrlPts,
                                       int nCtrl,
                                       double &area,
                                       double &volume);

} // namespace slimed
