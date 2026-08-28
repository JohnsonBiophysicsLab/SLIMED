/**
 * @file Mesh_legacy_volume.cpp
 * @brief TEMPORARY: reproduces the pre-fix volume accumulator for one release.
 *
 * Before the dot_row fix, dot_row() bounded its loop with the row count while
 * indexing columns, so on the 1 x 3 row vectors it is called with it returned
 * the x-component product alone. Every volume SLIMED reported was that x-only
 * integrand paired with the full-divergence 1/6 factor.
 *
 * This file recomputes that historical number, for reporting only, so results
 * from earlier runs can be mapped onto the corrected scale. It is deliberately
 * isolated in its own translation unit: nothing here feeds energy, force, or
 * the volume constraint, and deleting the file plus its two call sites removes
 * the whole of it.
 *
 * @note Step 4 of docs/volume_functional_split.md. Delete after one release.
 *
 * @warning Do not promote any of this into a runtime "legacy volume" mode. The
 * plan's first non-goal is explicit: a switch between a correct functional and
 * a buggy one guarantees both live forever and some future route picks the
 * wrong default. This is a diagnostic print, not a mode.
 */

#include "mesh/Mesh.hpp"

namespace
{
/**
 * The literal the accumulator carried before it was named. Kept verbatim
 * rather than replaced with 1.0 / 6.0 so the reported number reproduces what
 * earlier runs actually printed, to the last digit.
 */
constexpr double kLegacyVolumeQuadratureLiteral = 0.16666666666;
} // namespace

double Mesh::enumerate_legacy_x_only_volume(const Matrix &matOneRingVertex)
{
    double volume = 0.0;
#pragma omp parallel for reduction(+ \
                                   : volume)
    for (int j = 0; j < param.gaussQuadratureCoeff.nrow(); j++)
    {
        Matrix &sf = param.shapeFunctions[j];
        Matrix x = sf.get_row(0) * matOneRingVertex;
        Matrix a_3 = cross_row(sf.get_row(1) * matOneRingVertex,
                               sf.get_row(2) * matOneRingVertex);
        const double coeff = param.gaussQuadratureCoeff(j, 0);

        // This is what the old dot_row(x, a_3) returned: one term, not three.
        volume += kLegacyVolumeQuadratureLiteral * coeff * x(0, 0) * a_3(0, 0);
    }
    return volume;
}

double Mesh::sum_legacy_x_only_volume()
{
    const auto &subMat = param.subMatrix;
    const Matrix &M = subMat.irregM;
    const Matrix &M1 = subMat.irregM1;
    const Matrix &M2 = subMat.irregM2;
    const Matrix &M3 = subMat.irregM3;
    const Matrix &M4 = subMat.irregM4;

    double volume = 0.0;
    for (const Face &face : faces)
    {
        if (face.isGhost)
            continue;

        // Mirrors calculate_element_area_volume(): the 12-control regular case
        // and the 11-control subdivision recursion, with no other valences.
        switch (static_cast<int>(face.oneRingVertices.size()))
        {
        case 12:
            volume += enumerate_legacy_x_only_volume(get_one_ring_vertex_matrix(face));
            break;

        case 11:
        {
            Matrix matOrigOneRingVertex = get_one_ring_vertex_matrix(face);
            Matrix matNewNodes17;
            for (int j = 0; j < param.subDivideTimes; j++)
            {
                matNewNodes17 = M * matOrigOneRingVertex;
                volume += enumerate_legacy_x_only_volume(M1 * matNewNodes17);
                volume += enumerate_legacy_x_only_volume(M2 * matNewNodes17);
                volume += enumerate_legacy_x_only_volume(M3 * matNewNodes17);
                matOrigOneRingVertex = M4 * matNewNodes17;
            }
        }
        break;
        }
    }
    return volume;
}

void Mesh::report_volume_rebaseline()
{
    const double corrected = param.vol0;
    const double legacy = sum_legacy_x_only_volume();

    std::cout << "[Volume rebaseline] corrected vol0 = " << corrected
              << " ; legacy (pre-fix, x-only) vol0 = " << legacy << std::endl;

    if (legacy != 0.0)
    {
        std::cout << "[Volume rebaseline] corrected / legacy = " << corrected / legacy
                  << std::endl;
        std::cout << "[Volume rebaseline] On a closed surface this ratio is exactly 3 and "
                     "converts any vol0 or volume recorded before the fix. On an open "
                     "surface neither value is a volume, the ratio is configuration "
                     "dependent, and it does NOT convert a whole trajectory."
                  << std::endl;
    }

    std::cout << "[Volume rebaseline] This reporting is temporary and will be removed. "
                 "See docs/volume_functional_split.md."
              << std::endl;
}
