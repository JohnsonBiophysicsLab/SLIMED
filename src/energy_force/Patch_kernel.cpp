/**
 * @file Patch_kernel.cpp
 * @brief Holds the patch kernel's constants to the mesh's.
 *
 * The kernel itself is header-only -- see the note at the foot of
 * Patch_kernel.hpp. What needs a translation unit is the agreement between the
 * device-side constants and the mesh-side originals, which is checked below.
 */
#include "energy_force/Patch_kernel.hpp"

#include "mesh/Irregular_patch_rows.hpp"
#include "mesh/Subdivision_matrices.hpp"

// Patch_kernel.hpp cannot include the headers these come from -- they pull in
// GSL, which has no place in code nvcc compiles for a device -- so it carries
// its own copies. This is the one host translation unit that sees both, and so
// the one place the two can be held together.
static_assert(slimed::kMinPatchValence == kMinIrregularValence,
              "slimed::kMinPatchValence has drifted from kMinIrregularValence");
static_assert(slimed::kMaxPatchValence == kMaxIrregularValence,
              "slimed::kMaxPatchValence has drifted from kMaxIrregularValence");
static_assert(slimed::kChildrenPerSubdivisionStep == kRegularChildrenPerStep,
              "slimed::kChildrenPerSubdivisionStep has drifted from kRegularChildrenPerStep");
static_assert(kMaxIrregularValence + 6 <= slimed::kMaxControlPoints,
              "the widest patch no longer fits the kernel's stack buffers");
