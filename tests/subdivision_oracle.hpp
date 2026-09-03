/**
 * @file subdivision_oracle.hpp
 * @brief The literal valence-5 Loop subdivision matrices, owned by the tests.
 *
 * These were production code -- Mesh's `get_subdivision_matrices()`, a direct
 * port of the reference POC -- until WP5 replaced every caller with the
 * generated row table. They are the only independently-derived subdivision
 * data in the tree, so they are kept here verbatim as the permanent oracle the
 * generator is checked against.
 *
 * Do not regenerate these from build_subdivision_matrices(). The whole value
 * of an oracle is that it was derived some other way.
 *
 * Shapes: M(17,11), M1..M3(12,17), M4(11,17).
 *
 * @note WP2/WP5 of docs/irregular_patch_valence_4_to_8_plan.md.
 */
#pragma once

#include "linalg/Linear_algebra.hpp"

void get_subdivision_matrices_oracle(Matrix &mat,
                                     Matrix &subMat1,
                                     Matrix &subMat2,
                                     Matrix &subMat3,
                                     Matrix &subMat4);
