#include "Dynamics.hpp"

using namespace std;

// calculate mesh2surface matrix
// surface = mesh2surface * mesh

DynamicMesh::DynamicMesh(Param &srcParam) :
Mesh(srcParam)
{
}

void DynamicMesh::setup_flat() {
    // Call the superclass's setup_flat()
    Mesh::setup_flat();

    // The dense path stores M^-1, computed once. An edge flip changes four
    // rows of M, and there is no way to update a stored inverse for that short
    // of recomputing it -- O(N^3), per accepted flip. A fluid run configured
    // this way would either be unusably slow or, worse, quietly keep using an
    // inverse that no longer describes its mesh.
    if (param.edgeFlipEnabled && param.surfaceSolver != "iterative")
    {
        throw std::runtime_error(
            "[DynamicMesh::setup_flat] edgeFlipEnabled = true needs surfaceSolver = iterative. "
            "The dense conversion stores an inverse of the limit mask, and an edge flip changes "
            "the mask it was built from. See docs/edge_flip_plan.md section 3.7.");
    }

    // A flip on a triangulation of valence-6 vertices leaves two faces with
    // extraordinary corners at both ends of the new edge, so a fluid mesh is
    // full of faces with more than one of them -- and DeviceMeshLayout cannot
    // build those yet. Left alone, the run would start, flip, and then throw
    // out of the layout builder partway through, with a mesh already changed.
    // Say so before the first step instead.
    if (param.edgeFlipEnabled && param.forceBackend == "gpu")
    {
        throw std::runtime_error(
            "[DynamicMesh::setup_flat] edgeFlipEnabled = true needs forceBackend = cpu. A flip "
            "creates faces with more than one extraordinary corner, which the device layout "
            "cannot represent; only the CPU kernel evaluates those. See "
            "docs/edge_flip_plan.md work package 5.");
    }

    if (param.edgeFlipEnabled)
    {
        report_edge_flip_feasibility();
    }

    // Assign mesh2surface and surface2mesh
    //
    // Skipped entirely under the iterative solver, which is the whole point of
    // it: the dense path allocates two N x N matrices and inverts one of them,
    // which on a 300 nm sheet (4331 vertices) is 55 seconds and 300 MB before
    // the first step runs. The sparse structure it is replaced by is built
    // lazily by ensure_surface_solver(), in a quarter of a millisecond.
    if (param.surfaceSolver != "iterative")
    {
        if (param.VERBOSE_MODE)
        {
            std::cout << "[DynamicMesh] Assigning conversion matrices." << std::endl;
        }
        assign_mesh2surface();
        if (param.VERBOSE_MODE)
        {
            std::cout << "[DynamicMesh] mesh2surface matrix : " << std::endl
                      << mesh2surface << std::endl;
        }
        mesh2surface.get_inverted(surface2mesh);
        if (param.VERBOSE_MODE)
        {
            std::cout << "[DynamicMesh] surface2mesh matrix : " << std::endl
                      << surface2mesh << std::endl;
        }
    }
    matMesh = mat_calloc(vertices.size(), 3);
    matSurface = mat_calloc(vertices.size(), 3);
    mark_slaved_periodic_vertices();
}

void DynamicMesh::report_edge_flip_feasibility()
{
    // A flip on a rhombus of two equilateral triangles of side l replaces the
    // short diagonal l by the long one, l * sqrt(3). With a harmonic tether at
    // rest length l the move therefore has to climb
    //
    //     dE = (k / 2) (sqrt(3) - 1)^2 l^2
    //
    // and Metropolis accepts it with probability exp(-dE / kT). That barrier
    // grows with the square of the mesh spacing, so a stiffness that is
    // reasonable at one lFace can freeze the membrane solid at another, and
    // nothing downstream says so: the run proceeds, reports its attempts, and
    // accepts none of them. This is why the dynamically triangulated surface
    // literature uses a flat-bottomed tether rather than a spring -- inside
    // the allowed range a flip costs nothing.
    if (param.edgeSpringEnabled)
    {
        const double restLength =
            (param.edgeSpringRestLength > 0.0) ? param.edgeSpringRestLength : param.lFace;
        const double reach = std::sqrt(3.0) - 1.0;
        const double barrier = 0.5 * param.edgeSpringConstant * reach * reach * restLength *
                               restLength;
        const double inKT = (param.KBT > 0.0) ? barrier / param.KBT : 0.0;
        std::cout << "[DynamicMesh] Edge-flip barrier from the edge spring: " << barrier
                  << " pN.nm = " << inKT << " kT (k = " << param.edgeSpringConstant
                  << " pN/nm, l0 = " << restLength << " nm)." << std::endl;
        if (inKT > 10.0)
        {
            std::cout << "[DynamicMesh] WARNING: at " << inKT
                      << " kT that barrier accepts roughly exp(-" << inKT
                      << ") of the flips offered, so the membrane will not be fluid. The "
                         "acceptance is what edgeFlipAttemptRate is calibrated against, and it "
                         "will read as zero however high the rate is set. Lower "
                         "edgeSpringConstant -- a few kT of barrier means k of order "
                      << 10.0 * param.KBT / (0.5 * reach * reach * restLength * restLength)
                      << " pN/nm at this lFace." << std::endl;
        }
    }
    else
    {
        std::cout
            << "[DynamicMesh] WARNING: edgeFlipEnabled = true with edgeSpringEnabled = false. "
               "The regularization then in force measures each face against its own edges in "
               "the reference configuration, which an edge a flip has just created never had. "
               "Those faces carry an energy that is not a function of the mesh, and the "
               "Metropolis sweep samples it: measured on a 100 nm sheet, every accepted flip "
               "reported about -1200 pN.nm of it. Set edgeSpringEnabled = true. See "
               "docs/edge_flip_plan.md section 3.7."
            << std::endl;
    }
}

void DynamicMesh::mark_slaved_periodic_vertices()
{
    isSlavedPeriodic.assign(vertices.size(), 0);
    if (param.boundaryCondition != BoundaryType::Periodic)
    {
        return;
    }
    int nSlaved = 0;
    for (int i = 0; i < static_cast<int>(vertices.size()); i++)
    {
        if (vertices[i].isGhost)
        {
            continue;
        }
        if (get_relative_pt_periodic(i, param.nFaceX, param.nFaceY) != 0)
        {
            isSlavedPeriodic[i] = 1;
            nSlaved++;
        }
    }
    std::cout << "[DynamicMesh::mark_slaved_periodic_vertices] "
              << nSlaved << " of " << vertices.size()
              << " vertices are periodic duplicates; they are "
              << (param.integratePeriodicDuplicates ? "integrated anyway (legacy)"
                                                    : "left to their partners")
              << "." << std::endl;

    // A periodic duplicate is not an independent coordinate:
    // postprocess_ghost_periodic() overwrites it from its partner every step.
    // Flipping an edge around one would change its connectivity while its
    // partner's stayed put, so the two would no longer describe the same
    // patch. Freeze them for the flip move, and re-settle every edge that
    // touches one.
    flipFrozenVertex.assign(vertices.size(), 0);
    for (int i = 0; i < static_cast<int>(vertices.size()); i++)
    {
        flipFrozenVertex[i] = isSlavedPeriodic.empty() ? 0 : isSlavedPeriodic[i];
    }
    for (int iEdge = 0; iEdge < static_cast<int>(edges.size()); iEdge++)
    {
        refresh_edge_flippability(iEdge);
    }
}

void DynamicMesh::assign_mesh2surface()
{
    // Initialize mesh2surface matrix
    mesh2surface = mat_calloc(vertices.size(), vertices.size());
    
    // Iterate through all vertices
    for (int i = 0; i < vertices.size(); i++)
    {
        // find adjacent face that is not ghost
        int indexAdjFace = -1; // -1 means only ghost
        for (int jAdjFace = 0; jAdjFace < vertices[i].adjacentFaces.size(); jAdjFace++)
        {
            const Face& adjF = faces[vertices[i].adjacentFaces[jAdjFace]];
            if (!adjF.isGhost && !adjF.isBoundary)
            {
                // assign non ghost / boundary face index to indexAdjFace
                indexAdjFace = vertices[i].adjacentFaces[jAdjFace];
            }
        }
        if (indexAdjFace < 0)
        {
            // If only ghost or boundary faces are adjacent to vertex, shape function = (1, 0, 0, ...)
            mesh2surface.set(i, i, 1.0);
        }
        else
        {
            // If there exists non-ghost / boundary face
            // find vwu of vertex on face
            //! vwu assumed to be same sequence as defined in adjacent face
            Face &faceAdj = faces[indexAdjFace];
            std::vector<int> &adjVertices = faceAdj.adjacentVertices;
            int vwuInd = -1;
            for (int j = 0; j < adjVertices.size(); j++)
            {
                if ((adjVertices)[j] == i)
                {
                    vwuInd = j;
                }
            }
            Matrix vwuAdj = mat_calloc(1, 3); // get_shapefunction takes row vector
            vwuAdj.set(0, vwuInd, 1.0);
            // std::cout<< i << " ," << vwuAdj << endl;

            // use vwu to get sf - irregular patch -currently not used!
            // @todo implement support for irregular patch
            // Matrix sf(7, 12);
            // get_shapefunction(vwuAdj, sf);
            // vector<double> sfAdj = determine_ShapeFunctions(vwuAdj)[0]; // transposed sf

            // set mat(AB) (mesh2surface) to sf
            //@TODO: currently due to vwu not in order in adjacent vertex of face
            // using 0.5 / 0.0833333 directly for regular patches
            std::vector<int> &iAdjVertices = vertices[i].adjacentVertices;

            // Loop's limit mask: half on the vertex, the other half shared
            // equally among its neighbours. The denominator used to be a
            // literal 6 whatever the vertex's valence was, which is right on a
            // regular mesh and wrong everywhere else -- and a fluid membrane
            // is mostly not valence 6. On the workloads this tree has run,
            // every vertex that reaches here is at valence 6, so the
            // correction changes nothing there.
            const int valence = static_cast<int>(iAdjVertices.size());
            mesh2surface.set(i, i, 0.5);
            if (valence > 0)
            {
                for (const int &iAdj : iAdjVertices)
                {
                    mesh2surface.set(i, iAdj, 0.5 / valence);
                }
            }
        }
        // check if ghost / boundary faces are correctly recongnized
        // std::cout << "vertex: " << i << ", indexAdjFace: " << indexAdjFace << endl;
    }
    // The Brownian update maps the nodal force onto the limit-surface
    // coordinates with M^-1 in place of M^-T, which is only the same thing
    // when M is symmetric.  It very nearly is: the interior rows are the
    // valence-6 limit mask, which is symmetric.  It is not exactly, because a
    // vertex whose adjacent faces are all ghost or boundary gets the identity
    // row above while its neighbours' rows still reference it.  Report the
    // asymmetry rather than assume it away -- a non-reciprocal mobility does
    // not just bias the sampled distribution, it can do net work on the
    // membrane, and that shows up as an amplitude that never settles.
    double maxAsymmetry = 0.0;
    for (int i = 0; i < static_cast<int>(vertices.size()); i++)
    {
        for (int j = i + 1; j < static_cast<int>(vertices.size()); j++)
        {
            maxAsymmetry = std::max(maxAsymmetry,
                                    std::abs(mesh2surface(i, j) - mesh2surface(j, i)));
        }
    }
    std::cout << "[DynamicMesh::assign_mesh2surface] mesh2surface built; "
              << "max |M - M^T| = " << maxAsymmetry << std::endl;

    if (param.VERBOSE_MODE)
    {
        std::cout << mesh2surface << std::endl;
    }
}

void DynamicMesh::ensure_surface_solver()
{
    if (surfaceSolver.empty() || surfaceSolver.topologyVersion != topologyVersion)
    {
        surfaceSolver.build(*this);
        if (param.VERBOSE_MODE)
        {
            std::cout << "[DynamicMesh::ensure_surface_solver] built for topology version "
                      << topologyVersion << ": " << surfaceSolver.nFree() << " of "
                      << surfaceSolver.nVertices() << " vertices are degrees of freedom."
                      << std::endl;
        }
    }
}

void DynamicMesh::apply_mesh_to_surface()
{
    if (param.surfaceSolver == "iterative")
    {
        ensure_surface_solver();
        surfaceSolver.mesh_to_surface(matMesh, matSurface);
        return;
    }
    matSurface = mesh2surface * matMesh;
}

void DynamicMesh::apply_surface_to_mesh()
{
    if (param.surfaceSolver == "iterative")
    {
        ensure_surface_solver();
        // The previous step's control net is a very good starting point --
        // the step moved the surface by a Brownian increment, not by much --
        // and warm starting cuts the iteration count severalfold.
        surfaceSolver.surface_to_mesh(matSurface, matMesh, &matMesh);
        return;
    }
    matMesh = surface2mesh * matSurface;
}

void DynamicMesh::apply_nodal_force_to_surface(const Matrix &nodalForce, Matrix &surfaceForce)
{
    if (param.surfaceSolver == "iterative")
    {
        ensure_surface_solver();
        surfaceSolver.nodal_force_to_surface(nodalForce, surfaceForce);
        return;
    }
    // The historical path: M^-1 where M^-T belongs. Identical wherever the
    // mask is symmetric, which is every interior valence-6 row; the difference
    // lives in the ghost band, whose vertices the step does not integrate.
    surfaceForce = surface2mesh * nodalForce;
}

void DynamicMesh::update_vertices_mat_with_vector()
{
    for (int i = 0; i < vertices.size(); i++){
        for (int j = 0; j < 3; j++){
            matMesh.set(i, j, vertices[i].coord(j,0));
        }
    }
}

void DynamicMesh::update_vertices_vector_with_mat()
{
    for (int i = 0; i < vertices.size(); i++){
        for (int j = 0; j < 3; j++){
            vertices[i].coord.set(j, 0, matMesh(i, j));
        }
    }
}

