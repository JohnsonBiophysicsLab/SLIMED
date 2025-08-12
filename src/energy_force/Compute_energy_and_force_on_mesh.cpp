#include "mesh/Mesh.hpp"
#pragma omp declare reduction(vector_plus                                                                                             \
                              : std::vector <Force>                                                                                   \
                              : std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus <Force>())) \
    initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))

using namespace std;

/**
 *
 * @brief Computes the energy and force on each vertex and face of the mesh.
 *
 * This function performs the following steps:
 *       Calculates the area and volume of each element triangle and sums up the total area and volume of the membrane.
 *       Iterates through faces and calculates the energy and force on each triangular patch.
 *       Regularizes the force and energy.
 *       Sums up the energy and force on each vertex and face.
 *       Calculates the energy due to area constraint.
 *       Calculates the energy due to volume constraint.
 *       Calculates the energy due to scaffolding.
 *
 * @return void
 */
void Mesh::Compute_Energy_And_Force()
{

    // Step 1.
    // Calculate the area and volume of each element triangle
    calculate_element_area_volume();

    // Sum up the total area and volume of the membrane and sync with parameter
    sum_membrane_area_and_volume(param.area, param.vol);

    // Reset the force on vertices and energy on faces
    clear_force_on_vertices_and_energy_on_faces();

    // Initialize forces to zeros
    for (Vertex &vertex : vertices)
    {
        vertex.force = Force();
    }

    // Step 2.
    // Iterate through faces and calculate forces
#pragma omp parallel for
    for (Face &face : faces)
    {
        // std::cout << "Face index: " << i << endl;
        if (face.isBoundary)
            continue;

        // Get number of one ring vertices
        int nOneRingVertices = face.oneRingVertices.size();
//cout << "CEAF 54" << endl;
        // Get coord of one ring vertices
        std::vector<Matrix> coordOneRingVertices(nOneRingVertices);
        std::transform(face.oneRingVertices.begin(), face.oneRingVertices.end(),
                       coordOneRingVertices.begin(),
                       [this](int iVertex)
                       { return vertices[iVertex].coord; });

        double spontCurv = face.spontCurvature;         // spontaneous curvature of each patch,
        double eBend = 0.0;                             // curvature Energy of this element;
        double meanCurv = 0.0;                          // mean curvature of this element;
        Matrix fBend = mat_calloc(nOneRingVertices, 3); // bending or curvature term
        Matrix fArea = mat_calloc(nOneRingVertices, 3); // area term
        Matrix fVol = mat_calloc(nOneRingVertices, 3);  // volume term
        face.normVector.free();                         // reinitialize empty normal vector
        face.normVector = mat_calloc(3, 1);             // normal vector
                                                        // Calculate energy and force on the given triangular patch
//cout << "CEAF 71" << endl;
        // regular patch
        if (nOneRingVertices == 12)
        {
//cout << "CEAF 75" << endl;
            element_energy_force_regular(coordOneRingVertices,
                                         face,
                                         spontCurv,
                                         meanCurv,
                                         face.normVector,
                                         eBend,
                                         fBend,
                                         fArea,
                                         fVol);

            // irregular patch
        }
        else if (nOneRingVertices == 11)
        {
            //@todo energy force irregular
            // element_energy_force_irregular(coordOneRingVertices, param, spontCurv, meanCurv, normVector, eBend, fBend, fArea, fVol, GaussQuadratureCoeff, ShapeFunctions, subMatrix);
            element_energy_force_regular(coordOneRingVertices,
                                         face,
                                         spontCurv,
                                         meanCurv,
                                         face.normVector,
                                         eBend,
                                         fBend,
                                         fArea,
                                         fVol);
        }
        face.energy.energyCurvature = eBend; ///< store curvature energy in face object

        
        
        face.meanCurvature = meanCurv;       ///< store mean curvature in face object
        /*
        if (eBend > 0.1){

            std::cout << "energyCurvature: " << face.index << " = " << face.energy.energyCurvature << std::endl;
            std::cout << "meanCurvature: " << face.index << " = " << face.meanCurvature << std::endl;
            std::cout << "area: " << face.index << " = " << face.elementArea << std::endl;
            double estimate_ebend = 2.0 * 83.4 * face.elementArea * meanCurv * meanCurv;
            std::cout << "estimateEbend = " << estimate_ebend << std::endl;
        }*/

        for (int j = 0; j < nOneRingVertices; j++)
        {
            int iVertex = face.oneRingVertices[j];
            Force forceTmp;
            forceTmp.forceCurvature = fBend.get_row(j).get_transposed(); ///< store curvature force
            forceTmp.forceArea = fArea.get_row(j).get_transposed();      ///< store area force
            forceTmp.forceVolume = fVol.get_row(j).get_transposed();     ///< store volume force
            forceTmp.calculate_total_force();
            this->vertices[iVertex].force += forceTmp; ///< add the calculated forces to new_forces array
            // std::cout << "Force at node: " << fBend[j][0] << fBend[j][1]<< fBend[j][2] << F_consA[j][0] << F_consV[j][0] << endl;
        }

        // Free allocation
        fBend.free();
        fArea.free();
        fVol.free();

    }

    // Step 3.
    energy_force_regularization(); // regularization Force and Energy

    // Step 4.
    // Summing up energy and force
    // On each vertex, the total Force is the sum of internal, constraint and regularization Force.
#pragma omp parallel for
    for (Vertex &vertex : vertices)
    {
        vertex.force.calculate_total_force();
        //std::cout << "F@v = " << vertex.index << "," << vertex.force.forceTotal.get(2,0)<< std::endl;
    }

    // On each face, the total Energy is the sum of internal, constraint and regularization Energy.
    Energy sumEnergy;
//#pragma omp parallel for
    for (Face &face : faces)
    {
        face.energy.calculateTotalEnergy();
        sumEnergy += face.energy;
    }

    // Step 5. Calculate area constraint energy
    if (param.isGlobalConstraint)
    {
        if (param.area0 == 0.0)
        {
            sumEnergy.energyArea = 0.0;
        }
        else
        {
            sumEnergy.energyArea = 0.5 * param.uSurf / param.area0 * pow(param.area - param.area0, 2.0);
        }
    }

    // Step 6. Calculate volume constraint energy - volume energy is always global
    if (param.vol0 == 0.0)
    {
        sumEnergy.energyVolume = 0.0;
    }
    else
    {
        sumEnergy.energyVolume = 0.5 * param.uVol / param.vol0 * pow(param.vol - param.vol0, 2.0);
    }
    sumEnergy.calculateTotalEnergy();
    param.energy = sumEnergy;

    // Step 7. Calculate scaffolding energy
    if (param.isEnergyHarmonicBondIncluded)
    {
        param.energy.energyHarmonicBond = this->calculate_scaffolding_energy_force(false);
        param.energy.calculateTotalEnergy();
    }

    // Step 8. Boundary conditions - manage ghost vertices
    manage_force_for_boundary_ghost_vertex();

    if (false && param.VERBOSE_MODE)
    {
        for (int i = 0; i < vertices.size(); i++)
        {
            std::cout << vertices[i].force.to_string_full_at_vertex(i) << std::endl;
        }
    }
}

void Mesh::element_energy_force_regular(const std::vector<Matrix> &coordOneRingVertices,
                                        Face& face,
                                        const double spontCurv,
                                        double &meanCurv,
                                        Matrix &normVector,
                                        double &eBend,
                                        Matrix &fBend,
                                        Matrix &fArea,
                                        Matrix &fVolume)
{
    // fBend is the Force related to the curvature
    // fArea is the Force related to the area-constraint
    // fVolume is the Force related to the volume-constraint.
    // All these forces must be 0 as input.
    // Normal vector must be 0 as input
    // initialize output parameters
//cout << "EEFR 202" << endl;
    eBend = 0.0;
    meanCurv = 0.0;
    double halfGaussQuadratureCoeff = 0.0;
    //////////////////////////////////////////////////////////////
    double kCurv = param.kCurv; //test output
    double uSurfPerArea = param.uSurf / param.area0;
    int nFaces = this->faces.size();
    double area0 = param.area0;
    double area = param.area;
    double uVol = param.uVol / param.vol0;
    double vol0 = param.vol0;
    double vol = param.vol;
    /////////////////////////////////////////////////////////////
    // double definition
    double sqa = 0.0;
    double sqa_1 = 0.0;
    double sqa_2 = 0.0;

    double H_curv = 0.0;
    double eBend_tmp = 0.0;
    // Matrix definition

    Matrix sfDotOneRingV(7, 3);
    Matrix x(3, 1);
    Matrix a_1(3, 1);
    Matrix a_2(3, 1);
    Matrix a_11(3, 1);
    Matrix a_22(3, 1);
    Matrix a_12(3, 1);
    Matrix a_21(3, 1);
    Matrix xa(3, 1);

    Matrix xa_1(3, 1);
    Matrix xa_2(3, 1);

    Matrix a_3(3, 1);
    Matrix a_31(3, 1);
    Matrix a_32(3, 1);
    // Matrix d(3, 1);
    // Matrix d_1(3, 1);
    // Matrix d_2(3, 1);
    Matrix a1(3, 1);
    Matrix a2(3, 1);
    Matrix a11(3, 1);
    Matrix a12(3, 1);
    Matrix a21(3, 1);
    Matrix a22(3, 1);

    Matrix n1_be(3, 1);
    Matrix n2_be(3, 1);
    Matrix m1_be(3, 1);
    Matrix m2_be(3, 1);
    Matrix n1_cons(3, 1);
    Matrix n2_cons(3, 1);
    Matrix n1_conv(3, 1);
    Matrix n2_conv(3, 1);

    // 2d vectors
    Matrix f_be(12, 3);
    Matrix f_cons(12, 3);
    Matrix f_conv(12, 3);
    Matrix da1(3, 3);
    Matrix da2(3, 3);

    // TEMP vectors / double
    Matrix tmp_f(3, 1);
    Matrix tmp_l(3, 1);
    Matrix tmp_sum(3, 1);
    Matrix tmp_a2x3(3, 1);
    Matrix tmp_a3x1(3, 1);
    Matrix tempf_be(3, 1);
    Matrix tempf_cons(3, 1);
    Matrix tempf_conv(3, 1);
    Matrix tmp_kron_prod(3, 3);

    double tmp_const_f = 0.0;
    double tmp_const_l = 0.0;
    double tmp_evol = 0.0;
    double tmp_sqa_sqrinv = 0.0;
////cout << "EEFR 279" << endl;
    // convert vector of matrices to matrix
    Matrix matOneRingVertices(coordOneRingVertices.size(), 3);
    for (int i = 0; i < coordOneRingVertices.size(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            matOneRingVertices.set(i, j, coordOneRingVertices[i](j, 0));
        }
    }

    // Gaussian quadrature, second-order or 3 points.
    //std::cout << param.shapeFunctions.size() << std::endl;
    for (int i = 0; i < param.shapeFunctions.size(); i++)
    {
        halfGaussQuadratureCoeff = 0.5 * param.gaussQuadratureCoeff(i, 0);
        //std::cout << halfGaussQuadratureCoeff << std::endl;
        Matrix &sf = param.shapeFunctions[i];

        multiplication(sf, matOneRingVertices, sfDotOneRingV); //< (7, 12) . (12, 3) = (7, 3)
////cout << "EEFR 297" << endl;
        /*
         * The (7, 3) matrix represents the coordinate and derivatives of the point
         * on the limit surface:
         *   x : point on the limit surface; shapefunction(0,:), shape functions;
         *   a_1 : dx/du; shapefunction(1,:), differential to v;
         *   a_2 : dx/dv; shapefunction(2,:), differential to w;
         *   a_11 : d^2x/du^2; shapefunction(3,:), double differential to v;
         *   a_22 : d^2x/dv^2; shapefunction(4,:), double differential to w;
         *   a_12 : d^2x/dudv; shapefunction(5,:), differential to v and w;
         *   a_21 : d^2x/dvdu; shapefunction(6,:), differential to w and v;
         */

        assign_rowVec_to_colVec(sfDotOneRingV.get_row(0), x);
        assign_rowVec_to_colVec(sfDotOneRingV.get_row(1), a_1);
        assign_rowVec_to_colVec(sfDotOneRingV.get_row(2), a_2);
        assign_rowVec_to_colVec(sfDotOneRingV.get_row(3), a_11);
        assign_rowVec_to_colVec(sfDotOneRingV.get_row(4), a_22);
        assign_rowVec_to_colVec(sfDotOneRingV.get_row(5), a_12);
        assign_rowVec_to_colVec(sfDotOneRingV.get_row(6), a_21);

        cross(a_1, a_2, xa); // a_1 x a_2, perpendicular to the tangent plane
        sqa = xa.calculate_norm();
        tmp_sqa_sqrinv = 1.0 / sqa / sqa; // square inverse of sqa

        // Matrix xa_1 = cross(a_11,a_2) + cross(a_1,a_21);
        a_cross_b_plus_c_cross_d(a_11, a_2, a_1, a_21, tmp_f, tmp_l, xa_1);
        // Matrix xa_2 = cross(a_12,a_2) + cross(a_1,a_22);
        a_cross_b_plus_c_cross_d(a_12, a_2, a_1, a_22, tmp_f, tmp_l, xa_2);

        sqa_1 = dot_col(xa, xa_1) / sqa;
        sqa_2 = dot_col(xa, xa_2) / sqa;
        // Matrix a_3 = xa/sqa;
        const_division(xa, sqa, a_3);
        // Matrix a_31 = 1.0/sqa/sqa * (xa_1*sqa - xa*sqa_1);
        const_multiplication(xa_1, sqa, tmp_f);
        const_multiplication(xa, sqa_1, tmp_l);
        subtraction(tmp_f, tmp_l, tmp_sum);
        const_multiplication(tmp_sum, tmp_sqa_sqrinv, a_31);
        // Matrix a_32 = 1.0/sqa/sqa * (xa_2*sqa - xa*sqa_2);
        const_multiplication(xa_2, sqa, tmp_f);
        const_multiplication(xa, sqa_2, tmp_l);
        subtraction(tmp_f, tmp_l, tmp_sum);
        const_multiplication(tmp_sum, tmp_sqa_sqrinv, a_32);
////cout << "EEFR 341" << endl;
        // d = a_3;
        // d_1 = a_31;
        // d_2 = a_32;

        // tmp matrices to store scross result
        cross(a_2, a_3, tmp_a2x3);
        cross(a_3, a_1, tmp_a3x1);

        // Matrix a1 = cross(a_2,a_3)/sqa;
        const_division(tmp_a2x3, sqa, a1);

        // Matrix a2 = cross(a_3,a_1)/sqa;
        const_division(tmp_a3x1, sqa, a2);

        // Matrix a11 = 1.0/sqa/sqa *( (cross(a_21,a_3)+cross(a_2,a_31))*sqa - cross(a_2,a_3)*sqa_1 );
        a_cross_b_plus_c_cross_d(a_21, a_3, a_2, a_31, tmp_f, tmp_l, tmp_sum);
        const_multiplication(tmp_a2x3, sqa_1, tmp_l);
        const_multiplication(tmp_sum, sqa, tmp_f);
        subtraction(tmp_f, tmp_l, tmp_sum);
        const_multiplication(tmp_sum, tmp_sqa_sqrinv, a11);

        // Matrix a12 = 1.0/sqa/sqa *( (cross(a_22,a_3)+cross(a_2,a_32))*sqa - cross(a_2,a_3)*sqa_2 );
        a_cross_b_plus_c_cross_d(a_22, a_3, a_2, a_32, tmp_f, tmp_l, tmp_sum);
        const_multiplication(tmp_a2x3, sqa_2, tmp_l);
        const_multiplication(tmp_sum, sqa, tmp_f);
        subtraction(tmp_f, tmp_l, tmp_sum);
        const_multiplication(tmp_sum, tmp_sqa_sqrinv, a12);

        // Matrix a21 = 1.0/sqa/sqa *( (cross(a_31,a_1)+cross(a_3,a_11))*sqa - cross(a_3,a_1)*sqa_1 );
        a_cross_b_plus_c_cross_d(a_31, a_1, a_3, a_11, tmp_f, tmp_l, tmp_sum);
        const_multiplication(tmp_a3x1, sqa_1, tmp_l);
        const_multiplication(tmp_sum, sqa, tmp_f);
        subtraction(tmp_f, tmp_l, tmp_sum);
        const_multiplication(tmp_sum, tmp_sqa_sqrinv, a21);

        // Matrix a22 = 1.0/sqa/sqa *( (cross(a_32,a_1)+cross(a_3,a_12))*sqa - cross(a_3,a_1)*sqa_2 );
        a_cross_b_plus_c_cross_d(a_32, a_1, a_3, a_12, tmp_f, tmp_l, tmp_sum);
        const_multiplication(tmp_a3x1, sqa_2, tmp_l);
        const_multiplication(tmp_sum, sqa, tmp_f);
        subtraction(tmp_f, tmp_l, tmp_sum);
        const_multiplication(tmp_sum, 1.0 / sqa / sqa, a22);

        H_curv = 0.5 * (dot_col(a1, a_31) + dot_col(a2, a_32));

        // Matrix n1_be = -kCurv*(2.0*H_curv-spontCurv)*(dot_col(a1,a1)*d_1 + dot_col(a1,a2)*d_2) + kCurv*0.5*pow(2.0*H_curv-spontCurv,2.0)*a1;
        tmp_const_f = -kCurv * (2.0 * H_curv - spontCurv);
        tmp_const_l = kCurv * 0.5 * pow(2.0 * H_curv - spontCurv, 2);
        const_multiplication(a_31, dot_col(a1, a1), tmp_f);
        const_multiplication(a_32, dot_col(a1, a2), tmp_l);
        addition(tmp_f, tmp_l, tmp_sum);
        const_multiplication(tmp_sum, tmp_const_f, tmp_f);
        const_multiplication(a1, tmp_const_l, tmp_l);
        addition(tmp_f, tmp_l, n1_be);

        // Matrix n2_be = -kCurv*(2.0*H_curv-spontCurv)*(dot_col(a2,a1)*d_1 + dot_col(a2,a2)*d_2) + kCurv*0.5*pow(2.0*H_curv-spontCurv,2.0)*a2;
        const_multiplication(a_31, dot_col(a2, a1), tmp_f);
        const_multiplication(a_32, dot_col(a2, a2), tmp_l);
        addition(tmp_f, tmp_l, tmp_sum);
        const_multiplication(tmp_sum, tmp_const_f, tmp_f);
        const_multiplication(a2, tmp_const_l, tmp_l);
        addition(tmp_f, tmp_l, n2_be);

        // Matrix m1_be = kCurv*(2.0*H_curv-spontCurv)*a1;
        const_multiplication(a1, -tmp_const_f, m1_be);

        // Matrix m2_be = kCurv*(2.0*H_curv-spontCurv)*a2;
        const_multiplication(a2, -tmp_const_f, m2_be);

        // Matrix n1_cons = uSurfPerArea * (area-area0)*a1;
        // Matrix n2_cons = uSurfPerArea * (area-area0)*a2;

        if (param.isGlobalConstraint){
            tmp_const_f = uSurfPerArea * (area - area0); // Global mode
        } else {
            tmp_const_f = 0.5 * uSurfPerArea * nFaces * (nFaces * face.elementArea / area0 - 1); // Local mode
        }

        tmp_const_f = uSurfPerArea * (area - area0);
        const_multiplication(a1, tmp_const_f, n1_cons);
        const_multiplication(a2, tmp_const_f, n2_cons);

        // Matrix n1_conv = 1.0/3.0*uVol*(vol-vol0)*(dot_col(x,d)*a1 - dot_col(x,a1)*d);
        // Matrix n2_conv = 1.0/3.0*uVol*(vol-vol0)*(dot_col(x,d)*a2 - dot_col(x,a2)*d);
        tmp_evol = uVol * (vol - vol0) / 3.0;
        const_multiplication(a1, dot_col(x, a_3), tmp_f);
        const_multiplication(a_3, dot_col(x, a1), tmp_l);
        subtraction(tmp_f, tmp_l, tmp_sum);
        const_multiplication(tmp_sum, tmp_evol, n1_conv);

        const_multiplication(a2, dot_col(x, a_3), tmp_f);
        const_multiplication(a_3, dot_col(x, a2), tmp_l);
        subtraction(tmp_f, tmp_l, tmp_sum);
        const_multiplication(tmp_sum, tmp_evol, n2_conv);
////cout << "EEFR 429" << endl;
        //std::cout << sqa << std::endl;
        eBend_tmp = 0.5 * kCurv * sqa * pow(2.0 * H_curv - spontCurv, 2); // bending Energy
        // std::cout << "(2H, spontCurv, 2h-spontCurv, ebe)" << 2.0*H_curv << ", " << spontCurv << ", " << (2.0*H_curv-spontCurv) <<", " << ebe << endl;
        for (int j = 0; j < 12; j++)
        {   
            //da1 = -sf(3, j) * kron(a1, a_3)
            //- sf(1, j) * kron(a11, a_3)
            //- sf(1, j) * kron(a1, a_31)
            //- sf(6, j) * kron(a2, a_3)
            //- sf(2, j) * kron(a21, a_3)
            //- sf(2, j) * kron(a2, a_31);
            da1.set_all(0.0);
            kron(a1, a_3, tmp_kron_prod);
            tmp_kron_prod *= -sf(3,j);
            da1 += tmp_kron_prod;
            kron(a11, a_3, tmp_kron_prod);
            tmp_kron_prod *= -sf(1,j);
            da1 += tmp_kron_prod;
            kron(a1, a_31, tmp_kron_prod);
            tmp_kron_prod *= -sf(1,j);
            da1 += tmp_kron_prod;
            kron(a2, a_3, tmp_kron_prod);
            tmp_kron_prod *= -sf(6,j);
            da1 += tmp_kron_prod;
            kron(a21, a_3, tmp_kron_prod);
            tmp_kron_prod *= -sf(2,j);
            da1 += tmp_kron_prod;
            kron(a2, a_31, tmp_kron_prod);
            tmp_kron_prod *= -sf(2,j);
            da1 += tmp_kron_prod;


////cout << "EEFR 436" << endl;
            
            //da2 = -sf(5, j) * kron(a1, a_3)
            //- sf(1, j) * kron(a12, a_3)
            //- sf(1, j) * kron(a1, a_32)
            //- sf(4, j) * kron(a2, a_3)
            //- sf(2, j) * kron(a22, a_3)
            //- sf(2, j) * kron(a2, a_32);
            da2.set_all(0.0);// //cout << "EEFR 471" << endl;
            kron(a1, a_3, tmp_kron_prod);////cout << "EEFR 472" << endl;
            tmp_kron_prod *= -sf(5,j);////cout << "EEFR 473" << endl;
            da2 += tmp_kron_prod;////cout << "EEFR 474" << endl;
            kron(a12, a_3, tmp_kron_prod);////cout << "EEFR 475" << endl;
            tmp_kron_prod *= -sf(1,j);////cout << "EEFR 476" << endl;
            da2 += tmp_kron_prod;////cout << "EEFR 477" << endl;
            kron(a1, a_32, tmp_kron_prod);////cout << "EEFR 478" << endl;
            tmp_kron_prod *= -sf(1,j);////cout << "EEFR 479" << endl;
            da2 += tmp_kron_prod;////cout << "EEFR 480" << endl;
            kron(a2, a_3, tmp_kron_prod);////cout << "EEFR 481" << endl;
            tmp_kron_prod *= -sf(4,j);////cout << "EEFR 482" << endl;
            da2 += tmp_kron_prod;////cout << "EEFR 483" << endl;
            kron(a22, a_3, tmp_kron_prod);//cout << "EEFR 484" << endl;
            tmp_kron_prod *= -sf(2,j);//cout << "EEFR 485" << endl;
            da2 += tmp_kron_prod;////cout << "EEFR 486" << endl;
            kron(a2, a_32, tmp_kron_prod);////cout << "EEFR 487" << endl;
            tmp_kron_prod *= -sf(2,j);////cout << "EEFR 488" << endl;
            da2 += tmp_kron_prod;////cout << "EEFR 489" << endl;

////cout << "EEFR 438" << endl;
            // Matrix tempf_be = n1_be * sf(1, j) + m1_be * da1 + n2_be * sf(2, j) + m2_be * da2;
            tempf_be.set_all(0.0);
            colvec_matrix_multiplication(m1_be, da1, tmp_f); // (1, 3)(3, 3)->(1, 3)--> (3, 1),(3, 3) -> (3, 1)
            tempf_be += tmp_f;
            colvec_matrix_multiplication(m2_be, da2, tmp_l);
            tempf_be += tmp_l;
            const_multiplication(n1_be, sf(1, j), tmp_f);
            const_multiplication(n2_be, sf(2, j), tmp_l);
            tempf_be += tmp_f;
            tempf_be += tmp_l;
            tempf_be *= -sqa;
            f_be.set_row_from_col(j, tempf_be, 0); // the Force is the negative derivative
////cout << "EEFR 445" << endl;
            //tempf_cons = n1_cons * sf(1, j) + n2_cons * sf(2, j);
            const_multiplication(n1_cons, sf(1,j), tmp_f);
            const_multiplication(n2_cons, sf(2,j), tmp_l);
            addition(tmp_f, tmp_l, tempf_cons);
            tempf_cons *= -sqa;
            f_cons.set_row_from_col(j, tempf_cons, 0); // the Force is the negative derivative
////cout << "EEFR 448" << endl;
            //tempf_conv = n1_conv * sf(1, j)
            //+ n2_conv * sf(2, j)
            //+ tmp_evol * sf(0, j) * a_3;
            tempf_conv.set_all(0.0);
            const_multiplication(n1_conv, sf(1,j), tmp_f);
            tempf_conv += tmp_f;
            const_multiplication(n2_conv, sf(2,j), tmp_l);
            tempf_conv += tmp_l;
            const_multiplication(a_3, tmp_evol * sf(0,j), tmp_l);
            tempf_conv += tmp_l;
            tempf_conv *= -sqa;
            f_conv.set_row_from_col(j, tempf_conv, 0);
        }
////cout << "EEFR 451" << endl;
        meanCurv += halfGaussQuadratureCoeff * H_curv;
        eBend += halfGaussQuadratureCoeff * eBend_tmp;
        fBend += halfGaussQuadratureCoeff * f_be;
        fArea += halfGaussQuadratureCoeff * f_cons;
        fVolume += halfGaussQuadratureCoeff * f_conv;
        normVector += halfGaussQuadratureCoeff * a_3;
    }
    get_unit_vector(normVector);
}

/**
 * @brief Calculates the regularization energy and force for each face
 *
 * This function calculates the regularization energy and force
 * for each face in the mesh. The regularization energy is
 * calculated based on the deformation of the face's shape and
 * area relative to an equilateral triangle with the same side length.
 * The regularization force is then calculated based on the energy
 * and applied to the vertices of the face. The function also updates
 * the deformation counts for shapes and areas.
 *
 * @return void
 *
 */
void Mesh::energy_force_regularization()
{
    double kCurv = param.kCurv;

    // shape deformation meausures the shape degression of element triangles from an equilateral triangle
    int shapeDeformCount = 0;

    // area deformation measures the area difference from the target (equilibrium) element triangle area
    int areaDeformCount = 0;

    vector<Matrix> fReg(vertices.size(), mat_calloc(3, 1)); ///< Regularization force

//@todo Do I need reduction here?
// #pragma omp parallel for reduction(+ : fReg.data()[:vertices.size()])
#pragma omp parallel for
    for (Face &face : faces)
    {
        double eReg = 0.0; ///< Regularization energy

        // Indices of three vertices of this face (element triangle)
        int iVertex0 = face.adjacentVertices[0];
        int iVertex1 = face.adjacentVertices[1];
        int iVertex2 = face.adjacentVertices[2];

        // Vectors representing three sides of this face and magnitude
        Matrix vector10 = vertices[iVertex0].coord - vertices[iVertex1].coord;
        double length10 = vector10.calculate_norm();
        Matrix vector21 = vertices[iVertex1].coord - vertices[iVertex2].coord;
        double length21 = vector21.calculate_norm();
        Matrix vector02 = vertices[iVertex2].coord - vertices[iVertex0].coord;
        double length02 = vector02.calculate_norm();

        // semi-perimeter of the triangle
        double semiperi = (length10 + length21 + length02) / 2.0;

        // Area deformation
        //   s = 0.5(a + b + c)
        //   a = sqrt(s(s - a)(s - b)(s - c))
        double area = sqrt(semiperi * (semiperi - length10) * (semiperi - length21) * (semiperi - length02));

        // Shape deformation
        //   l = lmean = 1/3(a + b + c)
        //   gama = ((a - l)^2 + (b - l)^2 + (c - l)^2)/l^2
        // This calculate digression of triangle shape from equilateral
        // by calculating variance of three sides
        double meanSideLength = (length10 + length21 + length02) / 3.0;
        double gama = (pow(length10 - meanSideLength, 2.0) +
                       pow(length21 - meanSideLength, 2.0) +
                       pow(length02 - meanSideLength, 2.0)) /
                      pow(meanSideLength, 2.0);

        Matrix refVector10 = vertices[iVertex0].coordRef - vertices[iVertex1].coordRef;
        double refLength10 = refVector10.calculate_norm();
        Matrix refVector21 = vertices[iVertex1].coordRef - vertices[iVertex2].coordRef;
        double refLength21 = refVector21.calculate_norm();
        Matrix refVector02 = vertices[iVertex2].coordRef - vertices[iVertex0].coordRef;
        double refLength02 = refVector02.calculate_norm();
        // std::cout<<"v0"<<vector10.size()<<"::"<<refVector10.size()<<endl;

        // Here semipri = semi-perimeter of the reference triangle
        semiperi = (refLength10 + refLength21 + refLength02) / 2.0;

        // Calculate reference area
        double refArea = sqrt(semiperi * (semiperi - refLength10) * (semiperi - refLength21) * (semiperi - refLength02)); // double refArea = S0/face.n_rows;

        bool isDeformShape = (gama > param.gamaShape && param.usingRpi);
        bool isDeformArea = (abs(area - refArea) / refArea >= param.gamaArea && param.usingRpi);

        // convert vector10, vector21, vector02 unit vector
        vector10 /= length10;
        vector21 /= length21;
        vector02 /= length02;

        /*
         * Use a switch statement to check the values of
         * isDeformShape and isDeformArea , then performs the appropriate calculations based on each case.
         * Additionally, it may be possible to optimize the calculation of meanSideLength
         * by using a constant value or pre-calculated value instead of computing it at runtime.
         *
         * isDeformShape << isDeformArea:
         * -------------------------------
         * false, false: 0
         * false, true: 1
         * true, x: default
         */
        switch (isDeformShape << isDeformArea)
        {
        case 0:
            eReg = kCurv / 2.0 * (pow(length10 - refLength10, 2.0) + pow(length21 - refLength21, 2.0) + pow(length02 - refLength02, 2.0));
            fReg[iVertex0] += kCurv * ((refLength10 - length10) * (vector10) + (length02 - refLength02) * (vector02));
            fReg[iVertex1] += kCurv * ((refLength21 - length21) * (vector21) + (length10 - refLength10) * (vector10));
            fReg[iVertex2] += kCurv * ((refLength02 - length02) * (vector02) + (length21 - refLength21) * (vector21));
            break;
        case 1:
            {
                areaDeformCount++;
                double meanSideLengthRef = sqrt(4.0 * refArea / sqrt(3.0));
                eReg = kCurv / 2.0 * (pow(length10 - meanSideLengthRef, 2.0) + pow(length21 - meanSideLengthRef, 2.0) + pow(length02 - meanSideLengthRef, 2.0));
                fReg[iVertex0] += kCurv * ((meanSideLengthRef - length10) * (vector10) + (length02 - meanSideLengthRef) * (vector02));
                fReg[iVertex1] += kCurv * ((meanSideLengthRef - length21) * (vector21) + (length10 - meanSideLengthRef) * (vector10));
                fReg[iVertex2] += kCurv * ((meanSideLengthRef - length02) * (vector02) + (length21 - meanSideLengthRef) * (vector21));
            }
            break;
        default:
            {
                shapeDeformCount++;
                double meanSideLengthOld = sqrt(4.0 * area / sqrt(3.0));
                eReg = kCurv / 2.0 * (pow(length10 - meanSideLengthOld, 2.0) + pow(length21 - meanSideLengthOld, 2.0) + pow(length02 - meanSideLengthOld, 2.0));
                fReg[iVertex0] += kCurv * ((meanSideLengthOld - length10) * (vector10) + (length02 - meanSideLengthOld) * (vector02));
                fReg[iVertex1] += kCurv * ((meanSideLengthOld - length21) * (vector21) + (length10 - meanSideLengthOld) * (vector10));
                fReg[iVertex2] += kCurv * ((meanSideLengthOld - length02) * (vector02) + (length21 - meanSideLengthOld) * (vector21));
            }
            break;
        }

        // Store energy on faces
        face.energy.energyRegularization = eReg;

        vector10.free();
        vector21.free();
        vector02.free();
        refVector10.free();
        refVector21.free();
        refVector02.free();
    }

    // Store Force on vertices
#pragma omp parallel for
    for (int i = 0; i < vertices.size(); i++)
    {
        vertices[i].force.forceRegularization.free();
        vertices[i].force.forceRegularization = fReg[i];
    }

    // Update deformation count
    param.deformationCount.shapeDeformCount = shapeDeformCount;
    param.deformationCount.areaDeformCount = areaDeformCount;
    param.deformationCount.noDeformCount = faces.size() - shapeDeformCount - areaDeformCount;
}

/**
 * @brief Manages forces depending on different boundary conditions.
 *
 * This method sets the force of nodes that are part of the mesh's boundaries and ghost vertices to zero, based
 * on the type of boundary condition specified in the input parameters.
 *
 * @note See enum BoundaryType in Parameters.hpp
 *
 */
void Mesh::manage_force_for_boundary_ghost_vertex()
{

    switch (param.boundaryCondition)
    {
    // If the boundary condition is fixed, set the force of nodes on the boundary faces to zero.
    case BoundaryType::Fixed:
    {
        Force zeroForce;
#pragma omp parallel for
        for (Face &face : faces)
        {
            if (!face.isBoundary)
                continue;
            int node1 = face.adjacentVertices[0];
            int node2 = face.adjacentVertices[1];
            int node3 = face.adjacentVertices[2];
            vertices[node1].force = zeroForce;
            vertices[node2].force = zeroForce;
            vertices[node3].force = zeroForce;
        }
        break;
    }

    // If the boundary condition is periodic, set the force of ghost vertices to zero.
    case BoundaryType::Periodic:
    {
        Force zeroForce;
#pragma omp parallel for
        for (Vertex &vertex : vertices)
        {
            if (vertex.isGhost || vertex.isBoundary)
            {
                vertex.force = zeroForce;
            }
        }
        break;
    }

    // If the boundary condition is free, set the force of nodes on the outermost rows and columns of the
    // mesh to zero.
    case BoundaryType::Free:
    {
        Force zeroForce;
        int nFaceX = param.nFaceX; ///< x axis division
        int nFaceY = param.nFaceY; ///< y axis division
#pragma omp parallel for
        for (int i = 0; i < nFaceX + 1; i++)
        {
            int j = 0;
            int index = (nFaceX + 1) * j + i;
            vertices[index].force = zeroForce;
            index = (nFaceX + 1) * nFaceY + i;
            vertices[index].force = zeroForce;
        }
#pragma omp parallel for
        for (int j = 0; j < nFaceY + 1; j++)
        {
            int i = 0;
            int index = (nFaceX + 1) * j + i;
            vertices[index].force = zeroForce;
            index = (nFaceX + 1) * j + nFaceX;
            vertices[index].force = zeroForce;
        }
        break;
    }
    }
}

/**
 * @brief Get the max force scale of vertices.force.get_total_force_magnitude()
 *
 * @return double max force magnitude
 */
double Mesh::get_max_force_magnitude()
{
    double max_magnitude = 0;
    for (Vertex &vh : vertices)
    {
        double magnitude = vh.force.get_total_force_magnitude();
        if (magnitude > max_magnitude)
        {
            max_magnitude = magnitude;
        }
    }
    return max_magnitude;
}

double Mesh::calculate_mean_force()
{
    double sum = 0.0;
#pragma omp parallel for
    for (Vertex &vertex : vertices)
    {
        sum += vertex.force.get_total_force_magnitude();
    }
    return sum / vertices.size();
}
