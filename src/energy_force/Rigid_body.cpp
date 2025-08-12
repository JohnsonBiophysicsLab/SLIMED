/*#include "energy_force.hpp"

//using namespace std;

void calculate_rigid_force_and_torque(
        const Param& param,
        const vector<vector<double>>& scaffoldingPoints,
        const vector<Vertex>& vertices,
        const vector<int> scaffoldingPoints_correspondingVertexIndex,
        vector<vector<double>>& rigidBodyForces,



    //1. Initialize center of mass (COM) of scaffolding lattice
    vector<double> splinePointsCOM(3, 0.0);
    //2. Get a vector of forces based on vertices and correspondence bwt
    //spline points and vertices
    vector<vector<double>> rigidBodyRawForces(scaffoldingPoints.size(), vector<double>(3));
    for (int i = 0; i < scaffoldingPoints.size(); i++){
        //copy spline force to rigid body raw forces
        rigidBodyRawForces[i] = vertices[scaffoldingPoints_correspondingVertexIndex[i]].force.ForceSpline;
        //Calculate center of mass (COM) of scaffolding lattice - sum up coord first
        splinePointsCOM += scaffoldingPoints[i];
    }
    //Calculate COM - divide by number of spline points
    splinePointsCOM = splinePointsCOM / scaffoldingPoints.size();

    //3. Separate raw force into force in the direction of COM->scaffolding pt (for force) and
    //force perpendicular to COM->scaffolding pt (for torque):
    //f_tot = f_r + f_h
    vector<double> vecCOMScf(3, 0.0);
    for (int i = 0; i < scaffoldingPoints.size(); i++){
        //3.1 Calculate vec COM->scaffodling pt
        vecCOMScf = scaffoldingPoints[i] - splinePointsCOM;

        //3.2 Calculate force in the direction of COM->scaffolding(xvec):
        //f_r = (f_tot dot xvec) / (xvec dot xvec) * xvec
        rigidBodyForces[i] = ((dot(rigidBodyRawForces[i], vecCOMScf))/(dot(vecCOMScf, vecCOMScf))) * vecCOMScf;
        rigidBodyMeanForce = rigidBodyMeanForce + rigidBodyForces[i];
        //f_h = f_tot - f_r
        //torque = xvec x f_h
        rigidBodyTorques[i] = cross(vecCOMScf, (rigidBodyRawForces[i] - rigidBodyForces[i]));
        rigidBodyMeanTorque = rigidBodyMeanTorque + rigidBodyTorques[i];

        //3.3 Subtract COM from all spline points for rotation
        //splinePoints (COM zeroed) Value confirmed with test below
        splinePoints[i] += - splinePointsCOM;
        //@TEST
        //cout << "-COM:" << splinePoints[i][0] << "," << splinePoints[i][1] << "," << splinePoints[i][2] << endl;
    }
    
    //4. calculate mean force and torque and apply scale
    rigidBodyMeanForce = rigidBodyMeanForce / splinePoints.size() * forceDispScaleConst;
    rigidBodyMeanTorque = rigidBodyMeanTorque / splinePoints.size() * torqueDispScaleConst;
    cout << "MF: " << rigidBodyMeanForce[0] << ", " << rigidBodyMeanForce[1] << ", " << rigidBodyMeanForce[2] << endl;
    cout << "MT: " << rigidBodyMeanTorque[0] << ", " << rigidBodyMeanTorque[1] << ", " << rigidBodyMeanTorque[2] << endl;

    //5. move spline points based on forces and torques
    //Values confirmed with TEST below and rotation matrix math confirmed with
    //Wolfram alpha and ChatGPT
    gsl_matrix* rotMatrix = get_rotation_matrix(rigidBodyMeanTorque);
    //@TEST
    //print_mat_contents(rotMatrix);
    for (int i = 0; i < splinePoints.size(); i++){
        gsl_matrix* coordMatrix = vec_to_col_matrix(splinePoints[i]);
        
        // rotate
        gsl_matrix* coordMatrixRotated = gsl_matrix_calloc(3,1);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, rotMatrix, coordMatrix, 0.0, coordMatrixRotated);
        //@TEST
        //print_mat_contents(coordMatrixRotated);
        
        // convert back
        splinePoints[i] = col_matrix_to_vec(coordMatrixRotated);
        
        // add COM and force
        splinePoints[i] += splinePointsCOM + rigidBodyMeanForce;
        //@TEST
        //cout << "-COM:" << splinePoints[i][0] << "," << splinePoints[i][1] << "," << splinePoints[i][2] << endl;
        
    }
    
}
*/