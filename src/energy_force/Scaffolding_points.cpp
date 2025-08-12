/**
 * @file Scaffolding_points.cpp
 * @author Y Ying (yying7@jh.edu)
 * @brief This file defines function in Mesh class that deals with
 * harmonic bond energy and force between the membrane and the
 * scaffolding.
 * @date 2023-04-07
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#include "mesh/Mesh.hpp"
#include "io/io.hpp"

using namespace std;

// --- 1. Find closest rcap from predefined list ---
int Mesh::findClosestRcap(double membraneCapRadius) {
    const std::vector<int> rcap_array = {10, 20, 25, 30, 35, 40};

    auto closest = std::min_element(rcap_array.begin(), rcap_array.end(),
        [membraneCapRadius](int a, int b) {
            return std::abs(a - membraneCapRadius) < std::abs(b - membraneCapRadius);
        });

    return *closest;
}

// --- 2. Interpolation ---
double Mesh::interpolateHeight(double r, const std::vector<double>& r_vals, const std::vector<double>& h_vals) {
    if (r <= r_vals.front()) return h_vals.front();
    if (r >= r_vals.back()) return h_vals.back();

    auto it = std::lower_bound(r_vals.begin(), r_vals.end(), r);
    size_t idx = std::distance(r_vals.begin(), it) - 1;

    double r0 = r_vals[idx], r1 = r_vals[idx + 1];
    double h0 = h_vals[idx], h1 = h_vals[idx + 1];

    return h0 + (r - r0) * (h1 - h0) / (r1 - r0);
}

// --- 3. Load and use in case3 ---
void Mesh::loadMembraneShapeProfile(double membraneCapRadius,
                              std::vector<double>& r_vals,
                              std::vector<double>& h_vals)
{
    int closest_rcap = findClosestRcap(membraneCapRadius);

    std::ostringstream filename;
    filename << "membrane_shape_" << closest_rcap << ".csv";

    std::vector<std::vector<double>> membrane_data = read_data_from_csv<double>(filename.str());

    r_vals.clear();
    h_vals.clear();
    for (const auto& row : membrane_data) {
        if (row.size() >= 2) {
            r_vals.push_back(row[0]);
            h_vals.push_back(row[1]);
        }
    }
}

/**
 * @brief
 * This method takes in the vector of spline point and calculate the
 * average coordinates. Based on the difference between spline points
 * and mesh vertices, a difference vector is calculated and compared to
 * the target bond length. (Supposed only in Z direction). Afterwards,
 * all the mesh points are moved in the direction of the target difference
 * vector.
 * 
 * @note
 * Refer to docs to see a graphic representation of the geometric variables
 *
 * @param fixDir default to true; fix move vector to (0, 0, 50) if set to
 * true; otherwise move vector is determined based on average scaffolding points
 * @return \code{true} if success
 */
bool Mesh::move_vertices_based_on_scaffolding(bool fixDir)
{
    // If fixed direction, all vertices will be moved upwards by bond length
    // (lbond defined in Param class)
    if (fixDir)
    {
        for (Vertex& vertex : vertices)
        { 
            vertex.coord.set(2, 0, vertex.coord(2, 0) + param.lbond);
        }
        return true;
    }
    else
    {   // If not fixed direction, the vertices will be moved based on the
        // input scaffolding
        std::cout << "[Mesh::move_vertices_based_on_scaffolding] Custom move vertices based on scaffolding."
            << std::endl;
        // Here we move the vertices so that the membrane cap wraps around
        // the gag cap perfectly. 
        bool useDefault = false; // For potential switch
        this->centerScaffoldingSphere = this->find_center_of_scaffolding_sphere(useDefault);

        double gagSphereRadius = param.scaffoldingSphereRaidus;  //R 
        double gagCapRadius = this->approximate_scaffolding_cap_radius(useDefault); std::cout << "50"<< std::endl;//r
        double lbond = param.lbond; //l
        double membraneSphereRadius = gagSphereRadius + lbond; // R+l
        double membraneCapRadiusSquared =  gagCapRadius * gagCapRadius / gagSphereRadius / gagSphereRadius
                                * membraneSphereRadius * membraneSphereRadius;//r^2/R^2*(R+l)^2
        double membraneCapRadius = sqrt(membraneCapRadiusSquared);
        // h0 = sqrt(Rs^2 - Rc^2)
        double h0 = sqrt(membraneSphereRadius * membraneSphereRadius - membraneCapRadiusSquared);
        std:: cout << "[Mesh::move_vertices_based_on_scaffolding] Membrane sphere radius approximated : "
            << membraneSphereRadius << std::endl;
        // iterate thru all vertices
        // if vertices is with cap region (x^2 + y^2 <= r^2 / R^2 * (R + l)^2)
        // then move the vertices up
        // z - z0 = sqrt((R + l)^2 - (x - x0)^2 - (y - y0)^2)
        // otherwise use default zShift defined as z - z0 according to the equation above
        // but instead xyQuadrance = r^2/R^2*(R+l)^2 = membraneCapRadiusSquared
        double zShiftRelaxStart = sqrt(membraneSphereRadius * membraneSphereRadius - membraneCapRadiusSquared);

        // Different modes for geometry model
        int mode = 2;
        switch (mode)
        {
            // Mode 0: move vertices on cap, ignore the transition region
            case 0:
                for (Vertex& vertex : vertices)
                {
                    double xShift = vertex.coord(0, 0) - this->centerScaffoldingSphere(0, 0); //x-x0
                    double yShift = vertex.coord(1, 0) - this->centerScaffoldingSphere(1, 0); //y-y0
                    double xyQuadrance = xShift * xShift + yShift * yShift; // (x - x0)^2 + (y - y0)^2
                    //examine if vertex is with in cap region
                    if (xyQuadrance <= membraneCapRadiusSquared)
                    {
                        // z - z0 = sqrt((R + l)^2 - (x - x0)^2 - (y - y0)^2)
                        double zShift = sqrt(membraneSphereRadius * membraneSphereRadius - xyQuadrance);
                        vertex.coord.set(2, 0, this->centerScaffoldingSphere(2, 0) + zShift);
                    } else {
                        vertex.coord.set(2, 0, this->centerScaffoldingSphere(2, 0) + zShiftRelaxStart);
                    }
                }
                break;
            // Mode 1: use a naive inversed arc to approximate transition region
            case 1:
                for (Vertex& vertex : vertices)
                {
                    // end of relaxation : z = 2h0 - Rs
                    double zShiftRelaxEnd = 2 * h0 - membraneSphereRadius;
                    double xShift = vertex.coord(0, 0) - this->centerScaffoldingSphere(0, 0); //x-x0
                    double yShift = vertex.coord(1, 0) - this->centerScaffoldingSphere(1, 0); //y-y0
                    double xyQuadrance = xShift * xShift + yShift * yShift; // (x - x0)^2 + (y - y0)^2
                    //examine if vertex is with in cap region
                    if (xyQuadrance <= membraneCapRadiusSquared)
                    {
                        // z - z0 = sqrt((R + l)^2 - (x - x0)^2 - (y - y0)^2)
                        double zShift = sqrt(membraneSphereRadius * membraneSphereRadius - xyQuadrance);
                        vertex.coord.set(2, 0, this->centerScaffoldingSphere(2, 0) + zShift);
                    } else if (xyQuadrance <= membraneCapRadiusSquared * 4) {
                        double rMinusTwoRc = sqrt(xyQuadrance) - 2 * membraneCapRadius;
                        // z = 2h0 - sqrt(- (r - 2Rc) ^ 2 + Rs^2)
                        double zShift = 2 * h0 - sqrt( - rMinusTwoRc * rMinusTwoRc
                        + membraneSphereRadius * membraneSphereRadius);
                        vertex.coord.set(2, 0, this->centerScaffoldingSphere(2, 0) + zShift);
                    } else {
                        vertex.coord.set(2, 0, this->centerScaffoldingSphere(2, 0) + zShiftRelaxEnd);
                    }
                }
                break;
            // Mode 2: use tensionless approximation with relaxation length = membrane cap radius
            // Eqn: h(r) = c1 * r^2 + c3 ln(r) + c4
            // c1 = r0/(2rR (2r0 + rR)) * tanTheta
            // c3 = - tanTheta * r0 (r0 + rR) ^ 2 / (rR(2r0 + rR))
            case 2:
                for (Vertex& vertex : vertices)
                {
                    // constants
                    double r0 = membraneCapRadius;
                    double rR = r0 + r0 * param.relaxLengthRatioApproximation;
                    double r0PlusrR = r0 + rR;
                    double twor0PlusrR = 2 * r0 + rR;
                    // tangent theta (sector angle)
                    double tangentSectorAngle = gagCapRadius / sqrt(gagSphereRadius * gagSphereRadius - gagCapRadius * gagCapRadius);
                    // c1 = r0/(2rR (2r0 + rR)) * tanTheta
                    double c1 = r0 / (2 * rR * twor0PlusrR) * tangentSectorAngle;
                    // c3 = - t * r0 (r0 + rR) ^ 2 / (rR(2r0 + rR))
                    double c3 = - r0 * r0PlusrR * r0PlusrR / rR / twor0PlusrR * tangentSectorAngle;
                    // c4 = -  c1 r0^2 - c3 ln(r0)
                    double c4 = - c1 * r0 * r0 - c3 * log(r0);

                    // Define geometry shift from origin
                    double xShift = vertex.coord(0, 0) - this->centerScaffoldingSphere(0, 0); //x-x0
                    double yShift = vertex.coord(1, 0) - this->centerScaffoldingSphere(1, 0); //y-y0
                    double xyQuadrance = xShift * xShift + yShift * yShift; // (x - x0)^2 + (y - y0)^2

                    // Use shape function to calcualte end of relaxation
                    double zShiftRelaxEnd = c1 * rR * rR + c3 * log(rR) + c4;

                    //examine if vertex is with in cap region
                    if (xyQuadrance <= membraneCapRadiusSquared)
                    {
                        // z - z0 = sqrt((R + l)^2 - (x - x0)^2 - (y - y0)^2)
                        double zShift = sqrt(membraneSphereRadius * membraneSphereRadius - xyQuadrance);
                        vertex.coord.set(2, 0, this->centerScaffoldingSphere(2, 0) + zShift);
                    } else if (xyQuadrance <= membraneCapRadiusSquared * 4) {
                        // Use shape function
                        
                        double r = sqrt(xyQuadrance);
                        double zShift = c1 * r * r + c3 * log(r) + c4;
                        // Rule out boundary point
                        if (vertex.isBoundary)
                        {
                            vertex.coord.set(2, 0, this->centerScaffoldingSphere(2, 0) + zShiftRelaxEnd + zShiftRelaxStart);
                        }
                        else
                        {
                            vertex.coord.set(2, 0, this->centerScaffoldingSphere(2, 0) + zShift + zShiftRelaxStart);
                        }
                    } else {
                        vertex.coord.set(2, 0, this->centerScaffoldingSphere(2, 0) + zShiftRelaxEnd + zShiftRelaxStart);
                    }
                }
                break;
            case 3:
                // Load membrane shape data
                std::vector<double> r_vals, h_vals;
                this->loadMembraneShapeProfile(membraneCapRadius, r_vals, h_vals);

                for (Vertex& vertex : vertices)
                {
                    // constants
                    double r0 = membraneCapRadius;
                    double rR = 120.0;
                    double r0PlusrR = r0 + rR;
                    double twor0PlusrR = 2 * r0 + rR;
                    // tangent theta (sector angle)
                    double tangentSectorAngle = gagCapRadius / sqrt(gagSphereRadius * gagSphereRadius - gagCapRadius * gagCapRadius);

                    // Define geometry shift from origin
                    double xShift = vertex.coord(0, 0) - this->centerScaffoldingSphere(0, 0); //x-x0
                    double yShift = vertex.coord(1, 0) - this->centerScaffoldingSphere(1, 0); //y-y0
                    double xyQuadrance = xShift * xShift + yShift * yShift; // (x - x0)^2 + (y - y0)^2

                    // Use shape function to calcualte end of relaxation
                    double zShiftRelaxEnd = this->interpolateHeight(rR, r_vals, h_vals);;

                    //examine if vertex is with in cap region
                    if (xyQuadrance <= membraneCapRadiusSquared)
                    {
                        // z - z0 = sqrt((R + l)^2 - (x - x0)^2 - (y - y0)^2)
                        double zShift = sqrt(membraneSphereRadius * membraneSphereRadius - xyQuadrance);
                        vertex.coord.set(2, 0, this->centerScaffoldingSphere(2, 0) + zShift);
                    } else if (xyQuadrance <= membraneCapRadiusSquared * 4) {
                        // Use shape function
                        
                        double r = sqrt(xyQuadrance);
                        double zShift = this->interpolateHeight(r, r_vals, h_vals);
                        // Rule out boundary point
                        if (vertex.isBoundary)
                        {
                            vertex.coord.set(2, 0, this->centerScaffoldingSphere(2, 0) + zShiftRelaxEnd + zShiftRelaxStart);
                        }
                        else
                        {
                            vertex.coord.set(2, 0, this->centerScaffoldingSphere(2, 0) + zShift + zShiftRelaxStart);
                        }
                    } else {
                        vertex.coord.set(2, 0, this->centerScaffoldingSphere(2, 0) + zShiftRelaxEnd + zShiftRelaxStart);
                    }
                }

                // Continue with the rest of your shape model...
                // e.g., use zShift to compute the membrane height at a point (x, y)
                break;
        }

        
        return true;
    }
}

/**
 * @brief A function that finds the center of the scaffolding sphere.
 * Currently the function assumes that the spherical cap is oriented
 * in a way that the point with biggest z-coord corresponds to the
 * x,y of the spherical center.
 * 
 * @note this function does NOT find the radius of the sphere. It takes
 * in the gag sphere radius from the parameter structure in mesh (this->param).
 * 
 * @param use_default if set to true, the function with return the defautl
 * center of scaffolding sphere (0, 0, 0)
 * @return instance of Matrix(3, 1) that represents the center of scaffolding
 * sphere.
 */
Matrix& Mesh::find_center_of_scaffolding_sphere(bool use_default){
    // if enabled, the default center of scaffolding sphere (0, 0, 0) will be used
    if (use_default){
        // return a matrix of (0, 0, 0)
        this->centerScaffoldingSphere = mat_calloc(3, 1);
        return this->centerScaffoldingSphere;
    }

    // Implement a simple function that determines center x-y based on max(z)
    // Move the center x-y by radius to determine the center coord
    double max_zcoord = vertices[0].coord(2, 0);
    int index_with_max_zcoord = 0;
    // Iterate thru all vertices and compare to find maximum
    for (int i = 0; i < param.scaffoldingPoints.size(); i++) {
        if (param.scaffoldingPoints[i](2, 0) > max_zcoord) {
            max_zcoord = param.scaffoldingPoints[i](2, 0);
            index_with_max_zcoord = i;
        }
    }
    // Generate a new copy of matrix representing the center of scaffolding sphere
    this->centerScaffoldingSphere = Matrix(param.scaffoldingPoints[index_with_max_zcoord]);
    this->centerScaffoldingSphere.set(2, 0, this->centerScaffoldingSphere(2, 0) - this->param.scaffoldingSphereRaidus);
    std::cout << "Center of scaffolding sphere" << this->centerScaffoldingSphere << std::endl;
    // Return the matrix
    return this->centerScaffoldingSphere;
}

/**
 * @brief Approximates the radius of the scaffolding cap.
 *
 * This function approximates the radius of the scaffolding cap based on the furthest point from the center of the scaffolding sphere in terms of x and y coordinates.
 * 
 * @param use_default Flag indicating whether to use the default center of scaffolding sphere.
 *                    If set to true, the default center (0,0,0) will be used.
 *                    If set to false, the center will be determined based on the maximum z-coordinate of the scaffolding points.
 * 
 * @return The radius of the scaffolding cap.
 */
double Mesh::approximate_scaffolding_cap_radius(bool use_default){
    // Get center of scaffolding sphere using the member function
    this->find_center_of_scaffolding_sphere(use_default);
    // Iterate thru all gags to find the gag that is furthest away in terms of x,y coord
    // Use quadrance (distance ^ 2)
    double max_quad_2d = 0.0;
    int index_with_max_2d_quad = 0;
    for (int i = 0; i < param.scaffoldingPoints.size(); i++) {
        // Calculate quadrance in x and y direction
        double dx = param.scaffoldingPoints[i](0, 0) - this->centerScaffoldingSphere(0, 0);
        double dy = param.scaffoldingPoints[i](1, 0) - this->centerScaffoldingSphere(1, 0);
        double quad_2d = dx * dx + dy * dy;
        // Get maximum quadrance
        if (quad_2d> max_quad_2d) {
            index_with_max_2d_quad = i;
            max_quad_2d = quad_2d;
        }
    }
    // Cap radius is the assumed to be distance in x,y direction
    double capRadius = sqrt(max_quad_2d);
    return capRadius;
}

/**
 * @brief
 * Used in getClosestVertexIndex
 * Calculates the squared distance between a point represented by a 3x1 matrix
 * and a Vertex object. Returns the squared distance.
 */
double Mesh::get_squared_distance_sp_and_v(const Matrix &scaffoldingPoint, const Vertex &vertex)
{
    // get difference vector
    Matrix diffVec = vertex.coord - scaffoldingPoint;
    return std::pow(diffVec(0, 0), 2) + std::pow(diffVec(1, 0), 2) + std::pow(diffVec(2, 0), 2);
}

/**
 * @brief
 * Gets a vector of indexes of vertices that are closest to
 * the `param.scaffoldingPoints` vector provided. Then sets
 * `param.scaffoldingPoints_correspondingVertexIndex` to represent
 * the vertices bonded with each scaffolding point.
 * 
 * @note
 * This section might be accelerated with ckdtree
 */
void Mesh::set_scaffolding_vertices_correspondence()
{
    // initialize
    vector<int> scaffoldingPoints_correspondingVertexIndex(param.scaffoldingPoints.size(), 0);
    // loop over spline points and find closest vertices for each
    for (int i = 0; i < param.scaffoldingPoints.size(); i++)
    {
        // initialize with the first vertex
        double minDistance2 = get_squared_distance_sp_and_v(param.scaffoldingPoints[i], vertices[0]);
        double newDistance2 = minDistance2;
        int minDistanceVertexIndex = 0;
        // loop over vertices to search for min distance
        for (int j = 1; j < vertices.size(); j++)
        {
            newDistance2 = get_squared_distance_sp_and_v(param.scaffoldingPoints[i], vertices[j]);
            // overide min and index if new d < min
            if (newDistance2 < minDistance2)
            {
                minDistance2 = newDistance2;
                minDistanceVertexIndex = j;
            }
        }
        // set corresponding vertex index to minDistanceVertexIndex
        scaffoldingPoints_correspondingVertexIndex[i] = minDistanceVertexIndex;
    }
    if (param.VERBOSE_MODE)
    {
        for (int i = 0; i < scaffoldingPoints_correspondingVertexIndex.size(); i++)
        {
            for (int j = 0; j < i; j++)
            {
                if (scaffoldingPoints_correspondingVertexIndex[i] == scaffoldingPoints_correspondingVertexIndex[j])
                    std::cout << "Repeated Vertex Index: " << scaffoldingPoints_correspondingVertexIndex[i] << endl;
            }
        }
    }
    param.scaffoldingPoints_correspondingVertexIndex = scaffoldingPoints_correspondingVertexIndex;
}

/*
 * Used in getCloesestVertexIndex
 * Calculate the distance between two points denoted by (3) vector
 * and Vertex respectively.
 * Return the distance
 */
double getDistance(Matrix &scaffoldingPoint, Vertex &vertex)
{
    // get difference vector
    bool singleSide = true; // set r to be negative when z_diffVec is negative
    // (membrane is lower than gag)
    Matrix diffVec = vertex.coord - scaffoldingPoint;
    // std::cout << "Mz: " << vertex.Coord[2] << "; Gz: " << scaffoldingPoint[2] << endl;
    double distance = diffVec.calculate_norm();
    if (singleSide && diffVec(2, 0) < 0.0)
    {
        return (-distance);
    }
    else
        return distance;
}

/**
 *
 * @brief Calculates the energy and force due to the harmonic bond between scaffold points and membrane vertices.
 * @param doLocalSearch Flag indicating whether or not to perform a local search for each vertex (default=false).
 * @return double The total energy of the system due to the scaffold-membrane interactions.
 *
 * The function iterates over each scaffold point, calculates the distance between the point and its corresponding vertex,
 * and calculates the energy and force due to the harmonic bond between the two. If the energy flag is not set to include
 * harmonic bonding, the function returns 0.
 *
 * @param[in] doLocalSearch Flag to indicate whether a local search should be performed for each vertex (default=false).
 *
 *
 * @note The energy is calculated according to the formula: E = 0.5 * k * (r - l)^2, where k is the spring constant, r is the
 * distance between the scaffold point and vertex, and l is the resting length of the bond.
 * @note The force is calculated according to the formula: F = -k * (r - l) * (unit vector pointing from vertex to scaffold point).
 * @note If the distance between the scaffold point and vertex is negative, the force is multiplied by -1 in order to prevent
 * overlapping between the scaffold and membrane.
 * @note If the verbose flag is set, the function prints out information about each vertex's force due to the scaffold-membrane bond.
 */
double Mesh::calculate_scaffolding_energy_force(bool doLocalSearch)
{
    if (!param.isEnergyHarmonicBondIncluded)
    {
        return 0;
    }

    // initialize total energy to zero; reset total force
    double totalEnergy = 0.0;
    forceTotalOnScaffolding = mat_calloc(3, 1);

    // iterate over spline points
    for (int i = 0; i < param.scaffoldingPoints.size(); i++)
    {
        int index = param.scaffoldingPoints_correspondingVertexIndex[i];
        double distance = getDistance(param.scaffoldingPoints[i], vertices[index]);
        // @test get absolute distance
        // distance = abs(distance);

        // calculate energy = 0.5 * k * (r - l)^2
        totalEnergy += 0.5 * param.springConst * pow(distance - param.lbond, 2);

        // vertices[index].energy.EnergySpline = 0.5 * springConst * (distance - lbond) * (distance - lbond);
        // vertices[index].energy.EnergyTotal += vertices[index].energy.EnergySpline;
        // calculate force = - k * (r - l) * normalize(vertexvec- splinepointvec)
        Matrix r_unit(3, 1);
        get_unit_vector(vertices[index].coord - param.scaffoldingPoints[i], r_unit);

        // If distance < 0.0, this means that the scaffolding points are above the plane
        // If using standard harmonic bonding equation, it would push the membrane further up
        // However because we do not want any overlapping between the membrane and the scaffolding
        // This is set to negative of original harmonic bond force in this case.
        if (distance > 0.0)
        {
            vertices[index].force.forceHarmonicBond = -param.springConst * (distance - param.lbond) * r_unit;
        }
        else
        {
            //adapt negative distance using:
            vertices[index].force.forceHarmonicBond = param.springConst * (distance - param.lbond) * r_unit;
            //vertices[index].force.forceHarmonicBond = param.springConst * (distance + param.lbond) * r_unit;
        }
        // Add harmonic bond force to total nodal force
        vertices[index].force.forceTotal += vertices[index].force.forceHarmonicBond;
        // Force exerted by membrane on scaffolding = reverse of force exerted by scaffolding on membrane
        forceTotalOnScaffolding -= vertices[index].force.forceHarmonicBond;
    }
    // Commandline print if enabled
    if (param.VERBOSE_MODE)
    {
    std::cout << "[Mesh::calculate_scaffolding_energy_force] Total force exerted on scaffolding = ("
    << forceTotalOnScaffolding(0, 0) << ", "
    << forceTotalOnScaffolding(1, 0) << ", "
    << forceTotalOnScaffolding(2, 0) << ") " << std::endl;
    }
    return totalEnergy;
}

bool Mesh::propagate_scaffolding()
{
    double scaffolding_propagation_stepsize = 1.0e-7;
    int max_scaffolding_propagation_nsteps = param.propagateScaffoldingNstep;
    for (int i = 0; i < max_scaffolding_propagation_nsteps; i++){
        // recalculate scaffolding energy and force
        this->calculate_scaffolding_energy_force(false);
        // move scaffolding based on calculated force
        Matrix displacement_vector = scaffolding_propagation_stepsize * forceTotalOnScaffolding;
        // iterate thru scaffolding and apply displacement vector
        for (int j = 0; j < param.scaffoldingPoints.size(); j++)
        {
            param.scaffoldingPoints[j] += displacement_vector;
        }
        // record movement in scaffoldingMovementVector
        scaffoldingMovementVector += displacement_vector;
    }
    std::cout << "[Mesh::propagate_scaffolding()] Total movement : ("
    << scaffoldingMovementVector(0, 0) << ", "
    << scaffoldingMovementVector(1, 0) << ", "
    << scaffoldingMovementVector(2, 0) << ") " << std::endl;
    return true;
}