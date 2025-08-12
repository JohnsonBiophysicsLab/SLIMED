#include "mesh/Vertex.hpp"

/****************************************************************/
/**************************Constructors**************************/
/****************************************************************/

Vertex::Vertex() : coord(Matrix(3, 1))
{}

Vertex::Vertex(const int index, const double x, const double y, const double z) : index(index),
                                                          coord(Matrix(3, 1))
{
    coord.set(0, 0, x);
    coord.set(1, 0, y);
    coord.set(2, 0, z);
}

Vertex::Vertex(const int index, const Matrix& coord) : index(index),
                                               coord(coord)
{
}

/**
 * @brief Update the current vertex coordinate with the previous coordinate by
 * copying the values of the previous coordinate matrix to the current
 * coordinate. 
 * 
 */
void Vertex::update_coord_with_prev_coord()
{
    coord = Matrix(coordPrev);
}

/**
 * @brief Update the current vertex coordinate with the previous coordinate by
 * copying the values of the previous coordinate matrix to the current
 * coordinate.
 *
 */
void Vertex::approximate_unit_normal_vector(std::vector<Face> &faces)
{   
    
    Matrix totalNormVector = mat_calloc(3, 1);
    // iterate through all surrounding faces and add up normal vector
    //std::cout << "==========================================================" << std::endl;
    //std::cout << "==========================================================" << std::endl;
    //std::cout << totalNormVector << std::endl;


    for (int faceIndex : adjacentFaces){
        if (faces[faceIndex].normVector.get(2,0) > 0){
            std::cout<<"REVERSE NORM VEC ON " << faceIndex <<std::endl;
            negative(faces[faceIndex].normVector, faces[faceIndex].normVector);
        }
        //std::cout << "-------------------------------" << std::endl;
        //std::cout << "Face: " << faceIndex << std::endl;
        //std::cout << faces[faceIndex].normVector << std::endl;
        totalNormVector += faces[faceIndex].normVector;
        //std::cout << totalNormVector << std::endl;
        //std::cout << "-------------------------------" << std::endl;
    }

    // calcualte the unit vector and assign it to the member normal vector
    this->normVector = mat_calloc(3, 1); 
    //std::cout << "-------------------------------" << std::endl;
    //std::cout << totalNormVector << std::endl;
    if (totalNormVector.calculate_quad() > 1.0e-8){
        get_unit_vector(totalNormVector, this->normVector);
    }
    //std::cout << this->normVector << std::endl;
    //std::cout << "==========================================================" << std::endl;
    //std::cout << "==========================================================" << std::endl;

}

// Output operator for Matrix objects
std::ostream &operator<<(std::ostream &stream, const Vertex& vertex)
{
    stream << vertex.coord << std::endl;
    return stream;
}