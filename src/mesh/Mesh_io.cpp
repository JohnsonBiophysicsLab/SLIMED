#include "mesh/Mesh.hpp"

void Mesh::write_faces_csv(const std::string &outfile_name)
{
    // Generate ofstream based on given string
    std::ofstream outfile(outfile_name);
    std::cout << "[Mesh::write_faces_csv] Creating output file : " << outfile_name << std::endl;
    // Output csv file
    for (Face& face : faces)
    {
        outfile << face.adjacentVertices[0] << ","
                << face.adjacentVertices[1] << ","
                << face.adjacentVertices[2] << "\n";
    }
    outfile.close();
}

/**
 * @brief Writes a csv file containing the coordinates of each vertex in the mesh
 *
 * @param outfile_name The name of the output csv file
 */
void Mesh::write_vertices_csv(const std::string &outfile_name)
{
    // Generate ofstream based on given string
    std::ofstream outfile(outfile_name);
    std::cout << "[Mesh::write_vertices_csv] Creating output file : " << outfile_name << std::endl;
    // Output csv file
    for (Vertex& vertex : vertices)
    {
        outfile << setprecision(16)
                << vertex.coord(0, 0) << ", "
                << vertex.coord(1, 0) << ", "
                << vertex.coord(2, 0) << std::endl;
        
    }
    outfile.close();
}

/**
 * @brief Writes a csv file containing the coordinates and types of each vertex in the mesh
 *
 * @param outfile_name The name of the output csv file
 */
void Mesh::write_vertices_csv_with_type(const std::string &outfile_name)
{
    // Generate ofstream based on given string
    std::ofstream outfile(outfile_name);
    std::cout << "[Mesh::write_vertices_csv] Creating output file : " << outfile_name << std::endl;
    // Output csv file
    for (Vertex& vertex : vertices)
    {
        outfile << setprecision(16)
                << vertex.coord(0, 0) << ", "
                << vertex.coord(1, 0) << ", "
                << vertex.coord(2, 0) << ", "
                << static_cast<int>(vertex.type) << ", "
                << vertex.reflectiveVertexIndex << ", " << std::endl;
        
    }
    outfile.close();
}