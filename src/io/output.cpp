#include "io/io.hpp"

using namespace std;

/**
 * @brief Writes the vertex data for the mesh to a CSV file.
 *
 * This function writes the coordinates of each vertex in the mesh to a CSV file named "vertex[nFile].csv", where nFile is
 * equal to iteration/100. The output is formatted as follows:
 *
 *     x1,y1,z1
 *     x2,y2,z2
 *     ...
 *
 * where xi, yi, and zi are the coordinates of the ith vertex in the mesh.
 *
 * @note Step size and file name are left in definition in case customization is needed.
 *
 * @param mesh The mesh whose vertex data should be written to a CSV file.
 * @param iteration The current iteration number of the optimization algorithm.
 */
void write_vertex_data_to_csv(const Mesh &mesh, const int iteration)
{
    const int STEP_SIZE = 100;
    const std::string VERTEX_FILENAME_PREFIX = "vertex";
    const std::string VERTEX_FILENAME_SUFFIX = ".csv";
    const int nFile = iteration / STEP_SIZE;
    const std::string filename = VERTEX_FILENAME_PREFIX + std::to_string(nFile) + VERTEX_FILENAME_SUFFIX;

    std::ofstream outfile(filename);
    for (const auto &vertex : mesh.vertices)
    {
        outfile << vertex.coord(0, 0) << ',' << vertex.coord(1, 0) << ',' << vertex.coord(2, 0) << '\n';
    }
    outfile.close();
}

/**
 * @brief Writes the energy force data to a CSV file.
 *
 * This function writes the energy and force data for each iteration in the following format:
 *
 * "E_Curvature, E_Area, E_Regularization, E_Total ((pN.nm)), Mean Force (pN)"
 *
 * @note File name is left in definition in case customization is needed.
 *
 * @param model The model that cotains a mesh whose vertex data should be written to a CSV file.
 */
void write_energy_force_data_to_csv(const Model &model)
{
    // Output Energy and meanforce
    std::string ENERGY_FORCE_FILENAME = "EnergyForce.csv";
    ofstream outfileEF(ENERGY_FORCE_FILENAME);
    for (int j = 0; j <= model.iteration; j++)
    {
        if (j == 0)
        {
            outfileEF << "E_Curvature, E_Area, E_Regularization, E_Total ((pN.nm)), Mean Force (pN)" << '\n';
        }
        outfileEF << model.record.energyVec[j].energyCurvature << ", "
                  << model.record.energyVec[j].energyArea << ", "
                  << model.record.energyVec[j].energyRegularization << ", "
                  << model.record.energyVec[j].energyTotal << ", "
                  << model.record.meanForce[j] << '\n';
    }
    outfileEF.close();
}

/**
 * @brief Writes the current iteration element face energy data to a CSV file.
 *
 * This function writes the energy data for each element face of
 *  the current iteration in the following format:
 *
 * "E_Curvature, E_Area, E_Regularization, E_Total"
 *
 * @note @note File name is left in definition in case customization is needed.
 *
 * @param model The model that cotains a mesh whose vertex data should be written to a CSV file.
 *
 * 
 */
void write_element_face_energy_to_csv(const Model &model)
{
    // Output face index and then its corresponding element face energies
    std::string ELEMENT_FACE_ENERGY_FILENAME = "ElementFaceEnergy.csv";
    ofstream outfileEFE(ELEMENT_FACE_ENERGY_FILENAME);
    
    // Add header
    outfileEFE << "Face_index, E_Curvature, E_Area, E_Regularization, E_Total" << "\n";

    // iterate thru face and append to the output file
    for (Face& face: model.mesh.faces)
    {
        outfileEFE << face.index << ","
                << face.energy.energyCurvature << ","
                << face.energy.energyArea << ","
                << face.energy.energyTotal << "\n";
    }
    outfileEFE.close();
}

/**
 * @brief Writes the last vertex data for the mesh to a CSV file.
 *
 * This function writes the coordinates of each vertex in the last step
 * in the mesh to a CSV file named "vertexfinal.csv". The output is formatted as follows:
 *
 *     x1,y1,z1
 *     x2,y2,z2
 *     ...
 *
 * where xi, yi, and zi are the coordinates of the ith vertex in the mesh.
 *
 * @note file name are left in definition in case customization is needed.
 *
 * @param mesh The mesh whose vertex data should be written to a CSV file.
 */
void write_final_vertex_data_to_csv(const Mesh &mesh)
{
    const std::string FINAL_VERTEX_FILENANME = "vertexfinal.csv";

    ofstream outfile_final(FINAL_VERTEX_FILENANME);
    for (int j = 0; j < mesh.vertices.size(); j++)
    {
        outfile_final << setprecision(16) << mesh.vertices[j].coordPrev(0, 0) << ',' << mesh.vertices[j].coordPrev(1, 0) << ',' << mesh.vertices[j].coordPrev(2, 0) << '\n';
    }
    outfile_final.close();
}

void dynamics_create_trajectory_files(DynamicMesh &mesh, const std::string &filename)
{
        ofstream surfacepoint_csv("surfacepoint" + filename + ".csv");
        ofstream meshpoint_csv("meshpoint" + filename + ".csv");

        for (int i = 0; i < mesh.vertices.size(); i++)
        {
            meshpoint_csv << setprecision(8) << mesh.matMesh(i, 0) << ',' << mesh.matMesh(i, 1) << ',' << mesh.matMesh(i, 2) << ',';
            surfacepoint_csv << setprecision(8) << mesh.matSurface(i, 0) << ',' << mesh.matSurface(i, 1) << ',' << mesh.matSurface(i, 2) << ',';
        }
        surfacepoint_csv << '\n';
        meshpoint_csv << '\n';
        surfacepoint_csv.close();
        meshpoint_csv.close();
}

void output_trajectory_files(Mesh &mesh, const std::string &input_filename)
{
    ofstream meshpoint_csv_temp;
    meshpoint_csv_temp.open("meshpoint" + input_filename + ".csv", std::ios::app);

    for (int i = 0; i < mesh.vertices.size(); i++) {
        meshpoint_csv_temp << setprecision(8) << mesh.vertices[i].coord(0, 0) << ","
        << mesh.vertices[i].coord(1, 0) << ","
        << mesh.vertices[i].coord(2, 0) << ",";
    }
    meshpoint_csv_temp << '\n';
    meshpoint_csv_temp.close();
}

void dynamics_output_trajectory_files(DynamicMesh &mesh, const std::string &filename)
{
    if (mesh.param.meshpointOutput) {
        ofstream surfacepoint_csv_temp;
        ofstream meshpoint_csv_temp;
        surfacepoint_csv_temp.open("surfacepoint" + filename + ".csv", std::ios::app);
        meshpoint_csv_temp.open("meshpoint" + filename + ".csv", std::ios::app);
        
        for (int i = 0; i < mesh.vertices.size(); i++) {
            meshpoint_csv_temp << setprecision(8) << mesh.matMesh(i, 0) << ',' << mesh.matMesh(i, 1) << ',' << mesh.matMesh(i, 2) << ',';
            surfacepoint_csv_temp << setprecision(8) << mesh.matSurface(i, 0) << ',' << mesh.matSurface(i, 1) << ',' << mesh.matSurface(i, 2) << ',';
        }
        meshpoint_csv_temp << '\n';
        surfacepoint_csv_temp << '\n';
        meshpoint_csv_temp.close();
        surfacepoint_csv_temp.close();
    }
}

/*
std::string format_spaces_decimals(double number) {
    // Limit the values if they exceed the maximum/minimum possible value.
    if (number > 9999) {
        number = 9999.99999;
    }
    if (number < -9999) {
        number = -9999.99999;
    }

    // Format the number to 5 decimal places and 13 characters long.
    std::stringstream ss;
    ss.precision(5);
    ss << std::fixed << std::setw(13) << number;
    return ss.str();
}

// Function to convert a DataFrame to an XYZ string.
std::string to_xyz_string(const std::vector<std::vector<double>>& df) {
    std::stringstream xyz_string;

    // Loop over each row in the DataFrame.
    for (size_t row = 0; row < df.size(); row++) {
        // Calculate the number of vertices per iteration.
        size_t num_vertices = df[row].size() / 3;
        xyz_string << num_vertices << "\n";

        // Write the iteration number.
        xyz_string << "iteration: " << row << "\n";

        // Write the coordinates for each vertex.
        for (size_t i = 0; i < num_vertices; i++) {
            xyz_string << "  do";
            for (size_t j = 0; j < 3; j++) {
                double coord = df[row][i * 3 + j];
                xyz_string << format_spaces_decimals(coord);
            }
            xyz_string << "\n";
        }

        std::cout << "complete:" << row << " out of " << df.size() << std::endl;
    }

    return xyz_string.str();
}

// Main function.
int main() {
    // Read the input DataFrame from a CSV file.
    std::ifstream infile("meshpoint.csv");
    std::vector<std::vector<double>> df;
    std::string line;
    while (std::getline(infile, line)) {
        std::vector<double> row;
        std::istringstream iss(line);
        double val;
        while (iss >> val) {
            row.push_back(val);
            if (iss.peek() == ',')
                iss.ignore();
        }
        df.push_back(row);
    }

    // Convert the DataFrame to an XYZ string.
    std::string xyz_str = to_xyz_string(df);

    // Write the XYZ string to an output file.
    std::ofstream outfile("out_trajectory.xyz");
    outfile << xyz_str;
    outfile.close();

    return 0;
}
*/