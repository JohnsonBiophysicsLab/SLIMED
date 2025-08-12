#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

// Function to format a number to 5 decimal places and 13 characters long.
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