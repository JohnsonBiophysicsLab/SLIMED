#include "mesh/Mesh.hpp"

namespace
{
Matrix current_interface_coord_for_output(const Mesh::GagSubunit &subunit,
                                          const std::string &interfaceName,
                                          const std::vector<Matrix> &coms)
{
    return coms[subunit.pointIndex] + subunit.rotation * subunit.localInterfaces.at(interfaceName);
}
}

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

void Mesh::write_gag_scaffolding_state_dat(const std::string &outfile_name) const
{
    if ((!param.isGagScaffoldingEnergyIncluded && !param.isIdealizedProteinLatticeEnergyIncluded) ||
        gagSubunits.empty())
    {
        return;
    }

    std::ofstream outfile(outfile_name);
    std::cout << "[Mesh::write_gag_scaffolding_state_dat] Creating output file : "
              << outfile_name << std::endl;

    std::vector<const GagSubunit *> orderedSubunits;
    orderedSubunits.reserve(gagSubunits.size());
    for (const GagSubunit &subunit : gagSubunits)
    {
        if (subunit.moleculeId >= 0 && subunit.pointIndex >= 0 &&
            subunit.pointIndex < static_cast<int>(param.scaffoldingPoints.size()))
        {
            orderedSubunits.push_back(&subunit);
        }
    }
    std::sort(orderedSubunits.begin(), orderedSubunits.end(),
              [](const GagSubunit *lhs, const GagSubunit *rhs) {
                  return lhs->moleculeId < rhs->moleculeId;
              });

    outfile << "{\n";
    outfile << "  \"molecules\": {\n";
    for (size_t subunitIndex = 0; subunitIndex < orderedSubunits.size(); ++subunitIndex)
    {
        const GagSubunit &subunit = *orderedSubunits[subunitIndex];
        const Matrix &com = param.scaffoldingPoints[subunit.pointIndex];
        outfile << "    \"" << subunit.moleculeId << "\": {\n";
        outfile << "      \"com\": [\n";
        outfile << "        " << std::setprecision(16) << com(0, 0) << ",\n";
        outfile << "        " << std::setprecision(16) << com(1, 0) << ",\n";
        outfile << "        " << std::setprecision(16) << com(2, 0) << "\n";
        outfile << "      ],\n";
        outfile << "      \"interfaces\": {\n";

        std::vector<std::string> interfaceNames;
        interfaceNames.reserve(subunit.localInterfaces.size());
        for (const auto &entry : subunit.localInterfaces)
        {
            interfaceNames.push_back(entry.first);
        }
        std::sort(interfaceNames.begin(), interfaceNames.end());

        for (size_t interfaceIndex = 0; interfaceIndex < interfaceNames.size(); ++interfaceIndex)
        {
            const std::string &interfaceName = interfaceNames[interfaceIndex];
            const Matrix interfaceCoord =
                current_interface_coord_for_output(subunit, interfaceName, param.scaffoldingPoints);
            outfile << "        \"" << interfaceName << "\": {\n";
            outfile << "          \"coord\": [\n";
            outfile << "            " << std::setprecision(16) << interfaceCoord(0, 0) << ",\n";
            outfile << "            " << std::setprecision(16) << interfaceCoord(1, 0) << ",\n";
            outfile << "            " << std::setprecision(16) << interfaceCoord(2, 0) << "\n";
            outfile << "          ]\n";
            outfile << "        }";
            outfile << (interfaceIndex + 1 < interfaceNames.size() ? ",\n" : "\n");
        }

        outfile << "      }\n";
        outfile << "    }";
        outfile << (subunitIndex + 1 < orderedSubunits.size() ? ",\n" : "\n");
    }
    outfile << "  },\n";
    outfile << "  \"bound_states\": [\n";

    for (size_t interactionIndex = 0; interactionIndex < gagInteractions.size(); ++interactionIndex)
    {
        const GagInteraction &interaction = gagInteractions[interactionIndex];
        const GagSubunit &subunit1 = gagSubunits[interaction.subunit1];
        const GagSubunit &subunit2 = gagSubunits[interaction.subunit2];

        outfile << "    {\n";
        outfile << "      \"mol1\": {\n";
        outfile << "        \"molecule_index\": " << subunit1.moleculeId << ",\n";
        outfile << "        \"interface_name\": \"" << interaction.interface1 << "\",\n";
        outfile << "        \"interface_index\": 0\n";
        outfile << "      },\n";
        outfile << "      \"mol2\": {\n";
        outfile << "        \"molecule_index\": " << subunit2.moleculeId << ",\n";
        outfile << "        \"interface_name\": \"" << interaction.interface2 << "\",\n";
        outfile << "        \"interface_index\": 0\n";
        outfile << "      }\n";
        outfile << "    }";
        outfile << (interactionIndex + 1 < gagInteractions.size() ? ",\n" : "\n");
    }

    outfile << "  ]\n";
    outfile << "}\n";
    outfile.close();
}
