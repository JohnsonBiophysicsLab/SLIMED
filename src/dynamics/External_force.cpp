/**
 * @file External_force.cpp
 * @brief The external force on a region of the membrane: which vertices it
 * acts on, when, and the file that records the region.
 *
 * The force itself is added in DynamicModel::next_step(), to the surface
 * coordinate the Brownian step displaces. See docs/external_force.md.
 */
#include "Dynamics.hpp"

#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>

void DynamicModel::select_external_force_vertices()
{
    const Param &param = mesh.param;
    externalForceMask.clear();
    externalForceVertexCount = 0;
    if (!param.externalForceEnabled)
    {
        return;
    }

    const bool disc = (param.externalForceShape == "disc");
    const bool rectangle = (param.externalForceShape == "rectangle");
    if (!disc && !rectangle)
    {
        throw std::runtime_error(
            "[DynamicModel::select_external_force_vertices] externalForceShape must be "
            "\"disc\" or \"rectangle\"; got \"" +
            param.externalForceShape + "\"");
    }
    if (disc && !(param.externalForceRadius > 0.0))
    {
        throw std::runtime_error(
            "[DynamicModel::select_external_force_vertices] externalForceRadius must be "
            "positive; got " +
            std::to_string(param.externalForceRadius));
    }
    if (rectangle && !(param.externalForceHalfWidthX > 0.0 && param.externalForceHalfWidthY > 0.0))
    {
        throw std::runtime_error(
            "[DynamicModel::select_external_force_vertices] externalForceHalfWidthX and "
            "externalForceHalfWidthY must be positive; got " +
            std::to_string(param.externalForceHalfWidthX) + " and " +
            std::to_string(param.externalForceHalfWidthY));
    }

    const int nVertices = static_cast<int>(mesh.vertices.size());
    externalForceMask.assign(nVertices, 0);
    for (int i = 0; i < nVertices; i++)
    {
        // Only a coordinate the step advances can be pushed. A fixed vertex,
        // a ghost, a periodic image or the global mode's slaved duplicate is
        // skipped by next_step(), so selecting it would only misreport the
        // count.
        if (!mesh.is_independent_vertex(i))
        {
            continue;
        }
        if (!mesh.isSlavedPeriodic.empty() && mesh.isSlavedPeriodic[i] != 0)
        {
            continue;
        }
        const double dx = mesh.vertices[i].coord(0, 0) - param.externalForceCenterX;
        const double dy = mesh.vertices[i].coord(1, 0) - param.externalForceCenterY;
        const bool inside =
            disc ? (dx * dx + dy * dy <= param.externalForceRadius * param.externalForceRadius)
                 : (std::fabs(dx) <= param.externalForceHalfWidthX &&
                    std::fabs(dy) <= param.externalForceHalfWidthY);
        if (inside)
        {
            externalForceMask[i] = 1;
            externalForceVertexCount++;
        }
    }

    if (externalForceVertexCount == 0)
    {
        throw std::runtime_error(
            "[DynamicModel::select_external_force_vertices] the external force region selects "
            "no vertex: a " +
            param.externalForceShape + " centred at (" + std::to_string(param.externalForceCenterX) +
            ", " + std::to_string(param.externalForceCenterY) +
            ") nm contains no free vertex. The generated sheet is centred on the origin; "
            "check externalForceCenterX/Y and the region size against sideX/sideY.");
    }

    std::cout << "[DynamicModel] External force (" << param.externalForceX << ", "
              << param.externalForceY << ", " << param.externalForceZ << ") pN per vertex on "
              << externalForceVertexCount << " vertices in a " << param.externalForceShape
              << " centred at (" << param.externalForceCenterX << ", " << param.externalForceCenterY
              << ") nm, from iteration " << param.externalForceStartIteration << " to ";
    if (param.externalForceStopIteration < 0)
    {
        std::cout << "the end of the run";
    }
    else
    {
        std::cout << param.externalForceStopIteration << " (exclusive)";
    }
    std::cout << "." << std::endl;
}

bool DynamicModel::external_force_active(int iteration) const
{
    const Param &param = mesh.param;
    if (!param.externalForceEnabled)
    {
        return false;
    }
    if (iteration < param.externalForceStartIteration)
    {
        return false;
    }
    if (param.externalForceStopIteration >= 0 && iteration >= param.externalForceStopIteration)
    {
        return false;
    }
    return true;
}

void DynamicModel::write_external_force_vertices_csv(const std::string &filepath) const
{
    std::ofstream out(filepath);
    if (!out)
    {
        throw std::runtime_error(
            "[DynamicModel::write_external_force_vertices_csv] cannot open " + filepath);
    }
    for (std::size_t i = 0; i < externalForceMask.size(); i++)
    {
        if (externalForceMask[i] != 0)
        {
            out << i << '\n';
        }
    }
}
