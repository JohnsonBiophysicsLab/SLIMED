/**
 * @file Divergence_guard.cpp
 * @brief The run's divergence check (include/dynamics/Divergence_guard.hpp).
 */
#include "dynamics/Divergence_guard.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

DivergenceGuard::DivergenceGuard(const Mesh &mesh)
{
    std::array<double, 3> low;
    std::array<double, 3> high;
    low.fill(std::numeric_limits<double>::infinity());
    high.fill(-std::numeric_limits<double>::infinity());
    for (const Vertex &vertex : mesh.vertices)
    {
        for (int axis = 0; axis < 3; ++axis)
        {
            const double x = vertex.coord(axis, 0);
            low[axis] = std::min(low[axis], x);
            high[axis] = std::max(high[axis], x);
        }
    }
    // The size of the starting mesh: its largest side, and never less than an
    // edge, so a mesh that starts as a point or a line still gets a box.
    double size = std::max(mesh.param.lFace, 1.0);
    for (int axis = 0; axis < 3; ++axis)
    {
        if (low[axis] > high[axis])
        {
            low[axis] = high[axis] = 0.0;   // no vertices
        }
        size = std::max(size, high[axis] - low[axis]);
    }
    for (int axis = 0; axis < 3; ++axis)
    {
        lower_[axis] = low[axis] - kCoordinateReach * size;
        upper_[axis] = high[axis] + kCoordinateReach * size;
    }

    const double thermal = static_cast<double>(mesh.vertices.size()) * std::fabs(mesh.param.KBT);
    double scale = std::max(thermal, 1.0);
    if (std::isfinite(mesh.param.energy.energyTotal))
    {
        scale = std::max(scale, std::fabs(mesh.param.energy.energyTotal));
    }
    energyBound_ = kEnergyMultiple * scale;
}

std::string DivergenceGuard::check(const Mesh &mesh) const
{
    const double energy = mesh.param.energy.energyTotal;
    std::ostringstream reason;
    if (!std::isfinite(energy))
    {
        reason << "The total energy is " << energy;
        return reason.str();
    }

    static const char *const names[3] = {"x", "y", "z"};
    for (std::size_t i = 0; i < mesh.vertices.size(); ++i)
    {
        for (int axis = 0; axis < 3; ++axis)
        {
            const double x = mesh.vertices[i].coord(axis, 0);
            // Written so that a NaN fails it.
            if (!(x >= lower_[axis] && x <= upper_[axis]))
            {
                reason << "Vertex " << i << " has " << names[axis] << " = " << x
                       << " nm, outside [" << lower_[axis] << ", " << upper_[axis]
                       << "] (the starting mesh widened by " << kCoordinateReach
                       << " times its size)";
                return reason.str();
            }
        }
    }

    if (std::fabs(energy) > energyBound_)
    {
        reason << "The total energy is " << energy << " pN.nm, beyond " << energyBound_
               << " (" << kEnergyMultiple << " times its starting scale)";
        return reason.str();
    }
    return std::string();
}
