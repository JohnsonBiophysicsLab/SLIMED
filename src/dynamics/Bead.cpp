/**
 * @file Bead.cpp
 * @brief A rigid spherical bead pushed through the membrane: where it is,
 * the harmonic wall that keeps the control points out of it, and the log.
 *
 * The force is added to the Brownian drive in DynamicModel::next_step().
 * See docs/bead_push_benchmark.md.
 */
#include "Dynamics.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>

void DynamicModel::setup_bead()
{
    const Param &param = mesh.param;
    beadForce = mat_calloc(mesh.vertices.size(), 3);
    beadCenter = {{param.beadCenterX, param.beadCenterY, param.beadCenterZ}};
    beadEnergy = 0.0;
    beadContacts = 0;
    beadMaxPenetration = 0.0;
    beadAdhesion = 0.0;
    beadAdhered = 0;
    beadReaction = {{0.0, 0.0, 0.0}};
    if (!param.beadEnabled)
    {
        return;
    }
    if (!(param.beadRadius > 0.0))
    {
        throw std::runtime_error("[DynamicModel::setup_bead] beadRadius must be positive; got " +
                                 std::to_string(param.beadRadius));
    }
    if (!(param.beadStiffness > 0.0))
    {
        throw std::runtime_error("[DynamicModel::setup_bead] beadStiffness must be positive; got " +
                                 std::to_string(param.beadStiffness));
    }
    std::cout << "[DynamicModel] Bead of radius " << param.beadRadius << " nm, wall "
              << param.beadStiffness << " pN/nm, from (" << param.beadCenterX << ", "
              << param.beadCenterY << ", " << param.beadCenterZ << ") nm at (" << param.beadVelocityX
              << ", " << param.beadVelocityY << ", " << param.beadVelocityZ
              << ") nm/us, moving from iteration " << param.beadStartIteration << " to ";
    if (param.beadStopIteration < 0)
    {
        std::cout << "the end of the run";
    }
    else
    {
        const std::array<double, 3> last = bead_center_at(param.beadStopIteration);
        std::cout << param.beadStopIteration << ", where it holds at (" << last[0] << ", " << last[1]
                  << ", " << last[2] << ") nm";
    }
    std::cout << "." << std::endl;
    if (param.beadAdhesionEnergy != 0.0)
    {
        if (!(param.beadAdhesionRange > 0.0))
        {
            throw std::runtime_error("[DynamicModel::setup_bead] beadAdhesionRange must be positive; got " +
                                     std::to_string(param.beadAdhesionRange));
        }
        std::cout << "[DynamicModel] Bead adhesion " << param.beadAdhesionEnergy << " pN.nm per vertex ("
                  << param.beadAdhesionEnergy / param.KBT << " kT) over " << param.beadAdhesionRange
                  << " nm outside the wall." << std::endl;
    }
}

std::array<double, 3> DynamicModel::bead_center_at(int iteration) const
{
    const Param &param = mesh.param;
    int moving = std::max(0, iteration - param.beadStartIteration);
    if (param.beadStopIteration >= 0)
    {
        moving = std::min(moving, std::max(0, param.beadStopIteration - param.beadStartIteration));
    }
    const double t = moving * param.timeStep;
    return {{param.beadCenterX + param.beadVelocityX * t, param.beadCenterY + param.beadVelocityY * t,
             param.beadCenterZ + param.beadVelocityZ * t}};
}

void DynamicModel::update_bead(int iteration)
{
    beadCenter = bead_center_at(iteration);
    refresh_bead_contact();
}

void DynamicModel::refresh_bead_contact()
{
    const Param &param = mesh.param;
    const int nVertices = static_cast<int>(mesh.vertices.size());
    const double radius = param.beadRadius;
    const double stiffness = param.beadStiffness;
    beadEnergy = 0.0;
    beadContacts = 0;
    beadMaxPenetration = 0.0;
    beadAdhesion = 0.0;
    beadAdhered = 0;
    beadReaction = {{0.0, 0.0, 0.0}};
    const double well = param.beadAdhesionEnergy;
    const double range = param.beadAdhesionRange;
    for (int i = 0; i < nVertices; i++)
    {
        beadForce.set(i, 0, 0.0);
        beadForce.set(i, 1, 0.0);
        beadForce.set(i, 2, 0.0);
        // Only a coordinate the step advances can be pushed; the rest are
        // held or slaved and a force on them would go nowhere.
        if (!mesh.is_independent_vertex(i))
        {
            continue;
        }
        if (!mesh.isSlavedPeriodic.empty() && mesh.isSlavedPeriodic[i] != 0)
        {
            continue;
        }
        const Matrix &coord = mesh.vertices[i].coord;
        const double dx = coord(0, 0) - beadCenter[0];
        const double dy = coord(1, 0) - beadCenter[1];
        const double dz = coord(2, 0) - beadCenter[2];
        const double d = std::sqrt(dx * dx + dy * dy + dz * dz);
        if (d >= radius)
        {
            // Outside the wall: the adhesion well, if there is one. A smooth
            // step from full depth at the surface to nothing at the range.
            if (well != 0.0 && range > 0.0 && d < radius + range)
            {
                const double x = (d - radius) / range;
                const double s = 1.0 - 3.0 * x * x + 2.0 * x * x * x;
                const double dsdx = -6.0 * x + 6.0 * x * x;
                // E = -w s(x); the radial force -dE/dd = w s'(x) / delta, which
                // is negative: toward the centre.
                const double f = well * dsdx / range;
                const double ux = dx / d;
                const double uy = dy / d;
                const double uz = dz / d;
                beadForce.set(i, 0, f * ux);
                beadForce.set(i, 1, f * uy);
                beadForce.set(i, 2, f * uz);
                beadAdhesion -= well * s;
                beadAdhered++;
                beadReaction[0] -= f * ux;
                beadReaction[1] -= f * uy;
                beadReaction[2] -= f * uz;
            }
            continue;
        }
        // Inside: the wall, and the full depth of the well.
        beadAdhesion -= well;
        beadAdhered++;
        const double depth = radius - d;
        // A vertex at the very centre has no outward direction; send it up.
        double ux = 0.0;
        double uy = 0.0;
        double uz = 1.0;
        if (d > 1e-12)
        {
            ux = dx / d;
            uy = dy / d;
            uz = dz / d;
        }
        const double f = stiffness * depth;
        beadForce.set(i, 0, f * ux);
        beadForce.set(i, 1, f * uy);
        beadForce.set(i, 2, f * uz);
        beadEnergy += 0.5 * stiffness * depth * depth;
        beadContacts++;
        beadMaxPenetration = std::max(beadMaxPenetration, depth);
        beadReaction[0] -= f * ux;
        beadReaction[1] -= f * uy;
        beadReaction[2] -= f * uz;
    }
}

void DynamicModel::append_bead_log(const std::string &filepath, int iteration)
{
    std::ofstream out(filepath, std::ios::app);
    if (!out)
    {
        throw std::runtime_error("[DynamicModel::append_bead_log] cannot open " + filepath);
    }
    if (!beadLogHeaderWritten)
    {
        out << "iteration,time_us,centerX,centerY,centerZ,energy_pNnm,contacts,maxPenetration_nm,"
               "forceOnBeadX,forceOnBeadY,forceOnBeadZ,adhesion_pNnm,adhered\n";
        beadLogHeaderWritten = true;
    }
    out << iteration << ',' << std::setprecision(10) << iteration * mesh.param.timeStep << ','
        << beadCenter[0] << ',' << beadCenter[1] << ',' << beadCenter[2] << ',' << beadEnergy << ','
        << beadContacts << ',' << beadMaxPenetration << ',' << beadReaction[0] << ','
        << beadReaction[1] << ',' << beadReaction[2] << ',' << beadAdhesion << ',' << beadAdhered << '\n';
}
