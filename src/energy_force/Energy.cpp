#include "energy_force/Energy.hpp"

Energy::Energy():
    energyCurvature(0.0),
    energyArea(0.0),
    energyVolume(0.0),
    energyThickness(0.0),
    energyTilt(0.0),
    energyRegularization(0.0),
    energyHarmonicBond(0.0),
    energyTotal(0.0)
{
}

Energy::Energy(const Energy& energy):
    energyCurvature (energy.energyCurvature),
    energyArea(energy.energyArea),
    energyVolume(energy.energyVolume),
    energyThickness(energy.energyThickness),
    energyTilt(energy.energyTilt),
    energyRegularization(energy.energyRegularization),
    energyHarmonicBond(energy.energyHarmonicBond),
    energyTotal(energy.energyTotal)
{
}

double Energy::calculateTotalEnergy() {
    energyTotal = energyCurvature + energyArea + energyVolume + energyThickness + energyTilt + energyRegularization + energyHarmonicBond;
    return energyTotal;
}

std::ostream& operator<<(std::ostream& stream, const Energy& energy) {
    stream << std::to_string(energy.energyTotal);
    return stream;
}

Energy& operator+=(Energy& e1, const Energy& e2) {
    e1.energyCurvature += e2.energyCurvature;
    e1.energyArea += e2.energyArea;
    e1.energyVolume += e2.energyVolume;
    e1.energyThickness += e2.energyThickness;
    e1.energyTilt += e2.energyTilt;
    e1.energyRegularization += e2.energyRegularization;
    e1.energyHarmonicBond += e2.energyHarmonicBond;
    e1.energyTotal += e2.energyTotal;

    return e1;
}