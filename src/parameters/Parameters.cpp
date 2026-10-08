#include "Parameters.hpp"
#include <cctype>

std::ostream& operator<<(std::ostream& os, const Param& param)
{
    os << "VERBOSE_MODE: " << param.VERBOSE_MODE << std::endl;
    os << "kCurv: " << param.kCurv << std::endl;
    os << "uSurf: " << param.uSurf << std::endl;
    os << "uVol: " << param.uVol << std::endl;
    os << "kReg: " << param.kReg << std::endl;
    os << "kSpring: " << param.kSpring << std::endl;
    os << "area0: " << param.area0 << std::endl;
    os << "area: " << param.area << std::endl;
    os << "vol0: " << param.vol0 << std::endl;
    os << "vol: " << param.vol << std::endl;
    os << "insertCurv: " << param.insertCurv << std::endl;
    os << "spontCurv: " << param.spontCurv << std::endl;
    os << "sideX: " << param.sideX << std::endl;
    os << "sideY: " << param.sideY << std::endl;
    os << "radius: " << param.radius << std::endl;
    os << "lFace: " << param.lFace << std::endl;
    os << "nFaceX: " << param.nFaceX << std::endl;
    os << "nFaceY: " << param.nFaceY << std::endl;
    os << "dFaceX: " << param.dFaceX << std::endl;
    os << "dFaceY: " << param.dFaceY << std::endl;
    os << "meanL: " << param.meanL << std::endl;
    os << "sigma: " << param.sigma << std::endl;
    os << "isInsertionAreaConstraint: " << param.isInsertionAreaConstraint << std::endl;
    os << "isAdditiveScheme: " << param.isAdditiveScheme << std::endl;
    os << "isGlobalConstraint: " << param.isGlobalConstraint << std::endl;
    os << "elementTriangleArea0: " << param.elementTriangleArea0 << std::endl;
    os << "gaussQuadratureN: " << param.gaussQuadratureN << std::endl;
    os << "VWU: " << std::endl << param.VWU << std::endl;
    os << "gaussQuadratureCoeff: " << std::endl << param.gaussQuadratureCoeff << std::endl;
    os << "shapeFunctions: " << std::endl;
    for (const auto& sf : param.shapeFunctions)
    {
        os << sf << std::endl;
    }
    os << "boundaryCondition: " << (int)(param.boundaryCondition) << std::endl;
    os << "boundaryConditionX: " << boundary_type_name(param.boundaryConditionX) << std::endl;
    os << "boundaryConditionY: " << boundary_type_name(param.boundaryConditionY) << std::endl;
    os << "fixedBoundaryRings: " << param.fixedBoundaryRings << std::endl;
    os << "meshVerticesFile: " << param.meshVerticesFile << std::endl;
    os << "meshFacesFile: " << param.meshFacesFile << std::endl;
    os << "usingNCG: " << param.usingNCG << std::endl;
    os << "isNCGstuck: " << param.isNCGstuck << std::endl;
    os << "gamaShape: " << param.gamaShape << std::endl;
    os << "gamaArea: " << param.gamaArea << std::endl;
    os << "usingRpi: " << param.usingRpi << std::endl;
    os << "xyzOutput: " << param.xyzOutput << std::endl;
    os << "meshpointOutput: " << param.meshpointOutput << std::endl;
    os << "isEnergyHarmonicBondIncluded: " << param.isEnergyHarmonicBondIncluded << std::endl;
    os << std::endl << "springConst: " << param.springConst << std::endl;
    os << "springConstScaler: " << param.springConstScaler << std::endl;
    os << "springConstScalingInterv: " << param.springConstScalingInterv << std::endl;
    os << "springConstUpperBound: " << param.springConstUpperBound << std::endl;
    os << "splinePointsZcoordScaling: " << param.splinePointsZcoordScaling << std::endl;
    os << "lbond: " << param.lbond << std::endl;
    os << "scaffoldingZeroPlaneZ: " << param.scaffoldingZeroPlaneZ << std::endl;
    os << "scaffoldingFileName: " << param.scaffoldingFileName << std::endl;
    os << "isGagScaffoldingEnergyIncluded: " << param.isGagScaffoldingEnergyIncluded << std::endl;
    os << "gagReferenceStateFileName: " << param.gagReferenceStateFileName << std::endl;
    os << "gagReactionFileName: " << param.gagReactionFileName << std::endl;
    os << "isIdealizedProteinLatticeEnergyIncluded: " << param.isIdealizedProteinLatticeEnergyIncluded << std::endl;
    os << "idealizedProteinLatticeFileName: " << param.idealizedProteinLatticeFileName << std::endl;
    os << "gagKsigma: " << param.gagKsigma << std::endl;
    os << "gagKtheta: " << param.gagKtheta << std::endl;
    os << "gagKphi: " << param.gagKphi << std::endl;
    os << "gagKomega: " << param.gagKomega << std::endl;
    os << "gagFiniteDifferenceStep: " << param.gagFiniteDifferenceStep << std::endl;
    os << "gagPropagationStepSize: " << param.gagPropagationStepSize << std::endl;
    os << "gagPreRelaxSteps: " << param.gagPreRelaxSteps << std::endl;
    os << "numIterations: " << param.maxIterations << std::endl;
    os << "timeStep: " << param.timeStep << std::endl;
    os << "diffConst: " << param.diffConst << std::endl;
    os << "KBT: " << param.KBT << std::endl;
    os << "thermalFluctuationEnabled: " << param.thermalFluctuationEnabled << std::endl;
    os << "thermalFluctuationPureMMC: " << param.thermalFluctuationPureMMC << std::endl;
    os << "thermalFluctuationInterval: " << param.thermalFluctuationInterval << std::endl;
    os << "thermalFluctuationTemperatureKelvin: " << param.thermalFluctuationTemperatureKelvin << std::endl;
    os << "thermalFluctuationMinTemperatureKelvin: " << param.thermalFluctuationMinTemperatureKelvin << std::endl;
    os << "thermalFluctuationCoolingRate: " << param.thermalFluctuationCoolingRate << std::endl;
    os << "thermalFluctuationStepScale: " << param.thermalFluctuationStepScale << std::endl;
    os << "inPlaneDynamicsEnabled: " << param.inPlaneDynamicsEnabled << std::endl;
    os << "edgeSpringEnabled: " << param.edgeSpringEnabled << std::endl;
    os << "edgeSpringConstant: " << param.edgeSpringConstant << std::endl;
    os << "edgeSpringRestLength: " << param.edgeSpringRestLength << std::endl;
    os << "edgeTetherShape: " << param.edgeTetherShape << std::endl;
    os << "edgeTetherMinRatio: " << param.edgeTetherMinRatio << std::endl;
    os << "edgeTetherMaxRatio: " << param.edgeTetherMaxRatio << std::endl;
    os << "triangleShapeEnabled: " << param.triangleShapeEnabled << std::endl;
    os << "triangleShapeMinAltitudeRatio: " << param.triangleShapeMinAltitudeRatio << std::endl;
    os << "triangleShapeConstant: " << param.triangleShapeConstant << std::endl;
    os << "creaseWallEnabled: " << param.creaseWallEnabled << std::endl;
    os << "creaseWallAngle: " << param.creaseWallAngle << std::endl;
    os << "creaseWallConstant: " << param.creaseWallConstant << std::endl;
    os << "surfaceSolver: " << param.surfaceSolver << std::endl;
    os << "irregularPatchDepthScale: " << param.irregularPatchDepthScale << std::endl;
    os << "edgeFlipEnabled: " << param.edgeFlipEnabled << std::endl;
    os << "edgeFlipAttemptRate: " << param.edgeFlipAttemptRate << std::endl;
    os << "edgeFlipInterval: " << param.edgeFlipInterval << std::endl;
    os << "edgeFlipMinValence: " << param.edgeFlipMinValence << std::endl;
    os << "edgeFlipMaxValence: " << param.edgeFlipMaxValence << std::endl;
    os << "meshWallType: " << param.meshWallType << std::endl;
    os << "vertexMoveEnabled: " << param.vertexMoveEnabled << std::endl;
    os << "vertexMoveStepRatio: " << param.vertexMoveStepRatio << std::endl;
    os << "vertexMoveAttemptsPerVertex: " << param.vertexMoveAttemptsPerVertex << std::endl;
    os << "vertexMoveInterval: " << param.vertexMoveInterval << std::endl;
    os << "brownianStepEnabled: " << param.brownianStepEnabled << std::endl;
    os << "externalForceEnabled: " << param.externalForceEnabled << std::endl;
    os << "externalForceX: " << param.externalForceX << std::endl;
    os << "externalForceY: " << param.externalForceY << std::endl;
    os << "externalForceZ: " << param.externalForceZ << std::endl;
    os << "externalForceShape: " << param.externalForceShape << std::endl;
    os << "externalForceCenterX: " << param.externalForceCenterX << std::endl;
    os << "externalForceCenterY: " << param.externalForceCenterY << std::endl;
    os << "externalForceRadius: " << param.externalForceRadius << std::endl;
    os << "externalForceHalfWidthX: " << param.externalForceHalfWidthX << std::endl;
    os << "externalForceHalfWidthY: " << param.externalForceHalfWidthY << std::endl;
    os << "externalForceStartIteration: " << param.externalForceStartIteration << std::endl;
    os << "externalForceStopIteration: " << param.externalForceStopIteration << std::endl;
    os << "beadEnabled: " << param.beadEnabled << std::endl;
    os << "beadRadius: " << param.beadRadius << std::endl;
    os << "beadStiffness: " << param.beadStiffness << std::endl;
    os << "beadCenterX: " << param.beadCenterX << std::endl;
    os << "beadCenterY: " << param.beadCenterY << std::endl;
    os << "beadCenterZ: " << param.beadCenterZ << std::endl;
    os << "beadVelocityX: " << param.beadVelocityX << std::endl;
    os << "beadVelocityY: " << param.beadVelocityY << std::endl;
    os << "beadVelocityZ: " << param.beadVelocityZ << std::endl;
    os << "beadStartIteration: " << param.beadStartIteration << std::endl;
    os << "beadStopIteration: " << param.beadStopIteration << std::endl;
    os << "beadAdhesionEnergy: " << param.beadAdhesionEnergy << std::endl;
    os << "beadAdhesionRange: " << param.beadAdhesionRange << std::endl;
    os << "elasticBoxEnabled: " << param.elasticBoxEnabled << std::endl;
    os << "elasticBoxTension: " << param.elasticBoxTension << std::endl;
    os << "elasticBoxMobility: " << param.elasticBoxMobility << std::endl;
    os << "elasticBoxThermal: " << param.elasticBoxThermal << std::endl;
    os << "elasticBoxCoupling: " << param.elasticBoxCoupling << std::endl;
    os << "elasticBoxMaxStrainPerStep: " << param.elasticBoxMaxStrainPerStep << std::endl;
    os << "elasticBoxVertexGasTerm: " << param.elasticBoxVertexGasTerm << std::endl;
    os << "lipidFlowMode: " << param.lipidFlowMode << std::endl;
    os << "lipidFlowDiffusivity: " << param.lipidFlowDiffusivity << std::endl;
    os << "lipidFlowSubsteps: " << param.lipidFlowSubsteps << std::endl;
    os << "lipidFlowThermal: " << param.lipidFlowThermal << std::endl;
    os << "lipidFlowLinearStretching: " << param.lipidFlowLinearStretching << std::endl;
    os << "lipidReservoirRegion: " << param.lipidReservoirRegion << std::endl;
    os << "lipidReservoirTension: " << param.lipidReservoirTension << std::endl;
    os << "lipidReservoirRate: " << param.lipidReservoirRate << std::endl;
    os << "lipidFlipResampling: " << param.lipidFlipResampling << std::endl;
    os << "localOutput: " << param.localOutput << std::endl;
    os << "energy: " << param.energy << std::endl;
    os << "energyPrev: " << param.energyPrev << std::endl;
    return os;
}

const char *boundary_type_name(BoundaryType type)
{
    switch (type)
    {
    case BoundaryType::Fixed:
        return "Fixed";
    case BoundaryType::Periodic:
        return "Periodic";
    case BoundaryType::Free:
        return "Free";
    case BoundaryType::Mixed:
        return "Mixed";
    }
    return "unknown";
}

bool parse_boundary_type(const std::string &text, BoundaryType &type)
{
    std::string lowered;
    lowered.reserve(text.size());
    for (char c : text)
    {
        if (c == ' ' || c == '\t' || c == '\r')
        {
            continue;
        }
        lowered.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
    }
    if (lowered == "fixed")
    {
        type = BoundaryType::Fixed;
        return true;
    }
    if (lowered == "periodic")
    {
        type = BoundaryType::Periodic;
        return true;
    }
    if (lowered == "free")
    {
        type = BoundaryType::Free;
        return true;
    }
    if (lowered == "mixed" || lowered == "pervertex")
    {
        type = BoundaryType::Mixed;
        return true;
    }
    return false;
}
