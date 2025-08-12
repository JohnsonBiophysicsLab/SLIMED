#include "Parameters.hpp"

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
    os << "subDivideTimes: " << param.subDivideTimes << std::endl;
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
    os << "subMatrix.irregM: " << std::endl << param.subMatrix.irregM << std::endl;
    os << "subMatrix.irregM1: " << std::endl << param.subMatrix.irregM1 << std::endl;
    os << "subMatrix.irregM2: " << std::endl << param.subMatrix.irregM2 << std::endl;
    os << "subMatrix.irregM3: " << std::endl << param.subMatrix.irregM3 << std::endl;
    os << "subMatrix.irregM4: " << std::endl << param.subMatrix.irregM4 << std::endl;
    os << "boundaryCondition: " << (int)(param.boundaryCondition) << std::endl;
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
    os << "scaffoldingFileName: " << param.scaffoldingFileName << std::endl;
    os << "numIterations: " << param.maxIterations << std::endl;
    os << "timeStep: " << param.timeStep << std::endl;
    os << "diffConst: " << param.diffConst << std::endl;
    os << "KBT: " << param.KBT << std::endl;
    os << "energy: " << param.energy << std::endl;
    os << "energyPrev: " << param.energyPrev << std::endl;
    return os;
}