#include "io/io.hpp"
#include <cctype>
#include <cmath>
#include <sstream>
#include <stdexcept>

using namespace std;

/*
 * pop all spaces, tabs & carriage returns from the given string
 * (\s, \t and \r)
 * return popped string
 * the input string is unchanged
 *
 * getline() splits on '\n', so a CRLF-terminated line keeps a trailing
 * '\r'. Treating it as whitespace means no exact comparison downstream can be
 * defeated by an invisible character -- notably boundaryType, where 'Periodic\r'
 * matched neither "Periodic" nor "periodic" and fell through to
 * BoundaryType::Fixed: wrong physics, and no error to show for it.
 */
std::string pop_space(std::string rawString)
{
	std::string poppedString = rawString;
	while (poppedString.find(" ") != std::string::npos)
	{
		poppedString.erase(poppedString.find(" "), 1);
	}
	while (poppedString.find("\t") != std::string::npos)
	{
		poppedString.erase(poppedString.find("\t"), 1);
	}
	while (poppedString.find("\r") != std::string::npos)
	{
		poppedString.erase(poppedString.find("\r"), 1);
	}
	return poppedString;
}

std::string trim_trailing_semicolon(std::string rawString)
{
	if (!rawString.empty() && rawString.back() == ';')
	{
		rawString.pop_back();
	}
	return rawString;
}

/*
 * pop trailing carriage returns from the given string
 * return popped string
 * the input string is unchanged
 *
 * getline() splits on '\n', so a CRLF-terminated line keeps a trailing
 * '\r'. Strip it before parsing so an invisible character cannot defeat the
 * exact value comparisons below (e.g. "cpu", "true", "Periodic").
 * .gitattributes keeps checkouts LF; this is the backstop for files that were
 * hand-edited on Windows.
 */
std::string trim_trailing_cr(std::string rawString)
{
	while (!rawString.empty() && rawString.back() == '\r')
	{
		rawString.pop_back();
	}
	return rawString;
}

/*
 * import parameters from key-value strings
 * store value in given Param object
 * pop space before match
 */
bool import_kv_string(std::string variableNameStr, std::string variableValueStr, Param &param)
{
	if (variableNameStr.compare("boundaryType") == 0)
	{
		if (variableValueStr.compare("Periodic") == 0 ||
			variableValueStr.compare("periodic") == 0)
		{
			param.boundaryCondition = BoundaryType::Periodic;
		}
		else if (variableValueStr.compare("Free") == 0 ||
				 variableValueStr.compare("free") == 0)
		{
			param.boundaryCondition = BoundaryType::Free;
		}
		else if (variableValueStr.compare("Mixed") == 0 ||
				 variableValueStr.compare("mixed") == 0)
		{
			param.boundaryCondition = BoundaryType::Mixed;
		}
		else
		{
			param.boundaryCondition = BoundaryType::Fixed;
		}
		std::cout << "BOUNDARY_TYPE set to : " << variableValueStr
				  << std::endl;
		return true;
	}
	else if (variableNameStr.compare("boundaryTypeX") == 0 ||
			 variableNameStr.compare("boundaryTypeY") == 0)
	{
		// The per-axis boundary of the generated sheet under boundaryType =
		// Mixed. Validated here, like forceBackend: a typo should stop the run
		// before it builds a sheet with the wrong edges.
		BoundaryType axisType = BoundaryType::Periodic;
		if (!parse_boundary_type(variableValueStr, axisType) || axisType == BoundaryType::Mixed)
		{
			throw std::runtime_error("[read_param_file] " + variableNameStr +
									 " must be Periodic, Free or Fixed; got '" +
									 variableValueStr + "'");
		}
		if (variableNameStr.back() == 'X')
		{
			param.boundaryConditionX = axisType;
		}
		else
		{
			param.boundaryConditionY = axisType;
		}
		std::cout << variableNameStr << " set to: " << boundary_type_name(axisType) << std::endl;
		return true;
	}
	else if (variableNameStr.compare("fixedBoundaryRings") == 0)
	{
		param.fixedBoundaryRings = std::stoi(variableValueStr);
		if (param.fixedBoundaryRings < 1)
		{
			param.fixedBoundaryRings = 1;
		}
		std::cout << "fixedBoundaryRings set to: " << param.fixedBoundaryRings << std::endl;
		return true;
	}
	else if (variableNameStr.compare("meshVerticesFile") == 0)
	{
		param.meshVerticesFile = variableValueStr;
		std::cout << "meshVerticesFile set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("meshFacesFile") == 0)
	{
		param.meshFacesFile = variableValueStr;
		std::cout << "meshFacesFile set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("maxIterations") == 0)
	{
		param.maxIterations = std::stoi(variableValueStr);
		std::cout << "MAXITERATIONS set to: " << variableValueStr
				  << std::endl;
		return true;
	}
	else if (variableNameStr.compare("restartInputFile") == 0 ||
			 variableNameStr.compare("restartFrom") == 0)
	{
		param.restartInputFile = variableValueStr;
		std::cout << "restartInputFile set to: " << variableValueStr
				  << std::endl;
		return true;
	}
	else if (variableNameStr.compare("checkpointOutputFile") == 0 ||
			 variableNameStr.compare("checkpointFile") == 0)
	{
		param.checkpointOutputFile = variableValueStr;
		std::cout << "checkpointOutputFile set to: " << variableValueStr
				  << std::endl;
		return true;
	}
	else if (variableNameStr.compare("checkpointOutputInterval") == 0 ||
			 variableNameStr.compare("checkpointInterval") == 0)
	{
		param.checkpointOutputInterval = std::stoi(variableValueStr);
		std::cout << "checkpointOutputInterval set to: " << variableValueStr
				  << std::endl;
		return true;
	}
	else if (variableNameStr.compare("meshpointOutput") == 0)
	{
		param.meshpointOutput = (variableValueStr.compare("true") == 0);
		std::cout << "MESHPOINTOUTPUT set to: " << variableValueStr
				  << std::endl;
		return true;
	}
	else if (variableNameStr.compare("meshpointOutputInterval") == 0)
	{
		param.meshpointOutputInterval = std::stoi(variableValueStr);
		if (param.meshpointOutputInterval < 1)
		{
			param.meshpointOutputInterval = 1;
		}
		std::cout << "meshpointOutputInterval set to: " << param.meshpointOutputInterval
				  << std::endl;
		return true;
	}
	else if (variableNameStr.compare("xyzOutput") == 0)
	{
		param.xyzOutput = (variableValueStr.compare("true") == 0);
		std::cout << "XYZOUTPUT set to: " << variableValueStr
				  << std::endl;
		return true;
	}
	else if (variableNameStr.compare("deltaEnergyConverge") == 0)
	{
		param.deltaEnergyConverge = std::stod(variableValueStr);
		std::cout << "deltaEnergyConverge set to: " << variableValueStr
				  << std::endl;
		return true;
	}
	else if (variableNameStr.compare("deltaForceScaleConverge") == 0)
	{
		param.deltaForceScaleConverge = std::stod(variableValueStr);
		std::cout << "deltaForceScaleConverge set to: " << variableValueStr
				  << std::endl;
		return true;
	}
	else if (variableNameStr.compare("usingRpi") == 0)
	{
		// Validated, like forceBackend: this switches which energy the
		// regularization is, so "True" must not quietly mean false.
		if (variableValueStr.compare("true") != 0 && variableValueStr.compare("false") != 0)
		{
			throw std::runtime_error("[read_param_file] usingRpi must be true or false; got '" +
									 variableValueStr + "'");
		}
		param.usingRpi = (variableValueStr.compare("true") == 0);
		std::cout << "usingRpi set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("gamaShape") == 0 ||
			 variableNameStr.compare("gamaArea") == 0)
	{
		// The R-pi thresholds. Zero is legal (every face that deviates at all
		// is in the adaptive case); a negative or non-finite threshold is not
		// a threshold.
		const double threshold = std::stod(variableValueStr);
		if (!std::isfinite(threshold) || threshold < 0.0)
		{
			throw std::runtime_error("[read_param_file] " + variableNameStr +
									 " must be a finite number >= 0; got '" +
									 variableValueStr + "'");
		}
		if (variableNameStr.compare("gamaShape") == 0)
		{
			param.gamaShape = threshold;
		}
		else
		{
			param.gamaArea = threshold;
		}
		std::cout << variableNameStr << " set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("sideX") == 0)
	{
		param.sideX = std::stod(variableValueStr);
		std::cout << "SideX set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("sideY") == 0)
	{
		param.sideY = std::stod(variableValueStr);
		std::cout << "SideY set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("lFace") == 0)
	{
		param.lFace = std::stod(variableValueStr);
		std::cout << "LMESHSIDE set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("c0Insertion") == 0)
	{
		param.insertCurv = std::stod(variableValueStr);
		std::cout << "C0INSERTION set to: " << variableValueStr
				  << std::endl;
		return true;
	}
	else if (variableNameStr.compare("c0Membrane") == 0)
	{
		param.spontCurv = std::stod(variableValueStr);
		std::cout << "C0MEMBRANE set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("kcMembraneBending") == 0)
	{
		param.kCurv = std::stod(variableValueStr);
		std::cout << "KCMEMBRANEBENDING set to: " << variableValueStr
				  << std::endl;
		return true;
	}
	else if (variableNameStr.compare("usMembraneStretching") == 0)
	{
		param.uSurf = std::stod(variableValueStr);
		std::cout << "USMEMBRANESTRETCHING set to: " << variableValueStr
				  << std::endl;
		return true;
	}
	else if (variableNameStr.compare("uvVolumeConstraint") == 0)
	{
		param.uVol = std::stod(variableValueStr);
		std::cout << "UVVOLUMECONSTRAINT set to: " << variableValueStr
				  << std::endl;
		return true;
	}
	else if (variableNameStr.compare("forceBackend") == 0)
	{
		// Validated here rather than at first use: a typo should stop the run
		// before it spends an hour on the wrong backend.
		if (variableValueStr.compare("cpu") == 0 || variableValueStr.compare("gpu") == 0 ||
			variableValueStr.compare("auto") == 0)
		{
			param.forceBackend = variableValueStr;
		}
		else
		{
			throw std::runtime_error("[read_param_file] forceBackend must be cpu, gpu or auto; got '" +
									 variableValueStr + "'");
		}
		std::cout << "FORCE_BACKEND set to : " << variableValueStr
				  << std::endl;
		return true;
	}
	else if (variableNameStr.compare("isGlobalConstraint") == 0)
	{
		if (variableValueStr.compare("true") == 0)
		{
			param.isGlobalConstraint = true;
		}
		else
		{
			param.isGlobalConstraint = false;
		}
		std::cout << "GLOBAL CONSTRAINT for area and volume : " << variableValueStr
				  << std::endl;
		return true;
	}
	else if (variableNameStr.compare("timeStep") == 0)
	{
		param.timeStep = std::stod(variableValueStr);
		std::cout << "TIMESTEP set to: " << variableValueStr
				  << std::endl;
		return true;
	}
	else if (variableNameStr.compare("diffConst") == 0)
	{
		param.diffConst = std::stod(variableValueStr);
		std::cout << "DIFFCONST set to: " << variableValueStr
				  << std::endl;
		return true;
	}
	else if (variableNameStr.compare("KBT") == 0)
	{
		param.KBT = std::stod(variableValueStr);
		std::cout << "KBT set to: " << variableValueStr
				  << std::endl;
		return true;
	}
	else if (variableNameStr.compare("integratePeriodicDuplicates") == 0)
	{
		param.integratePeriodicDuplicates = (variableValueStr.compare("true") == 0);
		std::cout << "integratePeriodicDuplicates set to: " << variableValueStr
				  << std::endl;
		return true;
	}
	else if (variableNameStr.compare("fdtConsistentSurfaceUpdate") == 0)
	{
		param.fdtConsistentSurfaceUpdate = (variableValueStr.compare("true") == 0);
		std::cout << "fdtConsistentSurfaceUpdate set to: " << variableValueStr
				  << std::endl;
		return true;
	}
	else if (variableNameStr.compare("meshForceSubsteps") == 0)
	{
		param.meshForceSubsteps = std::stoi(variableValueStr);
		std::cout << "meshForceSubsteps set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("meshForceSubstepNoise") == 0)
	{
		if (variableValueStr != "once" && variableValueStr != "split")
		{
			throw std::invalid_argument("meshForceSubstepNoise must be once or split, got " +
			                            variableValueStr);
		}
		param.meshForceSubstepNoise = variableValueStr;
		std::cout << "meshForceSubstepNoise set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("inPlaneDynamicsEnabled") == 0)
	{
		param.inPlaneDynamicsEnabled = (variableValueStr.compare("true") == 0);
		std::cout << "inPlaneDynamicsEnabled set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("edgeSpringEnabled") == 0)
	{
		param.edgeSpringEnabled = (variableValueStr.compare("true") == 0);
		std::cout << "edgeSpringEnabled set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("edgeSpringConstant") == 0)
	{
		param.edgeSpringConstant = std::stod(variableValueStr);
		std::cout << "edgeSpringConstant set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("edgeSpringRestLength") == 0)
	{
		param.edgeSpringRestLength = std::stod(variableValueStr);
		std::cout << "edgeSpringRestLength set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("irregularPatchDepths") == 0)
	{
		// Four depths, for valences 4, 5, 7 and 8. pop_space() has already
		// taken the spaces out, so commas separate them; brackets are allowed
		// so `[12, 8, 7, 7]` reads as well.
		std::string list = variableValueStr;
		for (char &c : list)
		{
			if (c == '[' || c == ']' || c == ';')
			{
				c = ',';
			}
		}
		std::vector<int> depths;
		std::stringstream stream(list);
		std::string item;
		while (std::getline(stream, item, ','))
		{
			if (item.empty())
			{
				continue;
			}
			std::size_t used = 0;
			int depth = 0;
			try
			{
				depth = std::stoi(item, &used);
			}
			catch (const std::exception &)
			{
				used = 0;
			}
			if (used != item.size() || depth < 1 || depth > 24)
			{
				throw std::runtime_error(
					"[read_param_file] irregularPatchDepths entries must be integers in [1, 24]; got '" +
					item + "'");
			}
			depths.push_back(depth);
		}
		if (depths.size() != param.irregularPatchDepths.size())
		{
			throw std::runtime_error(
				"[read_param_file] irregularPatchDepths needs four depths, for valences 4, 5, 7 "
				"and 8; got '" + variableValueStr + "'");
		}
		for (std::size_t i = 0; i < depths.size(); i++)
		{
			param.irregularPatchDepths[i] = depths[i];
		}
		std::cout << "irregularPatchDepths set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("irregularPatchTailCorrection") == 0)
	{
		// Validated: it changes the Hamiltonian, so "True" must not quietly
		// mean false.
		if (variableValueStr.compare("true") != 0 && variableValueStr.compare("false") != 0)
		{
			throw std::runtime_error(
				"[read_param_file] irregularPatchTailCorrection must be true or false; got '" +
				variableValueStr + "'");
		}
		param.irregularPatchTailCorrection = (variableValueStr.compare("true") == 0);
		std::cout << "irregularPatchTailCorrection set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("irregularPatchDepthScale") == 0)
	{
		param.irregularPatchDepthScale = std::stod(variableValueStr);
		std::cout << "irregularPatchDepthScale set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("edgeTetherShape") == 0)
	{
		param.edgeTetherShape = variableValueStr;
		std::cout << "edgeTetherShape set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("edgeTetherMinRatio") == 0)
	{
		param.edgeTetherMinRatio = std::stod(variableValueStr);
		std::cout << "edgeTetherMinRatio set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("edgeTetherMaxRatio") == 0)
	{
		param.edgeTetherMaxRatio = std::stod(variableValueStr);
		std::cout << "edgeTetherMaxRatio set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("triangleShapeEnabled") == 0)
	{
		param.triangleShapeEnabled = (variableValueStr.compare("true") == 0);
		std::cout << "triangleShapeEnabled set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("triangleShapeMinAltitudeRatio") == 0)
	{
		param.triangleShapeMinAltitudeRatio = std::stod(variableValueStr);
		std::cout << "triangleShapeMinAltitudeRatio set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("triangleShapeConstant") == 0)
	{
		param.triangleShapeConstant = std::stod(variableValueStr);
		std::cout << "triangleShapeConstant set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("creaseWallEnabled") == 0)
	{
		param.creaseWallEnabled = (variableValueStr.compare("true") == 0);
		std::cout << "creaseWallEnabled set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("creaseWallAngle") == 0)
	{
		param.creaseWallAngle = std::stod(variableValueStr);
		std::cout << "creaseWallAngle set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("creaseWallConstant") == 0)
	{
		param.creaseWallConstant = std::stod(variableValueStr);
		std::cout << "creaseWallConstant set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("meshWallType") == 0)
	{
		param.meshWallType = variableValueStr;
		std::cout << "meshWallType set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("vertexMoveEnabled") == 0)
	{
		param.vertexMoveEnabled = (variableValueStr.compare("true") == 0);
		std::cout << "vertexMoveEnabled set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("vertexMoveStepRatio") == 0)
	{
		param.vertexMoveStepRatio = std::stod(variableValueStr);
		std::cout << "vertexMoveStepRatio set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("vertexMoveAttemptsPerVertex") == 0)
	{
		param.vertexMoveAttemptsPerVertex = std::stod(variableValueStr);
		std::cout << "vertexMoveAttemptsPerVertex set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("vertexMoveInterval") == 0)
	{
		param.vertexMoveInterval = std::stoi(variableValueStr);
		std::cout << "vertexMoveInterval set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("brownianStepEnabled") == 0)
	{
		param.brownianStepEnabled = (variableValueStr.compare("true") == 0);
		std::cout << "brownianStepEnabled set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("surfaceSolver") == 0)
	{
		param.surfaceSolver = variableValueStr;
		std::cout << "surfaceSolver set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("edgeFlipEnabled") == 0)
	{
		param.edgeFlipEnabled = (variableValueStr.compare("true") == 0);
		std::cout << "edgeFlipEnabled set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("edgeFlipAttemptRate") == 0)
	{
		param.edgeFlipAttemptRate = std::stod(variableValueStr);
		std::cout << "edgeFlipAttemptRate set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("edgeFlipInterval") == 0)
	{
		param.edgeFlipInterval = std::stoi(variableValueStr);
		std::cout << "edgeFlipInterval set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("edgeFlipMinValence") == 0)
	{
		param.edgeFlipMinValence = std::stoi(variableValueStr);
		std::cout << "edgeFlipMinValence set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("edgeFlipMaxValence") == 0)
	{
		param.edgeFlipMaxValence = std::stoi(variableValueStr);
		std::cout << "edgeFlipMaxValence set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("edgeFlipIncrementalForces") == 0)
	{
		param.edgeFlipIncrementalForces = (variableValueStr.compare("true") == 0);
		std::cout << "edgeFlipIncrementalForces set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("randomSeed") == 0)
	{
		param.randomSeed = static_cast<unsigned int>(std::stoul(variableValueStr));
		std::cout << "randomSeed set to: " << variableValueStr
				  << std::endl;
		return true;
	}
	else if (variableNameStr.compare("thermalFluctuationEnabled") == 0)
	{
		param.thermalFluctuationEnabled = (variableValueStr.compare("true") == 0);
		std::cout << "thermalFluctuationEnabled set to: " << variableValueStr
				  << std::endl;
		return true;
	}
	else if (variableNameStr.compare("thermalFluctuationPureMMC") == 0)
	{
		param.thermalFluctuationPureMMC = (variableValueStr.compare("true") == 0);
		std::cout << "thermalFluctuationPureMMC set to: " << variableValueStr
				  << std::endl;
		return true;
	}
	else if (variableNameStr.compare("thermalFluctuationInterval") == 0)
	{
		param.thermalFluctuationInterval = std::stoi(variableValueStr);
		std::cout << "thermalFluctuationInterval set to: " << variableValueStr
				  << std::endl;
		return true;
	}
	else if (variableNameStr.compare("thermalFluctuationTemperatureKelvin") == 0)
	{
		param.thermalFluctuationTemperatureKelvin = std::stod(variableValueStr);
		std::cout << "thermalFluctuationTemperatureKelvin set to: " << variableValueStr
				  << std::endl;
		return true;
	}
	else if (variableNameStr.compare("thermalFluctuationMinTemperatureKelvin") == 0)
	{
		param.thermalFluctuationMinTemperatureKelvin = std::stod(variableValueStr);
		std::cout << "thermalFluctuationMinTemperatureKelvin set to: " << variableValueStr
				  << std::endl;
		return true;
	}
	else if (variableNameStr.compare("thermalFluctuationCoolingRate") == 0)
	{
		param.thermalFluctuationCoolingRate = std::stod(variableValueStr);
		std::cout << "thermalFluctuationCoolingRate set to: " << variableValueStr
				  << std::endl;
		return true;
	}
	else if (variableNameStr.compare("thermalFluctuationStepScale") == 0)
	{
		param.thermalFluctuationStepScale = std::stod(variableValueStr);
		std::cout << "thermalFluctuationStepScale set to: " << variableValueStr
				  << std::endl;
		return true;
	}
	else if (variableNameStr.compare("setRelaxAreaToDefault") == 0)
	{
		if (variableValueStr.compare("true") == 0)
		{
			param.setRelaxAreaToDefault = true;
		}
		else
		{
			param.setRelaxAreaToDefault = false;
		}
		std::cout << "setRelaxAreaToDefault set to : " << variableValueStr
				  << std::endl;
		return true;
	}
	else if (variableNameStr.compare("relaxArea") == 0)
	{
		param.area0 = std::stod(variableValueStr);
		std::cout << "relaxArea set to: " << variableValueStr
				  << std::endl;
		return true;
	}
	else if (variableNameStr.compare("scaffoldingFileName") == 0)
	{
		param.scaffoldingFileName = variableValueStr;
		std::cout << "scaffoldingFileName set to: " << variableValueStr
				  << std::endl;
		return true;
	}
	else if (variableNameStr.compare("isEnergyHarmonicBondIncluded") == 0)
	{
		if (variableValueStr.compare("true") == 0)
		{
			param.isEnergyHarmonicBondIncluded = true;
		}
		else
		{
			param.isEnergyHarmonicBondIncluded = false;
		}
		std::cout << "HARMONIC BOND included : " << variableValueStr
				  << std::endl;
		return true;
	}
	else if (variableNameStr.compare("relaxLengthRatioApproximation") == 0)
	{
		param.relaxLengthRatioApproximation = std::stod(variableValueStr);
		std::cout << "relaxLengthRatioApproximation set to: " << variableValueStr
				  << std::endl;
		return true;
	}
	else if (variableNameStr.compare("scaffoldingZeroPlaneZ") == 0)
	{
		param.scaffoldingZeroPlaneZ = std::stod(variableValueStr);
		std::cout << "scaffoldingZeroPlaneZ set to: " << variableValueStr
				  << std::endl;
		return true;
	}
	else if (variableNameStr.compare("isGagScaffoldingEnergyIncluded") == 0)
	{
		param.isGagScaffoldingEnergyIncluded = (variableValueStr.compare("true") == 0);
		std::cout << "GAG SCAFFOLDING ENERGY included : " << variableValueStr
				  << std::endl;
		return true;
	}
	else if (variableNameStr.compare("gagReferenceStateFileName") == 0)
	{
		param.gagReferenceStateFileName = variableValueStr;
		std::cout << "gagReferenceStateFileName set to: " << variableValueStr
				  << std::endl;
		return true;
	}
	else if (variableNameStr.compare("gagReactionFileName") == 0)
	{
		param.gagReactionFileName = variableValueStr;
		std::cout << "gagReactionFileName set to: " << variableValueStr
				  << std::endl;
		return true;
	}
	else if (variableNameStr.compare("isIdealizedProteinLatticeEnergyIncluded") == 0)
	{
		param.isIdealizedProteinLatticeEnergyIncluded = (variableValueStr.compare("true") == 0);
		std::cout << "IDEALIZED PROTEIN LATTICE ENERGY included : " << variableValueStr
				  << std::endl;
		return true;
	}
	else if (variableNameStr.compare("idealizedProteinLatticeFileName") == 0)
	{
		param.idealizedProteinLatticeFileName = variableValueStr;
		std::cout << "idealizedProteinLatticeFileName set to: " << variableValueStr
				  << std::endl;
		return true;
	}
	else if (variableNameStr.compare("gagKsigma") == 0)
	{
		param.gagKsigma = std::stod(variableValueStr);
		std::cout << "gagKsigma set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("gagKtheta") == 0)
	{
		param.gagKtheta = std::stod(variableValueStr);
		std::cout << "gagKtheta set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("gagKphi") == 0)
	{
		param.gagKphi = std::stod(variableValueStr);
		std::cout << "gagKphi set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("gagKomega") == 0)
	{
		param.gagKomega = std::stod(variableValueStr);
		std::cout << "gagKomega set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("gagFiniteDifferenceStep") == 0)
	{
		param.gagFiniteDifferenceStep = std::stod(variableValueStr);
		std::cout << "gagFiniteDifferenceStep set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("gagPropagationStepSize") == 0)
	{
		param.gagPropagationStepSize = std::stod(variableValueStr);
		std::cout << "gagPropagationStepSize set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("gagPreRelaxSteps") == 0)
	{
		param.gagPreRelaxSteps = std::stoi(variableValueStr);
		std::cout << "gagPreRelaxSteps set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("propagateScaffoldingInterv") == 0)
	{
		param.propagateScaffoldingInterv = std::stoi(variableValueStr);
		std::cout << "propagateScaffoldingInterv set to: " << variableValueStr
				  << std::endl;
		return true;
	}
	else if (variableNameStr.compare("propagateScaffoldingNstep") == 0)
	{
		param.propagateScaffoldingNstep = std::stoi(variableValueStr);
		std::cout << "propagateScaffoldingNstep set to: " << variableValueStr
				  << std::endl;
		return true;
	}
	// External force on a region of the membrane; see docs/external_force.md.
	else if (variableNameStr.compare("externalForceEnabled") == 0)
	{
		param.externalForceEnabled = (variableValueStr.compare("true") == 0);
		std::cout << "externalForceEnabled set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("externalForceX") == 0)
	{
		param.externalForceX = std::stod(variableValueStr);
		std::cout << "externalForceX set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("externalForceY") == 0)
	{
		param.externalForceY = std::stod(variableValueStr);
		std::cout << "externalForceY set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("externalForceZ") == 0)
	{
		param.externalForceZ = std::stod(variableValueStr);
		std::cout << "externalForceZ set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("externalForceShape") == 0)
	{
		param.externalForceShape = variableValueStr;
		std::cout << "externalForceShape set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("externalForceCenterX") == 0)
	{
		param.externalForceCenterX = std::stod(variableValueStr);
		std::cout << "externalForceCenterX set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("externalForceCenterY") == 0)
	{
		param.externalForceCenterY = std::stod(variableValueStr);
		std::cout << "externalForceCenterY set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("externalForceRadius") == 0)
	{
		param.externalForceRadius = std::stod(variableValueStr);
		std::cout << "externalForceRadius set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("externalForceHalfWidthX") == 0)
	{
		param.externalForceHalfWidthX = std::stod(variableValueStr);
		std::cout << "externalForceHalfWidthX set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("externalForceHalfWidthY") == 0)
	{
		param.externalForceHalfWidthY = std::stod(variableValueStr);
		std::cout << "externalForceHalfWidthY set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("externalForceStartIteration") == 0)
	{
		param.externalForceStartIteration = std::stoi(variableValueStr);
		std::cout << "externalForceStartIteration set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("externalForceStopIteration") == 0)
	{
		param.externalForceStopIteration = std::stoi(variableValueStr);
		std::cout << "externalForceStopIteration set to: " << variableValueStr << std::endl;
		return true;
	}
	// A rigid bead pushed through the membrane; see docs/bead_push_benchmark.md.
	else if (variableNameStr.compare("beadEnabled") == 0)
	{
		param.beadEnabled = (variableValueStr.compare("true") == 0);
		std::cout << "beadEnabled set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("beadRadius") == 0)
	{
		param.beadRadius = std::stod(variableValueStr);
		std::cout << "beadRadius set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("beadStiffness") == 0)
	{
		param.beadStiffness = std::stod(variableValueStr);
		std::cout << "beadStiffness set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("beadCenterX") == 0)
	{
		param.beadCenterX = std::stod(variableValueStr);
		std::cout << "beadCenterX set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("beadCenterY") == 0)
	{
		param.beadCenterY = std::stod(variableValueStr);
		std::cout << "beadCenterY set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("beadCenterZ") == 0)
	{
		param.beadCenterZ = std::stod(variableValueStr);
		std::cout << "beadCenterZ set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("beadVelocityX") == 0)
	{
		param.beadVelocityX = std::stod(variableValueStr);
		std::cout << "beadVelocityX set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("beadVelocityY") == 0)
	{
		param.beadVelocityY = std::stod(variableValueStr);
		std::cout << "beadVelocityY set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("beadVelocityZ") == 0)
	{
		param.beadVelocityZ = std::stod(variableValueStr);
		std::cout << "beadVelocityZ set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("beadStartIteration") == 0)
	{
		param.beadStartIteration = std::stoi(variableValueStr);
		std::cout << "beadStartIteration set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("beadStopIteration") == 0)
	{
		param.beadStopIteration = std::stoi(variableValueStr);
		std::cout << "beadStopIteration set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("beadAdhesionEnergy") == 0)
	{
		param.beadAdhesionEnergy = std::stod(variableValueStr);
		std::cout << "beadAdhesionEnergy set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("beadAdhesionRange") == 0)
	{
		param.beadAdhesionRange = std::stod(variableValueStr);
		std::cout << "beadAdhesionRange set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("beadPotential") == 0)
	{
		param.beadPotential = variableValueStr;
		std::cout << "beadPotential set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("beadCaptureRange") == 0)
	{
		param.beadCaptureRange = std::stod(variableValueStr);
		std::cout << "beadCaptureRange set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("beadBreakRange") == 0)
	{
		param.beadBreakRange = std::stod(variableValueStr);
		std::cout << "beadBreakRange set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("noiseTruncation") == 0)
	{
		// 0 is off. Below 0.1 the rejection loop that draws the truncated
		// normal accepts fewer than 8 % of its attempts, and no cutoff that
		// small is a model of thermal noise anyone means.
		const double cutoff = std::stod(variableValueStr);
		if (!std::isfinite(cutoff) || cutoff < 0.0 || (cutoff > 0.0 && cutoff < 0.1))
		{
			throw std::runtime_error("[read_param_file] noiseTruncation must be 0 (off) or a finite "
									 "number >= 0.1; got '" + variableValueStr + "'");
		}
		param.noiseTruncation = cutoff;
		std::cout << "noiseTruncation set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("noiseTruncationRenormalize") == 0)
	{
		param.noiseTruncationRenormalize = (variableValueStr.compare("true") == 0);
		std::cout << "noiseTruncationRenormalize set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("meanHeightPinned") == 0)
	{
		param.meanHeightPinned = (variableValueStr.compare("true") == 0);
		std::cout << "meanHeightPinned set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("elasticBoxEnabled") == 0)
	{
		param.elasticBoxEnabled = (variableValueStr.compare("true") == 0);
		std::cout << "elasticBoxEnabled set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("elasticBoxTension") == 0)
	{
		param.elasticBoxTension = std::stod(variableValueStr);
		std::cout << "elasticBoxTension set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("elasticBoxMobility") == 0)
	{
		param.elasticBoxMobility = std::stod(variableValueStr);
		std::cout << "elasticBoxMobility set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("elasticBoxThermal") == 0)
	{
		param.elasticBoxThermal = (variableValueStr.compare("true") == 0);
		std::cout << "elasticBoxThermal set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("elasticBoxCoupling") == 0)
	{
		param.elasticBoxCoupling = variableValueStr;
		std::cout << "elasticBoxCoupling set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("elasticBoxMaxStrainPerStep") == 0)
	{
		param.elasticBoxMaxStrainPerStep = std::stod(variableValueStr);
		std::cout << "elasticBoxMaxStrainPerStep set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("elasticBoxVertexGasTerm") == 0)
	{
		param.elasticBoxVertexGasTerm = (variableValueStr.compare("true") == 0);
		std::cout << "elasticBoxVertexGasTerm set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("lipidFlowMode") == 0)
	{
		param.lipidFlowMode = variableValueStr;
		std::cout << "lipidFlowMode set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("lipidFlowDiffusivity") == 0)
	{
		param.lipidFlowDiffusivity = std::stod(variableValueStr);
		std::cout << "lipidFlowDiffusivity set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("lipidFlowSubsteps") == 0)
	{
		param.lipidFlowSubsteps = std::stoi(variableValueStr);
		std::cout << "lipidFlowSubsteps set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("lipidFlowThermal") == 0)
	{
		param.lipidFlowThermal = (variableValueStr.compare("true") == 0);
		std::cout << "lipidFlowThermal set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("lipidFlowLinearStretching") == 0)
	{
		param.lipidFlowLinearStretching = (variableValueStr.compare("true") == 0);
		std::cout << "lipidFlowLinearStretching set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("lipidReservoirRegion") == 0)
	{
		param.lipidReservoirRegion = variableValueStr;
		std::cout << "lipidReservoirRegion set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("lipidFlipResampling") == 0)
	{
		param.lipidFlipResampling = variableValueStr;
		std::cout << "lipidFlipResampling set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("lipidReservoirTension") == 0)
	{
		param.lipidReservoirTension = std::stod(variableValueStr);
		std::cout << "lipidReservoirTension set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("lipidReservoirRate") == 0)
	{
		param.lipidReservoirRate = std::stod(variableValueStr);
		std::cout << "lipidReservoirRate set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("localOutput") == 0)
	{
		param.localOutput = (variableValueStr.compare("true") == 0);
		std::cout << "localOutput set to: " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("preRefineMesh") == 0)
	{
		param.isPreRefinementEnabled = (variableValueStr.compare("true") == 0);
		std::cout << "preRefineMesh set to : " << variableValueStr << std::endl;
		return true;
	}
	else if (variableNameStr.compare("VERBOSE_MODE") == 0)
	{
		if (variableValueStr.compare("true") == 0)
		{
			param.VERBOSE_MODE = true;
		}
		else
		{
			param.VERBOSE_MODE = false;
		}
		std::cout << "VERBOSED_MODE set to : " << variableValueStr
				  << std::endl;
		return true;
	}
	std::cout << "VARIABLE NOT SUPPORTED: " << variableNameStr << std::endl;
	return false;
}

/*
 * load parameter file
 * input file path
 */
bool import_param_file(Param &param, std::string filepath)
{
	// load in parameter file with the given filename
	std::ifstream ifile(filepath, std::ios::in);
	std::vector<std::string> parameters; // convert params file data to vector of rows

	// check to see that the file was opened correctly:
	if (!ifile.is_open())
	{
		std::cout << "There was a problem opening the parameter file!\n"
				  << endl;
		// exit(1); //exit or do additional error checking
	}

	std::string str = "";
	// keep reading line by line from the text file so long as data exists:
	//(delete rows starting with # )
	while (getline(ifile, str))
	{
		str = trim_trailing_cr(str);
		if (str[0] != '#')
		{
			parameters.push_back(str);
		}
	}

	// store kv pair
	std::string variableNameStr = "";
	std::string variableValueStr = "";

	// iterate over rows
	for (int i = 0; i < parameters.size(); ++i)
	{

		// use "=" and "#" as marks to find name-value pairs
		if (parameters[i].find("=") != std::string::npos)
		{

			variableNameStr = parameters[i].substr(0, parameters[i].find("=")); // raw variable name string! need to pop space!
			variableValueStr = parameters[i].substr(
				parameters[i].find("=") + 1); // need to delete comment
			// std::cout << "RAW VALUE: " << variableValueStr << std::endl; //for testing only
			// std::cout << variableValueStr.find("#") << std::endl;//for testing only
			if (variableValueStr.find("#") != std::string::npos)
			{
				variableValueStr = variableValueStr.substr(0,
														   variableValueStr.find("#")); // raw value string! need to pop space!
			}
			// std::cout << "UNCOMMENTED VALUE: " << variableValueStr << std::endl; //for testing only

			// pop space
			// std::cout << variableNameStr << "::" << variableValueStr << std::endl;//for testing only
			variableNameStr = pop_space(variableNameStr);
			variableValueStr = pop_space(variableValueStr);
			variableValueStr = trim_trailing_semicolon(variableValueStr);
			// std::cout << variableNameStr << "::" << variableValueStr << std::endl;//for testing only

			// import kv string and load values to local variables
			import_kv_string(variableNameStr, variableValueStr, param);
		}
	}

	// End of import
	std::cout
		<< "============================END OF INPUT============================"
		<< std::endl;
	return true;
}

/**
 * @brief Read a vertex file and a faces file into MeshFileData.
 *
 * The parsing half of import_mesh_from_vertices_faces(); see io.hpp for the
 * file format. Coordinates, an optional type and mirror per vertex, three
 * corners and an optional copy flag per face.
 * 
 * @param mesh The Mesh object to write vertices and faces data.
 * @param verticesFilepath The file path of the vertices file.
 * @param facesFilepath The file path of the faces file.
 * @return True if the mesh is successfully imported, false otherwise.
 */
namespace
{
/// Strip spaces, tabs and carriage returns from both ends.
std::string trim_field(const std::string &raw)
{
	const std::string blank = " \t\r";
	const std::size_t first = raw.find_first_not_of(blank);
	if (first == std::string::npos)
	{
		return "";
	}
	const std::size_t last = raw.find_last_not_of(blank);
	return raw.substr(first, last - first + 1);
}

/**
 * The rows of a mesh file: fields split on commas (or on whitespace when a
 * line has no comma), trimmed, with blank lines and '#' comment lines skipped
 * and a trailing empty field dropped -- the legacy vertex-type writer ends
 * every line with a comma. @p lineNumbers gets the file line of each row, for
 * error messages.
 */
std::vector<std::vector<std::string>> read_mesh_file_rows(const std::string &filepath,
														  std::vector<int> &lineNumbers)
{
	std::ifstream file(filepath);
	if (!file.is_open())
	{
		throw std::invalid_argument("[read_mesh_vertices_faces_files] Unable to open " + filepath);
	}
	std::vector<std::vector<std::string>> rows;
	std::string line;
	int lineNumber = 0;
	while (std::getline(file, line))
	{
		lineNumber++;
		const std::string stripped = trim_field(line);
		if (stripped.empty() || stripped[0] == '#')
		{
			continue;
		}
		std::vector<std::string> fields;
		std::istringstream stream(stripped);
		std::string field;
		if (stripped.find(',') == std::string::npos)
		{
			while (stream >> field)
			{
				fields.push_back(field);
			}
		}
		else
		{
			while (std::getline(stream, field, ','))
			{
				fields.push_back(trim_field(field));
			}
		}
		while (!fields.empty() && fields.back().empty())
		{
			fields.pop_back();
		}
		rows.push_back(fields);
		lineNumbers.push_back(lineNumber);
	}
	return rows;
}

std::string mesh_file_where(const std::string &filepath, int lineNumber)
{
	return filepath + ":" + std::to_string(lineNumber);
}

double parse_mesh_double(const std::string &text, const std::string &filepath, int lineNumber,
						 const char *what)
{
	try
	{
		std::size_t consumed = 0;
		const double value = std::stod(text, &consumed);
		if (consumed == text.size())
		{
			return value;
		}
	}
	catch (const std::exception &)
	{
	}
	throw std::invalid_argument("[read_mesh_vertices_faces_files] " +
								mesh_file_where(filepath, lineNumber) + ": cannot read " + what +
								" from '" + text + "'");
}

int parse_mesh_int(const std::string &text, const std::string &filepath, int lineNumber,
				   const char *what)
{
	try
	{
		std::size_t consumed = 0;
		const long value = std::stol(text, &consumed);
		if (consumed == text.size())
		{
			return static_cast<int>(value);
		}
	}
	catch (const std::exception &)
	{
	}
	throw std::invalid_argument("[read_mesh_vertices_faces_files] " +
								mesh_file_where(filepath, lineNumber) + ": cannot read " + what +
								" from '" + text + "'");
}
} // namespace

MeshFileData read_mesh_vertices_faces_files(const std::string &verticesFilepath,
											const std::string &facesFilepath)
{
	MeshFileData data;

	// Vertices: x, y, z and optionally a type and a mirror.
	std::vector<int> lineNumbers;
	const std::vector<std::vector<std::string>> vertexRows =
		read_mesh_file_rows(verticesFilepath, lineNumbers);
	const int nVertices = static_cast<int>(vertexRows.size());
	std::vector<VertexType> types(nVertices, VertexType::Free);
	std::vector<int> mirrors(nVertices, -1);
	bool anyType = false;
	data.vertices.reserve(nVertices);
	for (int row = 0; row < nVertices; row++)
	{
		const std::vector<std::string> &fields = vertexRows[row];
		const int lineNumber = lineNumbers[row];
		if (fields.size() < 3)
		{
			throw std::invalid_argument("[read_mesh_vertices_faces_files] " +
										mesh_file_where(verticesFilepath, lineNumber) +
										": a vertex needs x, y and z; got " +
										std::to_string(fields.size()) + " field(s)");
		}
		data.vertices.push_back({parse_mesh_double(fields[0], verticesFilepath, lineNumber, "x"),
								 parse_mesh_double(fields[1], verticesFilepath, lineNumber, "y"),
								 parse_mesh_double(fields[2], verticesFilepath, lineNumber, "z")});
		if (fields.size() >= 4)
		{
			anyType = true;
			if (!parse_vertex_type(fields[3], types[row]))
			{
				throw std::invalid_argument(
					"[read_mesh_vertices_faces_files] " + mesh_file_where(verticesFilepath, lineNumber) +
					": unknown vertex type '" + fields[3] +
					"'; expected free, fixed, periodic or ghost");
			}
			if (fields.size() >= 5)
			{
				mirrors[row] = parse_mesh_int(fields[4], verticesFilepath, lineNumber, "the mirror index");
			}
			if (types[row] == VertexType::Periodic && mirrors[row] < 0)
			{
				throw std::invalid_argument(
					"[read_mesh_vertices_faces_files] " + mesh_file_where(verticesFilepath, lineNumber) +
					": a periodic vertex needs the index of the vertex it mirrors in the fifth column");
			}
		}
	}
	if (anyType)
	{
		data.types = types;
		data.mirrors = mirrors;
	}

	// Faces: three corners and optionally a copy flag.
	lineNumbers.clear();
	const std::vector<std::vector<std::string>> faceRows =
		read_mesh_file_rows(facesFilepath, lineNumbers);
	const int nFaces = static_cast<int>(faceRows.size());
	std::vector<char> flags(nFaces, 0);
	bool anyFlag = false;
	data.faces.reserve(nFaces);
	for (int row = 0; row < nFaces; row++)
	{
		const std::vector<std::string> &fields = faceRows[row];
		const int lineNumber = lineNumbers[row];
		if (fields.size() < 3)
		{
			throw std::invalid_argument("[read_mesh_vertices_faces_files] " +
										mesh_file_where(facesFilepath, lineNumber) +
										": a face needs three vertex indices; got " +
										std::to_string(fields.size()) + " field(s)");
		}
		std::vector<int> corners(3);
		for (int k = 0; k < 3; k++)
		{
			corners[k] = parse_mesh_int(fields[k], facesFilepath, lineNumber, "a vertex index");
			if (corners[k] < 0 || corners[k] >= nVertices)
			{
				throw std::invalid_argument(
					"[read_mesh_vertices_faces_files] " + mesh_file_where(facesFilepath, lineNumber) +
					": vertex index " + std::to_string(corners[k]) + " is outside the " +
					std::to_string(nVertices) + " vertices of " + verticesFilepath);
			}
		}
		data.faces.push_back(corners);
		if (fields.size() >= 4)
		{
			anyFlag = true;
			std::string flag;
			for (char c : fields[3])
			{
				flag.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
			}
			if (flag == "0" || flag == "real" || flag == "false")
			{
				flags[row] = 0;
			}
			else if (flag == "1" || flag == "copy" || flag == "ghost" || flag == "image" ||
					 flag == "true")
			{
				flags[row] = 1;
			}
			else
			{
				throw std::invalid_argument(
					"[read_mesh_vertices_faces_files] " + mesh_file_where(facesFilepath, lineNumber) +
					": unknown face flag '" + fields[3] + "'; expected 0/real or 1/copy");
			}
		}
	}
	if (anyFlag)
	{
		data.faceIsCopy = flags;
	}

	std::cout << "[read_mesh_vertices_faces_files] Read " << nVertices << " vertices from "
			  << verticesFilepath << (anyType ? " (with boundary types)" : "") << " and " << nFaces
			  << " faces from " << facesFilepath << (anyFlag ? " (with copy flags)" : "")
			  << std::endl;
	return data;
}

/*
 * load in model mesh file for adhesion of the triangular mesh
 * to the model mesh via adding an extra energy term -- E adhesion_geometry
 * and save the mesh infomation (vector<vector<double>> (n,3)) in model_mesh
 *
 * input file path in the format of .csv (assuming no endline comma):
 * -- "x, y, z"
 *
 */
vector<Matrix> import_scaffolding_mesh(std::string filepath)
{
	// load in parameter file with the given filename
	std::ifstream ifile(filepath, std::ios::in);
	std::vector<std::string> meshdata; // convert params file data to vector of rows

	// check to see that the file was opened correctly:
	if (!ifile.is_open())
	{
		std::cerr << "There was a problem opening the parameter file!\n";
		exit(1); // exit or do additional error checking
	}

	std::string str = "";
	// keep reading line by line from the text file so long as data exists:
	// pop all spaces and tabs in the process
	//(delete rows starting with # )
	while (getline(ifile, str))
	{
		str = trim_trailing_cr(str);
		if (str[0] != '#' && str[0] != 'x')
		{
			str = pop_space(str);
			meshdata.push_back(str);
		}
	}

	// initialize vector of matrices
	vector<Matrix> model_mesh(meshdata.size(), Matrix(3, 1));

	// iterate over rows of meshdata
	// strings of mesh data values along x, y, z-axis
	std::string x_str = "0";
	std::string y_str = "0";
	std::string z_str = "0";
	// index of comma
	int comma_index = 0;
	for (int i = 0; i < meshdata.size(); ++i)
	{

		// search for comma (all spaces and tabs were popped) to find values
		// assuming no spaces and tabs in the data right now
		// stoi error if not enough values! (fewer than 3)
		// ignore the fourth value and all values afterwards

		// x
		comma_index = meshdata[i].find(",");
		x_str = meshdata[i].substr(0, comma_index);
		model_mesh[i].set(0, 0, std::stod(x_str));
		meshdata[i] = meshdata[i].substr(comma_index + 1);
		// y
		comma_index = meshdata[i].find(",");
		y_str = meshdata[i].substr(0, comma_index);
		model_mesh[i].set(1, 0, std::stod(y_str));
		// z
		z_str = meshdata[i].substr(comma_index + 1);
		model_mesh[i].set(2, 0, std::stod(z_str));
	}

	return model_mesh;
}
