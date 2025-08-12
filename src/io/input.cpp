#include "io/io.hpp"

using namespace std;

/*
 * pop all spaces & tabs from the given string (\s and \t)
 * return popped string
 * the input string is unchanged
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
	return poppedString;
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
		else
		{
			param.boundaryCondition = BoundaryType::Fixed;
		}
		std::cout << "BOUNDARY_TYPE set to : " << variableValueStr
				  << std::endl;
		return true;
	}
	else if (variableNameStr.compare("maxIterations") == 0)
	{
		param.maxIterations = std::stoi(variableValueStr);
		std::cout << "MAXITERATIONS set to: " << variableValueStr
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
 * @brief Import a mesh from separate files containing vertices and faces.
 *
 * This function reads vertices and faces from specified files and constructs a mesh.
 * 
 * @param mesh The Mesh object to write vertices and faces data.
 * @param verticesFilepath The file path of the vertices file.
 * @param facesFilepath The file path of the faces file.
 * @return True if the mesh is successfully imported, false otherwise.
 */
bool import_mesh_from_vertices_faces(Mesh& mesh, std::string verticesFilepath, std::string facesFilepath){
	char SEPARATOR = ',';
	std::vector<std::vector<int>> facesData = read_data_from_csv<int>(facesFilepath, SEPARATOR);
	std::vector<std::vector<double>> verticesData = read_data_from_csv<double>(verticesFilepath, SEPARATOR);

	// set up and link geometric compnents of the mesh
	mesh.setup_from_vertices_faces(verticesData, facesData);

	return true;
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