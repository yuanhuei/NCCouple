#include <iostream>
#include <algorithm>
#include "MOCMesh.h"
#include "MOCIndex.h"
#include "CFDMesh.h"
#include "Solver.h"
#include "Logger.h"
#include "MHT_polyhedron/PolyhedronSet.h"
#include<time.h>
#include "./FileConvertor.h"
#include "./ConfigurationFile.h"
#include <string>
#include <sstream>
#ifdef _WIN32
#include <process.h>
#else
#include <unistd.h>
#endif
#include "./MHT_common/SystemControl.h"
#include "./MHT_mesh/UnGridFactory.h"
#include "./MHT_mesh/RegionConnection.h"
#include "./MHT_mesh/Mesh.h"
#include "./MHT_field/Field.h"
#include "./MHT_IO/FieldIO.h"

int g_iProcessID = 0;


Scalar initialT(Scalar x, Scalar y, Scalar z)
{
	Vector PointCenter(x, y, z);
	Vector axisCenter(0.63, 0.63, 0.0);
	Vector axisNorm(0.0, 0.0, 1.0);
	Vector OP = PointCenter - axisCenter;
	Vector axicialLocation = (OP & axisNorm) * axisNorm;
	Scalar radius = (OP - axicialLocation).Mag();
	Scalar sigma = 1.0;
	Scalar temperature = (1000.0 / sqrt(2.0 * PI) / sigma) * exp(-radius * radius / 2 / pow(sigma, 2));
	return temperature;
}

Scalar initialRho(Scalar x, Scalar y, Scalar z)
{
	Scalar add = 100.0 * z * (5.0 - z);
	return 1000.0 + add;
}

//this example was designed for test of
//(1) rewritting a apl file
//(2) reading and writting of inp files
void MOC_APL_INP_FileTest()
{
	WarningContinue("MOC_APL_INP_FileTest");
	MOCMesh mocMesh("pin_c1.apl", "pin_c1.inp",MeshKernelType::MHT_KERNEL);
	mocMesh.WriteTecplotFile("H2O", "h2O.plt");
	mocMesh.WriteTecplotFile("Zr4", "zr4.plt");
	mocMesh.WriteTecplotFile("UO2", "u2o.plt");
	//mocMesh.InitMOCValue("pin_c1.inp","pin_c1.txt");
	//mocMesh.OutputStatus("pin_c1_out.inp");
	return;
}

#include "./MHT_IO/VTKIO.h"
#include <fstream>

void VTKReaderTest()
{
	MHTVTKReader reader("pipe.msh");
	std::vector<std::string> fileName;
	fileName.push_back("pipe_solid.vtk");
	fileName.push_back("pipe_fluid.vtk");
	std::vector<int> IDList;
	IDList.push_back(0);
	IDList.push_back(1);
	std::vector<std::string> fieldName;
	fieldName.push_back("temperature");
	fieldName.push_back("Pressure");
	reader.ReadVTKFile(fileName, IDList, fieldName);
	
	for (size_t i = 0; i < reader.GetFieldIOList().size(); i++)
	{
		reader.GetFieldIO(i).WriteTecplotField("pipe_"+std::to_string(i)+".plt");
		reader.GetFieldIO(i).WriteVTKField("pipe_" + std::to_string(i) + ".vtk");
	}

	Field<Scalar> temperatureField(reader.GetMeshListPtr()[0],0.0,"temperature");
	temperatureField.ReadVTK_Field("pipe_solid.vtk");
	temperatureField.WriteVTK_Field("test.vtk");
//	for (size_t i = 0; i < IDList.size(); i++)
//	{
//		int regionID = IDList[i];
//		FieldIO result = reader.GetFieldIO(regionID);
//		result.WriteTecplotField(std::to_string(regionID) + ".plt");
//	}
	return;
}



void EntranceOfCreateMapper(std::vector<std::string>& fileNames)
{
	if (3 != fileNames.size())
	{
		std::cout << "Please give 3 file names if you are intending to create solver, like this:" << std::endl;
		std::cout << "NCCouple createmapper (MOCMesh) (MOCField) (CFDMesh)" << std::endl;
		Logger::LogError("inccorrect number of file names");
	}
	CreateMapper(fileNames[0], fileNames[1], fileNames[2]);
	return;
}

void EntranceOfCFDToMOC(std::vector<std::string>& fileNames)
{
	if (0 != fileNames.size())
	{
		std::cout << "Please do not give any file name if you are intending to map CFD to MOC solver, like this:" << std::endl;
		std::cout << "NCCouple cfdtomoc" << std::endl;
		Logger::LogError("inccorrect number of file names");
	}
	CFDFieldsToMOC();
	return;
}

void EntranceOfMOCToCFD(std::vector<std::string>& fileNames)
{
	if (0 != fileNames.size())
	{
		std::cout << "Please do not give any file name if you are intending to map MOC to CFD solver, like this:" << std::endl;
		std::cout << "NCCouple moctocfd" << std::endl;
		Logger::LogError("inccorrect number of file names");
	}
	MOCFieldsToCFD();
	return;
}

void RunWithParameters(std::vector<std::string>& parameters)
{
	if (0 == parameters.size())
	{
		Logger::LogError("no parameter given in RunWithParameters()");
	}
	std::vector<std::string> filenames;
	filenames.resize(parameters.size() - 1);
	for (int i = 1;i < parameters.size();i++)
	{
		filenames[i - 1] = parameters[i];
	}
	if ("--help" == parameters[0])
	{
		if (0 != filenames.size())
		{
			Logger::LogError("no more parameter need to be given in --help");
		}
		DisplayHelpInfo();
	}
	else if ("clear" == parameters[0])
	{
		if (0 != filenames.size())
		{
			Logger::LogError("no more parameter need to be given in clear");
		}
		Logger::LogInfo("hello world");
	}
	else if ("createmapper" == parameters[0])
	{
		EntranceOfCreateMapper(filenames);
	}
	else if ("cfdtomoc" == parameters[0])
	{
		EntranceOfCFDToMOC(filenames);
	}
	else if ("moctocfd" == parameters[0])
	{
		EntranceOfMOCToCFD(filenames);
	}
	else
	{
		std::cout << "the following parameters are acceptabel: " << std::endl;
		std::cout << "(1) cfdtomoc" << std::endl;
		std::cout << "(2) clear" << std::endl;
		std::cout << "(3) createmapper" << std::endl;
		std::cout << "(4) moctocfd" << std::endl;
		std::cout << "(5) --help" << std::endl;
		Logger::LogError("invalid parameter given: " + parameters[0]);
	}
	return;
}

int main(int argc, char** argv)
{
	//get processor ID


	g_iProcessID = (int)getpid();
	if (argc == 1)
	{
		RenameFile("lingkong.txt", "lingkongxxx.txt");
		return 0;
	}
	else
	{
		std::vector<std::string> parameterList;
		parameterList.resize(argc - 1);
		for (int i = 1;i < argc;i++)
		{
			parameterList[i - 1] = argv[i];
		}
		RunWithParameters(parameterList);
	}

	VTKReaderTest();
	return 0;
}