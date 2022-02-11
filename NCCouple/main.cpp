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



#include "./MHT_IO/VTKIO.h"
int g_iProcessID = 0;

//this example was designed for test of
//(1) rewritting a apl file
//(2) reading and writting of inp files

void WriteTotecplot(const MOCMesh& mocmesh,const CFDMesh& cfdmesh,
	const std::vector<int> v_iCFD,const std::vector<SMocIndex>& vSMocindex,std::string fileName)
{
	std::string  mType;
	std::ofstream ofile(fileName);
	ofile << "TITLE =\"" << "polyhedron" << "\"" << endl;
	ofile << "VARIABLES = " << "\"x\"," << "\"y\"," << "\"z\"" << endl;
	for (int i = 0; i < vSMocindex.size(); i++)
	{

		//平移坐标,
		int iAssembly = vSMocindex[i].iAssemblyIndex;
		double x = mocmesh.m_vAssembly[iAssembly].pAssembly_type->vAssemblyType_LeftDownPoint.x_- mocmesh.m_vAssembly[iAssembly].vAssembly_LeftDownPoint.x_;
		double y= mocmesh.m_vAssembly[iAssembly].pAssembly_type->vAssemblyType_LeftDownPoint.y_-mocmesh.m_vAssembly[iAssembly].vAssembly_LeftDownPoint.y_;
		
		const MHTMocMeshPoint& mocPoint = dynamic_cast<const MHTMocMeshPoint&>(*mocmesh.GetMocMeshPointPtr(vSMocindex[i]));

		MHTMocMeshPoint meshPoint = mocPoint;
		Vector vPoint(-x, -y, 0);
		meshPoint.Move(vPoint);
		meshPoint.WriteTecplotZones(ofile);
	}
	for (int i = 0; i < v_iCFD.size(); i++)
	{
		const MHTMeshPoint& mhtPolyhedron = dynamic_cast<const MHTMeshPoint&>(*cfdmesh.GetMeshPointPtr(v_iCFD[i]));
		mhtPolyhedron.WriteTecplotZones(ofile);
	}

	ofile.close();
	return;
}

void MOC_APL_INP_FileTest()
{
	WarningContinue("MOC_APL_INP_FileTest");
	MOCMesh mocMesh("c5g72l.apl", "c5g7.inp", MeshKernelType::MHT_KERNEL);
	//mocMesh.WriteSurfaceTecplotFile("filename.plt");
	//MOCMesh mocMesh("pin_c1.apl", "pin_c1.inp", MeshKernelType::MHT_KERNEL); 
	mocMesh.WriteTecplotFile("mMOD", "c5g721_cell12_mMOD.plt");
	mocMesh.WriteTecplotFile("mHE", "c5g721_cell12_mHE.plt");
	mocMesh.WriteTecplotFile("mCD", "c5g721_cell12_mCD.plt");
	mocMesh.WriteTecplotFile("mUO2", "c5g721_cell12_mUO2.plt");
	//mocMesh.WriteTecplotFile("Zr4", "zr4.plt");
	//mocMesh.WriteTecplotFile("UO2", "u2o.plt");
	//mocMesh.InitMOCValue("pin_c1.inp","pin_c1.txt");
	//mocMesh.OutputStatus("pin_c1_out.inp");
	return;
}
void MapTest()
{
	//RegisterMapper("c5g72l.apl", "c5g72l_out.apl", "PIN9Coarse.msh");
	//MOCFieldsToCFD();

	std::string configFile = "MapFile_FileNames";
	std::vector<std::vector<std::string> > matches = GetMatchList(configFile);
	std::vector<std::string>& materialList = matches[0];
	std::vector<std::string>& regionList = matches[1];
	std::string mocMeshFile = GetFileName(configFile, "inputApl");
	std::string outMocMeshFile = GetFileName(configFile, "outputApl");
	std::string cfdMeshFile = GetFileName(configFile, "inputMsh");

	MOCMesh mocMesh(mocMeshFile, outMocMeshFile);

	MHTVTKReader reader(cfdMeshFile);
	for (size_t i = 0; i < regionList.size(); i++)
	{
		if (i > 0)break;
		int CFDMeshID = reader.GetIDOfRegion(regionList[i]);
		Mesh* pmesh = reader.GetMeshListPtr()[CFDMeshID];
		//read cfd mesh and create solver
		CFDMesh cfdMesh(pmesh, MeshKernelType::MHT_KERNEL, CFDMeshID);
		/*
		std::vector<SMocIndex> thh;
		for (int j = 0; j < 64; j++)
		{
			thh.push_back(SMocIndex(0, 51, j));
		}
		WriteTotecplot(mocMesh, cfdMesh, std::vector<int>{2559}, thh, "temp.plt");
		*/
		//std::string outfilename = "mesh_" + std::to_string(i) + ".plt";
		//cfdMesh.WriteTecplotFile(outfilename);
		Solver solverMapper(mocMesh, cfdMesh, materialList[i], true);
		solverMapper.CheckMappingWeights();
	}


	//InitCFDMeshValue(cfdMesh);
	//solver.CFDtoMOCinterception(ValueType::DENSITY);

	//mocMesh.OutputStatus("pin_c1.inp");

}
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

Scalar initialHeatPower(Scalar x, Scalar y, Scalar z)
{
	Vector PointCenter(x, y, z);
	Vector axisCenter(0.63, 0.63, 0.0);
	Vector axisNorm(0.0, 0.0, 1.0);
	Vector OP = PointCenter - axisCenter;
	Vector axicialLocation = (OP & axisNorm) * axisNorm;
	Scalar radius = (OP - axicialLocation).Mag();
	Scalar sigma = 1.0;
	Scalar hp = (1000.0 / sqrt(2.0 * PI) / sigma) * exp(-radius * radius / 2 / pow(sigma, 2));
	hp = hp * z * (5.0 - z) / 6.25;
	return hp;
}

#include "./MHT_IO/VTKIO.h"
#include <fstream>

void VTKReaderTest()
{
	MHTVTKReader reader("pinWR.msh");
	std::vector<std::string> fileName;
	fileName.push_back("H2O.vtk");
	fileName.push_back("Zr4.vtk");
	std::vector<int> IDList;
	IDList.push_back(0);
	IDList.push_back(1);
	std::vector<std::string> fieldName;
	fieldName.push_back("heatpower");
	reader.ReadVTKFile(fileName, IDList, fieldName);
	reader.GetFieldIO(0).WriteTecplotField("heatpower_0.plt");
	reader.GetFieldIO(1).WriteTecplotField("heatpower_1.plt");
	return;
}

void EntranceOfRegister(std::vector<std::string>& fileNames)
{
	if (3 != fileNames.size())
	{
		std::cout << "Please give 3 file names if you are intending to register solver, like this:" << std::endl;
		std::cout << "NCCouple register (MOCMesh) (MOCField) (CFDMesh)" << std::endl;
		Logger::LogError("inccorrect number of file names");
	}
	RegisterMapper(fileNames[0], fileNames[1], fileNames[2]);
	return;
}

void EntranceOfCreateMapper(std::vector<std::string>& fileNames)
{
	if (0 != fileNames.size())
	{
		std::cout << "Please do not give any file name if you are intending to create mapper, like this:" << std::endl;
		std::cout << "NCCouple createmapper" << std::endl;
		Logger::LogError("inccorrect number of file names");
	}
	CreateMapper();
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
		ClearMapFiles();
	}
	else if ("register" == parameters[0])
	{
		EntranceOfRegister(filenames);
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
		std::cout << "(5) register" << std::endl;
		std::cout << "(6) --help" << std::endl;
		Logger::LogError("invalid parameter given: " + parameters[0]);
	}
	return;
}

int main(int argc, char** argv)
{
	//get processor ID
	//g_iProcessID = (int)getpid();
	if (argc == 1)
	{
		//PolyhedronSet box(Vector(0, 0, 0), Vector(1, 1, 1));
		//box.MHT::Polyhedron::Display();
		//MOC_APL_INP_FileTest();
		//CFDFieldsToMOC();
		MapTest();
		//MOCMesh mocmesh = MOCMesh();
		//mocmesh.InitMOCFromInputFile("c5g7.inp");
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
	return 0;
}