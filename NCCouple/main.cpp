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

#include <vtkDelaunay3D.h>
#include <vtkNew.h>
#include <vtkSphereSource.h>
#include <vtkXMLPUnstructuredGridWriter.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkAppendFilter.h>
#include <vtkMultiBlockDataSet.h>
void VTK_Test()
{
	vtkObject::GlobalWarningDisplayOff();
	vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
	reader->SetFileName("500_0.vtk");
	reader->SetReadAllColorScalars(true);
	reader->SetReadAllFields(true);
	reader->SetReadAllScalars(true);
	reader->Update();

	vtkSmartPointer<vtkUnstructuredGrid> Grid;
	Grid = reader->GetOutput();

	vtkSmartPointer<vtkUnstructuredGridReader> reader_1 = vtkSmartPointer<vtkUnstructuredGridReader>::New();
	reader_1->SetFileName("500_1.vtk");
	reader_1->SetReadAllColorScalars(true);
	reader_1->SetReadAllFields(true);
	reader_1->SetReadAllScalars(true);
	reader_1->Update();

	vtkSmartPointer<vtkUnstructuredGrid> Grid_1;
	Grid_1 = reader_1->GetOutput();

	// create the append filter
	vtkSmartPointer<vtkAppendFilter> append =
		vtkSmartPointer<vtkAppendFilter>::New();

	// add each data set
	append->AddInputData(Grid);
	append->AddInputData(Grid_1);
	append->Update();


	vtkSmartPointer<vtkUnstructuredGridWriter> writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
	writer->SetInputData(append->GetOutput());

	writer->SetFileName("500_0_1.vtk");
	writer->Update();
}

#include "./MHT_IO/VTKIO.h"
#include <fstream>
void VTKReaderTest()
{

	std::vector<std::string> fieldName;


	fieldName.push_back("Pressure");
	fieldName.push_back("temperature");

	std::vector<std::string> vVTKName;
	vVTKName.push_back("pipe_solid.vtk");
	vVTKName.push_back("pipe_fluid.vtk");


	MHTVTKReader reader("pipe.msh");		//initialize with meshFile and vVTKFile
	reader.ReadVTKFile(vVTKName, fieldName);
//	MHTVTKReader reader("pinWR.msh", fieldName);						//initialize with meshFile and EmptyField
//	MHTVTKReader reader("pinWR.msh");									//initialize with meshFile

	for (size_t i = 0; i < reader.GetFieldIOList().size(); i++)
	{
		reader.GetFieldIO(i).WriteTecplotField("pipe_" + std::to_string(i) + ".plt");
	}
	reader.WriteDataFile("TestData.txt");


//	for (size_t i = 0; i < reader.GetMeshListPtr().size(); i++)
//	{
//		reader.GetMeshListPtr()[i]->WriteTecplotMesh("pinWR_"+std::to_string(i)+".plt");
//
//	}
	
	std::cout <<"success end"<<std::endl;
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
	if (4 != fileNames.size())
	{
		std::cout << "Please give 4 file names if you are intending to map CFD to MOC solver, like this:" << std::endl;
		std::cout << "NCCouple cfdtomoc (CFDMesh) (CFDField) (MOCMesh) (MOCField)" << std::endl;
		Logger::LogError("inccorrect number of file names");
	}
	CFDFieldsToMOC(fileNames[0], fileNames[1], fileNames[2], fileNames[3]);
	return;
}

void EntranceOfMOCToCFD(std::vector<std::string>& fileNames)
{
	if (5 != fileNames.size())
	{
		std::cout << "Please give 5 file names if you are intending to map MOC to CFD solver, like this:" << std::endl;
		std::cout << "NCCouple moctocfd (MOCMesh) (MOCField) (MOCPower) (CFDMesh) (CFDField)" << std::endl;
		Logger::LogError("inccorrect number of file names");
	}
	MOCFieldsToCFD(fileNames[0], fileNames[1], fileNames[2], fileNames[3], fileNames[4]);
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
		std::cout << "hello world" << std::endl;
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

//	VTKReaderTest();
	return 0;
}