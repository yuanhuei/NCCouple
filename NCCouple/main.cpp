#include <iostream>
#include <algorithm>
#include "MOCMesh.h"
#include "MOCIndex.h"
#include "CFDMesh.h"
#include "Solver.h"
#include "Logger.h"
#include "MHT_polyhedron/PolyhedronSet.h"
#include<time.h>
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
enum class Material
{
	H2O,
	Zr4,
	UO2
};

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


void InitCFDMeshValue(const GeneralMesh& cfdMesh)
{
	for (int i = 0; i < cfdMesh.GetMeshPointNum(); i++)
	{
		double x, y, z;
		std::tie(x, y, z) = cfdMesh.GetMeshPointPtr(i)->CentralCoordinate();
		double _x = x - 0.63;
		double _y = y - 0.63;
		double r = sqrt(_x * _x + _y * _y);
		double value = r + z + _x / r;
		cfdMesh.GetMeshPointPtr(i)->SetValue(value, ValueType::DENSITY);
	}
	return;
}

void ConservationValidation
(
	const GeneralMesh& sourceMesh,
	const GeneralMesh& targetMesh,
	ValueType vt
)
{
	double sourceIntegralValue = 0.0;
	double targetIntegralValue = 0.0;

	for (int i = 0; i < sourceMesh.GetMeshPointNum(); i++) {
		double sourceValue = sourceMesh.GetMeshPointPtr(i)->GetValue(vt);
		double pointVolume = sourceMesh.GetMeshPointPtr(i)->Volume();
		sourceIntegralValue += sourceValue * pointVolume;
	}

	for (int i = 0; i < targetMesh.GetMeshPointNum(); i++) {
		double targetValue = targetMesh.GetMeshPointPtr(i)->GetValue(vt);
		double pointVolume = targetMesh.GetMeshPointPtr(i)->Volume();
		targetIntegralValue += targetValue * pointVolume;
	}

	std::string sourceMeshName, targetMeshName;
	if (dynamic_cast<const CFDMesh*>(&sourceMesh)) {
		sourceMeshName = "CFD Mesh";
		targetMeshName = "MOC Mesh";
	}
	else {
		sourceMeshName = "MOC Mesh";
		targetMeshName = "CFD Mesh";
	}
	Logger::LogInfo(FormatStr("The Integral Value of Source %s : %.6lf", sourceMeshName.c_str(), sourceIntegralValue));
	Logger::LogInfo(FormatStr("The Integral Value of Target %s : %.6lf", targetMeshName.c_str(), targetIntegralValue));
	return;
}

//this example was designed for test of
//(1) reading a CFD mesh file
//(2) creating a field and read variables from a vkt file
void ReadCFDMeshAndFieldTest()
{
	WarningContinue("ReadCFDMeshAndFieldTest");
	//open a mesh file
	UnGridFactory meshFactoryCon("pinW.msh", UnGridFactory::ugtFluent);
	FluentMeshBlock* FluentPtrCon = dynamic_cast<FluentMeshBlock*>(meshFactoryCon.GetPtr());
	RegionConnection Bridges;
	FluentPtrCon->Decompose(Bridges);
	Mesh* pmesh = &(FluentPtrCon->v_regionGrid[0]);
	//create field and read from vtk file
	Field<Scalar> T(pmesh, 0.0, "T");
	Field<Scalar> rho(pmesh, 0.0, "Rho");
	T.ReadVTK_Field("pinW.vtk");
	rho.ReadVTK_Field("pinW.vtk");
	T.WriteTecplotField("T.plt");
	rho.WriteTecplotField("rho.plt");
	return;
}

//this example was designed for test of
//(1) creating mapping solvers for different regions
//(2) conservation validation of H2O region
void SolverCreatingTest()
{
	WarningContinue("SolverCreatingTest");
	//get processor ID
	g_iProcessID = (int)getpid();
	MOCMesh mocMesh("pin_c1.apl", MeshKernelType::MHT_KERNEL);
	//examples for writing tecplot files of each materials
	//Note: these file can be viewed by Tecplot
	mocMesh.WriteTecplotFile("H2O", "H2OMOCFile.plt");
	mocMesh.WriteTecplotFile("Zr4", "Zr4MOCFile.plt");
	mocMesh.WriteTecplotFile("UO2", "U2OMOCFile.plt");
	//create an index for fast searching
	MOCIndex mocIndex(mocMesh);
	//the following information should be given for a specified tube
	mocIndex.axisNorm = Vector(0.0, 0.0, 1.0);
	mocIndex.axisPoint = Vector(0.63, 0.63, 0.0);
	mocIndex.theetaStartNorm = Vector(1.0, 0.0, 0.0);
	mocIndex.circularCellNum = 8;
	mocIndex.axialCellNum = 5;
	mocIndex.axialCellSize = 1.0;
	std::vector<Scalar> radiusList;
	radiusList.push_back(0.1024);
	radiusList.push_back(0.2048);
	radiusList.push_back(0.3072);
	radiusList.push_back(0.4096);
	radiusList.push_back(0.475);
	mocIndex.SetRadial(radiusList);
	mocIndex.BuildUpIndex();
	//mapper solvers are created for each zone (H2O, Zr4 and U2O) 
	time_t start, end;
	start = time(NULL);
	//read cfd mesh
	UnGridFactory meshFactoryCon("pinW.msh", UnGridFactory::ugtFluent);
	FluentMeshBlock* FluentPtrCon = dynamic_cast<FluentMeshBlock*>(meshFactoryCon.GetPtr());
	RegionConnection Bridges;
	FluentPtrCon->Decompose(Bridges);
	Mesh* pmesh0 = &(FluentPtrCon->v_regionGrid[0]);
	Mesh* pmesh1 = &(FluentPtrCon->v_regionGrid[1]);
	Mesh* pmesh2 = &(FluentPtrCon->v_regionGrid[2]);
	//read cfd mesh and create solver
	CFDMesh H2OcfdMesh(pmesh0, MeshKernelType::MHT_KERNEL, int(Material::H2O));
	Solver H2OMapper(mocMesh, H2OcfdMesh, mocIndex, "H2O");
	H2OMapper.CheckMappingWeights();
	//read cfd mesh and create solver
	CFDMesh Zr4cfdMesh(pmesh1, MeshKernelType::MHT_KERNEL, int(Material::Zr4));
	Solver Zr4Mapper(mocMesh, Zr4cfdMesh, mocIndex, "Zr4");
	Zr4Mapper.CheckMappingWeights();
	//read cfd mesh and create solver
	CFDMesh U2OcfdMesh(pmesh2, MeshKernelType::MHT_KERNEL, int(Material::UO2));
	Solver U2OMapper(mocMesh, U2OcfdMesh, mocIndex, "UO2");
	U2OMapper.CheckMappingWeights();
	end = time(NULL);
	Logger::LogInfo(FormatStr("Time for caculatation: %d second", int(difftime(end, start))));
	//Conservation Validation in H2O region
	ConservationValidation(H2OcfdMesh, mocMesh, ValueType::DENSITY);
	ConservationValidation(mocMesh, H2OcfdMesh, ValueType::DENSITY);
	return;
}

//this example was designed for test of
//(1) creating a solver
//(2) reading field from a vtk file
//(3) interpolations of CFD to MOC and then MOC to CFD
void SolverCreatingAndMappingTest()
{
	WarningContinue("SolverCreatingAndMappingTest");
	//get processor ID
	g_iProcessID = (int)getpid();
	MOCMesh mocMesh("pin_c1.apl", MeshKernelType::MHT_KERNEL);
	//create an index for fast searching
	MOCIndex mocIndex(mocMesh);
	//the following information should be given for a specified tube
	mocIndex.axisNorm = Vector(0.0, 0.0, 1.0);
	mocIndex.axisPoint = Vector(0.63, 0.63, 0.0);
	mocIndex.theetaStartNorm = Vector(1.0, 0.0, 0.0);
	mocIndex.circularCellNum = 8;
	mocIndex.axialCellNum = 5;
	mocIndex.axialCellSize = 1.0;
	std::vector<Scalar> radiusList;
	radiusList.push_back(0.1024);
	radiusList.push_back(0.2048);
	radiusList.push_back(0.3072);
	radiusList.push_back(0.4096);
	radiusList.push_back(0.475);
	mocIndex.SetRadial(radiusList);
	mocIndex.BuildUpIndex();
	//create MHT mesh
	UnGridFactory meshFactoryCon("pinW.msh", UnGridFactory::ugtFluent);
	FluentMeshBlock* FluentPtrCon = dynamic_cast<FluentMeshBlock*>(meshFactoryCon.GetPtr());
	RegionConnection Bridges;
	FluentPtrCon->Decompose(Bridges);
	Mesh* pmesh = &(FluentPtrCon->v_regionGrid[0]);
	//create MHT field
	Field<Scalar> rho(pmesh, 0.0, "Rho");
	//read cfd mesh and create solver
	CFDMesh H2OcfdMesh(pmesh, MeshKernelType::MHT_KERNEL, int(Material::H2O));
	Solver H2OMapper(mocMesh, H2OcfdMesh, mocIndex, "H2O");
	H2OMapper.CheckMappingWeights();
	//read CFD Field From Field
	rho.ReadVTK_Field("pinW.vtk");
	H2OcfdMesh.SetValueVec(rho.elementField.v_value, ValueType::DENSITY);
	H2OMapper.CFDtoMOCinterception(ValueType::DENSITY);
	H2OMapper.MOCtoCFDinterception(ValueType::DENSITY);
	H2OcfdMesh.SetFieldValue(rho, ValueType::DENSITY);
	rho.WriteTecplotField("rho.plt");
	return;
}

//this example was designed for test of
//(1) rewritting a apl file
//(2) reading and writting of inp files
void MOC_APL_INP_FileTest()
{
	WarningContinue("MOC_APL_INP_FileTest");
	MOCMesh mocMesh("pin_c1.apl", MeshKernelType::MHT_KERNEL);
	mocMesh.InitMOCValue("pin_c1.inp");
	mocMesh.OutputStatus("pin_c1_out.inp");
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


	fieldName.push_back("T");
	fieldName.push_back("Rho");


	MHTVTKReader reader("pinWR.msh", "Field_Data.txt", fieldName);
	for (size_t i = 0; i < reader.GetFieldIOList().size(); i++)
	{
		reader.GetFieldIO(i).WriteTecplotField("pinWR_"+std::to_string(i)+".plt");

	}
	
	std::cout <<"success end"<<std::endl;
}
int main()
{
	//	ReadCFDMeshAndFieldTest();
	//	SolverCreatingTest();
	//	SolverCreatingAndMappingTest();
	//	MOC_APL_INP_FileTest();
	VTKReaderTest();
	return 0;


}