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
#include "./MHT_mesh/UnGridFactory.h"
#include "./MHT_mesh/RegionConnection.h"
#include "./MHT_mesh/Mesh.h"
#include "./MHT_field/Field.h"
int g_iProcessID = 0;
enum class Material
{
	H2O,
	Zr4,
	UO2
};
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

void ConservationValidation(const GeneralMesh& sourceMesh, const GeneralMesh& targetMesh, ValueType vt) {
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

void MOCCFDMapping()
{
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
	//read cfd mesh and create solver
	std::string CFDMeshFile = "pinW.msh";
	CFDMesh H2OcfdMesh(CFDMeshFile, MeshKernelType::MHT_KERNEL, int(Material::H2O));
	//CFDMesh H2OcfdMesh("outcfdtemp_cfd", MeshKernelType::MHT_KERNEL);
	Solver H2OMapper(mocMesh, H2OcfdMesh, mocIndex, "H2O");
	H2OMapper.CheckMappingWeights();
	//read cfd mesh and create solver
	CFDMesh Zr4cfdMesh(CFDMeshFile, MeshKernelType::MHT_KERNEL, int(Material::Zr4));
	Solver Zr4Mapper(mocMesh, Zr4cfdMesh, mocIndex, "Zr4");
	Zr4Mapper.CheckMappingWeights();
	//read cfd mesh and create solver
	CFDMesh U2OcfdMesh(CFDMeshFile, MeshKernelType::MHT_KERNEL, int(Material::UO2));
	Solver U2OMapper(mocMesh, U2OcfdMesh, mocIndex, "UO2");
	U2OMapper.CheckMappingWeights();
	end = time(NULL);
	Logger::LogInfo(FormatStr("Time for caculatation: %d second", int(difftime(end, start))));
	//initialization of a scalar field on CFD mesh at H2O region
	InitCFDMeshValue(H2OcfdMesh);
	//Mapping of the scalar field from CFD to MOC
	H2OMapper.CFDtoMOCinterception(ValueType::DENSITY);
	//writting input file
	mocMesh.OutputStatus("pin_c1.inp");
	ConservationValidation(H2OcfdMesh, mocMesh, ValueType::DENSITY);
	ConservationValidation(mocMesh, H2OcfdMesh, ValueType::DENSITY);
	return;
}

void ReadCFDMesh()
{
	UnGridFactory meshFactoryCon("pinW.msh", UnGridFactory::ugtFluent);
	FluentMeshBlock* FluentPtrCon = dynamic_cast<FluentMeshBlock*>(meshFactoryCon.GetPtr());
	RegionConnection Bridges;
	FluentPtrCon->Decompose(Bridges);
	Mesh* pmesh = &(FluentPtrCon->v_regionGrid[0]);

	Field<Scalar> T(pmesh, 300.0, "T");
	T.WriteTecplotField("T.plt");
	return;
}

int main()
{
	MOCCFDMapping();
	//ReadCFDMesh();
	return 0;
}