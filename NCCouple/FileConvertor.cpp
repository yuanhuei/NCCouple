#include <iostream>
#include <algorithm>
#include "MOCMesh.h"
#include "MOCIndex.h"
#include "CFDMesh.h"
#include "Solver.h"
#include "Logger.h"
#include "MHT_polyhedron/PolyhedronSet.h"
#include "FileConvertor.h"

#include "./MHT_common/SystemControl.h"
#include "./MHT_mesh/UnGridFactory.h"
#include "./MHT_mesh/RegionConnection.h"
#include "./MHT_mesh/Mesh.h"
#include "./MHT_field/Field.h"
#include "./MHT_IO/FieldIO.h"



void MOCFieldsToCFD
(
	std::string  strInput_aplFileName,
	std::string strInput_inpFileName,
	std::string strInput_txtFileName,
	std::string strInput_meshFileName,
	std::string strOutput_vtkFileName,
	bool bRenew
)
{
	WarningContinue("SolverCreatingAndMappingTest");
	MOCMesh mocMesh(strInput_aplFileName, strInput_inpFileName, MeshKernelType::MHT_KERNEL);
	mocMesh.InitMOCHeatPower(strInput_txtFileName);

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

	//mocMesh.InitMOCValue(strInput_inpFileName);
	//create MHT mesh
	UnGridFactory meshFactoryCon(strInput_meshFileName, UnGridFactory::ugtFluent);
	FluentMeshBlock* FluentPtrCon = dynamic_cast<FluentMeshBlock*>(meshFactoryCon.GetPtr());
	RegionConnection Bridges;
	FluentPtrCon->Decompose(Bridges);
	Mesh* pmesh = &(FluentPtrCon->v_regionGrid[0]);
	//create MHT field
	Field<Scalar> heatpower (pmesh, 0.0, "heatpower");
	//read cfd mesh and create solver
	CFDMesh H2OcfdMesh(pmesh, MeshKernelType::MHT_KERNEL, int(Material::H2O));
	Solver H2OMapper;
	if(bRenew)
		H2OMapper=Solver(mocMesh, H2OcfdMesh, "H2O"); 
	else
		H2OMapper=Solver(mocMesh, H2OcfdMesh, mocIndex, "H2O");

	H2OMapper.CheckMappingWeights();

	H2OMapper.MOCtoCFDinterception(ValueType::HEATPOWER);

	H2OcfdMesh.SetFieldValue(heatpower.elementField.v_value, ValueType::HEATPOWER);
	heatpower.WriteVTK_Field(strOutput_vtkFileName);

	ConservationValidation(H2OcfdMesh, mocMesh, ValueType::HEATPOWER);
	ConservationValidation(mocMesh, H2OcfdMesh, ValueType::HEATPOWER);
}

void CFDFieldsToMOC
(
	std::string strInput_meshFileName,
	std::string strInput_vtkFileName,
	std::string  strInput_aplFileName,
	std::string strInput_inpFileName,
	bool bRenew//判断是否是第一次做插值，false为第一次
)
{
	WarningContinue("SolverCreatingAndMappingTest");

	MOCMesh mocMesh(strInput_aplFileName, strInput_inpFileName,MeshKernelType::MHT_KERNEL);
	//mocMesh.InitMOCValue("pin_c1.inp","pin_c1.txt");
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
	UnGridFactory meshFactoryCon(strInput_meshFileName, UnGridFactory::ugtFluent);
	FluentMeshBlock* FluentPtrCon = dynamic_cast<FluentMeshBlock*>(meshFactoryCon.GetPtr());
	RegionConnection Bridges;
	FluentPtrCon->Decompose(Bridges);
	Mesh* pmesh = &(FluentPtrCon->v_regionGrid[0]);
	//create MHT field
	Field<Scalar> rho(pmesh, 0.0, "Rho");
	rho.ReadVTK_Field(strInput_vtkFileName);
	Field<Scalar> T(pmesh, 0.0, "T");
	T.ReadVTK_Field(strInput_vtkFileName);

	//read cfd mesh and create solver
	CFDMesh H2OcfdMesh(pmesh, MeshKernelType::MHT_KERNEL, int(Material::H2O));
	H2OcfdMesh.SetValueVec(rho.elementField.v_value, ValueType::DENSITY);
	H2OcfdMesh.SetValueVec(T.elementField.v_value, ValueType::TEMPERAURE);
	Solver H2OMapper;
	if (bRenew)
		H2OMapper = Solver(mocMesh, H2OcfdMesh, "H2O");
	else
		H2OMapper = Solver(mocMesh, H2OcfdMesh, mocIndex, "H2O");

	H2OMapper.CheckMappingWeights();

	H2OMapper.CFDtoMOCinterception(ValueType::DENSITY);
	H2OMapper.CFDtoMOCinterception(ValueType::TEMPERAURE);

	ConservationValidation(H2OcfdMesh, mocMesh, ValueType::DENSITY);
	ConservationValidation(mocMesh, H2OcfdMesh, ValueType::DENSITY);
	ConservationValidation(H2OcfdMesh, mocMesh, ValueType::TEMPERAURE);
	ConservationValidation(mocMesh, H2OcfdMesh, ValueType::TEMPERAURE);

	std::string strOutput_inpName = strInput_inpFileName.substr(0, strInput_inpFileName.find(".")) + "_out.inp";
	mocMesh.OutputStatus(strOutput_inpName);

	
}