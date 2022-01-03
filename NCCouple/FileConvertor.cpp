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
#include "./MHT_IO/VTKIO.h"


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

	MHTVTKReader reader("pinWR.msh");									//initialize with meshFile
	for (size_t i = 0; i < reader.GetMeshListPtr().size(); i++)
	{
		//reader.GetMeshListPtr()[i]->WriteTecplotMesh("pinWR_" + std::to_string(i) + ".plt");
		Mesh* pmesh = reader.GetMeshListPtr()[i];
		Field<Scalar> heatpower(pmesh, 0.0, "heatpower");
		//read cfd mesh and create solver
		CFDMesh cfdMesh(pmesh, MeshKernelType::MHT_KERNEL, i);
		Solver H2OMapper;
		if (bRenew)
			H2OMapper = Solver(mocMesh, cfdMesh, pmesh->st_meshName);
		else
			H2OMapper = Solver(mocMesh, cfdMesh, mocIndex, pmesh->st_meshName);

		H2OMapper.CheckMappingWeights();

		H2OMapper.MOCtoCFDinterception(ValueType::HEATPOWER);

		cfdMesh.SetFieldValue(heatpower.elementField.v_value, ValueType::HEATPOWER);
		std::string strOutput_inpName = strOutput_vtkFileName.substr(0, strOutput_vtkFileName.find(".")) + pmesh->st_meshName+ ".vtk";

		heatpower.WriteVTK_Field(strOutput_inpName);

		//create MHT field
		ConservationValidation(cfdMesh, mocMesh, ValueType::HEATPOWER,pmesh->st_meshName);
		ConservationValidation(mocMesh, cfdMesh, ValueType::HEATPOWER,pmesh->st_meshName);
	}

	//mocMesh.InitMOCValue(strInput_inpFileName);
	//create MHT mesh
	/*
	UnGridFactory meshFactoryCon(strInput_meshFileName, UnGridFactory::ugtFluent);
	FluentMeshBlock* FluentPtrCon = dynamic_cast<FluentMeshBlock*>(meshFactoryCon.GetPtr());
	RegionConnection Bridges;
	FluentPtrCon->Decompose(Bridges);
	Mesh* pmesh = &(FluentPtrCon->v_regionGrid[0]);
	*/
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

	std::vector<std::string> fieldName;
	fieldName.push_back("T");
	fieldName.push_back("Rho");
	MHTVTKReader reader(strInput_meshFileName,strInput_vtkFileName, fieldName);//initialize with meshFile
	for (size_t i = 0; i < reader.GetMeshListPtr().size(); i++)
	{
		//reader.GetMeshListPtr()[i]->WriteTecplotMesh("pinWR_" + std::to_string(i) + ".plt");
		Mesh* pmesh = reader.GetMeshListPtr()[i];
		//read cfd mesh and create solver
		CFDMesh cfdMesh(pmesh, MeshKernelType::MHT_KERNEL, i);
		for (int j = 0; j < reader.GetFieldList()[i].size(); j++)
		{
			const Field<Scalar>& field = reader.GetField(i,j);
			if(field.st_name=="T")
				cfdMesh.SetValueVec(field.elementField.v_value, ValueType::TEMPERAURE);
			else if(field.st_name=="Rho")
				cfdMesh.SetValueVec(field.elementField.v_value, ValueType::DENSITY);
		}
		Solver H2OMapper;
		if (bRenew)
			H2OMapper = Solver(mocMesh, cfdMesh, pmesh->st_meshName);
		else
			H2OMapper = Solver(mocMesh, cfdMesh, mocIndex, pmesh->st_meshName);

		H2OMapper.CheckMappingWeights();

		H2OMapper.CFDtoMOCinterception(ValueType::DENSITY);
		H2OMapper.CFDtoMOCinterception(ValueType::TEMPERAURE);

		ConservationValidation(cfdMesh, mocMesh, ValueType::DENSITY,pmesh->st_meshName);
		ConservationValidation(cfdMesh, mocMesh, ValueType::TEMPERAURE,pmesh->st_meshName);


	}
	std::string strOutput_inpName = strInput_inpFileName.substr(0, strInput_inpFileName.find(".")) + "_out.inp";
	mocMesh.OutputStatus(strOutput_inpName);

}