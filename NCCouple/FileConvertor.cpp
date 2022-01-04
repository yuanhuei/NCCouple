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

//integration of a specified value over a specified region
//Note: it can be used both for MOC and CFD
Scalar Integration
(
	const GeneralMesh& mesh,
	ValueType vt,
	std::string strZoneName
)
{
	bool bSourceMoc = true;
	std::string sourceMeshName, targetMeshName;
	if (dynamic_cast<const CFDMesh*>(&mesh))
	{
		bSourceMoc = false;
	}
	Scalar integration = 0.0;
	for (int i = 0; i < mesh.GetMeshPointNum(); i++)
	{
		double sourceValue = 0;
		if (bSourceMoc)
		{
			const MOCMeshPoint& mocPoint = dynamic_cast<const MOCMeshPoint&>(*mesh.GetMeshPointPtr(i));
			if (mocPoint.GetMaterialName() != strZoneName) continue;
			sourceValue = mesh.GetMeshPointPtr(i)->GetValue(vt);
		}
		else
		{
			sourceValue = mesh.GetMeshPointPtr(i)->GetValue(vt);
		}
		double pointVolume = mesh.GetMeshPointPtr(i)->Volume();
		integration += sourceValue * pointVolume;
	}
	return integration;
}

void ConservationValidation
(
	const GeneralMesh& sourceMesh,
	const GeneralMesh& targetMesh,
	ValueType vt,
	std::string strZoneName
)
{
	std::string valueName = NameOfValueType(vt);
	double sourceIntegralValue = Integration(sourceMesh, vt, strZoneName);
	double targetIntegralValue = Integration(targetMesh, vt, strZoneName);
	Logger::LogInfo(FormatStr("Integral of %s on region %s of source mesh: %.6lf", valueName, strZoneName, sourceIntegralValue));
	Logger::LogInfo(FormatStr("Integral of %s on region %s of target mesh: %.6lf", valueName, strZoneName, targetIntegralValue));
	return;
}

void MOCFieldsToCFD
(
	std::string strInput_aplFileName,
	std::string strInput_inpFileName,
	std::string strInput_txtFileName,
	std::string strInput_meshFileName,
	std::string strOutput_vtkFileName,
	bool bRenew
)
{
	Logger::LogInfo("****** Mapping from MOC to CFD ******");
	Logger::LogInfo("Reading MOC files: " + strInput_aplFileName + ", " + strInput_inpFileName);
	MOCMesh mocMesh(strInput_aplFileName, strInput_inpFileName, MeshKernelType::MHT_KERNEL);
	Logger::LogInfo("Reading MOC heat power file: " + strInput_txtFileName);
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
	Logger::LogInfo("Reading CFD mesh file: " + strInput_meshFileName);
	//initialize with meshFile
	MHTVTKReader reader(strInput_meshFileName);
	for (size_t i = 0; i < reader.GetMeshListPtr().size(); i++)
	{
		Mesh* pmesh = reader.GetMeshListPtr()[i];
		Field<Scalar> heatpower(pmesh, 0.0, "heatpower");
		//read cfd mesh and create solver
		CFDMesh cfdMesh(pmesh, MeshKernelType::MHT_KERNEL, i);
		Solver solverMapper;
		std::string strMocMeshName = strInput_aplFileName.substr(0, strInput_aplFileName.find(".")) + "_";
		if (bRenew)
		{
			solverMapper = Solver(mocMesh, cfdMesh, pmesh->st_meshName, strMocMeshName);
		}
		else
		{
			solverMapper = Solver(mocMesh, cfdMesh, mocIndex, pmesh->st_meshName, strMocMeshName);
		}

		solverMapper.CheckMappingWeights();

		solverMapper.MOCtoCFDinterception(ValueType::HEATPOWER);

		cfdMesh.SetFieldValue(heatpower.elementField.v_value, ValueType::HEATPOWER);
		std::string strOutput_inpName = strOutput_vtkFileName.substr(0, strOutput_vtkFileName.find(".")) + pmesh->st_meshName+ ".vtk";
		Logger::LogInfo("Writing CFD vtk file: " + strOutput_inpName);
		heatpower.WriteVTK_Field(strOutput_inpName);
		ConservationValidation(mocMesh, cfdMesh, ValueType::HEATPOWER,pmesh->st_meshName);
	}
	return;
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
		Solver solverMapper;
		std::string strMocMeshName = strInput_aplFileName.substr(0, strInput_aplFileName.find(".")) + "_";
		if (bRenew)
			solverMapper = Solver(mocMesh, cfdMesh, pmesh->st_meshName, strMocMeshName);
		else
			solverMapper = Solver(mocMesh, cfdMesh, mocIndex, pmesh->st_meshName, strMocMeshName);

		solverMapper.CheckMappingWeights();

		solverMapper.CFDtoMOCinterception(ValueType::DENSITY);
		solverMapper.CFDtoMOCinterception(ValueType::TEMPERAURE);

		ConservationValidation(cfdMesh, mocMesh, ValueType::DENSITY,pmesh->st_meshName);
		ConservationValidation(cfdMesh, mocMesh, ValueType::TEMPERAURE,pmesh->st_meshName);
	}
	std::string strOutput_inpName = strInput_inpFileName.substr(0, strInput_inpFileName.find(".")) + "_out.inp";
	mocMesh.OutputStatus(strOutput_inpName);

}