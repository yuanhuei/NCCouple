#include <iostream>
#include <algorithm>
#include "MOCMesh.h"
#include "MOCIndex.h"
#include "CFDMesh.h"
#include "Solver.h"
#include "Logger.h"
#include "./MHT_polyhedron/PolyhedronSet.h"
#include "./FileConvertor.h"
#include "./ConfigurationFile.h"
#include "./MHT_common/SystemControl.h"
#include "./MHT_mesh/UnGridFactory.h"
#include "./MHT_mesh/RegionConnection.h"
#include "./MHT_mesh/Mesh.h"
#include "./MHT_field/Field.h"
#include "./MHT_IO/FieldIO.h"
#include "./MHT_IO/VTKIO.h"

std::string configFile = "MapFile_FileNames";

//integration of a specified value over a specified region
//Note: it can be used both for MOC and CFD
Scalar Integration
(
	const GeneralMesh& mesh,
	ValueType vt,
	std::string strZoneName,
	Solver& solverMapper
)
{
	bool bSourceMoc = true;
	std::string sourceMeshName, targetMeshName;
	if (dynamic_cast<const CFDMesh*>(&mesh))
	{
		bSourceMoc = false;
	}
	Scalar integration = 0.0;
	double sourceValue = 0, pointVolume = 0;
	if (bSourceMoc)
	{
		 const MOCMesh& mocMesh = dynamic_cast<const MOCMesh&> (mesh);
		 std::vector< SMocIndex> vSMocIndex;
		 solverMapper.GetMocIndexByMapValue(vSMocIndex);

		for (int i = 0; i < vSMocIndex.size(); i++)
		{
			//const MHTMocMeshPoint& mocPoint = dynamic_cast<const MHTMocMeshPoint&>(*mocMesh.GetMocMeshPointPtr(vSMocIndex[i]));
			if (mocMesh.GetMaterialNameAtIndex(vSMocIndex[i]) != strZoneName) continue;
			sourceValue = mocMesh.GetValueAtIndex(vSMocIndex[i],vt);
			pointVolume = mocMesh.GetVolumeAtIndex(vSMocIndex[i]);
			integration += sourceValue * pointVolume;
		}
	}
	else
	{
		for (int i = 0; i < mesh.GetMeshPointNum(); i++)
		{
			sourceValue = mesh.GetMeshPointPtr(i)->GetValue(vt);
			pointVolume = mesh.GetMeshPointPtr(i)->Volume();
			integration += sourceValue * pointVolume;
		}
	}
	return integration;
}

void ConservationValidation
(
	const GeneralMesh& sourceMesh,
	const GeneralMesh& targetMesh,
	ValueType vt,
	std::string strZoneName,
	Solver& solverMapper
)
{
	std::string valueName = NameOfValueType(vt);
	double sourceIntegralValue = Integration(sourceMesh, vt, strZoneName, solverMapper);
	double targetIntegralValue = Integration(targetMesh, vt, strZoneName, solverMapper);
	Logger::LogInfo(FormatStr("Integral of %s on region %s of source mesh: %.6lf", valueName.c_str(), strZoneName.c_str(), sourceIntegralValue));
	Logger::LogInfo(FormatStr("Integral of %s on region %s of target mesh: %.6lf", valueName.c_str(), strZoneName.c_str(), targetIntegralValue));
	return;
}

void DisplayHelpInfo()
{
	std::cout << "Please input command like one list below:" << std::endl;
	std::cout << "1. NCCouple" << std::endl;
	std::cout << "2. NCCouple --help" << std::endl;
	std::cout << "3. NCCouple clear" << std::endl;
	std::cout << "4. NCCouple createmapper" << std::endl;
	std::cout << "5. NCCouple cfdtomoc" << std::endl;
	std::cout << "6. NCCouple moctocfd" << std::endl;
	std::cout << "7. NCCouple register (MOCMesh) (MOCMesh) (CFDMesh)" << std::endl;
	return;
}

void InsertWhenNotFound(std::vector<std::string>& nameList, std::string name)
{
	bool found = false;
	for (int j = 0;j < nameList.size();j++)
	{
		if (name == nameList[j])
		{
			found = true;
			break;
		}
	}
	if (false == found)
	{
		nameList.push_back(name);
	}
	return;
}

void RegisterMapper
(
	std::string strInput_aplFileName,
	std::string strOutput_aplFileName
)
{
	//checking file names
	if (strInput_aplFileName.find(".apl") == std::string::npos)
	{
		Logger::LogError("in CreateMapper, the input MOC mesh " + strInput_aplFileName + " is not an .apl file");
	}
	if (strOutput_aplFileName.find(".apl") == std::string::npos)
	{
		Logger::LogError("in CreateMapper, the output MOC mesh " + strOutput_aplFileName + " is not an .apl file");
	}
	if (strOutput_aplFileName == strInput_aplFileName)
	{
		Logger::LogError("in CreateMapper, the output MOC mesh " + strOutput_aplFileName + " is the same with the input one " + strInput_aplFileName);
	}
	MOCMesh mocMesh(strInput_aplFileName, strOutput_aplFileName, MeshKernelType::MHT_KERNEL);
	std::vector<std::string> MOCRegionList;
	std::vector<std::string> CFDRegionList;
	//collect region names from MOC mesh
	for (size_t i = 0;i < mocMesh.m_vSMocIndex.size(); i++)
	{
		const MOCMeshPoint& mocPoint = dynamic_cast<const MOCMeshPoint&>(*mocMesh.GetMocMeshPointPtr(mocMesh.m_vSMocIndex[i]));
		std::string thisName = mocPoint.GetMaterialName();
		InsertWhenNotFound(MOCRegionList, thisName);
	}
	WriteConfigurationFile(configFile, strInput_aplFileName, strOutput_aplFileName, MOCRegionList);
	WarningContinue("please fill the blanks in file " + configFile + " before createmapper");
	return;
}

void CreateMapper()
{
	std::vector<std::vector<std::string> > matches = GetMatchList(configFile);
	std::vector<std::string>& materialList = matches[0];
	std::vector<std::string>& inputVTKList = matches[1];
	std::string mocMeshFile = GetFileName(configFile, "inputApl");
	std::string outMocMeshFile = GetFileName(configFile, "outputApl");
	//checking file names
	if (mocMeshFile.find(".apl") == std::string::npos)
	{
		Logger::LogError("in MOCFieldsToCFD, " + mocMeshFile + " is not an .apl file");
	}
	if (outMocMeshFile.find(".apl") == std::string::npos)
	{
		Logger::LogError("in MOCFieldsToCFD, " + outMocMeshFile + " is not an .apl file");
	}
	for (int i = 0;i < inputVTKList.size();i++)
	{
		if (inputVTKList[i].find(".vtk") == std::string::npos)
		{
			Logger::LogError("in MOCFieldsToCFD, " + inputVTKList[i] + " is not a .vtk file");
		}
	}
	MOCMesh mocMesh(mocMeshFile, outMocMeshFile, MeshKernelType::MHT_KERNEL);
	//create an index for fast searching
	Logger::LogInfo("Reading CFD mesh from input vtk files");
	Scalar ratio = GetValueInFile(configFile, "scaleRatio");
	MHTVTKReader reader(inputVTKList, ratio);
	for (size_t i = 0; i < inputVTKList.size(); i++)
	{
		Mesh* pmesh = reader.GetMeshListPtr()[i];
		//read cfd mesh and create solver
		CFDMesh cfdMesh(pmesh, MeshKernelType::MHT_KERNEL, i);
		Solver solverMapper(mocMesh, cfdMesh, materialList[i],true);
		solverMapper.CheckMappingWeights();
	}
	return;
}

void MOCFieldsToCFD()
{
	std::string mocMeshFile = GetFileName(configFile, "inputApl");
	std::string outMocMeshFile = GetFileName(configFile, "outputApl");
	std::string mocPowerFile = GetFileName(configFile, "mocPower");
	std::vector<std::vector<std::string> > matches = GetMatchList(configFile);
	std::vector<std::string>& materialList = matches[0];
	std::vector<std::string>& inputVTKList = matches[1];
	std::vector<std::string>& outputVtkList = matches[2];
	//checking file names
	if (mocMeshFile.find(".apl") == std::string::npos)
	{
		Logger::LogError("in MOCFieldsToCFD, " + mocMeshFile + " is not an .apl file");
	}
	if (outMocMeshFile.find(".apl") == std::string::npos)
	{
		Logger::LogError("in MOCFieldsToCFD, " + outMocMeshFile + " is not an .apl file");
	}
	if (mocPowerFile.find(".txt") == std::string::npos)
	{
		Logger::LogError("in MOCFieldsToCFD, heat power file " + mocPowerFile + " is not a .txt file");
	}
	for (size_t i = 0;i < outputVtkList.size();i++)
	{
		if (inputVTKList[i].find(".vtk") == std::string::npos)
		{
			Logger::LogError("in MOCFieldsToCFD, " + inputVTKList[i] + " is not given as .vtk file");
		}
		if (outputVtkList[i].find(".vtk") == std::string::npos)
		{
			Logger::LogError("in MOCFieldsToCFD, " + outputVtkList[i] + " is not given as .vtk file");
		}
	}
	//MOCMesh mocMesh(mocMeshFile, outMocMeshFile, MeshKernelType::MHT_KERNEL);
	MOCMesh mocMesh(materialList);
	mocMesh.InitMOCHeatPower(mocPowerFile);
	//initialize with meshFile
	Scalar ratio = GetValueInFile(configFile, "scaleRatio");
	MHTVTKReader reader(inputVTKList, ratio);
	//read interpolation weights only in available regions
	for (size_t i = 0; i < inputVTKList.size(); i++)
	{
		Mesh* pmesh = reader.GetMeshListPtr()[i];
		Field<Scalar> heatpower(pmesh, 0.0, "heatpower");
		//read cfd mesh and create solver
		CFDMesh cfdMesh(pmesh, MeshKernelType::MHT_KERNEL, i);
		Solver solverMapper(mocMesh, cfdMesh, materialList[i]);
		solverMapper.MOCtoCFDinterception(ValueType::HEATPOWER);
		cfdMesh.SetFieldValue(heatpower.elementField.v_value, ValueType::HEATPOWER);
		std::string strOutput_inpName = outputVtkList[i];
		RenameFile(strOutput_inpName, GetFileNameOfPrevious(strOutput_inpName,"vtk"));
		Logger::LogInfo("Writing CFD vtk file: " + strOutput_inpName);
		heatpower.WriteVTK_Field(strOutput_inpName);
		ConservationValidation(mocMesh, cfdMesh, ValueType::HEATPOWER, materialList[i], solverMapper);
	}
	return;
}

void CFDFieldsToMOC()
{
	std::vector<std::vector<std::string> > matches = GetMatchList(configFile);
	std::vector<std::string>& materialList = matches[0];
	std::vector<std::string>& inputVTKList = matches[1];
	std::string cfdMeshFile = GetFileName(configFile, "inputMsh");
	std::string mocMeshFile = GetFileName(configFile, "inputApl");
	std::string outMocMeshFile = GetFileName(configFile, "outputApl");
	std::string mocFieldFile = GetFileName(configFile, "inputInp");
	std::string outMocFieldFile = GetFileName(configFile, "outputInp");

	//checking file names
	if (cfdMeshFile.find(".msh") == std::string::npos)
	{
		Logger::LogError("in CFDFieldsToMOC, " + cfdMeshFile + " is not a .msh file");
	}
	for (size_t i = 0;i < inputVTKList.size();i++)
	{
		if (inputVTKList[i].find(".vtk") == std::string::npos)
		{
			Logger::LogError("in CFDFieldsToMOC, " + inputVTKList[i] + " is not a .vtk file");
		}
	}
	if (mocMeshFile.find(".apl") == std::string::npos)
	{
		Logger::LogError("in CFDFieldsToMOC, " + mocMeshFile + " is not an .apl file");
	}
	if (mocFieldFile.find(".inp") == std::string::npos)
	{
		Logger::LogError("in CFDFieldsToMOC, " + mocFieldFile + " is not an .inp file");
	}
	if (outMocFieldFile.find(".inp") == std::string::npos)
	{
		Logger::LogError("in CFDFieldsToMOC, " + outMocFieldFile + " is not an .inp file");
	}
	MOCMesh mocMesh(materialList);
	mocMesh.InitMOCFromInputFile(mocFieldFile);
	std::vector<std::string> fieldName;
	fieldName.push_back("T");
	fieldName.push_back("Rho");
	//initialize with meshFile
	Scalar ratio = GetValueInFile(configFile, "scaleRatio");
	MHTVTKReader reader(inputVTKList, ratio);
	//read CFD field with given ID list
	reader.ReadVTKFile(inputVTKList, fieldName);
	//read interpolation weights only in available regions
	for (size_t i = 0; i < inputVTKList.size(); i++)
	{
		Mesh* pmesh = reader.GetMeshListPtr()[i];
		//read cfd mesh and create solver
		CFDMesh cfdMesh(pmesh, MeshKernelType::MHT_KERNEL, i);
		for (int j = 0; j < reader.GetFieldList()[i].size(); j++)
		{
			const Field<Scalar>& field = reader.GetField(i,j);
			if (field.st_name == "T")
			{
				cfdMesh.SetValueVec(field.elementField.v_value, ValueType::TEMPERAURE);
			}
			else if (field.st_name == "Rho")
			{
				cfdMesh.SetValueVec(field.elementField.v_value, ValueType::DENSITY);
			}
		}
		Solver solverMapper(mocMesh, cfdMesh, materialList[i]);
		solverMapper.CFDtoMOCinterception(ValueType::DENSITY);
		solverMapper.CFDtoMOCinterception(ValueType::TEMPERAURE);
		ConservationValidation(cfdMesh, mocMesh, ValueType::DENSITY, materialList[i], solverMapper);
		ConservationValidation(cfdMesh, mocMesh, ValueType::TEMPERAURE, materialList[i], solverMapper);
	}
	RenameFile(outMocFieldFile, GetFileNameOfPrevious(outMocFieldFile, "inp"));
	mocMesh.OutputStatus(outMocFieldFile);
	return;
}

void ClearMapFiles()
{
	std::vector<std::vector<std::string> > matches = GetMatchList(configFile);
	std::vector<std::string>& materialList = matches[0];
	for (int i = 0;i < materialList.size();i++)
	{
		RemoveFile("MapFile_" + materialList[i] + "_CFDtoMOC");
		RemoveFile("MapFile_" + materialList[i] + "_MOCtoCFD");
	}
	RemoveFile(configFile);
	return;
}