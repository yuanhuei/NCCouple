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
	double sourceValue = 0, pointVolume = 0,cfdTotalVolume=0,mocTotalVolume=0;
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
			std::vector<int>vCFDID;
			double dValue=solverMapper.GetMocMeshMapValue(vSMocIndex[i]);
			double dTotalValue = solverMapper.GetMocMeshMapTotalValue(vSMocIndex[i]);

			pointVolume = mocMesh.GetVolumeAtIndex(vSMocIndex[i])* dValue/ dTotalValue;
			integration += sourceValue * pointVolume;


			mocTotalVolume += pointVolume;
		}
	}
	else
	{
		for (int i = 0; i < mesh.GetMeshPointNum(); i++)
		{
			sourceValue = mesh.GetMeshPointPtr(i)->GetValue(vt);
			pointVolume = mesh.GetMeshPointPtr(i)->Volume();
			integration += sourceValue * pointVolume;

			cfdTotalVolume += pointVolume;
		}
	}
	//if(g_iMpiID==1)
		//WriteToLog(FormatStr("the total cfd vulume is :%.6lf the total moc volume is:%.6lf", cfdTotalVolume, mocTotalVolume));
	return integration;
}

std::pair<double,double> ConservationValidation
(
	const GeneralMesh& sourceMesh,
	const GeneralMesh& targetMesh,
	ValueType vt,
	std::string strZoneName,
	Solver& solverMapper
)
{
	std::string valueName = NameOfValueType(vt);
	double sourceIntegralValue = Integration(sourceMesh,vt, strZoneName, solverMapper);
	double targetIntegralValue = Integration(targetMesh,vt, strZoneName, solverMapper);
	Logger::LogInfo(FormatStr("Process %d,Integral of %s on region %s of source mesh: %.6lf",g_iMpiID, valueName.c_str(), strZoneName.c_str(), sourceIntegralValue));
	Logger::LogInfo(FormatStr("Process %d,Integral of %s on region %s of target mesh: %.6lf", g_iMpiID,valueName.c_str(), strZoneName.c_str(), targetIntegralValue));
	//WriteToLog(FormatStr("Process %d,Integral of %s on region %s of source mesh: %.6lf", g_iMpiID, valueName.c_str(), strZoneName.c_str(), sourceIntegralValue), g_iMpiID);
	//WriteToLog(FormatStr("Process %d,Integral of %s on region %s of target mesh: %.6lf", g_iMpiID, valueName.c_str(), strZoneName.c_str(), targetIntegralValue),g_iMpiID);

	return std::make_pair(sourceIntegralValue,targetIntegralValue);
}

void DisplayHelpInfo()
{
	std::cout << "Please input command like one list below:" << std::endl;
	std::cout << "1. NCCouple --help" << std::endl;
	std::cout << "2. NCCouple clear" << std::endl;
	std::cout << "3. mpirun -np (n) NCCouple createmapper" << std::endl;
	std::cout << "4. mpirun -np (n) NCCouple cfdtomoc" << std::endl;
	std::cout << "5. mpirun -np (n) NCCouple moctocfd" << std::endl;
	std::cout << "6. NCCouple register (MOCMesh) (MOCMesh_out)" << std::endl;
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
	//std::string strInput_meshFileName
//)
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
	for (size_t i = 0; i < mocMesh.m_vSMocIndex.size(); i++)
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
		int iPos = inputVTKList[i].find(".vtk");
		if (iPos == std::string::npos)
		{
			Logger::LogError("in MOCFieldsToCFD, " + inputVTKList[i] + " is not a .vtk file");
		}
		inputVTKList[i] = inputVTKList[i].substr(0, iPos) + "_" + std::to_string(g_iMpiID) + ".vtk";
	}
	MOCMesh mocMesh(mocMeshFile, outMocMeshFile, MeshKernelType::MHT_KERNEL,true);
	//create an index for fast searching
	Logger::LogInfo("Reading CFD mesh from input vtk files");
	Scalar ratio = GetValueInFile(configFile, "scaleRatio");
	MHTVTKReader reader(inputVTKList, ratio,true);
	for (size_t i = 0; i < inputVTKList.size(); i++)
	{
		Mesh* pmesh = reader.GetMeshListPtr()[i];
		//read cfd mesh and create solver
		CFDMesh cfdMesh(pmesh, MeshKernelType::MHT_KERNEL, i);
		Solver solverMapper(mocMesh, cfdMesh, materialList[i],true);
		solverMapper.CheckMappingWeights();

		//send to mainprocess
		std::vector<STRMocMapValue> vMocMapValue;
		solverMapper.SetMocMapValueToStruct(vMocMapValue);
		SendMocMapValueToMainProcess(vMocMapValue);
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
	for (size_t i = 0; i < outputVtkList.size(); i++)
	{
		if (inputVTKList[i].find(".vtk") == std::string::npos)
		{
			Logger::LogError("in MOCFieldsToCFD, " + inputVTKList[i] + " is not given as .vtk file");
		}
		if (outputVtkList[i].find(".vtk") == std::string::npos)
		{
			Logger::LogError("in MOCFieldsToCFD, " + outputVtkList[i] + " is not given as .vtk file");
		}
		int iPos = inputVTKList[i].find(".vtk");
		inputVTKList[i] = inputVTKList[i].substr(0, iPos) + "_" + std::to_string(g_iMpiID) + ".vtk";
		iPos = outputVtkList[i].find(".vtk");
		outputVtkList[i] = outputVtkList[i].substr(0, iPos) + "_" + std::to_string(g_iMpiID) + ".vtk";
	}

	//MOCMesh mocMesh(mocMeshFile, outMocMeshFile, MeshKernelType::MHT_KERNEL);
	MOCMesh mocMesh(materialList);
	mocMesh.InitMOCHeatPower(mocPowerFile,true);
	//initialize with meshFile
	Scalar ratio = GetValueInFile(configFile, "scaleRatio");
	MHTVTKReader reader(inputVTKList, ratio);
	//read interpolation weights only in available regions
	for (size_t i = 0; i < materialList.size(); i++)
	{
		Mesh* pmesh = reader.GetMeshListPtr()[i];
		Field<Scalar> heatpower(pmesh, 0.0, "heatpower");
		//read cfd mesh and create solver
		CFDMesh cfdMesh(pmesh, MeshKernelType::MHT_KERNEL, i);
		Solver solverMapper(mocMesh, cfdMesh, materialList[i]);
		Logger::LogInfo(" solverMapper  ");
		solverMapper.MOCtoCFDinterception(ValueType::HEATPOWER);
		cfdMesh.SetFieldValue(heatpower.elementField.v_value, ValueType::HEATPOWER);
		std::string strOutput_inpName = outputVtkList[i];
		RenameFile(strOutput_inpName, GetFileNameOfPrevious(strOutput_inpName,"vtk"));
		Logger::LogInfo("Writing CFD vtk file: " + strOutput_inpName);
		heatpower.WriteVTK_Field(strOutput_inpName);
		std::pair<double, double> pairHeatpower;
		pairHeatpower=ConservationValidation(mocMesh, cfdMesh, ValueType::HEATPOWER, materialList[i], solverMapper);

		ofstream oFile("CFD_HeatPower_" + materialList[i] + "_" + std::to_string(g_iMpiID));
		oFile << pairHeatpower.second << std::endl;
		oFile.close();
	}
	return;
}

void CFDFieldsToMOC()
{
	std::vector<std::vector<std::string> > matches = GetMatchList(configFile);
	std::vector<std::string>& materialList = matches[0];
	std::vector<std::string>& inputVTKList = matches[1];
	std::string mocMeshFile = GetFileName(configFile, "inputApl");
	std::string outMocMeshFile = GetFileName(configFile, "outputApl");
	std::string mocFieldFile = GetFileName(configFile, "inputInp");
	std::string outMocFieldFile = GetFileName(configFile, "outputInp");

	//checking file names
	for (size_t i = 0;i < inputVTKList.size();i++)
	{
		if (inputVTKList[i].find(".vtk") == std::string::npos)
		{
			Logger::LogError("in CFDFieldsToMOC, " + inputVTKList[i] + " is not a .vtk file");
		}
		int iPos = inputVTKList[i].find(".vtk");
		inputVTKList[i] = inputVTKList[i].substr(0, iPos) + "_" + std::to_string(g_iMpiID) + ".vtk";
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
	fieldName.push_back("temperature");
	fieldName.push_back("rho");
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
			if (field.st_name == "temperature")
			{
				cfdMesh.SetValueVec(field.elementField.v_value, ValueType::TEMPERAURE);
			}
			else if (field.st_name == "rho")
			{
				cfdMesh.SetValueVec(field.elementField.v_value, ValueType::DENSITY);
			}
		}
		Solver solverMapper(mocMesh, cfdMesh, materialList[i]);
		solverMapper.CFDtoMOCinterception(ValueType::DENSITY);
		solverMapper.CFDtoMOCinterception(ValueType::TEMPERAURE);
		std::pair<double, double> pairDensity, pairTemerature;
		pairDensity=ConservationValidation(cfdMesh, mocMesh, ValueType::DENSITY, materialList[i], solverMapper);
		pairTemerature=ConservationValidation(cfdMesh, mocMesh, ValueType::TEMPERAURE, materialList[i], solverMapper);

		ofstream oFile("CFD_VALUE_"+ materialList[i]+"_"+ std::to_string(g_iMpiID));
		//ofstream oFile2("CFD_TEPERATURE_VALUE_" + materialList[i] + "_" + std::to_string(g_iMpiID));
		oFile << pairDensity.first << " " << pairTemerature.first << std::endl;
		oFile.close();

		 oFile.open("MOC_VALUE_" + materialList[i] + "_" + std::to_string(g_iMpiID));
		//ofstream oFile2("CFD_TEPERATURE_VALUE_" + materialList[i] + "_" + std::to_string(g_iMpiID));
		oFile << pairDensity.second << " " << pairTemerature.second << std::endl;
		oFile.close();
	}
	//RenameFile(outMocFieldFile, GetFileNameOfPrevious(outMocFieldFile, "inp"));
	//mocMesh.OutputStatus(outMocFieldFile);
	// solution 2
	//mocMesh.WriteFieldInfortoFile();

	//solution 1
	std::vector<STRMocField> vMocField;
	mocMesh.SetFieldValueToMpiType(vMocField);
	SendMocValueToMainProcess(vMocField);
	return;
}

void ClearMapFiles()
{
	std::vector<std::vector<std::string> > matches = GetMatchList(configFile);
	std::vector<std::string>& materialList = matches[0];
	for (int i = 0;i < materialList.size();i++)
	{
		RemoveFile("./temp/MapFile_" + materialList[i] + "_CFDtoMOC");
		RemoveFile("./temp/MapFile_" + materialList[i] + "_MOCtoCFD");
	}
	RemoveFile(configFile);
	return;
}