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
	std::cout << "4. NCCouple createmapper (MOCMesh) (MOCMesh) (CFDMesh)" << std::endl;
	std::cout << "5. NCCouple cfdtomoc" << std::endl;
	std::cout << "6. NCCouple moctocfd" << std::endl;
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

void CreateMapper
(
	std::string strInput_aplFileName,
	std::string strOutput_aplFileName,
	std::string strInput_meshFileName
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
	if (strInput_meshFileName.find(".msh") == std::string::npos)
	{
		Logger::LogError("in CreateMapper, the CFD mesh " + strInput_meshFileName + " is not a .msh file");
	}
	MOCMesh mocMesh(strInput_aplFileName, strOutput_aplFileName, MeshKernelType::MHT_KERNEL);
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
	MHTVTKReader reader(strInput_meshFileName);
	std::vector<std::string> MOCRegionList;
	std::vector<std::string> CFDRegionList;
	std::vector<int> availableRegionID;
	std::vector<std::string> availableRegionName;
	//collect region names from MOC mesh
	for (size_t i = 0;i < mocMesh.GetMeshPointNum();i++)
	{
		MOCMeshPoint* mocPoint = dynamic_cast<MOCMeshPoint*>(mocMesh.GetMeshPointPtr(i));
		std::string thisName = mocPoint->GetMaterialName();
		InsertWhenNotFound(MOCRegionList, thisName);
	}
	//collect region names from CFD mesh
	for (size_t i = 0; i < reader.GetMeshListPtr().size(); i++)
	{
		Mesh* pmesh = reader.GetMeshListPtr()[i];
		std::string thisName = pmesh->st_meshName;
		CFDRegionList.push_back(thisName);
	}
	for (size_t i = 0;i < MOCRegionList.size();i++)
	{
		for (size_t j = 0;j < CFDRegionList.size();j++)
		{
			if (MOCRegionList[i] == CFDRegionList[j])
			{
				availableRegionID.push_back(j);
				availableRegionName.push_back(CFDRegionList[j]);
				Logger::LogInfo("available region found: " + CFDRegionList[j]);
			}
		}
	}
	if (0 == availableRegionID.size())
	{
		Logger::LogInfo("no available region found!");
		std::cout << "In MOC, we have regions: " << std::endl;
		for (size_t i = 0;i < MOCRegionList.size();i++)
		{
			std::cout << "(" << i+1 << ") " << MOCRegionList[i] << std::endl;
		}
		std::cout << "In CFD, we have regions: " << std::endl;
		for (size_t i = 0;i < CFDRegionList.size();i++)
		{
			std::cout << "(" << i+1 << ") " << CFDRegionList[i] << std::endl;
		}
		exit(EXIT_FAILURE);
	}
	else
	{
		WriteConfigurationFile(configFile, strInput_meshFileName, strInput_aplFileName, strOutput_aplFileName, availableRegionID, availableRegionName);
	}
	for (size_t i = 0; i < availableRegionID.size(); i++)
	{
		int CFDMeshID = availableRegionID[i];
		Mesh* pmesh = reader.GetMeshListPtr()[CFDMeshID];
		//read cfd mesh and create solver
		CFDMesh cfdMesh(pmesh, MeshKernelType::MHT_KERNEL, i);
		Solver solverMapper(mocMesh, cfdMesh, mocIndex, pmesh->st_meshName);
		solverMapper.CheckMappingWeights();
	}
	return;
}

void MOCFieldsToCFD()
{
	std::string mocMeshFile = GetFileName(configFile, "inputApl");
	std::string outMocMeshFile = GetFileName(configFile, "outputApl");
	std::string mocPowerFile = GetFileName(configFile, "mocPower");
	std::string cfdMeshFile = GetFileName(configFile, "inputMsh");
	std::vector<std::string> outputVtkList = GetFileNameList(configFile, "outputVtk");
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
	if (cfdMeshFile.find(".msh") == std::string::npos)
	{
		Logger::LogError("in MOCFieldsToCFD, " + cfdMeshFile + " is not a .msh file");
	}
	for (size_t i = 0;i < outputVtkList.size();i++)
	{
		if (outputVtkList[i].find(".vtk") == std::string::npos)
		{
			Logger::LogError("in MOCFieldsToCFD, " + outputVtkList[i] + " is not given as .vtk file");
		}
	}
	//reading available region IDs
	std::vector<int> avaiableRegionIDList = GetRegionIDList(configFile);
	if (avaiableRegionIDList.size() != outputVtkList.size())
	{
		std::cout << "In previously generated mapper, " << avaiableRegionIDList.size() << " regions are available, ";
		std::cout << "but you specified " << outputVtkList.size() << " vtk files as output" << std::endl;
		Logger::LogError("numbers of vtk files and avaiable regions are not consistent");
	}
	MOCMesh mocMesh(mocMeshFile, outMocMeshFile, MeshKernelType::MHT_KERNEL);
	mocMesh.InitMOCHeatPower(mocPowerFile);
	//initialize with meshFile
	MHTVTKReader reader(cfdMeshFile);
	//read interpolation weights only in available regions
	for (size_t i = 0; i < avaiableRegionIDList.size(); i++)
	{
		int RegionID = avaiableRegionIDList[i];
		Mesh* pmesh = reader.GetMeshListPtr()[RegionID];
		Field<Scalar> heatpower(pmesh, 0.0, "heatpower");
		//read cfd mesh and create solver
		CFDMesh cfdMesh(pmesh, MeshKernelType::MHT_KERNEL, RegionID);
		Solver solverMapper;
		solverMapper = Solver(mocMesh, cfdMesh, pmesh->st_meshName);
		solverMapper.MOCtoCFDinterception(ValueType::HEATPOWER);
		cfdMesh.SetFieldValue(heatpower.elementField.v_value, ValueType::HEATPOWER);
		std::string strOutput_inpName = outputVtkList[i];
		RenameFile(strOutput_inpName, GetFileNameOfPrevious(strOutput_inpName,"vtk"));
		Logger::LogInfo("Writing CFD vtk file: " + strOutput_inpName);
		heatpower.WriteVTK_Field(strOutput_inpName);
		ConservationValidation(mocMesh, cfdMesh, ValueType::HEATPOWER,pmesh->st_meshName);
	}
	return;
}

void CFDFieldsToMOC()
{

	std::string cfdMeshFile = GetFileName(configFile, "inputMsh");
	std::vector<std::string> vtkFileName = GetFileNameList(configFile, "inputVtk");
	std::string mocMeshFile = GetFileName(configFile, "inputApl");
	std::string outMocMeshFile = GetFileName(configFile, "outputApl");
	std::string mocFieldFile = GetFileName(configFile, "inputInp");
	std::string outMocFieldFile = GetFileName(configFile, "outputInp");

	//checking file names
	if (cfdMeshFile.find(".msh") == std::string::npos)
	{
		Logger::LogError("in CFDFieldsToMOC, " + cfdMeshFile + " is not a .msh file");
	}
	for (size_t i = 0;i < vtkFileName.size();i++)
	{
		if (vtkFileName[i].find(".vtk") == std::string::npos)
		{
			Logger::LogError("in CFDFieldsToMOC, " + vtkFileName[i] + " is not a .vtk file");
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
	//reading available region IDs
	std::vector<int> avaiableRegionIDList = GetRegionIDList(configFile);
	if (avaiableRegionIDList.size() != vtkFileName.size())
	{
		std::cout << "In previously generated mapper, " << avaiableRegionIDList.size() << " regions are available, ";
		std::cout << "but you specified " << vtkFileName.size() << " vtk files as input" << std::endl;
		Logger::LogError("numbers of vtk files and avaiable regions are not consistent");
	}
	MOCMesh mocMesh(mocMeshFile, outMocMeshFile,MeshKernelType::MHT_KERNEL);
	mocMesh.InitMOCFromInputFile(mocFieldFile);
	std::vector<std::string> fieldName;
	fieldName.push_back("T");
	fieldName.push_back("Rho");
	//initialize with meshFile
	MHTVTKReader reader(cfdMeshFile);
	//read CFD field with given ID list
	reader.ReadVTKFile(vtkFileName, avaiableRegionIDList, fieldName);
	//read interpolation weights only in available regions
	for (size_t i = 0; i < avaiableRegionIDList.size(); i++)
	{
		int RegionID = avaiableRegionIDList[i];
		Mesh* pmesh = reader.GetMeshListPtr()[RegionID];
		//read cfd mesh and create solver
		CFDMesh cfdMesh(pmesh, MeshKernelType::MHT_KERNEL, RegionID);
		for (int j = 0; j < reader.GetFieldList()[RegionID].size(); j++)
		{
			const Field<Scalar>& field = reader.GetField(RegionID,j);
			if (field.st_name == "T")
			{
				cfdMesh.SetValueVec(field.elementField.v_value, ValueType::TEMPERAURE);
			}
			else if (field.st_name == "Rho")
			{
				cfdMesh.SetValueVec(field.elementField.v_value, ValueType::DENSITY);
			}
		}
		Solver solverMapper;
		solverMapper = Solver(mocMesh, cfdMesh, pmesh->st_meshName);
		solverMapper.CFDtoMOCinterception(ValueType::DENSITY);
		solverMapper.CFDtoMOCinterception(ValueType::TEMPERAURE);
		ConservationValidation(cfdMesh, mocMesh, ValueType::DENSITY,pmesh->st_meshName);
		ConservationValidation(cfdMesh, mocMesh, ValueType::TEMPERAURE,pmesh->st_meshName);
	}
	RenameFile(outMocFieldFile, GetFileNameOfPrevious(outMocFieldFile, "inp"));
	mocMesh.OutputStatus(outMocFieldFile);
	return;
}

void ClearMapFiles()
{
	std::string cfdMeshFile = GetFileName(configFile, "inputMsh");
	std::vector<std::string> regionNameList = GetFileNameList(configFile, "regionName");
	for (int i = 0;i < regionNameList.size();i++)
	{
		RemoveFile("MapFile_" + regionNameList[i] + "_CFDtoMOC");
		RemoveFile("MapFile_" + regionNameList[i] + "_MOCtoCFD");
	}
	RemoveFile(configFile);
	return;
}