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

void DisplayHelpInfo()
{
	std::cout << "Please input command like one list below:" << std::endl;
	std::cout << "1. NCCouple" << std::endl;
	std::cout << "2. NCCouple --help" << std::endl;
	std::cout << "3. NCCouple clear" << std::endl;
	std::cout << "4. NCCouple createmapper (MOCMesh) (MOCField) (CFDMesh)" << std::endl;
	std::cout << "5. NCCouple cfdtomoc (CFDMesh) (CFDField) (MOCMesh) (MOCField)" << std::endl;
	std::cout << "6. NCCouple moctocfd (MOCMesh) (MOCField) (MOCPower) (CFDMesh) (CFDField)" << std::endl;
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

void ReadAvalaibleRegionIDs(std::vector<int>& IDList)
{
	std::string fileName = "MapFile_Regions";
	ifstream infile(fileName);
	if (!infile.is_open())
	{
		Logger::LogError("cannot find the region record file: " + fileName);
		return;
	}
	IDList.clear();
	int RegionNum = 0;
	infile >> RegionNum;
	for (size_t i = 0;i < RegionNum;i++)
	{
		int RegionID = 0;
		infile >> RegionID;
		IDList.push_back(RegionID);
	}
	infile.close();
	return;
}

void CreateMapper
(
	std::string strInput_aplFileName,
	std::string strInput_inpFileName,
	std::string strInput_meshFileName
)
{
	//checking file names
	if (strInput_aplFileName.find(".apl") == std::string::npos)
	{
		Logger::LogError("in CreateMapper, the 1st file " + strInput_aplFileName + " is not an .apl file");
	}
	if (strInput_inpFileName.find(".inp") == std::string::npos)
	{
		Logger::LogError("in CreateMapper, the 2nd file " + strInput_aplFileName + " is not an .inp file");
	}
	if (strInput_meshFileName.find(".msh") == std::string::npos)
	{
		Logger::LogError("in CreateMapper, the 3rd file " + strInput_meshFileName + " is not a .msh file");
	}
	Logger::LogInfo("Reading MOC files: " + strInput_aplFileName + ", " + strInput_inpFileName);
	MOCMesh mocMesh(strInput_aplFileName, strInput_inpFileName, MeshKernelType::MHT_KERNEL);
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
	std::vector<int> availableRegionList;
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
				availableRegionList.push_back(j);
				Logger::LogInfo("available region found: " + CFDRegionList[j]);
			}
		}
	}
	if (0 == availableRegionList.size())
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
		std::ofstream regionIDRecord("MapFile_Regions");
		regionIDRecord << availableRegionList.size() << std::endl;
		for (size_t i = 0;i < availableRegionList.size();i++)
		{
			regionIDRecord << availableRegionList[i] << "\t";
		}
		regionIDRecord.close();
	}

	for (size_t i = 0; i < availableRegionList.size(); i++)
	{
		int CFDMeshID = availableRegionList[i];
		Mesh* pmesh = reader.GetMeshListPtr()[CFDMeshID];
		//read cfd mesh and create solver
		CFDMesh cfdMesh(pmesh, MeshKernelType::MHT_KERNEL, i);
		Solver solverMapper;
		solverMapper = Solver(mocMesh, cfdMesh, mocIndex, pmesh->st_meshName);
		solverMapper.CheckMappingWeights();
	}
	return;
}

void MOCFieldsToCFD
(
	std::string strInput_aplFileName,
	std::string strInput_inpFileName,
	std::string strInput_txtFileName,
	std::string strInput_meshFileName,
	std::string strOutput_vtkFileName
)
{
	//checking file names
	if (strInput_aplFileName.find(".apl") == std::string::npos)
	{
		Logger::LogError("in MOCFieldsToCFD, the 1st file " + strInput_aplFileName + " is not an .apl file");
	}
	if (strInput_inpFileName.find(".inp") == std::string::npos)
	{
		Logger::LogError("in MOCFieldsToCFD, the 2nd file " + strInput_inpFileName + " is not an .inp file");
	}
	if (strInput_txtFileName.find(".txt") == std::string::npos)
	{
		Logger::LogError("in MOCFieldsToCFD, the 3rd file " + strInput_txtFileName + " is not a .txt file");
	}
	if (strInput_meshFileName.find(".msh") == std::string::npos)
	{
		Logger::LogError("in MOCFieldsToCFD, the 4th file " + strInput_meshFileName + " is not a .msh file");
	}
	if (strOutput_vtkFileName.find(".vtk") == std::string::npos)
	{
		Logger::LogError("in MOCFieldsToCFD, the 5th file " + strOutput_vtkFileName + " is not given as .vtk file");
	}
	Logger::LogInfo("****** Mapping from MOC to CFD ******");
	Logger::LogInfo("Reading MOC files: " + strInput_aplFileName + ", " + strInput_inpFileName);
	MOCMesh mocMesh(strInput_aplFileName, strInput_inpFileName, MeshKernelType::MHT_KERNEL);
	Logger::LogInfo("Reading MOC heat power file: " + strInput_txtFileName);
	mocMesh.InitMOCHeatPower(strInput_txtFileName);
	Logger::LogInfo("Reading CFD mesh file: " + strInput_meshFileName);
	//initialize with meshFile
	MHTVTKReader reader(strInput_meshFileName);
	//reading available region IDs
	std::vector<int> avaiableRegionIDList;
	ReadAvalaibleRegionIDs(avaiableRegionIDList);
	//read interpolation weights only in available regions
	for (size_t i = 0; i < avaiableRegionIDList.size(); i++)
	{
		int RegionID = avaiableRegionIDList[i];
		Mesh* pmesh = reader.GetMeshListPtr()[RegionID];
		Field<Scalar> heatpower(pmesh, 0.0, "heatpower");
		//read cfd mesh and create solver
		CFDMesh cfdMesh(pmesh, MeshKernelType::MHT_KERNEL, i);
		Solver solverMapper;
		solverMapper = Solver(mocMesh, cfdMesh, pmesh->st_meshName);
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
	std::string strInput_aplFileName,
	std::string strInput_inpFileName
)
{
	//checking file names
	if (strInput_meshFileName.find(".msh") == std::string::npos)
	{
		Logger::LogError("in CFDFieldsToMOC, the 1st file " + strInput_meshFileName + " is not a .msh file");
	}
	if (strInput_vtkFileName.find(".txt") == std::string::npos)
	{
		Logger::LogError("in CFDFieldsToMOC, the 2nd file " + strInput_vtkFileName + " is not a .txt file");
	}
	if (strInput_aplFileName.find(".apl") == std::string::npos)
	{
		Logger::LogError("in CFDFieldsToMOC, the 3rd file " + strInput_aplFileName + " is not an .apl file");
	}
	if (strInput_inpFileName.find(".inp") == std::string::npos)
	{
		Logger::LogError("in CFDFieldsToMOC, the 4th file " + strInput_inpFileName + " is not an .inp file");
	}
	MOCMesh mocMesh(strInput_aplFileName, strInput_inpFileName,MeshKernelType::MHT_KERNEL);
	std::vector<std::string> fieldName;
	fieldName.push_back("T");
	fieldName.push_back("Rho");
	MHTVTKReader reader(strInput_meshFileName,strInput_vtkFileName, fieldName);//initialize with meshFile
	//reading available region IDs
	std::vector<int> avaiableRegionIDList;
	ReadAvalaibleRegionIDs(avaiableRegionIDList);
	//read interpolation weights only in available regions
	for (size_t i = 0; i < avaiableRegionIDList.size(); i++)
	{
		int RegionID = avaiableRegionIDList[i];
		Mesh* pmesh = reader.GetMeshListPtr()[RegionID];
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
		solverMapper = Solver(mocMesh, cfdMesh, pmesh->st_meshName);
		solverMapper.CFDtoMOCinterception(ValueType::DENSITY);
		solverMapper.CFDtoMOCinterception(ValueType::TEMPERAURE);
		ConservationValidation(cfdMesh, mocMesh, ValueType::DENSITY,pmesh->st_meshName);
		ConservationValidation(cfdMesh, mocMesh, ValueType::TEMPERAURE,pmesh->st_meshName);
	}
	std::string strOutput_inpName = strInput_inpFileName.substr(0, strInput_inpFileName.find(".")) + "_out.inp";
	mocMesh.OutputStatus(strOutput_inpName);

}