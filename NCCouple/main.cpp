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
#include "./ConfigurationFile.h"

#include <sstream>
#include <mpi.h>
#ifdef _WIN32
#include <process.h>
#include <string>
#else
#include <unistd.h>
#include <string>
#include <string.h>
#endif
#include "./MHT_common/SystemControl.h"
#include "./MHT_mesh/UnGridFactory.h"
#include "./MHT_mesh/RegionConnection.h"
#include "./MHT_mesh/Mesh.h"
#include "./MHT_field/Field.h"
#include "./MHT_IO/FieldIO.h"
#include "./MHT_IO/VTKIO.h"

int g_iMpiID = -1;
int g_iNumProcs = 0;

//this example was designed for test of
//(1) rewriting a apl file
//(2) reading and writting of inp files

void WriteTotecplot( MOCMesh& mocmesh, CFDMesh& cfdmesh,
	const std::vector<int> &v_iCFD,const std::vector<SMocIndex>& vSMocindex,std::string fileName)
{
	std::string  mType;
	std::ofstream ofile(fileName);
	ofile << "TITLE =\"" << "polyhedron" << "\"" << endl;
	ofile << "VARIABLES = " << "\"x\"," << "\"y\"," << "\"z\"" << endl;
	for (int i = 0; i < vSMocindex.size(); i++)
	{

		//move coordinate
		/*
		int iAssembly = vSMocindex[i].iAssemblyIndex;
		double x = mocmesh.m_vAssembly[iAssembly].pAssembly_type->vAssemblyType_LeftDownPoint.x_- mocmesh.m_vAssembly[iAssembly].vAssembly_LeftDownPoint.x_;
		double y= mocmesh.m_vAssembly[iAssembly].pAssembly_type->vAssemblyType_LeftDownPoint.y_-mocmesh.m_vAssembly[iAssembly].vAssembly_LeftDownPoint.y_;
		
		const MHTMocMeshPoint& mocPoint = dynamic_cast<const MHTMocMeshPoint&>(*mocmesh.GetMocMeshPointPtr(vSMocindex[i]));
		*/
		shared_ptr<MHTMocMeshPoint> ptrMocMesh;
		mocmesh.MoveMeshToSysCoordinate(ptrMocMesh, vSMocindex[i]);
		MHTMocMeshPoint& mocPoint = *ptrMocMesh;

		//MHTMocMeshPoint meshPoint = mocPoint;
		//Vector vPoint(-x, -y, 0);
		//meshPoint.Move(vPoint);
		mocPoint.WriteTecplotZones(ofile);
	}
	for (int i = 0; i < v_iCFD.size(); i++)
	{
		const MHTMeshPoint& mhtPolyhedron = dynamic_cast<const MHTMeshPoint&>(*cfdmesh.GetMeshPointPtr(v_iCFD[i]));
		mhtPolyhedron.WriteTecplotZones(ofile);
	}

	ofile.close();
	return;
}

void MOC_APL_INP_FileTest()
{
	WarningContinue("MOC_APL_INP_FileTest");
	MOCMesh mocMesh("3_3cells.apl", "3_3cells_out.apl", MeshKernelType::MHT_KERNEL);
	//mocMesh.WriteSurfaceTecplotFile("filename.plt");
	//MOCMesh mocMesh("pin_c1.apl", "pin_c1.inp", MeshKernelType::MHT_KERNEL); 
	mocMesh.WriteTecplotFile("3_3cells_mMOD.plt","mMOD");
	mocMesh.WriteTecplotFile("3_3cells_mFUEL.plt","mFUEL");
	//mocMesh.WriteTecplotFile("3_3cells_mCLAD.plt", "mCLAD");
	//mocMesh.WriteTecplotFile("mUO2", "c5g721_cell12_mUO2.plt");
	//mocMesh.WriteTecplotFile("Zr4", "zr4.plt");
	//mocMesh.WriteTecplotFile("UO2", "u2o.plt");
	//mocMesh.InitMOCValue("pin_c1.inp","pin_c1.txt");
	//mocMesh.OutputStatus("pin_c1_out.inp");
	return;
}
void MapTest()
{
	//RegisterMapper("c5g72l.apl", "c5g72l_out.apl", "PIN9Coarse.msh");
	//MOCFieldsToCFD();

	std::string configFile = "MapFile_FileNames";
	std::vector<std::vector<std::string> > matches = GetMatchList(configFile);
	std::vector<std::string>& materialList = matches[0];
	std::vector<std::string>& regionList = matches[1];
	std::string mocMeshFile = GetFileName(configFile, "inputApl");
	std::string outMocMeshFile = GetFileName(configFile, "outputApl");

	MOCMesh mocMesh(mocMeshFile, outMocMeshFile);
	std::vector<std::string> inputVTKList;
	inputVTKList.push_back("fluid.vtk");
	inputVTKList.push_back("solid.vtk");
	Scalar ratio = GetValueInFile(configFile, "scaleRatio");
	MHTVTKReader reader(inputVTKList, ratio);
	for (size_t i = 0; i < inputVTKList.size(); i++)
	{
		Mesh* pmesh = reader.GetMeshListPtr()[i];
		//read cfd mesh and create solver
		CFDMesh cfdMesh(pmesh, MeshKernelType::MHT_KERNEL, i);
		cfdMesh.WriteTecplotFile("cfd_50_"+std::to_string(i)+".plt");
	}

	/*
	MHTVTKReader reader(cfdMeshFile, ratio);
	for (size_t i = 0; i < regionList.size(); i++)
	{
		//if (i > 0)break;
		int CFDMeshID = reader.GetIDOfRegion(regionList[i]);
		Mesh* pmesh = reader.GetMeshListPtr()[CFDMeshID];
		//read cfd mesh and create solver
		CFDMesh cfdMesh(pmesh, MeshKernelType::MHT_KERNEL, CFDMeshID);
		
		std::vector<SMocIndex> thh;
		for (int j = 0; j < 64; j++)
		{
			thh.push_back(SMocIndex(0, 51, j));
		}
		WriteTotecplot(mocMesh, cfdMesh, std::vector<int>{2559}, thh, "temp.plt");
		
		//std::string outfilename = "mesh_" + std::to_string(i) + ".plt";
		//cfdMesh.WriteTecplotFile(outfilename);
		Solver solverMapper(mocMesh, cfdMesh, materialList[i], true);
		solverMapper.CheckMappingWeights();
	}
	*/

	//InitCFDMeshValue(cfdMesh);
	//solver.CFDtoMOCinterception(ValueType::DENSITY);

	//mocMesh.OutputStatus("pin_c1.inp");

}
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

Scalar initialHeatPower(Scalar x, Scalar y, Scalar z)
{
	Vector PointCenter(x, y, z);
	Vector axisCenter(0.63, 0.63, 0.0);
	Vector axisNorm(0.0, 0.0, 1.0);
	Vector OP = PointCenter - axisCenter;
	Vector axicialLocation = (OP & axisNorm) * axisNorm;
	Scalar radius = (OP - axicialLocation).Mag();
	Scalar sigma = 1.0;
	Scalar hp = (1000.0 / sqrt(2.0 * PI) / sigma) * exp(-radius * radius / 2 / pow(sigma, 2));
	hp = hp * z * (5.0 - z) / 6.25;
	return hp;
}

#include "./MHT_IO/VTKIO.h"
#include <fstream>

void VTKReaderTest()
{
	MHTVTKReader reader("pinWR.msh",1.0);
	std::vector<std::string> fileName;
	fileName.push_back("H2O.vtk");
	fileName.push_back("Zr4.vtk");
	std::vector<int> IDList;
	IDList.push_back(0);
	IDList.push_back(1);
	std::vector<std::string> fieldName;
	fieldName.push_back("heatpower");
	reader.ReadVTKFile(fileName, IDList, fieldName);
	reader.GetFieldIO(0).WriteTecplotField("heatpower_0.plt");
	reader.GetFieldIO(1).WriteTecplotField("heatpower_1.plt");
	return;
}

void EntranceOfRegister(std::vector<std::string>& fileNames)
{
	if (2 != fileNames.size())
	{
		std::cout << "Please give 2 file names if you are intending to register solver, like this:" << std::endl;
		std::cout << "NCCouple register (MOCMesh) (outMOCMesh)" << std::endl;
		Logger::LogError("inccorrect number of file names");
	}
	RegisterMapper(fileNames[0], fileNames[1]);
	return;
}

void EntranceOfCreateMapper(std::vector<std::string>& fileNames)
{
	if (0 != fileNames.size())
	{
		std::cout << "Please do not give any file name if you are intending to create mapper, like this:" << std::endl;
std::cout << "NCCouple createmapper" << std::endl;
Logger::LogError("inccorrect number of file names");
	}
	CreateMapper();
	return;
}

void EntranceOfCFDToMOC(std::vector<std::string>& fileNames)
{
	if (0 != fileNames.size())
	{
		std::cout << "Please do not give any file name if you are intending to map CFD to MOC solver, like this:" << std::endl;
		std::cout << "NCCouple cfdtomoc" << std::endl;
		Logger::LogError("inccorrect number of file names");
	}
	CFDFieldsToMOC();
	return;
}

void EntranceOfMOCToCFD(std::vector<std::string>& fileNames)
{
	if (0 != fileNames.size())
	{
		std::cout << "Please do not give any file name if you are intending to map MOC to CFD solver, like this:" << std::endl;
		std::cout << "NCCouple moctocfd" << std::endl;
		Logger::LogError("inccorrect number of file names");
	}
	MOCFieldsToCFD();
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
	for (int i = 1; i < parameters.size(); i++)
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
		ClearMapFiles();
	}
	else if ("register" == parameters[0])
	{
		EntranceOfRegister(filenames);
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
		std::cout << "(5) register" << std::endl;
		std::cout << "(6) --help" << std::endl;
		Logger::LogError("invalid parameter given: " + parameters[0]);
	}
	return;
}
void writeheatpower()
{
	std::ofstream ofs("heatpower_1cm.txt");
	int iNumber_Assembly = 9;
	for (int i = 1; i <= iNumber_Assembly; i++)
	{
		int iNum = 48*50;
		double iValue = 1.9e+10;
		ofs << i << "_" << iNum << std::endl;
		for (int j = 0; j < 240; j++)
		{
			for(int k = 0; k < 10; k++)
			{
				ofs << iValue << "  ";
			}
			ofs << std::endl;
		}
		//for (int k = 0; k < (iNum-40); k++)
		//{
			//ofs << iValue << "  ";
		//}
		ofs << std::endl;
	}
}

void VTKReadMeshTest()
{
	std::vector<std::string> vVTKname;
	vVTKname.push_back("1_PART-FLUID_couple.vtk");
	vVTKname.push_back("2_solid.vtk");

	std::vector<std::string> vFieldName;
	vFieldName.push_back("temperature");
	MHTVTKReader mhtvtkreader(vVTKname, vFieldName, 1.0);
	std::cout<<"mhtvtkreader.GetMeshListPtr()[0]->v_vertice.size() :" << mhtvtkreader.GetMeshListPtr()[0]->v_vertice.size()<<std::endl;
	mhtvtkreader.GetFieldIO(0).WriteTecplotField("temperature_0.plt");
	mhtvtkreader.GetFieldIO(1).WriteTecplotField("temperature_1.plt");
}

void FinalCFDFieldsToMOC()
{
	std::vector<std::vector<std::string> > matches = GetMatchList(configFile);
	std::string mocFieldFile = GetFileName(configFile, "inputInp");
	std::string outMocFieldFile = GetFileName(configFile, "outputInp");

	//checking file names
	std::vector<std::string> vStrMocFileName;
	vStrMocFileName.resize(g_iNumProcs-1);
	for (int i = 0; i < g_iNumProcs-1; i++)
	{
		vStrMocFileName[i] = "MocField_" + std::to_string(i+1);
	}
		
	MOCMesh mocMesh(vStrMocFileName,true);
	mocMesh.InitMOCFromInputFile(mocFieldFile);
	RenameFile(outMocFieldFile, GetFileNameOfPrevious(outMocFieldFile, "inp"));
	mocMesh.OutputStatus(outMocFieldFile);
}
void CFDTOMOC_ValueValidation(MOCMesh& mocMesh, std::vector<std::string>& materialList)
{
	std::vector<Scalar> integrationDesityValue, integrationTemperatureValue;
	double sourceDesityValue = 0, sourceTemperatureValue = 0, pointVolume = 0;
	integrationDesityValue.resize(materialList.size());
	integrationTemperatureValue.resize(materialList.size());
	SMocIndex vSMocIndex;
	for (int i = 0; i < mocMesh.m_vSMocIndex.size(); i++)
	{
		vSMocIndex = mocMesh.m_vSMocIndex[i];
		int k = -1;
		for (int j = 0; j < materialList.size(); j++)
		{
			if (materialList[j] == mocMesh.GetMaterialNameAtIndex(vSMocIndex))
			{
				k = j;
				break;
			}
		}
		if (k == -1)
			continue;
		sourceDesityValue = mocMesh.GetValueAtIndex(vSMocIndex, ValueType::DENSITY);
		sourceTemperatureValue = mocMesh.GetValueAtIndex(vSMocIndex, ValueType::TEMPERAURE);
		pointVolume = mocMesh.GetVolumeAtIndex(vSMocIndex);
		integrationDesityValue[k] += sourceDesityValue * pointVolume;
		integrationTemperatureValue[k] += sourceTemperatureValue * pointVolume;
	}

	ifstream iFile;
	for (int i = 0; i < materialList.size(); i++)
	{
		double dTotalDensity = 0, dTotalTemperature = 0;
		for (int j = 1; j < g_iNumProcs; j++)
		{
			double dDensity = -1, dTemperature = -1;
			std::string filename = "CFD_VALUE_" + materialList[i] + "_" + std::to_string(j);
			iFile.open(filename);
			iFile >> dDensity >> dTemperature;
			std::cout << "desity " << dDensity << "  temperature " << dTemperature
				<< "from process: " << i << "in file " << filename << std::endl;
			dTotalDensity += dDensity;
			dTotalTemperature += dTemperature;
			iFile.close();
		}
		Logger::LogInfo(FormatStr("Moc total Integration of Desity for Material:%s is %.6lf", materialList[i].c_str(), integrationDesityValue[i]));
		Logger::LogInfo(FormatStr("CFD total Integration of Desity for Material:%s is %.6lf", materialList[i].c_str(), dTotalDensity));
		Logger::LogInfo(FormatStr("Moc total Integration of Temperature for Material:%s is %.6lf", materialList[i].c_str(), integrationTemperatureValue[i]));
		Logger::LogInfo(FormatStr("CFD total Integration of Temperature for Material:%s is %.6lf", materialList[i].c_str(), dTotalTemperature));
	}
}
void MOCTOCFD_ValueValidation(MOCMesh& mocMesh, std::vector<std::string>& materialList)
{
	std::vector<Scalar> integrationHeatpower;
	double sourceHeatpower = 0, pointVolume = 0;
	integrationHeatpower.resize(materialList.size());
	SMocIndex vSMocIndex;
	for (int i = 0; i < mocMesh.m_vSMocIndex.size(); i++)
	{
		vSMocIndex = mocMesh.m_vSMocIndex[i];
		int k = -1;
		for (int j = 0; j < materialList.size(); j++)
		{
			if (materialList[j] == mocMesh.GetMaterialNameAtIndex(vSMocIndex))
			{
				k = j;
				break;
			}
		}
		if (k == -1)
			continue;
		sourceHeatpower = mocMesh.GetValueAtIndex(vSMocIndex, ValueType::HEATPOWER);
		pointVolume = mocMesh.GetVolumeAtIndex(vSMocIndex);
		integrationHeatpower[k] += sourceHeatpower * pointVolume;
	}

	ifstream iFile;
	for (int i = 0; i < materialList.size(); i++)
	{
		double dTotalHeatPower = 0;
		for (int j = 1; j < g_iNumProcs; j++)
		{
			double dHeatPower = -1;
			std::string filename = "CFD_HeatPower_" + materialList[i] + "_" + std::to_string(j);
			iFile.open(filename);
			iFile >> dHeatPower;
			std::cout << "HeatPower " << dHeatPower 
				<< "from process: " << i << "in file " << filename << std::endl;
			dTotalHeatPower += dHeatPower;
			iFile.close();
		}
		Logger::LogInfo(FormatStr("Moc total Integration of HeatPower for Material:%s is %.6lf", materialList[i].c_str(), integrationHeatpower[i]));
		Logger::LogInfo(FormatStr("CFD total Integration of HeatPower for Material:%s is %.6lf", materialList[i].c_str(), dTotalHeatPower));
	}
}




void Test()
{
	ifstream iFile("3_3cells.apl");
	ofstream oFile("3_3cells_bin.apl",std::ios_base::binary);
	std::string line;
	while (getline(iFile, line))
	{
		oFile << line << std::endl;
	}
	iFile.close();
	oFile.close();
	int iFileSize = File_size("3_3cells.apl");

	char* cFile = new char[iFileSize+1];
	MPI_File fh;
	MPI_Status status;
	MPI_File_open(MPI_COMM_WORLD, "3_3cells.apl", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
	MPI_File_read_at_all(fh, 0, cFile, iFileSize, MPI_CHAR, &status);
	cFile[iFileSize] = '\0';
	stringstream ss;
	ss << cFile;
	oFile.open("3_3cells_bin_temp.apl");
	getline(ss, line);
	oFile << ss.str();
	oFile.close();
	MPI_File_close(&fh);
	delete[] cFile;

}
void CaculateMinCoordinate()
{
	// recvied all min Coordinate from all cfd vtk file,caculuat the min x y z and send back
	std::vector<Vector> vPoint;
	double xMin = 100000, yMin = 100000, zMin = 100000;
	for (int iSourceID = 1; iSourceID < g_iNumProcs; iSourceID++)
	{
		std::vector<double> vCoordinate;
		vCoordinate.resize(3);
		//std::cout << "MPI_Recv from process: " << g_iMpiID << std::endl;
		MPI_Recv(&vCoordinate[0], 3, MPI_DOUBLE, iSourceID, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		vPoint.push_back(Vector(vCoordinate[0], vCoordinate[1], vCoordinate[2]));
		std::cout << "receive CFD MIN coordinate:" << Vector(vCoordinate[0], vCoordinate[1], vCoordinate[2]) << "from proc,ess:" << iSourceID << std::endl;
		xMin = min(xMin, vCoordinate[0]);
		yMin = min(yMin, vCoordinate[1]);
		zMin = min(zMin, vCoordinate[2]);
	}
	std::vector<double> vCoordinate;
	vCoordinate.push_back(xMin);
	vCoordinate.push_back(yMin);
	vCoordinate.push_back(zMin);
	MPI_Bcast(&vCoordinate[0], 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	std::cout << "MPI_Bcast send " << Vector(vCoordinate[0], vCoordinate[1], vCoordinate[2])
		<< "from process : " << g_iMpiID << std::endl;

}

void ReceiveAndSetMocValue(MOCMesh &mocMesh)
{
	MPI_Datatype mpiMocField_type;
	InitMocFieldToMpiType(mpiMocField_type);
	std::vector<STRMocField> vReciveField;
	for (int iSourceID = 1; iSourceID < g_iNumProcs; iSourceID++)
	{
		int iSize;
		MPI_Recv(&iSize, 1, MPI_INT, iSourceID, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		std::cout << "Receive isize:" << iSize << "from process:" << iSourceID << std::endl;
		vReciveField.resize(iSize);
		MPI_Recv(&vReciveField[0], iSize, mpiMocField_type, iSourceID, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		//printf("MPI process %d received person:\n\t- iAssembly = %d\n\t- icell = %d\n\t- name = %s\n", my_rank, vReciveField[1].iAssemblyIndex,
			//vReciveField[1].iCellIndex, vReciveField[1].cMaterialName);
		std::cout << "Receive data from process: " << iSourceID << std::endl;
		mocMesh.SetFieldByMpiType(vReciveField);
		vReciveField.clear();
	}
	MPI_Type_free(&mpiMocField_type);
}
void CheckMocMappingWeights(MOCMesh& mocMesh, std::vector<std::string>& materialList)
{
	std::vector < std::vector < std::vector < std::unordered_map < std::pair<int, int>, double, pair_hash>>>> m_MOC_CFD_MapWithID;
	int iMax_iAssembly, iMax_iCell, iMax_iMoc;
	std::tuple<int, int, int>tupIndex = GetMaxIndexOfMoc();
	iMax_iAssembly = std::get<0>(tupIndex);
	iMax_iCell = std::get<1>(tupIndex);
	iMax_iMoc = std::get<2>(tupIndex);

	m_MOC_CFD_MapWithID.resize(iMax_iAssembly);
	for (int i = 0; i < iMax_iAssembly; i++)
	{
		m_MOC_CFD_MapWithID[i].resize(iMax_iCell);
		for (int j = 0; j < iMax_iCell; j++)
		{
			m_MOC_CFD_MapWithID[i][j].resize(iMax_iMoc);
		}
	}
	std::string strFileName,line;
	fstream iFile;
	for (int i = 0; i < materialList.size(); i++)
	{

		strFileName = "MapFile_" + materialList[i] + "_MOCtoCFD";
		iFile.open(strFileName);
		Logger::LogInfo("reading MOC to CFD map file in material: " + materialList[i]);
		while (getline(iFile, line))
		{
			int i, j, k, m, n;
			double dValue;
			stringstream stringline(line);
			stringline >> i >> j >> k >> m >> n >> dValue;
			//std::cout << i << " " << j << " " << k << " " << m << " " << n << " " << dValue << std::endl;
			m_MOC_CFD_MapWithID[i][j][k].emplace(std::make_pair(m, n), dValue);
			//for debug
			//if(i==8&&j==0&&k==2391)
				//Logger::LogInfo(FormatStr("8 0 2391 MOC dValue=%.6lf", dValue));
		}
		iFile.close();
	}
	double minSumWeight = 1;
	double maxSumWeight = 0;
	int sumOfZero = 0;
	for (int j = 0; j < mocMesh.m_vSMocIndex.size(); j++)
	{
		double sumSumWeight = 0.0;
		SMocIndex smocTemp= mocMesh.m_vSMocIndex[j];
		for (auto& iter : m_MOC_CFD_MapWithID[smocTemp.iAssemblyIndex][smocTemp.iCellIndex][smocTemp.iMocIndex])
		{
			sumSumWeight += iter.second;
		}
		//for debug
		//if (smocTemp == SMocIndex(8, 0, 2391))
			//Logger::LogInfo(FormatStr("8 0 2391 MOC sumSumWeight=%.6lf", sumSumWeight));

		minSumWeight = min(minSumWeight, sumSumWeight);
		maxSumWeight = max(maxSumWeight, sumSumWeight);
		if (sumSumWeight < SMALL)
		{
			std::stringstream streamTemp;
			streamTemp << "0 sumSumWeight moc id: " <<mocMesh.GetMaterialNameAtIndex(smocTemp)<<" " << smocTemp.iAssemblyIndex
				<< " " << smocTemp.iCellIndex << " " << smocTemp.iMocIndex << std::endl;
			Logger::LogInfo(streamTemp.str(), true);
			sumOfZero++;
		}
	}
	Logger::LogInfo(FormatStr("sum of CFD->MOC weights ranges from %.6lf to %6lf", minSumWeight, maxSumWeight));
	Logger::LogInfo(FormatStr("%d MOC cells have no weights from CFD cells", sumOfZero));
}
int main(int argc, char** argv)
{
	MPI_Status status;
	char message[100];
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &g_iMpiID);
	MPI_Comm_size(MPI_COMM_WORLD, &g_iNumProcs);
	//for test
	//addbytes();
	//for test
	std::vector<std::string> parameterList;
	if (argc > 1)
	{
		parameterList.resize(argc - 1);
		for (int i = 1; i < argc; i++)
			parameterList[i - 1] = argv[i];
	}
	else
	{
		if(g_iMpiID==0)
			DisplayHelpInfo();
		MPI_Finalize();
		return 0;
	}

	if (g_iNumProcs == 1)
	{
		if (parameterList[0] != "createmapper" && parameterList[0] != "cfdtomoc" && parameterList[0] != "moctocfd")
		{
			RunWithParameters(parameterList);
		}
		else
		{
			Logger::LogError("Wrong input parameter");
			DisplayHelpInfo();
		}
		MPI_Finalize();
		return 0;
	}

	//sub process
	if (g_iMpiID > 0) {
		std::cout << "this is child,pid=" << getpid() << std::endl;
		//for debug
		int num = 10;
		while (num == 11)
		{
		#ifdef _WIN32
			Sleep(10);
		#else
			sleep(10);
		#endif
		}
		std::string strMessage;
		if (parameterList[0] == "createmapper"&&argc==2)
		{
			CreateMapper();
			//notify to primary process
			strMessage = "Process " + std::to_string(g_iMpiID) + " createmapper finished";
			strcpy(message, strMessage.data());// "caclulation finished!");
			MPI_Send(message, strlen(message) + 1, MPI_CHAR, 0, 99,MPI_COMM_WORLD);
		}
		else if (parameterList[0] == "moctocfd" && argc == 2)
		{
			MOCFieldsToCFD();
			//notify to primary process
			strMessage = "Process " + std::to_string(g_iMpiID) + " moctocfd finished";
			strcpy(message, strMessage.data());// "caclulation finished!");
			MPI_Send(message, strlen(message) + 1, MPI_CHAR, 0, 99,MPI_COMM_WORLD);
		}
		else if (parameterList[0] == "cfdtomoc" && argc == 2)
		{
			CFDFieldsToMOC();

		}
		else
		{
			Logger::LogError("Wrong input parameter");
			DisplayHelpInfo();
		}
	}
	else if(g_iMpiID ==0&& g_iNumProcs>1)//primary process
	{
		std::cout << "this is 0 Process,pid=" << getpid() << std::endl;
		//for debug
		int num = 10;
		while (num == 11)
		{
		#ifdef _WIN32
			Sleep(10);
		#else
			sleep(10);
		#endif
		}
		std::vector<std::vector<std::string> > matches = GetMatchList(configFile);
		std::vector<std::string>& materialList = matches[0];
		std::string mocMeshFile = GetFileName(configFile, "inputApl");
		std::string outMocMeshFile = GetFileName(configFile, "outputApl");
		std::string mocPowerFile = GetFileName(configFile, "mocPower");
		std::string mocFieldFile = GetFileName(configFile, "inputInp");
		std::string outMocFieldFile = GetFileName(configFile, "outputInp");
		if(argc==2 && parameterList[0]=="createmapper")
		{
			WriteToLog("Createmapper start.");
			std::vector<std::vector<std::string> > matches = GetMatchList(configFile);
			std::vector<std::string>& materialList = matches[0];
			std::string mocMeshFile = GetFileName(configFile, "inputApl");
			stringstream strTemp;
			//MPI_File_open moc date file need all process involed 
			MPI_OpenFile_To_Stream(mocMeshFile, strTemp);

			// recvied all min Coordinate from all cfd vtk file,caculuat the min x y z and send back
			CaculateMinCoordinate();

			//solution 2
			stringstream  MOCtoCFD_MapFile_stream;
			MOCtoCFD_MapFile_stream << "";
			for(int i=0;i< materialList.size();i++)
				MPI_WriteStream_To_File(("MapFile_" + materialList[i] + "_MOCtoCFD"), MOCtoCFD_MapFile_stream);
			//solution end

			for (int iSourceID = 1; iSourceID < g_iNumProcs; iSourceID++) {
				MPI_Recv(message, 100, MPI_CHAR, iSourceID, 99, MPI_COMM_WORLD, &status);
				Logger::LogInfo(FormatStr("Main process received message from No.%d process: %s\n", iSourceID, message));
			}
			MOCMesh mocMesh(mocMeshFile, outMocMeshFile, MeshKernelType::MHT_KERNEL, false);
			CheckMocMappingWeights(mocMesh, materialList);
			//solution 1 
			//ConvergeMocMapInfor();
			Logger::LogInfo("All createmapper finished.");
			WriteToLog("Createmapper end.");
		}
		if (argc == 2 && parameterList[0] == "cfdtomoc")
		{
			WriteToLog("cfdtomoc start.");
			//MPI_File_read_at_all read MapFile_***_MOCtoCFD.   all process need to be  involved
			for (int i = 0; i < materialList.size(); i++)
			{
				stringstream strTemp;
				std::string fileName = "MapFile_" + materialList[i] + "_MOCtoCFD";
				MPI_OpenFile_To_Stream(fileName, strTemp);
			}
			//solution 2
			//FinalCFDFieldsToMOC();
			//solution 1
			MOCMesh mocMesh(mocMeshFile, outMocMeshFile, MeshKernelType::MHT_KERNEL,false);
			mocMesh.InitMOCFromInputFile(mocFieldFile);
			mocMesh.SetDesityAndTemperatureToZero();

			ReceiveAndSetMocValue(mocMesh);
			CFDTOMOC_ValueValidation(mocMesh, materialList);
			
			RenameFile(outMocFieldFile, GetFileNameOfPrevious(outMocFieldFile, "inp"));
			mocMesh.OutputStatus(outMocFieldFile);
			
			Logger::LogInfo("CFD to MOC finished.");
			WriteToLog("cfdtomoc end.");
		}
		else if (argc == 2 && parameterList[0] == "moctocfd")
		{
			WriteToLog("moctocfd start.");
			stringstream ifs;
			//MPI_File_open heatepower date file need all process involed 
			MPI_OpenFile_To_Stream(mocPowerFile, ifs);
			//MPI_File_open MapFile_***_MOCtoCFD file need all process involed 
			for (int i = 0; i < materialList.size(); i++)
			{
				stringstream strTemp;
				std::string fileName = "MapFile_" + materialList[i] + "_MOCtoCFD";
				MPI_OpenFile_To_Stream(fileName, strTemp);
			}
			MOCMesh mocMesh(mocMeshFile, outMocMeshFile, MeshKernelType::MHT_KERNEL,false);
			mocMesh.InitMOCHeatPower(mocPowerFile);

			for (int iSourceID = 1; iSourceID < g_iNumProcs; iSourceID++) {
				MPI_Recv(message, 100, MPI_CHAR, iSourceID, 99, MPI_COMM_WORLD, &status);
				Logger::LogInfo(FormatStr("Main process received message from No.%d process: %s\n", iSourceID, message));
			}
			MOCTOCFD_ValueValidation(mocMesh, materialList);
			Logger::LogInfo(" MOC to CFD finished.");
			WriteToLog("moctocfd end.");
		}
	}
	else
	{
		Logger::LogInfo(" Wrong input parameter");
		DisplayHelpInfo();
	}
	MPI_Finalize();
	return 0;
}
