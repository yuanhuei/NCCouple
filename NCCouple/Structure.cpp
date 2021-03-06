#include "Structure.h"
#include<string>

#include "Logger.h"
#include <future>
#include <mutex>
#include<algorithm>
#include "index.h"
#include "MOCMesh.h"
#include "ConfigurationFile.h"
#include<mpi.h>
#ifdef _WIN32
#include <string>
#else
#include <string>
#include <string.h>
#endif
#include <sys/stat.h>
void WriteToLog(std::string strInfo, int iMpiID)
{
	std::string strLogName = "log" + std::to_string(iMpiID) + ".txt";
	ofstream iFile(strLogName,std::ios::app);
	time_t timep;
	time(&timep);
	char* cTime = ctime(&timep);
	cTime[strlen(cTime) - 2] = '\0';
	iFile << "[" << cTime << "]" << strInfo << std::endl;
	iFile.close();
	return;
}
void ConvergeMocMapInfor()
{


	int iMax_iAssembly, iMax_iCell, iMax_iMoc;
	std::tuple<int, int, int>tupIndex = GetMaxIndexOfMoc();
	iMax_iAssembly = std::get<0>(tupIndex);
	iMax_iCell = std::get<1>(tupIndex);
	iMax_iMoc = std::get<2>(tupIndex);


	std::vector<std::vector<std::vector<std::unordered_map<std::pair<int,int>, double, pair_hash>>>> MOC_CFD_Map;

	MOC_CFD_Map.resize(iMax_iAssembly);
	for (int i = 0; i < iMax_iAssembly; i++)
	{
		MOC_CFD_Map[i].resize(iMax_iCell);
		for (int j = 0; j < iMax_iCell; j++)
		{
			MOC_CFD_Map[i][j].resize(iMax_iMoc);
		}
	}
	//configFile = "MapFile_FileNames";
	std::string mocMeshFile = GetFileName(configFile, "inputApl");
	std::string outMocMeshFile = GetFileName(configFile, "outputApl");
	std::string mocPowerFile = GetFileName(configFile, "mocPower");
	std::vector<std::vector<std::string> > matches = GetMatchList(configFile);
	std::vector<std::string>& materialList = matches[0];
	std::string fileName,line;
	ifstream infile;
	for (int j = 0; j < materialList.size(); j++ )
	{
		for (int i = 1; i < g_iNumProcs; i++)
		{
			fileName = "./temp/MapFile_" + materialList[j] + "_MOCtoCFD"+"_"+std::to_string(i);
			infile.open(fileName);
			if (!infile.is_open())
			{
				Logger::LogError("cannot find the MOC to CFD map file:" + fileName);
				exit(EXIT_FAILURE);
				return;
			}
			Logger::LogInfo("reading MOC to CFD map file in material: " + materialList[j]+"in the file "+ fileName);
			while (getline(infile, line))
			{
				int iAssembly , iCell, iMesh, n;
				double dValue;
				stringstream stringline(line);
				stringline >> iAssembly >> iCell >> iMesh >> n >> dValue;
				MOC_CFD_Map[iAssembly][iCell][iMesh].emplace(std::make_pair(i,n),dValue);
			}
			infile.close();

		}
		ofstream MOCtoCFD_MapFile("./temp/MapFile_" + materialList[j] + "_MOCtoCFD");
		for (int i = 0; i < MOC_CFD_Map.size(); i++)
		{
			for (int m = 0; m < MOC_CFD_Map[i].size(); m++)
			{
				for (int k = 0; k < MOC_CFD_Map[i][m].size(); k++)
				{
					//std::unordered_map<std::pair<int,int>, double, pair_hash>::iterator it;
					for (auto it = MOC_CFD_Map[i][m][k].begin(); it != MOC_CFD_Map[i][m][k].end();it++)
					{
							MOCtoCFD_MapFile << i << " " << m << " " << k << " " << it->first.first
								<< " " << it->first.second << " " << it->second << std::endl;

					}
					if(!MOC_CFD_Map[i][m][k].empty())
						MOC_CFD_Map[i][m][k].clear();
				}
			}
		}
		MOCtoCFD_MapFile.close();
	}
	Logger::LogInfo("Writing Moc finished." );
	return;
}


void SendMocValueToMainProcess(const std::vector<STRMocField>& vMocField)
{
	//STRMocField stField;
	MPI_Datatype STRMocFieldMPIType;
	InitMocFieldToMpiType(STRMocFieldMPIType);


	int iSize = vMocField.size();
	MPI_Send(&iSize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	MPI_Send(&vMocField[0], iSize, STRMocFieldMPIType, 0, 1, MPI_COMM_WORLD);

	MPI_Type_free(&STRMocFieldMPIType);

	return;
}

void InitMocFieldToMpiType(MPI_Datatype &mpiMocField_type)
{
	STRMocField stField;
	// Create the datatype

	int lengths[7] = { 1, 1,1,1,1,20,20 };
	MPI_Aint displacements[7];
	MPI_Aint base_address;
	MPI_Get_address(&stField, &base_address);
	MPI_Get_address(&stField.iAssemblyIndex, &displacements[0]);
	MPI_Get_address(&stField.iCellIndex, &displacements[1]);
	MPI_Get_address(&stField.iMeshIndex, &displacements[2]);
	MPI_Get_address(&stField.dTempValue, &displacements[3]);
	MPI_Get_address(&stField.dDensityValue, &displacements[4]);
	MPI_Get_address(&stField.cMaterialName[0], &displacements[5]);
	MPI_Get_address(&stField.cTempName[0], &displacements[6]);

	displacements[0] = MPI_Aint_diff(displacements[0], base_address);
	displacements[1] = MPI_Aint_diff(displacements[1], base_address);
	displacements[2] = MPI_Aint_diff(displacements[2], base_address);
	displacements[3] = MPI_Aint_diff(displacements[3], base_address);
	displacements[4] = MPI_Aint_diff(displacements[4], base_address);
	displacements[5] = MPI_Aint_diff(displacements[5], base_address);
	displacements[6] = MPI_Aint_diff(displacements[6], base_address);

	MPI_Datatype types[7] = { MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE,MPI_DOUBLE, MPI_CHAR,MPI_CHAR };
	MPI_Type_create_struct(7, lengths, displacements, types, &mpiMocField_type);
	MPI_Type_commit(&mpiMocField_type);
}
void SendMocMapValueToMainProcess(const std::vector<STRMocMapValue>& vMocMapValue)
{
	//STRMocField stField;
	MPI_Datatype STRMocFieldMPIType;
	InitMocMapValueToMpiType(STRMocFieldMPIType);


	int iSize = vMocMapValue.size();
	MPI_Send(&iSize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	MPI_Send(&vMocMapValue[0], iSize, STRMocFieldMPIType, 0, 1, MPI_COMM_WORLD);

	MPI_Type_free(&STRMocFieldMPIType);

	return;
}
void InitMocMapValueToMpiType(MPI_Datatype& mpiMocMapValue_type)
{
	STRMocMapValue stField;
	// Create the datatype

	int lengths[4] = { 1, 1,1,1 };
	MPI_Aint displacements[4];
	MPI_Aint base_address;
	MPI_Get_address(&stField, &base_address);
	MPI_Get_address(&stField.iAssemblyIndex, &displacements[0]);
	MPI_Get_address(&stField.iCellIndex, &displacements[1]);
	MPI_Get_address(&stField.iMeshIndex, &displacements[2]);
	MPI_Get_address(&stField.dMapValue, &displacements[3]);

	displacements[0] = MPI_Aint_diff(displacements[0], base_address);
	displacements[1] = MPI_Aint_diff(displacements[1], base_address);
	displacements[2] = MPI_Aint_diff(displacements[2], base_address);
	displacements[3] = MPI_Aint_diff(displacements[3], base_address);

	MPI_Datatype types[4] = { MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE };
	MPI_Type_create_struct(4, lengths, displacements, types, &mpiMocMapValue_type);
	MPI_Type_commit(&mpiMocMapValue_type);
}

Vector UpdateMinLocation(Vector vPoint)
{
	std::cout << "UpdateMinLocation called by:" << g_iMpiID << std::endl;
	std::vector<double> vCoordinate;
	vCoordinate.push_back(vPoint.x_);
	vCoordinate.push_back(vPoint.y_);
	vCoordinate.push_back(vPoint.z_);
	
	MPI_Send(&vCoordinate[0], 3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	std::cout << "MPI_Send " << vPoint << "from process : " << g_iMpiID << std::endl;
	std::vector<double> vCoordinate_new;
	vCoordinate_new.resize(3);
	//MPI_Recv(&vCoordinate_new[0], 3, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Bcast(&vCoordinate_new[0], 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	std::cout << "MPI_Bcast receive " << Vector(vCoordinate_new[0], vCoordinate_new[1], vCoordinate_new[2])
		<< "from process : " << g_iMpiID << std::endl;
	return Vector(vCoordinate_new[0], vCoordinate_new[1], vCoordinate_new[2]);
}

int File_size(const char* filename)//get size(byte) of file

{
	struct stat statbuf;
	int ret;
	ret = stat(filename, &statbuf);
	if (ret != 0) return -1;
	return statbuf.st_size;
}

void MPI_OpenFile_To_Stream(std::string strFileName,stringstream& strTempStream)
{

	int iFileSize = File_size(strFileName.c_str());
	char* cFile = new char[iFileSize + 1];
	MPI_File fh;
	MPI_Status status;
	MPI_File_open(MPI_COMM_WORLD, strFileName.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
	MPI_File_read_at_all(fh, 0, cFile, iFileSize, MPI_CHAR, &status);
	cFile[iFileSize] = '\0';
	strTempStream << cFile;
	MPI_File_close(&fh);
	delete[] cFile;
}
void MPI_WriteStream_To_File(std::string strFileName, stringstream& strTempStream)
{

	std::string strTemp = strTempStream.str().c_str();
	MPI_File fh;
	MPI_Status status;
	MPI_File_open(MPI_COMM_WORLD, strFileName.c_str(),MPI_MODE_CREATE+MPI_MODE_RDWR, MPI_INFO_NULL, &fh);
	MPI_File_write_ordered(fh,strTemp.c_str(), strTemp.length(), MPI_CHAR, &status);
	MPI_File_close(&fh);
}

std::tuple<int, int, int>GetMaxIndexOfMoc()
{
	std::string fileName;
	if(g_iMpiID==0)
		fileName = "./temp/MocMeshMaxSize_info_"+std::to_string(1);
	else
		fileName = "./temp/MocMeshMaxSize_info_" + std::to_string(g_iMpiID);

	ifstream infile_MocMeshMaxSize_info(fileName);
	if (!infile_MocMeshMaxSize_info.is_open())
	{
		Logger::LogError("cannot find the MocMeshMaxSize_info file:" + fileName);
		exit(EXIT_FAILURE);
	}
	int iAssembly = -1, iCell = -1, iMesh = -1;
	infile_MocMeshMaxSize_info >> iAssembly >> iCell >> iMesh;
	infile_MocMeshMaxSize_info.close();
	return std::make_tuple(iAssembly, iCell, iMesh);
}