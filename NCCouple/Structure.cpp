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

void ConvergeMocMapInfor()
{
	//m_CFD_MOC_Map.resize(m_cfdMeshPtr->GetMeshPointNum());
	std::string fileName = "MocMeshMaxSize_info";
	ifstream infile(fileName);
	std::string line;
	//getline(infile, line);
	//stringstream stringline(line);
	int iAssemby, iCell, iMesh;
	infile >> iAssemby >> iCell >> iMesh ;
	infile.close();
	std::vector<std::vector<std::vector<std::unordered_map<std::pair<int,int>, double, pair_hash>>>> MOC_CFD_Map;

	MOC_CFD_Map.resize(iAssemby);
	for (int i = 0; i <iAssemby; i++)
	{
		MOC_CFD_Map[i].resize(iCell);
		for (int j = 0; j < iCell; j++)
		{
			MOC_CFD_Map[i][j].resize(iMesh);
		}
	}
	//configFile = "MapFile_FileNames";
	std::string mocMeshFile = GetFileName(configFile, "inputApl");
	std::string outMocMeshFile = GetFileName(configFile, "outputApl");
	std::string mocPowerFile = GetFileName(configFile, "mocPower");
	std::vector<std::vector<std::string> > matches = GetMatchList(configFile);
	std::vector<std::string>& materialList = matches[0];

	for (int j = 0; j < materialList.size(); j++ )
	{
		for (int i = 1; i < g_iNumProcs; i++)
		{
			fileName = "MapFile_" + materialList[j] + "_MOCtoCFD"+"_"+std::to_string(i);
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
		ofstream MOCtoCFD_MapFile("MapFile_" + materialList[j] + "_MOCtoCFD");
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
	std::cout << "Writing Moc finished" << std::endl;
	return;
}


void SendFieldForTest(const std::vector<STRMocField>& vMocField)
{
	//STRMocField stField;
	MPI_Datatype STRMocFieldMPIType;
	InitMocFieldToMpiType(STRMocFieldMPIType);

	/*
	std::vector<STRMocField> vField;
	struct STRMocField buffer;
	buffer.iAssemblyIndex = 20;
	buffer.iCellIndex = 20;
	buffer.iMeshIndex = 20;
	buffer.dDensityValue = 20;
	buffer.dTempValue = 20;
	strncpy(buffer.cMaterialName, "Tom", 19);
	buffer.cMaterialName[19] = '\0';
	strncpy(buffer.cTempName, "Jone", 19);
	buffer.cTempName[19] = '\0';

	vField.push_back(buffer);
	buffer.iAssemblyIndex = 40;
	vField.push_back(buffer);
	printf("MPI process  sends messge\n");
	*/
	int iSize = vMocField.size();
	MPI_Send(&iSize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	MPI_Send(&vMocField[0], iSize, STRMocFieldMPIType, 0, 1, MPI_COMM_WORLD);

	MPI_Type_free(&STRMocFieldMPIType);

	return;
}

void InitMocFieldToMpiType(MPI_Datatype &mpiMocField_type)
{
	STRMocField stField;
	/*
	struct STRMocField
	{
		int iAssemblyIndex, iCellIndex, iMeshIndex;
		double dTempValue, dDensityValue;
		char cMaterialName[20], cTempName[20];
	};*/
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

