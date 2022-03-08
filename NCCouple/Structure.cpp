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
						if (it->second != -1)
						{
							MOCtoCFD_MapFile << i << " " << m << " " << k << " " << it->first.first
								<< " " << it->first.second << " " << it->second << std::endl;
							it->second = -1;
						}

						//MOC_CFD_Map.erase(*itit->second
						//it=MOC_CFD_Map[i][m][k].erase(it);
					}
					//if(!MOC_CFD_Map[i][m][k].empty())
						//MOC_CFD_Map[i][m][k].clear();
				}
			}
		}
		MOCtoCFD_MapFile.close();
	}
	std::cout << "Writing Moc finished" << std::endl;
	return;
}

void SetFieldValue(
	const std::vector<std::vector<std::vector<std::shared_ptr<MocMeshField>>>> &vAssemlbyField,
	std::vector<STRMocField>& vMocField)
{
	STRMocField  stMocTemp;
	for (int i = 0; i < vAssemlbyField.size(); i++)
	{
		for (int j = 0; j < vAssemlbyField[i].size(); j++)
		{
			for (int k = 0; k < vAssemlbyField[i][j].size(); k++)
			{
				if (vAssemlbyField[i][j][k])
				{
					stMocTemp.iAssemblyIndex = i;
					stMocTemp.iCellIndex = j;
					stMocTemp.iMeshIndex = k;
					stMocTemp.dDensityValue = vAssemlbyField[i][j][k]->GetValue(ValueType::DENSITY);
					stMocTemp.dTempValue= vAssemlbyField[i][j][k]->GetValue(ValueType::TEMPERAURE);
					strcpy(stMocTemp.cMaterialName , vAssemlbyField[i][j][k]->GetMaterialNameWithID().c_str());
					strcpy(stMocTemp.cTempName, vAssemlbyField[i][j][k]->GetTemperatureName().c_str());

				}

			}
		}

	}

}
void SendFieldForTest(int argc, char** argv)
{

	MPI_Init(&argc, &argv);

	// Get the number of processes and check only 2 processes are used
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	if (size != 2)
	{
		printf("This application is meant to be run with 2 processes.\n");
		MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
	}

	STRMocField stField;

	struct STRMocField
	{
		int iAssemblyIndex, iCellIndex, iMeshIndex;
		double dTempValue, dDensityValue;
		char cMaterialName[20], cTempName[20];
	};
	// Create the datatype
	MPI_Datatype person_type;
	int lengths[7] = { 1, 1,1,1,1,20,20 };

	// Calculate displacements
	// In C, by default padding can be inserted between fields. MPI_Get_address will allow
	// to get the address of each struct field and calculate the corresponding displacement
	// relative to that struct base address. The displacements thus calculated will therefore
	// include padding if any.
	MPI_Aint displacements[7];
	struct STRMocField dummy_person;
	MPI_Aint base_address;
	MPI_Get_address(&dummy_person, &base_address);
	MPI_Get_address(&dummy_person.iAssemblyIndex, &displacements[0]);
	MPI_Get_address(&dummy_person.iCellIndex, &displacements[1]);
	MPI_Get_address(&dummy_person.iMeshIndex, &displacements[2]);
	MPI_Get_address(&dummy_person.dTempValue, &displacements[3]);
	MPI_Get_address(&dummy_person.dDensityValue, &displacements[4]);
	MPI_Get_address(&dummy_person.cMaterialName[0], &displacements[5]);
	MPI_Get_address(&dummy_person.cTempName[0], &displacements[6]);

	displacements[0] = MPI_Aint_diff(displacements[0], base_address);
	displacements[1] = MPI_Aint_diff(displacements[1], base_address);
	displacements[2] = MPI_Aint_diff(displacements[2], base_address);
	displacements[3] = MPI_Aint_diff(displacements[3], base_address);
	displacements[4] = MPI_Aint_diff(displacements[4], base_address);
	displacements[5] = MPI_Aint_diff(displacements[5], base_address);
	displacements[6] = MPI_Aint_diff(displacements[6], base_address);

	MPI_Datatype types[7] = { MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE,MPI_DOUBLE, MPI_CHAR,MPI_CHAR };
	MPI_Type_create_struct(7, lengths, displacements, types, &person_type);
	MPI_Type_commit(&person_type);

	int blockLength[] = { 1,1,4,1 };
	MPI_Datatype oldTypes[] = { MPI_INT,MPI_DOUBLE,MPI_DOUBLE,MPI_FLOAT };

	//直接计算偏移量,每一个数据以double对齐  
	MPI_Aint addressOffsets[] = { 0,sizeof(double),2 * sizeof(double),6 * sizeof(double) };


	// Get my rank and do the corresponding job
	enum rank_roles { SENDER, RECEIVER };
	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	switch (my_rank)
	{
	case SENDER:
	{
		// Send the message
		struct STRMocField buffer;
		buffer.iAssemblyIndex = 20;
		buffer.iCellIndex = 20;
		buffer.iMeshIndex = 20;
		buffer.dDensityValue = 20;
		buffer.dTempValue = 20;

		strncpy(buffer.cMaterialName, "Tom", 19);
		buffer.cMaterialName[19] = '\0';
		strncpy(buffer.cTempName, "Jone", 19);
		buffer.cMaterialName[19] = '\0';

		printf("MPI process %d sends person");
		MPI_Send(&buffer, 1, person_type, RECEIVER, 0, MPI_COMM_WORLD);
		break;
	}
	case RECEIVER:
	{
		// Receive the message
		struct STRMocField received;
		MPI_Recv(&received, 1, person_type, SENDER, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		printf("MPI process %d received person:\n\t- iAssembly = %d\n\t- icell = %d\n\t- name = %s\n", my_rank, received.iAssemblyIndex,
			received.iCellIndex, received.cMaterialName);
		break;
	}
	}

	MPI_Finalize();

	return;
}

