#include "Structure.h"
#include<string>

#include "Logger.h"
#include <future>
#include <mutex>
#include<algorithm>
#include "index.h"
#include "MOCMesh.h"
#include "ConfigurationFile.h"


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




