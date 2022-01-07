#include "./ConfigurationFile.h"
#include <sstream>
#include <iostream>
#include <string>
#include "./MHT_common/StringTools.h"
#include "./Logger.h"

void WriteConfigurationFile
(
	std::string configFile,
	std::string& mshFile,
	std::string& aplFile,
	std::string& outAplFile,
	std::vector<int>& IDList,
	std::vector<std::string>& regionNameList
)
{
	if (IDList.size() != regionNameList.size())
	{
		Logger::LogError("in WriteConfigurationFile(), numbers of region IDs and region names are not consistent");
	}
	std::ofstream config(configFile);
	//write msh file name
	config << "------inputMsh------" << std::endl;
	config << "(" << mshFile << ")" << std::endl << std::endl;
	config << "------regionID------" << std::endl;
	config << "(" << std::endl;
	for (size_t i = 0;i < IDList.size();i++)
	{
		config << "\t" << IDList[i] << std::endl;
	}
	config << ")" << std::endl << std::endl;
	config << "------regionName------" << std::endl;
	config << "(" << std::endl;
	for (size_t i = 0;i < IDList.size();i++)
	{
		config << regionNameList[i] << std::endl;
	}
	config << ")" << std::endl << std::endl;
	//write apl file infomation
	config << "------inputApl------" << std::endl;
	config << "("<< aplFile << ")" << std::endl << std::endl;
	config << "------outputApl------" << std::endl;
	config << "(" << outAplFile << ")" << std::endl << std::endl;
	//write inp file information
	config << "------inputInp------" << std::endl;
	config << "(***.inp)" << std::endl << std::endl;
	config << "------outputInp------" << std::endl;
	config << "(***.inp)" << std::endl << std::endl;
	//write moc heat power file information
	config << "------mocPower------" << std::endl;
	config << "(***.txt)" << std::endl << std::endl;
	//write vtk
	config << "------inputVtk------" << std::endl;
	config << "(" << std::endl;
	for (size_t i = 0;i < IDList.size();i++)
	{
		config << "\t" << "***.vtk" << std::endl;
	}
	config << ")" << std::endl << std::endl;
	config << "------outputVtk------" << std::endl;
	config << "(" << std::endl;
	for (size_t i = 0;i < IDList.size();i++)
	{
		config << "\t" << "***.vtk" << std::endl;
	}
	config << ")" << std::endl << std::endl;
	config.close();
	return;
}

std::string GetFileName
(
	std::string configFile,
	std::string symbol
)
{
	std::ifstream file(configFile);
	if (!file.is_open())
	{
		Logger::LogError("configuration file MapFile_FileNames is not existing");
	}
	bool found = false;
	std::string result;
	while (!file.eof())
	{
		std::getline(file, result);
		if (std::string::npos != result.find(symbol))
		{
			found = true;
			break;
		}
	}
	if (found == true)
	{
		std::stringstream ss;
		GetBySymbol(file, ss, '(', ')');
		file.close();
		ss >> result;
	}
	else
	{
		Logger::LogError("in GetFileName(), symbol " + symbol + " is not found");
	}
	return result;
}

std::vector<std::string> GetFileNameList
(
	std::string configFile,
	std::string symbol
)
{
	std::vector<std::string> nameList;
	std::ifstream file(configFile);
	if (!file.is_open())
	{
		Logger::LogError("configuration file MapFile_FileNames is not existing");
	}
	bool found = false;
	std::string result;
	while (!file.eof())
	{
		std::getline(file, result);
		if (std::string::npos != result.find(symbol))
		{
			found = true;
			break;
		}
	}
	if (found == true)
	{
		std::stringstream ss;
		GetBySymbol(file, ss, '(', ')');
		while (!ss.eof())
		{
			std::string oneFileName;
			ss >> oneFileName;
			if (!oneFileName.empty())
			{
				nameList.push_back(oneFileName);
			}
		}
		file.close();
	}
	else
	{
		file.close();
		Logger::LogError("in GetFileNameList(), symbol " + symbol + " is not found");
	}
	return nameList;
}

std::vector<int> GetRegionIDList
(
	std::string configFile
)
{
	std::vector<int> IDList;
	std::vector<std::string> nameList;
	nameList = GetFileNameList(configFile, "regionID");
	for (size_t i = 0;i < nameList.size();i++)
	{
		IDList.push_back(std::stoi(nameList[i]));
	}
	return IDList;
}