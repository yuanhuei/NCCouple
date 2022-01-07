#ifndef _ConfigurationFile_
#define _ConfigurationFile_

#include <vector>
#include <string>

void WriteConfigurationFile
(
	std::string configFile,
	std::string& mshFile,
	std::string& aplFile,
	std::string& inpFile,
	std::vector<int>& IDList,
	std::vector<std::string>& regionNameList
);

std::string GetFileName
(
	std::string configFile,
	std::string symbol
);

std::vector<std::string> GetFileNameList
(
	std::string configFile,
	std::string symbol
);

std::vector<int> GetRegionIDList
(
	std::string configFile
);

#endif