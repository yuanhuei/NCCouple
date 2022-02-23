#ifndef _ConfigurationFile_
#define _ConfigurationFile_

#include <vector>
#include <string>

std::string GetFileNameOfPrevious
(
	std::string fileName,
	std::string fileType
);

void RemoveFile
(
	std::string nameStr
);

void RenameFile
(
	std::string oldname,
	std::string newname
);

void WriteConfigurationFile
(
	std::string configFile,
	//std::string& mshFile,
	std::string& aplFile,
	std::string& outAplFile,
	std::vector<std::string>& MOCMaterials
	//std::vector<std::string>& CFDRegions
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

std::vector<std::vector<std::string> > GetMatchList
(
	std::string configFile
);

#endif