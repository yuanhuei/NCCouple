#include "./ConfigurationFile.h"
#include <sstream>
#include <iostream>
#include <string>
#include <stdio.h>
#include <fcntl.h>
#include "./MHT_common/StringTools.h"
#include "./Logger.h"

std::string GetFileNameOfPrevious
(
	//original file name
	std::string fileName,
	//file type name (such as txt vtk apl inp)
	std::string fileType
)
{
	int nLength = fileName.size();
	std::string fileTypeWithPoint = "." + fileType;
	int nTypeName = fileTypeWithPoint.size();
	//check the given file name
	if (nLength <= nTypeName)
	{
		Logger::LogError(fileName + " is not a valid " + fileType + " file");
	}
	int clipLocation = nLength - nTypeName;
	std::string fileTypePart = fileName.substr(clipLocation, nLength);
	if (fileTypePart != fileTypeWithPoint)
	{
		Logger::LogError(fileName + " is not a valid " + fileType + " file");
	}
	//add _0 to the end of the file and append with file type name
	std::string extenedName = fileName.substr(0, clipLocation) + "_0" + fileTypeWithPoint;
	return extenedName;
}

void RemoveFile
(
	std::string nameStr
)
{
	const char* nameChar = nameStr.data();
	//delete if the new-name file is existing
	std::ifstream file(nameStr);
	if (file.is_open())
	{
		file.close();
		if (0 == remove(nameChar))
		{
			Logger::LogInfo(nameStr + " is removed");
		}
	}
	return;
}

void RenameFile
(
	std::string oldnameStr,
	std::string newnameStr
)
{
	const char* oldname = oldnameStr.data();
	const char* newname = newnameStr.data();
	//delete if the new-name file is existing
	std::ifstream file(newnameStr);
	if (file.is_open())
	{
		file.close();
		if (0 == remove(newname))
		{
			Logger::LogInfo(newnameStr + " is removed");
		}
	}
	//make sure the old file is existing
	file.open(oldnameStr);
	if (file.is_open())
	{
		file.close();
		if (0 == rename(oldname, newname))
		{
			Logger::LogInfo(oldnameStr + " is renamed as " + newnameStr);
		}
	}
	return;
}

void WriteConfigurationFile
(
	std::string configFile,
	std::string& aplFile,
	std::string& outAplFile,
	std::vector<std::string>& MOCMaterials
)
{
	std::ofstream config(configFile);
	config << "******MOC materials******" << std::endl;
	config << "(" << std::endl;
	for (size_t i = 0;i < MOCMaterials.size();i++)
	{
		config << MOCMaterials[i] << std::endl;
	}
	config << ")" << std::endl << std::endl;

	//write material-region match list
	config << "------scaleRatio------" << std::endl;
	config << "(" << std::endl;
	config << "\t1.0" << std::endl;
	config << ")" << std::endl << std::endl;

	//write material-region match list
	config << "------matchList------" << std::endl;
	config << "(" << std::endl;
	for (size_t i = 0;i < MOCMaterials.size();i++)
	{
		config << "(" << std::endl;
		config << "\t" << MOCMaterials[i] << std::endl;
		config << "\t[write input file].vtk" << std::endl;
		config << "\t[write output file].vtk" << std::endl;
		config << ")" << std::endl << std::endl;
	}
	config << ")" << std::endl << std::endl;
	//write apl file infomation
	config << "------inputApl------" << std::endl;
	config << "(" << aplFile << ")" << std::endl << std::endl;
	config << "------outputApl------" << std::endl;
	config << "(" << outAplFile << ")" << std::endl << std::endl;
	//write inp file infomation
	config << "------inputInp------" << std::endl;
	config << "([write MOC input data file].inp)" << std::endl << std::endl;
	config << "------outputInp------" << std::endl;
	config << "([write MOC output data file].inp)" << std::endl << std::endl;
	//write moc heat power file information
	config << "------mocPower------" << std::endl;
	config << "([write MOC heat power file].txt)" << std::endl << std::endl;
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

double GetValueInFile
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
	std::string oneLine;
	double value = 0.0;
	while (!file.eof())
	{
		std::getline(file, oneLine);
		if (std::string::npos != oneLine.find(symbol))
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
		ss >> value;
	}
	else
	{
		Logger::LogError("in GetFileName(), symbol " + symbol + " is not found");
	}
	return value;
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

std::vector<std::vector<std::string> > GetMatchList
(
	std::string configFile
)
{
	std::vector<std::vector<std::string> > matchList;
	//1. MOC material, 2. input vtk, 3. output vtk
	matchList.resize(3);
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
		if (std::string::npos != result.find("matchList"))
		{
			found = true;
			break;
		}
	}
	if (found == true)
	{
		std::stringstream fullList;
		//extract all match lists along with brackets
		GetBySymbol(file, fullList, '(', ')');
		while (!fullList.eof())
		{
			std::stringstream oneMatch;
			//extract one match with no brackets
			GetBySymbol(fullList, oneMatch, '(', ')');
			//members of one match shoule be 1. MOC material, 2. input vtk, 3. output vtk
			std::vector<std::string> strList;
			while (!oneMatch.eof())
			{
				std::string memberOfMatch;
				oneMatch >> memberOfMatch;
				if (!memberOfMatch.empty())
				{
					strList.push_back(memberOfMatch);
				}
			}
			if (0 == strList.size()) continue;
			if (3 != strList.size())
			{
				Logger::LogError("incorrect match list: " + oneMatch.str());
			}
			for (int i = 0;i < 3;i++)
			{
				matchList[i].push_back(strList[i]);
			}
		}
		file.close();
	}
	else
	{
		file.close();
		Logger::LogError("in GetFileNameList(), symbol matchList is not found");
	}
	return matchList;
}