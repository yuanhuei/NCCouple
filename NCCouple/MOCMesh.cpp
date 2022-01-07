#include "MOCMesh.h"
#include <iostream> 
#include <fstream>
#include <map>
#include <tuple>
#include <regex>
#include "Logger.h"
using namespace std;

#define PI 3.14159265358979323846
#define NA 6.022e23
#define BARN 1.0e-24

MOCMesh::MOCMesh(std::string meshFileName, std::string outAplFileName, MeshKernelType kernelType) 
{
	ofstream outFile(outAplFileName);
	ifstream infile(meshFileName);
	if (!infile.is_open())
	{
		Logger::LogError("cannot find the MOC mesh file:" + meshFileName);
		exit(EXIT_FAILURE);
	}
	string line;
	vector<string> meshMaterialNameTemperary;
	vector<string> meshTemperatureNameTemperary;
	vector<string> meshFaceTypeTemperary;
	vector<string> meshFaceTemperatureNameTemperary;
	vector<string> fileNameTemperary;
	vector<int> meshIDTemperary;
	int nFineMesh = kernelType == MeshKernelType::CGAL_KERNEL ? 4 : 1;  //fine mesh
	std::vector<Surface>allMeshFaces;   //all face objects
	std::vector<MOCEdge>allEdges;    //all edge objects
	while (getline(infile, line))  //read mesh data
	{
		outFile << line << endl;
		stringstream stringline(line);
		string token;
		while (stringline >> token)
		{
			if (token == "*")
			{
				stringline >> token;
				if (token == "Mesh")
				{
					getline(infile, line);
					outFile << line << endl;
					setMeshInformation(line);
					break;
				}
				if (token == "AXIAL")
				{
					getline(infile, line);
					outFile << line << endl;
					setAxialInformation(line);
					break;
				}
				if (token == "EDGE")
				{
					stringline >> token;
					int edgeIDTemperary = stod(token);
					string lineType;
					string linePosition;
					getline(infile, lineType);
					outFile << lineType << endl;
					getline(infile, line);
					outFile << line << endl;
					getline(infile, linePosition);
					outFile << linePosition << endl;
					setEdgeInformation(lineType, linePosition, edgeIDTemperary, allEdges, nFineMesh);
					break;
				}
				if (token == "Mesh_no_of_each_coarse_mesh")
				{
					string tokenMeshId = "";
					int out0 = 1;
					while (out0 != 0)
					{
						getline(infile, line);
						outFile << line << endl;
						stringstream stringlineMeshID(line);
						while (stringlineMeshID >> tokenMeshId)
						{
							if (tokenMeshId == "*")
							{
								out0 = 0;
								stringlineMeshID >> token;
								break;
							}
							meshIDTemperary.push_back(stod(tokenMeshId));
						}
					}
				}
				if (token == "Material_name_of_each_mesh")
				{
					int ID_order = 1;
					string tokenMaterialType = "";
					int out0 = 1;
					while (out0 != 0)
					{
						getline(infile, line);
						stringstream stringlineMaterialType(line);
						while (stringlineMaterialType >> tokenMaterialType)
						{
							if (tokenMaterialType == "*")
							{
								out0 = 0;
								stringlineMaterialType >> token;
								break;
							}
							meshMaterialNameTemperary.push_back(tokenMaterialType);
							outFile << tokenMaterialType << "_" << ID_order << "  ";
							ID_order++;
						}
						if (out0 != 0)
						{
							outFile << endl;
						}
						else
						{
							outFile << line << endl;
						}
					}
				}
				if (token == "Temperature_name_of_each_mesh")
				{
					int ID_order = 1;
					string tokenTemperatureName = "";
					int out0 = 1;
					while (out0 != 0)
					{
						getline(infile, line);
						stringstream stringlineTemperatureName(line);
						while (stringlineTemperatureName >> tokenTemperatureName)
						{
							if (tokenTemperatureName == "*")
							{
								out0 = 0;
								stringlineTemperatureName >> token;
								break;
							}
							meshTemperatureNameTemperary.push_back(tokenTemperatureName);
							outFile << tokenTemperatureName << "_" << ID_order << "  ";
							ID_order++;
						}
						if (out0 != 0)
						{
							outFile << endl;
						}
						else
						{
							outFile << line << endl;
						}
					}
				}
			}
			break;
		}
	}
	outFile.close();

	for (int j = 0; j < axialNum; j++)
	{
		for (int i = 0; i < layerMeshNum; i++)
		{
			int meshIDtemp_ = meshIDTemperary[i] + j * layerMeshNum - 1;
			meshFaceTypeTemperary.push_back(meshMaterialNameTemperary[meshIDtemp_]);
			meshFaceTemperatureNameTemperary.push_back(meshTemperatureNameTemperary[meshIDtemp_]);
		}
	}

	setMeshFaceInformation(meshIDTemperary, meshFaceTypeTemperary, meshFaceTemperatureNameTemperary, allMeshFaces, allEdges);
	ThreeDemMeshOutput(fileNameTemperary, allMeshFaces, meshFaceTypeTemperary, nFineMesh);
	m_meshPointPtrVec.resize(axialNum * layerMeshNum);
	for (int j = 0; j < axialNum; j++)
	{
		for (int i = 0; i < layerMeshNum; i++)
		{
			int meshIDtemp_ = meshIDTemperary[i] + j * layerMeshNum;
			int index = i + j * layerMeshNum;

			if (kernelType == MeshKernelType::MHT_KERNEL) {
				Vector point, norm;
				for (auto& edge : allMeshFaces[i].faceEdges) {
					if (edge.edgeType == 3) {
						point = edge.arcCenter;
						norm = edge.arcAxisDir;
						break;
					}
				}
				std::ifstream ifs(fileNameTemperary[index]);
				m_meshPointPtrVec[index] = std::make_shared<MHTMocMeshPoint>(
					meshIDtemp_, ifs, allMeshFaces[i].curveInfo, point, norm,
					meshFaceTypeTemperary[index], meshFaceTemperatureNameTemperary[index]);
			}


			const char* removeFile = fileNameTemperary[index].data();
			if (remove(removeFile)) //delete file
			{
				cout << "delete file fail" << endl;
			}
		}
	}
}

void MOCMesh::setMeshInformation(string line)
{
	stringstream stringline(line);
	string token;
	stringline >> token;
	layerMeshNum = stod(token);  //mesh node number
	stringline >> token;
	EdgeNum = stod(token);  //edge number
	stringline >> token;
	coarseMeshNum = stod(token);  //coarse mesh number
}

void MOCMesh::setAxialInformation(string line)
{
	stringstream stringline(line);
	string token;
	stringline >> token;
	stringline >> token;    //delete 0
	stringline >> token;
	axialNum = stod(token);  //layer number in axial direction
	for (int i = 0; i < axialNum; i++)
	{
		std::pair<int, double> temp0;
		stringline >> token;
		temp0.first = stod(token);  // mesh number in axial direction
		stringline >> token;
		temp0.second = stod(token);  // mesh height in axial direction
		axialInformation.push_back(temp0);
	}
}


//set all edge objects
void MOCMesh::setEdgeInformation(string lineType, string linePosition, int edgeIDTemperary, std::vector<MOCEdge>& allEdges, int nFineMesh)
{
	vector<int> meshIDTemperary;
	stringstream stringline(lineType);
	string token;
	stringline >> token;
	int edgeTypeTemperary = stod(token);  //edge type
	int AditionArcEdgeID = 0;
	if (edgeTypeTemperary == 1)  //line
	{
		stringline >> token;
		meshIDTemperary.push_back(stod(token));  //node in the right
		stringline >> token;
		meshIDTemperary.push_back(stod(token));  //node in the left
		stringstream stringline0(linePosition);
		stringline0 >> token;
		double beginPoint_x = stod(token);  //coordinate of first point
		stringline0 >> token;
		double beginPoint_y = stod(token);
		stringline0 >> token;
		double endPoint_x = stod(token) + beginPoint_x;   //coordinate of second point
		stringline0 >> token;
		double endPoint_y = stod(token) + beginPoint_y;
		std::array<double, 3> beginPoint{ beginPoint_x, beginPoint_y, 0.0 };
		std::array<double, 3> endPoint{ endPoint_x, endPoint_y, 0.0 };
		MOCEdge edge0 = MOCEdge(beginPoint, endPoint, meshIDTemperary, edgeIDTemperary, edgeTypeTemperary);  //edge object
		edge0.arcCenter = Vector(0.0, 0.0, 0.0);
		edge0.arcAxisDir = Vector(0.0, 0.0, 0.0);
		allEdges.push_back(edge0);  //conserve edge object
	}

	if (edgeTypeTemperary == 2)//circle 
	{

	}
	if (edgeTypeTemperary == 3)//arc
	{
		stringline >> token;
		meshIDTemperary.push_back(stod(token));  //node in the right side
		stringline >> token;
		meshIDTemperary.push_back(stod(token));  //node in the left side

		double Center_point_X, Center_point_Y, Radius, Angle, Delta_angle;
		stringstream stringline0(linePosition);
		stringline0 >> token;
		Center_point_X = stod(token);//circle center coordinate
		stringline0 >> token;
		Center_point_Y = stod(token);

		stringline0 >> token;
		Radius = stod(token);//radius

		stringline0 >> token;
		Angle = stod(token);//angle of first point

		stringline0 >> token;
		Delta_angle = stod(token);//included angle
		if (Delta_angle < 0)
		{
			Delta_angle = Delta_angle + 360;
		}
		double fine_Delta_angle = Delta_angle / nFineMesh;
		for (int i = 1; i <= nFineMesh; i++)
		{
			double beginPoint_x = Center_point_X + Radius * cos(Angle * PI / 180.0);  //first point coordinate
			double beginPoint_y = Center_point_Y + Radius * sin(Angle * PI / 180.0);
			double endPoint_x = Center_point_X + Radius * cos((Angle + fine_Delta_angle) * PI / 180.0);  //first second coordinate
			double endPoint_y = Center_point_Y + Radius * sin((Angle + fine_Delta_angle) * PI / 180.0);
			std::array<double, 3> beginPoint{ beginPoint_x, beginPoint_y, 0.0 };
			std::array<double, 3> endPoint{ endPoint_x, endPoint_y, 0.0 };
			MOCEdge edge0 = MOCEdge(beginPoint, endPoint, meshIDTemperary, edgeIDTemperary, edgeTypeTemperary);  //edge object
			if (i > 1)
			{
				AditionArcEdgeID++;//为了保存多出来的弧边的ID
				edge0 = MOCEdge(beginPoint, endPoint, meshIDTemperary, edgeIDTemperary + EdgeNum + AditionArcEdgeID, edgeTypeTemperary);  //edge object
			}
			edge0.arcCenter = Vector(Center_point_X, Center_point_Y, 0.0);
			edge0.arcAxisDir = Vector(0.0, 0.0, 1.0);
			allEdges.push_back(edge0);  //conserve edge object
			Angle = Angle + fine_Delta_angle;
		}
	}

}
//set all surface objects
void MOCMesh::setMeshFaceInformation(vector<int> meshIDTransfer, vector<string> meshFaceTypeTransfer, vector<string> meshFaceTemperatureNameTransfer, std::vector<Surface>& allMeshFaces, std::vector<MOCEdge>& allEdges)
{
	for (int i = 0; i < layerMeshNum; i++)
	{
		Surface face0 = Surface(meshIDTransfer[i], meshIDTransfer[i], allEdges, meshFaceTypeTransfer[i], meshFaceTemperatureNameTransfer[i]);
		allMeshFaces.push_back(face0);
	}
}

void MOCMesh::ThreeDemMeshOutput(std::vector<std::string>& fileNameTransfer, std::vector<Surface>& allMeshFaces, std::vector<std::string>& meshFaceTypeTransfer, int nFineMesh)
{
	string filename = "";
	int H0 = 0, H1 = 0;
	int index0 = 0;
	for (int j = 0; j < axialNum; j++)
	{
		if (j == 0)
		{
			H0 = 0;
			H1 = axialInformation[j].second;
		}
		else
		{
			H0 = H0 + axialInformation[j - 1].second;
			H1 = H1 + axialInformation[j].second;
		}
		for (int i = 0; i < layerMeshNum; i++)
		{
			stringstream ssID;
			stringstream ssaxialID;
			string sID;
			string axialID;
			ssID << allMeshFaces[i].faceID;
			ssID >> sID;
			ssaxialID << (j + 1);
			ssaxialID >> axialID;
			filename = nFineMesh + "_poly" + meshFaceTypeTransfer[index0] + "_" + sID + "_" + axialID;
			index0++;

			//filename = filename + ".off";		
			//g_iProcessID进程ID，在输出临时文件时加到文件名里面，不然MPI多进程跑起来会出错		
			filename = filename + "_" + std::to_string(g_iProcessID) + ".off";
			
			fileNameTransfer.push_back(filename);
			ofstream outFile(filename);
			int pointNumPerMesh = allMeshFaces[i].facePointPosition.size() - 1;
			int ArcNumber = 0;
			for (int j = 0; j < allMeshFaces[i].faceEdges.size(); j++)
			{
				if (allMeshFaces[i].faceEdges[j].edgeType == 3)
				{
					ArcNumber++;
				}
			}
			outFile << "OFF" << endl;//off
			outFile << pointNumPerMesh * 2 << "\t" << pointNumPerMesh - (ArcNumber / nFineMesh) * (nFineMesh - 1) + 2 << "\t" << 0 << endl;;//point numbers,face number,0
			for (int j = 0; j < pointNumPerMesh; j++)  //coordinates of bottom points；
			{
				outFile << allMeshFaces[i].facePointPosition[j][0] << "\t" << allMeshFaces[i].facePointPosition[j][1] << "\t" << allMeshFaces[i].facePointPosition[j][2] + H0 << endl;
			}
			for (int j = 0; j < pointNumPerMesh; j++)  //coordinates of top points；
			{
				outFile << allMeshFaces[i].facePointPosition[j][0] << "\t" << allMeshFaces[i].facePointPosition[j][1] << "\t" << allMeshFaces[i].facePointPosition[j][2] + H1 << endl;
			}
			outFile << pointNumPerMesh << "\t";//number of bottom points
			for (int j = 0; j < pointNumPerMesh; j++)  //point id of bottom points;
			{
				outFile << pointNumPerMesh - 1 - allMeshFaces[i].facePointID[j] << "\t";//clockwise
			}
			allMeshFaces[i].curveInfo.push_back(0);
			allMeshFaces[i].curveFaceCenter.push_back(Vector(0.0, 0.0, 0.0));
			allMeshFaces[i].curveFaceAxisDir.push_back(Vector(0.0, 0.0, 0.0));
			outFile << endl;
			outFile << pointNumPerMesh << "\t";//number of top points
			for (int j = 0; j < pointNumPerMesh; j++)  //point id of top points;
			{
				outFile << allMeshFaces[i].facePointID[j] + pointNumPerMesh << "\t";//anticlockwise
			}
			allMeshFaces[i].curveInfo.push_back(0);
			allMeshFaces[i].curveFaceCenter.push_back(Vector(0.0, 0.0, 0.0));
			allMeshFaces[i].curveFaceAxisDir.push_back(Vector(0.0, 0.0, 0.0));
			outFile << endl;

			for (int j = 0; j < pointNumPerMesh; j++)  //side faces；
			{
				int presentEdgeTpye = allMeshFaces[i].faceEdges[j].edgeType;
				if (presentEdgeTpye == 1)//直线
				{
					outFile << 4 << "\t";//points number
					outFile << allMeshFaces[i].facePointID[j] << "\t";//anticlockwise
					outFile << allMeshFaces[i].facePointID[j + 1] << "\t";
					outFile << allMeshFaces[i].facePointID[j + 1] + pointNumPerMesh << "\t";
					outFile << allMeshFaces[i].facePointID[j] + pointNumPerMesh << "\t";
					outFile << endl;
					allMeshFaces[i].curveInfo.push_back(0);
					allMeshFaces[i].curveFaceCenter.push_back(Vector(0.0, 0.0, 0.0));
					allMeshFaces[i].curveFaceAxisDir.push_back(Vector(0.0, 0.0, 0.0));
				}
				else if (presentEdgeTpye == 2)//圆
				{

				}
				else       //圆弧
				{
					outFile << nFineMesh * 2 + 2 << "\t";//points number
					for (int k = 0; k <= nFineMesh; k++)
					{
						outFile << allMeshFaces[i].facePointID[j + k] << "\t";//anticlockwise
					}
					for (int k = nFineMesh; k >= 0; k--)
					{
						outFile << allMeshFaces[i].facePointID[j + k] + pointNumPerMesh << "\t";//anticlockwise
					}
					outFile << endl;
					j = j + nFineMesh - 1;//跳过多加进来的弧边
					allMeshFaces[i].curveInfo.push_back(1);
					allMeshFaces[i].curveFaceCenter.push_back(allMeshFaces[i].faceEdges[j].arcCenter);
					allMeshFaces[i].curveFaceAxisDir.push_back(allMeshFaces[i].faceEdges[j].arcAxisDir);
				}
			}
			outFile.close();  //close iostream
		}
	}
}

void MOCMesh::OutputStatus(std::string outputFileName) const {
	std::ofstream ofs(outputFileName);

	ofs << m_preContext.str();
	ofs << std::endl << "\t\t*   material definition" << std::endl;
	for (auto meshPointPtr : m_meshPointPtrVec) {
		ofs << std::endl;
		std::shared_ptr<MOCMeshPoint> p_mocMeshPoint = std::dynamic_pointer_cast<MOCMeshPoint>(meshPointPtr);
		std::string materialName = p_mocMeshPoint->GetMaterialName();
		const Medium* p_medium = nullptr;
		auto iter1 = m_mediumMap.find(materialName);
		if (iter1 != m_mediumMap.end())
			p_medium = &(iter1->second);

		std::smatch m;
		if (!std::regex_search(materialName, m, std::regex(R"(_\d+)"))) {
			materialName += "_" + std::to_string(p_mocMeshPoint->PointID());
			if (!p_medium) {
				auto iter2 = m_mediumMap.find(materialName);
				if (iter2 != m_mediumMap.end())
					p_medium = &(iter2->second);
			}
		}
			
		ofs << FormatStr("\t\t\t '%s' = MAT (  /", materialName.c_str());
		for (size_t i = 0; i < p_medium->eleFlagVec.size(); i++) {
			int eleFlag = p_medium->eleFlagVec[i];
			double eleDensity = p_medium->eleDensCalcFunVec[i](p_mocMeshPoint->GetValue(ValueType::DENSITY));
			ofs << eleFlag << "," << FormatStr("%.6f", eleDensity);

			if (i != p_medium->eleFlagVec.size() - 1)
				ofs << ";" << std::endl << "\t\t\t\t\t";
			else
				ofs << ")" << std::endl;
		}
	}

	ofs << std::endl << "\t\t*   temperature definition" << std::endl;
	for (auto meshPointPtr : m_meshPointPtrVec) {
		ofs << std::endl;
		std::shared_ptr<MOCMeshPoint> p_mocMeshPoint = std::dynamic_pointer_cast<MOCMeshPoint>(meshPointPtr);
		std::string tempName = p_mocMeshPoint->GetTemperatureName();
		std::smatch m;
		if (!std::regex_search(tempName, m, std::regex(R"(_\d+)")))
			tempName += "_" + std::to_string(p_mocMeshPoint->PointID());
		
		ofs << FormatStr("\t\t\t '%s' = TEMP(%.6lf)", tempName.c_str(), p_mocMeshPoint->GetValue(ValueType::TEMPERAURE)) << std::endl;
	}

	ofs << std::endl << m_sufContext.str();

	ofs.close();

	return;
}

void MOCMesh::InitMOCFromInputFile(std::string inputFileName) {
	std::ifstream ifs(inputFileName);
	if (!ifs.is_open()) 
	{
		Logger::LogError("cannot find the MOC data file:" + inputFileName);
		exit(EXIT_FAILURE);
	}
	m_mediumMap.clear();
	m_preContext.clear();
	m_sufContext.clear();

	std::unordered_map<std::string, double> materialDensityMap, temperatureMap;
	enum TOKEN {
		MATERIAL_DEFINITION,
		TEMPERATURE_DEFINITION,
		OTHER_DEFINITION
	} token = OTHER_DEFINITION;
	std::string detailDef;
	std::stringstream* p_currentContext = &m_preContext;
	while (!ifs.eof()) {
		std::string line;
		std::getline(ifs, line);
		TOKEN nextToken = token;
		if (std::regex_search(line, std::regex(R"(\*\s*material\s+definition)")))
			nextToken = MATERIAL_DEFINITION;
		if (std::regex_search(line, std::regex(R"(\*\s*temperature\s+definition)")))
			nextToken = TEMPERATURE_DEFINITION;
		if (std::regex_search(line, std::regex(R"(MODULE)")))
			nextToken = OTHER_DEFINITION;

		if (nextToken != token) {
			p_currentContext = &m_sufContext;
			if (token == MATERIAL_DEFINITION) {
				std::smatch outer_match;
				while (std::regex_search(detailDef, outer_match, std::regex(R"('\w+')"))) {
					std::string materialName = outer_match.str(0).substr(1, outer_match.str(0).length() - 2);
					detailDef = outer_match.suffix().str();

					std::regex_search(detailDef, outer_match, std::regex(R"(MAT[^\)]+\))"));
					std::string matInfo = outer_match.str(0);
					std::smatch innerMatch;
					Medium& medium = m_mediumMap[materialName];
					while (std::regex_search(matInfo, innerMatch, std::regex(R"([\d\.]+)"))) {
						int eleFlag = std::stoi(innerMatch.str(0));
						medium.eleFlagVec.push_back(eleFlag);

						matInfo = innerMatch.suffix().str();
						std::regex_search(matInfo, innerMatch, std::regex(R"([\d\.]+)"));
						double eleDens = std::stod(innerMatch.str(0));
						medium.eleDensCalcFunVec.push_back([eleDens](double) {return eleDens; });

						matInfo = innerMatch.suffix().str();
					}

					if (materialName.find("H2O") != materialName.npos) {
						double slackH2ODensity = 1e3;
						for (size_t i = 0; i < medium.eleFlagVec.size(); i++) {
							int eleFlag = medium.eleFlagVec[i];
							if (eleFlag == 1001) {
								double eleDensity = medium.eleDensCalcFunVec[i](0.0);
								slackH2ODensity = eleDensity / (NA * BARN * 2) * 1800.0;
								medium.eleDensCalcFunVec[i] = [](double density) {
									return density / 1800.0 * NA * BARN * 2;
								};
							}
							else if (eleFlag == 8016) {
								medium.eleDensCalcFunVec[i] = [](double density) {
									return density / 1800.0 * NA * BARN;
								};
							}
						}
						for (size_t i = 0; i < medium.eleFlagVec.size(); i++) {
							int eleFlag = medium.eleFlagVec[i];
							if (eleFlag == 5000) {
								double eleDensity = medium.eleDensCalcFunVec[i](0.0);
								medium.eleDensCalcFunVec[i] = [eleDensity, slackH2ODensity](double density) {
									return density / slackH2ODensity * eleDensity;
								};
							}
						}
						materialDensityMap[materialName] = slackH2ODensity;
					}
					else
						materialDensityMap[materialName] = 0.0;

					detailDef = outer_match.suffix().str();
				}
			}
			else if (token == TEMPERATURE_DEFINITION) {
				std::smatch outer_match;
				while (std::regex_search(detailDef, outer_match, std::regex(R"('\w+')"))) {
					std::string tempName = outer_match.str(0).substr(1, outer_match.str(0).length() - 2);
					detailDef = outer_match.suffix().str();

					std::regex_search(detailDef, outer_match, std::regex(R"(TEMP[^\)]+\))"));
					std::string tempInfo = outer_match.str(0);
					std::smatch innerMatch;
					std::regex_search(tempInfo, innerMatch, std::regex(R"([\d\.]+)"));
					double tempValue = std::stod(innerMatch.str(0));

					temperatureMap[tempName] = tempValue;
					detailDef = outer_match.suffix().str();
				}
			}
			
			detailDef.clear();
			token = nextToken;
		}

		if (token == MATERIAL_DEFINITION ||
			token == TEMPERATURE_DEFINITION) {
			detailDef += line;
		}
		else
			*p_currentContext << line << std::endl;
	}

	ifs.close();

	for (auto meshPointPtr : m_meshPointPtrVec) {
		std::shared_ptr<MOCMeshPoint> p_mocMeshPoint = std::dynamic_pointer_cast<MOCMeshPoint>(meshPointPtr);
		if (p_mocMeshPoint) {
			{
				std::smatch m;
				std::string materialName = p_mocMeshPoint->GetMaterialName();
				auto iter1 = materialDensityMap.find(materialName);
				if (iter1 != materialDensityMap.end())
					p_mocMeshPoint->SetValue(iter1->second, ValueType::DENSITY);

				if (std::regex_search(materialName, m, std::regex(R"(_\d+)"))) {
					std::string materialMetaName = m.prefix().str();
					auto iter2 = materialDensityMap.find(materialMetaName);
					if (iter2 != materialDensityMap.end()) {
						p_mocMeshPoint->SetValue(iter2->second, ValueType::DENSITY);
					}
				}
			}

			{
				std::smatch m;
				std::string temperatureName = p_mocMeshPoint->GetTemperatureName();
				auto iter1 = temperatureMap.find(temperatureName);
				if (iter1 != temperatureMap.end())
					p_mocMeshPoint->SetValue(iter1->second, ValueType::TEMPERAURE);

				if (std::regex_search(temperatureName, m, std::regex(R"(_\d+)"))) {
					std::string tempMetaName = m.prefix().str();
					auto iter2 = temperatureMap.find(tempMetaName);
					if (iter2 != temperatureMap.end()) {
						p_mocMeshPoint->SetValue(iter2->second, ValueType::TEMPERAURE);
					}
				}
			}
		}
	}


	return;
}

void MOCMesh::InitMOCHeatPower(std::string heatPowerFileName) {
	std::ifstream ifs(heatPowerFileName);
	if (!ifs.is_open()) {
		Logger::LogError("cannot find the moc data file:" + heatPowerFileName);
		exit(EXIT_FAILURE);
	}
	std::vector<double> powerInput;
	while (!ifs.eof()) {
		std::string line;
		std::getline(ifs, line);
		if (line == "") continue;
		powerInput.push_back(std::stod(line));
	}
	for (auto p_meshPoint : m_meshPointPtrVec) {
		p_meshPoint->SetValue(powerInput.at(p_meshPoint->PointID() - 1) * 1e6, ValueType::HEATPOWER);
	}
	ifs.close();
}

void MOCMesh::WriteTecplotFile
(
	std::string  mType,
	std::string fileName
)
{
	std::ofstream ofile(fileName);
	ofile << "TITLE =\"" << "polyhedron" << "\"" << endl;
	ofile << "VARIABLES = " << "\"x\"," << "\"y\"," << "\"z\"" << endl;
	for (int i = 0; i < this->m_meshPointPtrVec.size(); i++)
	{
		const MHTMeshPoint& mhtPolyhedron = dynamic_cast<const MHTMeshPoint&>(*m_meshPointPtrVec[i]);
		const MOCMeshPoint& mocPoint = dynamic_cast<const MOCMeshPoint&>(*m_meshPointPtrVec[i]);
		if (mType != mocPoint.GetMaterialName()) continue;
		mhtPolyhedron.WriteTecplotZones(ofile);
	}
	ofile.close();
	return;
}

Surface::Surface()
{
	facePointPosition.clear();
	facePointID.clear();
	faceID = 0;
	faceEdges.clear();
	faceType = "";
}

Surface::Surface(int faceID0, int nodeID, vector<MOCEdge> allEdgesTransfer, string meshFaceTypeTransfer, string meshTemperatureNameTransfer)
{
	int edgeNum0 = 0;
	faceID = faceID0;
	edgeNum0 = allEdgesTransfer.size();
	faceType = meshFaceTypeTransfer;
	face_temperatureName = meshTemperatureNameTransfer;
	for (int i = 0; i < edgeNum0; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			if (allEdgesTransfer[i].sideMeshID[j] == nodeID)
			{
				faceEdges.push_back(allEdgesTransfer[i]);
			}
		}
	}
	faceEdgeOrder(nodeID);
}

void Surface::faceEdgeOrder(int nodeID)  //save the edge in the anticlockwise
{
	int face_edge_num = faceEdges.size();
	MOCEdge mostLeftEdge;
	double minx = 1.0e8;
	double RoundingError = 1.0e-12;  //Rounding error of computer
	vector <MOCEdge> edgeTemperary;
	for (int i = 0; i < face_edge_num; i++)  //the most left edge and the edge type is line
	{
		for (int j = 0; j < 2; j++)
		{
			if (faceEdges[i].edgePoints[j][0] < minx && faceEdges[i].edgeType == 1)
			{
				minx = faceEdges[i].edgePoints[j][0];
				mostLeftEdge = faceEdges[i];
			}
		}
	}
	edgeTemperary.push_back(mostLeftEdge);
	vector<MOCEdge>faceEdgesTemporary;    //temporary variables for edge object
	for (int i = 0; i < faceEdges.size(); i++)
	{
		if (mostLeftEdge.edgeID == faceEdges[i].edgeID)continue;
		faceEdgesTemporary.push_back(faceEdges[i]);
	}
	std::array<double, 3> connectPoint{ 0.0, 0.0, 0.0 };    //connection point for two edge
	int pointOrderIndex = 0;
	//If the target nodes in the left side of the edge, the direction is anticlockwise
	if (mostLeftEdge.sideMeshID[1] == nodeID)
	{
		facePointID.push_back(pointOrderIndex);
		pointOrderIndex++;
		facePointPosition.push_back(mostLeftEdge.edgePoints[0]);
		facePointID.push_back(pointOrderIndex);
		pointOrderIndex++;
		facePointPosition.push_back(mostLeftEdge.edgePoints[1]);
		connectPoint = mostLeftEdge.edgePoints[1];
	}
	else
	{
		facePointID.push_back(pointOrderIndex);
		pointOrderIndex++;
		facePointPosition.push_back(mostLeftEdge.edgePoints[1]);
		facePointID.push_back(pointOrderIndex);
		pointOrderIndex++;
		facePointPosition.push_back(mostLeftEdge.edgePoints[0]);
		connectPoint = mostLeftEdge.edgePoints[0];
	}
	//Sort other edges and points anticlockwise

	MOCEdge presentEdge = mostLeftEdge;
	for (int i = 1; i < face_edge_num; i++)
	{
		for (int j = 0; j < faceEdgesTemporary.size(); j++)
		{
			if (faceEdgesTemporary[j].edgeID == presentEdge.edgeID)
			{
				continue;
			}
			if (fabs(faceEdgesTemporary[j].edgePoints[0][0] - connectPoint[0]) < RoundingError && fabs(faceEdgesTemporary[j].edgePoints[0][1] - connectPoint[1]) < RoundingError)
			{
				presentEdge = faceEdgesTemporary[j];
				connectPoint = faceEdgesTemporary[j].edgePoints[1];
				if (i == face_edge_num - 1)   //If it's the last edge, connect to the first point
				{
					facePointID.push_back(facePointID[0]);
				}
				else
				{
					facePointID.push_back(pointOrderIndex);
				}
				pointOrderIndex++;
				facePointPosition.push_back(faceEdgesTemporary[j].edgePoints[1]);
				edgeTemperary.push_back(faceEdgesTemporary[j]);
				faceEdgesTemporary.erase(faceEdgesTemporary.begin() + j);
				break;
			}
			if (fabs(faceEdgesTemporary[j].edgePoints[1][0] - connectPoint[0]) < RoundingError && fabs(faceEdgesTemporary[j].edgePoints[1][1] - connectPoint[1]) < RoundingError)
			{
				presentEdge = faceEdgesTemporary[j];
				connectPoint = faceEdgesTemporary[j].edgePoints[0];
				if (i == face_edge_num - 1)
				{
					facePointID.push_back(facePointID[0]);
				}
				else
				{
					facePointID.push_back(pointOrderIndex);
				}
				pointOrderIndex++;
				facePointPosition.push_back(faceEdgesTemporary[j].edgePoints[0]);
				edgeTemperary.push_back(faceEdgesTemporary[j]);
				faceEdgesTemporary.erase(faceEdgesTemporary.begin() + j);
				break;
			}
		}
	}
	faceEdgesTemporary.clear();

	for (int i = 0; i < face_edge_num; i++)   //������ʱ�뷽���ÿ������ı߽�����������
	{
		faceEdges[i] = edgeTemperary[i];
	}
}

MOCEdge::MOCEdge()
{
	edgePoints.clear();
	sideMeshID.clear();
	edgeID = 0;
	edgeType = 0;
}

MOCEdge::MOCEdge(std::array<double, 3> beginPoint, std::array<double, 3> endPoint, vector<int> meshIDTransfer, int edgeIDTransfer, int edgeTypeTransfer)
{
	int sideMeshIDnum = 0;
	edgeID = edgeIDTransfer;
	edgePoints.push_back(beginPoint);
	edgePoints.push_back(endPoint);
	edgeType = edgeTypeTransfer;
	sideMeshIDnum = meshIDTransfer.size();
	for (int i = 0; i < sideMeshIDnum; i++)
	{
		sideMeshID.push_back(meshIDTransfer[i]);
	}
}