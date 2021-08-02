#include "MOCMesh.h"
#include <iostream> 
#include <fstream>
using namespace std;
#define PI 3.14159265358979323846

template <class... Args>
static std::string FormatStr(const std::string& fmtStr, Args... args) {
	auto size_buf = std::snprintf(nullptr, 0, fmtStr.c_str(), args ...) + 1;
	std::unique_ptr<char[]> buf(new(std::nothrow) char[size_buf]);

	if (!buf)
		return std::string("");

	std::snprintf(buf.get(), size_buf, fmtStr.c_str(), args ...);
	return std::string(buf.get(), buf.get() + size_buf - 1);
}

double DensityConversion(std::string mediumName, std::string element, double mediumRho)
{
	typedef struct
	{
		std::string mediumName;
		std::vector<std::string> element;
		std::vector<int> elementNumber;
		double relativeMolecularMass;
	} Medium;

	std::vector<Medium> mediumInformation{
		{ "H2O",{"H","O"},{2,1},18.0 },
		{ "UO2",{"U","O"},{1,2},270.0277 },
		{ "He",{"He"},{1},4.002602 },
		{ "Zr4",{"Zr"},{4},364.896 }
	};

	double NA = 6.022e23;
	double barn = 1.0e-24;
	double rhoNuclide = 0.0;
	for (int i = 0; i < mediumInformation.size(); i++)
	{
		if (mediumName == mediumInformation[i].mediumName)
		{
			rhoNuclide = NA * mediumRho / mediumInformation[i].relativeMolecularMass * barn;
			for (int j = 0; j < mediumInformation[i].element.size(); j++)
			{
				if (mediumInformation[i].element[j] == element)
				{
					rhoNuclide = rhoNuclide * mediumInformation[i].elementNumber[j];
				}
			}
		}
	}
	return rhoNuclide;
}

MOCMesh::MOCMesh(std::string meshFileName) {
	meshHighZ = 0.5;                       //自己临时设置的网格高度，这个后面还要改
	ifstream infile(meshFileName);
	string line;
	vector<string> meshFaceTypeTemperary;
	vector<string> meshFaceTemperatureNameTemperary;
	vector<string> fileNameTemperary;
	vector<int> meshIDTemperary;

	std::vector<Surface>allMeshFaces;   //all face objects
	std::vector<Edge>allEdges;    //all edge objects

	while (getline(infile, line))  //read mesh data
	{
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
					setMeshInformation(line);
					break;
				}
				if (token == "EDGE")
				{
					stringline >> token;
					int edgeIDTemperary = stod(token);
					string lineType;
					string linePosition;
					getline(infile, lineType);
					getline(infile, line);
					getline(infile, linePosition);
					setEdgeInformation(lineType, linePosition, edgeIDTemperary, allEdges);
					break;
				}
				if (token == "Mesh_no_of_each_coarse_mesh")
				{
					string tokenMeshId = "";
					int out0 = 1;
					while (out0 != 0)
					{
						getline(infile, line);
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
							meshFaceTypeTemperary.push_back(tokenMaterialType);
						}
					}
				}
				if (token == "Temperature_name_of_each_mesh")
				{
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
							meshFaceTemperatureNameTemperary.push_back(tokenTemperatureName);
						}
					}
				}
			}
			break;
		}
	}
	setMeshFaceInformation(meshIDTemperary, meshFaceTypeTemperary, meshFaceTemperatureNameTemperary, allMeshFaces, allEdges);
	ThreeDemMeshOutput(fileNameTemperary, allMeshFaces);
	m_meshPointPtrVec.resize(MeshNum);
	static std::unordered_map<std::string, MaterialType> materialNameTypeMap{
		{ "H2O", MaterialType::H2O},
		{ "UO2", MaterialType::UO2 },
		{ "He", MaterialType::He },
		{ "Zr4", MaterialType::Zr4 }
	};
	for (int i = 0; i < MeshNum; i++)
	{
		MaterialType faceType = MaterialType::UNKNOWN;
		auto iter = materialNameTypeMap.find(meshFaceTypeTemperary[i]);
		if (iter != materialNameTypeMap.end())
			faceType = iter->second;
		m_meshPointPtrVec[i] = std::make_shared<MOCMeshPoint>(meshIDTemperary[i], fileNameTemperary[i], faceType, meshFaceTemperatureNameTemperary[i]);
		const char* removeFile = fileNameTemperary[i].data();
		if (remove(removeFile)) //删除生成的文件
		{
			cout << "删除.off文件失败" << endl;
		}
	}
}

void MOCMesh::setMeshInformation(string line)
{
	stringstream stringline(line);
	string token;
	stringline >> token;
	MeshNum = stod(token);  //mesh node number
	stringline >> token;
	EdgeNum = stod(token);  //edge number
	stringline >> token;
	coarseMeshNum = stod(token);  //coarse mesh number
}
//set all edge objects
void MOCMesh::setEdgeInformation(string lineType, string linePosition, int edgeIDTemperary, std::vector<Edge>& allEdges)
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
		Edge edge0 = Edge(beginPoint, endPoint, meshIDTemperary, edgeIDTemperary, edgeTypeTemperary);  //edge object
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
			Edge edge0 = Edge(beginPoint, endPoint, meshIDTemperary, edgeIDTemperary, edgeTypeTemperary);  //edge object
			if (i > 1)
			{
				AditionArcEdgeID++;//为了保存多出来的弧边的ID
				edge0 = Edge(beginPoint, endPoint, meshIDTemperary, edgeIDTemperary + EdgeNum + AditionArcEdgeID, edgeTypeTemperary);  //edge object
			}
			allEdges.push_back(edge0);  //conserve edge object
			Angle = Angle + fine_Delta_angle;
		}
	}

}
//set all surface objects
void MOCMesh::setMeshFaceInformation(vector<int> meshIDTransfer, vector<string> meshFaceTypeTransfer, vector<string> meshFaceTemperatureNameTransfer, std::vector<Surface>& allMeshFaces, std::vector<Edge>& allEdges)
{
	for (int i = 0; i < MeshNum; i++)
	{
		Surface face0 = Surface(meshIDTransfer[i], meshIDTransfer[i], allEdges, meshFaceTypeTransfer[i], meshFaceTemperatureNameTransfer[i]);
		allMeshFaces.push_back(face0);
	}
}
//void MOCMesh::ThreeDemMeshOutput()
//{
//	string filename = "";
//	for (int i = 0; i < MeshNum; i++)
//	{
//		stringstream ssID;
//		string sID;
//		ssID << allMeshFaces[i].faceID;
//		ssID >> sID;
//		filename = "mesh/poly" + allMeshFaces[i].faceType +"_"+ sID;
//		filename = filename + ".txt";
//		ofstream outFile(filename);
//		int pointNumPerMesh = allMeshFaces[i].facePointPosition.size() - 1;
//		outFile << pointNumPerMesh * 2 << endl;;//point numbers
//		for (int j = 0; j < pointNumPerMesh; j++)  //coordinates of bottom points；
//		{
//			outFile << allMeshFaces[i].facePointPosition[j].x_ << "\t" << allMeshFaces[i].facePointPosition[j].y_ << "\t" << allMeshFaces[i].facePointPosition[j].z_ << endl;
//		}
//		for (int j = 0; j < pointNumPerMesh; j++)  //coordinates of top points；
//		{
//			outFile << allMeshFaces[i].facePointPosition[j].x_ << "\t" << allMeshFaces[i].facePointPosition[j].y_ << "\t" << allMeshFaces[i].facePointPosition[j].z_ + meshHighZ << endl;
//		}
//		int ArcNumber = 0;
//		for (int j = 0; j < allMeshFaces[i].faceEdges.size(); j++)
//		{
//			if (allMeshFaces[i].faceEdges[j].edgeType == 3)
//			{
//				ArcNumber++;
//			}
//		}
//		outFile << pointNumPerMesh - (ArcNumber / nFineMesh) * (nFineMesh - 1) + 2 << endl;;//face number
//
//		outFile << pointNumPerMesh << "\t";//number of bottom points
//		for (int j = 0; j < pointNumPerMesh; j++)  //point id of bottom points;
//		{
//			outFile << pointNumPerMesh - 1 - allMeshFaces[i].facePointID[j] << "\t";//clockwise
//		}
//		outFile << endl;
//		outFile << pointNumPerMesh << "\t";//number of top points
//		for (int j = 0; j < pointNumPerMesh; j++)  //point id of top points;
//		{
//			outFile << allMeshFaces[i].facePointID[j] + pointNumPerMesh << "\t";//anticlockwise
//		}
//		outFile << endl;
//
//		for (int j = 0; j < pointNumPerMesh; j++)  //side faces；
//		{
//			int presentEdgeTpye = allMeshFaces[i].faceEdges[j].edgeType;
//			if (presentEdgeTpye == 1)//直线
//			{
//				outFile << 4 << "\t";//points number
//				outFile << allMeshFaces[i].facePointID[j] << "\t";//anticlockwise
//				outFile << allMeshFaces[i].facePointID[j + 1] << "\t";
//				outFile << allMeshFaces[i].facePointID[j + 1] + pointNumPerMesh << "\t";
//				outFile << allMeshFaces[i].facePointID[j] + pointNumPerMesh << "\t";
//				outFile << endl;
//			}
//			else if (presentEdgeTpye == 2)//圆
//			{
//
//			}
//			else       //圆弧
//			{
//				outFile << nFineMesh * 2 + 2 << "\t";//points number
//				for (int k = 0; k <= nFineMesh; k++)
//				{
//					outFile << allMeshFaces[i].facePointID[j + k] << "\t";//anticlockwise
//				}
//				for (int k = nFineMesh; k >= 0; k--)
//				{
//					outFile << allMeshFaces[i].facePointID[j + k] + pointNumPerMesh << "\t";//anticlockwise
//				}
//				outFile << endl;
//				j = j + nFineMesh - 1;//跳过多加进来的弧边
//			}
//		}
//		outFile.close();  //close iostream
//	}
//}

void MOCMesh::ThreeDemMeshOutput(std::vector<std::string>& fileNameTransfer, std::vector<Surface>& allMeshFaces)
{
	string filename = "";
	for (int i = 0; i < MeshNum; i++)
	{
		stringstream ssID;
		string sID;
		ssID << allMeshFaces[i].faceID;
		ssID >> sID;
		filename = "poly" + allMeshFaces[i].faceType + "_" + sID;
		filename = filename + ".off";
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
			outFile << allMeshFaces[i].facePointPosition[j][0] << "\t" << allMeshFaces[i].facePointPosition[j][1] << "\t" << allMeshFaces[i].facePointPosition[j][2] << endl;
		}
		for (int j = 0; j < pointNumPerMesh; j++)  //coordinates of top points；
		{
			outFile << allMeshFaces[i].facePointPosition[j][0] << "\t" << allMeshFaces[i].facePointPosition[j][1] << "\t" << allMeshFaces[i].facePointPosition[j][2] + meshHighZ << endl;
		}
		outFile << pointNumPerMesh << "\t";//number of bottom points
		for (int j = 0; j < pointNumPerMesh; j++)  //point id of bottom points;
		{
			outFile << pointNumPerMesh - 1 - allMeshFaces[i].facePointID[j] << "\t";//clockwise
		}
		outFile << endl;
		outFile << pointNumPerMesh << "\t";//number of top points
		for (int j = 0; j < pointNumPerMesh; j++)  //point id of top points;
		{
			outFile << allMeshFaces[i].facePointID[j] + pointNumPerMesh << "\t";//anticlockwise
		}
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
			}
		}
		outFile.close();  //close iostream
	}
}

void MOCMesh::OutputStatus(std::string outputFileName) const {
	std::ofstream ofs(outputFileName);

	std::unordered_map<MaterialType, double> materailVolumeMap;
	std::unordered_map<MaterialType, double> materailMassMap;
	std::unordered_map<std::string, double> temperatureMassMap;
	std::unordered_map<std::string, double> temperatureEnergyMap;
	for (int i = 0; i < m_meshPointPtrVec.size(); i++) {
		std::shared_ptr<MOCMeshPoint> pointPtr = dynamic_pointer_cast<MOCMeshPoint>(m_meshPointPtrVec[i]);
		MaterialType pointMaterialType = pointPtr->GetMaterialType();
		materailVolumeMap[pointMaterialType] += pointPtr->Volume();
		materailMassMap[pointMaterialType] += pointPtr->Volume() * pointPtr->GetValue(ValueType::DENSITY);

		std::string temperatureName = pointPtr->GetTemperatureName();
		temperatureMassMap[temperatureName] += pointPtr->Volume() * pointPtr->GetValue(ValueType::DENSITY);
		temperatureEnergyMap[temperatureName] += pointPtr->Volume() * pointPtr->GetValue(ValueType::DENSITY) *
			pointPtr->GetValue(ValueType::TEMPERAURE);
	}

	std::string materialDef = "\t\t*   material definition\n";
	for (auto iter : materailVolumeMap) {
		double materialAverageDensity = materailMassMap[iter.first] / (iter.second + 1E-8);
		switch (iter.first)
		{
		case MaterialType::H2O:
		{
			std::string defTemplate = R"(
			'H2O' = MAT (  /1001, %lf;
							8016, %lf;
							5000, %lf)
			)";
			double H_density = DensityConversion("H2O", "H", materialAverageDensity);
			double O_density = DensityConversion("H2O", "O", materialAverageDensity);
			std::string h2oDef = FormatStr(defTemplate, H_density, H_density, O_density);
			materialDef += h2oDef;
			break;
		}
		case MaterialType::UO2:
		{
			std::string defTemplate = R"(
			'UO2' = MAT ( /92235, %lf;
							92238, %lf;
							8001,  %lf)
			)";
			double U_density = DensityConversion("U2O", "U", materialAverageDensity);
			double O_density = DensityConversion("U2O", "O", materialAverageDensity);
			std::string uo2Def = FormatStr(defTemplate, U_density, O_density, O_density);
			materialDef += uo2Def;
			break;
		}
		case MaterialType::He:
		{
			std::string defTemplate = R"(
			'He' = MAT (/1, %lf)
			)";
			double He_density = DensityConversion("He", "He", materialAverageDensity);
			std::string heDef = FormatStr(defTemplate, He_density);
			materialDef += heDef;
			break;
		}
		case MaterialType::Zr4: {
			std::string defTemplate = R"(
			'Zr4' = MAT (/  40000, %lf;
							26000, %lf;
							8016, %lf;
							6000, %lf;
							72178, %lf;
							73181, %lf;
							41093, %lf;
							24000, %lf)
			)";
			double Zr_density = DensityConversion("Zr4", "Zr", materialAverageDensity);
			std::string zrDef = FormatStr(defTemplate, Zr_density, Zr_density, Zr_density, Zr_density,
				Zr_density, Zr_density, Zr_density, Zr_density);
			materialDef += zrDef;
			break;
		}
		default:
			break;
		}
	}

	std::string temperatureDef = "\t\t*   temperature definition\n";
	for (auto iter : temperatureMassMap) {
		double averageTemp = temperatureEnergyMap[iter.first] / (iter.second + 1E-8);
		std::string curtempDef = FormatStr("\t\t\t\t'%s' = TEMP(%lf)\n", iter.first.c_str(), averageTemp);
		temperatureDef += curtempDef;
	}

	std::string prefix = R"(
!KYLIN-2 input file

#PROGRAM  pin.hdf5

	MODULE: GENERAL
		Title  
		Date   
		Case_Name program
                         
	MODULE: ATTRIBUTION     
	)";
	std::string suffix = R"(
	MODULE: RAY_TRACE
   
		Distance_Between_Rays              0.001 
		Azimuth_Angle_Num  10   
		Polar_Angle         1  3
      
      
	MODULE: LIBRARY
		Group_Num          45
		Lib_Name       lib\hy045n18g17a.dat
        
	MODULE: RESONANCE
		Nuclide_Categorization   9
     
	MODULE: FLUX
   
		Keff_Iterative_Error 1.0E-5
		Flux_Iterative_Error 1.0E-4
		Pn_Order              0
     
	MODULE: LEAKAGE
   
	MODULE: DEPLETION 
   
		Steps_Number   1
      
		Burnup_Type   1 1  
      
		Burnup_Differential_Steps  0 25   
      
      
		Power_Density 20.19  20.19
      
		Iterary_Mode P 
   
	MODULE: EDIT_OUTPUT_FOR_CORCA3D
		Output_Name pin4_edit_for_corca3d      
      
		Output_Group 2 6.25E-7
		Output_Burnup_Point 2  0 25
      
            
		Edit_Geometry RECTANGLAR   1 1      
		Edit_Geometry_Size_X 1.26
		Edit_Geometry_Size_Y 1.26
		Type ASSEMBLY WHOLE
		Output_Cell_Layout 
			1                           
              
#PROGRAM_END
	)";

	ofs << prefix << endl;
	ofs << materialDef << endl;
	ofs << temperatureDef << endl;
	ofs << suffix;

	ofs.close();
}

Surface::Surface()
{
	facePointPosition.clear();
	facePointID.clear();
	faceID = 0;
	faceEdges.clear();
	faceType = "";
}

Surface::Surface(int faceID0, int nodeID, vector<Edge> allEdgesTransfer, string meshFaceTypeTransfer, string meshTemperatureNameTransfer)
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
	Edge mostLeftEdge;
	double minx = 1.0e8;
	double RoundingError = 1.0e-12;  //Rounding error of computer
	vector <Edge> edgeTemperary;
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
	vector<Edge>faceEdgesTemporary;    //temporary variables for edge object
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

	Edge presentEdge = mostLeftEdge;
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

Edge::Edge()
{
	edgePoints.clear();
	sideMeshID.clear();
	edgeID = 0;
	edgeType = 0;
}

Edge::Edge(std::array<double, 3> beginPoint, std::array<double, 3> endPoint, vector<int> meshIDTransfer, int edgeIDTransfer, int edgeTypeTransfer)
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
