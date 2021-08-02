#ifndef MOCMESH_HEADER
#define MOCMESH_HEADER

#include "Mesh.h"
#include <algorithm>
#include <array>
#define nFineMesh 4   //��Ϊ��������

enum class MaterialType {
	H2O,
	UO2,
	He,
	Zr4,
	UNKNOWN
};

class MOCMeshPoint : public MeshPoint
{
public:
	MOCMeshPoint() = delete;
	MOCMeshPoint(int pointID, std::string polyFileName, MaterialType materialType, std::string temperatureName)
		:MeshPoint(pointID, polyFileName), m_materialType(materialType), m_temperatureName(temperatureName) {}

public:
	void SetValue(double value, ValueType vt) override {
		switch (vt)
		{
		case ValueType::TEMPERAURE:
			m_temperature = value;
			break;
		case ValueType::HEATPOWER:
			break;
		case ValueType::DENSITY:
			m_density = value;
			break;
		default:
			break;
		}

		return;
	}
	double GetValue(ValueType vt) const override {
		switch (vt)
		{
		case ValueType::TEMPERAURE:
			return m_temperature;
		case ValueType::HEATPOWER:
			break;
		case ValueType::DENSITY:
			return m_density;
		default:
			break;
		}
		return 0.0;
	}

	MaterialType GetMaterialType() const {
		return m_materialType;
	}
	std::string GetTemperatureName() const {
		return m_temperatureName;
	}

private:
	MaterialType m_materialType = MaterialType::H2O;
	std::string m_temperatureName;
	double m_density = 0.0;
	double m_temperature = 0.0;
};

class Edge
{
public:
	std::vector <std::array<double, 3>> edgePoints;
	std::vector <int> sideMeshID;
	int edgeID;
	int edgeType;
public:
	Edge();
	Edge(std::array<double, 3> beginPoint, std::array<double, 3> beginEnd, std::vector<int> meshIDTransfer, int edgeIDTransfer, int edgeTypeTransfer);
};

class Surface
{
public:
	std::vector <std::array<double, 3>> facePointPosition; 
	std::vector <int> facePointID;
	std::vector<Edge>faceEdges;
	int faceID;
	std::string faceType;
	std::string face_temperatureName;
public:
	Surface();
	Surface(int faceID0, int nodeID, std::vector<Edge> allEdgesTransfer, std::string meshFaceTypeTransfer, std::string meshFaceTemperatureNameTransfer);
	void faceEdgeOrder(int nodeID);
};

class MOCMesh : public Mesh
{
public:
	MOCMesh() = delete;
	MOCMesh(std::string meshFileName);
	void ThreeDemMeshOutput(std::vector<std::string>& fileNameTransfer, std::vector<Surface>& allMeshFaces);   //output 3D mesh

public:
	void OutputStatus(std::string outputFileName) const override;

private:
	void setMeshInformation(std::string line); //set mesh information
	void setEdgeInformation(std::string lineType, std::string linePosition, int edgeIDTemperary, std::vector<Edge>& allEdges);//set edge objects
	void setMeshFaceInformation(std::vector<int> meshIDTransfer, std::vector<std::string> meshFaceTypeTransfer, std::vector<std::string> meshFaceTemperatureNameTransfer, std::vector<Surface>& allMeshFaces, std::vector<Edge>& allEdges);  //set surface objects

private:
	std::vector<std::vector<int>> m_coarseMeshInfo;

	//std::vector<Surface>allMeshFaces;   //all face objects
	//std::vector<Edge>allEdges;    //all edge objects
	int MeshNum;          //fine mesh number
	int EdgeNum;            //edge number
	int coarseMeshNum;     //coarse mesh number
	double meshHighZ;   //mesh length in z direction
};
#endif
