#ifndef MOCMESH_HEADER
#define MOCMESH_HEADER

#include "GeneralMesh.h"
#include <algorithm>
#include <array>
#include <unordered_map>
#include <functional>

//#define nFineMesh 4

extern int g_iProcessID;
//enum class MaterialType {
//	H2O,
//	UO2,
//	He,
//	Zr4,
//	UNKNOWN
//};

class MOCMeshPoint : virtual public MeshPoint
{
public:
	MOCMeshPoint() = delete;
	MOCMeshPoint(std::string materialName, std::string temperatureName) :
		m_materialName(materialName), m_temperatureName(temperatureName) {}

public:
	void SetValue(double value, ValueType vt) override {
		switch (vt)
		{
		case ValueType::TEMPERAURE:
			m_temperature = value;
			break;
		case ValueType::HEATPOWER:
			m_heatPower = value;
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
			return m_heatPower;
			break;
		case ValueType::DENSITY:
			return m_density;
		default:
			break;
		}
		return 0.0;
	}

	std::string GetMaterialName() const {
		return m_materialName;
	}
	std::string GetTemperatureName() const {
		return m_temperatureName;
	}

private:
	std::string m_materialName;
	std::string m_temperatureName;
	double m_density = 0.0; //< Unit: kg/m3
	double m_temperature = 0.0;
	double m_heatPower = 0.0;
};


class MHTMocMeshPoint : public MOCMeshPoint, public MHTMeshPoint
{
public:
	MHTMocMeshPoint(
		int pointID,
		std::istream& isf,
		std::vector<int>& curveInfoVec,
		Vector axisPoint,
		Vector axisNorm,
		std::string materialName,
		std::string temperatureName) :
		MeshPoint(pointID), MOCMeshPoint(materialName, temperatureName), MHTMeshPoint(isf, curveInfoVec, axisPoint, axisNorm) {}
};

class MOCEdge
{
public:
	std::vector <std::array<double, 3>> edgePoints;
	std::vector <int> sideMeshID;
	int edgeID;
	int edgeType;
	Vector arcCenter;
	Vector arcAxisDir;
public:
	MOCEdge();
	MOCEdge(std::array<double, 3> beginPoint, std::array<double, 3> beginEnd, std::vector<int> meshIDTransfer, int edgeIDTransfer, int edgeTypeTransfer);
};

class Surface
{
public:
	std::vector <std::array<double, 3>> facePointPosition;
	std::vector <int> facePointID;
	std::vector<MOCEdge>faceEdges;
	std::vector <int> curveInfo;
	std::vector<Vector> curveFaceCenter;
	std::vector<Vector> curveFaceAxisDir;
	int faceID;
	std::string faceType;
	std::string face_temperatureName;
public:
	Surface();
	Surface(int faceID0, int nodeID, std::vector<MOCEdge> allEdgesTransfer, std::string meshFaceTypeTransfer, std::string meshFaceTemperatureNameTransfer);
	void faceEdgeOrder(int nodeID);
};

class MOCMesh : public GeneralMesh
{
public:
	struct Medium {
		std::string mediumName;
		std::vector<int> eleFlagVec;
		std::vector<std::function<double(double)>> eleDensCalcFunVec;
	};

public:
	MOCMesh() = delete;
	MOCMesh(std::string meshFileName, std::string inputFileName, MeshKernelType kernelType);
	void ThreeDemMeshOutput(std::vector<std::string>& fileNameTransfer, std::vector<Surface>& allMeshFaces, std::vector<std::string>& meshFaceTypeTransfer, int nFineMesh);   //output 3D mesh

public:
	void OutputStatus(std::string outputFileName) const override;
	void reOrganaziIndex();
	void WriteTecplotFile(std::string, std::string);
	void InitMOCHeatPower(std::string heatPowerFileName);

private:
	void setMeshInformation(std::string line); //set mesh information
	void setAxialInformation(std::string line); //set mesh information
	void setEdgeInformation(std::string lineType, std::string linePosition, int edgeIDTemperary, std::vector<MOCEdge>& allEdges, int nFineMesh);//set edge objects
	void setMeshFaceInformation(std::vector<int> meshIDTransfer, std::vector<std::string> meshFaceTypeTransfer, std::vector<std::string> meshFaceTemperatureNameTransfer, std::vector<Surface>& allMeshFaces, std::vector<MOCEdge>& allEdges);  //set surface objects
	void InitMOCFromInputFile(std::string inputFileName);

private:
	std::vector<std::vector<int>> m_coarseMeshInfo;
	int layerMeshNum;          //fine mesh number
	int EdgeNum;            //edge number
	int coarseMeshNum;     //coarse mesh number
	int axialNum;
	std::vector<std::pair<int, double>> axialInformation;

private:
	std::unordered_map<std::string, Medium> m_mediumMap;
	std::stringstream m_preContext, m_sufContext;
};
int CalMeshIndex(double x, double y);
int CalMeshIndexbyCFD(double x, double y);
#endif