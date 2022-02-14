#ifndef MOCMESH_HEADER
#define MOCMESH_HEADER

#include "GeneralMesh.h"
#include <algorithm>
#include <array>
#include <unordered_map>
#include <functional>
#include"Structure.h"

//#define nFineMesh 4
class Solver;
class AssemblyIndex;
//extern int g_iProcessID;
//enum class MaterialType {
//	H2O,
//	UO2,
//	He,
//	Zr4,
//	UNKNOWN
//};
class MocMeshField
{
public:
	MocMeshField() {};

public:
	void SetValue(double value, ValueType vt) {
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
	double GetValue(ValueType vt) const {
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
	double GetVolume() const{ return m_volume; };
	void SetVolume(double volume) { m_volume = volume; };
	std::string GetMaterialName() const{ return m_materialName; };
	void SetMaterialName(std::string strName) { m_materialName = strName; };
	std::string GetTemperatureName() const { return m_temperatureName; };
	int m_iPointID = -1;
	std::string m_temperatureName;
private:
	double m_density = 0.0; //< Unit: kg/m3
	double m_temperature = 0.0;
	double m_heatPower = 0.0;
	//add 
	double m_volume = 0;
	std::string  m_materialName;
	//std::string m_temperatureName;
	

};

struct Cell
{
	std::vector<std::shared_ptr<MeshPoint>> vMeshPointPtrVec; 
	Vector vCell_LeftDownPoint, vCell_RightUpPoint;
};
 struct Assembly_Type
{
	std::vector<Cell> v_Cell;
	int iAssemblyType;
	double xLength, yLength;
	Vector vAssemblyType_LeftDownPoint, vAssemblyType_RightUpPoint;

};
struct Assembly
{
	Assembly_Type* pAssembly_type;
	Vector vAssembly_LeftDownPoint, vAssembly_RightUpPoint;
	int iAssemblyID;
	int iAssemblyType;
    std::vector<std::vector<MocMeshField>> v_field;
};



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
	struct Medium 
	{
		std::string mediumName;
		std::vector<int> eleFlagVec;
		std::vector<std::function<double(double)>> eleDensCalcFunVec;
	};

	MOCMesh() = delete;
	MOCMesh(const std::vector<std::string>& vMaterailName);
	MOCMesh(std::string meshFileName, std::string outAplFileName,
		MeshKernelType kernelType = MeshKernelType::MHT_KERNEL);

	//void ThreeDemMeshOutput(std::vector<std::string>& fileNameTransfer, std::vector<Surface>& allMeshFaces, std::vector<std::string>& meshFaceTypeTransfer, int nFineMesh);   //output 3D mesh
	void ThreeDemMeshOutput(vector< shared_ptr<stringstream>>& vStreamTemperay, std::vector<Surface>& allMeshFaces, std::vector<std::string>& meshFaceTypeTransfer, int nFineMesh);   //output 3D mesh
	void InitMOCFromInputFile(std::string inputFileName);
	void InitMOCHeatPower(std::string heatPowerFileName);

	void OutputStatus(std::string outputFileName) const override;
	void WriteTecplotFile(std::string, std::string);
	void WriteSurfaceTecplotFile(std::string);
	void WriteHeatPowerTxtFile();
	std::pair<int, Scalar> GetAxialInformation();
	//
	SMocIndex getIndex(Vector vPoint);
	//move coppied mesh to system coordinate 
	void MoveMeshToSysCoordinate(shared_ptr<MHTMocMeshPoint>& ptrMocMesh, const SMocIndex& sMocIndex)
	{
		double x = m_vAssembly[sMocIndex.iAssemblyIndex].pAssembly_type->vAssemblyType_LeftDownPoint.x_-m_vAssembly[sMocIndex.iAssemblyIndex].vAssembly_LeftDownPoint.x_;
		double y= m_vAssembly[sMocIndex.iAssemblyIndex].pAssembly_type->vAssemblyType_LeftDownPoint.y_-m_vAssembly[sMocIndex.iAssemblyIndex].vAssembly_LeftDownPoint.y_;
		ptrMocMesh = std::make_shared<MHTMocMeshPoint>(*std::dynamic_pointer_cast<MHTMocMeshPoint> (GetMocMeshPointPtr(sMocIndex)));
		Vector vPoint(-x, -y, 0);
		ptrMocMesh->Move(vPoint);
	};
	//move point to system coordinate
	void MovePoint(Vector& vPoint, int iAssembly)
	{
		double x = m_vAssembly[iAssembly].pAssembly_type->vAssemblyType_LeftDownPoint.x_ - m_vAssembly[iAssembly].vAssembly_LeftDownPoint.x_;
		double y = m_vAssembly[iAssembly].pAssembly_type->vAssemblyType_LeftDownPoint.y_ - m_vAssembly[iAssembly].vAssembly_LeftDownPoint.y_;

		vPoint.x_ = vPoint.x_ - x;
		vPoint.y_ = vPoint.y_ - y;
	};
	void Display();

private:
	void setMeshInformation(std::string line); //set mesh information
	void setAxialInformation(std::string line); //set mesh information
	void setEdgeInformation(std::string lineType, std::string linePosition, int edgeIDTemperary, std::vector<MOCEdge>& allEdges, int nFineMesh);//set edge objects
	void setMeshFaceInformation(std::vector<int> meshIDTransfer, std::vector<std::string> meshFaceTypeTransfer, std::vector<std::string> meshFaceTemperatureNameTransfer, std::vector<Surface>& allMeshFaces, std::vector<MOCEdge>& allEdges);  //set surface objects

	std::vector<std::vector<int>> m_coarseMeshInfo;
	int layerMeshNum;          //fine mesh number
	int EdgeNum;            //edge number
	int coarseMeshNum;     //coarse mesh number
	int axialNum;
	std::vector<std::pair<int, double>> axialInformation;

	std::unordered_map<std::string, Medium> m_mediumMap;
	std::stringstream m_preContext, m_sufContext;

public://multi assembly multi cell 
	std::vector<Assembly_Type> m_vAssemblyType;
	std::vector<Assembly> m_vAssembly;
	std::shared_ptr<AssemblyIndex>  m_pAssemblyIndex;
	std::vector< SMocIndex> m_vSMocIndex;

	std::vector<std::vector<std::vector<MocMeshField>>>m_vAssemblyField;

	bool m_bSingleCell = true;
	bool m_firstCreated = true;
	MocMeshField*  GetFieldPointerAtIndex(const SMocIndex& sIndex)
	{
		return &m_vAssembly[sIndex.iAssemblyIndex].v_field[sIndex.iCellIndex][sIndex.iMocIndex];
	};
	int GetPointIDAtIndex(const SMocIndex& sIndex)const
	{
		if (m_firstCreated)
		{
			const MOCMeshPoint& mocPoint = dynamic_cast<const MOCMeshPoint&>(*GetMocMeshPointPtr(sIndex));
			return mocPoint.PointID();
		}
		else
			return m_vAssemblyField[sIndex.iAssemblyIndex][sIndex.iCellIndex][sIndex.iMocIndex].m_iPointID;

	}
	double GetVolumeAtIndex(const SMocIndex& sIndex)const
	{
		if (m_firstCreated)
		{
			const MOCMeshPoint& mocPoint = dynamic_cast<const MOCMeshPoint&>(*GetMocMeshPointPtr(sIndex));
			return mocPoint.Volume();
		}
		else
			return m_vAssemblyField[sIndex.iAssemblyIndex][sIndex.iCellIndex][sIndex.iMocIndex].GetVolume();

	}
	std::string GetMaterialNameAtIndex(const SMocIndex sIndex)const
	{
		if (m_firstCreated)
		{
			const MOCMeshPoint& mocPoint = dynamic_cast<const MOCMeshPoint&>(*GetMocMeshPointPtr(sIndex));
			return mocPoint.GetMaterialName();
		}
		else
			return m_vAssemblyField[sIndex.iAssemblyIndex][sIndex.iCellIndex][sIndex.iMocIndex].GetMaterialName();
	}
	std::string GetTemperatureNameAtIndex(const SMocIndex sIndex)const
	{
		if (m_firstCreated)
		{
			const MOCMeshPoint& mocPoint = dynamic_cast<const MOCMeshPoint&>(*GetMocMeshPointPtr(sIndex));
			return mocPoint.GetTemperatureName();
		}
		else
			return m_vAssemblyField[sIndex.iAssemblyIndex][sIndex.iCellIndex][sIndex.iMocIndex].GetTemperatureName();
	}
	double GetValueAtIndex(const SMocIndex& sIndex, ValueType vt) const
	{
		if (m_firstCreated)
			return m_vAssembly[sIndex.iAssemblyIndex].v_field[sIndex.iCellIndex][sIndex.iMocIndex].GetValue(vt);
		else
			return m_vAssemblyField[sIndex.iAssemblyIndex][sIndex.iCellIndex][sIndex.iMocIndex].GetValue(vt);

	};
	void SetValueAtIndex(const SMocIndex& sIndex, double value, ValueType vt)
	{
		if (m_firstCreated)
			m_vAssembly[sIndex.iAssemblyIndex].v_field[sIndex.iCellIndex][sIndex.iMocIndex].SetValue(value, vt);
		else
			m_vAssemblyField[sIndex.iAssemblyIndex][sIndex.iCellIndex][sIndex.iMocIndex].SetValue(value, vt);
	};
	std::shared_ptr<MeshPoint> GetMocMeshPointPtr(const SMocIndex pointID) const {
		return m_vAssembly[pointID.iAssemblyIndex].pAssembly_type->v_Cell[pointID.iCellIndex].vMeshPointPtrVec[pointID.iMocIndex];
	};
	
	void readMapFile(const std::vector<std::string>& materialList);
private:
	Assembly_Type* GetAssemblyTypePointer(int iAssemblyType);
	void InitAssembly();
	//put all index in m_vSMocIndex
	void GetAllMocIndex(std::vector< SMocIndex>& vSMocIndex);

};
#endif