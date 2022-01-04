#ifndef CFDMESH_HEADER
#define CFDMESH_HEADER

#include "GeneralMesh.h"
class Mesh;
template<class Type>
class Field;

class CFDMeshPoint : virtual public MeshPoint
{
public:
	CFDMeshPoint() {};

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
		case ValueType::DENSITY:
			return m_density;
		default:
			break;
		}
		return 0.0;
	}

private:
	double m_temperature = 0.0;
	double m_heatPower = 0.0;
	double m_density = 0.0;
};


class MHTCFDMeshPoint : public CFDMeshPoint, public MHTMeshPoint
{
public:
	MHTCFDMeshPoint(
		int pointID,
		std::istream& isf,
		std::vector<int>& curveInfoVec,
		Vector axisPoint,
		Vector axisNorm) :
		MeshPoint(pointID), CFDMeshPoint(), MHTMeshPoint(isf, curveInfoVec, axisPoint, axisNorm) {}
};

class CFDMesh : public GeneralMesh
{
public:
	CFDMesh() = delete;
	CFDMesh(Mesh* pmesh, MeshKernelType kernelType, int iMeshRegionZone);
	void WriteTecplotFile(std::string);
	void WriteTecplotFile(std::string,std::vector<int>&);
	void SetFieldValue(std::vector<double>& v_value, ValueType vt);
};

#endif