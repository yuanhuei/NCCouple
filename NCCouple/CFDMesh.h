#ifndef CFDMESH_HEADER
#define CFDMESH_HEADER

#include "Mesh.h"

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

class CGALCFDMeshPoint : public CFDMeshPoint, public CGALMeshPoint
{
public:
	CGALCFDMeshPoint(int pointID, std::istream& isf) :
		MeshPoint(pointID), CFDMeshPoint(), CGALMeshPoint(isf) {}
};

class LingCFDMeshPoint : public CFDMeshPoint, public LingMeshPoint
{
public:
	LingCFDMeshPoint(
		int pointID,
		std::istream& isf,
		std::vector<int>& curveInfoVec,
		Vector axisPoint,
		Vector axisNorm) :
		MeshPoint(pointID), CFDMeshPoint(), LingMeshPoint(isf, curveInfoVec, axisPoint, axisNorm) {}
};

class CFDMesh : public Mesh
{
public:
	CFDMesh() = delete;
	CFDMesh(std::string fileName, MeshKernelType kernelType);
};

#endif