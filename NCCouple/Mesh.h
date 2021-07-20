#ifndef MESHPOINT_HEADER
#define MESHPOINT_HEADER

#include "CAGLWraper.h"
#include <tuple>

enum class ValueType
{
	TEMPERAURE,
	HEATPOWER,
	DENSITY
};

class MeshPoint
{	
public:
	double Volume() const {
		return m_volume;
	}
	int PointID() const {
		return m_pointID;
	}
	std::tuple<double, double, double> CentralCoordinate() const {
		return { CGAL::to_double(m_centerPoint.x()),
			CGAL::to_double(m_centerPoint.y()),
			CGAL::to_double(m_centerPoint.z()) };
	}
	double IntersectedVolume(const MeshPoint& other) const;

public:
	virtual void SetValue(double value, ValueType vt) = 0;
	virtual double GetValue(ValueType vt) const = 0;

protected:
	MeshPoint() = delete;
	MeshPoint(int pointID, std::string polyFileName);
	MeshPoint(int pointID, std::istream& isf) {
		Init(pointID, isf);
	}
	virtual ~MeshPoint() {}

private:
	void Init(int pointID, std::istream& isf);

private:
	int m_pointID = 0;
	mutable Polyhedron m_poly;
	double m_volume = 0.0;
	Kernel::Point_3 m_centerPoint;
};

class Mesh
{
public:
	int GetMeshPointNum() const {
		return m_meshPointPtrVec.size();
	}
	MeshPoint* GetMeshPointPtr(int pointID) const {
		return m_meshPointPtrVec.at(pointID).get();
	}
	void SetValueVec(const std::vector<double>& valueVec, ValueType vt) {
		for (size_t i = 0; i < valueVec.size(); i++)
			m_meshPointPtrVec[i]->SetValue(valueVec[i], vt);
	}
	std::vector<double> GetValueVec(ValueType vt) const {
		std::vector<double> res(m_meshPointPtrVec.size());
		for (size_t i = 0; i < m_meshPointPtrVec.size(); i++)
			res[i] = m_meshPointPtrVec[i]->GetValue(vt);
		return res;
	}

protected:
	Mesh() = default;
	virtual ~Mesh() {}

protected:
	std::vector<std::shared_ptr<MeshPoint>> m_meshPointPtrVec;
};

#endif

