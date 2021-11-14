#ifndef MESHPOINT_HEADER
#define MESHPOINT_HEADER

#include "CAGLWraper.h"
#include "MHT_polyhedron/PolyhedronSet.h"
#include <tuple>

enum class ValueType
{
	TEMPERAURE,
	HEATPOWER,
	DENSITY
};

enum class MeshKernelType
{
	CGAL_KERNEL,
	LING_KERNEL
};

class MeshPoint
{
public:
	int PointID() const {
		return m_pointID;
	}

public:
	virtual double Volume() const = 0;
	virtual std::tuple<double, double, double> CentralCoordinate() const = 0;
	virtual double IntersectedVolume(const MeshPoint& other) const = 0;
	virtual int VerticesNum() const = 0;
	virtual std::tuple<double, double, double> VerticeCoordinate(int verticeID) const = 0;

public:
	virtual void SetValue(double value, ValueType vt) = 0;
	virtual double GetValue(ValueType vt) const = 0;

protected:
	MeshPoint() = delete;
	MeshPoint(int pointID) : m_pointID(pointID) {}
	virtual ~MeshPoint() {}

private:
	int m_pointID = 0;
};

class CGALMeshPoint :virtual public MeshPoint
{	
public:
	double Volume() const override {
		return m_volume;
	}
	std::tuple<double, double, double> CentralCoordinate() const override {
		return std::make_tuple(CGAL::to_double(m_centerPoint.x()),
			CGAL::to_double(m_centerPoint.y()),
			CGAL::to_double(m_centerPoint.z()));
	}
	double IntersectedVolume(const MeshPoint& other) const override;
	int VerticesNum() const override {
		return m_verticesVec.size();
	}
	std::tuple<double, double, double> VerticeCoordinate(int verticeID) const override {
		return m_verticesVec.at(verticeID);
	}
	void WriteToOFF(std::ostream& out) const {
		CGAL::write_off(out, m_poly);
	}

protected:
	CGALMeshPoint() = delete;
	CGALMeshPoint(std::string polyFileName) {
		std::ifstream ifs(polyFileName);
		ifs >> m_poly;
		ifs.close();
		Init();
	}
	CGALMeshPoint(std::istream& isf) {
		isf >> m_poly;
		Init();
	}
	virtual ~CGALMeshPoint() {}

private:
	void Init();

private:
	Polyhedron m_poly;
	double m_volume = 0.0;
	Kernel::Point_3 m_centerPoint;
	std::vector<std::tuple<double, double, double>> m_verticesVec;
};

class LingMeshPoint : virtual public MeshPoint
{
public:
	double Volume() const override {
		return m_poly.GetVolume();
	}
	std::tuple<double, double, double> CentralCoordinate() const override {
		Vector center = m_poly.GetCenter();
		return std::make_tuple(center.x_, center.y_, center.z_);
	}
	double IntersectedVolume(const MeshPoint& other) const override;
	int VerticesNum() const override {
		return m_poly.v_point.size();
	}
	std::tuple<double, double, double> VerticeCoordinate(int verticeID) const override {
		Vector vertice = m_poly.v_point.at(verticeID);
		return std::make_tuple(vertice.x_, vertice.y_, vertice.z_);
	}

protected:
	LingMeshPoint() = delete;
	LingMeshPoint(std::istream& isf, std::vector<int>& curveInfoVec, Vector axisPoint, Vector axisNorm) :
		m_poly(isf, curveInfoVec, axisPoint, axisNorm) {
		m_poly.CalculateVolume();
	}
	virtual ~LingMeshPoint() {}

private:
	PolyhedronSet m_poly;
};

class Mesh
{
public:
	virtual void OutputStatus(std::string outputFileName) const {
		return;
	}

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

