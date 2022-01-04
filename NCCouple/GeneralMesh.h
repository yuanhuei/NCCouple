#ifndef MESHPOINT_HEADER
#define MESHPOINT_HEADER

//#include "CAGLWraper.h"
#include "MHT_polyhedron/PolyhedronSet.h"
#include <tuple>
#include <memory>
enum class ValueType
{
	TEMPERAURE,
	HEATPOWER,
	DENSITY
};

std::string NameOfValueType(ValueType vt);

enum class MeshKernelType
{
	CGAL_KERNEL,
	MHT_KERNEL
};

class MeshPoint
{
public:
	int PointID() const {
		return m_pointID;
	}

public:
	virtual double Volume() const = 0;
	virtual Vector Center() const = 0;
	virtual std::tuple<double, double, double> CentralCoordinate() const = 0;
	virtual double IntersectedVolume(const MeshPoint& other) const = 0;
	virtual int VerticesNum() const = 0;
	virtual std::tuple<double, double, double> VerticeCoordinate(int verticeID) const = 0;

public:
	virtual void SetValue(double value, ValueType vt) = 0;
	virtual double GetValue(ValueType vt) const = 0;

protected:
	MeshPoint();
	MeshPoint(int pointID) : m_pointID(pointID) {}
	virtual ~MeshPoint() {}

private:
	int m_pointID = 0;
};

class MHTMeshPoint : virtual public MeshPoint
{
public:
	double Volume() const override {
		return m_poly.GetVolume();
	}

	Vector Center() const override {
		return m_poly.GetCenter();
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
	void WriteToTecplotFile(std::string fileName) const {
		m_poly.WriteTecplotFile(fileName);
	}

	void WriteTecplotHeader(std::ofstream& ofile) const {
		m_poly.WriteTecplotHeader(ofile);
	}

	void WriteTecplotZones(std::ofstream& ofile) const {
		m_poly.WriteTecplotZones(ofile);
	}

	Vector GetCenter()
	{
		return m_poly.GetCenter();
	}

protected:
	MHTMeshPoint() = delete;
	MHTMeshPoint(std::istream& isf, std::vector<int>& curveInfoVec, Vector axisPoint, Vector axisNorm)
		:m_poly(isf, curveInfoVec, axisPoint, axisNorm)
	{
		m_poly.CalculateVolume();
		m_poly.Check();
	}
	virtual ~MHTMeshPoint() {}

private:
	PolyhedronSet m_poly;
};

class GeneralMesh
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
	GeneralMesh() = default;
	virtual ~GeneralMesh() {}

protected:
	std::vector<std::shared_ptr<MeshPoint>> m_meshPointPtrVec;
};

#endif

