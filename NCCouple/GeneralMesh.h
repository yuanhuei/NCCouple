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
	virtual vector<MHT::Polygon> GetFacesOnBoxBoundary(Vector, Vector, Scalar) const = 0;
	virtual std::vector<Scalar> GetRadiusList() const = 0;
	virtual std::pair<bool, Vector> AxisCenter() const = 0;
	virtual double IntersectedVolume(const MeshPoint& other) const = 0;
	virtual int VerticesNum() const = 0;
	virtual Vector VerticeCoordinate(int verticeID) const = 0;

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

	Vector Center() const override
	{
		return m_poly.GetCenter();
	}

	vector<MHT::Polygon> GetFacesOnBoxBoundary(Vector nodeMin, Vector nodeMax, Scalar tolerance) const override
	{
		return m_poly.GetFacesOnBoxBoundary(nodeMin, nodeMax, tolerance);
	}

	std::vector<Scalar> GetRadiusList() const override
	{
		return m_poly.GetRaduisList();
	}

	std::pair<bool, Vector> AxisCenter() const override
	{
		return m_poly.GetAxisCenter();
	}

	double IntersectedVolume(const MeshPoint& other) const override;
	double IntersectedVolume(const PolyhedronSet& other) const;

	int VerticesNum() const override {
		return m_poly.v_point.size();
	}

	Vector VerticeCoordinate(int verticeID) const override {
		Vector vertice = m_poly.v_point.at(verticeID);
		return vertice;
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
	void Move(Vector& vPoint) {
		m_poly.Move(vPoint);
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

	void GetValueVec(ValueType vt, std::vector<double>&res) const
	{
		res.resize(m_meshPointPtrVec.size());
		for (size_t i = 0; i < m_meshPointPtrVec.size(); i++)
		{
			res[i] = m_meshPointPtrVec[i]->GetValue(vt);
		}
	}

protected:
	GeneralMesh() = default;
	virtual ~GeneralMesh() {}

protected:
	std::vector<std::shared_ptr<MeshPoint>> m_meshPointPtrVec;

};

#endif

