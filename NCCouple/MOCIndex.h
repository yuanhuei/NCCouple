#ifndef MOCINDEX_HEADER
#define MOCINDEX_HEADER

#include <vector>
#include <tuple>
#include <string>
#include "MHT_common/Vector.h"
//#include "MOCMesh.h"
class MOCMesh;
class MeshPoint;
struct Cell;
//An index build up for fast seaching when building up Solver
class MOCIndex
{
public:
	//pointer to the MOC mesh
	MOCMesh* pMOCMesh = nullptr;
	Cell* m_pCell = nullptr;
	//index for MOC IDs written in an array
	std::vector<std::vector<std::vector<int> > > v_MOCID;


	//axis center of the tube
	Vector axisPoint;
	//axis normal of the tube
	Vector axisNorm;
	//the vector denoting theeta = 0
	Vector theetaStartNorm;
	//number of cells in cicular direction
	int circularCellNum;
	//radius list in radial direction
	std::vector<double> v_radius;
	//number of cells in axial direction
	int axialCellNum;
	//cell size in axial direction
	Scalar axialCellSize;
	//tolerance for size
	Scalar scaleTolerance;
public:
	//Constructor with a given MOC mesh
	MOCIndex(MOCMesh&);
	MOCIndex(MOCMesh&, Cell&);
	//extract information from MOC mesh
	void Initialization();
	//estimate a tolerance
	void SetTolerance();
	//get radius list from MOCMesh
	void GetRadiusList();
	//get axial cell number and size
	void SetAxialInfo();
	//set theetastartNorm and number in circular direction
	void SetCircularInfo();
	//resize the index and write IDs according to MOC cell coordinates
	void BuildUpIndex();
	//convert a Cartesian coordinate (x, y, z) to a Cylindrical coordinate (theeta, r, z)
	std::tuple<Scalar, Scalar, Scalar> GetCylindricalCoordinate(Scalar x, Scalar y, Scalar z);
	//get the i,j,k of v_MOCID from a given Cartesian coordinate
	std::tuple<int, int, int> GetIJKWithPoint(Scalar x, Scalar y, Scalar z);
	//get the ID of the MOC cell containing a given point
	int GetMOCIDWithPoint(Scalar x, Scalar y, Scalar z);
	

	void CheckIndex();
	//a member function written for test
	void Display();
};

#endif