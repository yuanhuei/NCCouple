#ifndef MOCINDEX_HEADER
#define MOCINDEX_HEADER

#include "MOCMesh.h"
#include <vector>
#include <tuple>
#include <string>
#include "MHT_common/Vector.h"

//An index build up for fast seaching when building up Solver
class MOCIndex
{
public:
	//pointer to the MOC mesh
	MOCMesh* pMOCMesh = nullptr;
	//index for MOC IDs written in an array
	std::vector<std::vector<std::vector<int> > > v_MOCID;

	//四维数组：栅元mesh id 和xIndex,yIndex,zIndex获取网格ID
	std::vector<std::vector<std::vector<std::vector<int> > >> v_MOC_Mesh_ID;
	//四维数组：燃料堆 stackXindex,stackYindex,栅元meshxIndex,meshyIndex,获取栅元ID
	std::vector<std::vector<std::vector<std::vector<int>>>> v_Stack_Mesh_index;
	//堆的map,记录堆的索引到堆的坐标和宽度，N*N的大小
	std::map<std::pair<int, int>, std::tuple<Vector , double , int >> m_stackMap;
	//整个区域的信息
	std::tuple<Vector, double, int> m_Zone;

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
	
	//new added
	std::pair<int, int> MOCIndex::getIndex(std::tuple<Vector, double, int> tup_Zone, Vector vMeshpoint);
	std::pair<int,int> GetMOCIDWithPoint_(Scalar x, Scalar y, Scalar z);
	void BuildUpIndex_(Vector vLeftDown, int N_N, double iLength);

	void CheckIndex();
	//a member function written for test
	void Display()
	{
		std::cout << "axisPoint = " << axisPoint << std::endl;
		std::cout << "axisNorm = " << axisNorm << std::endl;
		std::cout << "theetaStartNorm = " << theetaStartNorm << std::endl;
		std::cout << "circularCellNum = " << circularCellNum << std::endl;
		std::cout << "axialCellSize = " << axialCellSize << std::endl;
		std::cout << "axialCellNum = " << axialCellNum << std::endl;
		std::cout << "scaleTolerance = " << scaleTolerance << std::endl;
		for (int i = 0;i < this->v_radius.size();i++)
		{
			std::cout << "v_radius[" << i << "] = " << v_radius[i] << std::endl;
		}
		for (int k = 0; k < axialCellNum;k++)
		{
			for (int i = 0;i < this->v_MOCID.size();i++)
			{
				for (int j = 0;j < this->v_MOCID[i].size();j++)
				{
					std::cout << i << "\t" << j << "\t" << k << "\t" << this->v_MOCID[i][j][k] << std::endl;
				}
			}
		}
		
		return;
	}
};

#endif