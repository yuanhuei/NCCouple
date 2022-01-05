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
	//axis center of the tube
	Vector axisPoint;
	//axis normal of the tube
	Vector axisNorm;
	//the vector denoting theeta = 0
	Vector theetaStartNorm;
	//number of cells in cicular direction
	int circularCellNum;
	//number of cells in axial direction
	int axialCellNum;
	//radius list in radial direction
	std::vector<double> v_radius;
	//cell size in axial direction
	Scalar axialCellSize;
public:
	//Constructor with a given MOC mesh
	MOCIndex(MOCMesh&);
	//get radius list from MOCMesh
	void GetRadiusList();
	//resize the index and write IDs according to MOC cell coordinates
	void BuildUpIndex();
	//specify the radius list
	void SetRadial(std::vector<Scalar>& radiusList);
	//convert a Cartesian coordinate (x, y, z) to a Cylindrical coordinate (theeta, r, z)
	std::tuple<Scalar, Scalar, Scalar> GetCylindricalCoordinate(Scalar x, Scalar y, Scalar z);
	//get the i,j,k of v_MOCID from a given Cartesian coordinate
	std::tuple<int, int, int> GetIJKWithPoint(Scalar x, Scalar y, Scalar z);
	//get the ID of the MOC cell containing a given point
	int GetMOCIDWithPoint(Scalar x, Scalar y, Scalar z);
	//a member function written for test
	void Display()
	{
		std::cout << "axisPoint = " << axisPoint << std::endl;
		std::cout << "axisNorm = " << axisNorm << std::endl;
		std::cout << "theetaStartNorm = " << theetaStartNorm << std::endl;
		std::cout << "circularCellNum = " << circularCellNum << std::endl;
		std::cout << "axialCellSize = " << axialCellSize << std::endl;
		std::cout << "axialCellNum = " << axialCellNum << std::endl;
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