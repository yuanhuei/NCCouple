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
	std::vector<std::vector<std::vector<int> > > v_MOCID;
	Vector axisPoint;
	Vector axisNorm;
	Vector theetaStartNorm;
	int circularCellNum;
	int axialCellNum;
	std::vector<double> v_radius;
	Scalar axialCellSize;
public:
	MOCIndex(MOCMesh&);
	void BuildUpIndex();
	void SetAxial(Vector startPoint,Vector norm);
	void SetCircular(int num);
	void SetRadial(std::vector<Scalar>& radiusList);
	std::tuple<Scalar, Scalar, Scalar> GetCylindricalCoordinate(Scalar x, Scalar y, Scalar z);
	std::tuple<int, int, int> GetIJKWithPoint(Scalar x, Scalar y, Scalar z);
	int GetMOCIDWithPoint(Scalar x, Scalar y, Scalar z);
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
		for (int i = 0;i < this->v_MOCID.size();i++)
		{
			for (int j = 0;j < this->v_MOCID[i].size();j++)
			{
				for (int k = 0; k < this->v_MOCID[i][j].size();k++)
				{
					std::cout << i << "\t" << j << "\t" << k << "\t" << this->v_MOCID[i][j][k] << std::endl;
				}
			}
		}
		return;
	}
};

#endif