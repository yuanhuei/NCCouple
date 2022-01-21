#pragma once
#include "MOCMesh.h"
#include <vector>
#include <tuple>
#include <string>
#include <set>
#include "MHT_common/Vector.h"
#include "MOCIndex.h"

class CellIndex {

public:
	CellIndex(Assembly& mocAssembly) :pAssembly(&mocAssembly) {};
	CellIndex() {};
	~CellIndex() {};
private:
	std::vector<double> m_x;
	std::vector<double> m_y;
	std::vector<std::vector<int>> m_cellIndex;
	Assembly* pAssembly = nullptr;
	std::vector<MOCIndex> v_MocIndex;

public:
	void buildIndex();
	int getCellIndex(Vector vPoint);
	int getMeshID(int iCellID,Vector vPoint) 
	{
		return v_MocIndex[iCellID].GetMOCIDWithPoint(vPoint.x_, vPoint.y_, vPoint.z_);
	};

};

class AssemblyIndex {

public:
	AssemblyIndex(MOCMesh& mocMesh);
	~AssemblyIndex() {};
private:
	std::vector<double> m_x;
	std::vector<double> m_y;
	std::vector<std::vector<int>> m_assemblyIndex;
	MOCMesh* pMOCMesh = nullptr;
	std::vector<CellIndex> v_CellIndex;
public:
	void buildIndex();
	int getAssemblyIndex(Vector vPoint);
	//根据坐标获取燃料组件，栅元，以及网格的id
	std::tuple<int, int, int> getIndex(Vector vPoint);
};

