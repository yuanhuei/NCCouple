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

	std::vector<double> m_x;
	std::vector<double> m_y;
	std::vector<std::vector<int>> m_cellIndex;
	Assembly* pAssembly = nullptr;
	std::vector<MOCIndex> v_MocIndex;

	void buildIndex();
	int getCellIndex(Vector vPoint);


};

class AssemblyIndex {

public:
	AssemblyIndex(MOCMesh& mocMesh);
	~AssemblyIndex() {};

	std::vector<double> m_x;
	std::vector<double> m_y;
	std::vector<std::vector<int>> m_assemblyIndex;
	MOCMesh* pMOCMesh = nullptr;
	std::vector<CellIndex> v_CellIndex;

	void buildIndex();
	int getAssemblyIndex(Vector vPoint);
	//根据坐标获取燃料组件，栅元，以及网格的id
	std::tuple<int, int, int> getIndex(Vector vPoint);
};

