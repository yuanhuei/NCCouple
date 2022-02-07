#ifndef INDEX_HEADER
#define INDEX_HEADER

#include <vector>
#include <tuple>
#include <string>
#include <set>
#include "MHT_common/Vector.h"
#include "MOCIndex.h"
class MOCMesh;
struct Assembly;
struct Assembly_Type;

class CellIndex {

public:
	CellIndex(Assembly_Type& mocAssemblyType, MOCMesh* mocMesh);

private:
	std::vector<double> m_x;
	std::vector<double> m_y;

	std::vector<std::vector<int>> m_cellIndex;
	Assembly_Type* pAssemblyType = nullptr;
	MOCMesh* pMOCMesh = nullptr;
	std::vector<std::shared_ptr<MOCIndex>> v_MocIndex;

public:
	int iAssembly_type;
	void buildIndex();
	int getCellIndex(Vector vPoint);
	int getMeshID(int iCellID, Vector vPoint);
	void checkCellIndex();
};

class AssemblyIndex {

public:
	AssemblyIndex(MOCMesh& mocMesh) :pMOCMesh(&mocMesh) {};
private:

	MOCMesh* pMOCMesh = nullptr;
	std::vector<std::shared_ptr<CellIndex>> v_CellIndex;
public:
	std::vector<double> m_x;
	std::vector<double> m_y;

	std::vector<std::vector<int>> m_assemblyIndex;
	void buildIndex();
	int getAssemblyIndex(Vector vPoint);
	int getAssemblyIndex(int xIndex, int yIndex)
	{
		return m_assemblyIndex[xIndex][yIndex];
	};

	//根据坐标获取燃料组件，栅元，以及网格的id
	std::tuple<int, int, int> getIndex(Vector vPoint);
	void checkAssemblyIndex();
};
#endif