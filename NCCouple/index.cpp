#include "index.h"
#include "Logger.h"
#include "./MHT_common/SystemControl.h"
//v_X已经排序
int getindex(double x, std::vector<double>v_X)
{
	int i = 0;
	for (; i < v_X.size() - 1; i++)
	{
		if (x >= v_X[i] && x < v_X[i + 1])
			break;
	}
	return i;
}
AssemblyIndex::AssemblyIndex(MOCMesh& mocMesh):pMOCMesh(&mocMesh)
{
}

void AssemblyIndex::buildIndex()
{
	std::vector<Assembly>& v_Assembly = pMOCMesh->m_meshAssembly;
	//初始化m_x,m_y
	for (int i == ; i < v_Assembly.size(); i++)
	{

		//遍历堆，获取左下角坐标x,y，加入一个set，最后将set中的x,y放入m_x,m_y vector中

	}
	//m_assemblyIndex初始化
	m_assemblyIndex.resize(m_x.size());
	int ySize = m_y.size();
	for (int i = 0; i < m_x.size(); i++)
	{
		m_assemblyIndex[i].resize(ySize);
		std::fill(m_assemblyIndex[i].begin(), m_assemblyIndex[i].end(), -1);
	}
	for (int i == ; i < v_Assembly.size(); i++)
	{

		//遍历堆，根据堆的左下角和右上角坐标x获取x>=m_x[i],且x<m_x[i+1],同样获取y，
		//根据获得的I,J给m_assemblyIndex赋值

	}
	//调用v_CellIndex buildindex
	v_CellIndex.resize(v_Assemblysize());
	for (int i == ; i < v_Assembly.size(); i++)
	{
		
		v_CellIndex[i] = CellIndex(v_Assembly[i]);
		v_CellIndex[i].buildIndex();
	}

}

int AssemblyIndex::getAssemblyIndex(Vector vPoint)
{
	int xIndex = getindex(vPoint.x_, m_x);
	int yIndex = getindex(vPoint.y_, m_y);
	rerurn m_assemblyIndex[xIndex][yIndex];
}
std::tuple<int, int, int> AssemblyIndex::getIndex(Vector vPoint)
{
	int iAssembly = getAssemblyIndex(vPoint);
	int iCell=v_CellIndex[iAssembly].getCellIndex(vPoint);
	int iMesh = v_CellIndex[iAssembly].v_MocIndex[iCell].GetMOCIDWithPoint(vPoint.x_, vPoint.y_, vPoint.z_);
	return std::make_tuple<int, int, int>(iAssembly, iCell, iMesh);
}


void CellIndex::buildIndex()
{
	std::vector<Cell>& v_Cell = pAssembly->v_Cell;
	//初始化m_x,m_y
	for (int i == ; i < v_Cell.size(); i++)
	{

		//遍历堆中的栅元，获取左下角，右上角坐标x,y，加入一个set，最后将set中的x,y放入m_x,m_y vector中

	}
	//m_assemblyIndex初始化
	m_cellIndex.resize(m_x.size());
	int ySize = m_y.size();
	for (int i = 0; i < m_x.size(); i++)
	{
		m_cellIndex[i].resize(ySize);
		std::fill(m_cellIndex[i].begin(), m_cellIndex[i].end(), -1);//初始值赋值为-1
	}
	for (int i == ; i < v_Cell.size(); i++)
	{

		//遍历堆中栅元，根据栅元的左下角坐标x获取索引 满足x>=m_x[i],且x<m_x[i+1],同样获取y，
		//根据获得的I,J给m_cellIndex赋值

	}
	//调用v_MocIndex buildindex
	v_MocIndex.resize(v_Cell.size());
	for (int i == ; i < v_Cell.size(); i++)
	{

		v_MocIndex[i] = MOCIndex(v_Cell[i].meshPointPtrVec);
		v_MocIndex[i].buildIndex();
	}


}
int CellIndex::getCellIndex(Vector vPoint)
{
	int xIndex = getindex(vPoint.x_, m_x);
	int yIndex = getindex(vPoint.y_, m_y);
	rerurn m_cellIndex[xIndex][yIndex];

}