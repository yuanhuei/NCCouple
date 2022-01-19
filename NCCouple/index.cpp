#include "index.h"
#include<iterator>
#include "Logger.h"
#include "./MHT_common/SystemControl.h"
//get  index in x vector according to x
int getindex(double x, std::vector<double>&v_X)
{
	int i = 0;
	for (; i < v_X.size() - 1; i++)
	{
		if (x >= v_X[i] && x <= v_X[i + 1])
			break;
	}
	return i;
}
//get ordered x y vector according to leftdownpoint and rightuppoint

void getOrderXY(const std::vector<Vector>& vLeftDownPoint, const std::vector<Vector>& vRightUpPoint, 
	std::vector<double>& v_X, std::vector<double>& v_Y)
{
	std::set<double> set_X, set_Y;
	for (int i = 0; i < vLeftDownPoint.size(); i++)
	{
		set_X.insert(vLeftDownPoint[i].x_);
		set_X.insert(vRightUpPoint[i].x_);
		set_Y.insert(vLeftDownPoint[i].y_);
		set_Y.insert(vRightUpPoint[i].y_);
	}
	v_X.resize(set_X.size());
	v_Y.resize(set_Y.size());
	int i = 0;
	for (auto it = set_X.begin(); it != set_X.end(); it++)
	{
		v_X[i] = *it;
		i++;
	}
	i = 0;
	for (auto it = set_Y.begin(); it != set_Y.end(); it++)
	{
		v_Y[i] = *it;
		i++;
	}
}
AssemblyIndex::AssemblyIndex(MOCMesh& mocMesh):pMOCMesh(&mocMesh)
{
}

void AssemblyIndex::buildIndex()
{
	std::vector<Assembly>& v_Assembly = pMOCMesh->m_meshAssembly;
	//initialize m_x,m_y
	std::vector<Vector> vLeftDownPoint, vRightUpPoint;
	vLeftDownPoint.resize(v_Assembly.size());
	vRightUpPoint.resize(v_Assembly.size());
	for (int i =0; i < v_Assembly.size(); i++)
	{
		vLeftDownPoint[i] = v_Assembly[i].vAssembly_LeftDownPoint;
		vRightUpPoint[i] = v_Assembly[i].vAssembly_RightUpPoint;		
	}
	getOrderXY(vLeftDownPoint, vRightUpPoint, m_x, m_y);


	//m_assemblyIndex初始化为-1,代表无燃料组件在里面
	m_assemblyIndex.resize(m_x.size());
	int ySize = m_y.size();
	for (int i = 0; i < m_x.size(); i++)
	{
		m_assemblyIndex[i].resize(ySize);
		std::fill(m_assemblyIndex[i].begin(), m_assemblyIndex[i].end(), -1);
	}
	for (int i =0 ; i < v_Assembly.size(); i++)
	{
		int xIndex = getindex(v_Assembly[i].vAssembly_LeftDownPoint.x_, m_x);
		int yIndex= getindex(v_Assembly[i].vAssembly_LeftDownPoint.y_, m_y);
		m_assemblyIndex[xIndex][yIndex] = i;

	}
	//call v_CellIndex buildindex
	v_CellIndex.resize(v_Assembly.size());
	for (int i =0 ; i < v_Assembly.size(); i++)
	{
		v_CellIndex[i] = CellIndex(v_Assembly[i]);
		v_CellIndex[i].buildIndex();
	}
}

int AssemblyIndex::getAssemblyIndex(Vector vPoint)
{
	int xIndex = getindex(vPoint.x_, m_x);
	int yIndex = getindex(vPoint.y_, m_y);
	return m_assemblyIndex[xIndex][yIndex];
}
std::tuple<int, int, int> AssemblyIndex::getIndex(Vector vPoint)
{
	int iAssembly = getAssemblyIndex(vPoint);
	int iCell=v_CellIndex[iAssembly].getCellIndex(vPoint);
	int iMesh = v_CellIndex[iAssembly].v_MocIndex[iCell].GetMOCIDWithPoint(vPoint.x_, vPoint.y_, vPoint.z_);
	return std::make_tuple(iAssembly, iCell, iMesh);
}


void CellIndex::buildIndex()
{
	std::vector<Cell>& v_Cell = pAssembly->v_Cell;
	std::vector<Vector> vLeftDownPoint, vRightUpPoint;
	vLeftDownPoint.resize(v_Cell.size());
	vRightUpPoint.resize(v_Cell.size());
	for (int i =0 ; i < v_Cell.size(); i++)
	{
		vLeftDownPoint[i] = v_Cell[i].vCell_LeftDownPoint;
		vRightUpPoint[i] = v_Cell[i].vCell_RightUpPoint;

	}
	getOrderXY(vLeftDownPoint, vRightUpPoint, m_x, m_y);

	m_cellIndex.resize(m_x.size());
	int ySize = m_y.size();
	for (int i = 0; i < m_x.size(); i++)
	{
		m_cellIndex[i].resize(ySize);
		std::fill(m_cellIndex[i].begin(), m_cellIndex[i].end(), -1);//初始值赋值为-1
	}
	for (int i =0 ; i < v_Cell.size(); i++)
	{
		//遍历堆中栅元，根据栅元的左下角坐标x获取索引 满足x>=m_x[i],且x<m_x[i+1],同样获取y，
		//根据获得的I,J给m_cellIndex赋值
		int xIndex = getindex(v_Cell[i].vCell_LeftDownPoint.x_, m_x);
		int yIndex = getindex(v_Cell[i].vCell_RightUpPoint.y_, m_y);
		m_cellIndex[xIndex][yIndex] = i;
	}
	//call v_MocIndex buildindex
	v_MocIndex.resize(v_Cell.size());
	for (int i =0 ; i < v_Cell.size(); i++)
	{

		v_MocIndex[i] = MOCIndex(v_Cell[i].meshPointPtrVec);
		v_MocIndex[i].BuildUpIndex();
	}
}
int CellIndex::getCellIndex(Vector vPoint)
{
	int xIndex = getindex(vPoint.x_, m_x);
	int yIndex = getindex(vPoint.y_, m_y);
	return m_cellIndex[xIndex][yIndex];

}