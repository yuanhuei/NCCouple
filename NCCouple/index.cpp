#include "index.h"
#include<iterator>
#include "Logger.h"
#include "./MHT_common/SystemControl.h"
#include "MOCMesh.h"
#include "MOCIndex.h"
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


void AssemblyIndex::buildIndex()
{
	std::vector<Assembly>& v_Assembly = pMOCMesh->m_vAssembly;
	std::vector<Assembly_Type>& v_AssemblyType = pMOCMesh->m_vAssemblyType;
	
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
	/*

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
	*/
	//call v_CellIndex buildindex
	checkAssemblyIndex();

	v_CellIndex.resize(v_AssemblyType.size());
	for (int i =0 ; i < v_AssemblyType.size(); i++)
	{
		v_CellIndex[i] = std::make_shared<CellIndex>(v_AssemblyType[i], pMOCMesh);
		v_CellIndex[i]->buildIndex();
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
	int iAssemblyIndex = getAssemblyIndex(vPoint);
	int iAsssemblyTypeIndex;
	for (int i = 0; i < v_CellIndex.size(); i++)
	{
		if (v_CellIndex[i]->iAssembly_type == pMOCMesh->m_vAssembly[iAssemblyIndex].iAssemblyType)
		{
			iAsssemblyTypeIndex = i;
			break;
		}
		else
		{
			Logger::LogError("wrong iAssemblyType ");
			exit(EXIT_FAILURE);
		}
	}
	//坐标需要经过两次转换
	Vector vPoin_on_AssemblyType = vPoint - pMOCMesh->m_vAssembly[iAssemblyIndex].vAssembly_LeftDownPoint
		+pMOCMesh->m_vAssembly[iAssemblyIndex].pAssembly_type->vAssemblyType_LeftDownPoint;
	int iCell=v_CellIndex[iAsssemblyTypeIndex]->getCellIndex(vPoin_on_AssemblyType);
	int iMesh = v_CellIndex[iAsssemblyTypeIndex]->getMeshID(iCell, vPoin_on_AssemblyType);
	return std::make_tuple(iAssemblyIndex, iCell, iMesh);
}
void AssemblyIndex::checkAssemblyIndex()
{
	std::vector<Assembly>& v_Assembly = pMOCMesh->m_vAssembly;
	std::cout << "AssemblyIndex" << std::endl;
	std::set<int> setAssembly;
	for (int i = 0; i < m_assemblyIndex.size(); i++)
	{
		for (int j = 0; j < m_assemblyIndex[i].size(); j++)
		{
			setAssembly.insert(m_assemblyIndex[i][j]);
			std::cout << "i=" << i << " j=" << j << " m_assemblyIndex=" << m_assemblyIndex[i][j] << std::endl;
		}
	}
	std::cout << "number of Assembly: " << v_Assembly.size() << std::endl;
	std::cout << "number of m_assemblyIndex:" << setAssembly.size() << std::endl;
	for (int i = 0; i < v_Assembly.size(); i++)
	{
		std::cout << "Assembly " << i << " leftdownpoint=" << v_Assembly[i].vAssembly_LeftDownPoint << std::endl;
		std::cout << "Assembly " << i << " rightuppoint=" << v_Assembly[i].vAssembly_RightUpPoint << std::endl;

	}
	if (v_Assembly.size() != setAssembly.size())
	{
		FatalError("wrong Assembly index");
	}



}


void CellIndex::buildIndex()
{
	std::vector<Cell>& v_Cell = pAssemblyType->v_Cell;
	std::vector<Vector> vLeftDownPoint, vRightUpPoint;
	vLeftDownPoint.resize(v_Cell.size());
	vRightUpPoint.resize(v_Cell.size());
	for (int i =0 ; i < v_Cell.size(); i++)
	{
		vLeftDownPoint[i] = v_Cell[i].vCell_LeftDownPoint;
		vRightUpPoint[i] = v_Cell[i].vCell_RightUpPoint;

	}
	getOrderXY(vLeftDownPoint, vRightUpPoint, m_x, m_y);

	m_cellIndex.resize(m_x.size()-1);
	int ySize = m_y.size()-1;
	for (int i = 0; i < m_x.size()-1; i++)
	{
		m_cellIndex[i].resize(ySize);
		std::fill(m_cellIndex[i].begin(), m_cellIndex[i].end(), -1);//初始值赋值为-1
	}
	for (int i =0 ; i < v_Cell.size(); i++)
	{
		//遍历堆中栅元，根据栅元中心坐标x获取索引 满足x>=m_x[i],且x<m_x[i+1],同样获取y，
		//根据获得的I,J给m_cellIndex赋值
		int xIndex = getindex((v_Cell[i].vCell_LeftDownPoint.x_+ v_Cell[i].vCell_RightUpPoint.x_)/2, m_x);
		int yIndex = getindex((v_Cell[i].vCell_LeftDownPoint.y_ + v_Cell[i].vCell_RightUpPoint.y_) / 2, m_y);
		m_cellIndex[xIndex][yIndex] = i;
	}
	//call v_MocIndex buildindex

	checkCellIndex();

	v_MocIndex.resize(v_Cell.size());
	for (int i =0 ; i < v_Cell.size(); i++)
	{
		v_MocIndex[i] = std::make_shared<MOCIndex>(*pMOCMesh,v_Cell[i]);
		v_MocIndex[i]->BuildUpIndex();
	}
}
void CellIndex::checkCellIndex()
{
	std::vector<Cell>& v_Cell = pAssemblyType->v_Cell;
	std::cout << "CellIndex" << std::endl;
	std::set<int> setCell;
	for (int i = 0; i < m_cellIndex.size(); i++)
	{
		for (int j = 0; j < m_cellIndex[i].size(); j++)
		{
			setCell.insert(m_cellIndex[i][j]);
			std::cout << "i=" << i << " j=" << j << " m_cellIndex=" << m_cellIndex[i][j] << std::endl;
		}
	}
	std::cout << "number of cell: " << v_Cell.size() << std::endl;
	std::cout << "number of m_cellIndex:" << setCell.size() << std::endl;
	if (v_Cell.size() != setCell.size())
	{
		FatalError("wrong cell index");
	}

}
CellIndex::CellIndex(Assembly_Type& mocAssemblyType, MOCMesh* mocMesh) :pAssemblyType(&mocAssemblyType), pMOCMesh(mocMesh)
{
	iAssembly_type = pAssemblyType->iAssemblyType;
};


int CellIndex::getCellIndex(Vector vPoint)
{
	int xIndex = getindex(vPoint.x_, m_x);
	int yIndex = getindex(vPoint.y_, m_y);
	return m_cellIndex[xIndex][yIndex];

}

int CellIndex::getMeshID(int iCellID, Vector vPoint)
{
	return v_MocIndex[iCellID]->GetMOCIDWithPoint(vPoint.x_, vPoint.y_, vPoint.z_);
}