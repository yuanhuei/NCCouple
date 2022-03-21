#pragma once

#include "MOCMesh.h"
#include "CFDMesh.h"
#include "MOCIndex.h"
#include <unordered_map>
#include "Structure.h"
class MOCMesh;
extern int g_iMpiID;
class Solver {
public:
	Solver() {};
	Solver(MOCMesh& mocMesh, CFDMesh& cfdMesh);
	//Solver(MOCMesh& mocMesh, CFDMesh& cfdMesh, MOCIndex& mocIndex, std::string mName);
	Solver(MOCMesh& mocMesh, CFDMesh& cfdMesh, std::string mName,bool bFirstCreated=false);
	void CheckMappingWeights();
	void GetMocIndexByMapValue(std::vector< SMocIndex>& vSMocIndex);

public:
	const MOCMesh* GetMOCMeshPtr() const {
		return m_mocMeshPtr;
	}
	const CFDMesh* GetCFDMeshPtr() const {
		return m_cfdMeshPtr;
	}
	void CFDtoMOCinterception(ValueType vt) {
		Interception_fromCFDToMOC(*m_cfdMeshPtr, *m_mocMeshPtr, vt);
		return;
	}
	void MOCtoCFDinterception(ValueType vt) {
		Interception_fromMocToCFD(*m_mocMeshPtr, *m_cfdMeshPtr, vt);
		return;
	}
	double GetMocMeshMapValue(SMocIndex mocIndex)
	{
		double dValue = 0;
		for (auto& iter : m_MOC_CFD_Map[mocIndex.iAssemblyIndex][mocIndex.iCellIndex][mocIndex.iMocIndex])
		{
			dValue += iter.second;
		}
		return dValue;
	}
	double GetMocMeshMapTotalValue(SMocIndex mocIndex)
	{
		double dValue = 0;
		for (auto& iter : m_MOC_CFD_MapWithID[mocIndex.iAssemblyIndex][mocIndex.iCellIndex][mocIndex.iMocIndex])
		{
			dValue += iter.second;
		}
		return dValue;
	}
	void SetMocMapValueToStruct(std::vector<STRMocMapValue>& vMocMapValue)
	{
		for (int i = 0; i < m_MOC_CFD_Map.size(); i++)
		{
			for (int j = 0; j < m_MOC_CFD_Map[i].size(); j++)
			{
				for (int k = 0; k < m_MOC_CFD_Map[i][j].size(); k++)
				{
					double dValue = 0;
					for (auto& iter : m_MOC_CFD_Map[i][j][k])
					{
						dValue += iter.second;
					}
					if (dValue != 0)
					{
						vMocMapValue.push_back(STRMocMapValue(i, j, k, dValue));
					}
				}
			}
		}
	}
private:
	//void Interception(const GeneralMesh* sourceMesh, GeneralMesh* targetMesh, ValueType vt);
	void Interception_fromMocToCFD
	(
		MOCMesh& sourceMesh,
		CFDMesh& targetMesh,
		ValueType vt
	);
	void Interception_fromCFDToMOC
	(
		CFDMesh& sourceMesh,
		MOCMesh& targetMesh,
		ValueType vt
	);
	void WriteMapInfortoFile();
	void ReadMapInfor();
	void DisplayEmptyMap();

private:
	MOCMesh* m_mocMeshPtr = nullptr;
	CFDMesh* m_cfdMeshPtr = nullptr;
	//a solver is created for the mapping of a specific type;
	std::string materialName;
	//std::vector<std::unordered_map<int, double>> m_CFD_MOC_Map;
	//std::vector<std::unordered_map<int, double>> m_MOC_CFD_Map;
	std::vector<std::vector<std::vector<std::unordered_map<int, double>>>> m_MOC_CFD_Map;
	std::vector<std::unordered_map<SMocIndex, double>> m_CFD_MOC_Map;

	std::vector < std::vector < std::vector < std::unordered_map < std::pair<int, int>, double, pair_hash>>>> m_MOC_CFD_MapWithID;

};