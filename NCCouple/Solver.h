#pragma once

#include "MOCMesh.h"
#include "CFDMesh.h"
#include "MOCIndex.h"
#include <unordered_map>
#include "Structure.h"
class MOCMesh;
class Solver {
public:
	Solver() {};
	Solver(MOCMesh& mocMesh, CFDMesh& cfdMesh);
	//Solver(MOCMesh& mocMesh, CFDMesh& cfdMesh, MOCIndex& mocIndex, std::string mName);
	Solver(MOCMesh& mocMesh, CFDMesh& cfdMesh, std::string mName,bool bFirstCreated=false);
	void CheckMappingWeights();

public:
	const MOCMesh* GetMOCMeshPtr() const {
		return m_mocMeshPtr;
	}
	const CFDMesh* GetCFDMeshPtr() const {
		return m_cfdMeshPtr;
	}
	void CFDtoMOCinterception(ValueType vt) {
		//Interception(m_cfdMeshPtr, m_mocMeshPtr, vt);
		Interception_fromCFDToMOC(*m_cfdMeshPtr, *m_mocMeshPtr, vt);
		return;
	}
	void MOCtoCFDinterception(ValueType vt) {
		//Interception(m_mocMeshPtr, m_cfdMeshPtr, vt);
		Interception_fromMocToCFD(*m_mocMeshPtr, *m_cfdMeshPtr, vt);
		return;
	}

private:
	void Interception(const GeneralMesh* sourceMesh, GeneralMesh* targetMesh, ValueType vt);
	void Solver::Interception_fromMocToCFD
	(
		MOCMesh& sourceMesh,
		CFDMesh& targetMesh,
		ValueType vt
	);
	void Solver::Interception_fromCFDToMOC
	(
		CFDMesh& sourceMesh,
		MOCMesh& targetMesh,
		ValueType vt
	);

	void writeMapInfortoFile();
	void readMapInfor();

private:
	MOCMesh* m_mocMeshPtr = nullptr;
	CFDMesh* m_cfdMeshPtr = nullptr;
	//a solver is created for the mapping of a specific type;
	std::string materialName;
	//std::vector<std::unordered_map<int, double>> m_CFD_MOC_Map;
	//std::vector<std::unordered_map<int, double>> m_MOC_CFD_Map;
	std::vector<std::vector<std::vector<std::unordered_map<int, double>>>> m_MOC_CFD_Map;
	std::vector<std::unordered_map<SMocIndex, double, hash_name>> m_CFD_MOC_Map;

};