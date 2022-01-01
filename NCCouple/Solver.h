#pragma once

#include "MOCMesh.h"
#include "CFDMesh.h"
#include "MOCIndex.h"
#include <unordered_map>
#include <memory>

class Solver {
public:
	Solver() {};
	Solver(MOCMesh& mocMesh, CFDMesh& cfdMesh);
	Solver(MOCMesh& mocMesh, CFDMesh& cfdMesh, MOCIndex& mocIndex, std::string mName);
	Solver(MOCMesh& mocMesh, CFDMesh& cfdMesh, std::string mName);
	void CheckMappingWeights();

public:
	const MOCMesh* GetMOCMeshPtr() const {
		return m_mocMeshPtr;
	}
	const CFDMesh* GetCFDMeshPtr() const {
		return m_cfdMeshPtr;
	}
	void CFDtoMOCinterception(ValueType vt) {
		Interception(m_cfdMeshPtr, m_mocMeshPtr, vt);
		return;
	}
	void MOCtoCFDinterception(ValueType vt) {
		Interception(m_mocMeshPtr, m_cfdMeshPtr, vt);
		return;
	}

private:
	void Interception(const GeneralMesh* sourceMesh, GeneralMesh* targetMesh, ValueType vt);
	void writeMapInfortoFile();
	int readMapInfor();

private:
	MOCMesh* m_mocMeshPtr = nullptr;
	CFDMesh* m_cfdMeshPtr = nullptr;
	//a solver is created for the mapping of a specific type;
	std::string materialName;
	std::vector<std::unordered_map<int, double>> m_CFD_MOC_Map;
	std::vector<std::unordered_map<int, double>> m_MOC_CFD_Map;
};