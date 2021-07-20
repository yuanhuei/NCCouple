#pragma once

#include "MOCMesh.h"
#include "CFDMesh.h"
#include <unordered_map>
#include <memory>

class Solver {
public:
	Solver() = delete;
	Solver(MOCMesh& mocMesh, CFDMesh& cfdMesh);

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
	void Interception(const Mesh* sourceMesh, Mesh* targetMesh, ValueType vt);

private:
	MOCMesh* m_mocMeshPtr = nullptr;
	CFDMesh* m_cfdMeshPtr = nullptr;
	std::vector<std::unordered_map<int, double>> m_CFD_MOC_Map;
	std::vector<std::unordered_map<int, double>> m_MOC_CFD_Map;
};