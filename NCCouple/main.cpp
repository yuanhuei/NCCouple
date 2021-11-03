#include <iostream>
#include <algorithm>
#include "MOCMesh.h"
#include "MOCIndex.h"
#include "CFDMesh.h"
#include "Solver.h"
#include "Logger.h"

void InitCFDMeshValue(const Mesh& cfdMesh) {
	for (int i = 0; i < cfdMesh.GetMeshPointNum(); i++) {
		double x, y, z;
		std::tie(x, y, z) = cfdMesh.GetMeshPointPtr(i)->CentralCoordinate();
		double _x = x - 0.63;
		double _y = y - 0.63;
		double r = sqrt(_x * _x + _y * _y);
		double value = r + z + _x / r;
		cfdMesh.GetMeshPointPtr(i)->SetValue(value, ValueType::DENSITY);
	}
	
	return;
}

void ConservationValidation(const Mesh& sourceMesh, const Mesh& targetMesh, ValueType vt) {
	double sourceIntegralValue = 0.0;
	double targetIntegralValue = 0.0;

	for (int i = 0; i < sourceMesh.GetMeshPointNum(); i++) {
		double sourceValue = sourceMesh.GetMeshPointPtr(i)->GetValue(vt);
		double pointVolume = sourceMesh.GetMeshPointPtr(i)->Volume();
		sourceIntegralValue += sourceValue * pointVolume;
	}

	for (int i = 0; i < targetMesh.GetMeshPointNum(); i++) {
		double targetValue = targetMesh.GetMeshPointPtr(i)->GetValue(vt);
		double pointVolume = targetMesh.GetMeshPointPtr(i)->Volume();
		targetIntegralValue += targetValue * pointVolume;
	}

	std::string sourceMeshName, targetMeshName;
	if (dynamic_cast<const CFDMesh*>(&sourceMesh)) {
		sourceMeshName = "CFD Mesh";
		targetMeshName = "MOC Mesh";
	}
	else {
		sourceMeshName = "MOC Mesh";
		targetMeshName = "CFD Mesh";
	}
	Logger::LogInfo(FormatStr("The Integral Value of Source %s : %.6lf", sourceMeshName.c_str(), sourceIntegralValue));
	Logger::LogInfo(FormatStr("The Integral Value of Target %s : %.6lf", targetMeshName.c_str(), targetIntegralValue));

	return;
}

int main()
{
	CFDMesh cfdMesh("CFDCELLSCoarse.txt");
	MOCMesh mocMesh("pin_c1.apl");
	MOCIndex mocIndex(mocMesh);

	mocIndex.axisNorm = Vector(0.0, 0.0, 1.0);
	mocIndex.axisPoint = Vector(0.63, 0.63, 0.0);
	mocIndex.theetaStartNorm = Vector(1.0, 0.0, 0.0);
	mocIndex.circularCellNum = 8;
	mocIndex.axialCellNum = 1;
	mocIndex.axialCellSize = 0.5;
	std::vector<Scalar> radiusList;
	radiusList.push_back(0.1024);
	radiusList.push_back(0.2048);
	radiusList.push_back(0.3072);
	radiusList.push_back(0.4096);
	radiusList.push_back(0.418);
	radiusList.push_back(0.475);
	mocIndex.SetRadial(radiusList);
	mocIndex.BuildUpIndex();
	/*
	std::cout << mocIndex.GetMOCIDWithPoint(1.0, 0.5, 0.25) << std::endl;
	std::cout << mocIndex.GetMOCIDWithPoint(1.0, 0.2, 0.25) << std::endl;
	std::cout << mocIndex.GetMOCIDWithPoint(0.1, 0.2, 0.25) << std::endl;
	std::cout << mocIndex.GetMOCIDWithPoint(0.5, 0.26, 0.25) << std::endl;
	std::cout << mocIndex.GetMOCIDWithPoint(0.25, 0.75, 0.4) << std::endl;
	*/
	Solver solver(mocMesh, cfdMesh, mocIndex);

	InitCFDMeshValue(cfdMesh);
	solver.CFDtoMOCinterception(ValueType::DENSITY);

	ConservationValidation(cfdMesh, mocMesh, ValueType::DENSITY);

	mocMesh.OutputStatus("pin_c1.inp");

	return 0;

}