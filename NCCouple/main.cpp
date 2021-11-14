#include <iostream>
#include <algorithm>
#include "MOCMesh.h"
#include "MOCIndex.h"
#include "CFDMesh.h"
#include "Solver.h"
#include "Logger.h"
#include "MHT_polyhedron/PolyhedronSet.h"
#include<time.h>

void InitCFDMeshValue(const Mesh& cfdMesh) 
{
	for (int i = 0; i < cfdMesh.GetMeshPointNum(); i++) 
	{
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
	/*
	std::ifstream testOFFFile("poly.off");
	std::istream& is = testOFFFile;
	std::vector<int> v_curvedList;
	v_curvedList.push_back(0);
	v_curvedList.push_back(0);
	v_curvedList.push_back(0);
	v_curvedList.push_back(0);
	v_curvedList.push_back(1);
	v_curvedList.push_back(0);
	Vector center = Vector(0.63, 0.63, 0.0);
	Vector norm = Vector(0.0, 0.0, 1.0);
	PolyhedronSet testPolyhedronSet(is, v_curvedList, center, norm);
	testOFFFile.close();
	std::cout << "Before clipping, the volume is " << testPolyhedronSet.volume << std::endl;
	testPolyhedronSet.ClipIntoSubPolygons(10.02);
	std::cout << "After clipping, the volume is " << testPolyhedronSet.volume << std::endl;
	std::ifstream anotherTestOFFFile("anotherPoly.off");
	std::istream& anotherIss = anotherTestOFFFile;
	PolyhedronSet anotherPolyhedron(anotherIss, v_curvedList, center, norm);
	anotherTestOFFFile.close();
	Scalar volume = testPolyhedronSet.IntersectionVolumeWithPolyhedronSet(anotherPolyhedron);
	std::cout << "volume = " << volume << std::endl;
	Scalar volume1 = anotherPolyhedron.IntersectionVolumeWithPolyhedronSet(testPolyhedronSet);
	std::cout << "volume1 = " << volume1 << std::endl;

	*/
	//return 1;
	time_t start, end;
	MOCMesh mocMesh("pin_c1.apl", MeshKernelType::LING_KERNEL);
	MOCIndex mocIndex(mocMesh);
	mocIndex.axisNorm = Vector(0.0, 0.0, 1.0);
	mocIndex.axisPoint = Vector(0.63, 0.63, 0.0);
	mocIndex.theetaStartNorm = Vector(1.0, 0.0, 0.0);
	mocIndex.circularCellNum = 8;
	mocIndex.axialCellNum = 5;
	mocIndex.axialCellSize = 1.0;
	std::vector<Scalar> radiusList;
	radiusList.push_back(0.1024);
	radiusList.push_back(0.2048);
	radiusList.push_back(0.3072);
	radiusList.push_back(0.4096);
	radiusList.push_back(0.475);
	mocIndex.SetRadial(radiusList);
	mocIndex.BuildUpIndex();

	CFDMesh cfdMesh("CFDCELLSCoarse.txt", MeshKernelType::LING_KERNEL);
	
	start = time(NULL);
	//Solver solver(mocMesh, cfdMesh);
	Solver solver(mocMesh, cfdMesh, mocIndex);
	end = time(NULL);

	InitCFDMeshValue(cfdMesh);
	solver.CFDtoMOCinterception(ValueType::DENSITY);
	Logger::LogInfo(FormatStr("Time for caculatation:%d second", int(difftime(end, start))));
	mocMesh.OutputStatus("pin_c1.inp");
	ConservationValidation(cfdMesh,mocMesh, ValueType::DENSITY);
	ConservationValidation(mocMesh,cfdMesh, ValueType::DENSITY);
	return 0;
}