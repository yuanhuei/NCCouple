#include <iostream>
#include <algorithm>
#include "MOCMesh.h"
#include "MOCIndex.h"
#include "CFDMesh.h"
#include "Solver.h"
#include "Logger.h"
#include "MHT_polyhedron/PolyhedronSet.h"
#include<time.h>

#include <stdio.h>
#include <string.h>
#include <mpi.h>

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

void  caculate(std::string strInputMOCfile, std::string strInputCFDfile, std::string strOutputFile)
{
	time_t start, end;
	MOCMesh mocMesh(strInputMOCfile, MeshKernelType::MHT_KERNEL);//"pin_c1.apl"
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

	CFDMesh cfdMesh(strInputCFDfile, MeshKernelType::MHT_KERNEL);//"CFDCELLS0.txt"
	
	start = time(NULL);
	//Solver solver(mocMesh, cfdMesh);
	Solver solver(mocMesh, cfdMesh, mocIndex);
	end = time(NULL);

	InitCFDMeshValue(cfdMesh);
	solver.CFDtoMOCinterception(ValueType::DENSITY);
	Logger::LogInfo(FormatStr("Time for caculatation:%d second", int(difftime(end, start))));
	mocMesh.OutputStatus(strOutputFile);// "pin_c1.inp"
	ConservationValidation(cfdMesh,mocMesh, ValueType::DENSITY);
	ConservationValidation(mocMesh,cfdMesh, ValueType::DENSITY);
}


int main(int argc, char* argv[])
{
	int numprocs, myid, source;
	MPI_Status status;
	char message[100];
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	if (myid != 0) {  //非0号进程发送消息
		std::string strInputMOCfile, strInputCFDfile, strOutputFile;
		strInputMOCfile = "pin_c1.apl";
		strInputCFDfile = "CFDCELLS0.txt";
		strOutputFile = "pin_c1.inp_"+std::to_string(myid);
		caculate(strInputMOCfile,strInputCFDfile,strOutputFile);

		std::string strMessage = "Process "+std::to_string(myid) + " caculation finished";
		strcpy(message, strMessage.data());// "caclulation finished!");
		MPI_Send(message, strlen(message) + 1, MPI_CHAR, 0, 99,
			MPI_COMM_WORLD);
	}
	else {   // myid == 0，即0号进程接收消息
		for (source = 1; source < numprocs; source++) {
			MPI_Recv(message, 100, MPI_CHAR, source, 99,
				MPI_COMM_WORLD, &status);
			//printf("接收到第%d号进程发送的消息：%s\n", source, message);
			Logger::LogInfo(FormatStr("Main process received message from No.%d process: %s\n", source, message));
		}
	}
	MPI_Finalize();
	return 0;
} /* end main */