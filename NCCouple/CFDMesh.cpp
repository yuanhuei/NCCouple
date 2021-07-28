#include "CFDMesh.h"
#include "Logger.h"
#include <fstream>
#include <future>

CFDMesh::CFDMesh(std::string fileName) {
	std::ifstream infile(fileName);
	int cellNum = 0;
	std::string cellNumStr;
	std::getline(infile, cellNumStr);
	cellNum = std::stoi(cellNumStr);
	m_meshPointPtrVec.resize(cellNum);
	
	Logger::LogInfo("reading CFD cells...");
	Logger::LogInfo(FormatStr("CFD cell number = %d", cellNum));

	std::mutex mtx;
	int currentConstructMeshNum = 0;
	std::vector<std::future<void>> futureVec;
	for (int i = 0; i < cellNum; i++)
	{
		std::vector<std::string> offFileLineVec;
		int verticesNum = 0;
		std::string verticesNumStr;
		std::getline(infile, verticesNumStr);
		verticesNum = std::stoi(verticesNumStr);
		for (int j = 0; j < verticesNum; j++) {
			std::string verticesCordinateStr;
			std::getline(infile, verticesCordinateStr);
			offFileLineVec.push_back(verticesCordinateStr);
		}
		int faceNum = 0;
		std::string faceNumStr;
		std::getline(infile, faceNumStr);
		faceNum = std::stoi(faceNumStr);
		for (int j = 0; j < verticesNum; j++) {
			std::string faceStr;
			std::getline(infile, faceStr);
			offFileLineVec.push_back(faceStr);
		}

		auto constructMeshFun = [this, i, offFileLineVec, &mtx, &currentConstructMeshNum, 
			verticesNum, faceNum, cellNum]() {
			std::string polyDesc = FormatStr("%d %d 0", verticesNum, faceNum);
			std::stringstream ss;
			ss << "OFF" << std::endl;
			ss << polyDesc << std::endl;
			for (auto& lineStr : offFileLineVec)
				ss << lineStr << std::endl;

			m_meshPointPtrVec[i] = std::make_shared<CFDMeshPoint>(i, ss);
			std::lock_guard<std::mutex> lg(mtx);
			currentConstructMeshNum++;
			if(currentConstructMeshNum % (cellNum / 10) == 0)
				Logger::LogInfo(FormatStr("%.2lf%% completed", currentConstructMeshNum * 100.0 / cellNum));
		};
		constructMeshFun();
		//futureVec.push_back(std::async(std::launch::async, constructMeshFun));
	}
	infile.close();

	for (size_t i = 0; i < futureVec.size(); i++)
		futureVec[i].get();
}