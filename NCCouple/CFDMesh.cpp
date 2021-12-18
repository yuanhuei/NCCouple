#include "CFDMesh.h"
#include "Logger.h"
#include <fstream>
#include <future>
/*
CFDMesh::CFDMesh(std::string fileName, MeshKernelType kernelType) {
	std::ifstream infile(fileName);
	if (!infile.is_open())
	{
		Logger::LogError("cannot find the cfd data file:" + fileName);
		exit(EXIT_FAILURE);
	}
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
		for (int j = 0; j < faceNum; j++) {
			std::string faceStr;
			std::getline(infile, faceStr);
			offFileLineVec.push_back(faceStr);
		}

		auto constructMeshFun = [this, i, offFileLineVec, &mtx, &currentConstructMeshNum, kernelType,
			verticesNum, faceNum, cellNum]() {
			std::string polyDesc = FormatStr("%d %d 0", verticesNum, faceNum);
			std::stringstream ss;
			ss << "OFF" << std::endl;
			ss << polyDesc << std::endl;
			for (auto& lineStr : offFileLineVec)
				ss << lineStr << std::endl;

			if (kernelType == MeshKernelType::MHT_KERNEL)
			{
				std::vector<int> curveInfo(faceNum, 0.0);
				Vector point, norm;
				m_meshPointPtrVec[i] = std::make_shared<MHTCFDMeshPoint>(i, ss,
					curveInfo, point, norm);
			}
				
			std::lock_guard<std::mutex> lg(mtx);
			currentConstructMeshNum++;
			if(currentConstructMeshNum % (cellNum / 10) == 0)
				Logger::LogInfo(FormatStr("%.2lf%% completed", currentConstructMeshNum * 100.0 / cellNum));
		};
		//constructMeshFun();
		futureVec.push_back(std::async(std::launch::async | std::launch::deferred, constructMeshFun));
	}
	infile.close();

	for (size_t i = 0; i < futureVec.size(); i++)
		futureVec[i].get();
}
*/
#include "./MHT_mesh/UnGridFactory.h"
#include "./MHT_mesh/RegionConnection.h"
#include "./MHT_mesh/Mesh.h"
#include "./MHT_field/Field.h"

CFDMesh::CFDMesh(std::string fileName, MeshKernelType kernelType) {
	std::ifstream infile(fileName);
	if (!infile.is_open())
	{
		Logger::LogError("cannot find the cfd data file:" + fileName);
		exit(EXIT_FAILURE);
	}

	UnGridFactory meshFactoryCon(fileName, UnGridFactory::ugtFluent);
	FluentMeshBlock* FluentPtrCon = dynamic_cast<FluentMeshBlock*>(meshFactoryCon.GetPtr());
	RegionConnection Bridges;
	FluentPtrCon->Decompose(Bridges);
	Mesh* pmesh = &(FluentPtrCon->v_regionGrid[0]);



	int cellNum = pmesh->n_elemNum;
	std::string cellNumStr;
	std::getline(infile, cellNumStr);
	
	m_meshPointPtrVec.resize(cellNum);

	Logger::LogInfo("reading CFD cells...");
	Logger::LogInfo(FormatStr("CFD cell number = %d", cellNum));

	std::mutex mtx;
	int currentConstructMeshNum = 0;
	std::vector<std::future<void>> futureVec;
	/*按照如下格式输入到字符流，根据字符流生成meshpoint的数据结构
	    OFF
		8 6 0
		0	0	0.05
		- 1.44958e-026	0.042	0.05
		0.0278931	0.0688267	0.0499999
		0.0296465	0.0292642	0.0499997
		0	0	0
		0	0.042	0
		0.0278886	0.0688195	0
		0.0296111	0.0292265	0
		4	3	2	1	0
		4	1	2	6	5
		4	7	3	0	4
		4	6	2	3	7
		4	5	6	7	4
		4	0	1	5	4
	*/
	for (int i = 0; i < cellNum; i++)
	{
		int verticesNum = pmesh->v_elem[i].v_nodeID.size();
		int faceNum = pmesh->v_elem[i].v_faceID.size();
		std::vector<std::string> offFileLineVec;
		std::map<int, int> nodeID_id_map;//nodeID到重新排序的的映射 每个cell中的节点重新设置一个从0开始的ID
		//生成点坐标行 如：0.0296111	0.0292265	0
		for (int j = 0; j < verticesNum; j++) {
			std::string verticesCordinateStr;
			int iNodeID	=pmesh->v_elem[i].v_nodeID[j];
			nodeID_id_map.insert(std::pair<int,int>(iNodeID, j));
			verticesCordinateStr = std::to_string(pmesh->v_node[iNodeID].x_) + " " + std::to_string(pmesh->v_node[iNodeID].y_) + " " + std::to_string(pmesh->v_node[iNodeID].z_);
			offFileLineVec.push_back(verticesCordinateStr);
		}
		//生成面形成的点id行 如： 4	3	2	1	0	
		for (int j = 0; j < faceNum; j++) {
			int iFaceID = pmesh->v_elem[i].v_faceID[j];
			int iNodeCount = pmesh->v_face[iFaceID].v_nodeID.size();
			std::string faceStr = std::to_string(iNodeCount) + " ";
			for (int k = 0; k < iNodeCount; k++)
			{
				int id = nodeID_id_map.at(pmesh->v_face[iFaceID].v_nodeID[k]);
				faceStr += std::to_string(id)+" ";
			}
			offFileLineVec.push_back(faceStr);
		}


		auto constructMeshFun = [this, i, offFileLineVec, &mtx, &currentConstructMeshNum, kernelType,
			verticesNum, faceNum, cellNum]() {
			std::string polyDesc = FormatStr("%d %d 0", verticesNum, faceNum);
			std::stringstream ss;
			ss << "OFF" << std::endl;
			ss << polyDesc << std::endl;
			for (auto& lineStr : offFileLineVec)
				ss << lineStr << std::endl;
			//std::ofstream of("outcfdtemp");
			//of << ss.rdbuf();
			//of.close();
			
			if (kernelType == MeshKernelType::MHT_KERNEL)
			{
				std::vector<int> curveInfo(faceNum, 0.0);
				Vector point, norm;
				m_meshPointPtrVec[i] = std::make_shared<MHTCFDMeshPoint>(i, ss,
					curveInfo, point, norm);
			}

			std::lock_guard<std::mutex> lg(mtx);
			currentConstructMeshNum++;
			if (currentConstructMeshNum % (cellNum / 10) == 0)
				Logger::LogInfo(FormatStr("%.2lf%% completed", currentConstructMeshNum * 100.0 / cellNum));
		};
		constructMeshFun();
		//futureVec.push_back(std::async(std::launch::async | std::launch::deferred, constructMeshFun));
	}
	infile.close();

	for (size_t i = 0; i < futureVec.size(); i++)
		futureVec[i].get();
}

void CFDMesh::WriteTecplotFile(std::string fileName)
{
	std::ofstream ofile(fileName);
	ofile << "TITLE =\"" << "polyhedron" << "\"" << endl;
	ofile << "VARIABLES = " << "\"x\"," << "\"y\"," << "\"z\"" << endl;
	for (int i = 0; i < this->m_meshPointPtrVec.size(); i++)
	{
		const MHTMeshPoint& mhtPolyhedron = dynamic_cast<const MHTMeshPoint&>(*m_meshPointPtrVec[i]);
		//const CFDMeshPoint& cfdPoint = dynamic_cast<const CFDMeshPoint&>(*m_meshPointPtrVec[i]);
		//if (mType != cfdPoint.GetMaterialName()) continue;
		mhtPolyhedron.WriteTecplotZones(ofile);
	}
	ofile.close();
	return;
}

void CFDMesh::WriteTecplotFile(std::string fileName, std::vector<int>& vMeshID)
{
	std::ofstream ofile(fileName);
	ofile << "TITLE =\"" << "polyhedron" << "\"" << endl;
	ofile << "VARIABLES = " << "\"x\"," << "\"y\"," << "\"z\"" << endl;
	for (int i = 0; i < vMeshID.size(); i++)
	{
		const MHTMeshPoint& mhtPolyhedron = dynamic_cast<const MHTMeshPoint&>(*m_meshPointPtrVec[vMeshID[i]]);
		//const CFDMeshPoint& cfdPoint = dynamic_cast<const CFDMeshPoint&>(*m_meshPointPtrVec[i]);
		//if (mType != cfdPoint.GetMaterialName()) continue;
		mhtPolyhedron.WriteTecplotZones(ofile);
	}
	ofile.close();
	return;
}