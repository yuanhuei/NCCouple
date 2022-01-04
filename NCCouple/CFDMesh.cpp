#include "CFDMesh.h"
#include "Logger.h"
#include <fstream>
#include <future>
#include "./MHT_mesh/UnGridFactory.h"
#include "./MHT_mesh/RegionConnection.h"
#include "./MHT_mesh/Mesh.h"
#include "./MHT_field/Field.h"

CFDMesh::CFDMesh(Mesh* pmesh, MeshKernelType kernelType,int iMeshRegionZone)
{
	if (pmesh == nullptr)
		return;
	if (iMeshRegionZone < 0)
	{
		Logger::LogError("wrong iRegionID input");
		exit(EXIT_FAILURE);
	}
	int cellNum = pmesh->n_elemNum;
	std::string regionName = pmesh->st_meshName;
	m_meshPointPtrVec.resize(cellNum);
	Logger::LogInfo("reading CFD cells in region: " + regionName);
	Logger::LogInfo(FormatStr("CFD cell number = %d", cellNum));
	std::mutex mtx;
	int currentConstructMeshNum = 0;
	std::vector<std::future<void>> futureVec;
	for (int i = 0; i < cellNum; i++)
	{
		int verticesNum = pmesh->v_elem[i].v_nodeID.size();
		int faceNum = pmesh->v_elem[i].v_faceID.size();
		std::vector<std::string> offFileLineVec;
		//mapping nodeIDs for each polyhedron
		std::map<int, int> nodeID_id_map;
		//write node coordinates in form like 0.0296111	0.0292265	0
		for (int j = 0; j < verticesNum; j++) {
			int iNodeID	=pmesh->v_elem[i].v_nodeID[j];
			nodeID_id_map.insert(std::pair<int,int>(iNodeID, j));
			std::stringstream strCor;
			strCor << pmesh->v_node[iNodeID].x_ << " " << pmesh->v_node[iNodeID].y_ << " " << pmesh->v_node[iNodeID].z_;
			offFileLineVec.push_back(strCor.str());
		}
		//create faceID list in form like  4	3	2	1	0	
		for (int j = 0; j < faceNum; j++) 
		{
			int iFaceID = pmesh->v_elem[i].v_faceID[j];
			//whether or not the element is the owner, initializad false
			bool bFaceOwner_Element = false;
			//the polyhedron saved is the owner
			if (pmesh->v_face[iFaceID].n_owner == i) bFaceOwner_Element = true;
			int iNodeCount = pmesh->v_face[iFaceID].v_nodeID.size();
			std::string faceStr = std::to_string(iNodeCount) + " ";
			for (int k = 0; k < iNodeCount; k++)
			{
				//the polyhedron saved is the owner, the original node ID list is used
				if (bFaceOwner_Element)
				{
					int id = nodeID_id_map.at(pmesh->v_face[iFaceID].v_nodeID[k]);
					faceStr += std::to_string(id) + " ";
				}
				//the polyhedron saved is the owner, the opposite list sequence is used
				else
				{
					int id = nodeID_id_map.at(pmesh->v_face[iFaceID].v_nodeID[iNodeCount-k-1]);
					faceStr += std::to_string(id) + " ";
				}
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
			if (kernelType == MeshKernelType::MHT_KERNEL)
			{
				std::vector<int> curveInfo(faceNum, 0.0);
				Vector point, norm;
				m_meshPointPtrVec[i] = std::make_shared<MHTCFDMeshPoint>(i, ss, curveInfo, point, norm);
			}

			std::lock_guard<std::mutex> lg(mtx);
			currentConstructMeshNum++;
			if (currentConstructMeshNum % (cellNum / 10) == 0)
				Logger::LogInfo(FormatStr("%.2lf%% completed", currentConstructMeshNum * 100.0 / cellNum));
		};
		constructMeshFun();
	}
	for (size_t i = 0; i < futureVec.size(); i++)
	{
		futureVec[i].get();
	}
}

void CFDMesh::WriteTecplotFile(std::string fileName)
{
	std::ofstream ofile(fileName);
	ofile << "TITLE =\"" << "polyhedron" << "\"" << endl;
	ofile << "VARIABLES = " << "\"x\"," << "\"y\"," << "\"z\"" << endl;
	for (int i = 0; i < this->m_meshPointPtrVec.size(); i++)
	{
		const MHTMeshPoint& mhtPolyhedron = dynamic_cast<const MHTMeshPoint&>(*m_meshPointPtrVec[i]);
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
		mhtPolyhedron.WriteTecplotZones(ofile);
	}
	ofile.close();
	return;
}

void CFDMesh::SetFieldValue(std::vector<double>& v_value, ValueType vt)
{
	for (int i = 0; i < m_meshPointPtrVec.size(); i++)
	{
		v_value[i]=m_meshPointPtrVec[i]->GetValue(vt);
	}
	return;
}