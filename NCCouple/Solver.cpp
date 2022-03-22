#include "Solver.h"
#include "Logger.h"
#include <future>
#include <mutex>
#include<algorithm>
#include "index.h"
#include "MOCMesh.h"

//extern int g_iProcessID;

#define INTERSECT_JUDGE_LIMIT 1e-10

Solver::Solver(MOCMesh& mocMesh, CFDMesh& cfdMesh)
	: m_mocMeshPtr(&mocMesh), m_cfdMeshPtr(&cfdMesh)
{
	std::mutex mtx;
	m_CFD_MOC_Map.resize(cfdMesh.GetMeshPointNum());
	m_MOC_CFD_Map.resize(mocMesh.GetMeshPointNum());
	m_MOC_CFD_Map.resize(mocMesh.m_vAssembly.size());
	for (int i = 0; i < cfdMesh.GetMeshPointNum(); i++)
	{
		std::vector<std::future<void>> futureVec;
		for (int j = 0; j < mocMesh.GetMeshPointNum(); j++) {
			auto fun = [this, &mtx, i, j]() {
				const CFDMeshPoint& cfdPoint = dynamic_cast<const CFDMeshPoint&>(*m_cfdMeshPtr->GetMeshPointPtr(i));
				const MOCMeshPoint& mocPoint = dynamic_cast<const MOCMeshPoint&>(*m_mocMeshPtr->GetMeshPointPtr(j));
				double cfdPointVolume = cfdPoint.Volume();
				double mocPointVolume = mocPoint.Volume();
				double intersectedVolume = 0.0;
				if (mocPoint.GetMaterialName().find("H2O") != std::string::npos)
					intersectedVolume = cfdPoint.IntersectedVolume(mocPoint);

				if (intersectedVolume > INTERSECT_JUDGE_LIMIT) {
					std::lock_guard<std::mutex> lg(mtx);
					//m_CFD_MOC_Map[i][j] = intersectedVolume / cfdPointVolume;
					//m_MOC_CFD_Map[j][i] = intersectedVolume / mocPointVolume;
				}
			};
			//fun();
			futureVec.push_back(std::async(std::launch::async | std::launch::deferred, fun));
		}
		for (size_t j = 0; j < futureVec.size(); j++)
			futureVec[j].get();

		if (i % 100 == 0 || i == cfdMesh.GetMeshPointNum())
			Logger::LogInfo(FormatStr("Solver Initialization: %.2lf%% Completed.", i * 100.0 / cfdMesh.GetMeshPointNum()));
	}
	WriteMapInfortoFile();
}

Solver::Solver
(
	MOCMesh& mocMesh, 
	CFDMesh& cfdMesh, 
	std::string mName,
	bool bFirstCreated 
)
	:
	m_mocMeshPtr(&mocMesh),
	m_cfdMeshPtr(&cfdMesh),
	materialName(mName)
{
	if (!bFirstCreated)
	{
		ReadMapInfor();
		return;
	}
	std::mutex mtx;
	m_CFD_MOC_Map.resize(cfdMesh.GetMeshPointNum());
	m_MOC_CFD_Map.resize(mocMesh.m_vAssembly.size());
	
	for (int i = 0; i < m_MOC_CFD_Map.size(); i++)
	{
		m_MOC_CFD_Map[i].resize(mocMesh.m_vAssembly[i].pAssembly_type->v_Cell.size());
		for (int j = 0; j < m_MOC_CFD_Map[i].size(); j++)
		{
			m_MOC_CFD_Map[i][j].resize(mocMesh.m_vAssembly[i].pAssembly_type->v_Cell[j].vMeshPointPtrVec.size());
		}
	}
	
	int iNum = 0;
	int nCFDNum = cfdMesh.GetMeshPointNum();

	for (int CFDID = 0; CFDID < nCFDNum; CFDID++)
	{
		std:vector< SMocIndex> vSMoc;
		//std::vector<int>vCFD;
		//vCFD.push_back(CFDID);

		Vector P = m_cfdMeshPtr->GetMeshPointPtr(CFDID)->Center();
		SMocIndex sFirstMocIndex = mocMesh.getIndex(P);
		vSMoc.push_back(sFirstMocIndex);
		//if (CFDID == 350)
			//WriteTotecplot(*m_mocMeshPtr, *m_cfdMeshPtr, vCFD, vSMoc, "debuf_350_0.plt");
		
		if (sFirstMocIndex == SMocIndex(-1, -1, -1)) 
		{
			std::cout << "The coordinate of CFD mesh is " << P << std::endl;
			Logger::LogError("can not find  moc index corresponding to cfd Coordinate");
			continue;
		}
		std::stringstream  streamOut;

		int iAssembly = sFirstMocIndex.iAssemblyIndex, iCell = sFirstMocIndex.iCellIndex, iMoc = sFirstMocIndex.iMocIndex;

		const MHTCFDMeshPoint& cfdPoint = dynamic_cast<const MHTCFDMeshPoint&>(*m_cfdMeshPtr->GetMeshPointPtr(CFDID));

		shared_ptr<MHTMocMeshPoint> ptrMocMesh;
		m_mocMeshPtr->MoveMeshToSysCoordinate(ptrMocMesh, sFirstMocIndex);
		MHTMocMeshPoint& mocPoint =*ptrMocMesh;

		double cfdPointVolume = cfdPoint.Volume();

		double mocPointVolume = mocPoint.Volume();
		double intersectedVolume = 0.0;
		if (mocPoint.GetMaterialName() == materialName)
			intersectedVolume = cfdPoint.IntersectedVolume(mocPoint);

		if (intersectedVolume > INTERSECT_JUDGE_LIMIT) {

			m_CFD_MOC_Map[CFDID][sFirstMocIndex]=intersectedVolume / cfdPointVolume;
			m_MOC_CFD_Map[iAssembly][iCell][iMoc][CFDID] = intersectedVolume / mocPointVolume;
		}

		if ((cfdPointVolume - intersectedVolume) <= INTERSECT_JUDGE_LIMIT)
		{
			iNum++;
			continue;
		}
		//Logger::LogInfo(streamOut.str());
		std::vector<std::future<void>> futureVec;
		//coordinate on the leftdown and rightup corner of cell
		Vector leftdownPoint = mocMesh.m_vAssembly[iAssembly].pAssembly_type->v_Cell[iCell].vCell_LeftDownPoint;
		Vector rightupPoint = mocMesh.m_vAssembly[iAssembly].pAssembly_type->v_Cell[iCell].vCell_RightUpPoint;
		// move system coordinate 
		mocMesh.MovePoint(leftdownPoint, iAssembly);
		mocMesh.MovePoint(rightupPoint, iAssembly);
		Vector middlePoint = (leftdownPoint + rightupPoint) / 2;
		//get the corner point which is nearest to CFD center
		Vector vNearestPoint;
		if (P.x_ >= middlePoint.x_ && P.y_ >= middlePoint.y_)
			vNearestPoint = rightupPoint;
		else if (P.x_ < middlePoint.x_ && P.y_ < middlePoint.y_)
			vNearestPoint = leftdownPoint;
		else if (P.x_ >= middlePoint.x_ && P.y_ < middlePoint.y_)
		{
			vNearestPoint.x_ = rightupPoint.x_;
			vNearestPoint.y_ = leftdownPoint.y_;
		}
		else
		{
			vNearestPoint.x_ = leftdownPoint.x_;
			vNearestPoint.y_ = rightupPoint.y_;
		}
		//get the nearby four cell
		std::vector<SMocIndex> vMocIndex(4);
		std::vector<Vector> vPoint{ Vector(vNearestPoint.x_ + 0.5, vNearestPoint.y_ + 0.5,0),
			Vector(vNearestPoint.x_ + 0.5, vNearestPoint.y_ - 0.5,0),
			Vector(vNearestPoint.x_ - 0.5, vNearestPoint.y_ + 0.5, 0),
			Vector(vNearestPoint.x_ - 0.5, vNearestPoint.y_ - 0.5, 0) };
		//search this four cell
		
		for(int i=0;i<vPoint.size();i++)
		{
			//get cell Index
			SMocIndex sMocIndex = mocMesh.getIndex(vPoint[i]);
			if (sMocIndex == SMocIndex(-1, -1, -1))
				continue;
			vMocIndex[i] = sMocIndex;

			Vector leftdownPoint = mocMesh.m_vAssembly[sMocIndex.iAssemblyIndex].pAssembly_type->v_Cell[sMocIndex.iCellIndex].vCell_LeftDownPoint;
			Vector rightupPoint = mocMesh.m_vAssembly[sMocIndex.iAssemblyIndex].pAssembly_type->v_Cell[sMocIndex.iCellIndex].vCell_RightUpPoint;
			mocMesh.MovePoint(leftdownPoint, sMocIndex.iAssemblyIndex);
			mocMesh.MovePoint(rightupPoint, sMocIndex.iAssemblyIndex);
			PolyhedronSet box(leftdownPoint, rightupPoint);
			if (cfdPoint.IntersectedVolume(box) <= INTERSECT_JUDGE_LIMIT)
				continue;

			std::vector<int> vMocID;
			mocMesh.m_pAssemblyIndex->getNearLayerMocID(vMocID, cfdPoint, sMocIndex.iAssemblyIndex,
				sMocIndex.iCellIndex);

			//search the nearby mesh in the cell
			for (int j = 0; j < vMocID.size(); j++)
			{
				sMocIndex.iMocIndex = vMocID[j];
				vSMoc.push_back(sMocIndex);
				if (sFirstMocIndex == sMocIndex)continue;
				const MOCMeshPoint& localMOCPoint = dynamic_cast<const MOCMeshPoint&>(*m_mocMeshPtr->GetMocMeshPointPtr(sMocIndex));
				if (materialName != localMOCPoint.GetMaterialName()) continue;
				auto fun = [this, &mtx, CFDID, sMocIndex]() 
				{
					const CFDMeshPoint& cfdPoint = dynamic_cast<const CFDMeshPoint&>(*m_cfdMeshPtr->GetMeshPointPtr(CFDID));
					shared_ptr<MHTMocMeshPoint> ptrMocMesh;
					m_mocMeshPtr->MoveMeshToSysCoordinate(ptrMocMesh, sMocIndex);
					MHTMocMeshPoint& mocPoint = *ptrMocMesh;

					double cfdPointVolume = cfdPoint.Volume();
					double mocPointVolume = mocPoint.Volume();
					double intersectedVolume = 0.0;
					intersectedVolume = cfdPoint.IntersectedVolume(mocPoint);
					if (intersectedVolume > INTERSECT_JUDGE_LIMIT)
					{
						std::lock_guard<std::mutex> lg(mtx);

						m_CFD_MOC_Map[CFDID][sMocIndex] = intersectedVolume / cfdPointVolume;
						m_MOC_CFD_Map[sMocIndex.iAssemblyIndex][sMocIndex.iCellIndex][sMocIndex.iMocIndex][CFDID] = intersectedVolume / mocPointVolume;
					}

				};
				//fun();
				futureVec.push_back(std::async(std::launch::async | std::launch::deferred, fun));
			}
		}

		//WriteTotecplot(*m_mocMeshPtr, *m_cfdMeshPtr, vCFD, vSMoc, "debuf.plt");
		for (size_t j = 0; j < futureVec.size(); j++)
			futureVec[j].get();

		if (CFDID % 100 == 0 || CFDID == cfdMesh.GetMeshPointNum())
		{
			Logger::LogInfo(FormatStr("Solver Initialization: %.2lf%% Completed.", CFDID * 100.0 / cfdMesh.GetMeshPointNum()));
		}
	}

	Logger::LogInfo(FormatStr("CFD cell number:%d, and %d of them are located inside one MOC cell, taking %.2lf percent", cfdMesh.GetMeshPointNum(), iNum, 100 * double(iNum) / cfdMesh.GetMeshPointNum()));
	WriteMapInfortoFile();

	DisplayEmptyMap();
}
void Solver::DisplayEmptyMap()
{
	std::stringstream  streamOutCFD,streamOutMOC;
	for (int i = 0; i < m_CFD_MOC_Map.size(); i++)
	{
		if (m_CFD_MOC_Map[i].empty())
		{
			std::stringstream  streamOut;
			streamOut << "empty cfd id:" << i << std::endl;

		}
	}
	Logger::LogInfo(streamOutCFD.str(),true);
	for (int i = 0; i < m_MOC_CFD_Map.size(); i++)
	{
		for (int j = 0; j < m_MOC_CFD_Map[i].size(); j++)
		{
			for (int k = 0; k < m_MOC_CFD_Map[i][j].size(); k++)
			{
				SMocIndex sTemp(i, j, k);
				const MOCMeshPoint& mocPoint = dynamic_cast<const MOCMeshPoint&>(*m_mocMeshPtr->GetMocMeshPointPtr(sTemp));

				if (m_MOC_CFD_Map[i][j][k].empty() && mocPoint.GetMaterialName() == materialName)
				{
					streamOutMOC << "empty moc id :" << i << " " << j << " " << k << std::endl;
				}

			}
		}
	}
	Logger::LogInfo(streamOutMOC.str(),true);
}

void Solver::GetMocIndexByMapValue(std::vector< SMocIndex>& vSMocIndex)
{

	for (int i = 0; i < m_MOC_CFD_Map.size(); i++)
	{
		for (int j = 0; j < m_MOC_CFD_Map[i].size(); j++)
		{
			for (int k = 0; k < m_MOC_CFD_Map[i][j].size(); k++)
			{
				std::unordered_map<int, double>::iterator it;
				if(!m_MOC_CFD_Map[i][j][k].empty())
					vSMocIndex.push_back(SMocIndex(i, j, k));
			}
		}
	}
}
void Solver::CheckMappingWeights()
{
	std::vector< SMocIndex> vSMocIndex;
	m_mocMeshPtr->GetMocIndexByMaterial(vSMocIndex,materialName);

	double totalMOCVolume = 0.0;
	for (int j = 0; j < vSMocIndex.size(); j++)
	{
		const MOCMeshPoint& mocPoint = dynamic_cast<const MOCMeshPoint&>(*m_mocMeshPtr->GetMocMeshPointPtr(vSMocIndex[j]));
		//if (this->materialName != mocPoint.GetMaterialName()) continue;
		totalMOCVolume += mocPoint.Volume();
	}
	double totalCFDVolume = 0.0;
	for (int j = 0; j < m_cfdMeshPtr->GetMeshPointNum(); j++)
	{
		const CFDMeshPoint& cfdPoint = dynamic_cast<const CFDMeshPoint&>(*m_cfdMeshPtr->GetMeshPointPtr(j));
		totalCFDVolume += cfdPoint.Volume();
	}
	Logger::LogInfo(FormatStr("Total cell number of MOC is %d,Total volume of MOC is %.6lf", vSMocIndex.size(),totalMOCVolume));
	Logger::LogInfo(FormatStr("Total cell number of CFD is %d,Total volume of CFD is %.6lf", m_cfdMeshPtr->GetMeshPointNum(), totalCFDVolume));
	//compute the maximum and the minimum of weights of CFD->MOC mapping
	double minSumWeight = 1;
	double maxSumWeight = 0;
	//they are exepected to be around 1


	int sumOfZero = 0;
	for (int j = 0; j < vSMocIndex.size(); j++)
	{
		//const MOCMeshPoint& mocPoint = dynamic_cast<const MOCMeshPoint&>(*m_mocMeshPtr->GetMocMeshPointPtr(vSMocIndex[j]));
		//if (mocPoint.GetMaterialName() != materialName) continue;
		double sumSumWeight = 0.0;
		for (auto& iter : m_MOC_CFD_Map[vSMocIndex[j].iAssemblyIndex][vSMocIndex[j].iCellIndex][vSMocIndex[j].iMocIndex])
		{
			sumSumWeight += iter.second;
		}
		minSumWeight = min(minSumWeight, sumSumWeight);
		maxSumWeight = max(maxSumWeight, sumSumWeight);
		if (sumSumWeight < SMALL)
		{
			std::stringstream streamTemp;
			streamTemp << "0 sumSumWeight moc id: " << vSMocIndex[j].iAssemblyIndex
				<<" "<< vSMocIndex[j].iCellIndex<<" "<< vSMocIndex[j].iMocIndex <<std::endl;
			//Logger::LogInfo(streamTemp.str(),true);
			sumOfZero++;
		}
	}
	Logger::LogInfo(FormatStr("sum of CFD->MOC weights ranges from %.6lf to %6lf", minSumWeight, maxSumWeight));
	Logger::LogInfo(FormatStr("%d MOC cells have no weights from CFD cells", sumOfZero));
	//compute the maximum and the minimum of weights of MOC->CFD mapping
	minSumWeight = 1;
	maxSumWeight = 0;
	//also, they are exepected to be around 1
	sumOfZero = 0;
	std::vector<int> topTen;
	int count = 0;
	for (int j = 0; j < m_cfdMeshPtr->GetMeshPointNum(); j++)
	{
		double sumSumWeight = 0.0;
		for (auto& iter : m_CFD_MOC_Map[j])
		{
			sumSumWeight += iter.second;
		}
		minSumWeight = min(minSumWeight, sumSumWeight);
		maxSumWeight = max(maxSumWeight, sumSumWeight);
		if (sumSumWeight < SMALL)
		{
			Logger::LogInfo(FormatStr("0 sumSumWeight cfd id  %d ", j),true);
			sumOfZero++;
			if (count < 10)
			{
				topTen.push_back(j);
				count++;
			}
		}
	}
	std::cout << "if any, top 10 CFD ID having no overlaping MOC cells are list below:" << std::endl;
	for (int j = 0;j < topTen.size();j++)
	{
		std::cout << topTen[j] << "\t";
	}
	Logger::LogInfo(FormatStr("sum of MOC->CFD weights ranges from %.6lf to %6lf", minSumWeight, maxSumWeight));
	Logger::LogInfo(FormatStr("%d CFD cells have no weights from MOC cells", sumOfZero));
	return;
}
/*
void Solver::Interception(const GeneralMesh* sourceMesh, GeneralMesh* targetMesh, ValueType vt)
{


	std::vector<std::unordered_map<int, double>>* interMap = nullptr;
	if (dynamic_cast<const CFDMesh*>(sourceMesh) && dynamic_cast<MOCMesh*>(targetMesh))
		interMap = &m_MOC_CFD_Map;
	else if (dynamic_cast<const MOCMesh*>(sourceMesh) && dynamic_cast<CFDMesh*>(targetMesh))
		interMap = &m_CFD_MOC_Map;
	else
		throw std::runtime_error("Mesh Cast Error!");

	std::vector<double> targetValueField = targetMesh->GetValueVec(vt);
	std::vector<double> sourceValueField = sourceMesh->GetValueVec(vt);
	for (int i = 0; i < targetMesh->GetMeshPointNum(); i++)
	{
		double interValue = 0.0;
		double selfInterCoe = 1.0;
		for (auto& pair : (*interMap)[i])
		{
			int MOCMeshPointID = pair.first;
			double interCoe = pair.second;
			interValue += interCoe * sourceValueField[MOCMeshPointID];
			selfInterCoe -= interCoe;
		}
		interValue += selfInterCoe * targetValueField[i];
		targetValueField[i] = interValue;
	}

	targetMesh->SetValueVec(targetValueField, vt);
	
	return;
}
*/
void Solver::Interception_fromMocToCFD
(
	MOCMesh& sourceMesh,
	CFDMesh& targetMesh, 
	ValueType vt
)
{

	std::vector<double> targetValueField;
	targetMesh.GetValueVec(vt,targetValueField);
	for (int i = 0; i < targetMesh.GetMeshPointNum(); i++) 
	{
		double interValue = 0.0;
		double selfInterCoe = 1.0;
		double dValue = 0;
		for (auto& pair : m_CFD_MOC_Map[i])
		{
			dValue += pair.second;
		}

		for (auto& pair : m_CFD_MOC_Map[i])
		{
			SMocIndex sMOCMeshPointID = pair.first;
			double interCoe = pair.second/ dValue;
			interValue += interCoe * sourceMesh.GetValueAtIndex(sMOCMeshPointID,vt);
			selfInterCoe -= interCoe;
		}
		if(dValue!=0)
			targetValueField[i] = interValue;
	}

	targetMesh.SetValueVec(targetValueField, vt);

	return;
}

void Solver::Interception_fromCFDToMOC
(
	CFDMesh& sourceMesh,
	MOCMesh& targetMesh,
	ValueType vt
)
{

	//std::vector< SMocIndex>& vSMocIndex = targetMesh.m_vSMocIndex;
	std::vector< SMocIndex> vSMocIndex;
	GetMocIndexByMapValue(vSMocIndex);
	std::vector<double> sourceValueField;
	sourceMesh.GetValueVec(vt, sourceValueField);
	double dMax = 0, dMin = 2;
	for (int i = 0; i < vSMocIndex.size(); i++)
	{
		double interValue = 0.0;
		double dValue= GetMocMeshMapTotalValue(vSMocIndex[i]);
		double iMin = 10000, iMax = 0;

		for (auto& pair : m_MOC_CFD_Map[vSMocIndex[i].iAssemblyIndex][vSMocIndex[i].iCellIndex][vSMocIndex[i].iMocIndex])
		{
			//if (pair.first.first != g_iMpiID)
				//continue;
			int  iCFDMeshPointID = pair.first;
			double interCoe = pair.second/ dValue;
			interValue += interCoe * sourceValueField[iCFDMeshPointID];

		}
		dMax = max(dMax, dValue);
		if(dValue!=0)
			dMin = min(dMin, dValue);
		if(dValue !=0)
			targetMesh.SetValueAtIndex(vSMocIndex[i], interValue, vt);
	}
	Logger::LogInfo(FormatStr("the max dValue:%.6lf the Min:%6.lf ", dMax, dMin));
	return;
}

void Solver::WriteMapInfortoFile()
{
	ofstream CFDtoMOC_MapFile("./temp/MapFile_"+materialName+"_CFDtoMOC"+ "_"+std::to_string(g_iMpiID));
	ofstream MOCtoCFD_MapFile("./temp/MapFile_"+materialName+"_MOCtoCFD"+ "_"+std::to_string(g_iMpiID));
	stringstream MOCtoCFD_MapFile_stream;

	CFDtoMOC_MapFile << m_CFD_MOC_Map.size() << std::endl;
	for (int i = 0; i < m_CFD_MOC_Map.size(); i++)
	{
		std::unordered_map<SMocIndex, double>::iterator it;
		for (it = m_CFD_MOC_Map[i].begin(); it != m_CFD_MOC_Map[i].end(); it++)
		{
			SMocIndex sTemp(it->first.iAssemblyIndex, it->first.iCellIndex, it->first.iMocIndex);
			//std::shared_ptr<MeshPoint>pMeshPoint = m_mocMeshPtr->GetMocMeshPointPtr(sTemp);
			const MOCMeshPoint& mocPoint = dynamic_cast<const MOCMeshPoint&>(*m_mocMeshPtr->GetMocMeshPointPtr(sTemp));
			
			CFDtoMOC_MapFile << i << " " << it->first.iAssemblyIndex
				<< " " << it->first.iCellIndex << " " << it->first.iMocIndex << " " << it->second 
				<<" "<< mocPoint.Volume()<<" "<<mocPoint.GetTemperatureName() 
				<<" "<< mocPoint.PointID() << " " << mocPoint.GetMaterialNameWithID() << std::endl;

		}
	}
	for (int i = 0; i < m_MOC_CFD_Map.size(); i++)
	{
		for (int j = 0; j < m_MOC_CFD_Map[i].size(); j++)
		{
			for (int k = 0; k < m_MOC_CFD_Map[i][j].size(); k++)
			{
				std::unordered_map<int, double>::iterator it;
				for (it = m_MOC_CFD_Map[i][j][k].begin(); it != m_MOC_CFD_Map[i][j][k].end(); it++)
				{
					MOCtoCFD_MapFile << i << " " << j << " " << k << " " << it->first << " " << it->second << std::endl;
					MOCtoCFD_MapFile_stream << i << " " << j << " " << k << " "
						<< g_iMpiID << " " << it->first << " " << it->second << std::endl;
				}
			}
		}
	}
	CFDtoMOC_MapFile.close();
	MOCtoCFD_MapFile.close();
	//MPI_WriteStream_To_File(("MapFile_" + materialName + "_MOCtoCFD"), MOCtoCFD_MapFile_stream);


}

void Solver::ReadMapInfor()
{
	m_CFD_MOC_Map.resize(m_cfdMeshPtr->GetMeshPointNum());
	std::string fileName = "./temp/MapFile_" + materialName + "_CFDtoMOC" + "_" + std::to_string(g_iMpiID);
	ifstream infile(fileName);
	if (!infile.is_open())
	{
		Logger::LogError("cannot find the CFD to MOC map file: " + fileName);
		exit(EXIT_FAILURE);
		return;
	}
	Logger::LogInfo("reading CFD to MOC map file in material: " + materialName);
	std::string line;
	getline(infile, line);
	while (getline(infile, line))
	{
		int i, j, k, m;
		double value;
		stringstream stringline(line);
		stringline >> i >> j >> k >> m >> value;
		SMocIndex sTemp;
		sTemp.iAssemblyIndex = j;
		sTemp.iCellIndex = k;
		sTemp.iMocIndex = m;
		m_CFD_MOC_Map[i][sTemp] = value;

	}
	infile.close();
	int iMax_iAssembly, iMax_iCell, iMax_iMoc;
	std::tuple<int, int, int>tupIndex = GetMaxIndexOfMoc();
	iMax_iAssembly = std::get<0>(tupIndex);
	iMax_iCell = std::get<1>(tupIndex);
	iMax_iMoc = std::get<2>(tupIndex);

	m_MOC_CFD_Map_Total.resize(iMax_iAssembly);
	m_MOC_CFD_Map.resize(iMax_iAssembly);
	for (int i = 0; i < iMax_iAssembly; i++)
	{
		m_MOC_CFD_Map_Total[i].resize(iMax_iCell);
		m_MOC_CFD_Map[i].resize(iMax_iCell);
		for (int j = 0; j < iMax_iCell; j++)
		{
			m_MOC_CFD_Map_Total[i][j].resize(iMax_iMoc);
			m_MOC_CFD_Map[i][j].resize(iMax_iMoc);
		}
	}

	stringstream strTemp;
	fileName = "./temp/MapFile_" + materialName + "_MOCtoCFD";
	MPI_OpenFile_To_Stream(fileName, strTemp);
	Logger::LogInfo("reading MOC to CFD map file in material: " + materialName);
	while (getline(strTemp, line))
	{
		int i, j, k, m, n;
		double dValue;
		stringstream stringline(line);
		stringline >> i >> j >> k >>dValue;
		//std::cout << i << " " << j << " " << k << " " << m << " " << n << " " << dValue << std::endl;
		m_MOC_CFD_Map_Total[i][j][k]=dValue;

	}

	fileName = "./temp/MapFile_" + materialName + "_MOCtoCFD" + "_" + std::to_string(g_iMpiID);;
	infile.open(fileName);
	if (!infile.is_open())
	{
		Logger::LogError("cannot find the MOC to CFD map file:" + fileName);
		exit(EXIT_FAILURE);
		return;
	}
	Logger::LogInfo("reading MOC to CFD_xx map file in material: " + materialName);

	while (getline(infile, line))
	{
		int i, j, k, m;
		double dValue;
		stringstream stringline(line);
		stringline >> i >> j >> k >> m >> dValue;
		//m_MOC_CFD_MapWithID[i][iCell][iMesh].emplace(std::make_pair(i, n), k);
		m_MOC_CFD_Map[i][j][k].emplace(m, dValue);
	}
	infile.close();
	/* for debug
	double dMax = 0, dMin = 2.0;
	for (int i = 0; i < m_MOC_CFD_Map_Total.size(); i++)
	{
		for (int j = 0; j < m_MOC_CFD_Map_Total[i].size(); j++)
		{
			for (int k = 0; k < m_MOC_CFD_Map_Total[i][j].size(); k++)
			{
				double dValue = m_MOC_CFD_Map_Total[i][j][k];

				dMax = max(dMax, dValue);
				dMin = min(dMin, dValue);
			}
		}
	}
	//WriteToLog(FormatStr("the max value is:%.6lf,the min value is:&%.6lf", dMax, dMin));
	*/
	return;
}

