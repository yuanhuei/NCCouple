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
	writeMapInfortoFile();
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
		readMapInfor();
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
		Vector P = m_cfdMeshPtr->GetMeshPointPtr(CFDID)->Center();
		SMocIndex sFirstMocIndex = mocMesh.getIndex(P);
		//Logger::LogError("can not find  moc index corresponding to cfd Coordinate");
		if (sFirstMocIndex == SMocIndex(-1, -1, -1)) 
		{
			std::cout << "The coordinate of CFD mesh is " << P << std::endl;
			Logger::LogError("can not find  moc index corresponding to cfd Coordinate");
			continue;
		}
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
		for (size_t j = 0; j < futureVec.size(); j++)
			futureVec[j].get();

		if (CFDID % 100 == 0 || CFDID == cfdMesh.GetMeshPointNum())
		{
			Logger::LogInfo(FormatStr("Solver Initialization: %.2lf%% Completed.", CFDID * 100.0 / cfdMesh.GetMeshPointNum()));
		}
	}

	Logger::LogInfo(FormatStr("CFD cell number:%d, and %d of them are located inside one MOC cell, taking %.2lf percent", cfdMesh.GetMeshPointNum(), iNum, 100 * double(iNum) / cfdMesh.GetMeshPointNum()));
	writeMapInfortoFile();
}
/*
Solver::Solver
(
	MOCMesh& mocMesh,
	CFDMesh& cfdMesh,
	MOCIndex& mocIndex,
	std::string mName
)
	:
	m_mocMeshPtr(&mocMesh),
	m_cfdMeshPtr(&cfdMesh),
	materialName(mName)
{
	std::mutex mtx;
	m_CFD_MOC_Map.resize(cfdMesh.GetMeshPointNum());
	m_MOC_CFD_Map.resize(mocMesh.GetMeshPointNum());
	int iNum = 0;
	//the code below was written for test
	int nCFDNum = cfdMesh.GetMeshPointNum();

	for (int CFDID = 0; CFDID < nCFDNum; CFDID++)
	{
		Vector P = m_cfdMeshPtr->GetMeshPointPtr(CFDID)->Center();
		int iMocIndex = mocIndex.GetMOCIDWithPoint(P.x_, P.y_, P.z_);
		const CFDMeshPoint& cfdPoint = dynamic_cast<const CFDMeshPoint&>(*m_cfdMeshPtr->GetMeshPointPtr(CFDID));
		const MOCMeshPoint& mocPoint = dynamic_cast<const MOCMeshPoint&>(*m_mocMeshPtr->GetMeshPointPtr(iMocIndex));
		double cfdPointVolume = cfdPoint.Volume();
		double mocPointVolume = mocPoint.Volume();
		double intersectedVolume = 0.0;
		if (mocPoint.GetMaterialName() == materialName)
			intersectedVolume = cfdPoint.IntersectedVolume(mocPoint);

		if (intersectedVolume > INTERSECT_JUDGE_LIMIT) {

			m_CFD_MOC_Map[CFDID][iMocIndex] = intersectedVolume / cfdPointVolume;
			m_MOC_CFD_Map[iMocIndex][CFDID] = intersectedVolume / mocPointVolume;
		}
		if ((cfdPointVolume - intersectedVolume) <= INTERSECT_JUDGE_LIMIT)
		{
			iNum++;
			continue;
		}
		std::vector<std::future<void>> futureVec;

		//compute the range of structured index in the axial direction 
		int nzMin = mocIndex.axialCellNum - 1;
		int nzMax = 0;
		int verticeNum = cfdPoint.VerticesNum();
		for (int j = 0; j < verticeNum; j++)
		{
			Vector cor = cfdPoint.VerticeCoordinate(j);
			std::tuple<int, int, int> structuredIJK = mocIndex.GetIJKWithPoint(cor.x_, cor.y_, cor.z_);
			int nz = std::get<2>(structuredIJK);
			nzMin = min(nz, nzMin);
			nzMax = max(nz, nzMax);
		}
		//loop over only a part of MOC cells in the range previously obtained
		for (int kk = max(0, nzMin); kk <= min(mocIndex.axialCellNum - 1, nzMax); kk++)
		{
			for (int ii = 0; ii < mocIndex.v_MOCID.size(); ii++)
			{
				for (int jj = 0; jj < mocIndex.v_MOCID[ii].size(); jj++)
				{
					int MOCID = mocIndex.v_MOCID[ii][jj][kk];
					const MOCMeshPoint& localMOCPoint = dynamic_cast<const MOCMeshPoint&>(*m_mocMeshPtr->GetMeshPointPtr(MOCID));
					if (materialName != localMOCPoint.GetMaterialName()) continue;
					auto fun = [this, &mtx, CFDID, MOCID]()
					{
						const CFDMeshPoint& cfdPoint = dynamic_cast<const CFDMeshPoint&>(*m_cfdMeshPtr->GetMeshPointPtr(CFDID));
						const MOCMeshPoint& mocPoint = dynamic_cast<const MOCMeshPoint&>(*m_mocMeshPtr->GetMeshPointPtr(MOCID));
						double cfdPointVolume = cfdPoint.Volume();
						double mocPointVolume = mocPoint.Volume();
						double intersectedVolume = 0.0;
						intersectedVolume = cfdPoint.IntersectedVolume(mocPoint);
						if (intersectedVolume > INTERSECT_JUDGE_LIMIT)
						{
							std::lock_guard<std::mutex> lg(mtx);
							m_CFD_MOC_Map[CFDID][MOCID] = intersectedVolume / cfdPointVolume;
							m_MOC_CFD_Map[MOCID][CFDID] = intersectedVolume / mocPointVolume;
						}
					};
					if (MOCID == iMocIndex) continue;
					futureVec.push_back(std::async(std::launch::async | std::launch::deferred, fun));
				}
			}
		}

		for (size_t j = 0; j < futureVec.size(); j++)
			futureVec[j].get();

		if (CFDID % 100 == 0 || CFDID == cfdMesh.GetMeshPointNum())
		{
			Logger::LogInfo(FormatStr("Solver Initialization: %.2lf%% Completed.", CFDID * 100.0 / cfdMesh.GetMeshPointNum()));
		}
	}

	Logger::LogInfo(FormatStr("CFD cell number:%d, and %d of them are located inside one MOC cell, taking %.2lf percent", cfdMesh.GetMeshPointNum(), iNum, 100 * double(iNum) / cfdMesh.GetMeshPointNum()));
	writeMapInfortoFile();
}
*/
void Solver::GetMocIndexByMapValue(std::vector< SMocIndex>& vSMocIndex)
{
	for (int i = 0; i < m_MOC_CFD_Map.size(); i++)
	{
		for (int j = 0; j < m_MOC_CFD_Map[i].size(); j++)
		{
			for (int k = 0; k < m_MOC_CFD_Map[i][j].size(); k++)
			{
				std::unordered_map<int, double>::iterator it;
				for (it = m_MOC_CFD_Map[i][j][k].begin(); it != m_MOC_CFD_Map[i][j][k].end(); it++)
				{
					vSMocIndex.push_back(SMocIndex(i, j, k));
					break;
				}
			}
		}
	}
}
void Solver::CheckMappingWeights()
{
	std::vector< SMocIndex> vSMocIndex;
	GetMocIndexByMapValue(vSMocIndex);

	double totalMOCVolume = 0.0;
	for (int j = 0; j < vSMocIndex.size(); j++)
	{
		const MOCMeshPoint& mocPoint = dynamic_cast<const MOCMeshPoint&>(*m_mocMeshPtr->GetMocMeshPointPtr(vSMocIndex[j]));
		if (this->materialName != mocPoint.GetMaterialName()) continue;
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
		const MOCMeshPoint& mocPoint = dynamic_cast<const MOCMeshPoint&>(*m_mocMeshPtr->GetMocMeshPointPtr(vSMocIndex[j]));
		if (mocPoint.GetMaterialName() != materialName) continue;
		double sumSumWeight = 0.0;
		for (auto& iter : m_MOC_CFD_Map[vSMocIndex[j].iAssemblyIndex][vSMocIndex[j].iCellIndex][vSMocIndex[j].iMocIndex])
		{
			sumSumWeight += iter.second;
		}
		minSumWeight = min(minSumWeight, sumSumWeight);
		maxSumWeight = max(maxSumWeight, sumSumWeight);
		if (sumSumWeight < SMALL)
		{
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

void Solver::Interception(const GeneralMesh* sourceMesh, GeneralMesh* targetMesh, ValueType vt)
{


	std::vector<std::unordered_map<int, double>>* interMap = nullptr;
	/*
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
	*/
	return;
}

void Solver::Interception_fromMocToCFD
(
	MOCMesh& sourceMesh,
	CFDMesh& targetMesh, 
	ValueType vt
)
{
	//std::vector< SMocIndex> vSMocIndex;
	//sourceMesh.GetAllMocIndex(vSMocIndex);

	std::vector<double> targetValueField = targetMesh.GetValueVec(vt);
	//std::vector<double> sourceValueField = sourceMesh->GetValueVec(vt);
	for (int i = 0; i < targetMesh.GetMeshPointNum(); i++) 
	{
		double interValue = 0.0;
		double selfInterCoe = 1.0;
		for (auto& pair : m_CFD_MOC_Map[i])
		{
			SMocIndex sMOCMeshPointID = pair.first;
			//MocMeshField* pMocField= sourceMesh.GetFieldPointerAtIndex(sMOCMeshPointID);
			double interCoe = pair.second;
			interValue += interCoe * sourceMesh.GetValueAtIndex(sMOCMeshPointID,vt);
			selfInterCoe -= interCoe;
		}
		interValue += selfInterCoe * targetValueField[i];
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

	std::vector< SMocIndex> vSMocIndex = targetMesh.m_vSMocIndex;
	//targetMesh.GetAllMocIndex(vSMocIndex);

	//std::vector<double> targetValueField = targetMesh.GetValueVec(vt);
	std::vector<double> sourceValueField = sourceMesh.GetValueVec(vt);
	for (int i = 0; i < vSMocIndex.size(); i++)
	{
		//MocMeshField mocField = *targetMesh.GetFieldPointerAtIndex(vSMocIndex[i]);
		double interValue = 0.0;
		double selfInterCoe = 1.0;
		for (auto& pair : m_MOC_CFD_Map[vSMocIndex[i].iAssemblyIndex][vSMocIndex[i].iCellIndex][vSMocIndex[i].iMocIndex])
		{
			int  iCFDMeshPointID = pair.first;
			//MocMeshField* pMocField = sourceMesh.GetValueAtIndex(sMOCMeshPointID);
			double interCoe = pair.second;
			interValue += interCoe * sourceValueField[iCFDMeshPointID];
			selfInterCoe -= interCoe;
		}
		interValue += selfInterCoe * targetMesh.GetValueAtIndex(vSMocIndex[i],vt);
		//mocField.SetValue(interValue,vt);
		targetMesh.SetValueAtIndex(vSMocIndex[i], interValue, vt);
	}

	//targetMesh.SetValueVec(targetValueField, vt);

	return;
}

void Solver::writeMapInfortoFile()
{
	ofstream CFDtoMOC_MapFile("MapFile_"+materialName+"_CFDtoMOC");
	ofstream MOCtoCFD_MapFile("MapFile_"+materialName+"_MOCtoCFD");

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
					MOCtoCFD_MapFile << i << " "<< j <<" "<< k <<" " << it->first << " " << it->second << std::endl;
			}
		}
	}
	CFDtoMOC_MapFile.close();
	MOCtoCFD_MapFile.close();

}

void Solver::readMapInfor()
{
	m_CFD_MOC_Map.resize(m_cfdMeshPtr->GetMeshPointNum());
	/*
	m_MOC_CFD_Map.resize(m_mocMeshPtr->m_vAssembly.size());
	for (int i = 0; i < m_MOC_CFD_Map.size(); i++)
	{
		m_MOC_CFD_Map[i].resize(m_mocMeshPtr->m_vAssembly[i].pAssembly_type->v_Cell.size());
		for (int j = 0; j < m_MOC_CFD_Map[i].size(); j++)
		{
			m_MOC_CFD_Map[i][j].resize(m_mocMeshPtr->m_vAssembly[i].pAssembly_type->v_Cell[j].vMeshPointPtrVec.size());
		}
	}
	*/
	m_MOC_CFD_Map.resize(m_mocMeshPtr->m_vAssemblyField.size());
	for (int i = 0; i < m_MOC_CFD_Map.size(); i++)
	{
		m_MOC_CFD_Map[i].resize(m_mocMeshPtr->m_vAssemblyField[i].size());
		for (int j = 0; j < m_MOC_CFD_Map[i].size(); j++)
		{
			m_MOC_CFD_Map[i][j].resize(m_mocMeshPtr->m_vAssemblyField[i][j].size());
		}
	}

	std::string fileName = "MapFile_" + materialName + "_CFDtoMOC";
	ifstream infile(fileName);
	if (!infile.is_open())
	{
		Logger::LogError("cannot find the CFD to MOC map file: "+ fileName);
		exit(EXIT_FAILURE);
		return;
	}
	Logger::LogInfo("reading CFD to MOC map file in material: " + materialName);
	std::string line;
	while (getline(infile, line))
	{
		int i, j,k,m;
		double value;
		stringstream stringline(line);
		stringline >> i >> j >> k>>m >> value ;
		SMocIndex sTemp;
		sTemp.iAssemblyIndex = j;
		sTemp.iCellIndex = k;
		sTemp.iMocIndex = m;
		m_CFD_MOC_Map[i][sTemp] = value;

	}
	infile.close();
	fileName = "MapFile_"  + materialName + "_MOCtoCFD";
	infile.open(fileName);
	if (!infile.is_open())
	{
		Logger::LogError("cannot find the MOC to CFD map file:" + fileName);
		exit(EXIT_FAILURE);
		return;
	}
	Logger::LogInfo("reading MOC to CFD map file in material: " + materialName);
	while (getline(infile, line))
	{
		int i, j,m,n;
		double k;
		stringstream stringline(line);
		stringline >> i >> j >>m >>n >> k;
		m_MOC_CFD_Map[i][j][m][n] = k;
	}
	infile.close();
	return;
}

/*
Solver::Solver(MOCMesh& mocMesh, CFDMesh& cfdMesh, std::string mName)
	:
	m_mocMeshPtr(&mocMesh),
	m_cfdMeshPtr(&cfdMesh),
	materialName(mName)
{
	readMapInfor();
}
*/