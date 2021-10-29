#include "Solver.h"
#include "Logger.h"
#include <future>
#include <mutex>

#define INTERSECT_JUDGE_LIMIT 1e-10
Solver::Solver(MOCMesh& mocMesh, CFDMesh& cfdMesh) : m_mocMeshPtr(&mocMesh), m_cfdMeshPtr(&cfdMesh) {
	std::mutex mtx;
	m_CFD_MOC_Map.resize(cfdMesh.GetMeshPointNum());
	m_MOC_CFD_Map.resize(mocMesh.GetMeshPointNum());
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
				if (mocPoint.GetMaterialType() == MaterialType::H2O)
					intersectedVolume = cfdPoint.IntersectedVolume(mocPoint);

				if (intersectedVolume > INTERSECT_JUDGE_LIMIT) {
					std::lock_guard<std::mutex> lg(mtx);
					m_CFD_MOC_Map[i][j] = intersectedVolume / cfdPointVolume;
					m_MOC_CFD_Map[j][i] = intersectedVolume / mocPointVolume;
				}
			};
			fun();
			//futureVec.push_back(std::async(std::launch::async, fun));
		}
		for (size_t j = 0; j < futureVec.size(); j++)
			futureVec[j].get();

		if (i % 100 == 0 || i == cfdMesh.GetMeshPointNum())
			Logger::LogInfo(FormatStr("Solver Initialization: %.2lf%% Completed.", i * 100.0 / cfdMesh.GetMeshPointNum()));
	}
}

Solver::Solver(MOCMesh& mocMesh, CFDMesh& cfdMesh,MOCIndex& mocIndex) : m_mocMeshPtr(&mocMesh), m_cfdMeshPtr(&cfdMesh) {
	std::mutex mtx;
	m_CFD_MOC_Map.resize(cfdMesh.GetMeshPointNum());
	m_MOC_CFD_Map.resize(mocMesh.GetMeshPointNum());
	int iNum = 0;
	for (int i = 0; i < cfdMesh.GetMeshPointNum(); i++)
	{

		double x = std::get<0>(m_cfdMeshPtr->GetMeshPointPtr(i)->CentralCoordinate());
		double y = std::get<1>(m_cfdMeshPtr->GetMeshPointPtr(i)->CentralCoordinate());
		double z=  std::get<2>(m_cfdMeshPtr->GetMeshPointPtr(i)->CentralCoordinate());
		int iMocIndex = mocIndex.GetMOCIDWithPoint(x,y,z);

		const CFDMeshPoint& cfdPoint = dynamic_cast<const CFDMeshPoint&>(*m_cfdMeshPtr->GetMeshPointPtr(i));
		const MOCMeshPoint& mocPoint = dynamic_cast<const MOCMeshPoint&>(*m_mocMeshPtr->GetMeshPointPtr(iMocIndex));

		double cfdPointVolume = cfdPoint.Volume();
		double mocPointVolume = mocPoint.Volume();
		double intersectedVolume = 0.0;
		if (mocPoint.GetMaterialType() == MaterialType::H2O)
			intersectedVolume = cfdPoint.IntersectedVolume(mocPoint);


		if ((cfdPointVolume - intersectedVolume) <= INTERSECT_JUDGE_LIMIT)
		{
			m_CFD_MOC_Map[i][iMocIndex] = intersectedVolume / cfdPointVolume;
			m_MOC_CFD_Map[iMocIndex][i] = intersectedVolume / mocPointVolume;

			iNum++;
			

			continue;
		}
		std::vector<std::future<void>> futureVec;
		for (int j = 0; j < mocMesh.GetMeshPointNum(); j++) {
			auto fun = [this, &mtx, i, j]() {
				const CFDMeshPoint& cfdPoint = dynamic_cast<const CFDMeshPoint&>(*m_cfdMeshPtr->GetMeshPointPtr(i));
				const MOCMeshPoint& mocPoint = dynamic_cast<const MOCMeshPoint&>(*m_mocMeshPtr->GetMeshPointPtr(j));

				double cfdPointVolume = cfdPoint.Volume();
				double mocPointVolume = mocPoint.Volume();
				double intersectedVolume = 0.0;
				if (mocPoint.GetMaterialType() == MaterialType::H2O)
					intersectedVolume = cfdPoint.IntersectedVolume(mocPoint);

				if (intersectedVolume > INTERSECT_JUDGE_LIMIT) {
					std::lock_guard<std::mutex> lg(mtx);
					m_CFD_MOC_Map[i][j] = intersectedVolume / cfdPointVolume;
					m_MOC_CFD_Map[j][i] = intersectedVolume / mocPointVolume;
				}
			};
			if (j == iMocIndex)
				continue;
			fun();
			//futureVec.push_back(std::async(std::launch::async, fun));
		}
		for (size_t j = 0; j < futureVec.size(); j++)
			futureVec[j].get();

		if (i % 100 == 0 || i == cfdMesh.GetMeshPointNum())
			Logger::LogInfo(FormatStr("Solver Initialization: %.2lf%% Completed.", i * 100.0 / cfdMesh.GetMeshPointNum()));
		if (/*cfdPoint.PointID()*/i == 2762)
		{
			double value = 0.0;
			for (auto& iter : m_CFD_MOC_Map[i]) {
				value += iter.second;
				Logger::LogInfo(FormatStr("插值权系数:%.6lf", iter.second));
			}
			Logger::LogInfo(FormatStr("CFD网格编号:%d,插值权系数总和:%.6lf/n", i, value));
		}
	}
	Logger::LogInfo(FormatStr("CFD总数量是:%d. 被包含在MOC中的CFD数量是 %d,占比是:百分之%.2lf ", cfdMesh.GetMeshPointNum(), iNum,100*double(iNum)/ cfdMesh.GetMeshPointNum()));
	for (int j = 0; j < mocMesh.GetMeshPointNum(); j++) {
		if (/*cfdPoint.PointID()*/j == 3){
			double value = 0.0;
			for (auto& iter : m_MOC_CFD_Map[j]) {
				value += iter.second;
				Logger::LogInfo(FormatStr("插值权系数:%.6lf", iter.second));
			}
			Logger::LogInfo(FormatStr("MOC网格编号:%d,插值权系数总和:%.6lf/n", j, value));
			//std::cout << value << std::endl;
		}
	}
}
/*
Solver::Solver(MOCMesh& mocMesh, CFDMesh& cfdMesh, MOCIndex& mocIndex)
	: m_mocMeshPtr(&mocMesh), m_cfdMeshPtr(&cfdMesh) 
{
	std::cout << mocIndex.GetMOCIDWithPoint(1.0, 0.5, 0.25) << std::endl;
	std::cout << mocIndex.GetMOCIDWithPoint(1.0, 0.2, 0.25) << std::endl;
	std::cout << mocIndex.GetMOCIDWithPoint(0.1, 0.2, 0.25) << std::endl;
	std::cout << mocIndex.GetMOCIDWithPoint(0.5, 0.26, 0.25) << std::endl;
	std::cout << mocIndex.GetMOCIDWithPoint(0.25, 0.75, 0.4) << std::endl;
}
*/
void Solver::Interception(const Mesh* sourceMesh, Mesh* targetMesh, ValueType vt) {
	std::vector<std::unordered_map<int, double>>* interMap = nullptr;
	if (dynamic_cast<const CFDMesh*>(sourceMesh) && dynamic_cast<MOCMesh*>(targetMesh))
		interMap = &m_MOC_CFD_Map;
	else if (dynamic_cast<const MOCMesh*>(sourceMesh) && dynamic_cast<CFDMesh*>(targetMesh))
		interMap = &m_CFD_MOC_Map;
	else
		throw std::runtime_error("Mesh Cast Error!");

	std::vector<double> targetValueField = targetMesh->GetValueVec(vt);
	std::vector<double> sourceValueField = sourceMesh->GetValueVec(vt);

	for (int i = 0; i < targetMesh->GetMeshPointNum(); i++) {
		double interValue = 0.0;
		double selfInterCoe = 1.0;
		for (auto& pair : (*interMap)[i]) {
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