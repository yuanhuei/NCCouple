#include "MOCIndex.h"
#include "Logger.h"
#include "./MHT_common/SystemControl.h"
#include "MOCMesh.h"

MOCIndex::MOCIndex(MOCMesh& mocMesh)
	:pMOCMesh(&mocMesh)
{
	Initialization();
	//BuildUpIndex();
	//CheckIndex();
}
MOCIndex::MOCIndex(MOCMesh&pMocMesh, Cell& pCell) :m_pCell(&pCell), pMOCMesh(&pMocMesh)
{ 
	Initialization();
}


void MOCIndex::Initialization()
{
	this->SetTolerance();
	this->SetAxialInfo();
	this->SetCircularInfo();
	this->GetRadiusList();
	return;
}

void MOCIndex::BuildUpIndex()
{
	AllocateDim(this->v_MOCID, this->circularCellNum, this->v_radius.size(), this->axialCellNum);
	//for (int i = 0; i < pMOCMesh->GetMeshPointNum(); i++)
	for (int i = 0; i < m_pCell->vMeshPointPtrVec.size(); i++)

	{
		const MOCMeshPoint& mocPoint = dynamic_cast<const MOCMeshPoint&>(*m_pCell->vMeshPointPtrVec[i]);
		int verticeNum = mocPoint.VerticesNum();
		double RSum = 0.0;
		double ZSum = 0.0;
		double COSTheetaSum = 0.0;
		Vector OPSum(0.0, 0.0, 0.0);
		double count1 = 0.0;
		double count2 = 0.0;
		for (int j = 0; j < verticeNum; j++)
		{
			Vector P = mocPoint.VerticeCoordinate(j);
			Vector OP = P - this->axisPoint;
			Vector axialProjection = (OP & this->axisNorm) * this->axisNorm;
			Vector radialProjection = OP - axialProjection;
			RSum += radialProjection.Mag();
			ZSum += axialProjection.Mag();
			OPSum = OPSum + OP;
			count1 += 1.0;
			if (radialProjection.Mag() > scaleTolerance)
			{
				COSTheetaSum += radialProjection.GetNormal() & this->theetaStartNorm;
				count2 += 1.0;
			}
		}
		Scalar radius = RSum / count1;
		Scalar height = ZSum / count1;
		Scalar cosTheeta = COSTheetaSum / count2;
		Scalar theeta = acos(cosTheeta);
		Vector OP = OPSum / count1;
		if (((theetaStartNorm ^ OP) & axisNorm) < 0.0)
		{
			theeta = 2.0 * PI - theeta;
		}
		int IndexI = int((double)circularCellNum * theeta / (2.0 * PI));
		int IndexJ = 0;
		for (int j = 0; j < this->v_radius.size(); j++)
		{
			if (v_radius[j] > radius)
			{
				break;
			}
			IndexJ = j;
		}
		int IndexK = int(height / axialCellSize);
		this->v_MOCID[IndexI][IndexJ][IndexK] = i;
	}
	//CheckIndex();
	return;
}

//find the largest radius to estimate a tolerance
void MOCIndex::SetTolerance()
{
	Scalar maxRadius = 0.0;
	for (int i = 0; i < m_pCell->vMeshPointPtrVec.size(); i++)
	{
		const MOCMeshPoint& mocPoint = dynamic_cast<const MOCMeshPoint&>(*m_pCell->vMeshPointPtrVec[i]);
		std::vector<Scalar> radiusList = mocPoint.GetRadiusList();
		for (size_t j = 0;j < radiusList.size();j++)
		{
			maxRadius = max(maxRadius, radiusList[j]);
		}
	}
	this->scaleTolerance = 1e-4 * maxRadius;
	return;
}

void MOCIndex::GetRadiusList()
{
	//insert into v_radius list
	this->v_radius.clear();
	for (int i = 0; i < m_pCell->vMeshPointPtrVec.size(); i++)
	{
		const MOCMeshPoint& mocPoint = dynamic_cast<const MOCMeshPoint&>(*m_pCell->vMeshPointPtrVec[i]);
		std::vector<Scalar> radiusList = mocPoint.GetRadiusList();
		//for each radius in this cell
		for (size_t j = 0;j < radiusList.size();j++)
		{
			//loop over the existing vadius list,
			bool found = false;
			for (size_t k = 0;k < v_radius.size();k++)
			{
				if (fabs(v_radius[k] - radiusList[j]) < scaleTolerance)
				{
					found = true;
					break;
				}
			}
			//insert if not found within given tolerance
			if (false == found)
			{
				v_radius.push_back(radiusList[j]);
			}
		}
	}
	v_radius.push_back(0.0);
	//re-order v_radius
	int vadiusNum = v_radius.size();
	for (size_t i = 0;i < vadiusNum;i++)
	{
		for (int j = 0;j < vadiusNum - 1;j++)
		{
			if (v_radius[j] > v_radius[j + 1])
			{
				Scalar temp = v_radius[j];
				v_radius[j] = v_radius[j + 1];
				v_radius[j + 1] = temp;
			}
		}
	}
	return;
}

void MOCIndex::SetAxialInfo()
{
	//axial cell number and cell size
	std::pair<int, Scalar> info = pMOCMesh->GetAxialInformation();
	this->axialCellNum = info.first;
	this->axialCellSize = info.second;
	//axis center
	Vector numerator(0.0,0.0,0.0);
	int denominator = 0;
	for (int i = 0; i < m_pCell->vMeshPointPtrVec.size(); i++)
	{
		const MOCMeshPoint& mocPoint = dynamic_cast<const MOCMeshPoint&>(*m_pCell->vMeshPointPtrVec[i]);
		std::pair<bool, Vector> axisCenterInfo = mocPoint.AxisCenter();
		if (axisCenterInfo.first)
		{
			numerator += axisCenterInfo.second;
			denominator++;
		}
	}
	if (0 == denominator)
	{
		Logger::LogError("in MOCIndex::SetAxialInfo(), no curved cell is found");
	}
	this->axisPoint = numerator / Scalar(denominator);
	this->axisPoint.z_ = 0.0;
	//axis norm
	this->axisNorm = Vector(0.0, 0.0, 1.0);
	return;
}

void MOCIndex::SetCircularInfo()
{
	this->theetaStartNorm = Vector(1.0, 0.0, 0.0);
	this->circularCellNum = 8;
	return;
}

std::tuple<Scalar, Scalar, Scalar> MOCIndex::GetCylindricalCoordinate
(
	Scalar x,
	Scalar y,
	Scalar z
)
{
	Vector P(x, y, z);
	Vector OP = P - this->axisPoint;
	Vector axialProjection = (OP & this->axisNorm) * this->axisNorm;
	Vector radialProjection = OP - axialProjection;
	Scalar radius = radialProjection.Mag();
	Scalar height = axialProjection.Mag();
	if (radialProjection.Mag() < SMALL)
	{
		return std::tuple<Scalar, Scalar, Scalar>(0.0, radius, height);
	}
	Scalar cosTheeta = radialProjection.GetNormal() & this->theetaStartNorm;
	Scalar theeta = 0.0;
	//it can be avoid theeta resulted as 2PI when cosTheeta = 1.0
	if (fabs(cosTheeta - 1.0) > 10.0 * SMALL)
	{
		theeta = acos(cosTheeta);
		if (((theetaStartNorm ^ radialProjection) & axisNorm) < 0.0)
		{
			theeta = 2.0 * PI - theeta;
		}
	}
	return std::tuple<Scalar, Scalar, Scalar>(theeta, radius, height);
}


std::tuple<int, int, int> MOCIndex::GetIJKWithPoint
(
	Scalar x,
	Scalar y,
	Scalar z
)
{
	std::tuple<Scalar,Scalar,Scalar> theetaRZ = GetCylindricalCoordinate(x, y, z);
	int circularIndex = int((double)circularCellNum * std::get<0>(theetaRZ) / (2.0 * PI));
	int radialIndex = 0;
	for (int j = 0; j < this->v_radius.size();j++)
	{
		if (v_radius[j] > std::get<1>(theetaRZ))
		{
			break;
		}
		radialIndex = j;
	}
	int axialIndex = int(std::get<2>(theetaRZ) / axialCellSize);
	return std::tuple<int, int, int>(circularIndex, radialIndex, axialIndex);
}

int MOCIndex::GetMOCIDWithPoint
(
	Scalar x,
	Scalar y,
	Scalar z
)
{
	std::tuple<int, int, int> indexes = GetIJKWithPoint(x,y,z);
	int i = std::get<0>(indexes);
	int j = std::get<1>(indexes);
	int k = std::get<2>(indexes);
	return this->v_MOCID[i][j][k];
}


void MOCIndex::CheckIndex()
{
	int meshNum = m_pCell->vMeshPointPtrVec.size();
	std::vector<bool> v_registered;
	v_registered.resize(meshNum);
	for (size_t i = 0;i < meshNum;i++)
	{
		v_registered[i] = false;
	}

	for (int i = 0;i < this->v_MOCID.size();i++)
	{
		for (int j = 0;j < this->v_MOCID[i].size();j++)
		{
			for (int k = 0;k < this->v_MOCID[i][j].size();k++)
			{
				int ID = this->v_MOCID[i][j][k];
				v_registered[ID] = true;
			}
		}
	}
	std::stringstream msg;
	int unregisteredNum = 0;
	for (int i = 0; i < meshNum; i++)
	{
		if (false == v_registered[i])
		{
			const MOCMeshPoint& mocPoint = dynamic_cast<const MOCMeshPoint&>(*m_pCell->vMeshPointPtrVec[i]);
			Vector meshCenter = mocPoint.Center();
			msg << "mesh #" << i << " is not registered, centering at " << meshCenter << std::endl;
			unregisteredNum++;
		}
	}
	if (0 != unregisteredNum)
	{
		FatalError(msg.str());
	}
	return;
}
void MOCIndex::Display()
{
	std::cout << "axisPoint = " << axisPoint << std::endl;
	std::cout << "axisNorm = " << axisNorm << std::endl;
	std::cout << "theetaStartNorm = " << theetaStartNorm << std::endl;
	std::cout << "circularCellNum = " << circularCellNum << std::endl;
	std::cout << "axialCellSize = " << axialCellSize << std::endl;
	std::cout << "axialCellNum = " << axialCellNum << std::endl;
	std::cout << "scaleTolerance = " << scaleTolerance << std::endl;
	for (int i = 0; i < this->v_radius.size(); i++)
	{
		std::cout << "v_radius[" << i << "] = " << v_radius[i] << std::endl;
	}
	for (int k = 0; k < axialCellNum; k++)
	{
		for (int i = 0; i < this->v_MOCID.size(); i++)
		{
			for (int j = 0; j < this->v_MOCID[i].size(); j++)
			{
				std::cout << i << "\t" << j << "\t" << k << "\t" << this->v_MOCID[i][j][k] << std::endl;
			}
		}
	}
	return;
}
