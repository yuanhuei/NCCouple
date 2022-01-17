#include "MOCIndex.h"
#include "Logger.h"
#include "./MHT_common/SystemControl.h"

MOCIndex::MOCIndex(MOCMesh& mocMesh)
	:pMOCMesh(&mocMesh)
{}

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
	for (int i = 0; i < pMOCMesh->GetMeshPointNum(); i++)
	{
		const MOCMeshPoint& mocPoint = dynamic_cast<const MOCMeshPoint&>(*pMOCMesh->GetMeshPointPtr(i));
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
	return;
}
//根据燃料堆左下角坐标，栅元中心坐标返回x和y向索引
std::pair<int,int> MOCIndex::getIndex(std::tuple<Vector, double, int> tup_Zone,Vector vMeshpoint)
{

	Vector vStack_LeftDownPoint = std::get<0>(tup_Zone);
	int N_N= std::get<2>(tup_Zone);
	int iLength = std::get<1>(tup_Zone);;
	double xLength = vMeshpoint.x_ - vStack_LeftDownPoint.x_;
	double yLength = vMeshpoint.y_ - vStack_LeftDownPoint.y_;
	int xIndex = xLength / (iLength / N_N);
	int yIndex = yLength / (iLength / N_N);
	return std::pair<int, int>(xIndex,yIndex);

}
//输入区域左下角左边，
void MOCIndex::BuildUpIndex_(Vector vLeftDown, int N_N,double iLength)
{

	m_Zone = std::make_tuple(vLeftDown, iLength, N_N);
	v_MOC_Mesh_ID.resize(pMOCMesh->GetMeshNum());
	for(int k = 0; k < pMOCMesh->GetMeshNum(); k++)
	{ 
		AllocateDim(this->v_MOC_Mesh_ID[k], this->circularCellNum, this->v_radius.size(), this->axialCellNum);

		meshPoint& mesh = pMOCMesh->m_meshPoint[k];

		//堆所在位置索引
		std::pair<int, int> iStackIndex = getIndex(std::make_tuple(vLeftDown, iLength,N_N), (mesh.vStack_LeftDownPoint + mesh.vStack_RightUpPoint) / 2);

        //栅元所在燃料堆中索引
		std::pair<int, int> iMeshIndex = getIndex(std::make_tuple(mesh.vStack_LeftDownPoint,
			(mesh.vStack_RightUpPoint.x_-mesh.vStack_LeftDownPoint.x_), mesh.iN_N),mesh.vPoint);

		v_Stack_Mesh_index[iStackIndex.first][iStackIndex.second][iMeshIndex.first][iMeshIndex.second] = k;

		m_stackMap[std::make_pair(iStackIndex.first, iStackIndex.second)] = std::make_tuple(mesh.vStack_LeftDownPoint, (mesh.vStack_RightUpPoint.x_ - mesh.vStack_LeftDownPoint.x_) / 2, mesh.iN_N);
		for (int i = 0; i < mesh.meshPointPtrVec.size(); i++)
		{
			const MOCMeshPoint& mocPoint = dynamic_cast<const MOCMeshPoint&>(*pMOCMesh->m_meshPoint[k].meshPointPtrVec[i]);
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
			this->v_MOC_Mesh_ID[k][IndexI][IndexJ][IndexK] = i;
		}
	}
	return;
}


//find the largest radius to estimate a tolerance
void MOCIndex::SetTolerance()
{
	Scalar maxRadius = 0.0;
	for (int i = 0; i < pMOCMesh->GetMeshPointNum(); i++)
	{
		const MOCMeshPoint& mocPoint = dynamic_cast<const MOCMeshPoint&>(*pMOCMesh->GetMeshPointPtr(i));
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
	for (int i = 0; i < pMOCMesh->GetMeshPointNum(); i++)
	{
		const MOCMeshPoint& mocPoint = dynamic_cast<const MOCMeshPoint&>(*pMOCMesh->GetMeshPointPtr(i));
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
	for (int i = 0; i < pMOCMesh->GetMeshPointNum(); i++)
	{
		const MOCMeshPoint& mocPoint = dynamic_cast<const MOCMeshPoint&>(*pMOCMesh->GetMeshPointPtr(i));
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

std::pair<int, int> MOCIndex::GetMOCIDWithPoint_(Scalar x, Scalar y, Scalar z)
{
	
	//输入点所在的堆的位置索引
	std::pair<int, int> iStackIndex = getIndex(m_Zone, Vector(x,y,z));

	//m_stackMap[iStackIndex]
	//栅元所在燃料堆中索引
	std::pair<int, int> iMeshIndex = getIndex(m_stackMap[iStackIndex], Vector(x, y, z));
	//栅元ID
	int index=v_Stack_Mesh_index[iStackIndex.first][iStackIndex.second][iMeshIndex.first][iMeshIndex.second];

	std::tuple<int, int, int> indexes = GetIJKWithPoint(x, y, z);
	int i = std::get<0>(indexes);
	int j = std::get<1>(indexes);
	int k = std::get<2>(indexes);
	return std::make_pair(index,this->v_MOC_Mesh_ID[index][i][j][k]);


}

void MOCIndex::CheckIndex()
{
	int meshNum = pMOCMesh->GetMeshPointNum();
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
			const MOCMeshPoint& mocPoint = dynamic_cast<const MOCMeshPoint&>(*pMOCMesh->GetMeshPointPtr(i));
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