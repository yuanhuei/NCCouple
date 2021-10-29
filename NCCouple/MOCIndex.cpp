#include "MOCIndex.h"
#include"Logger.h"

MOCIndex::MOCIndex(MOCMesh& mocMesh)
	:pMOCMesh(&mocMesh)
{
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
		for (int j = 0;j < verticeNum; j++)
		{
			std::tuple<double, double, double> cor = mocPoint.VerticeCoordinate(j);
			Vector P(std::get<0>(cor), std::get<1>(cor), std::get<2>(cor));
			Vector OP = P - this->axisPoint;
			Vector axialProjection = (OP & this->axisNorm) * this->axisNorm;
			Vector radialProjection = OP - axialProjection;
			RSum += radialProjection.Mag();
			ZSum += axialProjection.Mag();
			OPSum = OPSum + OP;
			count1 += 1.0;
			if (radialProjection.Mag() > SMALL)
			{
				COSTheetaSum += radialProjection.GetNormal()& this->theetaStartNorm;
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
		for (int j = 0; j < this->v_radius.size();j++)
		{
			if (v_radius[j] > radius)
			{
				break;
			}
			IndexJ = j;
		}
		int IndexK = int(height / axialCellSize);
		this->v_MOCID[IndexI][IndexJ][IndexK] = i;
		//Logger::LogInfo(FormatStr("IndexI,IndexJ,IndexK is :%d,%d,%d; MocIndex is %d :", IndexI, IndexJ, IndexK, i));
	}
	return;
}

void MOCIndex::SetRadial
(
	std::vector<Scalar>& radiusList
)
{
	this->v_radius.resize(radiusList.size() + 1);
	this->v_radius[0] = 0.0;
	for (int i = 1;i < v_radius.size();i++)
	{
		v_radius[i] = radiusList[i - 1];
	}
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
	Scalar theeta = acos(cosTheeta);
	if (((theetaStartNorm ^ radialProjection) & axisNorm) < 0.0)
	{
		theeta = 2.0 * PI - theeta;
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