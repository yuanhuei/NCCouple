#include "GeneralMesh.h"
#include "./MHT_common/SystemControl.h"
#include <fstream>

std::string NameOfValueType(ValueType vt)
{
	if (vt == ValueType::DENSITY)
	{
		return "density";
	}
	else if (vt == ValueType::HEATPOWER)
	{
		return "heat power";
	}
	else if (vt == ValueType::TEMPERAURE)
	{
		return "temperature";
	}
	else
	{
		FatalError("false ValueType given in NameOfValueType()");
		return "";
	}
}

double MHTMeshPoint::IntersectedVolume(const MeshPoint& other) const {
	const MHTMeshPoint* p_other = dynamic_cast<const MHTMeshPoint*>(&other);
	if (!p_other)
		throw std::runtime_error("the type of input MeshPoint is not same with current MeshPoint!");

	return m_poly.IntersectionVolumeWithPolyhedronSet(p_other->m_poly);
}