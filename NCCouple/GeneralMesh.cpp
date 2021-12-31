#include "GeneralMesh.h"
#include <fstream>

double MHTMeshPoint::IntersectedVolume(const MeshPoint& other) const {
	const MHTMeshPoint* p_other = dynamic_cast<const MHTMeshPoint*>(&other);
	if (!p_other)
		throw std::runtime_error("the type of input MeshPoint is not same with current MeshPoint!");

	return m_poly.IntersectionVolumeWithPolyhedronSet(p_other->m_poly);
}