#include "Mesh.h"
#include <fstream>

void MeshPoint::Init() {
	CGAL::Polygon_mesh_processing::triangulate_faces(m_poly);

	//m_poly = Nef_polyhedron(poly);
	m_volume = CGAL::to_double(CGAL::Polygon_mesh_processing::volume(m_poly));
	m_centerPoint = CGAL::centroid(m_poly.points_begin(), m_poly.points_end());
}

double MeshPoint::IntersectedVolume(const MeshPoint& other) const{
	Polyhedron poly1(m_poly), poly2(other.m_poly), interPoly;
	try {
		if (CGAL::Polygon_mesh_processing::corefine_and_compute_intersection(poly1, poly2, interPoly)) {
			return CGAL::to_double(CGAL::Polygon_mesh_processing::volume(interPoly));
		}
	}
	catch (std::exception & err) {
		return 0.0;
	}
		
	return 0.0;
}