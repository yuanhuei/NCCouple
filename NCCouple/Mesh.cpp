#include "Mesh.h"
#include <fstream>

void MeshPoint::Init() {
	CGAL::Polygon_mesh_processing::triangulate_faces(m_poly);

	//m_poly = Nef_polyhedron(poly);
	m_volume = CGAL::to_double(CGAL::Polygon_mesh_processing::volume(m_poly));
	m_centerPoint = CGAL::centroid(m_poly.points_begin(), m_poly.points_end());
	for (auto iter = m_poly.vertices_begin(); iter != m_poly.vertices_end(); iter++) {
		m_verticesVec.push_back(std::tuple<double, double, double>{
			CGAL::to_double(iter->point().x()),
				CGAL::to_double(iter->point().y()),
				CGAL::to_double(iter->point().z())
		});
	}
}

double MeshPoint::IntersectedVolume(const MeshPoint& other) const{
	Polyhedron poly1(m_poly), poly2(other.m_poly);
	CGAL::Surface_mesh<Point> mesh1, mesh2, interMesh;
	CGAL::Polygon_mesh_processing::triangulate_faces(poly1);
	CGAL::copy_face_graph(poly1, mesh1);
	CGAL::Polygon_mesh_processing::triangulate_faces(poly2);
	CGAL::copy_face_graph(poly2, mesh2);
	try {
		if (CGAL::Polygon_mesh_processing::corefine_and_compute_intersection(mesh1, mesh2, interMesh)) {
			return CGAL::to_double(CGAL::Polygon_mesh_processing::volume(interMesh));
		}
	}
	catch (std::exception & err) {
		return 0.0;
	}
		
	return 0.0;
}