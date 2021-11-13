#include "Mesh.h"
#include <fstream>

void CGALMeshPoint::Init() {
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

double CGALMeshPoint::IntersectedVolume(const MeshPoint& other) const{
	const CGALMeshPoint* p_other = dynamic_cast<const CGALMeshPoint*>(&other);
	if (!p_other)
		throw std::runtime_error("the type of input MeshPoint is not same with current MeshPoint!");

	Polyhedron poly1(m_poly), poly2(p_other->m_poly);
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

double LingMeshPoint::IntersectedVolume(const MeshPoint& other) const {
	const LingMeshPoint* p_other = dynamic_cast<const LingMeshPoint*>(&other);
	if (!p_other)
		throw std::runtime_error("the type of input MeshPoint is not same with current MeshPoint!");

	return m_poly.IntersectionVolumeWithPolyhedronSet(p_other->m_poly);
}