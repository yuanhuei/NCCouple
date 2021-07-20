#include <CGAL/Polyhedron_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/boost/graph/convert_nef_polyhedron_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>

typedef CGAL::Cartesian<double> Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

void fill_cube(Polyhedron& poly)
{
	std::string input = R"(
OFF
5 5 0
1.10339 0.590945    0.0384615
1.10338 0.590775    0
1.105   0.63    0
1.105   0.63    0.0384615
1.1361  0.609828    0
4   3   2   1   0
3   2   4   1
3   2   3   4
3   3   0   4
3   0   1   4
)";
	std::stringstream ss;
	ss << input;
	ss >> poly;
}

int main() {
	Polyhedron poly1, poly2, interPoly;
	fill_cube(poly1), fill_cube(poly2);
	if (CGAL::Polygon_mesh_processing::corefine_and_compute_intersection(poly1, poly2, interPoly)) {
		double volume = CGAL::Polygon_mesh_processing::volume(interPoly);
	}

	return 0;
}