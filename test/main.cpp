#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/centroid.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Exact_predicates_exact_constructions_kernel Other_kernel;
typedef Other_kernel::Point_3                             Point;

void fill_cube1(Polyhedron& poly)
{
    std::string input = R"(
OFF
5 5 0
0.63	1.105	0.5
0.63	1.105	0.461538
0.590945	1.10339	0.461538
0.590775	1.10338	0.5
0.609828	1.1361	0.5
4	3	2	1	0
3	0	4	3
3	2	3	4
3	0	1	4
3	1	2	4
)";
    std::stringstream ss;
    ss << input;
    ss >> poly;
}

void fill_cube2(Polyhedron& poly)
{
    std::string input = R"(
OFF
14      6       0
1.26    1.26    0
0.63    1.26    0
0.63    1.105   0
0.722668        1.09587 0
0.811775        1.06884 0
0.893896        1.02495 0
0.965876        0.965876        0
1.26    1.26    0.5
0.63    1.26    0.5
0.63    1.105   0.5
0.722668        1.09587 0.5
0.811775        1.06884 0.5
0.893896        1.02495 0.5
0.965876        0.965876        0.5
7       6       5       4       3       2       1       0
7       7       8       9       10      11      12      13
4       0       1       8       7
4       1       2       9       8
10      2       3       4       5       6       13      12      11      10      9
4       6       0       7       13
)";
    std::stringstream ss;
    ss << input;
    ss >> poly;
}


int main()
{
	Polyhedron poly1, poly2;
	fill_cube1(poly1);
	fill_cube2(poly2);
    CGAL::Surface_mesh<Point> mesh1, mesh2, interMesh;
	CGAL::Polygon_mesh_processing::triangulate_faces(poly1);
    CGAL::copy_face_graph(poly1, mesh1);
	CGAL::Polygon_mesh_processing::triangulate_faces(poly2);
    CGAL::copy_face_graph(poly2, mesh2);
	CGAL::Polygon_mesh_processing::corefine_and_compute_intersection(mesh1, mesh2, interMesh);
    double volume = CGAL::to_double(CGAL::Polygon_mesh_processing::volume(interMesh));
    double volume2 = CGAL::to_double(CGAL::Polygon_mesh_processing::volume(poly1));
    Kernel::Point_3 centerPoint = CGAL::centroid(poly1.points_begin(), poly1.points_end());
    double a = CGAL::to_double(centerPoint.x());

	return 0;
}
