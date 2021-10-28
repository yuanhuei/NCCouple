
#include "GeometryTools.h"

//Calculating center point and length of a line segment.
void CalculateLineSegment
(
Point3D<Scalar>& P1,
Point3D<Scalar>& P2,
Point3D<Scalar>& center,
Point3D<Scalar>& length
)
{
	center = (P1 + P2) / 2.0;
	length = (P2 - P1) ^ Point3D<Scalar>(0, 0, 1);
}

//Calculating center point and area of a triangle.
void CalculateTriangle
(
Point3D<Scalar>& P1,
Point3D<Scalar>& P2,
Point3D<Scalar>& P3,
Point3D<Scalar>& center,
Point3D<Scalar>& area
)
{
	center = (P1 + P2 + P3) / 3.0;
	area = ((P2 - P1) ^ (P3 - P1)) / 2;
}

//Calculating center point and area of a polygon.
//it also supports:
//(1) line segment, if 2 vertices are provided;
//(2) triangle, if 3 vertices are provided;
void GetCenterAndArea
(
vector<Point3D<Scalar> >& vertice,
Point3D<Scalar>& center,
Point3D<Scalar>& area
)
{
	int n = (int)vertice.size();

	if (n < 2)
	{
		std::cout << "in GetCenterAndArea(), at least 2 vertices are required.";
		std::cout << "You have only " << n << " vertices" << endl;
		system("pause");
		exit(1);
	}

	if (2 == n)
	{
		CalculateLineSegment(vertice[0], vertice[1], center, area);
	}
	else if (3 == n)
	{
		CalculateTriangle(vertice[0], vertice[1], vertice[2], center, area);
	}
	else
	{
		Point3D<Scalar> xG(0, 0, 0);
		for (int i = 0; i < n; i++)
		{
			xG = xG + vertice[i];
		}
		xG = xG / Scalar(n);

		Point3D<Scalar> numerator(0, 0, 0);
		Scalar dominator = 0;
		area = Point3D<Scalar>(0, 0, 0);
		for (int i = 0; i < n; i++)
		{
			Point3D<Scalar> triangleCenter, triangleArea;
			CalculateTriangle(vertice[i % n], vertice[(i + 1) % n], xG, triangleCenter, triangleArea);
			area += triangleArea;
			numerator += triangleArea.Mag() * triangleCenter;
			dominator += triangleArea.Mag();
		}
		center = (dominator > 10 * SMALL) ? (numerator / dominator) : xG;
	}
}

//Calculating center point and volume of a phyramid.
void GetCenterAndVolume
(
Point3D<Scalar>& apex,
vector<Point3D<Scalar> >& base,
Point3D<Scalar>& center,
Scalar& volume
)
{
	Point3D<Scalar> baseCenter, baseArea;
	GetCenterAndArea(base, baseCenter, baseArea);
	center = 0.75*baseCenter + 0.25*apex;
	Point3D<Scalar> norm = baseArea;
	norm.Normalize();
	Scalar height = abs(norm & (apex - baseCenter));
	volume = height*baseArea.Mag() / 3.0;
}

//Calculating center point and volume of a polyhedron (from vertices).
void GetCenterAndVolume
(
vector<vector<Point3D<Scalar> > >& polygonVertice,
Point3D<Scalar>& center,
Scalar& volume
)
{
	int n = (int)polygonVertice.size();

	Point3D<Scalar> xG(0, 0, 0);
	for (int i = 0; i < n; i++)
	{
		Point3D<Scalar> faceCenter, faceArea;
		GetCenterAndArea(polygonVertice[i], faceCenter, faceArea);
		xG = xG + faceCenter;
	}
	xG = xG / Scalar(n);

	Point3D<Scalar> numerator(0, 0, 0);
	Scalar dominator = 0;
	volume = 0;
	for (int i = 0; i < n; i++)
	{
		Point3D<Scalar> pyramidCenter;
		Scalar pyramidVolume;
		GetCenterAndVolume(xG, polygonVertice[i], pyramidCenter, pyramidVolume);
		numerator += pyramidVolume*pyramidCenter;
		dominator += pyramidVolume;
	}

	center = (dominator > 10 * SMALL) ? (numerator / dominator) : xG;
	volume = dominator;
}

//Calculating center point and volume of a polyhedron 
//(from center points and areas of its surrounding faces)
void GetCenterAndVolume
(
vector<Point3D<Scalar> >& faceCenter,
vector<Point3D<Scalar> >& faceArea,
Point3D<Scalar>& center,
Scalar& volume
)
{
	int n = (int)faceCenter.size();
	if (n != faceArea.size())
	{
		std::cout << "Error found in GetCenterAndVolume(): face numbers are not consistent" << std::endl;
		exit(1);
	}

	Point3D<Scalar> xG(0, 0, 0);
	for (int i = 0; i < n; i++)
	{
		xG = xG + faceCenter[i];
	}
	xG = xG / Scalar(n);

	Point3D<Scalar> numerator(0, 0, 0);
	Scalar dominator = 0;

	for (int i = 0; i < n; i++)
	{
		Point3D<Scalar> pyramidCenter = 0.75*faceCenter[i] + 0.25*xG;
		Point3D<Scalar> norm = faceArea[i];
		if (norm.Mag()<SMALL) continue;
		norm.Normalize();
		Scalar height = abs(norm & (xG - faceCenter[i]));
		Scalar pyramidVolume = height*faceArea[i].Mag() / 3.0;
		numerator += pyramidVolume*pyramidCenter;
		dominator += pyramidVolume;

	}
	center = (dominator > 10 * SMALL) ? (numerator / dominator) : xG;
	volume = dominator;
}