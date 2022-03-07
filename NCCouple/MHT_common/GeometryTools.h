/*---------------------------------------------------------------------------*\
File Name:
	GeometryTools.h

Description:
	Definations and Calculations on basic geometries: 
	line segment, polygon and polyhedron

	Author:		Kong Ling
	Date: 2016-09-30
\*---------------------------------------------------------------------------*/
#pragma once

#ifndef _GeometryTools_
#define _GeometryTools_

#include <iostream>
#include <vector>

#include "Point3D.h"
#include "Vector.h"

//Calculating center point and length of a line segment.
void CalculateLineSegment
(
Point3D<Scalar>& P1,
Point3D<Scalar>& P2,
Point3D<Scalar>& center,
Point3D<Scalar>& length
);

//Calculating center point and area of a triangle.
void CalculateTriangle
(
Point3D<Scalar>& P1,
Point3D<Scalar>& P2,
Point3D<Scalar>& P3,
Point3D<Scalar>& center,
Point3D<Scalar>& area
);

//Calculating center point and area of a polygon.
//it also supports:
//(1) line segment, if 2 vertices are provided;
//(2) triangle, if 3 vertices are provided;
void GetCenterAndArea
(
std::vector<Point3D<Scalar> >& vertice,
Point3D<Scalar>& center,
Point3D<Scalar>& area
);

//Calculating center point and volume of a phyramid.
void GetCenterAndVolume
(
Point3D<Scalar>& apex,
std::vector<Point3D<Scalar> >& base,
Point3D<Scalar>& center,
Scalar& volume
);

//Calculating center point and volume of a polyhedron (from vertices).
void GetCenterAndVolume
(
std::vector<std::vector<Point3D<Scalar> > >& polygonVertice,
Point3D<Scalar>& center,
Scalar& volume
);

//Calculating center point and volume of a polyhedron 
//(from center points and areas of its surrounding faces)
void GetCenterAndVolume
(
std::vector<Point3D<Scalar> >& faceCenter,
std::vector<Point3D<Scalar> >& faceArea,
Point3D<Scalar>& center,
Scalar& volume
);

//rotate a given point around a given axis at a given angle 
//https://en.wikipedia.org/wiki/Rotation_matrix#cite_note-5
Vector RotatePoint
(
	const Vector& GivenPoint,
	const Vector& AxisPoint,
	Vector axis,
	const Scalar& theeta
);

//added by Kong Ling
#include "../MHT_common/Tensor.h"
Tensor RotationTensor
(
	Vector axis,
	Scalar theeta
);

#endif
