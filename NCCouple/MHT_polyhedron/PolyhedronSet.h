#pragma once

#ifndef _PolyhedronSet_
#define _PolyhedronSet_

#include "Polyhedron.h"
#define maxDegree 20.01

class PolyhedronSet :public MHT::Polyhedron
{
public:
	bool clipped = false;
	std::vector<std::pair<int, Scalar> > v_curvedFace;
	Vector axisCenter;
	Vector axisNorm;
	std::vector<MHT::Polyhedron> v_subPolyhedron;
public:
	PolyhedronSet();

	PolyhedronSet(std::istream&,std::vector<int>&, Vector, Vector);

	PolyhedronSet(string);

	void ReadCurveFaces(ifstream& infile);

	void CalculateRadius();

	void Display();

	PolyhedronSet ClipByPlane(Vector pointOnPlane, Vector planeNorm);

	void ClipIntoSubPolygons(Scalar);

	void CalculateVolume();

	bool IsContaining(Vector&);

	Scalar IntersectionVolumeWithPolyhedron(MHT::Polyhedron&);

	Scalar IntersectionVolumeWithPolyhedronSet(PolyhedronSet&);

	void WriteTecplotFile(std::string);

};
#endif