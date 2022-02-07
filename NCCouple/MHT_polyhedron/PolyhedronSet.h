#pragma once

#ifndef _PolyhedronSet_
#define _PolyhedronSet_

#include "Polyhedron.h"
#define maxDegree 22.6

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

	void Display() const;

	PolyhedronSet ClipByPlane(Vector pointOnPlane, Vector planeNorm);

	void ClipIntoSubPolygons(Scalar);

	void CalculateVolume();

	bool IsContaining(Vector&) const;

	void Move(Vector&);

	PolyhedronSet Copy(Vector&) const;

	Scalar IntersectionVolumeWithPolyhedron(const MHT::Polyhedron&) const;

	Scalar IntersectionVolumeWithPolyhedronSet(const PolyhedronSet&) const;

	std::vector<Scalar> GetRaduisList() const;

	std::pair<bool, Vector> GetAxisCenter() const;

	std::vector<MHT::Polygon> GetFacesOnBoxBoundary(Vector, Vector, Scalar) const;

	void WriteTecplotFile(std::string) const;

	void WriteTecplotFile(ofstream&) const;

	void WriteTecplotHeader(ofstream&) const;

	void WriteTecplotZones(ofstream&) const;
};
#endif