/*---------------------------------------------------------------------------*\
Class: Polyhedron
Description: definition of polyhedron
Author:		Shan Fan, Yifan Wang
Date: 2018-05-26
\*---------------------------------------------------------------------------*/

#pragma once

#ifndef _Polyhedron_
#define _Polyhedron_

#include "../MHT_common/Vector.h"
#include <map>
#include <iomanip>

namespace MHT
{
	class Polygon
	{
	public:
		std::vector <Vector> v_point;
		Vector area;
		Vector center;
	public:
		Polygon() {}
	};


	class Polyhedron
	{
	public:
		std::vector<Vector> v_point;
		std::vector<std::vector<int> > v_facePointID;
		std::vector<Vector> v_faceArea;
		std::vector<Vector> v_faceCenter;
		bool geometryCalculated;
		Scalar volume;
	private:
		Vector center;
	public:

		Polyhedron();

		Polyhedron(string filename);

		Polyhedron(ifstream& infile);

		void Display() const;

		Scalar GetVolume() const;

		Vector GetCenter() const;

		void CalculateVolume();

		void Clear();

		void WriteDataFile(ofstream&) const;

		void WriteTecplotFile(const string&) const;

		bool CutWhenNeeded(std::vector<bool>&);

		void CreateOldVerticeMap(Vector, Vector, Polyhedron&, std::vector<bool>&, std::vector<int>&);

		bool IsContaining(Vector&) const;

		bool IsExisting() const;

		Polyhedron ClipByPlane(Vector, Vector, std::vector<int>&);

		Polygon SearchCutPosition(Vector norm, Scalar volume, Scalar tolerance);

		Polygon Reconstruction(Vector norm, Scalar volume, Scalar tolerance);

	private:

		void ReadGeometry(ifstream& infile);

		void CalculateFaceGeometry();

		void CutFaceIntoTriangle(int);

	};

	Polyhedron operator && (const Polyhedron& left, const Polyhedron& right);
}

#endif