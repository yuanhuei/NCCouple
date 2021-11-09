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
		void Display()
		{
			for (int i = 0; i < (int)v_point.size(); i++)
			{
				cout << v_point[i] << endl;
			}
		}
	};


	class Polyhedron
	{
	public:
		vector <Vector> v_point;
		vector<vector<int> > v_facePointID;
		vector<Vector> v_faceArea;
		vector<Vector> v_faceCenter;
		Scalar volume;
		Vector center;

		Polyhedron();

		Polyhedron(string filename);

		Polyhedron(ifstream& infile);

		void ReadGeometry(ifstream& infile);

		void Display();

		void CalculateFaceGeometry();

		void CalculateVolume();

		void Clear();

		void WriteDataFile(ofstream&);

		void WriteTecplotFile(const string&);

		void CutFaceIntoTriangle(int);

		bool CutWhenNeeded(vector<bool>&);

		bool IsContaining(Vector&);

		bool IsExisting();

		void CreateOldVerticeMap(Vector, Vector, Polyhedron&, vector<bool>&, vector<int>&);

		Polyhedron ClipByPlane(Vector, Vector, vector<int>&);

		pair<Scalar, Scalar> CutByPlane(Vector, Vector);

		Polygon SearchCutPosition(Vector norm, Scalar volume, Scalar tolerance);

		Polygon Reconstruction(Vector norm, Scalar volume, Scalar tolerance);

	};

	Polyhedron operator && (Polyhedron& left, Polyhedron& right);
}

#endif