#include "Polyhedron.h"
#include "../MHT_common/GeometryTools.h"
#include <algorithm>

namespace MHT
{
	Polyhedron::Polyhedron()
		:volume(0.0),geometryCalculated(false)
	{
		v_point.clear();
		v_facePointID.clear();
		v_faceArea.clear();
		v_faceCenter.clear();
	}

	Polyhedron::Polyhedron
	(
		string filename
	)
	{
		int numPoints;
		ifstream infile(filename);
		this->ReadGeometry(infile);
		infile.close();
		return;
	}

	Polyhedron::Polyhedron
	(
		ifstream& infile
	)
		:volume(0.0), geometryCalculated(false)
	{
		this->ReadGeometry(infile);
		return;
	}

	void Polyhedron::ReadGeometry(ifstream& infile)
	{
		int numPoints;
		infile >> numPoints;
		for (int i = 0; i < numPoints; i++)
		{
			Scalar x, y, z;
			infile >> x >> y >> z;
			v_point.push_back(Vector(x, y, z));
		}
		int faceNum;
		infile >> faceNum;
		this->v_facePointID.resize(faceNum);
		this->v_faceArea.resize(faceNum);
		this->v_faceCenter.resize(faceNum);
		for (int i = 0; i < faceNum; i++)
		{
			int numSides;
			infile >> numSides;
			this->v_facePointID[i].resize(numSides);
			for (int j = 0; j < numSides; j++)
			{
				int pointID;
				infile >> pointID;
				this->v_facePointID[i][j] = pointID;
			}
			this->v_faceArea[i] = Vector(0.0, 0.0, 0.0);
			this->v_faceCenter[i] = Vector(0.0, 0.0, 0.0);
		}
		return;
	}

	void Polyhedron::Clear()
	{
		this->v_point.clear();
		this->v_facePointID.clear();
		this->v_faceCenter.clear();
		this->v_faceArea.clear();
		this->volume = 0.0;
		this->geometryCalculated = false;
		return;
	}

	void Polyhedron::CalculateFaceGeometry()
	{
		int faceNum = (int)this->v_facePointID.size();
		this->v_faceCenter.resize(faceNum);
		this->v_faceArea.resize(faceNum);
		//calculate face centers and norms
		Vector centerTemp(0.0, 0.0, 0.0);
		Vector areaTemp(0.0, 0.0, 0.0);
		vector<Vector> verticeList;
		for (int faceID = 0; faceID < faceNum; faceID++)
		{
			verticeList.clear();
			for (int i = 0; i < (int)this->v_facePointID[faceID].size(); i++)
			{
				int pointID = v_facePointID[faceID][i];
				verticeList.push_back(this->v_point[pointID]);
			}
			GetCenterAndArea(verticeList, centerTemp, areaTemp);
			this->v_faceCenter[faceID] = centerTemp;
			this->v_faceArea[faceID] = areaTemp;
		}
		return;
	}

	void Polyhedron::CalculateVolume()
	{
		//calculate the volume of this polyhedron
		if (0 == this->v_point.size() || 0 == this->v_facePointID.size())
		{
			this->volume = 0.0;
			return;
		}
		this->CalculateFaceGeometry();
		GetCenterAndVolume(this->v_faceCenter, this->v_faceArea, this->center, this->volume);
		this->geometryCalculated = true;
		return;
	}

	void Polyhedron::Display() const
	{
		cout << "This polyhedron is composed of points:" << endl;
		for (int i = 0; i < (int)v_point.size(); i++)
		{
			cout << i << ": " << v_point[i] << endl;
		}
		cout << "and faces:" << endl;
		for (int i = 0; i < (int)v_facePointID.size(); i++)
		{
			cout << i << ": ";
			for (int j = 0; j < (int)v_facePointID[i].size(); j++)
			{
				cout << v_facePointID[i][j] << " ";
			}
			cout << endl;
		}
		return;
	}

	void Polyhedron::WriteDataFile(ofstream& outfile) const
	{
		outfile << this->v_point.size() << std::endl;
		for (int i = 0; i < this->v_point.size(); i++)
		{
			Scalar x = this->v_point[i].x_;
			Scalar y = this->v_point[i].y_;
			Scalar z = this->v_point[i].z_;
			outfile << x << "\t" << y << "\t" << z << std::endl;
		}
		outfile << this->v_facePointID.size() << std::endl;
		for (int i = 0; i < this->v_facePointID.size(); i++)
		{
			outfile << this->v_facePointID[i].size() << "\t";
			for (int j = 0; j < this->v_facePointID[i].size(); j++)
			{
				outfile << this->v_facePointID[i][j] << "\t";
			}
			outfile << std::endl;
		}
		return;
	}

	void Polyhedron::WriteTecplotFile(const string& filename) const
	{
		int nCount(0);
		int faceNum = (int)this->v_facePointID.size();
		for (int i = 0; i < faceNum; i++)
		{
			nCount += (int)this->v_facePointID[i].size();
		}

		ofstream	outFile(filename.c_str());

		outFile << "TITLE =\"" << "polyhedron" << "\"" << endl;
		outFile << "VARIABLES = " << "\"x\"," << "\"y\"," << "\"z\"" << endl;
		outFile << "ZONE T = " << "\"" << "polyhedron" << "\"," << " DATAPACKING = POINT, N = " << v_point.size() << ", E = "
			<< nCount << ", ZONETYPE = FELINESEG" << std::endl;

		for (unsigned int i = 0; i < v_point.size(); i++)
		{
			outFile << setprecision(8) << setiosflags(ios::scientific) << v_point[i].x_ << " " << v_point[i].y_ << " " << v_point[i].z_ << endl;
		}

		for (int i = 0; i < faceNum; i++)
		{
			for (int n = 0; n < (int)this->v_facePointID[i].size(); n++)
			{
				if (n == (int)this->v_facePointID[i].size() - 1)
				{
					outFile << v_facePointID[i][n] + 1 << "\t" << v_facePointID[i][0] + 1 << std::endl;
				}
				else
				{
					outFile << v_facePointID[i][n] + 1 << "\t" << v_facePointID[i][n + 1] + 1 << std::endl;
				}
			}
		}
		return;
	}

	void Polyhedron::CutFaceIntoTriangle(int nface)
	{
		vector<int> mirrorFace = this->v_facePointID[nface];
		int numVertices = (int)mirrorFace.size();
		if (3 == numVertices)
		{
			return;
		}
		//calculate the center point of the original face
		Vector center(0.0, 0.0, 0.0);
		for (int i = 0; i < numVertices; i++)
		{
			center += this->v_point[mirrorFace[i]];
		}
		center = center / (Scalar)numVertices;
		int newPointID = (int)this->v_point.size();
		this->v_point.push_back(center);
		//new triangular faces created
		for (int i = 1; i < numVertices; i++)
		{
			vector<int> newFace;
			newFace.push_back(mirrorFace[i]);
			newFace.push_back(mirrorFace[(i + 1) % numVertices]);
			newFace.push_back(newPointID);
			this->v_facePointID.push_back(newFace);
		}
		//modify the original face
		this->v_facePointID[nface][2] = newPointID;
		for (int i = 0; i < numVertices - 3; i++)
		{
			this->v_facePointID[nface].pop_back();
		}
		return;
	}

	void Polyhedron::CreateOldVerticeMap(Vector pointOnPlane, Vector norm, Polyhedron& target, vector<bool>& inside, vector<int>& map)
	{
		target.Clear();
		inside.resize(this->v_point.size());
		map.resize(this->v_point.size());
		//for all points, decide which side they are located in
		for (int i = 0; i < (int)this->v_point.size(); i++)
		{
			if (((pointOnPlane - v_point[i]) & norm) > 0)
			{
				inside[i] = true;
				map[i] = (int)target.v_point.size();
				target.v_point.push_back(this->v_point[i]);
			}
			else
			{
				map[i] = -1;
				inside[i] = false;
			}
		}
		return;
	}

	bool Polyhedron::CutWhenNeeded(vector<bool>& inside)
	{
		bool cutDone = false;
		vector<int> faceToBeCut;
		for (int faceID = 0; faceID < (int)this->v_facePointID.size(); faceID++)
		{
			vector<int>& thisFace = this->v_facePointID[faceID];
			int numVertices = (int)thisFace.size();
			int numIntersections = 0;
			for (int j = 0; j < (int)thisFace.size(); j++)
			{
				int startID = thisFace[j % numVertices];
				int endID = thisFace[(j + 1) % numVertices];
				if (inside[startID] != inside[endID])
				{
					numIntersections++;
				}
			}
			if (numIntersections > 2)
			{
				faceToBeCut.push_back(faceID);
				cutDone = true;
			}
		}
		for (int i = 0; i < (int)faceToBeCut.size(); i++)
		{
			//cout << "cut done" << endl;
			this->CutFaceIntoTriangle(faceToBeCut[i]);
		}
		return cutDone;
	}

	bool Polyhedron::IsContaining(Vector& point) const
	{
		bool contain = true;
		for (int i = 0; i < this->v_faceCenter.size(); i++)
		{
			Vector faceCenter = this->v_faceCenter[i];
			Vector pointToFace = faceCenter - point;
			if ((pointToFace & this->v_faceArea[i]) < 0)
			{
				contain = false;
				break;
			}
		}
		return contain;
	}

	bool Polyhedron::IsExisting() const
	{
		if (this->v_point.size() > 0)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	Polyhedron Polyhedron::ClipByPlane
	(
		Vector pointOnPlane, 
		Vector planeNorm, 
		vector<int>& newFacelist
	)
	{
		Vector clipNorm = planeNorm.GetNormal();
		newFacelist.clear();
		//create temperal data needed in the algorithm
		vector<bool> v_inside;
		vector<int> v_newPointID;
		//create a new polygon without any data
		Polyhedron targetPoly;

		CreateOldVerticeMap(pointOnPlane, clipNorm, targetPoly, v_inside, v_newPointID);

		//before creating the new polygon, check intersections in each face, and cut into triangles when needed
		if (true == CutWhenNeeded(v_inside))
		{
			//if cut occurs, operations on vertices need to be re-proceeded
			CreateOldVerticeMap(pointOnPlane, clipNorm, targetPoly, v_inside, v_newPointID);
		}

		int numOldPoints = (int)targetPoly.v_point.size();

		vector<int> v_oldFaceID;
		vector<bool> v_faceCut;

		//visit all faces to see whether they need to be created in the new polyhedron
		for (int i = 0; i < (int)this->v_facePointID.size(); i++)
		{
			vector<int>& thisFace = this->v_facePointID[i];
			for (int j = 0; j < (int)this->v_facePointID[i].size(); j++)
			{
				int pointID = thisFace[j];
				if (true == v_inside[pointID])
				{
					vector<int> newface;
					targetPoly.v_facePointID.push_back(newface);
					v_oldFaceID.push_back(i);
					v_faceCut.push_back(false);
					break;
				}
			}
		}

		//Note: the three integers in v_intersection refer to
		//(1) the smaller and (2) the greater point ID of the new polyhedron and
		//(3) the point ID in the new polyhedron
		map<pair<int, int>, int> v_intersection;
		//visit all faces and insert points and faces into the new polyhedron
		for (int i = 0; i < (int)v_oldFaceID.size(); i++)
		{
			int refFaceID = v_oldFaceID[i];
			vector<int>& refFace = this->v_facePointID[refFaceID];
			int numSides = (int)refFace.size();
			for (int j = 0; j < numSides; j++)
			{
				int thisID = refFace[j];
				int nextID = refFace[(j + 1) % numSides];
				if (v_inside[thisID])
				{
					targetPoly.v_facePointID[i].push_back(v_newPointID[thisID]);
				}
				if (v_inside[thisID] != v_inside[nextID])
				{
					v_faceCut[i] = true;
					//search
					pair<int, int> key = pair<int, int>(Min(thisID, nextID), Max(thisID, nextID));
					map<pair<int, int>, int>::iterator it = v_intersection.find(key);
					if (v_intersection.end() == it)
					{
						//if not found, insert
						Vector A = this->v_point[thisID];
						Vector B = this->v_point[nextID];
						Vector AToB = B - A;
						Scalar denominator = AToB & clipNorm;
						Scalar lambda = 0.5;
						if (fabs(denominator) > SMALL * 10)
						{
							lambda = ((pointOnPlane - A) & clipNorm) / denominator;
						}
						Vector intersection = A + lambda * AToB;
						int IDInNew = (int)targetPoly.v_point.size();
						targetPoly.v_point.push_back(intersection);
						targetPoly.v_facePointID[i].push_back(IDInNew);
						v_intersection.insert(pair<pair<int, int>, int>(key, IDInNew));
					}
					else
					{
						targetPoly.v_facePointID[i].push_back(it->second);
					}
				}
			}
		}

		//create the last polygon face
		vector<pair<int, int> > v_sideArror;
		for (int i = 0; i < (int)targetPoly.v_facePointID.size(); i++)
		{
			if (false == v_faceCut[i])
			{
				continue;
			}
			vector<int>& thisFace = targetPoly.v_facePointID[i];
			int numVertices = (int)thisFace.size();
			//start from an old point
			int jstart = 0;
			for (int j = 0; j < numVertices; j++)
			{
				if (thisFace[j] < numOldPoints)
				{
					jstart = j;
					break;
				}
			}
			for (int j = jstart; j < jstart + numVertices; j++)
			{
				int startID = thisFace[j % numVertices];
				if (startID >= numOldPoints)
				{
					int endID = thisFace[(j + 1) % numVertices];
					v_sideArror.push_back(pair<int, int>(startID, endID));
					break;
				}
			}
		}

		//insert the new face coming from a part of the clip plane
		if (0 == v_sideArror.size())
		{
			return targetPoly;
		}

		for (int i = 0; i < (int)v_sideArror.size() - 1; i++)
		{
			int desiredNext = i + 1;
			for (int j = i + 1; j < (int)v_sideArror.size(); j++)
			{
				if (v_sideArror[i].second == v_sideArror[j].first)
				{
					desiredNext = j;
					break;
				}
			}
			if (desiredNext != i + 1)
			{
				//swap
				pair<int, int> tempSideArror = v_sideArror[i + 1];
				v_sideArror[i + 1] = v_sideArror[desiredNext];
				v_sideArror[desiredNext] = tempSideArror;
			}
		}

		vector<vector<pair<int, int> > > v_groupedArror;

		v_groupedArror.resize(1);
		v_groupedArror[0].resize(1);
		v_groupedArror[0][0] = v_sideArror[0];
		for (int i = 1; i < (int)v_sideArror.size(); i++)
		{
			if (v_sideArror[i].first != v_sideArror[i - 1].second)
			{
				vector<pair<int, int> > newgroup;
				v_groupedArror.push_back(newgroup);
			}
			v_groupedArror[v_groupedArror.size() - 1].push_back(v_sideArror[i]);
		}

		for (int i = 0; i < (int)v_groupedArror.size(); i++)
		{
			vector<int> newface;
			for (int j = (int)v_groupedArror[i].size() - 1; j >= 0; j--)
			{
				newface.push_back(v_groupedArror[i][j].first);
			}
			targetPoly.v_facePointID.push_back(newface);
			newFacelist.push_back((int)targetPoly.v_facePointID.size() - 1);
		}

		return targetPoly;

	}

	Polyhedron operator && (const Polyhedron& left, const Polyhedron& right)
	{
		Polyhedron result = left;
		if (false == right.geometryCalculated)
		{
			std::cout << "Fatal error: in Polyhedron operator &&, the geometry of the rhs is not calculated" << std::endl;
			exit(1);
			return result;
		}
		vector<int> cutfaceID;
		for (int i = 0; i < (int)right.v_facePointID.size(); i++)
		{
			result = result.ClipByPlane(right.v_faceCenter[i], right.v_faceArea[i], cutfaceID);
		}
		return result;
	}

	Polygon Polyhedron::SearchCutPosition(Vector normal, Scalar volumeFraction, Scalar tolerance)
	{
		Vector norm = normal.GetNormal();
		Scalar volumeDesired = volumeFraction * this->volume;
		Vector refPoint = this->v_point[0];
		int imin(0), imax(0);
		Scalar tempmin(0.0), tempmax(0.0);
		for (int i = 0; i < (int)this->v_point.size(); i++)
		{
			Vector refVector = this->v_point[i] - refPoint;
			Scalar candidate = refVector & norm;
			if (candidate > tempmax)
			{
				tempmax = candidate;
				imax = i;
			}
			if (candidate < tempmin)
			{
				tempmin = candidate;
				imin = i;
			}
		}
		Scalar totalLength = (this->v_point[imax] - this->v_point[imin]) & norm;
		refPoint = this->v_point[imin];
		Scalar alphaLeft = 0.0;
		Scalar alphaRight = totalLength;
		Scalar alphaMiddle = volumeFraction * totalLength;
		Polyhedron cutResult;
		vector<int> v_cutFaceID;
		for (int i = 0; i < 100; i++)
		{
			cout << "i = " << i << endl;
			Vector pointToTry = refPoint + alphaMiddle * norm;
			cutResult = this->ClipByPlane(pointToTry, norm, v_cutFaceID);
			cutResult.CalculateVolume();
			Scalar cutVolume = cutResult.volume;
			Scalar cutArea = 0.0;
			for (int i = 0; i < (int)v_cutFaceID.size(); i++)
			{
				int faceID = (int)v_cutFaceID[i];
				vector<Vector> verticeList;
				for (int i = 0; i < (int)cutResult.v_facePointID[faceID].size(); i++)
				{
					int pointID = cutResult.v_facePointID[faceID][i];
					verticeList.push_back(cutResult.v_point[pointID]);
				}
				Vector center;
				Vector area;
				GetCenterAndArea(verticeList, center, area);
				cutArea += area.Mag();
			}

			if (fabs(cutVolume / this->volume - volumeFraction) < tolerance)
			{
				break;
			}
			if (cutArea < SMALL)
			{
				if (cutVolume > volumeDesired)
				{
					alphaRight = alphaMiddle;
				}
				else
				{
					alphaLeft = alphaMiddle;
				}
				alphaMiddle = 0.5 * (alphaLeft + alphaRight);
				continue;
			}
			Scalar alphaNext = alphaMiddle + (volumeDesired - cutVolume) / cutArea;
			if (alphaNext > alphaRight)
			{
				alphaLeft = alphaMiddle;
				alphaMiddle = 0.5 * (alphaLeft + alphaRight);
				continue;
			}
			if (alphaNext < alphaLeft)
			{
				alphaRight = alphaMiddle;
				alphaMiddle = 0.5 * (alphaLeft + alphaRight);
				continue;
			}
			alphaMiddle = alphaNext;
		}

		cutResult.WriteTecplotFile("finalPoly.plt");
		Polygon cutFace;
		if (0 == v_cutFaceID.size())
		{
			return cutFace;
		}
		int faceID = v_cutFaceID[0];
		for (int i = 0; i < (int)cutResult.v_facePointID[faceID].size(); i++)
		{
			int pointID = cutResult.v_facePointID[faceID][i];
			cutFace.v_point.push_back(cutResult.v_point[pointID]);
		}
		return cutFace;
	}

	Polygon Polyhedron::Reconstruction(Vector normal, Scalar volumeFraction, Scalar tolerance)
	{
		Polyhedron standardPoly = (*this);
		Scalar alpha = pow(this->volume, 1.0 / 3.0);
		for (int i = 0; i < (int)standardPoly.v_point.size(); i++)
		{
			standardPoly.v_point[i] = (standardPoly.v_point[i] - this->center) / alpha;
		}
		standardPoly.center = Vector(0.0, 0.0, 0.0);
		standardPoly.volume = 1.0;
		Polygon resultPolygon = standardPoly.SearchCutPosition(normal, volumeFraction, tolerance);
		for (int i = 0; i < (int)resultPolygon.v_point.size(); i++)
		{
			resultPolygon.v_point[i] = alpha * resultPolygon.v_point[i] + this->center;
		}
		return resultPolygon;
	}
}
