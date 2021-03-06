#include "Polyhedron.h"
#include "../MHT_common/GeometryTools.h"
#include "../MHT_common/SystemControl.h"
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
		this->CalculateVolume();
		return;
	}

	Polyhedron::Polyhedron
	(
		std::istream& is
	)
	{
		this->ReadGeometry(is);
		this->CalculateVolume();
		return;
	}

	Polyhedron::Polyhedron(Vector minVertice, Vector maxVertice)
	{
		Scalar xmin = minVertice.x_;
		Scalar xmax = maxVertice.x_;
		Scalar ymin = minVertice.y_;
		Scalar ymax = maxVertice.y_;
		Scalar zmin = minVertice.z_;
		Scalar zmax = maxVertice.z_;
		if (xmax <= xmin || ymax <= ymin || zmax <= zmin)
		{
			stringstream ss;
			ss << "invalid vertice: " << minVertice << " and " << maxVertice;
			FatalError(ss.str());
		}
		//create the 8 vertices
		this->v_point.resize(8);
		v_point[0] = Vector(xmin, ymin, zmin);
		v_point[1] = Vector(xmax, ymin, zmin);
		v_point[2] = Vector(xmax, ymax, zmin);
		v_point[3] = Vector(xmin, ymax, zmin);
		v_point[4] = Vector(xmin, ymin, zmax);
		v_point[5] = Vector(xmax, ymin, zmax);
		v_point[6] = Vector(xmax, ymax, zmax);
		v_point[7] = Vector(xmin, ymax, zmax);
		//create the 6 faces
		this->v_facePointID.resize(6);
		v_facePointID[0].push_back(0);
		v_facePointID[0].push_back(3);
		v_facePointID[0].push_back(2);
		v_facePointID[0].push_back(1);
		v_facePointID[1].push_back(0);
		v_facePointID[1].push_back(1);
		v_facePointID[1].push_back(5);
		v_facePointID[1].push_back(4);
		v_facePointID[2].push_back(0);
		v_facePointID[2].push_back(4);
		v_facePointID[2].push_back(7);
		v_facePointID[2].push_back(3);
		v_facePointID[3].push_back(2);
		v_facePointID[3].push_back(3);
		v_facePointID[3].push_back(7);
		v_facePointID[3].push_back(6);
		v_facePointID[4].push_back(1);
		v_facePointID[4].push_back(2);
		v_facePointID[4].push_back(6);
		v_facePointID[4].push_back(5);
		v_facePointID[5].push_back(4);
		v_facePointID[5].push_back(5);
		v_facePointID[5].push_back(6);
		v_facePointID[5].push_back(7);
		//calculate geomety
		this->CalculateVolume();
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

	void Polyhedron::ReadGeometry(std::istream& is)
	{
		std::string oneline;
		std::getline(is, oneline);
		std::getline(is, oneline);
		std::stringstream ss(oneline);
		int nodeNum = 0;
		int faceNum = 0;
		ss >> nodeNum >> faceNum;
		this->v_point.resize(nodeNum);
		this->v_facePointID.resize(faceNum);
		for (int nodeID = 0; nodeID < nodeNum; nodeID++)
		{
			std::getline(is, oneline);
			std::stringstream nodeData(oneline);
			Vector node;
			nodeData >> node.x_ >> node.y_ >> node.z_;
			this->v_point[nodeID] = node;
		}
		for (int faceID = 0; faceID < faceNum; faceID++)
		{
			std::getline(is, oneline);
			std::stringstream faceData(oneline);
			int nodeNumInFace = 0;
			faceData >> nodeNumInFace;
			this->v_facePointID[faceID].resize(nodeNumInFace);
			for (int nodeCount = 0;nodeCount < nodeNumInFace; nodeCount++)
			{
				int nodeID;
				faceData >> nodeID;
				this->v_facePointID[faceID][nodeCount] = nodeID;
			}
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
		std::vector<Vector> verticeList;
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
			this->geometryCalculated = true;
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
			if (geometryCalculated)
			{
				cout << " center" << this->v_faceCenter[i] << ", area" << this->v_faceArea[i];
			}
			cout << endl;
		}
		if (geometryCalculated)
		{
			cout << "polyhedron center: " << this->center << endl;
			cout << "polyhedron volume: " << this->volume << endl;
		}
		return;
	}

	void Polyhedron::Check() const
	{
		//if the polyhedron is not existint, it is meaningless to check
		if (false == this->IsExisting())
		{
			return;
		}
		std::stringstream errorMsg;
		//error 1: there is no face
		if (0 == this->v_faceArea.size())
		{
			FatalError("Polyhedron Check(): no face is found");
			return;
		}
		bool pass = true;
		Vector totalArea(0.0, 0.0, 0.0);
		Scalar denominator = 0.0;
		for (int faceID = 0;faceID < this->v_faceArea.size(); faceID++)
		{
			totalArea += this->v_faceArea[faceID];
			denominator += this->v_faceArea[faceID].Mag();
		}
		//error 2: there are faces, but no face has area
		if (denominator < SMALL)
		{
			errorMsg << "Polyhedron Check(): no face has any area" << std::endl;
			this->Display();
			pass = false;
		}
		else
		{
			Scalar faceAreaError = totalArea.Mag() / denominator;
			//error 3: face areas are not conservative
			if (faceAreaError > SMALL)
			{
				errorMsg << "Polyhedron Check(): the face areas are not conservative in this polyhedron:" << std::endl;
				for (int faceID = 0;faceID < this->v_faceArea.size(); faceID++)
				{
					errorMsg << "face #" << faceID << ": " << this->v_faceArea[faceID] << std::endl;
				}
				errorMsg << "total area:" << totalArea << std::endl;
				pass = false;
			}
		}
		//error 4: wrong face direction
		for (int faceID = 0;faceID < this->v_faceArea.size(); faceID++)
		{
			Vector faceArea = this->v_faceArea[faceID];
			Vector faceCenter = this->v_faceCenter[faceID];
			Vector OP = faceCenter - this->center;
			Scalar OPdotSf = OP & faceArea;
			if (OPdotSf < 0)
			{
				errorMsg << "Polyhedron Check(): direction of face #" << faceID << " is incorrect" << std::endl;
				pass = false;
			}
		}
		//error 5: negative volume
		if (volume < 0)
		{
			errorMsg << "Polyhedron Check(): the volume is negative: " << this->volume << std::endl;
		}
		if (!pass)
		{
			FatalError(errorMsg.str());
		}
		return;
	}

	Scalar Polyhedron::GetVolume() const
	{
		if (false == this->geometryCalculated)
		{
			std::cout << "can not get the volume of this polyhedron, the geomety is not yet calculated" << std::endl;
			exit(1);
			return 0.0;
		}
		else
		{
			return this->volume;
		}
	}

	Vector Polyhedron::GetCenter() const
	{
		if (false == this->geometryCalculated)
		{
			std::cout << "can not get the center of this polyhedron, the geomety is not yet calculated" << std::endl;
			exit(1);
			return Vector(0.0, 0.0, 0.0);
		}
		else
		{
			return this->center;
		}
	}

	Polygon Polyhedron::GetFace(int faceID) const
	{
		if (faceID < 0 || faceID >= this->v_facePointID.size())
		{
			FatalError("in Polyhedron::GetFace(), face ID " + std::to_string(faceID) + " is out of range.");
		}
		Polygon result;
		for (size_t i = 0;i < v_facePointID[faceID].size();i++)
		{
			int nodeID = v_facePointID[faceID][i];
			Vector vertice = v_point[nodeID];
			result.v_point.push_back(vertice);
		}
		result.center = this->v_faceCenter[faceID];
		result.area = this->v_faceArea[faceID];
		return result;
	}

	//collect all faces on boundaries of a box
	std::vector<Polygon> Polyhedron::GetFacesOnBoxBoundary
	(
		Vector verticeMin,//xmin, ymin, zmin
		Vector verticeMax,//xmax, ymax, zmax
		Scalar tolerance
	) const
	{
		std::vector<Polygon> faceList;
		for (size_t i = 0;i < v_faceCenter.size();i++)
		{
			Scalar xCenter = v_faceCenter[i].x_;
			Scalar yCenter = v_faceCenter[i].y_;
			Scalar zCenter = v_faceCenter[i].z_;
			bool isOnBoundary = false;
			if (xCenter - verticeMin.x_ < -tolerance) continue;
			if (xCenter - verticeMax.x_ > tolerance) continue;
			if (yCenter - verticeMin.y_ < -tolerance) continue;
			if (yCenter - verticeMax.y_ > tolerance) continue;
			if (zCenter - verticeMin.z_ < -tolerance) continue;
			if (zCenter - verticeMax.z_ > tolerance) continue;
			if (fabs(xCenter - verticeMin.x_) < tolerance || fabs(xCenter - verticeMax.x_) < tolerance) isOnBoundary = true;
			if (fabs(yCenter - verticeMin.y_) < tolerance || fabs(yCenter - verticeMax.y_) < tolerance) isOnBoundary = true;
			if (fabs(zCenter - verticeMin.z_) < tolerance || fabs(zCenter - verticeMax.z_) < tolerance) isOnBoundary = true;
			if (isOnBoundary)
			{
				//std::cout << "face center = " << v_faceCenter[i] << ", box range = " << verticeMin << " ~ " << verticeMax << std::endl;
				Polygon face = this->GetFace(i);
				//for (int j = 0;j < face.v_point.size();j++)
				//{
				//	std::cout << "point coordinate = " << face.v_point[j] << std::endl;
				//}
				faceList.push_back(face);
				//system("pause");
			}
		}
		return faceList;
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
		std::vector<int> mirrorFace = this->v_facePointID[nface];
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
			std::vector<int> newFace;
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

	void Polyhedron::CreateOldVerticeMap
	(
		Vector pointOnPlane, Vector norm, 
		Polyhedron& target, 
		std::vector<bool>& inside,
		std::vector<int>& map
	)
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

	bool Polyhedron::CutWhenNeeded(std::vector<bool>& inside)
	{
		bool cutDone = false;
		std::vector<int> faceToBeCut;
		for (int faceID = 0; faceID < (int)this->v_facePointID.size(); faceID++)
		{
			std::vector<int>& thisFace = this->v_facePointID[faceID];
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

	void Polyhedron::Move(Vector& translation)
	{
		//move vertices
		for (size_t i = 0;i < this->v_point.size();i++)
		{
			this->v_point[i] += translation;
		}
		if (geometryCalculated)
		{
			//move face centers
			for (size_t i = 0;i < this->v_faceCenter.size();i++)
			{
				this->v_faceCenter[i] += translation;
			}
			//move polyhedron center
			this->center += translation;
		}
		return;
	}

	Polyhedron Polyhedron::Copy(Vector& translation) const
	{
		Polyhedron result = *(this);
		result.Move(translation);
		return result;
	}

	Polyhedron Polyhedron::ClipByPlane
	(
		Vector pointOnPlane, 
		Vector planeNorm, 
		std::vector<int>& newFacelist
	)
	{
		Vector clipNorm = planeNorm.GetNormal();
		newFacelist.clear();
		//create temperal data needed in the algorithm
		std::vector<bool> v_inside;
		std::vector<int> v_newPointID;
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

		std::vector<int> v_oldFaceID;
		std::vector<bool> v_faceCut;

		//visit all faces to see whether they need to be created in the new polyhedron
		for (int i = 0; i < (int)this->v_facePointID.size(); i++)
		{
			std::vector<int>& thisFace = this->v_facePointID[i];
			for (int j = 0; j < (int)this->v_facePointID[i].size(); j++)
			{
				int pointID = thisFace[j];
				if (true == v_inside[pointID])
				{
					std::vector<int> newface;
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
			std::vector<int>& refFace = this->v_facePointID[refFaceID];
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
		std::vector<std::pair<int, int> > v_sideArror;
		for (int i = 0; i < (int)targetPoly.v_facePointID.size(); i++)
		{
			if (false == v_faceCut[i])
			{
				continue;
			}
			std::vector<int>& thisFace = targetPoly.v_facePointID[i];
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

		std::vector<std::vector<pair<int, int> > > v_groupedArror;

		v_groupedArror.resize(1);
		v_groupedArror[0].resize(1);
		v_groupedArror[0][0] = v_sideArror[0];
		for (int i = 1; i < (int)v_sideArror.size(); i++)
		{
			if (v_sideArror[i].first != v_sideArror[i - 1].second)
			{
				std::vector<pair<int, int> > newgroup;
				v_groupedArror.push_back(newgroup);
			}
			v_groupedArror[v_groupedArror.size() - 1].push_back(v_sideArror[i]);
		}

		for (int i = 0; i < (int)v_groupedArror.size(); i++)
		{
			std::vector<int> newface;
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
		std::vector<int> cutfaceID;
		for (int i = 0; i < (int)right.v_facePointID.size(); i++)
		{
			result = result.ClipByPlane(right.v_faceCenter[i], right.v_faceArea[i], cutfaceID);
		}
		result.CalculateVolume();
		return result;
	}

	Polygon Polyhedron::SearchCutPosition
	(
		Vector normal,
		Scalar volumeFraction,
		Scalar tolerance
	)
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
		std::vector<int> v_cutFaceID;
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
				std::vector<Vector> verticeList;
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

	Polygon Polyhedron::Reconstruction
	(
		Vector normal,
		Scalar volumeFraction,
		Scalar tolerance
	)
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


	void Polyhedron::WriteTecplotZones(std::ofstream& outFile) const
	{
		int nCount(0);
		int faceNum = (int)this->v_facePointID.size();
		for (int i = 0; i < faceNum; i++)
		{
			nCount += (int)this->v_facePointID[i].size();
		}

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
}
