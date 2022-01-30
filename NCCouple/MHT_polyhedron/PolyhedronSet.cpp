#include "PolyhedronSet.h"
#include "../MHT_common/Tensor.h"

PolyhedronSet::PolyhedronSet()
	:MHT::Polyhedron()
{}

PolyhedronSet::PolyhedronSet(std::string filename)
	: MHT::Polyhedron()
{
	ifstream infile(filename);
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
	this->ReadCurveFaces(infile);
	this->CalculateRadius();
	infile.close();
	this->ClipIntoSubPolygons(maxDegree);
	this->CalculateVolume();
	return;
}

PolyhedronSet::PolyhedronSet
(
	//OFF stream describing the entire polyhedron
	std::istream& OFFstream,
	//curved marked saved in a list, 0 for planar face and 1 for curcular face
	std::vector<int>& v_curvedFaceMark,
	//one point on the tube axis
	Vector axisPoint,
	//normal direction of the tube axis
	Vector axisNorm
)
	: MHT::Polyhedron(OFFstream)
{
	int faceNum = this->v_facePointID.size();
	this->v_curvedFace.resize(faceNum);
	for (int i = 0; i < faceNum; i++)
	{
		this->v_curvedFace[i].first = v_curvedFaceMark[i];
	}
	this->axisCenter = axisPoint;
	if (axisNorm.Mag() > 10.0 * SMALL)
	{
		this->axisNorm = axisNorm.GetNormal();
	}
	else
	{
		this->axisNorm = Vector(0.0, 0.0, 1.0);
	}
	this->CalculateRadius();
	//see whether or not the polyhedron has curved faces
	bool havingCurvedFace = false;
	for (int i = 0;i < this->v_curvedFace.size();i++)
	{
		if (1 == this->v_curvedFace[i].first)
		{
			havingCurvedFace = true;
			break;
		}
	}
	//if having curved faces, the polyhedron should be clipped into sub-polyhedrons
	if (havingCurvedFace)
	{
		this->ClipIntoSubPolygons(maxDegree);
	}
}

void PolyhedronSet::ReadCurveFaces(ifstream& infile)
{
	int faceNum = this->v_facePointID.size();
	this->v_curvedFace.resize(faceNum);
	for (int i = 0; i < faceNum; i++)
	{
		infile >> v_curvedFace[i].first;
	}
	infile >> axisCenter.x_ >> axisCenter.y_ >> axisCenter.z_;
	infile >> axisNorm.x_ >> axisNorm.y_ >> axisNorm.z_;
	axisNorm.Normalize();
	return;
}

void PolyhedronSet::CalculateRadius()
{
	for (int i = 0; i < v_curvedFace.size(); i++)
	{
		if (1 == v_curvedFace[i].first)
		{
			Scalar denominator = 0.0;
			Scalar numerator = 0.0;
			for (int j = 0; j < v_facePointID[i].size(); j++)
			{
				int verticeID = v_facePointID[i][j];
				Vector OP = v_point[verticeID] - axisCenter;
				Scalar radius = (OP - (OP&axisNorm)*axisNorm).Mag();
				denominator += 1.0;
				numerator += radius;
			}
			Scalar averageRadius = numerator / denominator;
			v_curvedFace[i].second = averageRadius;
		}
	}
	return;
}

void PolyhedronSet::Display() const
{
	std::cout << "this is a MOC polyhedron, axis center and norm are" << std::endl;
	std::cout << this->axisCenter << std::endl;
	std::cout << this->axisNorm << std::endl;
	cout << "This MOC polyhedron is composed of points:" << endl;
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
	std::cout << "if any, curved faces are list below" << std::endl;
	for (int i = 0; i < this->v_curvedFace.size(); i++)
	{
		if (this->v_curvedFace[i].first)
		{
			std::cout << "radius of face " << i << " is " << v_curvedFace[i].second << std::endl;
		}
	}
	return;
}

PolyhedronSet PolyhedronSet::ClipByPlane(Vector pointOnPlane, Vector planeNorm)
{
	Vector clipNorm = planeNorm.GetNormal();
	vector<int> newFacelist;
	//create temperal data needed in the algorithm
	vector<bool> v_inside;
	vector<int> v_newPointID;
	//create a new polygon without any data
	PolyhedronSet targetMOCPoly;
	targetMOCPoly.axisCenter = this->axisCenter;
	targetMOCPoly.axisNorm = this->axisNorm;
	CreateOldVerticeMap(pointOnPlane, clipNorm, targetMOCPoly, v_inside, v_newPointID);

	//before creating the new polygon, check intersections in each face, and cut into triangles when needed
	if (true == CutWhenNeeded(v_inside))
	{
		//if cut occurs, operations on vertices need to be re-proceeded
		CreateOldVerticeMap(pointOnPlane, clipNorm, targetMOCPoly, v_inside, v_newPointID);
	}
	int numOldPoints = (int)targetMOCPoly.v_point.size();

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
				targetMOCPoly.v_facePointID.push_back(newface);
				v_oldFaceID.push_back(i);
				v_faceCut.push_back(false);
				break;
			}
		}
	}

	//re-order oldFaceID, curved face should be cut first
	std::vector<int> reOrderedOldFaceID;
	for (int i = 0; i < (int)v_oldFaceID.size(); i++)
	{
		if (1 == this->v_curvedFace[v_oldFaceID[i]].first) reOrderedOldFaceID.push_back(v_oldFaceID[i]);
	}
	for (int i = 0; i < (int)v_oldFaceID.size(); i++)
	{
		if (0 == this->v_curvedFace[v_oldFaceID[i]].first) reOrderedOldFaceID.push_back(v_oldFaceID[i]);
	}

	//Note: the three integers in v_intersection refer to
	//(1) the smaller and (2) the greater point ID of the new polyhedron and
	//(3) the point ID in the new polyhedron
	map<pair<int, int>, int> v_intersection;
	//visit all faces and insert points and faces into the new polyhedron
	for (int i = 0; i < (int)reOrderedOldFaceID.size(); i++)
	{
		int refFaceID = reOrderedOldFaceID[i];
		int faceIsCurved = this->v_curvedFace[refFaceID].first;
		Scalar faceRadius = this->v_curvedFace[refFaceID].second;
		std::vector<int>& refFace = this->v_facePointID[refFaceID];
		int numSides = (int)refFace.size();
		for (int j = 0; j < numSides; j++)
		{
			int thisID = refFace[j];
			int nextID = refFace[(j + 1) % numSides];
			if (v_inside[thisID])
			{
				targetMOCPoly.v_facePointID[i].push_back(v_newPointID[thisID]);
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
					
					if (1 == faceIsCurved)
					{
						//when the edge is on a curved face, all vertice need to be re-computed
						Vector OP = intersection - axisCenter;
						Vector OR = (OP & axisNorm) * axisNorm;
						Vector RCenter = axisCenter + OR;
						intersection = RCenter + faceRadius * (intersection - RCenter).GetNormal();
					}
					//end of correction
					int IDInNew = (int)targetMOCPoly.v_point.size();
					targetMOCPoly.v_point.push_back(intersection);
					targetMOCPoly.v_facePointID[i].push_back(IDInNew);
					v_intersection.insert(pair<pair<int, int>, int>(key, IDInNew));
				}
				else
				{
					targetMOCPoly.v_facePointID[i].push_back(it->second);
				}
			}
		}
		targetMOCPoly.v_curvedFace.push_back(std::pair<int, Scalar>(faceIsCurved, faceRadius));
	}

	//create the last polygon face
	vector<pair<int, int> > v_sideArror;
	for (int i = 0; i < (int)targetMOCPoly.v_facePointID.size(); i++)
	{
		if (false == v_faceCut[i])
		{
			continue;
		}
		vector<int>& thisFace = targetMOCPoly.v_facePointID[i];
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
		return targetMOCPoly;
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
		targetMOCPoly.v_facePointID.push_back(newface);
		targetMOCPoly.v_curvedFace.push_back(std::pair<int, Scalar>(0, 0));//the last face cannot be curved
		newFacelist.push_back((int)targetMOCPoly.v_facePointID.size() - 1);
	}
	return targetMOCPoly;
}

Tensor RotateTensor(Vector& axis, Scalar theeta)
{
	Tensor result(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	Scalar nx = axis.x_;
	Scalar ny = axis.y_;
	Scalar nz = axis.z_;
	Scalar costheeta = cos(theeta);
	Scalar sintheeta = sin(theeta);
	result.mat_[0][0] = costheeta + nx*nx*(1.0 - costheeta);
	result.mat_[0][1] = nx*ny*(1.0 - costheeta) - nz*sintheeta;
	result.mat_[0][2] = nx*nz*(1.0 - costheeta) + ny*sintheeta;
	result.mat_[1][0] = nx*ny*(1.0 - costheeta) + nz*sintheeta;
	result.mat_[1][1] = costheeta + ny*ny*(1.0 - costheeta);
	result.mat_[1][2] = ny*nz*(1 - costheeta) - nx*sintheeta;
	result.mat_[2][0] = nx*nz*(1.0 - costheeta) - ny*sintheeta;
	result.mat_[2][1] = ny*nz*(1 - costheeta) + nx*sintheeta;
	result.mat_[2][2] = costheeta + nz*nz*(1.0 - costheeta);
	return result;
}

void PolyhedronSet::ClipIntoSubPolygons(Scalar maxAngleInDegree)
{
	this->v_subPolyhedron.clear();
	if (false == this->geometryCalculated)
	{
		std::cout << "we cannot clip the polyhedron into sub polyhedrons, the geometry is not calculated" << std::endl;
	}
	//calculating the bounds of clipping normals
	Vector globalCenter = this->GetCenter();
	Vector temp = globalCenter - this->axisCenter;
	Vector nCenter = (temp - (temp&axisNorm)*axisNorm).GetNormal();
	Scalar valueMax = 0.0;
	Scalar valueMin = 0.0;
	Vector nMax = nCenter;
	Vector nMin = nCenter;
	for (int i = 0; i < this->v_point.size(); i++)
	{
		Vector temp = this->v_point[i] - this->axisCenter;
		Vector radialVector = temp - (temp & axisNorm) * axisNorm;
		//a special case when the point is located on the axis, and it should be skipped
		if (radialVector.Mag() < 10.0 * SMALL)
		{
			continue;
		}
		Vector nCandidate = (temp - (temp&axisNorm)*axisNorm).GetNormal();
		Scalar valueToCompare = (nCenter^nCandidate)&axisNorm;
		if (valueToCompare > valueMax)
		{
			valueMax = valueToCompare;
			nMax = nCandidate;
		}
		if (valueToCompare < valueMin)
		{
			valueMin = valueToCompare;
			nMin = nCandidate;
		}
	}
	Scalar theeta = acos(nMin & nMax);
	int number = (theeta / (PI*maxAngleInDegree / 180.0) + 1);
	Tensor RMatrix = RotateTensor(axisNorm, theeta / Scalar(number));
	//clipping starts from here
	Vector clipNorm = axisNorm ^ nMin;
	PolyhedronSet tempMOCPolygon = *this;
	for (int i = 1; i < number; i++)
	{
		clipNorm = RMatrix*clipNorm;
		MHT::Polyhedron sub = tempMOCPolygon.ClipByPlane(axisCenter, clipNorm);
		sub.MHT::Polyhedron::CalculateVolume();
		this->v_subPolyhedron.push_back(sub);
		tempMOCPolygon = tempMOCPolygon.ClipByPlane(axisCenter, -clipNorm);
	}
	tempMOCPolygon.MHT::Polyhedron::CalculateVolume();
	this->v_subPolyhedron.push_back((MHT::Polyhedron)tempMOCPolygon);
	this->CalculateVolume();
	this->clipped = true;
	return;
}

void PolyhedronSet::CalculateVolume()
{
	if (this->clipped)
	{
		this->volume = 0.0;
		for (int i = 0; i < this->v_subPolyhedron.size(); i++)
		{
			this->volume += v_subPolyhedron[i].GetVolume();
		}
	}
	return;
}

bool PolyhedronSet::IsContaining(Vector& point) const
{
	bool contain = true;
	for (int i = 0; i < this->v_faceCenter.size(); i++)
	{
		if (0 == this->v_curvedFace[i].first)
		{
			Vector faceCenter = this->v_faceCenter[i];
			Vector pointToFace = faceCenter - point;
			if ((pointToFace & this->v_faceArea[i]) < 0)
			{
				return false;
			}
		}
		else
		{
			Vector faceCenter = this->v_faceCenter[i];
			Scalar temp = (faceCenter - axisCenter)&this->v_faceArea[i];
			Vector OP = point - axisCenter;
			Scalar distanceToAxis = (OP - (OP&axisNorm)*axisNorm).Mag();
			if (temp > 0 && distanceToAxis > this->v_curvedFace[i].second)
			{
				return false;
			}
			if (temp < 0 && distanceToAxis < this->v_curvedFace[i].second)
			{
				return false;
			}
		}
	}
	return contain;
}

void PolyhedronSet::Move(Vector& translation)
{
	this->Polyhedron::Move(translation);
	//move the axis center
	this->axisCenter += translation;
	return;
}

PolyhedronSet PolyhedronSet::Copy(Vector& translation) const
{
	PolyhedronSet result = *this;
	result.Move(translation);
	if (this->clipped)
	{
		for (size_t i = 0;i < result.v_subPolyhedron.size();i++)
		{
			result.v_subPolyhedron[i].Move(translation);
		}
	}
	return result;
}

Scalar PolyhedronSet::IntersectionVolumeWithPolyhedron
(
	const MHT::Polyhedron& poly
) const
{
	Scalar volume = 0.0;
	if (0 == this->v_subPolyhedron.size())
	{
		std::cout << "In PolyhedronSet::IntersectionVolumeWithPolyhedron(), PolyhedronSet has not been clipped into polyhedrons" << std::endl;
		exit(1);
	}
	for (int i = 0; i < this->v_subPolyhedron.size(); i++)
	{
		MHT::Polyhedron temp = poly&&this->v_subPolyhedron[i];
		volume += temp.GetVolume();
	}
	return volume;
}

Scalar PolyhedronSet::IntersectionVolumeWithPolyhedronSet(const PolyhedronSet& another) const
{
	Scalar totalVolume = 0.0;
	if (false == this->clipped)
	{
		if (false == another.clipped)
		{
			MHT::Polyhedron result = (*this) && another;
			totalVolume = result.GetVolume();
		}
		else
		{
			for (int i = 0;i < another.v_subPolyhedron.size();i++)
			{
				MHT::Polyhedron result = (*this) && another.v_subPolyhedron[i];
				totalVolume += result.GetVolume();
			}
		}
	}
	else
	{
		if (false == another.clipped)
		{
			for (int i = 0; i < this->v_subPolyhedron.size(); i++)
			{
				MHT::Polyhedron result = this->v_subPolyhedron[i] && another;
				totalVolume += result.GetVolume();
			}
		}
		else
		{
			for (int i = 0; i < this->v_subPolyhedron.size(); i++)
			{
				for (int j = 0; j < another.v_subPolyhedron.size();j++)
				{
					MHT::Polyhedron result = this->v_subPolyhedron[i] && this->v_subPolyhedron[j];
					totalVolume += result.GetVolume();
				}
			}
		}
	}
	return totalVolume;
}

std::vector<Scalar> PolyhedronSet::GetRaduisList() const
{
	std::vector<Scalar> radiusList;
	for (size_t i = 0;i < this->v_curvedFace.size();i++)
	{
		if (this->v_curvedFace[i].first == 1)
		{
			radiusList.push_back(this->v_curvedFace[i].second);
		}
	}
	return radiusList;
}

std::pair<bool, Vector> PolyhedronSet::GetAxisCenter() const
{
	std::pair<bool, Vector> info;
	info.second = this->axisCenter;
	bool hasCurvdFace = false;
	for (size_t i = 0;i < this->v_curvedFace.size();i++)
	{
		if (1 == this->v_curvedFace[i].first)
		{
			hasCurvdFace = true;
			break;
		}
	}
	info.first = hasCurvdFace;
	return info;
}

std::vector<MHT::Polygon> PolyhedronSet::GetFacesOnBoxBoundary
(
	Vector verticeMin,//xmin, ymin, zmin
	Vector verticeMax,//xmax, ymax, zmax
	Scalar tolerance
) const
{
	if (!clipped)
	{
		return Polyhedron::GetFacesOnBoxBoundary(verticeMin, verticeMax, tolerance);
	}
	else
	{
		std::vector<MHT::Polygon> faceList;
		for (size_t i = 0;i < this->v_subPolyhedron.size();i++)
		{
			std::vector<MHT::Polygon> subFaceList = this->v_subPolyhedron[i].GetFacesOnBoxBoundary(verticeMin, verticeMax, tolerance);
			faceList.insert(faceList.end(), subFaceList.begin(), subFaceList.end());
		}
		return faceList;
	}
}

void PolyhedronSet::WriteTecplotFile(std::string filename) const
{
	if (false == this->clipped)
	{
		this->MHT::Polyhedron::WriteTecplotFile(filename);
		return;
	}
	else
	{
		ofstream outFile(filename);
		this->WriteTecplotHeader(outFile);
		this->WriteTecplotZones(outFile);
		outFile.close();
		return;
	}
}

void PolyhedronSet::WriteTecplotHeader(std::ofstream& ofile) const
{
	ofile << "TITLE =\"" << "polyhedron" << "\"" << endl;
	ofile << "VARIABLES = " << "\"x\"," << "\"y\"," << "\"z\"" << endl;
	return;
}

void PolyhedronSet::WriteTecplotZones(std::ofstream& ofile) const
{
	if (true == this->clipped)
	{
		for (int polyID = 0; polyID < this->v_subPolyhedron.size(); polyID++)
		{
			int nCount(0);
			int faceNum = (int)this->v_subPolyhedron[polyID].v_facePointID.size();
			for (int i = 0; i < faceNum; i++)
			{
				nCount += (int)this->v_subPolyhedron[polyID].v_facePointID[i].size();
			}
			ofile << "ZONE T = " << "\"" << "polyhedron" << "\",";
			ofile << " DATAPACKING = POINT, N = " << v_subPolyhedron[polyID].v_point.size();
			ofile << ", E = " << nCount << ", ZONETYPE = FELINESEG" << std::endl;

			for (unsigned int i = 0; i < v_subPolyhedron[polyID].v_point.size(); i++)
			{
				Vector vertice = v_subPolyhedron[polyID].v_point[i];
				ofile << setprecision(8) << setiosflags(ios::scientific) << vertice.x_ << " " << vertice.y_ << " " << vertice.z_ << endl;
			}

			for (int i = 0; i < faceNum; i++)
			{
				for (int n = 0; n < (int)this->v_subPolyhedron[polyID].v_facePointID[i].size(); n++)
				{
					if (n == (int)this->v_subPolyhedron[polyID].v_facePointID[i].size() - 1)
					{
						ofile << v_subPolyhedron[polyID].v_facePointID[i][n] + 1 << "\t" << v_subPolyhedron[polyID].v_facePointID[i][0] + 1 << std::endl;
					}
					else
					{
						ofile << v_subPolyhedron[polyID].v_facePointID[i][n] + 1 << "\t" << v_subPolyhedron[polyID].v_facePointID[i][n + 1] + 1 << std::endl;
					}
				}
			}
		}
	}
	else
	{
		this->MHT::Polyhedron::WriteTecplotZones(ofile);
	}
	return;

}
/*
void PolyhedronSet::WriteTecplotZones(std::ofstream& ofile) const
{
	for (int polyID = 0; polyID < this->v_subPolyhedron.size(); polyID++)
	{
		int nCount(0);
		int faceNum = (int)this->v_subPolyhedron[polyID].v_facePointID.size();
		for (int i = 0; i < faceNum; i++)
		{
			nCount += (int)this->v_subPolyhedron[polyID].v_facePointID[i].size();
		}
		ofile << "ZONE T = " << "\"" << "polyhedron" << "\",";
		ofile << " DATAPACKING = POINT, N = " << v_subPolyhedron[polyID].v_point.size();
		ofile << ", E = " << nCount << ", ZONETYPE = FELINESEG" << std::endl;

		for (unsigned int i = 0; i < v_subPolyhedron[polyID].v_point.size(); i++)
		{
			Vector vertice = v_subPolyhedron[polyID].v_point[i];
			ofile << setprecision(8) << setiosflags(ios::scientific) << vertice.x_ << " " << vertice.y_ << " " << vertice.z_ << endl;
		}

		for (int i = 0; i < faceNum; i++)
		{
			for (int n = 0; n < (int)this->v_subPolyhedron[polyID].v_facePointID[i].size(); n++)
			{
				if (n == (int)this->v_subPolyhedron[polyID].v_facePointID[i].size() - 1)
				{
					ofile << v_subPolyhedron[polyID].v_facePointID[i][n] + 1 << "\t" << v_subPolyhedron[polyID].v_facePointID[i][0] + 1 << std::endl;
				}
				else
				{
					ofile << v_subPolyhedron[polyID].v_facePointID[i][n] + 1 << "\t" << v_subPolyhedron[polyID].v_facePointID[i][n + 1] + 1 << std::endl;
				}
			}
		}
	}
	return;

}
*/