/*---------------------------------------------------------------------------*\
Class:
    Mesh
File Name:
    Mesh.h
Description:
    Definition of standard grid structure

    Author:		Shuai Zhang, Kong Ling
    Date: once upon a time (I can't remember)

Revised:
    Description:
    1. Display functions for element- and face-zones were revised to standard forms

    Revisor:		Kong Ling
    Modified Date:	2017-01-14

Revised:
    Description:
    1. Convert "non const" input parameters to const;
    2. Convert "non const" function to const function;

    Revisor:		Shuai Zhang
    Modified Date:	2018-11-28
\*---------------------------------------------------------------------------*/

#include <sstream>
#include "../MHT_mesh/Mesh.h"
#include "../MHT_common/LogLevel.h"
#include "../MHT_common/SystemControl.h"
#include "../MHT_common/ListTools.h"
#include "../MHT_common/GeometryTools.h"

Mesh::Mesh()
	:
	n_elemNum(0),
	n_faceNum(0),
	n_nodeNum(0),
	md_meshDim(mdNoType),
	est_shapeType(Element::estNoType)
{}

Mesh::Mesh(const std::string& MshFileName)
	:st_fileName(MshFileName),
	n_elemNum(0),
	n_faceNum(0),
	n_nodeNum(0),
	md_meshDim(mdNoType),
	est_shapeType(Element::estNoType)
{}

void Mesh::EstablishVertice()
{
	int nodeNum = (int)this->v_node.size();
    this->v_vertice.resize(nodeNum);
	for (int i = 0; i < (int)this->v_elem.size(); i++)
    {
		//visiting one element
		for (int j = 0; j < (int)v_elem[i].v_nodeID.size(); j++)
        {
            int verticeNum = v_elem[i].v_nodeID[j];
			this->v_vertice[verticeNum].v_elemID.push_back(i);
		}
    }
	for (int i = 0; i < (int)this->v_face.size(); i++)
	{
		//visiting one face
		for (int j = 0; j < (int)v_face[i].v_nodeID.size(); j++)
		{
			int verticeNum = v_face[i].v_nodeID[j];
			this->v_vertice[verticeNum].v_faceID.push_back(i);
		}
	}
}

void Mesh::CalculateParameters()
{
	this->n_elemNum = (int)this->v_elem.size();
	this->n_faceNum = (int)this->v_face.size();
	this->n_nodeNum = (int)this->v_node.size();
}

void Mesh::CalculateCenterAndVolume()
{
	//Step 1: Get centers and areas for faces;
	for (int i = 0; i < this->n_faceNum; i++)
	{
		this->v_face[i].CalculateCenterAndArea(this->v_node);
	}
	//Step 2: Get centers and volumes for elements;
	if (this->md_meshDim == md2D)
	{
		for (int i = 0; i < this->n_elemNum; i++)
		{
			this->v_elem[i].CalculateCenterAndVolume(this->v_node);
		}
	}
	else if (this->md_meshDim == md3D)
	{
		for (int i = 0; i < this->n_elemNum; i++)
		{
			this->v_elem[i].CalculateCenterAndVolume(this->v_face);
		}
	}
	//Step 3: Correct directions for face areas
	for (int i = 0; i < this->n_faceNum; i++)
	{
		std::pair<int, int> OandN = this->SearchOwnerNeighbor(i);
		int ownerID = OandN.first;
		int neighborID = OandN.second;
		Node ownerCenter = this->v_elem[ownerID].center;
		Node neighborCenter(0.0, 0.0, 0.0);
		if (-1 == neighborID)
		{
			neighborCenter = v_face[i].center;
		}
		else
		{
			neighborCenter = this->v_elem[neighborID].center;
		}
		Node OtoN = neighborCenter - ownerCenter;
		if ((v_face[i].area & OtoN) < 0.0)
		{
			v_face[i].area = -v_face[i].area;
		}
	}
}

std::pair<int, int> Mesh::SearchOwnerNeighbor(int faceID) const
{
	if (faceID<0 || faceID>this->n_faceNum)
	{
		FatalError("face ID out of range in Mesh::SearchOwnerNeighbor()");
	}
	int ownerID = v_face[faceID].n_owner;
	int neighborID = v_face[faceID].n_neighbor;
    return std::pair<int, int>(ownerID, neighborID);
}

//Search the IDs of elements in the same element zone sharing faces with a element whose ID is given
std::vector<int> Mesh::SearchElementNeighbor(int EID) const
{
	if (EID<0 || EID>=this->n_elemNum)
	{
		FatalError("element ID out of range in Mesh::SearchElementNeighbor!");
	}
	int zoneID = v_elem[EID].n_ElementZoneID;
    std::vector<int> neighborList;
	for (int i = 0; i < (int)v_elem[EID].v_faceID.size(); i++)
	//visiting all faces bounding the element
	{
		int faceID = v_elem[EID].v_faceID[i];
		std::pair<int, int> oAndN = SearchOwnerNeighbor(faceID);
		if ((-1 != oAndN.first) && (EID != oAndN.first))
		{
			if (zoneID == v_elem[oAndN.first].n_ElementZoneID)
			{
				neighborList.push_back(oAndN.first);
			}
		}
		if ((-1 != oAndN.second) && (EID != oAndN.second))
		{
			if (zoneID == v_elem[oAndN.second].n_ElementZoneID)
			{
				neighborList.push_back(oAndN.second);
			}
		}
	}
	if (0 == (int)neighborList.size())
	{
        std::stringstream ss;
		ss << "element " << EID << " is found to be isolated by Mesh::SearchElementNeighbor";
		WarningContinue(ss.str());
	}
	return neighborList;
}

std::vector<int> Mesh::SearchBounaryFaceNeighbor(int faceID) const
{
	if (faceID<0 || faceID>=this->n_faceNum)
	{
		FatalError("face ID out of range in Mesh::SearchBounaryFaceNeighbor!");
	}
	int bondaryID = this->v_face[faceID].n_FaceZoneID;
	if (-1 == bondaryID)
	{
		FatalError("Cannot proceed Mesh::SearchBounaryFaceNeighbor, the given face is an interior one!");
	}
    std::vector<int> faceIDList;
	for (int i = 0; i < (int)this->v_face[faceID].v_nodeID.size(); i++)
	{
		int nodeID = this->v_face[faceID].v_nodeID[i];
		for (int j = 0; j < (int)this->v_vertice[nodeID].v_faceID.size(); j++)
		{
			int candiFaceID = this->v_vertice[nodeID].v_faceID[j];
			if (candiFaceID == faceID)
			{
				continue;
			}
			if (bondaryID != this->v_face[candiFaceID].n_FaceZoneID)
			{
				continue;
			}
			faceIDList.push_back(candiFaceID);
		}
	}
	return GetNonRepeatedList(faceIDList);
}

std::vector<int> Mesh::SearchNearbyElements(int EID, Scalar dis) const
{
	if (EID<0 || EID >= this->n_elemNum)
	{
		FatalError("element ID out of range in Mesh::SearchElementNeighbor!");
	}
    std::map<int, Scalar> elemAndDis;
    std::vector<int> v_front;
    std::vector<int> v_directNb = this->SearchElementNeighbor(EID);
	for (int i = 0; i < (int)v_directNb.size(); i++)
	{
		int thisID = v_directNb[i];
		Scalar disTemp = (this->v_elem[thisID].center - this->v_elem[EID].center).Mag();
		if (disTemp < dis)
		{
            elemAndDis.insert(std::pair<int, Scalar>(thisID, disTemp));
			v_front.push_back(thisID);
		}
	}
	while (!v_front.empty())
	{
        std::vector<int> v_newFront;
		for (int i = 0; i < (int)v_front.size(); i++)
		{
            std::vector<int> nbIDs = this->SearchElementNeighbor(v_front[i]);
			for (int j = 0; j < (int)nbIDs.size(); j++)
			{
				int thisID = nbIDs[j];
				if (thisID == EID) continue;
				if (elemAndDis.find(nbIDs[j]) == elemAndDis.end())
				{
					Scalar disTemp = (this->v_elem[thisID].center - this->v_elem[EID].center).Mag();
					if (disTemp < dis)
					{
                        elemAndDis.insert(std::pair<int, Scalar>(thisID, disTemp));
						v_newFront.push_back(thisID);
					}
				}
			}
		}
		v_front = GetNonRepeatedList(v_newFront);
	}
    std::vector<int> v_result;
	v_result.resize(elemAndDis.size());
    std::map<int, Scalar>::iterator it = elemAndDis.begin();
	int i = 0;
	while (it != elemAndDis.end())
	{
		v_result[i] = it->first;
		i++;
		it++;
	}
	return v_result;
}

std::vector<int> Mesh::SearchNearbyElements(int EID, int layerNum) const
{
	if (EID<0 || EID >= this->n_elemNum)
	{
		FatalError("element ID out of range in Mesh::SearchElementNeighbor(int,int)!");
	}
	if (layerNum < 1)
	{
		FatalError("layer number should be at least 1 in Mesh::SearchElementNeighbor(int,int)!");
	}
	if (1 == layerNum)
	{
        std::vector<int> nbElemIDList;
        for (int i = 0; i < (int)this->v_elem[EID].v_nodeID.size(); i++)
		{
			int nodeID = this->v_elem[EID].v_nodeID[i];
            for (int j = 0; j < (int)this->v_vertice[nodeID].v_elemID.size(); j++)
			{
				int nbElemID = this->v_vertice[nodeID].v_elemID[j];
				if (nbElemID != EID) nbElemIDList.push_back(nbElemID);
			}
		}
		return GetNonRepeatedList(nbElemIDList);
	}
	else
	{
		Scalar gridScale = this->GetCharacteristicLength();
		return SearchNearbyElements(EID, (Scalar)layerNum*gridScale);
	}
}

//calculate to get the characteristic length of this mesh
Scalar Mesh::GetCharacteristicLength() const
{
	Scalar averageVol = 0.0;
	for (int i = 0; i < (int)this->v_elem.size(); i++)
	{
		averageVol += this->v_elem[i].volume;
	}
	averageVol = averageVol / (Scalar)this->n_elemNum;
	Scalar gridScale = 0.0;
	if (this->md_meshDim == Mesh::md2D)
	{
		gridScale = sqrt(averageVol);
	}
	else if (this->md_meshDim == Mesh::md3D)
	{
		gridScale = pow(averageVol, 1.0 / 3.0);
	}
	return gridScale;
}

//Translate the entire mesh with a given vector
void Mesh::Translate(Vector transVector)
{
	for (int i = 0; i < this->n_nodeNum; i++)
	{
		this->v_node[i] = this->v_node[i]+transVector;
	}
	this->CalculateCenterAndVolume();
}

//Scaling of the entrire grid
void Mesh::Scale(Scalar factor)
{
	for (int i = 0; i < this->n_nodeNum; i++)
	{
		this->v_node[i] *= factor;
	}
	this->CalculateCenterAndVolume();
}

//Rotate the entire mesh with a given center, an axis and a rotating angle
void Mesh::Rotate(
	Vector AxisPoint,
	Vector axis,
	Scalar theeta
	)
{
	for (int i = 0; i < this->n_nodeNum; i++)
	{
		this->v_node[i] = RotatePoint(this->v_node[i], AxisPoint, axis, theeta);
	}
	this->CalculateCenterAndVolume();
}

void Mesh::SetBoundaryType(const std::string& patchName, FaceZone::BCType bc_Type)
{
	int found = -1;
	for (int i = 0; i < (int)this->v_boundaryFaceZone.size(); i++)
	{
		if (patchName == this->v_boundaryFaceZone[i].name)
		{
			found = i;
			break;
		}
	}
	if (-1 == found)
	{
		WarningContinue(patchName + " not found in boundary list");
	}
	else
	{
		this->v_boundaryFaceZone[found].bc_Type = bc_Type;
		std::cout << "Boundary condition successfully given on boundary: " << patchName << "(" << FaceZone::IntToBCDescription((int)this->v_boundaryFaceZone[found].bc_Type) << ")" << std::endl;
	}
}

int Mesh::GetBoundaryFaceZoneID(const std::string& patchName) const
{
	int found = -1;
	for (int i = 0; i < (int)this->v_boundaryFaceZone.size(); i++)
	{
		if (patchName == this->v_boundaryFaceZone[i].name)
		{
			found = i;
			break;
		}
	}
	if (-1 == found)
	{
		WarningContinue(patchName + " not found in boundary list");
	}
	return found;
}

Scalar Mesh::GetBoundaryArea(const std::string& patchName) const
{
	int boundaryID = this->GetBoundaryFaceZoneID(patchName);
	if (-1 == boundaryID)
	{
		FatalError("Cannot calculate boundary area for not finding boudnary "+patchName);
	}
	Scalar result = 0.0;
	for (unsigned int i = 0; i < this->v_boundaryFaceZone[boundaryID].v_faceID.size(); i++)
	{
		int faceID = this->v_boundaryFaceZone[boundaryID].v_faceID[i];
		result += this->v_face[faceID].area.Mag();
	}
	return result;
}

void Mesh::CheckBoundaryCondition() const
{
	std::cout << "Checking physical boundary conditions ..." << std::endl;

	for (int i = 0; i < (int)this->v_boundaryFaceZone.size(); i++)
	{
		if (this->v_boundaryFaceZone[i].bc_Type == FaceZone::bcNoType)
		{
            FatalError(std::string("Physical boundary condition not given on External face: ") + this->v_boundaryFaceZone[i].name);
		}
	}

	std::cout << LogLevel(LogLevel::mlOK, "Mesh:") << this->st_meshName << ", physical boundary condition OK." << std::endl;
}
