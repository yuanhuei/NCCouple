
#include "../MHT_mesh/MeshSupport.h"
#include "../MHT_common/LogLevel.h"
#include "../MHT_common/GeometryTools.h"

//--------------------------------Member functions of class Face-----------------------------------
//Constructor
Face::Face(int Node1, int Node2, int Node3)
{
	this->v_nodeID.push_back(Node1);
	this->v_nodeID.push_back(Node2);
	this->v_nodeID.push_back(Node3);
}

//Calculate center and area with a given set of nodes	
void Face::CalculateCenterAndArea(std::vector<Node>& gridNode)
{
    std::vector<Node> nodeList;
	for (int i = 0; i < (int)v_nodeID.size(); i++)
	{
		nodeList.push_back(gridNode[v_nodeID[i]]);
	}
	GetCenterAndArea(nodeList, center, area);
}

//correct face sides in case of incorrect indications for owner and neighbor cells
void Face::CorrectSide()
{
	if (-1 == this->n_owner)
	{
		//1. exchange IDs of owner and neighbor elements
		this->n_owner = this->n_neighbor;
		this->n_neighbor = -1;
		//2. reverse surface area (Note: it is a vector rather than a scalar)
		this->area = -(this->area);
		//3. reverse sequence of nodeID
        std::vector<int> temp = this->v_nodeID;
		this->v_nodeID.assign(temp.rbegin(), temp.rend());
	}
	return;
}

//--------------------------------Member functions of class Element-----------------------------------

//Calculate Center And Volume of a 2D element using polygon Geometry
void Element::CalculateCenterAndVolume(std::vector<Node>& gridNode)
{
    std::vector<Node> verticeList;
	//building up vertice list of this 2D element;
	for (int j = 0; j < (int)v_nodeID.size(); j++)
	{
		verticeList.push_back(gridNode[v_nodeID[j]]);
	}
	//calculating center point and volume of this 2D elment;
	Vector elementArea;
	GetCenterAndArea(verticeList, center, elementArea);
	volume = elementArea.Mag();
}

//Calculate Center And Volume of a 3D element using polyhedron Geometry
void Element::CalculateCenterAndVolume(std::vector<Face>& gridFace)
{
    std::vector<Node> faceCenterList, faceAreaList;
	//building up face list of this 3D element;
	for (int j = 0; j < (int)v_faceID.size(); j++)
	{
		faceCenterList.push_back(gridFace[v_faceID[j]].center);
		faceAreaList.push_back(gridFace[v_faceID[j]].area);
	}
	//calculating center point and volume of this 3D elment;
	GetCenterAndVolume(faceCenterList, faceAreaList, center, volume);
}

//--------------------------------Member functions of class FaceZone-----------------------------------
std::string FaceZone::IntToBCDescription(int id)
{
	switch (id)
	{
	case 2:
		return "Interior";
		break;
	case 3:
		return "Wall";
		break;
	case 4:
		return "Pressure Inlet";
		break;
	case 5:
		return "Pressure Outlet";
		break;
	case 7:
		return "Symmetry";
		break;
	case 8:
		return "Periodic Shadow";
		break;
	case 9:
		return "Pressure far field";
		break;
	case 10:
		return "Velocity Inlet";
		break;
	case 12:
		return "Periodic";
		break;
	case 14:
		return "Fan";
		break;
	case 20:
		return "Mass Flow Inlet";
		break;
	case 24:
		return "Interface";
		break;
	case 31:
		return "Parent";
		break;
	case 36:
		return "Outflow";
		break;
	case 37:
		return "Axis";
		break;
	default:
		std::cout << LogLevel(LogLevel::mlError, "Error in input parameter in function IntToBCDescription(), no type ") << id << std::endl;
		return "No Boundary Type";
		break;
	}
}

//--------------------------------External Functions-----------------------------------

//Calculating the difference between two given 2D faces according to their centers and areas
Scalar ErrorBetween2DFaces(Face& face1, Face&face2)
{
	Scalar denominator = face1.area.Mag();
	Scalar errorAreas = (face1.area - face2.area).Mag()/denominator;
	Scalar errorCenters = (face1.center - face2.center).Mag() / denominator;
	return errorAreas + errorCenters;
}

//Calculating the difference between two given 3D faces according to their centers and areas
Scalar ErrorBetween3DFaces(Face& face1, Face& face2)
{
	Scalar denominator = face1.area.Mag();
	Scalar errorAreas = (face1.area - face2.area).Mag() / denominator;
	Scalar centerDistance = (face1.center - face2.center).Mag();
	Scalar errorCenters = centerDistance * centerDistance / denominator;
	return errorAreas + errorCenters;
}
