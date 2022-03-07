
#ifndef _MeshSupport_
#define _MeshSupport_

#include <vector>
#include <string>

#include "../MHT_common/Configuration.h"
#include "../MHT_common/Vector.h"
#include "../MHT_common/Tensor.h"

typedef Point3D<Scalar> Node;

class Element;

//topological information for nodes
class Vertice
{
public:
	//IDs of elements associated with this vertice
    std::vector<int> v_elemID;
	//IDs of faces associated with this vertice
    std::vector<int> v_faceID;
	//constructor
	Vertice(){}
};

class Face
{
public:
	enum FaceType
	{
		ftMixed = 0,
		ftLinear = 2,
		ftTrangular = 3,
		ftQuadrilateral = 4,
		ftPolygon = 5,
		ftNoType
	};

public:

	FaceType	ft_faceType;

	//indicating which boundary this face is located on, valued as -1 for interior faces
	int n_FaceZoneID;

	//Face Center;
	Vector center;

	//Area of this face (Note: the "area" is a vector pointing to the normal direction)
	Vector area;

    std::vector<int> v_nodeID;

	std::vector<Vector> v_node;

	//element ID on one side of the face 
	int		n_owner;
	//element ID on the other side of the face, equaling to -1 for a boundary face
	int		n_neighbor;

public:
	//Constructor
	Face(){}

	//Constructor
	Face(int Node1, int Node2, int Node3);

	//operator comparing whether given two face are the same one
	inline friend bool operator == (const Face& f1, const Face& f2)
	{
		if (f1.v_nodeID.size() != f2.v_nodeID.size())
		{
			return false;
		}
		for (int i = 0; i < (int)f1.v_nodeID.size(); i++)
		{
			if (f1.v_nodeID[i] != f2.v_nodeID[i])
			{
				return false;
			}
		}
		return true;
	}

	//Calculate center and area with a given set of nodes
    void CalculateCenterAndArea(std::vector<Node>&);

	//correct face sides in case of incorrect indications for owner and neighbor cells
	void CorrectSide();
};

class Element
{
public:
	enum ElemType
	{
		real,
		ghost
	};

	enum ElemShapeType
	{
		estMixed			= 0,
		estTriangular		= 1,
		estTetrahedral		= 2,
		estQuadrilateral	= 3,
		estHexahedral		= 4,
		estPyramid			= 5,
		estWedge			= 6,
		estPolyhedron		= 7,
		estPolygon = 8,
		estNoType
	};

public:

	ElemShapeType	est_shapeType;

	//Element type, real or ghost
	ElemType	et_type;

	int n_ElementZoneID;

	//Element Center;
	Vector center;

	//Volume of this element. 
	//The area will be considered as the "volume" for two-dimensional element;
	Scalar volume;

	//Face ID in Element
    std::vector<int> v_faceID;

	//Node ID in Element
    std::vector<int> v_nodeID;

	//Calculate Center And Volume of a 2D element using polygon Geometry
    void CalculateCenterAndVolume(std::vector<Node>& gridNode);

	//Calculate Center And Volume of a 3D element using polyhedron Geometry
    void CalculateCenterAndVolume(std::vector<Face>& gridFace);
};

class ElementZone
{
public:
    std::string	name;
    std::vector<int> v_elementID;
    std::vector<int> v_subFaceZoneID;
	ElementZone()
	{
        std::vector<int>().swap(this->v_elementID);
	}
};

class FaceZone
{
public:
	//Enum for face zone type
	enum FZType
	{
		fztInterior = 0,
		fztExternal = 1,
		fztInternal = 2,
		fztMixed = 3,
		fztNoType
	};

	//Enum for physical BC
	enum BCType
	{
		bcInterior = 2,
		bcWall = 3,
		bcPressureInlet = 4,
		bcPressureOutlet = 5,
		bcSymmetry = 7,
		bcPeriodicShadow = 8,
		bcPressureFarfield = 9,
		bcVelocityInlet = 10,
		bcPeriodic = 12,
		bcFan = 14,
		bcMassFlowInlet = 20,
		bcInterface = 24,
		bcParent = 31,
		bcOutflow = 36,
		bcAxis = 37,
		bcNoType = 0
	};

	int n_patchID;
    std::string	name;
    std::vector<int> v_faceID;
	BCType bc_Type;
	FZType faceZoneType;

	FaceZone()
		:
		n_patchID(0),
		bc_Type(bcNoType),
		faceZoneType(fztNoType)
	{}

	BCType IntToBCType(int id);
    static std::string IntToBCDescription(int id);
};

//Calculating the difference between two given 2D faces according to their centers and areas
Scalar ErrorBetween2DFaces(Face&,Face&);

//Calculating the difference between two given 3D faces according to their centers and areas
Scalar ErrorBetween3DFaces(Face&, Face&);

#endif
