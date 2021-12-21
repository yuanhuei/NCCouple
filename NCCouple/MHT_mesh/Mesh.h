#pragma once
#ifndef _Mesh_
#define _Mesh_

#include <string>
#include "../MHT_common/Configuration.h"
#include "../MHT_mesh/MeshSupport.h"

class Mesh
{
public:

	enum MeshDim
	{
		md2D	= 2,
		md3D	= 3,
		mdNoType
	};

    std::string	st_fileName;								//File name
    std::string	st_meshName;								//Mesh name

	int	n_elemNum;										//Total element number
	int	n_faceNum;										//Total face number
	int	n_nodeNum;										//Total node number

    MeshDim	md_meshDim;

    Element::ElemShapeType est_shapeType;

    std::vector<Element>	v_elem;								//Vector of Element
    std::vector<Face>	v_face;								//Vector of Face
    std::vector<Node>	v_node;								//Vector of Node
    std::vector<Vertice> v_vertice;							//Vector of Vertice 

	FaceZone fz_interiorFaceZone;						//information of interior faces 
    std::vector<FaceZone> v_boundaryFaceZone;				//information of boundary faces

public:

	Mesh();

    Mesh(const std::string& MshFileName);

    virtual void ReadMesh(std::ifstream& inFile) = 0;

    virtual void WriteTecplotMesh(const std::string& outMshFileName) = 0;
	
	//Write IDs of elements and faces associated to a vertice such that they can be visited directly
	void EstablishVertice();

	//Calculate numbers of nodes, faces and elements
	void CalculateParameters();

	//Calculate centers of faces & elements, areas of faces, and volumes of elements 
	void CalculateCenterAndVolume();

    std::pair<int, int> SearchOwnerNeighbor(int faceID) const;

	//Search the IDs of elements in the same element zone sharing faces with a element whose ID is given
    std::vector<int> SearchElementNeighbor(int EID) const;

	//given the ID of a boundary face, search IDs of its neighbor located on the same boundary
    std::vector<int> SearchBounaryFaceNeighbor(int) const;

	//search of the IDs of element within a given distance to a specific element 
    std::vector<int> SearchNearbyElements(int, Scalar) const;

	//search of the IDs of element within a given number of layers
    std::vector<int> SearchNearbyElements(int, int) const;

	//calculate to get the characteristic length of this mesh
    Scalar GetCharacteristicLength() const;

	//Scaling of the entrire grid
	void Scale(Scalar);

	//Translate the entire mesh with a given vector
	void Translate(Vector);

	//Rotate the entire mesh with a given center, an axis and a rotating angle
	void Rotate(Vector, Vector, Scalar);

	//Set a boundary type on a faceZone with a specific name
    void SetBoundaryType(const std::string& patchName, FaceZone::BCType bc_Type);

	//Search a faceZone with a given name
    int GetBoundaryFaceZoneID(const std::string& patchName) const;

    Scalar GetBoundaryArea(const std::string&) const;

    void CheckBoundaryCondition() const;

	void CorrectSide();

};

#endif
