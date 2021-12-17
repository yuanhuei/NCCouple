/*---------------------------------------------------------------------------*\
File Name:
	FluentMeshBlock.h

Description:
	Derivation of Mesh to make it suited to FLUENT mesh file

	Author:		Shuai Zhang, Kong Ling
	Date: once upon a time (I can't remember)

Revised:
	Description:
	1. Display functions were removed;
	2. Styles of information display during mesh reading were entirely checked and revised,
	where logLevels and systemControl were applied.

	Revisor:		Kong Ling, Shuai Zhang
	Modified Date:	2017-01-14
\*---------------------------------------------------------------------------*/

#ifndef _FluentMeshBlock_
#define _FluentMeshBlock_

#include <sstream>

#include "../MHT_common/Configuration.h"
#include "../MHT_mesh/Mesh.h"

class StandardMeshBlock;
class RegionConnection;

class MshElementZone
{
public:

    std::string st_name;

	int	n_zoneID;
	int n_elemNum;
	int n_begID, n_endID;

	int n_active;
	int n_elemShape;

    Element::ElemShapeType est_shapeType;

    MshElementZone()
        :n_zoneID(0),
          n_elemNum(0),
          n_begID(0),
          n_endID(0),
          n_active(0),
          n_elemShape(0),
          est_shapeType(Element::estNoType)
    {}
};

class MshFaceZone
{
public:
    std::string st_name;
	int n_zoneID;
	int n_faceNum;
	int n_begID, n_endID;

	int n_BCType;
	int n_faceShape;

	MshFaceZone() : n_zoneID(0), n_begID(0), n_endID(0){}

	Face::FaceType IntToFaceType(int id);
};

class FluentMeshFileSection
{
public:
	int n_index;
    std::vector<int> v_parameters;
    std::stringstream ss_details;
public:
    FluentMeshFileSection();

    int ReadNew(std::ifstream& inFile);
};

class FluentMeshBlock :public Mesh
{
public:
	// These two parameters is for multiple region mesh
    std::vector<ElementZone> v_elementZone;					//Vector of element zone
    std::vector<FaceZone> v_faceZone;						//vector of face zones

	//composing region-grid, created when decomposing grids
    std::vector<StandardMeshBlock> v_regionGrid;

	//blocks of information in a msh file
    std::vector<MshFaceZone> v_FaceZoneInfo;
    std::vector<MshElementZone> v_ElementZoneInfo;

    std::vector<int> v_parameters;
public:

    FluentMeshBlock(const std::string& MshFileName);

    void ReadMesh(std::ifstream& inFile);
	void ReadMesh(FILE *F_FilePointer);

	int GetIndex(FILE *F_FilePointer);
    int GetData(FILE *F_FilePointer, std::string &s_InputData);
    Scalar GetData(FILE *F_FilePointer, std::string &s_InputData, int n_DataType);
	void ReadNew(FILE *F_FilePointer, int n_Index);

	void ReadComment(FluentMeshFileSection& inSection);
	void ReadComment(FILE *F_FilePointer);
	void ReadDimensionNumber(FluentMeshFileSection& inSection);
	void ReadDimensionNumber(FILE *F_FilePointer);
	void ReadNode(FluentMeshFileSection& inSection);
	void ReadNode(FILE *F_FilePointer);
	void ReadFace(FluentMeshFileSection& inSection);
	void ReadFace(FILE *F_FilePointer);
	void ReadElement(FluentMeshFileSection& inSection);
	void ReadElement(FILE *F_FilePointer);
	void ReadZone(FluentMeshFileSection& inSection);
	void ReadZone(FILE *F_FilePointer);

	//Build up element and face zones in standard forms according to the Element- and Face- Zone Infos.
	void EstablishZones();

	//Convert Origin Label (1:n) To (0:n-1)
	void RefineArrayIndex();
	 
	//Get Cell Node
	void PopulateCellNodes();
	void PopulateTriangleCell(int i);
	void PopulateTetraCell(int i);
	void PopulateQuadCell(int i);
	void PopulateHexahedronCell(int i);
	void PopulatePyramidCell(int i);
	void PopulateWedgeCell(int i);
	void PopulatePolyhedronCell(int i);

    void WriteTecplotMesh(const std::string& outMshFileName);

	//Decompose the global grid into regions, and write connection information into given RegionConnection
	void Decompose(RegionConnection&);
	/*===============Member of Decompose===================================*/
	//2. get zone number
	void GetZoneNumber(int &zoneNum);
	//3. distribute elements
	void DistributeElements(std::vector<std::pair<int, int> > &globalElementIndex);
	//4. distribute faces
	void DistributeFaces(std::vector<std::vector<std::pair<int, int> > > &globalFaceIndex);
	//5. distribute nodes & vertices
	void DistributeNodes_Vertices(std::vector<std::vector<std::pair<int, int> > > &globalNodeIndex);
	//6. rewrite topological information
	void RewriteTopologicalInformation(std::vector<std::vector<std::pair<int, int> > > &globalNodeIndex, std::vector<std::vector<std::pair<int, int> > > &globalFaceIndex, std::vector<std::pair<int, int> > &globalElementIndex);
	//7. distribute interior face zones
	void DistributeInteriorFaceZones(std::vector<std::vector<std::pair<int, int> > > &globalFaceIndex);
	//8. split external (one-side) face zones
	void SplitExternalFaceZones(std::vector<std::vector<std::pair<int, int> > > &globalFaceIndex);
	//9. split internal (two-side) face zones
	void SplitInternalFaceZones(std::vector<std::vector<std::pair<int, int> > > &globalFaceIndex, std::vector<std::pair<int, int> > &globalElementIndex, RegionConnection& regionConnectionInfo);
	//10. side check & correction for newly created faces at boundary
	void SideCheck();
	//11. write n_FaceZoneID for faces in each region grid
	void WriteFaceID();
	/*=====================================================================*/


	//Clear Global Fluent mesh
	void ClearGlobalFluentMesh();
};


#endif
