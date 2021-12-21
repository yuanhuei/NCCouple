/*---------------------------------------------------------------------------*\
File Name:
StandardMeshBlock.h

Description:
	1.Derivation of Mesh to make it suited to StandardMeshBlock mesh file.
	2.This Class is made for our own CFD code.
	3.This Class is single mesh for unstruct grid (eg. Fluent mesh with multiple could decompose into several UnGridBlock). 

Author:		Shuai Zhang, Kong Ling
Date: 2017-02-16

\*---------------------------------------------------------------------------*/

#include "../MHT_mesh/StandardMeshBlock.h"
#include <iomanip>
#include "../MHT_common/LogLevel.h"
#include "../MHT_common/SystemControl.h"

StandardMeshBlock::StandardMeshBlock(const std::string& MshFileName)
	:Mesh(MshFileName)
{
    std::ifstream	inFile(MshFileName.c_str());

	if (inFile.fail())
	{
        FatalError(std::string("Can not find Mesh file:\n") + MshFileName);
	}

	ReadMesh(inFile);
}

void StandardMeshBlock::ReadMesh(std::ifstream& inFile)
{
    std::string stTemp;

	// Read Head
	inFile >> this->n_nodeNum;
	this->v_node.resize(this->n_nodeNum);

	inFile >> this->n_faceNum;
	this->v_face.resize(this->n_faceNum);

	inFile >> this->n_elemNum;
	this->v_elem.resize(this->n_elemNum);

	inFile >> this->nFaceNodeNum;

    std::cout << this->v_node.size() << std::endl;
    std::cout << this->v_face.size() << std::endl;
    std::cout << this->v_elem.size() << std::endl;
    std::cout << this->nFaceNodeNum << std::endl;

	// Read node X
	for (int i = 0; i < (int) this->v_node.size(); i++)
	{
		inFile >> this->v_node[i].x_;
	}

	// Read node Y
	for (int i = 0; i < (int) this->v_node.size(); i++)
	{
		inFile >> this->v_node[i].y_;
	}

	// Read node Z
	for (int i = 0; i < (int) this->v_node.size(); i++)
	{
		inFile >> this->v_node[i].z_;
	}

	// Read face node number
	getline(inFile,stTemp);
	getline(inFile, stTemp);
    if (stTemp == std::string("# node count per face"))
	{
		std::cout << LogLevel(LogLevel::mlOK, "Reading \"" + stTemp + "\"") << std::endl;
	}
	else
	{
        FatalError(std::string("Node number and entities are not the same. Please check!"));
	}

	int nTemp(0);
	for (int i = 0; i < (int) this->v_face.size(); i++)
	{
		inFile >> nTemp;
		this->v_face[i].v_nodeID.resize(nTemp);
	}

	// Read face nodes and check face nodes number
	getline(inFile, stTemp);
	getline(inFile, stTemp);

    if (stTemp == std::string("# face nodes"))
	{
		std::cout << LogLevel(LogLevel::mlOK, "Reading \"" + stTemp + "\"") << std::endl;
	}
	else
	{
        FatalError(std::string("Face number and entities are not the same. Please check!"));
	}

	int nFaceNodeCheck(0);
	for (int i = 0; i < (int)this->v_face.size(); i++)
	{
		for (int j = 0; j < (int)this->v_face[i].v_nodeID.size(); j++)
		{
			inFile >> this->v_face[i].v_nodeID[j];
			nFaceNodeCheck++;
		}
	}

	if (nFaceNodeCheck != this->nFaceNodeNum)
	{
        FatalError(std::string("Face nodes number and entities are not the same. Please check!"));
	}

	// Read left elements(Owner)
	getline(inFile, stTemp);
	getline(inFile, stTemp);
    if (stTemp == std::string("# left elements"))
	{
		std::cout << LogLevel(LogLevel::mlOK, "Reading \"" + stTemp + "\"") << std::endl;
	}
	else
	{
        FatalError(std::string("Face nodes number and entities are not the same. Please check!"));
	}

	for (int i = 0; i < (int) this->v_face.size(); i++)
	{
		inFile >> this->v_face[i].n_owner;
	}

	// Read right elements(Neighbor)
	getline(inFile, stTemp);
	getline(inFile, stTemp);
    if (stTemp == std::string("# right elements"))
	{
		std::cout << LogLevel(LogLevel::mlOK, "Reading \"" + stTemp + "\"") << std::endl;
	}
	else
	{
        FatalError(std::string("Left elements and entities are not the same. Please check!"));
	}

	for (int i = 0; i < (int) this->v_face.size(); i++)
	{
		inFile >> this->v_face[i].n_neighbor;
	}
}

void StandardMeshBlock::WriteTecplotMesh(const std::string& outMshFileName)
{
    std::ofstream	outFile(outMshFileName.c_str());

	outFile << "TITLE =\"" << "UnGrid Output" << "\"" << std::endl;

	if (this->md_meshDim == Mesh::md2D)
	{
		outFile << "VARIABLES = " << "\"x [m]\"," << "\"y [m]\"" << std::endl;

		if (this->est_shapeType == Element::estTriangular)
		{
			outFile << "ZONE T = " << "\"" << "Tri_mesh" << "\"," << " DATAPACKING = POINT, N = " << this->n_nodeNum << ", E = " << this->n_elemNum
				<< ", ZONETYPE = FETRIANGLE" << std::endl;
		}
		else if (this->est_shapeType == Element::estQuadrilateral)
		{
			outFile << "ZONE T = " << "\"" << "Quad_mesh" << "\"," << " DATAPACKING = POINT, N = " << this->n_nodeNum << ", E = " << this->n_elemNum
				<< ", ZONETYPE = FEQUADRILATERAL" << std::endl;
		}
		else if (this->est_shapeType == Element::estMixed)
		{
			outFile << "ZONE T = " << "\"" << "Mix_mesh" << "\"," << " DATAPACKING = POINT, N = " << this->n_nodeNum << ", E = " << this->n_elemNum
				<< ", ZONETYPE = FEQUADRILATERAL" << std::endl;
		}

		for (int i = 0; i < (int)this->n_nodeNum; i++)
		{
			outFile << this->v_node[i].x_ << " " << this->v_node[i].y_ << std::endl;
		}

		for (int i = 0; i < (int)this->n_elemNum; i++)
		{
			if (this->est_shapeType == Element::estTriangular)
			{
				outFile << this->v_elem[i].v_nodeID[0] + 1 << " " << this->v_elem[i].v_nodeID[1] + 1 << " " << this->v_elem[i].v_nodeID[2] + 1 << std::endl;
			}
			else if (this->est_shapeType == Element::estQuadrilateral)
			{
				outFile << this->v_elem[i].v_nodeID[0] + 1 << " " << this->v_elem[i].v_nodeID[1] + 1 << " " << this->v_elem[i].v_nodeID[2] + 1 << " " << this->v_elem[i].v_nodeID[3] + 1 << std::endl;
			}
			else if (this->est_shapeType == Element::estMixed)
			{
				if (this->v_elem[i].est_shapeType == Element::estTriangular)
				{
					outFile << this->v_elem[i].v_nodeID[0] + 1 << " " << this->v_elem[i].v_nodeID[1] + 1 << " " << this->v_elem[i].v_nodeID[2] + 1 << " " << this->v_elem[i].v_nodeID[2] + 1 << std::endl;
				}
				else if (this->v_elem[i].est_shapeType == Element::estQuadrilateral)
				{
					outFile << this->v_elem[i].v_nodeID[0] + 1 << " " << this->v_elem[i].v_nodeID[1] + 1 << " " << this->v_elem[i].v_nodeID[2] + 1 << " " << this->v_elem[i].v_nodeID[3] + 1 << std::endl;
				}
			}
		}

		// Write BC
		for (int i = 0; i < (int)this->v_boundaryFaceZone.size(); i++)
		{
			outFile << "ZONE T = " << "\"" << this->v_boundaryFaceZone[i].name << "\"," << " DATAPACKING = POINT, N = " << this->n_nodeNum << ", E = "
				<< this->v_boundaryFaceZone[i].v_faceID.size() << ", ZONETYPE = FELINESEG, " << "VARSHARELIST = ([1,2]=1,[1,2])" << std::endl;

			for (int j = 0; j < (int)this->v_boundaryFaceZone[i].v_faceID.size(); j++)
			{
				int nFaceID = this->v_boundaryFaceZone[i].v_faceID[j];

				for (int n = 0; n < (int)this->v_face[nFaceID].v_nodeID.size(); n++)
				{
					int nNodeID = this->v_face[nFaceID].v_nodeID[n];
					outFile << this->v_face[nFaceID].v_nodeID[n] + 1 << "\t";
				}

				outFile << std::endl;
			}
		}
	}
	else if (this->md_meshDim == Mesh::md3D)
	{
		if (this->est_shapeType == Element::estPolyhedron)
		{
			int nFaceNodeNum(0);
			for (int i = 0; i < (int)this->v_face.size(); i++)
			{
				for (int j = 0; j < (int)this->v_face[i].v_nodeID.size(); j++)
				{
					nFaceNodeNum++;
				}
			}

			outFile << "VARIABLES = " << std::endl << "\"X\"" << std::endl << "\"Y\"" << std::endl << "\"Z\"" << std::endl;
			outFile << "DATASETAUXDATA Common.VectorVarsAreVelocity = \"TRUE\"" << std::endl;
			outFile << "ZONE T = \"fluid\"" << std::endl;
			outFile << "STRANDID=0, SOLUTIONTIME=0" << std::endl;
			outFile << "Nodes=" << this->n_nodeNum << ", Faces=" << this->n_faceNum << ", Elements=" << this->n_elemNum << ", ZONETYPE=FEPolyhedron" << std::endl;
			outFile << "DATAPACKING=BLOCK" << std::endl;
			outFile << "TotalNumFaceNodes=" << nFaceNodeNum << ", NumConnectedBoundaryFaces=0, TotalNumBoundaryConnections=0" << std::endl;
			outFile << "DT=(SINGLE SINGLE SINGLE )";

			// Write node X
			for (int i = 0; i < (int) this->v_node.size(); i++)
			{
				if (i % 5 == 0) outFile << std::endl;
                outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << this->v_node[i].x_ << " ";
			}

			// Write node Y
			for (int i = 0; i < (int) this->v_node.size(); i++)
			{
				if (i % 5 == 0) outFile << std::endl;
                outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << this->v_node[i].y_ << " ";
			}

			// Write node Z
			for (int i = 0; i < (int) this->v_node.size(); i++)
			{
				if (i % 5 == 0) outFile << std::endl;
                outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << this->v_node[i].z_ << " ";
			}

			outFile << std::endl << "# node count per face";
			for (int i = 0; i < (int) this->v_face.size(); i++)
			{
				if (i % 10 == 0) outFile << std::endl;
				outFile << this->v_face[i].v_nodeID.size() << " ";
			}

			int nFaceNodeCheck(0);
			outFile << std::endl << "# face nodes";
			for (int i = 0; i < (int) this->v_face.size(); i++)
			{
				for (int j = 0; j < (int)this->v_face[i].v_nodeID.size(); j++)
				{
					if (nFaceNodeCheck % 10 == 0) outFile << std::endl;
					outFile << this->v_face[i].v_nodeID[j] + 1 << " ";
					nFaceNodeCheck++;
				}
			}

			outFile << std::endl << "# left elements";
			for (int i = 0; i < (int) this->v_face.size(); i++)
			{
				if (i % 10 == 0) outFile << std::endl;
				outFile << this->v_face[i].n_owner + 1 << " ";
			}

			outFile << std::endl << "# node count per face";
			for (int i = 0; i < (int) this->v_face.size(); i++)
			{
				if (i % 10 == 0) outFile << std::endl;
				outFile << this->v_face[i].n_neighbor + 1 << " ";
			}

			//// Write BC
			//for (int i = 0; i < (int)this->v_boundaryFaceZone.size(); i++)
			//{
			//	outFile << "ZONE T = " << "\"" << this->v_boundaryFaceZone[i].name << std::endl;
			//	outFile << "PARENTZONE=1" << std::endl;
			//	outFile << "STRANDID=0, SOLUTIONTIME=0" << std::endl;
			//	outFile << "Nodes=" << this->n_nodeNum << ", Faces=" << this->n_faceNum << ", Elements=" << this->n_elemNum << ", ZONETYPE=FEPolyhedron" << std::endl;
			//	outFile << "DATAPACKING=BLOCK" << std::endl;
			//	outFile << "TotalNumFaceNodes=" << nFaceNodeNum << ", NumConnectedBoundaryFaces=0, TotalNumBoundaryConnections=0" << std::endl;
			//	outFile << "DT=(SINGLE SINGLE SINGLE )";

			//	<< "\"," << " DATAPACKING = POINT, N = " << this->n_nodeNum << ", E = "
			//		<< this->v_boundaryFaceZone[i].v_faceID.size() << ", ZONETYPE = FEQUADRILATERAL, " << "VARSHARELIST = ([1,2,3]=1,[1,2,3])" << std::endl;

			//	for (int j = 0; j < (int)this->v_boundaryFaceZone[i].v_faceID.size(); j++)
			//	{
			//		int nFaceID = this->v_boundaryFaceZone[i].v_faceID[j];

			//		for (int n = 0; n < (int)this->v_face[nFaceID].v_nodeID.size(); n++)
			//		{
			//			int nNodeID = this->v_face[nFaceID].v_nodeID[n];
			//			outFile << this->v_face[nFaceID].v_nodeID[n] + 1 << "\t";
			//		}

			//		if (this->v_face[nFaceID].ft_faceType == Face::ftTrangular)
			//		{
			//			// ZONETYPE = FEQUADRILATERAL need four Node to display trangle, here output the last again
			//			outFile << this->v_face[nFaceID].v_nodeID[(int)this->v_face[nFaceID].v_nodeID.size() - 1] + 1 << "\t";
			//		}
			//		outFile << std::endl;
			//	}
			//}
		}
		else
		{
			outFile << "VARIABLES = " << "\"x [m]\"," << "\"y [m]\"," << "\"z [m]\"" << std::endl;

			if (this->est_shapeType == Element::estTetrahedral)
			{
				outFile << "ZONE T = " << "\"" << "Tetrahedral_mesh" << "\"," << " DATAPACKING = POINT, N = " << this->n_nodeNum << ", E = " << this->n_elemNum
					<< ", ZONETYPE = FETETRAHEDRON" << std::endl;
			}
			else if (this->est_shapeType == Element::estHexahedral)
			{
				outFile << "ZONE T = " << "\"" << "Hexahedral_mesh" << "\"," << " DATAPACKING = POINT, N = " << this->n_nodeNum << ", E = " << this->n_elemNum
					<< ", ZONETYPE = FEBRICK" << std::endl;
			}
			else if (this->est_shapeType == Element::estMixed)
			{
				outFile << "ZONE T = " << "\"" << "Tet_mesh" << "\"," << " DATAPACKING = POINT, N = " << this->n_nodeNum << ", E = " << this->n_elemNum
					<< ", ZONETYPE = FEBRICK" << std::endl;
			}
			else if (this->est_shapeType == Element::estWedge)
			{
				outFile << "ZONE T = " << "\"" << "Wedge_mesh" << "\"," << " DATAPACKING = POINT, N = " << this->n_nodeNum << ", E = " << this->n_elemNum
					<< ", ZONETYPE = FEBRICK" << std::endl;
			}

			for (int i = 0; i < this->n_nodeNum; i++)
			{
                outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << this->v_node[i].x_ << " " << this->v_node[i].y_ << " " << this->v_node[i].z_ << std::endl;
			}

			if (this->est_shapeType == Element::estTetrahedral)
			{
				for (int i = 0; i < (int)this->v_elem.size(); i++)
				{
					outFile << this->v_elem[i].v_nodeID[0] + 1 << " " << this->v_elem[i].v_nodeID[1] + 1 << " " << this->v_elem[i].v_nodeID[2] + 1 << " " << this->v_elem[i].v_nodeID[3] + 1 << std::endl;
				}
			}
			else if (this->est_shapeType == Element::estHexahedral)
			{
				for (int i = 0; i < (int)this->v_elem.size(); i++)
				{
					outFile << this->v_elem[i].v_nodeID[0] + 1 << " " << this->v_elem[i].v_nodeID[1] + 1 << " " << this->v_elem[i].v_nodeID[2] + 1 << " " << this->v_elem[i].v_nodeID[3] + 1 << " "
						<< this->v_elem[i].v_nodeID[4] + 1 << " " << this->v_elem[i].v_nodeID[5] + 1 << " " << this->v_elem[i].v_nodeID[6] + 1 << " " << this->v_elem[i].v_nodeID[7] + 1 << std::endl;
				}
			}
			else if (this->est_shapeType == Element::estMixed)
			{
				for (int i = 0; i < (int)this->v_elem.size(); i++)
				{
					if (this->v_elem[i].est_shapeType == Element::estTetrahedral)
					{
						outFile << this->v_elem[i].v_nodeID[0] + 1 << " " << this->v_elem[i].v_nodeID[1] + 1 << " " << this->v_elem[i].v_nodeID[2] + 1 << " " << this->v_elem[i].v_nodeID[2] + 1 << " "
							<< this->v_elem[i].v_nodeID[3] + 1 << " " << this->v_elem[i].v_nodeID[3] + 1 << " " << this->v_elem[i].v_nodeID[3] + 1 << " " << this->v_elem[i].v_nodeID[3] + 1 << std::endl;
					}
					else if (this->v_elem[i].est_shapeType == Element::estPyramid)
					{
						outFile << this->v_elem[i].v_nodeID[0] + 1 << " " << this->v_elem[i].v_nodeID[1] + 1 << " " << this->v_elem[i].v_nodeID[2] + 1 << " " << this->v_elem[i].v_nodeID[3] + 1 << " "
							<< this->v_elem[i].v_nodeID[4] + 1 << " " << this->v_elem[i].v_nodeID[4] + 1 << " " << this->v_elem[i].v_nodeID[4] + 1 << " " << this->v_elem[i].v_nodeID[4] + 1 << std::endl;
					}
					else if (this->v_elem[i].est_shapeType == Element::estWedge)
					{
						outFile << this->v_elem[i].v_nodeID[0] + 1 << " " << this->v_elem[i].v_nodeID[1] + 1 << " " << this->v_elem[i].v_nodeID[4] + 1 << " " << this->v_elem[i].v_nodeID[3] + 1 << " "
							<< this->v_elem[i].v_nodeID[2] + 1 << " " << this->v_elem[i].v_nodeID[2] + 1 << " " << this->v_elem[i].v_nodeID[5] + 1 << " " << this->v_elem[i].v_nodeID[5] + 1 << std::endl;
					}
					else if (this->v_elem[i].est_shapeType == Element::estHexahedral)
					{
						outFile << this->v_elem[i].v_nodeID[0] + 1 << " " << this->v_elem[i].v_nodeID[1] + 1 << " " << this->v_elem[i].v_nodeID[2] + 1 << " " << this->v_elem[i].v_nodeID[3] + 1 << " "
							<< this->v_elem[i].v_nodeID[4] + 1 << " " << this->v_elem[i].v_nodeID[5] + 1 << " " << this->v_elem[i].v_nodeID[6] + 1 << " " << this->v_elem[i].v_nodeID[7] + 1 << std::endl;
					}
				}
			}
			else if (this->est_shapeType == Element::estWedge)
			{
				for (int i = 0; i < (int)this->v_elem.size(); i++)
				{
					outFile << this->v_elem[i].v_nodeID[0] + 1 << " " << this->v_elem[i].v_nodeID[1] + 1 << " " << this->v_elem[i].v_nodeID[5] + 1 << " " << this->v_elem[i].v_nodeID[3] + 1 << " "
						<< this->v_elem[i].v_nodeID[2] + 1 << " " << this->v_elem[i].v_nodeID[2] + 1 << " " << this->v_elem[i].v_nodeID[4] + 1 << " " << this->v_elem[i].v_nodeID[4] + 1 << std::endl;
				}
			}

			// Write BC
			for (int i = 0; i < (int)this->v_boundaryFaceZone.size(); i++)
			{
				outFile << "ZONE T = " << "\"" << this->v_boundaryFaceZone[i].name << "\"," << " DATAPACKING = POINT, N = " << this->n_nodeNum << ", E = "
					<< this->v_boundaryFaceZone[i].v_faceID.size() << ", ZONETYPE = FEQUADRILATERAL, " << "VARSHARELIST = ([1,2,3]=1,[1,2,3])" << std::endl;

				for (int j = 0; j < (int)this->v_boundaryFaceZone[i].v_faceID.size(); j++)
				{
					int nFaceID = this->v_boundaryFaceZone[i].v_faceID[j];

					for (int n = 0; n < (int)this->v_face[nFaceID].v_nodeID.size(); n++)
					{
						int nNodeID = this->v_face[nFaceID].v_nodeID[n];
						outFile << this->v_face[nFaceID].v_nodeID[n] + 1 << "\t";
					}

					if (this->v_face[nFaceID].ft_faceType == Face::ftTrangular)
					{
						// ZONETYPE = FEQUADRILATERAL need four Node to display trangle, here output the last again
						outFile << this->v_face[nFaceID].v_nodeID[(int)this->v_face[nFaceID].v_nodeID.size() - 1] + 1 << "\t";
					}
					outFile << std::endl;
				}
			}
		}
	}

	outFile.close();

	//system(outMshFileName.c_str());
}

void StandardMeshBlock::WriteTecplotMeshOneElem(const std::string& outMshFileName, int elemID)
{
    std::ofstream	outFile(outMshFileName.c_str());

	outFile << "TITLE =\"" << "UnGrid Output" << "\"" << std::endl;

	if (this->md_meshDim == Mesh::md2D)
	{
		outFile << "VARIABLES = " << "\"x [m]\"," << "\"y [m]\"" << std::endl;

		if (this->est_shapeType == Element::estTriangular)
		{
			outFile << "ZONE T = " << "\"" << "Tri_mesh" << "\"," << " DATAPACKING = POINT, N = " << this->n_nodeNum << ", E = " << "1"
				<< ", ZONETYPE = FETRIANGLE" << std::endl;
		}
		else if (this->est_shapeType == Element::estQuadrilateral)
		{
			outFile << "ZONE T = " << "\"" << "Quad_mesh" << "\"," << " DATAPACKING = POINT, N = " << this->n_nodeNum << ", E = " << "1"
				<< ", ZONETYPE = FEQUADRILATERAL" << std::endl;
		}
		else if (this->est_shapeType == Element::estMixed)
		{
			outFile << "ZONE T = " << "\"" << "Mix_mesh" << "\"," << " DATAPACKING = POINT, N = " << this->n_nodeNum << ", E = " << "1"
				<< ", ZONETYPE = FEQUADRILATERAL" << std::endl;
		}

		for (int i = 0; i < (int)this->n_nodeNum; i++)
		{
			outFile << this->v_node[i].x_ << " " << this->v_node[i].y_ << std::endl;
		}

		if (this->est_shapeType == Element::estTriangular)
		{
			outFile << this->v_elem[elemID].v_nodeID[0] + 1 << " " << this->v_elem[elemID].v_nodeID[1] + 1 << " " << this->v_elem[elemID].v_nodeID[2] + 1 << std::endl;
		}
		else if (this->est_shapeType == Element::estQuadrilateral)
		{
			outFile << this->v_elem[elemID].v_nodeID[0] + 1 << " " << this->v_elem[elemID].v_nodeID[1] + 1 << " " << this->v_elem[elemID].v_nodeID[2] + 1 << " " << this->v_elem[elemID].v_nodeID[3] + 1 << std::endl;
		}
		else if (this->est_shapeType == Element::estMixed)
		{
			if (this->v_elem[elemID].est_shapeType == Element::estTriangular)
			{
				outFile << this->v_elem[elemID].v_nodeID[0] + 1 << " " << this->v_elem[elemID].v_nodeID[1] + 1 << " " << this->v_elem[elemID].v_nodeID[2] + 1 << " " << this->v_elem[elemID].v_nodeID[2] + 1 << std::endl;
			}
			else if (this->v_elem[elemID].est_shapeType == Element::estQuadrilateral)
			{
				outFile << this->v_elem[elemID].v_nodeID[0] + 1 << " " << this->v_elem[elemID].v_nodeID[1] + 1 << " " << this->v_elem[elemID].v_nodeID[2] + 1 << " " << this->v_elem[elemID].v_nodeID[3] + 1 << std::endl;
			}
		}
	}
	else if (this->md_meshDim == Mesh::md3D)
	{
		outFile << "VARIABLES = " << "\"x [m]\"," << "\"y [m]\"," << "\"z [m]\"" << std::endl;

		if (this->est_shapeType == Element::estTetrahedral)
		{
			outFile << "ZONE T = " << "\"" << "Tetrahedral_mesh" << "\"," << " DATAPACKING = POINT, N = " << this->n_nodeNum << ", E = " << "1"
				<< ", ZONETYPE = FETETRAHEDRON" << std::endl;
		}
		else if (this->est_shapeType == Element::estHexahedral)
		{
			outFile << "ZONE T = " << "\"" << "Hexahedral_mesh" << "\"," << " DATAPACKING = POINT, N = " << this->n_nodeNum << ", E = " << "1"
				<< ", ZONETYPE = FEBRICK" << std::endl;
		}
		else if (this->est_shapeType == Element::estMixed)
		{
			outFile << "ZONE T = " << "\"" << "Tet_mesh" << "\"," << " DATAPACKING = POINT, N = " << this->n_nodeNum << ", E = " << "1"
				<< ", ZONETYPE = FEBRICK" << std::endl;
		}
		else if (this->est_shapeType == Element::estWedge)
		{
			outFile << "ZONE T = " << "\"" << "Wedge_mesh" << "\"," << " DATAPACKING = POINT, N = " << this->n_nodeNum << ", E = " << "1"
				<< ", ZONETYPE = FEBRICK" << std::endl;
		}

		for (int i = 0; i < this->n_nodeNum; i++)
		{
            outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << this->v_node[i].x_ << " " << this->v_node[i].y_ << " " << this->v_node[i].z_ << std::endl;
		}

		if (this->est_shapeType == Element::estTetrahedral)
		{
			outFile << this->v_elem[elemID].v_nodeID[0] + 1 << " " << this->v_elem[elemID].v_nodeID[1] + 1 << " " << this->v_elem[elemID].v_nodeID[2] + 1 << " " << this->v_elem[elemID].v_nodeID[3] + 1 << std::endl;
		}
		else if (this->est_shapeType == Element::estHexahedral)
		{
			outFile << this->v_elem[elemID].v_nodeID[0] + 1 << " " << this->v_elem[elemID].v_nodeID[1] + 1 << " " << this->v_elem[elemID].v_nodeID[2] + 1 << " " << this->v_elem[elemID].v_nodeID[3] + 1 << " "
				<< this->v_elem[elemID].v_nodeID[4] + 1 << " " << this->v_elem[elemID].v_nodeID[5] + 1 << " " << this->v_elem[elemID].v_nodeID[6] + 1 << " " << this->v_elem[elemID].v_nodeID[7] + 1 << std::endl;
		}
		else if (this->est_shapeType == Element::estMixed)
		{
			if (this->v_elem[elemID].est_shapeType == Element::estTetrahedral)
			{
				outFile << this->v_elem[elemID].v_nodeID[0] + 1 << " " << this->v_elem[elemID].v_nodeID[1] + 1 << " " << this->v_elem[elemID].v_nodeID[2] + 1 << " " << this->v_elem[elemID].v_nodeID[2] + 1 << " "
					<< this->v_elem[elemID].v_nodeID[3] + 1 << " " << this->v_elem[elemID].v_nodeID[3] + 1 << " " << this->v_elem[elemID].v_nodeID[3] + 1 << " " << this->v_elem[elemID].v_nodeID[3] + 1 << std::endl;
			}
			else if (this->v_elem[elemID].est_shapeType == Element::estPyramid)
			{
				outFile << this->v_elem[elemID].v_nodeID[0] + 1 << " " << this->v_elem[elemID].v_nodeID[1] + 1 << " " << this->v_elem[elemID].v_nodeID[2] + 1 << " " << this->v_elem[elemID].v_nodeID[3] + 1 << " "
					<< this->v_elem[elemID].v_nodeID[4] + 1 << " " << this->v_elem[elemID].v_nodeID[4] + 1 << " " << this->v_elem[elemID].v_nodeID[4] + 1 << " " << this->v_elem[elemID].v_nodeID[4] + 1 << std::endl;
			}
			else if (this->v_elem[elemID].est_shapeType == Element::estWedge)
			{
				outFile << this->v_elem[elemID].v_nodeID[0] + 1 << " " << this->v_elem[elemID].v_nodeID[1] + 1 << " " << this->v_elem[elemID].v_nodeID[4] + 1 << " " << this->v_elem[elemID].v_nodeID[3] + 1 << " "
					<< this->v_elem[elemID].v_nodeID[2] + 1 << " " << this->v_elem[elemID].v_nodeID[2] + 1 << " " << this->v_elem[elemID].v_nodeID[5] + 1 << " " << this->v_elem[elemID].v_nodeID[5] + 1 << std::endl;
			}
			else if (this->v_elem[elemID].est_shapeType == Element::estHexahedral)
			{
				outFile << this->v_elem[elemID].v_nodeID[0] + 1 << " " << this->v_elem[elemID].v_nodeID[1] + 1 << " " << this->v_elem[elemID].v_nodeID[2] + 1 << " " << this->v_elem[elemID].v_nodeID[3] + 1 << " "
					<< this->v_elem[elemID].v_nodeID[4] + 1 << " " << this->v_elem[elemID].v_nodeID[5] + 1 << " " << this->v_elem[elemID].v_nodeID[6] + 1 << " " << this->v_elem[elemID].v_nodeID[7] + 1 << std::endl;
			}
		}
		else if (this->est_shapeType == Element::estWedge)
		{
			outFile << this->v_elem[elemID].v_nodeID[0] + 1 << " " << this->v_elem[elemID].v_nodeID[1] + 1 << " " << this->v_elem[elemID].v_nodeID[5] + 1 << " " << this->v_elem[elemID].v_nodeID[3] + 1 << " "
				<< this->v_elem[elemID].v_nodeID[2] + 1 << " " << this->v_elem[elemID].v_nodeID[2] + 1 << " " << this->v_elem[elemID].v_nodeID[4] + 1 << " " << this->v_elem[elemID].v_nodeID[4] + 1 << std::endl;
		}	
	}

	outFile.close();
}
