/*---------------------------------------------------------------------------*\
File Name:
	FieldIO.cpp

Description:
	1. Pack Field and make IO for Field;

	Author:		Shuai Zhang
	Date: 		2017-12-11

	
Revised:
	Description:
	1. Convert "non const" input parameters to const;
	2. Convert "non const" function to const function;
	
	Revisor:		Shuai Zhang
	Modified Date:	2018-11-28
\*---------------------------------------------------------------------------*/

#include "../MHT_IO/FieldIO.h"
#include <map>
#include <iomanip>
#include "../MHT_mesh/Mesh.h"
#include "../MHT_field/Field.h"

FieldIO::FieldIO()
	:
	p_blockMesh(NULL)
{
}

void FieldIO::push_backScalarField(Field<Scalar>& ScalarField)
{
	if (p_blockMesh == NULL)
	{
		p_blockMesh = ScalarField.p_blockMesh;
	}

	v_scalarField.push_back(&ScalarField);
}

void FieldIO::push_backVectorField(Field<Vector>& VectorField)
{
	if (p_blockMesh == NULL)
	{
		p_blockMesh = VectorField.p_blockMesh;
	}

	v_vectorField.push_back(&VectorField);
}

void FieldIO::WriteTecplotField(const std::string& outMshFileName)
{
	this->FieldElementToNode();
    std::ofstream	outFile(outMshFileName.c_str());

    outFile << "TITLE =\"" << "NHT-CFD-Solver" << "\"" << std::endl;

	if (p_blockMesh->md_meshDim == Mesh::md2D)
	{
		int nValNum = 2 + (int)v_scalarField.size() + (int)v_vectorField.size() * 2;

		outFile << "VARIABLES = " << "\"x [m]\"," << "\"y [m]\",";
		if ((int)v_scalarField.size() > 0)
		{
			for (int i = 0; i < (int)v_scalarField.size() - 1; i++)
			{
				Field<Scalar>& CurField = *v_scalarField[i];
				outFile << "\"" << CurField.st_name << "\",";
			}
			outFile << "\"" << v_scalarField[(int)v_scalarField.size() - 1]->st_name << "\"";
		}

		if ((int)v_vectorField.size() > 0)
		{
			outFile << ",";
			for (int i = 0; i < (int)v_vectorField.size() - 1; i++)
			{
				Field<Vector>& CurField = *v_vectorField[i];
				outFile << "\"" << CurField.st_name << ".x\"," << "\"" << CurField.st_name << ".y\",";
			}
			outFile << "\"" << v_vectorField[(int)v_vectorField.size() - 1]->st_name << ".x\"," << "\"" << v_vectorField[(int)v_vectorField.size() - 1]->st_name << ".y\"";
		}
		outFile << std::endl;

		if (p_blockMesh->est_shapeType == Element::estTriangular)
		{
			outFile << "ZONE T = " << "\"" << "Tri_mesh" << "\"," << " DATAPACKING = POINT, N = " << p_blockMesh->n_nodeNum << ", E = " << p_blockMesh->n_elemNum
                << ", ZONETYPE = FETRIANGLE" << std::endl;
		}
		else if (p_blockMesh->est_shapeType == Element::estQuadrilateral)
		{
			outFile << "ZONE T = " << "\"" << "Quad_mesh" << "\"," << " DATAPACKING = POINT, N = " << p_blockMesh->n_nodeNum << ", E = " << p_blockMesh->n_elemNum
                << ", ZONETYPE = FEQUADRILATERAL" << std::endl;
		}
		else if (p_blockMesh->est_shapeType == Element::estMixed)
		{
			outFile << "ZONE T = " << "\"" << "Mix_mesh" << "\"," << " DATAPACKING = POINT, N = " << p_blockMesh->n_nodeNum << ", E = " << p_blockMesh->n_elemNum
                << ", ZONETYPE = FEQUADRILATERAL" << std::endl;
		}

		for (int i = 0; i < (int)p_blockMesh->n_nodeNum; i++)
		{
            outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << p_blockMesh->v_node[i].x_ << " " << p_blockMesh->v_node[i].y_ << " ";
			for (int j = 0; j < (int)v_scalarField.size() - 1; j++)
			{
				Field<Scalar>& CurField = *v_scalarField[j];
				outFile << CurField.nodeField.GetValue(i) << " ";
			}
			if ((int)v_scalarField.size() > 0)
			{
				outFile << (*v_scalarField[(int)v_scalarField.size() - 1]).nodeField.GetValue(i) << " ";
			}

			for (int j = 0; j < (int)v_vectorField.size() - 1; j++)
			{
				Field<Vector>& CurField = *v_vectorField[j];
				outFile << CurField.nodeField.GetValue(i).x_ << " " << CurField.nodeField.GetValue(i).y_;
			}
			if ((int)v_vectorField.size() > 0)
			{
				outFile << (*v_vectorField[(int)v_vectorField.size() - 1]).nodeField.GetValue(i).x_ << " " << (*v_vectorField[(int)v_vectorField.size() - 1]).nodeField.GetValue(i).y_;
			}
			outFile << std::endl;
		}

		for (int i = 0; i < (int)p_blockMesh->n_elemNum; i++)
		{
			if (p_blockMesh->est_shapeType == Element::estTriangular)
			{
                outFile << p_blockMesh->v_elem[i].v_nodeID[0] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[1] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[2] + 1 << std::endl;
			}
			else if (p_blockMesh->est_shapeType == Element::estQuadrilateral)
			{
                outFile << p_blockMesh->v_elem[i].v_nodeID[0] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[1] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[2] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[3] + 1 << std::endl;
			}
			else if (p_blockMesh->est_shapeType == Element::estMixed)
			{
				if (p_blockMesh->v_elem[i].est_shapeType == Element::estTriangular)
				{
                    outFile << p_blockMesh->v_elem[i].v_nodeID[0] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[1] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[2] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[2] + 1 << std::endl;
				}
				else if (p_blockMesh->v_elem[i].est_shapeType == Element::estQuadrilateral)
				{
                    outFile << p_blockMesh->v_elem[i].v_nodeID[0] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[1] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[2] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[3] + 1 << std::endl;
				}
			}
		}

		// Write BC 
		for (int i = 0; i < (int)this->p_blockMesh->v_boundaryFaceZone.size(); i++)
		{
			outFile << "ZONE T = " << "\"" << this->p_blockMesh->v_boundaryFaceZone[i].name << "\"," << " DATAPACKING = POINT, N = " << this->p_blockMesh->n_nodeNum << ", E = "
				<< this->p_blockMesh->v_boundaryFaceZone[i].v_faceID.size() << ", ZONETYPE = FELINESEG, " << "VARSHARELIST = ([";

			for (int j = 0; j < nValNum - 1; j++)
			{
				outFile << j + 1 << ",";
			}

			outFile << nValNum << "]=1,[1,2])" << std::endl;

			for (int j = 0; j < (int)this->p_blockMesh->v_boundaryFaceZone[i].v_faceID.size(); j++)
			{
				int nFaceID = this->p_blockMesh->v_boundaryFaceZone[i].v_faceID[j];

				for (int n = 0; n < (int)this->p_blockMesh->v_face[nFaceID].v_nodeID.size(); n++)
				{
					outFile << this->p_blockMesh->v_face[nFaceID].v_nodeID[n] + 1 << "\t";
				}

				outFile << std::endl;
			}
		}
	}
	else if (this->p_blockMesh->md_meshDim == Mesh::md3D)
	{
		if (this->p_blockMesh->est_shapeType == Element::estPolyhedron)
		{
			//Not implimented
		}
		else
		{
			int nValNum = 3 + (int)v_scalarField.size() + (int)v_vectorField.size() * 3;

			outFile << "VARIABLES = " << "\"x [m]\"," << "\"y [m]\"," << "\"z [m]\",";
			if ((int)v_scalarField.size() > 0)
			{
				for (int i = 0; i < (int)v_scalarField.size() - 1; i++)
				{
					Field<Scalar>& CurField = *v_scalarField[i];
					outFile << "\"" << CurField.st_name << "\",";
				}
				outFile << "\"" << v_scalarField[(int)v_scalarField.size() - 1]->st_name << "\"";
			}

			if ((int)v_vectorField.size() > 0)
			{
				outFile << ",";
				for (int i = 0; i < (int)v_vectorField.size() - 1; i++)
				{
					Field<Vector>& CurField = *v_vectorField[i];
					outFile << "\"" << CurField.st_name << ".x\"," << "\"" << CurField.st_name << ".y\"," << "\"" << CurField.st_name << ".z\",";
				}
				outFile << "\"" << v_vectorField[(int)v_vectorField.size() - 1]->st_name << ".x\"," << "\"" << v_vectorField[(int)v_vectorField.size() - 1]->st_name << ".y\""
					<< "\"" << v_vectorField[(int)v_vectorField.size() - 1]->st_name << ".z\"";
			}
			outFile << std::endl;

			if (p_blockMesh->est_shapeType == Element::estTetrahedral)
			{
				outFile << "ZONE T = " << "\"" << "Tetrahedral_mesh" << "\"," << " DATAPACKING = POINT, N = " << p_blockMesh->n_nodeNum << ", E = " << p_blockMesh->n_elemNum
                    << ", ZONETYPE = FETETRAHEDRON" << std::endl;
			}
			else if (p_blockMesh->est_shapeType == Element::estHexahedral)
			{
				outFile << "ZONE T = " << "\"" << "Hexahedral_mesh" << "\"," << " DATAPACKING = POINT, N = " << p_blockMesh->n_nodeNum << ", E = " << p_blockMesh->n_elemNum
                    << ", ZONETYPE = FEBRICK" << std::endl;
			}
			else if (p_blockMesh->est_shapeType == Element::estMixed)
			{
				outFile << "ZONE T = " << "\"" << "Tet_mesh" << "\"," << " DATAPACKING = POINT, N = " << p_blockMesh->n_nodeNum << ", E = " << p_blockMesh->n_elemNum
                    << ", ZONETYPE = FEBRICK" << std::endl;
			}
			else if (p_blockMesh->est_shapeType == Element::estWedge)
			{
				outFile << "ZONE T = " << "\"" << "Wedge_mesh" << "\"," << " DATAPACKING = POINT, N = " << p_blockMesh->n_nodeNum << ", E = " << p_blockMesh->n_elemNum
                    << ", ZONETYPE = FEBRICK" << std::endl;
			}

			for (int i = 0; i < (int)p_blockMesh->n_nodeNum; i++)
			{
                outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << p_blockMesh->v_node[i].x_ << " " << p_blockMesh->v_node[i].y_ << " " << p_blockMesh->v_node[i].z_ << " ";
				for (int j = 0; j < (int)v_scalarField.size() - 1; j++)
				{
					Field<Scalar>& CurField = *v_scalarField[j];
					outFile << CurField.nodeField.GetValue(i) << " ";
				}
				if ((int)v_scalarField.size() > 0)
				{
					outFile << (*v_scalarField[(int)v_scalarField.size() - 1]).nodeField.GetValue(i) << " ";
				}

				for (int j = 0; j < (int)v_vectorField.size() - 1; j++)
				{
					Field<Vector>& CurField = *v_vectorField[j];
					outFile << CurField.nodeField.GetValue(i).x_ << " " << CurField.nodeField.GetValue(i).y_ << " " << CurField.nodeField.GetValue(i).z_;
				}
				if ((int)v_vectorField.size() > 0)
				{ 
					outFile << (*v_vectorField[(int)v_vectorField.size() - 1]).nodeField.GetValue(i).x_ << " " << (*v_vectorField[(int)v_vectorField.size() - 1]).nodeField.GetValue(i).y_
						<< " " << (*v_vectorField[(int)v_vectorField.size() - 1]).nodeField.GetValue(i).z_;
				}
				outFile << std::endl;
			}

			if (p_blockMesh->est_shapeType == Element::estTetrahedral)
			{
				for (int i = 0; i < (int)p_blockMesh->v_elem.size(); i++)
				{
                    outFile << p_blockMesh->v_elem[i].v_nodeID[0] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[1] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[2] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[3] + 1 << std::endl;
				}
			}
			else if (p_blockMesh->est_shapeType == Element::estHexahedral)
			{
				for (int i = 0; i < (int)p_blockMesh->v_elem.size(); i++)
				{
					outFile << p_blockMesh->v_elem[i].v_nodeID[0] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[1] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[2] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[3] + 1 << " "
                        << p_blockMesh->v_elem[i].v_nodeID[4] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[5] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[6] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[7] + 1 << std::endl;
				}
			}
			else if (p_blockMesh->est_shapeType == Element::estMixed)
			{
				for (int i = 0; i < (int)p_blockMesh->v_elem.size(); i++)
				{
					if (p_blockMesh->v_elem[i].est_shapeType == Element::estTetrahedral)
					{
						outFile << p_blockMesh->v_elem[i].v_nodeID[0] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[1] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[2] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[2] + 1 << " "
                            << p_blockMesh->v_elem[i].v_nodeID[3] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[3] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[3] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[3] + 1 << std::endl;
					}
					else if (p_blockMesh->v_elem[i].est_shapeType == Element::estPyramid)
					{
						outFile << p_blockMesh->v_elem[i].v_nodeID[0] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[1] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[2] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[3] + 1 << " "
                            << p_blockMesh->v_elem[i].v_nodeID[4] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[4] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[4] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[4] + 1 << std::endl;
					}
					else if (p_blockMesh->v_elem[i].est_shapeType == Element::estWedge)
					{
						outFile << p_blockMesh->v_elem[i].v_nodeID[0] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[1] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[4] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[3] + 1 << " "
                            << p_blockMesh->v_elem[i].v_nodeID[2] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[2] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[5] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[5] + 1 << std::endl;
					}
					else if (p_blockMesh->v_elem[i].est_shapeType == Element::estHexahedral)
					{
						outFile << p_blockMesh->v_elem[i].v_nodeID[0] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[1] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[2] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[3] + 1 << " "
                            << p_blockMesh->v_elem[i].v_nodeID[4] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[5] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[6] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[7] + 1 << std::endl;
					}
				}
			}
			else if (p_blockMesh->est_shapeType == Element::estWedge)
			{
				for (int i = 0; i < (int)p_blockMesh->v_elem.size(); i++)
				{
					outFile << p_blockMesh->v_elem[i].v_nodeID[0] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[1] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[5] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[3] + 1 << " "
                        << p_blockMesh->v_elem[i].v_nodeID[2] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[2] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[4] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[4] + 1 << std::endl;
				}
			}

			// Write BC
			for (int i = 0; i < (int)this->p_blockMesh->v_boundaryFaceZone.size(); i++)
			{
				outFile << "ZONE T = " << "\"" << this->p_blockMesh->v_boundaryFaceZone[i].name << "\"," << " DATAPACKING = POINT, N = " << this->p_blockMesh->n_nodeNum << ", E = "
					<< this->p_blockMesh->v_boundaryFaceZone[i].v_faceID.size() << ", ZONETYPE = FEQUADRILATERAL, " << "VARSHARELIST = ([";
				for (int j = 0; j < nValNum - 1; j++)
				{
					outFile << j + 1 << ",";
				}
				outFile << nValNum << "]=1,[1,2,3])" << std::endl;

				for (int j = 0; j < (int)this->p_blockMesh->v_boundaryFaceZone[i].v_faceID.size(); j++)
				{
					int nFaceID = this->p_blockMesh->v_boundaryFaceZone[i].v_faceID[j];

					for (int n = 0; n < (int)this->p_blockMesh->v_face[nFaceID].v_nodeID.size(); n++)
					{
						outFile << this->p_blockMesh->v_face[nFaceID].v_nodeID[n] + 1 << "\t";
					}

					if (this->p_blockMesh->v_face[nFaceID].ft_faceType == Face::ftTrangular)
					{
						// ZONETYPE = FEQUADRILATERAL need four Node to display trangle, here output the last again
						outFile << this->p_blockMesh->v_face[nFaceID].v_nodeID[(int)this->p_blockMesh->v_face[nFaceID].v_nodeID.size() - 1] + 1 << "\t";
					}
					outFile << std::endl;
				}
			}
		}
	}

	outFile.close();
}

void FieldIO::WriteTecplotField(const std::string& outMshFileName, Scalar Time)
{
	this->FieldElementToNode();
	std::ofstream	outFile(outMshFileName.c_str());

	outFile << "TITLE =\"" << "NHT-CFD-Solver" << "\"" << std::endl;

	if (p_blockMesh->md_meshDim == Mesh::md2D)
	{
		int nValNum = 2 + (int)v_scalarField.size() + (int)v_vectorField.size() * 2;

		outFile << "VARIABLES = " << "\"x [m]\"," << "\"y [m]\",";
		if ((int)v_scalarField.size() > 0)
		{
			for (int i = 0; i < (int)v_scalarField.size() - 1; i++)
			{
				Field<Scalar>& CurField = *v_scalarField[i];
				outFile << "\"" << CurField.st_name << "\",";
			}
			outFile << "\"" << v_scalarField[(int)v_scalarField.size() - 1]->st_name << "\"";
		}

		if ((int)v_vectorField.size() > 0)
		{
			outFile << ",";
			for (int i = 0; i < (int)v_vectorField.size() - 1; i++)
			{
				Field<Vector>& CurField = *v_vectorField[i];
				outFile << "\"" << CurField.st_name << ".x\"," << "\"" << CurField.st_name << ".y\",";
			}
			outFile << "\"" << v_vectorField[(int)v_vectorField.size() - 1]->st_name << ".x\"," << "\"" << v_vectorField[(int)v_vectorField.size() - 1]->st_name << ".y\"";
		}
		outFile << std::endl;

		if (p_blockMesh->est_shapeType == Element::estTriangular)
		{
			outFile << "ZONE T = " << "\"" << "Tri_mesh" << "\"," << "STRANDID = 1, SOLUTIONTIME = "<< Time << " DATAPACKING = POINT, N = " << p_blockMesh->n_nodeNum << ", E = " << p_blockMesh->n_elemNum
				<< ", ZONETYPE = FETRIANGLE" << std::endl;
		}
		else if (p_blockMesh->est_shapeType == Element::estQuadrilateral)
		{
			outFile << "ZONE T = " << "\"" << "Quad_mesh" << "\"," << "STRANDID = 1, SOLUTIONTIME = " << Time << " DATAPACKING = POINT, N = " << p_blockMesh->n_nodeNum << ", E = " << p_blockMesh->n_elemNum
				<< ", ZONETYPE = FEQUADRILATERAL" << std::endl;
		}
		else if (p_blockMesh->est_shapeType == Element::estMixed)
		{
			outFile << "ZONE T = " << "\"" << "Mix_mesh" << "\"," << "STRANDID = 1, SOLUTIONTIME = " << Time << " DATAPACKING = POINT, N = " << p_blockMesh->n_nodeNum << ", E = " << p_blockMesh->n_elemNum
				<< ", ZONETYPE = FEQUADRILATERAL" << std::endl;
		}

		for (int i = 0; i < (int)p_blockMesh->n_nodeNum; i++)
		{
			outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << p_blockMesh->v_node[i].x_ << " " << p_blockMesh->v_node[i].y_ << " ";
			for (int j = 0; j < (int)v_scalarField.size() - 1; j++)
			{
				Field<Scalar>& CurField = *v_scalarField[j];
				outFile << CurField.nodeField.GetValue(i) << " ";
			}
			if ((int)v_scalarField.size() > 0)
			{
				outFile << (*v_scalarField[(int)v_scalarField.size() - 1]).nodeField.GetValue(i) << " ";
			}

			for (int j = 0; j < (int)v_vectorField.size() - 1; j++)
			{
				Field<Vector>& CurField = *v_vectorField[j];
				outFile << CurField.nodeField.GetValue(i).x_ << " " << CurField.nodeField.GetValue(i).y_;
			}
			if ((int)v_vectorField.size() > 0)
			{
				outFile << (*v_vectorField[(int)v_vectorField.size() - 1]).nodeField.GetValue(i).x_ << " " << (*v_vectorField[(int)v_vectorField.size() - 1]).nodeField.GetValue(i).y_;
			}
			outFile << std::endl;
		}

		for (int i = 0; i < (int)p_blockMesh->n_elemNum; i++)
		{
			if (p_blockMesh->est_shapeType == Element::estTriangular)
			{
				outFile << p_blockMesh->v_elem[i].v_nodeID[0] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[1] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[2] + 1 << std::endl;
			}
			else if (p_blockMesh->est_shapeType == Element::estQuadrilateral)
			{
				outFile << p_blockMesh->v_elem[i].v_nodeID[0] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[1] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[2] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[3] + 1 << std::endl;
			}
			else if (p_blockMesh->est_shapeType == Element::estMixed)
			{
				if (p_blockMesh->v_elem[i].est_shapeType == Element::estTriangular)
				{
					outFile << p_blockMesh->v_elem[i].v_nodeID[0] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[1] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[2] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[2] + 1 << std::endl;
				}
				else if (p_blockMesh->v_elem[i].est_shapeType == Element::estQuadrilateral)
				{
					outFile << p_blockMesh->v_elem[i].v_nodeID[0] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[1] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[2] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[3] + 1 << std::endl;
				}
			}
		}

		// Write BC 
		for (int i = 0; i < (int)this->p_blockMesh->v_boundaryFaceZone.size(); i++)
		{
			outFile << "ZONE T = " << "\"" << this->p_blockMesh->v_boundaryFaceZone[i].name << "\"," << " DATAPACKING = POINT, N = " << this->p_blockMesh->n_nodeNum << ", E = "
				<< this->p_blockMesh->v_boundaryFaceZone[i].v_faceID.size() << ", ZONETYPE = FELINESEG, " << "VARSHARELIST = ([";

			for (int j = 0; j < nValNum - 1; j++)
			{
				outFile << j + 1 << ",";
			}

			outFile << nValNum << "]=1,[1,2])" << std::endl;

			for (int j = 0; j < (int)this->p_blockMesh->v_boundaryFaceZone[i].v_faceID.size(); j++)
			{
				int nFaceID = this->p_blockMesh->v_boundaryFaceZone[i].v_faceID[j];

				for (int n = 0; n < (int)this->p_blockMesh->v_face[nFaceID].v_nodeID.size(); n++)
				{
					outFile << this->p_blockMesh->v_face[nFaceID].v_nodeID[n] + 1 << "\t";
				}

				outFile << std::endl;
			}
		}
	}
	else if (this->p_blockMesh->md_meshDim == Mesh::md3D)
	{
		if (this->p_blockMesh->est_shapeType == Element::estPolyhedron)
		{
			//Not implimented
		}
		else
		{
			int nValNum = 3 + (int)v_scalarField.size() + (int)v_vectorField.size() * 3;

			outFile << "VARIABLES = " << "\"x [m]\"," << "\"y [m]\"," << "\"z [m]\",";
			if ((int)v_scalarField.size() > 0)
			{
				for (int i = 0; i < (int)v_scalarField.size() - 1; i++)
				{
					Field<Scalar>& CurField = *v_scalarField[i];
					outFile << "\"" << CurField.st_name << "\",";
				}
				outFile << "\"" << v_scalarField[(int)v_scalarField.size() - 1]->st_name << "\"";
			}

			if ((int)v_vectorField.size() > 0)
			{
				outFile << ",";
				for (int i = 0; i < (int)v_vectorField.size() - 1; i++)
				{
					Field<Vector>& CurField = *v_vectorField[i];
					outFile << "\"" << CurField.st_name << ".x\"," << "\"" << CurField.st_name << ".y\"," << "\"" << CurField.st_name << ".z\",";
				}
				outFile << "\"" << v_vectorField[(int)v_vectorField.size() - 1]->st_name << ".x\"," << "\"" << v_vectorField[(int)v_vectorField.size() - 1]->st_name << ".y\""
					<< "\"" << v_vectorField[(int)v_vectorField.size() - 1]->st_name << ".z\"";
			}
			outFile << std::endl;

			if (p_blockMesh->est_shapeType == Element::estTetrahedral)
			{
				outFile << "ZONE T = " << "\"" << "Tetrahedral_mesh" << "\"," << "STRANDID = 1, SOLUTIONTIME = " << Time << " DATAPACKING = POINT, N = " << p_blockMesh->n_nodeNum << ", E = " << p_blockMesh->n_elemNum
					<< ", ZONETYPE = FETETRAHEDRON" << std::endl;
			}
			else if (p_blockMesh->est_shapeType == Element::estHexahedral)
			{
				outFile << "ZONE T = " << "\"" << "Hexahedral_mesh" << "\"," << "STRANDID = 1, SOLUTIONTIME = " << Time << " DATAPACKING = POINT, N = " << p_blockMesh->n_nodeNum << ", E = " << p_blockMesh->n_elemNum
					<< ", ZONETYPE = FEBRICK" << std::endl;
			}
			else if (p_blockMesh->est_shapeType == Element::estMixed)
			{
				outFile << "ZONE T = " << "\"" << "Tet_mesh" << "\"," << "STRANDID = 1, SOLUTIONTIME = " << Time << " DATAPACKING = POINT, N = " << p_blockMesh->n_nodeNum << ", E = " << p_blockMesh->n_elemNum
					<< ", ZONETYPE = FEBRICK" << std::endl;
			}
			else if (p_blockMesh->est_shapeType == Element::estWedge)
			{
				outFile << "ZONE T = " << "\"" << "Wedge_mesh" << "\"," << "STRANDID = 1, SOLUTIONTIME = " << Time << " DATAPACKING = POINT, N = " << p_blockMesh->n_nodeNum << ", E = " << p_blockMesh->n_elemNum
					<< ", ZONETYPE = FEBRICK" << std::endl;
			}

			for (int i = 0; i < (int)p_blockMesh->n_nodeNum; i++)
			{
				outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << p_blockMesh->v_node[i].x_ << " " << p_blockMesh->v_node[i].y_ << " " << p_blockMesh->v_node[i].z_ << " ";
				for (int j = 0; j < (int)v_scalarField.size() - 1; j++)
				{
					Field<Scalar>& CurField = *v_scalarField[j];
					outFile << CurField.nodeField.GetValue(i) << " ";
				}
				if ((int)v_scalarField.size() > 0)
				{
					outFile << (*v_scalarField[(int)v_scalarField.size() - 1]).nodeField.GetValue(i) << " ";
				}

				for (int j = 0; j < (int)v_vectorField.size() - 1; j++)
				{
					Field<Vector>& CurField = *v_vectorField[j];
					outFile << CurField.nodeField.GetValue(i).x_ << " " << CurField.nodeField.GetValue(i).y_ << " " << CurField.nodeField.GetValue(i).z_;
				}
				if ((int)v_vectorField.size() > 0)
				{
					outFile << (*v_vectorField[(int)v_vectorField.size() - 1]).nodeField.GetValue(i).x_ << " " << (*v_vectorField[(int)v_vectorField.size() - 1]).nodeField.GetValue(i).y_
						<< " " << (*v_vectorField[(int)v_vectorField.size() - 1]).nodeField.GetValue(i).z_;
				}
				outFile << std::endl;
			}

			if (p_blockMesh->est_shapeType == Element::estTetrahedral)
			{
				for (int i = 0; i < (int)p_blockMesh->v_elem.size(); i++)
				{
					outFile << p_blockMesh->v_elem[i].v_nodeID[0] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[1] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[2] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[3] + 1 << std::endl;
				}
			}
			else if (p_blockMesh->est_shapeType == Element::estHexahedral)
			{
				for (int i = 0; i < (int)p_blockMesh->v_elem.size(); i++)
				{
					outFile << p_blockMesh->v_elem[i].v_nodeID[0] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[1] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[2] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[3] + 1 << " "
						<< p_blockMesh->v_elem[i].v_nodeID[4] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[5] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[6] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[7] + 1 << std::endl;
				}
			}
			else if (p_blockMesh->est_shapeType == Element::estMixed)
			{
				for (int i = 0; i < (int)p_blockMesh->v_elem.size(); i++)
				{
					if (p_blockMesh->v_elem[i].est_shapeType == Element::estTetrahedral)
					{
						outFile << p_blockMesh->v_elem[i].v_nodeID[0] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[1] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[2] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[2] + 1 << " "
							<< p_blockMesh->v_elem[i].v_nodeID[3] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[3] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[3] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[3] + 1 << std::endl;
					}
					else if (p_blockMesh->v_elem[i].est_shapeType == Element::estPyramid)
					{
						outFile << p_blockMesh->v_elem[i].v_nodeID[0] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[1] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[2] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[3] + 1 << " "
							<< p_blockMesh->v_elem[i].v_nodeID[4] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[4] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[4] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[4] + 1 << std::endl;
					}
					else if (p_blockMesh->v_elem[i].est_shapeType == Element::estWedge)
					{
						outFile << p_blockMesh->v_elem[i].v_nodeID[0] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[1] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[4] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[3] + 1 << " "
							<< p_blockMesh->v_elem[i].v_nodeID[2] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[2] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[5] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[5] + 1 << std::endl;
					}
					else if (p_blockMesh->v_elem[i].est_shapeType == Element::estHexahedral)
					{
						outFile << p_blockMesh->v_elem[i].v_nodeID[0] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[1] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[2] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[3] + 1 << " "
							<< p_blockMesh->v_elem[i].v_nodeID[4] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[5] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[6] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[7] + 1 << std::endl;
					}
				}
			}
			else if (p_blockMesh->est_shapeType == Element::estWedge)
			{
				for (int i = 0; i < (int)p_blockMesh->v_elem.size(); i++)
				{
					outFile << p_blockMesh->v_elem[i].v_nodeID[0] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[1] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[5] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[3] + 1 << " "
						<< p_blockMesh->v_elem[i].v_nodeID[2] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[2] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[4] + 1 << " " << p_blockMesh->v_elem[i].v_nodeID[4] + 1 << std::endl;
				}
			}

			// Write BC
			for (int i = 0; i < (int)this->p_blockMesh->v_boundaryFaceZone.size(); i++)
			{
				outFile << "ZONE T = " << "\"" << this->p_blockMesh->v_boundaryFaceZone[i].name << "\"," << " DATAPACKING = POINT, N = " << this->p_blockMesh->n_nodeNum << ", E = "
					<< this->p_blockMesh->v_boundaryFaceZone[i].v_faceID.size() << ", ZONETYPE = FEQUADRILATERAL, " << "VARSHARELIST = ([";
				for (int j = 0; j < nValNum - 1; j++)
				{
					outFile << j + 1 << ",";
				}
				outFile << nValNum << "]=1,[1,2,3])" << std::endl;

				for (int j = 0; j < (int)this->p_blockMesh->v_boundaryFaceZone[i].v_faceID.size(); j++)
				{
					int nFaceID = this->p_blockMesh->v_boundaryFaceZone[i].v_faceID[j];

					for (int n = 0; n < (int)this->p_blockMesh->v_face[nFaceID].v_nodeID.size(); n++)
					{
						outFile << this->p_blockMesh->v_face[nFaceID].v_nodeID[n] + 1 << "\t";
					}

					if (this->p_blockMesh->v_face[nFaceID].ft_faceType == Face::ftTrangular)
					{
						// ZONETYPE = FEQUADRILATERAL need four Node to display trangle, here output the last again
						outFile << this->p_blockMesh->v_face[nFaceID].v_nodeID[(int)this->p_blockMesh->v_face[nFaceID].v_nodeID.size() - 1] + 1 << "\t";
					}
					outFile << std::endl;
				}
			}
		}
	}

	outFile.close();
}

void FieldIO::FieldElementToNode()
{
	for (int i = 0; i < (int)v_scalarField.size(); i++)
	{
		v_scalarField[i]->ElementToNode();
	}

	for (int i = 0; i < (int)v_vectorField.size(); i++)
	{
		v_vectorField[i]->ElementToNode();
	}
}


void FieldIO::WriteVTKField(const std::string& outFileName)
{
	std::ofstream outFile(outFileName.c_str());
	outFile << "# vtk DataFile Version 2.0" << std::endl;
	outFile << "Unstructured Grid VTK" << std::endl;
	outFile << "ASCII" << std::endl;
	outFile << "DATASET UNSTRUCTURED_GRID" << std::endl;
	outFile << "POINTS " << p_blockMesh->n_nodeNum << " float" << std::endl;
	// Write Points
	for (unsigned int i = 0; (int)i < p_blockMesh->n_nodeNum; i++)
	{
		outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << p_blockMesh->v_node[i].x_ << " " << p_blockMesh->v_node[i].y_ << " " << p_blockMesh->v_node[i].z_ << std::endl;
	}

	int nCount(p_blockMesh->n_elemNum);
	for (int i = 0; i < (int)p_blockMesh->v_elem.size(); i++)
	{
		nCount += (int)(p_blockMesh->v_elem[i].v_nodeID.size());
	}
	outFile << "CELLS " << p_blockMesh->n_elemNum << " " << nCount << std::endl;

	// Write element composiation
	for (int i = 0; i < (int)p_blockMesh->v_elem.size(); i++)
	{
		if (p_blockMesh->v_elem[i].est_shapeType == Element::estTriangular)
		{
			outFile << "3 " << p_blockMesh->v_elem[i].v_nodeID[0] << " " << p_blockMesh->v_elem[i].v_nodeID[1] << " " << p_blockMesh->v_elem[i].v_nodeID[2] << std::endl;
		}
		else if (p_blockMesh->v_elem[i].est_shapeType == Element::estQuadrilateral)
		{
			outFile << "4 " << p_blockMesh->v_elem[i].v_nodeID[0] << " " << p_blockMesh->v_elem[i].v_nodeID[1] << " " << p_blockMesh->v_elem[i].v_nodeID[2] << " " << p_blockMesh->v_elem[i].v_nodeID[3] << std::endl;
		}
		else if (p_blockMesh->v_elem[i].est_shapeType == Element::estTetrahedral)
		{
			outFile << "4 " << p_blockMesh->v_elem[i].v_nodeID[0] << " " << p_blockMesh->v_elem[i].v_nodeID[1] << " " << p_blockMesh->v_elem[i].v_nodeID[2] << " " << p_blockMesh->v_elem[i].v_nodeID[3] << std::endl;
		}
		else if (p_blockMesh->v_elem[i].est_shapeType == Element::estPyramid)
		{
			outFile << "5 " << p_blockMesh->v_elem[i].v_nodeID[0] << " " << p_blockMesh->v_elem[i].v_nodeID[1] << " " << p_blockMesh->v_elem[i].v_nodeID[2] << " " << p_blockMesh->v_elem[i].v_nodeID[3] << " "
				<< p_blockMesh->v_elem[i].v_nodeID[4] << std::endl;
		}
		else if (p_blockMesh->v_elem[i].est_shapeType == Element::estWedge)
		{
			outFile << "6 " << p_blockMesh->v_elem[i].v_nodeID[0] << " " << p_blockMesh->v_elem[i].v_nodeID[1] << " " << p_blockMesh->v_elem[i].v_nodeID[2] << " " << p_blockMesh->v_elem[i].v_nodeID[3] << " "
				<< p_blockMesh->v_elem[i].v_nodeID[4] << " " << p_blockMesh->v_elem[i].v_nodeID[5] << std::endl;
		}
		else if (p_blockMesh->v_elem[i].est_shapeType == Element::estHexahedral)
		{
			outFile << "8 " << p_blockMesh->v_elem[i].v_nodeID[0] << " " << p_blockMesh->v_elem[i].v_nodeID[1] << " " << p_blockMesh->v_elem[i].v_nodeID[2] << " " << p_blockMesh->v_elem[i].v_nodeID[3] << " "
				<< p_blockMesh->v_elem[i].v_nodeID[4] << " " << p_blockMesh->v_elem[i].v_nodeID[5] << " " << p_blockMesh->v_elem[i].v_nodeID[6] << " " << p_blockMesh->v_elem[i].v_nodeID[7] << std::endl;
		}
	}

	// Write element type
	outFile << "CELL_TYPES " << p_blockMesh->n_elemNum << std::endl;
	for (int i = 0; i < (int)p_blockMesh->v_elem.size(); i++)
	{
		if (p_blockMesh->v_elem[i].est_shapeType == Element::estTriangular)
		{
			outFile << "5" << std::endl;
		}
		else if (p_blockMesh->v_elem[i].est_shapeType == Element::estQuadrilateral)
		{
			outFile << "9" << std::endl;
		}
		else if (p_blockMesh->v_elem[i].est_shapeType == Element::estTetrahedral)
		{
			outFile << "10" << std::endl;
		}
		else if (p_blockMesh->v_elem[i].est_shapeType == Element::estPyramid)
		{
			outFile << "14" << std::endl;
		}
		else if (p_blockMesh->v_elem[i].est_shapeType == Element::estWedge)
		{
			outFile << "13" << std::endl;
		}
		else if (p_blockMesh->v_elem[i].est_shapeType == Element::estHexahedral)
		{
			outFile << "12" << std::endl;
		}
	}

	// Write field 
	outFile << "CELL_DATA " << p_blockMesh->v_elem.size() << std::endl;

	for (unsigned int sField = 0 ; sField < (int)v_scalarField.size(); sField++)
	{
		outFile << "SCALARS"<< " "<< (*v_scalarField[sField]).st_name << " " << "float " << std::endl;

		outFile << "LOOKUP_TABLE " << "default" << std::endl;

		for (int i = 0; i < (int)p_blockMesh->v_elem.size(); i++)
		{
			outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) <<
				(*v_scalarField[sField]).elementField.GetValue(i) <<  std::endl;

		}
	}


	for (unsigned int vField = 0; vField < (int)v_vectorField.size(); vField++)
	{
		outFile << "VECTORS"<< " " << (*v_vectorField[vField]).st_name << " " << "float "<< std::endl;


		for (int i = 0; i < (int)p_blockMesh->n_nodeNum; i++)
		{
			outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) <<
				(*v_vectorField[vField]).nodeField.GetValue(i).x_<<" "
				<< (*v_vectorField[vField]).nodeField.GetValue(i).y_ << " "
				<< (*v_vectorField[vField]).nodeField.GetValue(i).z_ << std::endl;

		}
	}
	outFile.close();
}

