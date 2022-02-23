/*---------------------------------------------------------------------------*\
Class:
Field
	File Name:
Field.h
	Description:
Definition and manipulations on field

	Author:		Kong Ling
	Date: 2016-11-24

Revised:
	Description:
	1.Made modification in Construction function(changed "const std::auto_ptr<Mesh>&" to "Mesh*")
	2.Made correction in bugs which caulse by modification(1.)

	Revisor:	Shuai Zhang
	Modified Date: 2017-01-12

Revised:
	Description:
	0. Happy National Day!
	1. Adding st_name as a member data of Field<Type>
	2. Adding an sting-type input parameter in constructor functions for the naming of a created field.
	
	Revisor:	Kong Ling
	Modified Date: 2017-10-01
	
Revised:
	Description:
	1. Convert "non const" input parameters to const;
	2. Convert "non const" function to const function;
	3. Correct UpdateBoundary::fbtCaclulated else from "v_BoundaryFaceField[i].KeepCurrentValue(faceField);" to 
	"v_BoundaryFaceField[i].CopyFromElement(elementField, faceField);"
	
	Revisor:		Shuai Zhang
	Modified Date:	2018-11-28
	
Revised:
	Description:
	1. Convert "static_cast<Field<Type> > (rhs).elementField.Assignment()" To
	"rhs.elementField.Assignment()";
	2. "boundaryLeft" option have been made in Field<Type>::UpdateBoundary;
	
	Revisor:		Shuai Zhang, Kong Ling
	Modified Date:	2018-11-28
\*---------------------------------------------------------------------------*/

#include "../MHT_field/Field.h"

#include <map>
#include <iomanip>
#include <sstream>
#include "../MHT_common/SystemControl.h"
#include "../MHT_common/InterpolationTools.h"
#include "../MHT_mesh/Mesh.h"
#include "../MHT_field/BoundaryCalculationUtility.h"

// =======================================Scalar Type Field=============================================

//creating boundary fields according to mesh information
template<>
void Field<Scalar>::CreateBoundaryField()
{
    for (int i = 0; i < (int)p_blockMesh->v_boundaryFaceZone.size(); i++)
    {
        FaceZone* pfz = &(p_blockMesh->v_boundaryFaceZone[i]);
        this->v_BoundaryFaceField.push_back(BoundaryFaceZoneField<Scalar>(p_blockMesh, pfz));
    }
}

//constructor
template<>
Field<Scalar>::Field
(
Mesh* UGBPointer,
const std::string& fieldName
):
p_blockMesh(UGBPointer),
st_name(fieldName),
elementField(UGBPointer),
elementField0(UGBPointer),
faceField(UGBPointer),
nodeField(UGBPointer),
interiorFaceField(UGBPointer)
{
    this->CreateBoundaryField();
}

//constructor
template<>
Field<Scalar>::Field
(
Mesh* UGBPointer,
Scalar fieldValue,
const std::string& fieldName
)
:
p_blockMesh(UGBPointer),
st_name(fieldName),
elementField(UGBPointer),
elementField0(UGBPointer),
faceField(UGBPointer),
nodeField(UGBPointer),
interiorFaceField(UGBPointer)
{
	elementField.Initialize(fieldValue);
    this->CreateBoundaryField();
}

//constructor
template<>
Field<Scalar>::Field
(
Mesh* UGBPointer,
const std::vector<Scalar>& valueList,
const std::string& fieldName
):
p_blockMesh(UGBPointer),
st_name(fieldName),
elementField(UGBPointer),
elementField0(UGBPointer),
faceField(UGBPointer),
nodeField(UGBPointer),
interiorFaceField(UGBPointer)
{
	elementField.Initialize(valueList);
    this->CreateBoundaryField();
}

//constructor
template<>
Field<Scalar>::Field
(
Mesh* UGBPointer,
Scalar(udf)(Scalar, Scalar, Scalar),
const std::string& fieldName
)
:
p_blockMesh(UGBPointer),
st_name(fieldName),
elementField(UGBPointer),
elementField0(UGBPointer),
faceField(UGBPointer),
nodeField(UGBPointer),
interiorFaceField(UGBPointer)
{
	elementField.Initialize(udf);
    this->CreateBoundaryField();
}

//specify a special boundary condition on a boundary face zone
template<>
void Field<Scalar>::SetBoundaryCondition
(
const std::string& patchName,
const std::string& boundaryType
)
{
	int faceZoneID = this->p_blockMesh->GetBoundaryFaceZoneID(patchName);
	if (-1 != faceZoneID)
	{
		if ("boundaryLeft" == boundaryType)
		{
			this->v_BoundaryFaceField[faceZoneID].SetSpecialBC(BoundaryFaceZoneField<Scalar>::fbtCaclulated);
            std::cout << "boundaryLeft condition successfully given on boundary: " << patchName << " of " << this->st_name << std::endl;
		}
		else if ("copy" == boundaryType)
		{
			this->v_BoundaryFaceField[faceZoneID].SetSpecialBC(BoundaryFaceZoneField<Scalar>::fbtCopy);
           // std::cout << "Copy condition successfully given on boundary: " << patchName << " of " << this->st_name << std::endl;
		}
		else if ("extrapolate" == boundaryType)
		{
			this->v_BoundaryFaceField[faceZoneID].SetSpecialBC(BoundaryFaceZoneField<Scalar>::fbtExtrapolate);
            std::cout << "Extrapolate condition successfully given on boundary: " << patchName << " of " << this->st_name << std::endl;
		}
		else
		{
			WarningContinue("unknown boundary type:" + boundaryType);
		}
	}
	return;
}

//creating v_external with given boundary condition according to external face zones in corresponding mesh
template<>
void Field<Scalar>::SetBoundaryCondition
(
const std::string& patchName,
Scalar bc_a,
Scalar bc_b,
Scalar bc_c
)
{
	int faceZoneID = this->p_blockMesh->GetBoundaryFaceZoneID(patchName);
	if (-1 != faceZoneID)
	{
		this->v_BoundaryFaceField[faceZoneID].SetNumerBC(bc_a, bc_b, bc_c);
        std::cout << "Boundary condition successfully given on boundary: " << patchName << " of " << this->st_name << std::endl;
	}
	return;
}

//specify a uniform general-type (ABC-form) boundary condition on a boundary face zone
template<>
void Field<Scalar>::SetBoundaryCondition
(
const std::string& patchName,
Scalar bc_a,
Scalar bc_b,
Scalar(*udf)(Scalar x, Scalar y, Scalar z)
)
{
	int faceZoneID = this->p_blockMesh->GetBoundaryFaceZoneID(patchName);
	if (-1 != faceZoneID)
	{
		this->v_BoundaryFaceField[faceZoneID].SetNumerBC(bc_a, bc_b, udf);
        std::cout << "Boundary condition successfully given on boundary: " << patchName << " of " << this->st_name << std::endl;
	}
	return;
}

//specify a uniform general-type (ABC-form) boundary condition on a boundary face zone
template<>
void Field<Scalar>::SetBoundaryCondition
(
const std::string& patchName,
Scalar bc_a,
Scalar bc_b,
std::vector<Scalar > bc_c_list
)
{
	int faceZoneID = this->p_blockMesh->GetBoundaryFaceZoneID(patchName);
	if (-1 != faceZoneID)
	{
        std::vector<NumerBC<Scalar> > abcList;
		for (int i = 0; i < (int)bc_c_list.size(); i++)
		{
			abcList.push_back(NumerBC<Scalar>(bc_a, bc_b, bc_c_list[i]));
		}
		this->v_BoundaryFaceField[faceZoneID].SetNumerBC(abcList);
        std::cout << "Boundary condition successfully given on boundary: " << patchName << std::endl;
	}
	return;
}

//check wether all boundaries have been correctly specified
template<>
void Field<Scalar>::CheckBoundaryCondition() const
{
    std::vector<std::string> v_facesWithoutBC;
	//collect names of faces with no boundary condtion given on;
	for (int i = 0; i < (int)v_BoundaryFaceField.size(); i++)
	{


		if (!(v_BoundaryFaceField[i].BCExist()))
		{
			v_facesWithoutBC.push_back(v_BoundaryFaceField[i].p_faceZone->name);
		}
	}
	int NumOfFacesWithoutBC = (int)v_facesWithoutBC.size();
	//if there exist face without boundary condition, display fatal error with a list of no-BC face names included.
	if (0 != NumOfFacesWithoutBC)
	{
        std::stringstream message;
        message << "For field " + this->st_name + ", boundary conidtions are not found on faces:" << std::endl;
		for (int i = 0; i < (int)v_facesWithoutBC.size(); i++)
		{
            message << "(" << i + 1 << ") " << v_facesWithoutBC[i] << std::endl;
		}
		FatalError(message.str());
	}
}

//copy boundary conditions to another field
template<>
void Field<Scalar>::CopyBoundaryConditionTo(Field<Scalar>& targetField)
{
	if (targetField.p_blockMesh != this->p_blockMesh)
	{
        std::stringstream message;
		message << "The host field " << this->st_name << " and the target one ";
		message << targetField.st_name << " are not built on the same mesh.";
		FatalError(message.str());
	}
    for (int i = 0; i < (int)this->v_BoundaryFaceField.size(); i++)
	{
		this->v_BoundaryFaceField[i].CopyBoundaryConditionTo(targetField.v_BoundaryFaceField[i]);
	}
}

//copying elementField to elementField0
template<>
void Field<Scalar>::SaveOld()
{
	if (false == this->elementField.Assignment())
	{
		FatalError("Cannot save old for not passing Assignment Check");
	}
	this->elementField0.v_value.assign(this->elementField.v_value.begin(), this->elementField.v_value.end());
	this->elementField0.fs_status = BaseField<Scalar>::fsAssigned;
}

//calculate face values on boundary
template<>
void Field<Scalar>::UpdateBoundary(const std::string& option)
{
	if (false == this->faceField.Assignment() || false == this->elementField.Assignment())
	{
		FatalError("Cannot UpdateBoundary for not passing Assignment Check, element- and face-field both need be available");
		return;
	}
	if ("boundaryLeft" == option)
	{
		for (int i = 0; i < (int)v_BoundaryFaceField.size(); i++)
		{
			v_BoundaryFaceField[i].KeepCurrentValue(faceField);
		}
	}
	else if ("copy" == option)
	{
		for (int i = 0; i < (int)v_BoundaryFaceField.size(); i++)
		{
			v_BoundaryFaceField[i].CopyFromElement(elementField, faceField);
		}
	}
	else
	{
		//Step 1: set initial value for at boundaries
		for (int i = 0; i < (int)v_BoundaryFaceField.size(); i++)
		{
			if (BoundaryFaceZoneField<Scalar>::fbtCaclulated== v_BoundaryFaceField[i].fbt_BCType)
			{
				v_BoundaryFaceField[i].KeepCurrentValue(faceField);
			}
			else
			{
				v_BoundaryFaceField[i].CopyFromElement(elementField, faceField);
			}
		}
		//Step 2: calculate the values according to the given boundary conditions
		for (int i = 0; i < (int)v_BoundaryFaceField.size(); i++)
		{
			v_BoundaryFaceField[i].CalculateValue(elementField, faceField, option);
		}
	}
	return;
}


void ElementToFaceWithAllParameters
(
Field<Scalar>& phi,
const std::string& interiorApproach,
const std::string& boundaryApproach
)
{
	//step 3: interpolate values on interior faces
	phi.interiorFaceField.CalculateValue(phi.elementField, phi.faceField, interiorApproach);
	//step 4: calculate values on boundary according to selected method
	phi.UpdateBoundary(boundaryApproach);
}

//interpolate from element to face. On boundaris, values are calculated by the given option
template<>
void Field<Scalar>::ElementToFace(const std::string& approach)
{
	ElementToFaceCheck(this->elementField,this->faceField);
    std::string interiorParameter = "linear";
    std::string boundaryParameter = "default";
	int parameter = ElementToFaceParameterCheck(approach);
	if (-1 == parameter)
	{
		FatalError("illegal parameter in ElementToFace: " + approach);
	}
	else if (0 == parameter)
	{
		interiorParameter = approach;
	}
	else if (1 == parameter)
	{
		boundaryParameter = approach;
	}
	ElementToFaceWithAllParameters(*this, interiorParameter, boundaryParameter);
	return;
}


//interpolate from element to face. On boundaris, values are calculated by the given option
template<>
void Field<Scalar>::ElementToFace(const std::string& option1, const std::string& option2)
{
	ElementToFaceCheck(this->elementField, this->faceField);
    std::string interiorParameter = "linear";
    std::string boundaryParameter = "default";
    std::vector<std::string> v_option;
	v_option.push_back(option1);
	v_option.push_back(option2);
	for (int i = 0; i < (int)v_option.size(); i++)
	{
		int parameter = ElementToFaceParameterCheck(v_option[i]);
		if (-1 == parameter)
		{
			FatalError("illegal parameter in ElementToFace: " + v_option[i]);
		}
		else if (0 == parameter)
		{
			interiorParameter = v_option[i];
		}
		else if (1 == parameter)
		{
			boundaryParameter = v_option[i];
		}
	}
	ElementToFaceWithAllParameters(*this, interiorParameter, boundaryParameter);
	return;
}

template<>
void Field<Scalar>::NodeToElement()
{
	//step 1: check whether the element field has values on it
	if (false == this->nodeField.Assignment())
	{
		FatalError("Cannot proceed NodeToElement for not passing Assignment Check");
	}
	//step 2: create elementfield if not existing
	if (BaseField<Scalar>::fsNotExist == nodeField.fs_status)
	{
		elementField.Initialize();
	}
	//step 3: visit all nodes in mesh data
	for (int i = 0; i < p_blockMesh->n_elemNum; i++)
	{
        std::vector<Node> relativePos;
        std::vector<Scalar> value;
		Vector cellCenter = p_blockMesh->v_elem[i].center;
		for (int j = 0;j < p_blockMesh->v_elem[i].v_nodeID.size();j++)
		{
			int nodeID = p_blockMesh->v_elem[i].v_nodeID[j];
			Scalar nodeValue = this->nodeField.GetValue(nodeID);
			relativePos.push_back(p_blockMesh->v_node[nodeID] - cellCenter);
			value.push_back(nodeValue);
		}
		Scalar interpolatedValue = Interpolation(relativePos, value);
		elementField.SetValue(i, interpolatedValue);
	}
	elementField.fs_status = BaseField<Scalar>::fsAssigned;
	return;
}

template<>
void Field<Scalar>::ElementToNode()
{
	//step 1: check whether the element field has values on it
	if (false == this->elementField.Assignment())
	{
		FatalError("Cannot proceed ElementToNode for not passing Assignment Check");
	}
	//step 2: create node field if not existing
	if (BaseField<Scalar>::fsNotExist == nodeField.fs_status)
	{
		nodeField.Initialize();
	}
	//step 3: visit all nodes in mesh data

	for (int i = 0; i < (int)p_blockMesh->v_vertice.size(); i++)
	{
		std::vector<Node> relativePos;
		std::vector<Scalar> value;
		std::vector<bool> isBoundary;
		//step 3.1: collect topological information from element field
		for (int j = 0; j < (int)p_blockMesh->v_vertice[i].v_elemID.size(); j++)
		{
			int elemID = p_blockMesh->v_vertice[i].v_elemID[j];
			Node elementCenter = p_blockMesh->v_elem[elemID].center;
			relativePos.push_back(elementCenter - p_blockMesh->v_node[i]);
			value.push_back(elementField.GetValue(elemID));
			//an element center is of course not on the boundary; 
			isBoundary.push_back(false);
		}
		std::cout << i << " , " << p_blockMesh->v_vertice[i].v_elemID.size() << std::endl;
		if (relativePos.size() == 0)
		{
			std::cout<< i <<" , "<< p_blockMesh->v_vertice[i].v_elemID.size() << std::endl;
			system("pause");
		}

		//step 3.2: if existing, collect toplogical information from boundary face field
		if (BaseField<Scalar>::fsAssigned == faceField.fs_status)
		{
			for (int j = 0; j < (int)p_blockMesh->v_vertice[i].v_faceID.size(); j++)
			{
				int faceID = p_blockMesh->v_vertice[i].v_faceID[j];
				//check wether this face is on a boundray and push back values correspondingly;
				if (-1 == p_blockMesh->SearchOwnerNeighbor(faceID).second)
				{
					Node faceCenter = p_blockMesh->v_face[faceID].center;
					relativePos.push_back(faceCenter - p_blockMesh->v_node[i]);
					value.push_back(faceField.GetValue(faceID));
					isBoundary.push_back(true);
				}
			}
		}
		//step 3.3: check wether we have interplating points located on boundary.
		//If having, delete those interpolating points at interior;
		bool hasBoundary = false;
		for (int j = 0; j < (int)isBoundary.size(); j++)
		{
			if (true == isBoundary[j])
			{
				hasBoundary = true;
				break;
			}
		}
		if (true == hasBoundary)
		{
			std::vector<Node> relativePos_backup;
			std::vector<Scalar> value_backup;
			for (int j = 0; j < (int)isBoundary.size(); j++)
			{
				relativePos_backup.push_back(relativePos[j]);
				value_backup.push_back(value[j]);
			}
			relativePos.clear();
			value.clear();
			for (int j = 0; j < (int)isBoundary.size(); j++)
			{
				if (true == isBoundary[j])
				{
					relativePos.push_back(relativePos_backup[j]);
					value.push_back(value_backup[j]);
				}
			}
		}
		//step 3.4: use interpolation tools and write values
		Scalar interpolatedValue = Interpolation(relativePos, value);
		nodeField.SetValue(i, interpolatedValue);
	}

	nodeField.fs_status = BaseField<Scalar>::fsAssigned;
}

//Kong Ling supplemented on 2019/12/23
template<>
std::vector<Scalar> Field<Scalar>::ExtractBoundary(const std::string& boundaryName)
{
	return this->faceField.ExtractBoundary(boundaryName);
}

//Kong Ling modified on 2019/12/23
//load values on a specific boundary with a ordered value list
template<>
void Field<Scalar>::LoadBoundary(const std::string& boundaryName,const std::vector<Scalar>& valueList)
{
	this->faceField.LoadBoundary(boundaryName, valueList);
	return;
}

template<>
Field<Scalar>& Field<Scalar>::operator = (const Field<Scalar>& rhs)
{
    if (true == rhs.nodeField.Assignment())
    {
        this->nodeField = rhs.nodeField;
    }
    if (true == rhs.faceField.Assignment())
    {
        this->faceField = rhs.faceField;
    }
    if (true == rhs.elementField.Assignment())
    {
        this->elementField = rhs.elementField;
    }
    return *this;
}

//write output file for a scalar field
template<>
void Field<Scalar>::WriteTecplotField(const std::string& outMshFileName)
{
	this->ElementToNode();
    std::ofstream	outFile(outMshFileName.c_str());

    outFile << "TITLE =\"" << this->st_name << "\"" << std::endl;

	if (p_blockMesh->md_meshDim == Mesh::md2D)
	{
		if (this->st_name!="")
		{
			outFile << "VARIABLES = " << "\"x [m]\"," << "\"y [m]\"," << "\""<<this->st_name << "\"" << std::endl;
		}
		else
		{
			outFile << "VARIABLES = " << "\"x [m]\"," << "\"y [m]\"," << "\"Field Var\"" << std::endl;
		}
 

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
            outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << p_blockMesh->v_node[i].x_ << " " << p_blockMesh->v_node[i].y_ << " " << this->nodeField.GetValue(i) << std::endl;
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
				<< this->p_blockMesh->v_boundaryFaceZone[i].v_faceID.size() << ", ZONETYPE = FELINESEG, " << "VARSHARELIST = ([1,2,3]=1,[1,2])" << std::endl;

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
			int nFaceNodeNum(0);
			for (int i = 0; i < (int)this->p_blockMesh->v_face.size(); i++)
			{
				for (int j = 0; j < (int)this->p_blockMesh->v_face[i].v_nodeID.size(); j++)
				{
					nFaceNodeNum++;
				}
			}

			outFile << "VARIABLES = " << std::endl << "\"X\"" << std::endl << "\"Y\"" << std::endl << "\"Z\"" << std::endl << "\"Phi\"" << std::endl;
			outFile << "DATASETAUXDATA Common.VectorVarsAreVelocity = \"TRUE\"" << std::endl;
			outFile << "ZONE T = \"fluid\"" << std::endl;
			outFile << "STRANDID=0, SOLUTIONTIME=0" << std::endl;
			outFile << "Nodes=" << this->p_blockMesh->n_nodeNum << ", Faces=" << this->p_blockMesh->n_faceNum << ", Elements=" << this->p_blockMesh->n_elemNum << ", ZONETYPE=FEPolyhedron" << std::endl;
			outFile << "DATAPACKING=BLOCK" << std::endl;
			outFile << "TotalNumFaceNodes=" << nFaceNodeNum << ", NumConnectedBoundaryFaces=0, TotalNumBoundaryConnections=0" << std::endl;
			outFile << "DT=(SINGLE SINGLE SINGLE )";

			// Write node X
			for (int i = 0; i < (int) this->p_blockMesh->v_node.size(); i++)
			{
				if (i % 5 == 0) outFile << std::endl;
                outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << this->p_blockMesh->v_node[i].x_ << " ";
			}

			// Write node Y
			for (int i = 0; i < (int) this->p_blockMesh->v_node.size(); i++)
			{
				if (i % 5 == 0) outFile << std::endl;
                outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << this->p_blockMesh->v_node[i].y_ << " ";
			}

			// Write node Z
			for (int i = 0; i < (int) this->p_blockMesh->v_node.size(); i++)
			{
				if (i % 5 == 0) outFile << std::endl;
                outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << this->p_blockMesh->v_node[i].z_ << " ";
			}

			// Write node Phi
			for (int i = 0; i < (int) this->p_blockMesh->v_node.size(); i++)
			{
				if (i % 5 == 0) outFile << std::endl;
                outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << this->nodeField.v_value[i] << " ";
			}

			outFile << std::endl << "# node count per face";
			for (int i = 0; i < (int) this->p_blockMesh->v_face.size(); i++)
			{
				if (i % 10 == 0) outFile << std::endl;
				outFile << this->p_blockMesh->v_face[i].v_nodeID.size() << " ";
			}

			int nFaceNodeCheck(0);
			outFile << std::endl << "# face nodes";
			for (int i = 0; i < (int) this->p_blockMesh->v_face.size(); i++)
			{
				for (int j = 0; j < (int)this->p_blockMesh->v_face[i].v_nodeID.size(); j++)
				{
					if (nFaceNodeCheck % 10 == 0) outFile << std::endl;
					outFile << this->p_blockMesh->v_face[i].v_nodeID[j] + 1 << " ";
					nFaceNodeCheck++;
				}
			}

			outFile << std::endl << "# left elements";
			for (int i = 0; i < (int) this->p_blockMesh->v_face.size(); i++)
			{
				if (i % 10 == 0) outFile << std::endl;
				outFile << this->p_blockMesh->v_face[i].n_owner + 1 << " ";
			}

			outFile << std::endl << "# node count per face";
			for (int i = 0; i < (int) this->p_blockMesh->v_face.size(); i++)
			{
				if (i % 10 == 0) outFile << std::endl;
				outFile << this->p_blockMesh->v_face[i].n_neighbor + 1 << " ";
			}
		}
		else
		{
			if (this->st_name!="")
			{
				outFile << "VARIABLES = " << "\"x [m]\"," << "\"y [m]\"," << "\"z [m]\"" << "\"" << this->st_name << "\"" << std::endl;
			}
			else
			{
				outFile << "VARIABLES = " << "\"x [m]\"," << "\"y [m]\"," << "\"z [m]\"" << "\"Field Var\"" << std::endl;
			}

            

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
                outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << p_blockMesh->v_node[i].x_ << " " << p_blockMesh->v_node[i].y_ << " " << p_blockMesh->v_node[i].z_ << " " << this->nodeField.GetValue(i) << std::endl;
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
					<< this->p_blockMesh->v_boundaryFaceZone[i].v_faceID.size() << ", ZONETYPE = FEQUADRILATERAL, " << "VARSHARELIST = ([1,2,3,4]=1,[1,2,3])" << std::endl;

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

//Writing BC Tecplot file for visualization
template<>
void Field<Scalar>::WriteTecplotField(const std::string& patchName, const std::string& outMshFileName)
{
	//Set v_nodeField number
	this->ElementToNode();

	int found = -1;
	for (int i = 0; i < (int)this->v_BoundaryFaceField.size(); i++)
	{
		if (patchName == this->v_BoundaryFaceField[i].p_faceZone->name)
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
		//Face ID list in current BC
        std::vector<int>& vFaceID = this->v_BoundaryFaceField[found].p_faceZone->v_faceID;
        std::vector<Face> vFaceOutput;
		vFaceOutput.resize((int)vFaceID.size());

        std::vector<Node> vNodeOutput;
        std::vector<Scalar> vNodeValue;
        std::map<int, int> NodeID_LocalID;
        std::map<int, int>::iterator iter;

		for (int i = 0; i < (int)vFaceID.size(); i++)
		{
			int curFaceID = vFaceID[i];
			Face& curFace = this->p_blockMesh->v_face[curFaceID];

			vFaceOutput[i].ft_faceType = curFace.ft_faceType;
			vFaceOutput[i].v_nodeID.resize((int)curFace.v_nodeID.size());

			for (int j = 0; j < (int)curFace.v_nodeID.size(); j++)
			{
				int curNodeID = curFace.v_nodeID[j];

				iter = NodeID_LocalID.begin();
				iter = NodeID_LocalID.find(curNodeID);

				//find global Node ID in map::NodeID_LocalID
				if (iter != NodeID_LocalID.end())
				{
					vFaceOutput[i].v_nodeID[j] = iter->second;
				}
				else
				{
					vFaceOutput[i].v_nodeID[j] = (int)vNodeOutput.size();

					vNodeOutput.push_back(this->p_blockMesh->v_node[curNodeID]);
					vNodeValue.push_back(this->nodeField.GetValue(curNodeID));
					NodeID_LocalID.insert(std::pair<int, int>(curNodeID, vFaceOutput[i].v_nodeID[j]));
				}
			}
		}

        std::ofstream	outFile(outMshFileName.c_str());

		if (this->p_blockMesh->md_meshDim == Mesh::md2D)
		{
            outFile << "TITLE =\"" << this->v_BoundaryFaceField[found].p_faceZone->name << "\"" << std::endl;
            outFile << "VARIABLES = " << "\"x [m]\"," << "\"y [m]\"," << "\"value\"" << std::endl;

			outFile << "ZONE T = " << "\"" << this->v_BoundaryFaceField[found].p_faceZone->name << "\"," << " DATAPACKING = POINT, N = " << vNodeOutput.size() << ", E = "
				<< vFaceOutput.size() << ", ZONETYPE = FELINESEG, " << std::endl;

			for (unsigned int i = 0; i < vNodeOutput.size(); i++)
			{
                outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << vNodeOutput[i].x_ << " " << vNodeOutput[i].y_ << " " << vNodeValue[i] << std::endl;
			}

			for (int i = 0; i < (int)vFaceOutput.size(); i++)
			{
				for (int n = 0; n < (int)vFaceOutput[i].v_nodeID.size(); n++)
				{
					outFile << vFaceOutput[i].v_nodeID[n] + 1 << "\t";
				}

				outFile << std::endl;
			}
		}
		else if (this->p_blockMesh->md_meshDim == Mesh::md3D)
		{
            outFile << "TITLE =\"" << this->v_BoundaryFaceField[found].p_faceZone->name << "\"" << std::endl;
            outFile << "VARIABLES = " << "\"x [m]\"," << "\"y [m]\"," << "\"z [m]\"," << "\"value\"" << std::endl;

			outFile << "ZONE T = " << "\"" << this->v_BoundaryFaceField[found].p_faceZone->name << "\"," << " DATAPACKING = POINT, N = " << vNodeOutput.size() << ", E = "
				<< vFaceOutput.size() << ", ZONETYPE = FEQUADRILATERAL, " << std::endl;

			for (unsigned int i = 0; i < vNodeOutput.size(); i++)
			{
                outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << vNodeOutput[i].x_ << " " << vNodeOutput[i].y_ << " " << vNodeOutput[i].z_ << " " << vNodeValue[i] << std::endl;
			}

			for (int i = 0; i < (int)vFaceOutput.size(); i++)
			{
				for (int n = 0; n < (int)vFaceOutput[i].v_nodeID.size(); n++)
				{
					outFile << vFaceOutput[i].v_nodeID[n] + 1 << "\t";
				}

				if (vFaceOutput[i].ft_faceType == Face::ftTrangular)
				{
					outFile << vFaceOutput[i].v_nodeID[(int)vFaceOutput[i].v_nodeID.size() - 1] + 1 << "\t";
				}
				outFile << std::endl;
			}
		}

		outFile.close();
	}
}

//write output file for a scalar field
template<>
void Field<Scalar>::WriteTecplotField(const std::string& outMshFileName, outType outFileType)
{
	if (outFileType == otCellNode)
	{
		this->WriteTecplotField(outMshFileName.c_str());
	}
	else if (outFileType == otCellCenter)
	{
        std::ofstream	outFile(outMshFileName.c_str());

        outFile << "TITLE =\"" << this->st_name << "\"" << std::endl;

		if (p_blockMesh->md_meshDim == Mesh::md2D)
		{
            outFile << "VARIABLES = " << "\"x [m]\"," << "\"y [m]\"," << "\"Field Var\"" << std::endl;

			if (p_blockMesh->est_shapeType == Element::estTriangular)
			{
				outFile << "ZONE T = " << "\"" << "Tri_mesh" << "\"," << " DATAPACKING = BLOCK, N = " << p_blockMesh->n_nodeNum << ", E = " << p_blockMesh->n_elemNum
                    << ", ZONETYPE = FETRIANGLE" << std::endl;
                outFile << "VARLOCATION = ([1-2]=NODAL, [3] = CELLCENTERED)" << std::endl;
			}
			else if (p_blockMesh->est_shapeType == Element::estQuadrilateral)
			{
				outFile << "ZONE T = " << "\"" << "Quad_mesh" << "\"," << " DATAPACKING = BLOCK, N = " << p_blockMesh->n_nodeNum << ", E = " << p_blockMesh->n_elemNum
                    << ", ZONETYPE = FEQUADRILATERAL" << std::endl;
                outFile << "VARLOCATION = ([1-2]=NODAL, [3] = CELLCENTERED)" << std::endl;
			}
			else if (p_blockMesh->est_shapeType == Element::estMixed)
			{
				outFile << "ZONE T = " << "\"" << "Mix_mesh" << "\"," << " DATAPACKING = BLOCK, N = " << p_blockMesh->n_nodeNum << ", E = " << p_blockMesh->n_elemNum
                    << ", ZONETYPE = FEQUADRILATERAL" << std::endl;
                outFile << "VARLOCATION = ([1-2]=NODAL, [3] = CELLCENTERED)" << std::endl;
			}

			for (int i = 0; i < (int)p_blockMesh->n_nodeNum; i++)
			{
                outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << p_blockMesh->v_node[i].x_ << std::endl;
			}

			for (int i = 0; i < (int)p_blockMesh->n_nodeNum; i++)
			{
                outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << p_blockMesh->v_node[i].y_ << std::endl;
			}

			for (int i = 0; i < (int)p_blockMesh->n_elemNum; i++)
			{
                outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << this->elementField.GetValue(i) << std::endl;
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
		}
		else if (p_blockMesh->md_meshDim == Mesh::md3D)
		{
            outFile << "VARIABLES = " << "\"x [m]\"," << "\"y [m]\"," << "\"z [m]\"" << "\"Field Var\"" << std::endl;

			if (p_blockMesh->est_shapeType == Element::estTetrahedral)
			{
				outFile << "ZONE T = " << "\"" << "Tetrahedral_mesh" << "\"," << " DATAPACKING = BLOCK, N = " << p_blockMesh->n_nodeNum << ", E = " << p_blockMesh->n_elemNum
                    << ", ZONETYPE = FETETRAHEDRON" << std::endl;
                outFile << "VARLOCATION = ([1-3]=NODAL, [4] = CELLCENTERED)" << std::endl;
			}
			else if (p_blockMesh->est_shapeType == Element::estHexahedral)
			{
				outFile << "ZONE T = " << "\"" << "Hexahedral_mesh" << "\"," << " DATAPACKING = BLOCK, N = " << p_blockMesh->n_nodeNum << ", E = " << p_blockMesh->n_elemNum
                    << ", ZONETYPE = FEBRICK" << std::endl;
                outFile << "VARLOCATION = ([1-3]=NODAL, [4] = CELLCENTERED)" << std::endl;
			}
			else if (p_blockMesh->est_shapeType == Element::estMixed)
			{
				outFile << "ZONE T = " << "\"" << "Tet_mesh" << "\"," << " DATAPACKING = BLOCK, N = " << p_blockMesh->n_nodeNum << ", E = " << p_blockMesh->n_elemNum
                    << ", ZONETYPE = FEBRICK" << std::endl;
                outFile << "VARLOCATION = ([1-3]=NODAL, [4] = CELLCENTERED)" << std::endl;
			}
			else if (p_blockMesh->est_shapeType == Element::estWedge)
			{
				outFile << "ZONE T = " << "\"" << "Wedge_mesh" << "\"," << " DATAPACKING = BLOCK, N = " << p_blockMesh->n_nodeNum << ", E = " << p_blockMesh->n_elemNum
                    << ", ZONETYPE = FEBRICK" << std::endl;
                outFile << "VARLOCATION = ([1-3]=NODAL, [4] = CELLCENTERED)" << std::endl;
			}

			for (int i = 0; i < (int)p_blockMesh->n_nodeNum; i++)
			{
                outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << p_blockMesh->v_node[i].x_ << std::endl;
			}

			for (int i = 0; i < (int)p_blockMesh->n_nodeNum; i++)
			{
                outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << p_blockMesh->v_node[i].y_ << std::endl;
			}

			for (int i = 0; i < (int)p_blockMesh->n_nodeNum; i++)
			{
                outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << p_blockMesh->v_node[i].z_ << std::endl;
			}

			for (int i = 0; i < (int)p_blockMesh->n_elemNum; i++)
			{
                outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << this->elementField.GetValue(i) << std::endl;
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
		}
		outFile.close();
	}
}

//write output file for a scalar field
template<>
void Field<Scalar>::WriteVTK_Field(const std::string& outMshFileName)
{
	this->ElementToNode();
    std::ofstream	outFile(outMshFileName.c_str());

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
    outFile << "SCALARS " <<this->st_name  << " float " << std::endl;
    outFile << "LOOKUP_TABLE " << "default" << std::endl;
	for (int i = 0; i < (int)p_blockMesh->v_elem.size(); i++)
	{
        outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << this->elementField.GetValue(i) << std::endl;
	}

	outFile.close();
}
template<>
void Field<Scalar>::ReadVTK_Field(const std::string& inVTKFileName)
{
	std::cout << "start read vtk mesh data" << std::endl;
	/****read vtk file mesh********/
	std::ifstream inFile(inVTKFileName);

	std::string comment;
	getline(inFile, comment);


	std::string WriterName;
	getline(inFile, WriterName);


	std::string codingFormat;
	inFile >> codingFormat;

	std::string DataSet;
	inFile >> DataSet;
	inFile >> DataSet;

	while (true)
	{
		std::string dataComment;
		inFile >> dataComment;

		if (dataComment == "POINTS")
		{
			int pointNum;
			std::string  pointNumType;

			inFile >> pointNum >> pointNumType;

			for (size_t i = 0; i < pointNum; i++)
			{
				Scalar point_x, point_y, point_z;
				inFile >> point_x >> point_y >> point_z;
			}
		}
		else if (dataComment == "CELLS")
		{
			int nCellNum, nCellTotalNum;
			inFile >> nCellNum >> nCellTotalNum;

			for (size_t i = 0; i < nCellNum; i++)
			{
				int nCellNodeNum, nNodeID;
				inFile >> nCellNodeNum;
				for (size_t j = 0; j < nCellNodeNum; j++)
				{
					inFile >> nNodeID;
				}
			}
		}
		else if (dataComment == "CELL_TYPES")
		{
			int nCellTypeNum, nCellType;
			inFile >> nCellTypeNum;
			for (size_t i = 0; i < nCellTypeNum; i++)
			{
				inFile >> nCellType;
			}
			break;
		}
		else
		{
			break;
		}

	}
	/** **read field data*************************/
	ReadVTKGridField(inFile);
}


template<>
void Field<Scalar>::ReadVTKGridField(std::ifstream& inFile)
{
	std::cout<<"start reading field " << this->st_name << " from vtk file" << std::endl;
	std::string comment;
	inFile >> comment;

	int DataNumber;
	inFile >> DataNumber;

	while(!inFile.eof())
	{
		std::string sDataType, sDataName, sFormatType;
		inFile >> sDataType >> sDataName >> sFormatType;
		std::string sLookuptable, sLookuptableModel;
		inFile >> sLookuptable >> sLookuptableModel;
		Scalar nData;
		if (sDataType == "SCALARS" && sDataName == this->st_name)
		{
			for (size_t i = 0; i < DataNumber; i++)
			{
				inFile >> nData;
				elementField.SetValue(i, nData);
			}
			inFile.close();
			return;
		}
		else
		{
			for (size_t i = 0; i < DataNumber; i++)
			{
				inFile >> nData;
			}
		}
	}
	FatalError("this vtk file has no field named " + this->st_name);
	return;
}


// =======================================Vector Type Field=============================================

//creating boundary fields according to mesh information
template<>
void Field<Vector>::CreateBoundaryField()
{
    for (int i = 0; i < (int)p_blockMesh->v_boundaryFaceZone.size(); i++)
    {
        FaceZone* pfz = &(p_blockMesh->v_boundaryFaceZone[i]);
        this->v_BoundaryFaceField.push_back(BoundaryFaceZoneField<Vector>(p_blockMesh, pfz));
    }
}

//constructor
template<>
Field<Vector>::Field
(
Mesh* UGBPointer,
const std::string& fieldName
)
:
p_blockMesh(UGBPointer),
st_name(fieldName),
elementField(UGBPointer),
elementField0(UGBPointer),
faceField(UGBPointer),
nodeField(UGBPointer),
interiorFaceField(UGBPointer)
{
	this->CreateBoundaryField();
}

//constructor
template<>
Field<Vector>::Field
(
Mesh* UGBPointer,
Vector fieldValue,
const std::string& fieldName
)
:
p_blockMesh(UGBPointer),
st_name(fieldName),
elementField(UGBPointer),
elementField0(UGBPointer),
faceField(UGBPointer),
nodeField(UGBPointer),
interiorFaceField(UGBPointer)
{
	elementField.Initialize(fieldValue);
	this->CreateBoundaryField();
}

//constructor
template<>
Field<Vector>::Field
(
Mesh* UGBPointer,
const std::vector<Vector>& valueList,
const std::string& fieldName
)
:
p_blockMesh(UGBPointer),
st_name(fieldName),
elementField(UGBPointer),
elementField0(UGBPointer),
faceField(UGBPointer),
nodeField(UGBPointer),
interiorFaceField(UGBPointer)
{
	elementField.Initialize(valueList);
	this->CreateBoundaryField();
}

//constructor
template<>
Field<Vector>::Field
(
Mesh* UGBPointer,
Vector(udf)(Scalar, Scalar, Scalar),
const std::string& fieldName
)
:
p_blockMesh(UGBPointer),
st_name(fieldName),
elementField(UGBPointer),
elementField0(UGBPointer),
faceField(UGBPointer),
nodeField(UGBPointer),
interiorFaceField(UGBPointer)
{
	elementField.Initialize(udf);
	this->CreateBoundaryField();
}

//specify a special boundary condition on a boundary face zone
template<>
void Field<Vector>::SetBoundaryCondition
(
const std::string& patchName,
const std::string& boundaryType
)
{
	int faceZoneID = this->p_blockMesh->GetBoundaryFaceZoneID(patchName);
	if (-1 != faceZoneID)
	{
		if ("boundaryLeft" == boundaryType)
		{
			this->v_BoundaryFaceField[faceZoneID].SetSpecialBC(BoundaryFaceZoneField<Vector>::fbtCaclulated);
            std::cout << "boundaryLeft condition successfully given on boundary: " << patchName << " of " << this->st_name << std::endl;
		}
		else if ("copy" == boundaryType)
		{
			this->v_BoundaryFaceField[faceZoneID].SetSpecialBC(BoundaryFaceZoneField<Vector>::fbtCopy);
            std::cout << "Copy condition successfully given on boundary: " << patchName << " of " << this->st_name << std::endl;
		}
		else if ("extrapolate" == boundaryType)
		{
			this->v_BoundaryFaceField[faceZoneID].SetSpecialBC(BoundaryFaceZoneField<Vector>::fbtExtrapolate);
            std::cout << "Extrapolate condition successfully given on boundary: " << patchName << " of " << this->st_name << std::endl;
		}
		else
		{
			WarningContinue("unknown boundary type:" + boundaryType);
		}
	}
	return;
}

//specify a uniform general-type (ABC-form) boundary condition on a boundary face zone
template<>
void Field<Vector>::SetBoundaryCondition
(
const std::string& patchName,
Scalar bc_a,
Scalar bc_b,
Vector bc_c
)
{
	int faceZoneID = this->p_blockMesh->GetBoundaryFaceZoneID(patchName);
	if (-1 != faceZoneID)
	{
		this->v_BoundaryFaceField[faceZoneID].SetNumerBC(bc_a, bc_b, bc_c);
        std::cout << "Boundary condition successfully given on boundary: " << patchName << " of " << this->st_name << std::endl;
	}
	return;
}

//specify a non-uniform general-type (ABC-form) boundary condition on a boundary face zone
template<>
void Field<Vector>::SetBoundaryCondition
(
const std::string& patchName,
Scalar bc_a,
Scalar bc_b,
Vector(*udf)(Scalar x, Scalar y, Scalar z)
)
{
	int faceZoneID = this->p_blockMesh->GetBoundaryFaceZoneID(patchName);
	if (-1 != faceZoneID)
	{
		this->v_BoundaryFaceField[faceZoneID].SetNumerBC(bc_a, bc_b, udf);
        std::cout << "Boundary condition successfully given on boundary: " << patchName << " of " << this->st_name << std::endl;
	}
	return;
}

//specify a uniform general-type (ABC-form) boundary condition on a boundary face zone
template<>
void Field<Vector>::SetBoundaryCondition
(
const std::string& patchName,
Scalar bc_a,
Scalar bc_b,
std::vector<Vector > bc_c_list
)
{
	int faceZoneID = this->p_blockMesh->GetBoundaryFaceZoneID(patchName);
	if (-1 != faceZoneID)
	{
        std::vector<NumerBC<Vector> > abcList;
		for (int i = 0; i < (int)bc_c_list.size(); i++)
		{
			abcList.push_back(NumerBC<Vector>(bc_a, bc_b, bc_c_list[i]));
		}
		this->v_BoundaryFaceField[faceZoneID].SetNumerBC(abcList);
        std::cout << "Boundary condition successfully given on boundary: " << patchName << std::endl;
	}
	return;
}

//check wether all boundaries have been correctly specified
template<>
void Field<Vector>::CheckBoundaryCondition() const
{
    std::vector<std::string> v_facesWithoutBC;
	//collect names of faces with no boundary condtion given on;
	for (int i = 0; i < (int)v_BoundaryFaceField.size(); i++)
	{
		if (!(v_BoundaryFaceField[i].BCExist()))
		{
			v_facesWithoutBC.push_back(v_BoundaryFaceField[i].p_faceZone->name);
		}
	}
	int NumOfFacesWithoutBC = (int)v_facesWithoutBC.size();
	//if there exist face without boundary condition, display fatal error with a list of no-BC face names included.
	if (0 != NumOfFacesWithoutBC)
	{
        std::stringstream message;
        message << "For field " + this->st_name + ", boundary conidtions are not found on faces:" << std::endl;
		for (int i = 0; i < (int)v_facesWithoutBC.size(); i++)
		{
            message << "(" << i + 1 << ") " << v_facesWithoutBC[i] << std::endl;
		}
		FatalError(message.str());
	}
}

//copy boundary conditions to another field
template<>
void Field<Vector>::CopyBoundaryConditionTo(Field<Vector>& targetField)
{
	if (targetField.p_blockMesh != this->p_blockMesh)
	{
        std::stringstream message;
		message << "The host field " << this->st_name << " and the target one ";
		message << targetField.st_name << " are not built on the same mesh.";
		FatalError(message.str());
	}
    for (int i = 0; i < (int)this->v_BoundaryFaceField.size(); i++)
	{
		this->v_BoundaryFaceField[i].CopyBoundaryConditionTo(targetField.v_BoundaryFaceField[i]);
	}
}

//copying elementField to elementField0
template<>
void Field<Vector>::SaveOld()
{
	if (false == this->elementField.Assignment())
	{
		FatalError("Cannot save old for not passing Assignment Check");
	}
	this->elementField0.v_value.assign(this->elementField.v_value.begin(), this->elementField.v_value.end());
	this->elementField0.fs_status = BaseField<Vector>::fsAssigned;
}

//calculate face values on boundary
template<>
void Field<Vector>::UpdateBoundary(const std::string& option)
{
	if (false == this->faceField.Assignment() || false == this->elementField.Assignment())
	{
		FatalError("Cannot UpdateBoundary for not passing Assignment Check, element- and face-field both need be available");
		return;
	}
	if ("boundaryLeft" == option)
	{
		for (int i = 0; i < (int)v_BoundaryFaceField.size(); i++)
		{
			v_BoundaryFaceField[i].KeepCurrentValue(faceField);
		}
	}
	else if ("copy" == option)
	{
		for (int i = 0; i < (int)v_BoundaryFaceField.size(); i++)
		{
			v_BoundaryFaceField[i].CopyFromElement(elementField, faceField);
		}
	}
	else
	{
		//Step 1: set initial value for at boundaries
		for (int i = 0; i < (int)v_BoundaryFaceField.size(); i++)
		{
			if (BoundaryFaceZoneField<Vector>::fbtCaclulated == v_BoundaryFaceField[i].fbt_BCType)
			{
				v_BoundaryFaceField[i].KeepCurrentValue(faceField);
			}
			else
			{
				v_BoundaryFaceField[i].CopyFromElement(elementField, faceField);
			}
		}
		//Step 2: calculate the values according to the given boundary conditions
		for (int i = 0; i < (int)v_BoundaryFaceField.size(); i++)
		{
			v_BoundaryFaceField[i].CalculateValue(elementField, faceField, option);
		}
	}
	return;
}

void ElementToFaceWithAllParameters
(
Field<Vector>& phi,
const std::string& interiorApproach,
const std::string& boundaryApproach
)
{
	//step 3: interpolate values on interior faces
	phi.interiorFaceField.CalculateValue(phi.elementField, phi.faceField, interiorApproach);
	//step 4: calculate values on boundary according to selected method
	phi.UpdateBoundary(boundaryApproach);
}

//interpolate from element to face. On boundaris, values are calculated by the given option
template<>
void Field<Vector>::ElementToFace(const std::string& approach)
{
	ElementToFaceCheck(this->elementField, this->faceField);
    std::string interiorParameter = "linear";
    std::string boundaryParameter = "default";
	int parameter = ElementToFaceParameterCheck(approach);
	if (-1 == parameter)
	{
		FatalError("illegal parameter in ElementToFace: " + approach);
	}
	else if (0 == parameter)
	{
		interiorParameter = approach;
	}
	else if (1 == parameter)
	{
		boundaryParameter = approach;
	}
	ElementToFaceWithAllParameters(*this, interiorParameter, boundaryParameter);
	return;
}


//interpolate from element to face. On boundaris, values are calculated by the given option
template<>
void Field<Vector>::ElementToFace(const std::string& option1, const std::string& option2)
{
	ElementToFaceCheck(this->elementField, this->faceField);
    std::string interiorParameter = "linear";
    std::string boundaryParameter = "default";
    std::vector<std::string> v_option;
	v_option.push_back(option1);
	v_option.push_back(option2);
	for (int i = 0; i < (int)v_option.size(); i++)
	{
		int parameter = ElementToFaceParameterCheck(v_option[i]);
		if (-1 == parameter)
		{
			FatalError("illegal parameter in ElementToFace: " + v_option[i]);
		}
		else if (0 == parameter)
		{
			interiorParameter = v_option[i];
		}
		else if (1 == parameter)
		{
			boundaryParameter = v_option[i];
		}
	}
	ElementToFaceWithAllParameters(*this, interiorParameter, boundaryParameter);
	return;
}

template<>
void Field<Vector>::ElementToNode()
{
	//step 1: check whether the element field has values on it
	if (false == this->elementField.Assignment())
	{
		FatalError("Cannot proceed ElementToNode for not passing Assignment Check");
	}
	//step 2: create node field if not existing
	if (BaseField<Vector>::fsNotExist == nodeField.fs_status)
	{
		nodeField.Initialize();
	}
	//step 3: visit all nodes in mesh data
	for (int i = 0; i < (int)p_blockMesh->v_vertice.size(); i++)
	{
        std::vector<Node> relativePos;
        std::vector<Vector> value;
        std::vector<bool> isBoundary;
		//step 3.1: collect topological information from element field
		for (int j = 0; j < (int)p_blockMesh->v_vertice[i].v_elemID.size(); j++)
		{
			int elemID = p_blockMesh->v_vertice[i].v_elemID[j];
			Node elementCenter = p_blockMesh->v_elem[elemID].center;
			relativePos.push_back(elementCenter - p_blockMesh->v_node[i]);
			value.push_back(elementField.GetValue(elemID));
			//an element center is of course not on the boundary; 
			isBoundary.push_back(false);
		}
		//step 3.2: if existing, collect toplogical information from boundary face field
		if (BaseField<Vector>::fsAssigned == faceField.fs_status)
		{
			for (int j = 0; j < (int)p_blockMesh->v_vertice[i].v_faceID.size(); j++)
			{
				int faceID = p_blockMesh->v_vertice[i].v_faceID[j];
				//check wether this face is on a boundray and push back values correspondingly;
				if (-1 == p_blockMesh->SearchOwnerNeighbor(faceID).second)
				{
					Node faceCenter = p_blockMesh->v_face[faceID].center;
					relativePos.push_back(faceCenter - p_blockMesh->v_node[i]);
					value.push_back(faceField.GetValue(faceID));
					isBoundary.push_back(true);
				}
			}
		}
		//step 3.3: check wether we have interplating points located on boundary.
		//If having, delete those interpolating points at interior;
		bool hasBoundary = false;
		for (int j = 0; j < (int)isBoundary.size(); j++)
		{
			if (true == isBoundary[j])
			{
				hasBoundary = true;
				break;
			}
		}
		if (true == hasBoundary)
		{
            std::vector<Node> relativePos_backup;
            std::vector<Vector> value_backup;
			for (int j = 0; j < (int)isBoundary.size(); j++)
			{
				relativePos_backup.push_back(relativePos[j]);
				value_backup.push_back(value[j]);
			}
			relativePos.clear();
			value.clear();
			for (int j = 0; j < (int)isBoundary.size(); j++)
			{
				if (true == isBoundary[j])
				{
					relativePos.push_back(relativePos_backup[j]);
					value.push_back(value_backup[j]);
				}
			}
		}
		//step 3.4: use interpolation tools and write values
		Vector interpolatedValue = Interpolation(relativePos, value);
		nodeField.SetValue(i, interpolatedValue);
	}

	nodeField.fs_status = BaseField<Vector>::fsAssigned;
}

template<>
void Field<Vector>::NodeToElement()
{
	//step 1: check whether the element field has values on it
	if (false == this->nodeField.Assignment())
	{
		FatalError("Cannot proceed NodeToElement for not passing Assignment Check");
	}
	//step 2: create elementfield if not existing
	if (BaseField<Vector>::fsNotExist == nodeField.fs_status)
	{
		elementField.Initialize();
	}
	//step 3: visit all nodes in mesh data
	for (int i = 0; i < p_blockMesh->n_elemNum; i++)
	{
		std::vector<Node> relativePos;
		std::vector<Vector> value;
		Vector cellCenter = p_blockMesh->v_elem[i].center;
		for (int j = 0;j < p_blockMesh->v_elem[i].v_nodeID.size();j++)
		{
			int nodeID = p_blockMesh->v_elem[i].v_nodeID[j];
			Vector nodeValue = this->nodeField.GetValue(nodeID);
			relativePos.push_back(p_blockMesh->v_node[nodeID] - cellCenter);
			value.push_back(nodeValue);
		}
		Vector interpolatedValue = Interpolation(relativePos, value);
		elementField.SetValue(i, interpolatedValue);
	}
	elementField.fs_status = BaseField<Vector>::fsAssigned;
	return;
}

//Kong Ling supplemented on 2019/12/23
template<>
std::vector<Vector> Field<Vector>::ExtractBoundary(const std::string& boundaryName)
{
	return this->faceField.ExtractBoundary(boundaryName);
}

//Kong Ling modified on 2019/12/23
//load values on a specific boundary with a ordered value list
template<>
void Field<Vector>::LoadBoundary(const std::string& boundaryName, const std::vector<Vector>& valueList)
{
	this->faceField.LoadBoundary(boundaryName, valueList);
	return;
}

template<>
Field<Vector>& Field<Vector>::operator = (const Field<Vector>& rhs)
{
    if (true == rhs.nodeField.Assignment())
    {
        this->nodeField = rhs.nodeField;
    }
    if (true == rhs.faceField.Assignment())
    {
        this->faceField = rhs.faceField;
    }
    if (true == rhs.elementField.Assignment())
    {
        this->elementField = rhs.elementField;
    }
    return *this;
}

//write output file for a scalar field in vtk form
template<>
void Field<Vector>::WriteVTK_Field(const std::string& outMshFileName)
{

	this->ElementToNode();
	std::ofstream	outFile(outMshFileName.c_str());

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
	outFile << "POINT_DATA " << p_blockMesh->n_nodeNum << std::endl;
	outFile << "VECTORS vectors " << "float "  << std::endl;
	//outFile << "LOOKUP_TABLE " << "default" << std::endl;
	for (int i = 0; i < (int)p_blockMesh->n_nodeNum; i++)
	{
		outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) <<
			this->nodeField.GetValue(i).x_ << " "
		<< this->nodeField.GetValue(i).y_ <<" "
		<< this->nodeField.GetValue(i).z_ <<std::endl;
	}

	outFile.close();

}

//write output file for a scalar field
template<>
void Field<Vector>::WriteTecplotField(const std::string& outMshFileName)
{
	//Set v_nodeField number
	this->ElementToNode();

    std::ofstream	outFile(outMshFileName.c_str());

    outFile << "TITLE =\"" << this->st_name << "\"" << std::endl;

	if (p_blockMesh->md_meshDim == Mesh::md2D)
	{
        outFile << "VARIABLES = " << "\"x [m]\"," << "\"y [m]\"," << "\"u\"," << "\"v\"" << std::endl;

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
            outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << p_blockMesh->v_node[i].x_ << " " << p_blockMesh->v_node[i].y_ << " " << this->nodeField.GetValue(i).x_ << " " << this->nodeField.GetValue(i).y_ << std::endl;
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
				<< this->p_blockMesh->v_boundaryFaceZone[i].v_faceID.size() << ", ZONETYPE = FELINESEG, " << "VARSHARELIST = ([1,2,3,4]=1,[1,2])" << std::endl;

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
	else if (p_blockMesh->md_meshDim == Mesh::md3D)
	{
		if (this->p_blockMesh->est_shapeType == Element::estPolyhedron)
		{
			int nFaceNodeNum(0);
			for (int i = 0; i < (int)this->p_blockMesh->v_face.size(); i++)
			{
				for (int j = 0; j < (int)this->p_blockMesh->v_face[i].v_nodeID.size(); j++)
				{
					nFaceNodeNum++;
				}
			}

			outFile << "VARIABLES = " << std::endl << "\"X\"" << std::endl << "\"Y\"" << std::endl << "\"Z\"" << std::endl << "\"u\"," << "\"v\"," << "\"w\"" << std::endl;
			outFile << "DATASETAUXDATA Common.VectorVarsAreVelocity = \"TRUE\"" << std::endl;
			outFile << "ZONE T = \"fluid\"" << std::endl;
			outFile << "STRANDID=0, SOLUTIONTIME=0" << std::endl;
			outFile << "Nodes=" << this->p_blockMesh->n_nodeNum << ", Faces=" << this->p_blockMesh->n_faceNum << ", Elements=" << this->p_blockMesh->n_elemNum << ", ZONETYPE=FEPolyhedron" << std::endl;
			outFile << "DATAPACKING=BLOCK" << std::endl;
			outFile << "TotalNumFaceNodes=" << nFaceNodeNum << ", NumConnectedBoundaryFaces=0, TotalNumBoundaryConnections=0" << std::endl;
			outFile << "DT=(SINGLE SINGLE SINGLE )";

			// Write node X
			for (int i = 0; i < (int) this->p_blockMesh->v_node.size(); i++)
			{
				if (i % 5 == 0) outFile << std::endl;
                outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << this->p_blockMesh->v_node[i].x_ << " ";
			}

			// Write node Y
			for (int i = 0; i < (int) this->p_blockMesh->v_node.size(); i++)
			{
				if (i % 5 == 0) outFile << std::endl;
                outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << this->p_blockMesh->v_node[i].y_ << " ";
			}

			// Write node Z
			for (int i = 0; i < (int) this->p_blockMesh->v_node.size(); i++)
			{
				if (i % 5 == 0) outFile << std::endl;
                outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << this->p_blockMesh->v_node[i].z_ << " ";
			}

			// Write node u
			for (int i = 0; i < (int) this->p_blockMesh->v_node.size(); i++)
			{
				if (i % 5 == 0) outFile << std::endl;
                outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << this->nodeField.v_value[i].x_ << " ";
			}

			// Write node v
			for (int i = 0; i < (int) this->p_blockMesh->v_node.size(); i++)
			{
				if (i % 5 == 0) outFile << std::endl;
                outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << this->nodeField.v_value[i].y_ << " ";
			}

			// Write node w
			for (int i = 0; i < (int) this->p_blockMesh->v_node.size(); i++)
			{
				if (i % 5 == 0) outFile << std::endl;
                outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << this->nodeField.v_value[i].z_ << " ";
			}

			outFile << std::endl << "# node count per face";
			for (int i = 0; i < (int) this->p_blockMesh->v_face.size(); i++)
			{
				if (i % 10 == 0) outFile << std::endl;
				outFile << this->p_blockMesh->v_face[i].v_nodeID.size() << " ";
			}

			int nFaceNodeCheck(0);
			outFile << std::endl << "# face nodes";
			for (int i = 0; i < (int) this->p_blockMesh->v_face.size(); i++)
			{
				for (int j = 0; j < (int)this->p_blockMesh->v_face[i].v_nodeID.size(); j++)
				{
					if (nFaceNodeCheck % 10 == 0) outFile << std::endl;
					outFile << this->p_blockMesh->v_face[i].v_nodeID[j] + 1 << " ";
					nFaceNodeCheck++;
				}
			}

			outFile << std::endl << "# left elements";
			for (int i = 0; i < (int) this->p_blockMesh->v_face.size(); i++)
			{
				if (i % 10 == 0) outFile << std::endl;
				outFile << this->p_blockMesh->v_face[i].n_owner + 1 << " ";
			}

			outFile << std::endl << "# node count per face";
			for (int i = 0; i < (int) this->p_blockMesh->v_face.size(); i++)
			{
				if (i % 10 == 0) outFile << std::endl;
				outFile << this->p_blockMesh->v_face[i].n_neighbor + 1 << " ";
			}
		}
		else
		{
            outFile << "VARIABLES = " << "\"x [m]\"," << "\"y [m]\"," << "\"z [m]\"" << "\"u\"," << "\"v\"," << "\"w\"" << std::endl;

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
                outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << p_blockMesh->v_node[i].x_ << " " << p_blockMesh->v_node[i].y_ << " " << p_blockMesh->v_node[i].z_ << " " << this->nodeField.GetValue(i).x_ << " " << this->nodeField.GetValue(i).y_ << " " << this->nodeField.GetValue(i).z_ << std::endl;
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
					<< this->p_blockMesh->v_boundaryFaceZone[i].v_faceID.size() << ", ZONETYPE = FEQUADRILATERAL, " << "VARSHARELIST = ([1,2,3,4,5,6]=1,[1,2,3])" << std::endl;

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

	//system(outMshFileName.c_str());
}

//Writing BC Tecplot file for visualization
template<>
void Field<Vector>::WriteTecplotField(const std::string& patchName, const std::string& outMshFileName)
{
	//Set v_nodeField number
	this->ElementToNode();

	int found = -1;
	for (int i = 0; i < (int)this->v_BoundaryFaceField.size(); i++)
	{
		if (patchName == this->v_BoundaryFaceField[i].p_faceZone->name)
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
		//Face ID list in current BC
        std::vector<int>& vFaceID = this->v_BoundaryFaceField[found].p_faceZone->v_faceID;
        std::vector<Face> vFaceOutput;
		vFaceOutput.resize((int)vFaceID.size());

        std::vector<Node> vNodeOutput;
        std::vector<Vector> vNodeValue;
        std::map<int, int> NodeID_LocalID;
        std::map<int, int>::iterator iter;

		for (int i = 0; i < (int)vFaceID.size(); i++)
		{
			int curFaceID = vFaceID[i];
			Face& curFace = this->p_blockMesh->v_face[curFaceID];

			vFaceOutput[i].ft_faceType = curFace.ft_faceType;
			vFaceOutput[i].v_nodeID.resize((int)curFace.v_nodeID.size());

			for (int j = 0; j < (int)curFace.v_nodeID.size(); j++)
			{
				int curNodeID = curFace.v_nodeID[j];

				iter = NodeID_LocalID.begin();
				iter = NodeID_LocalID.find(curNodeID);

				//find global Node ID in map::NodeID_LocalID
				if (iter != NodeID_LocalID.end())
				{
					vFaceOutput[i].v_nodeID[j] = iter->second;
				}
				else
				{
					vFaceOutput[i].v_nodeID[j] = (int)vNodeOutput.size();

					vNodeOutput.push_back(this->p_blockMesh->v_node[curNodeID]);
					vNodeValue.push_back(this->nodeField.GetValue(curNodeID));
					NodeID_LocalID.insert(std::pair<int, int>(curNodeID, vFaceOutput[i].v_nodeID[j]));
				}
			}
		}

        std::ofstream	outFile(outMshFileName.c_str());

		if (this->p_blockMesh->md_meshDim == Mesh::md2D)
		{
            outFile << "TITLE =\"" << this->v_BoundaryFaceField[found].p_faceZone->name << "\"" << std::endl;
            outFile << "VARIABLES = " << "\"x [m]\"," << "\"y [m]\"," << "\"u\"," << "\"v\"" << std::endl;

			outFile << "ZONE T = " << "\"" << this->v_BoundaryFaceField[found].p_faceZone->name << "\"," << " DATAPACKING = POINT, N = " << vNodeOutput.size() << ", E = "
				<< vFaceOutput.size() << ", ZONETYPE = FELINESEG, " << std::endl;

			for (unsigned int i = 0; i < vNodeOutput.size(); i++)
			{
                outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << vNodeOutput[i].x_ << " " << vNodeOutput[i].y_ << " " << vNodeValue[i].x_ << " " << vNodeValue[i].y_ << std::endl;
			}

			for (int i = 0; i < (int)vFaceOutput.size(); i++)
			{
				for (int n = 0; n < (int)vFaceOutput[i].v_nodeID.size(); n++)
				{
					outFile << vFaceOutput[i].v_nodeID[n] + 1 << "\t";
				}

				outFile << std::endl;
			}
		}
		else if (this->p_blockMesh->md_meshDim == Mesh::md3D)
		{
            outFile << "TITLE =\"" << this->v_BoundaryFaceField[found].p_faceZone->name << "\"" << std::endl;
            outFile << "VARIABLES = " << "\"x [m]\"," << "\"y [m]\"," << "\"z [m]\"," << "\"u\"," << "\"v\"," << "\"w\"" << std::endl;

			outFile << "ZONE T = " << "\"" << this->v_BoundaryFaceField[found].p_faceZone->name << "\"," << " DATAPACKING = POINT, N = " << vNodeOutput.size() << ", E = "
				<< vFaceOutput.size() << ", ZONETYPE = FEQUADRILATERAL, " << std::endl;

			for (unsigned int i = 0; i < vNodeOutput.size(); i++)
			{
                outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << vNodeOutput[i].x_ << " " << vNodeOutput[i].y_ << " " << vNodeOutput[i].z_ << " "
                    << vNodeValue[i].x_ << " " << vNodeValue[i].y_ << " " << vNodeValue[i].z_ << std::endl;
			}

			for (int i = 0; i < (int)vFaceOutput.size(); i++)
			{
				for (int n = 0; n < (int)vFaceOutput[i].v_nodeID.size(); n++)
				{
					outFile << vFaceOutput[i].v_nodeID[n] + 1 << "\t";
				}

				if (vFaceOutput[i].ft_faceType == Face::ftTrangular)
				{
					outFile << vFaceOutput[i].v_nodeID[(int)vFaceOutput[i].v_nodeID.size() - 1] + 1 << "\t";
				}
				outFile << std::endl;
			}
		}

		outFile.close();
	}
}

//write output file for a scalar field
template<>
void Field<Vector>::WriteTecplotField(const std::string& outMshFileName, outType outFileType)
{
	if (outFileType == otCellNode)
	{
		this->WriteTecplotField(outMshFileName.c_str());
	}
	else if (outFileType == otCellCenter)
	{
        std::ofstream	outFile(outMshFileName.c_str());

        outFile << "TITLE =\"" << this->st_name << "\"" << std::endl;

		if (p_blockMesh->md_meshDim == Mesh::md2D)
		{
            outFile << "VARIABLES = " << "\"x [m]\"," << "\"y [m]\"," << "\"u\"," << "\"v\"" << std::endl;

			if (p_blockMesh->est_shapeType == Element::estTriangular)
			{
				outFile << "ZONE T = " << "\"" << "Tri_mesh" << "\"," << " DATAPACKING = BLOCK, N = " << p_blockMesh->n_nodeNum << ", E = " << p_blockMesh->n_elemNum
                    << ", ZONETYPE = FETRIANGLE" << std::endl;
                outFile << "VARLOCATION = ([1-2]=NODAL, [3-4] = CELLCENTERED)" << std::endl;
			}
			else if (p_blockMesh->est_shapeType == Element::estQuadrilateral)
			{
				outFile << "ZONE T = " << "\"" << "Quad_mesh" << "\"," << " DATAPACKING = BLOCK, N = " << p_blockMesh->n_nodeNum << ", E = " << p_blockMesh->n_elemNum
                    << ", ZONETYPE = FEQUADRILATERAL" << std::endl;
                outFile << "VARLOCATION = ([1-2]=NODAL, [3-4] = CELLCENTERED)" << std::endl;
			}
			else if (p_blockMesh->est_shapeType == Element::estMixed)
			{
				outFile << "ZONE T = " << "\"" << "Mix_mesh" << "\"," << " DATAPACKING = BLOCK, N = " << p_blockMesh->n_nodeNum << ", E = " << p_blockMesh->n_elemNum
                    << ", ZONETYPE = FEQUADRILATERAL" << std::endl;
                outFile << "VARLOCATION = ([1-2]=NODAL, [3-4] = CELLCENTERED)" << std::endl;
			}

			for (int i = 0; i < (int)p_blockMesh->n_nodeNum; i++)
			{
                outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << p_blockMesh->v_node[i].x_ << std::endl;
			}

			for (int i = 0; i < (int)p_blockMesh->n_nodeNum; i++)
			{
                outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << p_blockMesh->v_node[i].y_ << std::endl;
			}

			for (int i = 0; i < (int)p_blockMesh->n_elemNum; i++)
			{
                outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << this->elementField.GetValue(i).x_ << std::endl;
			}

			for (int i = 0; i < (int)p_blockMesh->n_elemNum; i++)
			{
                outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << this->elementField.GetValue(i).y_ << std::endl;
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
		}
		else if (p_blockMesh->md_meshDim == Mesh::md3D)
		{
            outFile << "VARIABLES = " << "\"x [m]\"," << "\"y [m]\"," << "\"z [m]\"" << "\"u\"," << "\"v\"," << "\"w\"" << std::endl;

			if (p_blockMesh->est_shapeType == Element::estTetrahedral)
			{
				outFile << "ZONE T = " << "\"" << "Tetrahedral_mesh" << "\"," << " DATAPACKING = BLOCK, N = " << p_blockMesh->n_nodeNum << ", E = " << p_blockMesh->n_elemNum
                    << ", ZONETYPE = FETETRAHEDRON" << std::endl;
                outFile << "VARLOCATION = ([1-3]=NODAL, [4-6] = CELLCENTERED)" << std::endl;
			}
			else if (p_blockMesh->est_shapeType == Element::estHexahedral)
			{
				outFile << "ZONE T = " << "\"" << "Hexahedral_mesh" << "\"," << " DATAPACKING = BLOCK, N = " << p_blockMesh->n_nodeNum << ", E = " << p_blockMesh->n_elemNum
                    << ", ZONETYPE = FEBRICK" << std::endl;
                outFile << "VARLOCATION = ([1-3]=NODAL, [4-6] = CELLCENTERED)" << std::endl;
			}
			else if (p_blockMesh->est_shapeType == Element::estMixed)
			{
				outFile << "ZONE T = " << "\"" << "Tet_mesh" << "\"," << " DATAPACKING = BLOCK, N = " << p_blockMesh->n_nodeNum << ", E = " << p_blockMesh->n_elemNum
                    << ", ZONETYPE = FEBRICK" << std::endl;
                outFile << "VARLOCATION = ([1-3]=NODAL, [4-6] = CELLCENTERED)" << std::endl;
			}
			else if (p_blockMesh->est_shapeType == Element::estWedge)
			{
				outFile << "ZONE T = " << "\"" << "Wedge_mesh" << "\"," << " DATAPACKING = BLOCK, N = " << p_blockMesh->n_nodeNum << ", E = " << p_blockMesh->n_elemNum
                    << ", ZONETYPE = FEBRICK" << std::endl;
                outFile << "VARLOCATION = ([1-3]=NODAL, [4-6] = CELLCENTERED)" << std::endl;
			}

			for (int i = 0; i < (int)p_blockMesh->n_nodeNum; i++)
			{
                outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << p_blockMesh->v_node[i].x_ << std::endl;
			}

			for (int i = 0; i < (int)p_blockMesh->n_nodeNum; i++)
			{
                outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << p_blockMesh->v_node[i].y_ << std::endl;
			}

			for (int i = 0; i < (int)p_blockMesh->n_nodeNum; i++)
			{
                outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << p_blockMesh->v_node[i].z_ << std::endl;
			}

			for (int i = 0; i < (int)p_blockMesh->n_elemNum; i++)
			{
                outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << this->elementField.GetValue(i).x_ << std::endl;
			}

			for (int i = 0; i < (int)p_blockMesh->n_elemNum; i++)
			{
                outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << this->elementField.GetValue(i).y_ << std::endl;
			}

			for (int i = 0; i < (int)p_blockMesh->n_elemNum; i++)
			{
                outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << this->elementField.GetValue(i).z_ << std::endl;
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
		}
		outFile.close();
	}
}



template<>
void Field<Vector>::ReadVTK_Field(const std::string& inVTKFileName)
{
	
}

template<>
void Field<Vector>::ReadVTKGridField(std::ifstream& inFile)
{

}