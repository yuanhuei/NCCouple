/*---------------------------------------------------------------------------*\
File Name:
	BoundaryFaceZoneField.cpp

Description:
	defining classes BoundaryFaceZoneField

	Author:		Kong Ling
	Date: 2017-12-10

	
Revised:
	Description:
	1. Convert "non const" input parameters to const;
	2. Convert "non const" function to const function;
	
	Revisor:		Shuai Zhang
	Modified Date:	2018-11-28
\*---------------------------------------------------------------------------*/

#include "../MHT_field/BoundaryFaceZoneField.h"
#include "../MHT_common/SystemControl.h"
#include "../MHT_common/InterpolationTools.h"

#include "../MHT_mesh/Mesh.h"
#include "../MHT_mesh/MeshSupport.h"
#include "../MHT_mesh/LagrangeTools.h"

#include "../MHT_field/FaceField.h"
#include "../MHT_field/ElementField.h"
#include "../MHT_field/NumerBC.h"

// =======================================Scalar Type NodeField=============================================
//constructor
template<>
BoundaryFaceZoneField<Scalar>::BoundaryFaceZoneField(Mesh* UGB, FaceZone* FZ)
:p_blockMesh(UGB),
	p_faceZone(FZ),
	fbt_BCType(fbtNotSet)
{
	this->v_NumerBC.resize(this->p_faceZone->v_faceID.size());
}

//return true if boundary condition is imposed otherwise return false
template<>
bool BoundaryFaceZoneField<Scalar>::BCExist() const
{
	if (fbtNotSet == this->fbt_BCType)
	{
		return false;
	}
	else
	{
		return true;
	}
}

//set an uniform numerical boundary condition
template<>
void BoundaryFaceZoneField<Scalar>::SetNumerBC(Scalar a_bc, Scalar b_bc, Scalar c_bc)
{
	for (int i = 0; i < (int)this->p_faceZone->v_faceID.size(); i++)
	{
		this->v_NumerBC[i] = NumerBC<Scalar>(a_bc, b_bc, c_bc);
	}
	this->fbt_BCType = fbtGeneral;
	return;
}

//set an uniform numerical boundary condition
template<>
void BoundaryFaceZoneField<Scalar>::SetNumerBC(Scalar a_bc, Scalar b_bc, Scalar(*udf)(Scalar,Scalar,Scalar))
{
	for (int i = 0; i < (int)this->p_faceZone->v_faceID.size(); i++)
	{
		int faceID = this->p_faceZone->v_faceID[i];
		Vector faceCenter = this->p_blockMesh->v_face[faceID].center;
		this->v_NumerBC[i] = NumerBC<Scalar>(a_bc, b_bc, udf(faceCenter.x_, faceCenter.y_, faceCenter.z_));
	}
	this->fbt_BCType = fbtGeneral;
	return;
}

template<>
void BoundaryFaceZoneField<Scalar>::SetNumerBC(std::vector<NumerBC<Scalar> >& givenABC)
{
	if (givenABC.size() != this->p_faceZone->v_faceID.size())
	{
		FatalError("in BoundaryFaceZoneField::SetNumerBC, the number of the given boundary condition list is incorrect");
	}
	for (int i = 0; i < (int)this->p_faceZone->v_faceID.size(); i++)
	{
		this->v_NumerBC[i] = givenABC[i];
	}
	this->fbt_BCType = fbtGeneral;
	return;
}

//set special-type boundary condition like copy and extrapolate
template<>
void BoundaryFaceZoneField<Scalar>::SetSpecialBC(FieldBoundaryType bcType)
{
	this->v_NumerBC.resize(this->p_faceZone->v_faceID.size());
	this->fbt_BCType = bcType;
	return;
}

//Copy boundary condition to a target boundary for a specific field
template<>
void BoundaryFaceZoneField<Scalar>::CopyBoundaryConditionTo(BoundaryFaceZoneField& targetBoundary)
{
	if (targetBoundary.p_faceZone != this->p_faceZone)
	{
		FatalError("in BoundaryFaceZoneField:CopyBoundaryConditionTo(), boundary pointers are not consistent.");
	}
	targetBoundary.fbt_BCType = this->fbt_BCType;
	for (int i = 0; i < (int)this->v_NumerBC.size(); i++)
	{
		targetBoundary.v_NumerBC[i] = this->v_NumerBC[i];
	}
	return;
}

//keep current value and write corresponding coefficients
template<>
void BoundaryFaceZoneField<Scalar>::KeepCurrentValue(FaceField<Scalar>& FF)
{
	for (int i = 0; i < (int)this->p_faceZone->v_faceID.size(); i++)
	{
		int faceID = this->p_faceZone->v_faceID[i];
		this->v_NumerBC[i].con1 = 0.0;
		this->v_NumerBC[i].con2 = FF.GetValue(faceID);
	}
	return;
}

//copy values from owner elements
template<>
void BoundaryFaceZoneField<Scalar>::CopyFromElement(ElementField<Scalar>& EF, FaceField<Scalar>& FF)
{
	for (int i = 0; i < (int)this->p_faceZone->v_faceID.size(); i++)
	{
		int faceID = this->p_faceZone->v_faceID[i];
		std::pair<int, int> elementIDs = this->p_blockMesh->SearchOwnerNeighbor(faceID);
		int ownerID = elementIDs.first;
		FF.SetValue(faceID, EF.GetValue(ownerID));
		this->v_NumerBC[i].con1 = 1.0;
		this->v_NumerBC[i].con2 = 0.0;
	}
}

//extrapolate values linearly from elements
template<>
void BoundaryFaceZoneField<Scalar>::ExtrapolateFromElement(ElementField<Scalar>& EF, FaceField<Scalar>& FF)
{
	Mesh* pmesh = EF.p_blockMesh;
    Vector (*LeastSquareGradient)(const std::vector<Vector>&, const std::vector<Scalar>&) = NULL;
	if (Mesh::md2D == pmesh->md_meshDim)
	{
		LeastSquareGradient = &LeastSquareGradient2D;
	}
	else if (Mesh::md3D == pmesh->md_meshDim)
	{
		LeastSquareGradient = &LeastSquareGradient3D;
	}

	for (int i = 0; i < (int)this->p_faceZone->v_faceID.size(); i++)
	{
		int faceID = (int)this->p_faceZone->v_faceID[i];
		int elemID = (this->p_blockMesh->SearchOwnerNeighbor(faceID)).first;
		Scalar phiC = EF.GetValue(elemID);
        std::vector<int> nbIDList = pmesh->SearchNearbyElements(elemID, 1);
        std::vector<Scalar> deltaPhi;
        std::vector<Vector> deltaR;
		deltaPhi.resize(nbIDList.size());
		deltaR.resize(nbIDList.size());

		for (int k = 0; k < (int)nbIDList.size(); k++)
		{
			int nbID = nbIDList[k];
			deltaPhi[k] = EF.GetValue(nbID) - phiC;
			deltaR[k] = pmesh->v_elem[nbID].center - pmesh->v_elem[elemID].center;
		}
		Vector grad = LeastSquareGradient(deltaR, deltaPhi);
		Vector CToF = pmesh->v_face[faceID].center - pmesh->v_elem[elemID].center;
		Scalar correction = CToF & grad;
		FF.SetValue(faceID, phiC + correction);
		this->v_NumerBC[i].con1 = 1.0;
		this->v_NumerBC[i].con2 = correction;
	}
}

//Calclate field values at boundary in case of general-type boundary conditon
template<>
void BoundaryFaceZoneField<Scalar>::CalculateByABC(ElementField<Scalar>& EF, FaceField<Scalar>& FF)
{
	if (fbtGeneral != this->fbt_BCType)
	{
		FatalError("Cannot calculate boundary values by CalculateByABC(), boundary type is not given as a general ABC form.");
	}
	for (int i = 0; i < (int)this->p_faceZone->v_faceID.size(); i++)
	{
		int faceID = this->p_faceZone->v_faceID[i];
		Face& patchFace = this->p_blockMesh->v_face[faceID];
		int ownerID = patchFace.n_owner;
		Element& thisElement = this->p_blockMesh->v_elem[ownerID];
        std::vector<Scalar> v_FaceValue;
        std::vector<Node> v_FaceArea;
		//collecting informations for calculation
		for (int j = 0; j < (int)thisElement.v_faceID.size(); j++)
		{
			int thisFaceID = thisElement.v_faceID[j];
			Face& thisFace = this->p_blockMesh->v_face[thisFaceID];
			if (((thisFace.center - thisElement.center) & thisFace.area) > 0)
			{
				v_FaceArea.push_back(thisFace.area);
			}
			else
			{
				v_FaceArea.push_back(-(thisFace.area));
			}
			v_FaceValue.push_back(FF.GetValue(thisFaceID));
		}
		Vector OToF = patchFace.center - thisElement.center;
		Vector efUnit = OToF.GetNormal();
		Vector nf;
		if ((efUnit&patchFace.area) > 0)
		{
			nf = patchFace.area.GetNormal();
		}
		else
		{
			nf = -patchFace.area.GetNormal();
		}
		Vector ef = efUnit / (efUnit & nf);
		Vector tf = nf - ef;
		Scalar temp = (v_FaceArea[0] & tf)*v_FaceValue[0];
		for (int j = 1; j < (int)v_FaceArea.size(); j++)
		{
			temp = temp + (v_FaceArea[j] & tf)*v_FaceValue[j];
		}
		Scalar tempScalar = v_NumerBC[i].b*ef.Mag() / OToF.Mag();
		Scalar c1 = tempScalar / (v_NumerBC[i].a + tempScalar);
		Scalar c2 = (v_NumerBC[i].c - v_NumerBC[i].b*temp / thisElement.volume) / (v_NumerBC[i].a + tempScalar);
		v_NumerBC[i].con1 = c1;
		v_NumerBC[i].con2 = c2;
		Scalar updatedValue = v_NumerBC[i].con1*EF.GetValue(ownerID) + v_NumerBC[i].con2;
		FF.SetValue(faceID, updatedValue);
	}
}

//calculating face values on the patch accroding to boundary conditions 
template<>
void BoundaryFaceZoneField<Scalar>::CalculateValue(ElementField<Scalar>& EF, FaceField<Scalar>& FF, const std::string& option)
{
	//Decide which calculation method to use from the preset boundary condtion and the given option
	FieldBoundaryType presetType = this->fbt_BCType;
	if ("boundaryLeft" == option)
	{
        presetType = fbtCaclulated;
	}
	else if ("copy" == option)
	{
        presetType = fbtCopy;
	}
	else if ("extrapolate" == option)
	{
        presetType = fbtExtrapolate;
	}
    else if ("boundaryCondition" == option && presetType != fbtGeneral)
	{
        presetType = fbtCopy;
	}
	else if ("default" != option)
	{
		FatalError("unsupported option " + option + " in BoundaryFaceZoneField::CalculateValue()");
	}
	//Select the corresponding calculation method according to the decision.
    if (fbtCaclulated == presetType)
	{
		this->KeepCurrentValue(FF);
	}
    else if (fbtCopy == presetType || fbtNotSet == presetType)
	{
		this->CopyFromElement(EF, FF);
	}
    else if (fbtExtrapolate == presetType)
	{
		this->ExtrapolateFromElement(EF, FF);
	}
    else if (fbtGeneral == presetType)
	{
		this->CalculateByABC(EF, FF);
	}
	return;
}

// =======================================Vector Type NodeField=============================================
//constructor
template<>
BoundaryFaceZoneField<Vector>::BoundaryFaceZoneField(Mesh* UGB, FaceZone* FZ)
	:p_blockMesh(UGB),
	p_faceZone(FZ),
	fbt_BCType(fbtNotSet)
{
	this->v_NumerBC.resize(this->p_faceZone->v_faceID.size());
}

//return true if boundary condition is imposed otherwise return false
template<>
bool BoundaryFaceZoneField<Vector>::BCExist() const
{
	if (fbtNotSet == this->fbt_BCType)
	{
		return false;
	}
	else
	{
		return true;
	}
}

//set an uniform numerical boundary condition
template<>
void BoundaryFaceZoneField<Vector>::SetNumerBC(Scalar a_bc, Scalar b_bc, Vector c_bc)
{
	for (int i = 0; i < (int)this->p_faceZone->v_faceID.size(); i++)
	{
		this->v_NumerBC[i] = NumerBC<Vector>(a_bc, b_bc, c_bc);
	}
	this->fbt_BCType = fbtGeneral;
	return;
}

//set an uniform numerical boundary condition
template<>
void BoundaryFaceZoneField<Vector>::SetNumerBC(Scalar a_bc, Scalar b_bc, Vector(*udf)(Scalar, Scalar, Scalar))
{
	for (int i = 0; i < (int)this->p_faceZone->v_faceID.size(); i++)
	{
		int faceID = this->p_faceZone->v_faceID[i];
		Vector faceCenter = this->p_blockMesh->v_face[faceID].center;
		this->v_NumerBC[i] = NumerBC<Vector>(a_bc, b_bc, udf(faceCenter.x_, faceCenter.y_, faceCenter.z_));
	}
	this->fbt_BCType = fbtGeneral;
	return;
}

template<>
void BoundaryFaceZoneField<Vector>::SetNumerBC(std::vector<NumerBC<Vector> >& givenABC)
{
	if (givenABC.size() != this->p_faceZone->v_faceID.size())
	{
		FatalError("in BoundaryFaceZoneField::SetNumerBC, the number of the given boundary condition list is incorrect");
	}
	for (int i = 0; i < (int)this->p_faceZone->v_faceID.size(); i++)
	{
		this->v_NumerBC[i] = givenABC[i];
	}
	this->fbt_BCType = fbtGeneral;
	return;
}

//set special-type boundary condition like copy and extrapolate
template<>
void BoundaryFaceZoneField<Vector>::SetSpecialBC(FieldBoundaryType bcType)
{
	this->v_NumerBC.resize(this->p_faceZone->v_faceID.size());
	this->fbt_BCType = bcType;
	return;
}

//Copy boundary condition to a target boundary for a specific field
template<>
void BoundaryFaceZoneField<Vector>::CopyBoundaryConditionTo(BoundaryFaceZoneField& targetBoundary)
{
	if (targetBoundary.p_faceZone != this->p_faceZone)
	{
		FatalError("in BoundaryFaceZoneField:CopyBoundaryConditionTo(), boundary pointers are not consistent.");
	}
	targetBoundary.fbt_BCType = this->fbt_BCType;
	for (int i = 0; i < (int)this->v_NumerBC.size(); i++)
	{
		targetBoundary.v_NumerBC[i] = this->v_NumerBC[i];
	}
	return;
}

//keep current value and write corresponding coefficients
template<>
void BoundaryFaceZoneField<Vector>::KeepCurrentValue(FaceField<Vector>& FF)
{
	for (int i = 0; i < (int)this->p_faceZone->v_faceID.size(); i++)
	{
		int faceID = this->p_faceZone->v_faceID[i];
		this->v_NumerBC[i].con1 = 0.0;
		this->v_NumerBC[i].con2 = FF.GetValue(faceID);
	}
	return;
}

//copy values from owner elements
template<>
void BoundaryFaceZoneField<Vector>::CopyFromElement(ElementField<Vector>& EF, FaceField<Vector>& FF)
{
	for (int i = 0; i < (int)this->p_faceZone->v_faceID.size(); i++)
	{
		int faceID = this->p_faceZone->v_faceID[i];
		std::pair<int, int> elementIDs = this->p_blockMesh->SearchOwnerNeighbor(faceID);
		int ownerID = elementIDs.first;
		FF.SetValue(faceID, EF.GetValue(ownerID));
		this->v_NumerBC[i].con1 = 1.0;
		this->v_NumerBC[i].con2 = Vector(0.0, 0.0, 0.0);
	}
}

//extrapolate values linearly from elements
template<>
void BoundaryFaceZoneField<Vector>::ExtrapolateFromElement(ElementField<Vector>& EF, FaceField<Vector>& FF)
{
	Mesh* pmesh = EF.p_blockMesh;
    Tensor(*LeastSquareGradient)(const std::vector<Vector>&, const std::vector<Vector>&) = NULL;
	if (Mesh::md2D == pmesh->md_meshDim)
	{
		LeastSquareGradient = &LeastSquareGradient2D;
	}
	else if (Mesh::md3D == pmesh->md_meshDim)
	{
		LeastSquareGradient = &LeastSquareGradient3D;
	}

	for (int i = 0; i < (int)this->p_faceZone->v_faceID.size(); i++)
	{
		int faceID = (int)this->p_faceZone->v_faceID[i];
		int elemID = (this->p_blockMesh->SearchOwnerNeighbor(faceID)).first;
		Vector phiC = EF.GetValue(elemID);
        std::vector<int> nbIDList = pmesh->SearchNearbyElements(elemID, 1);
        std::vector<Vector> deltaPhi;
        std::vector<Vector> deltaR;
		deltaPhi.resize(nbIDList.size());
		deltaR.resize(nbIDList.size());

		for (int k = 0; k < (int)nbIDList.size(); k++)
		{
			int nbID = nbIDList[k];
			deltaPhi[k] = EF.GetValue(nbID) - phiC;
			deltaR[k] = pmesh->v_elem[nbID].center - pmesh->v_elem[elemID].center;
		}
		Tensor grad = LeastSquareGradient(deltaR, deltaPhi);
		Vector CToF = pmesh->v_face[faceID].center - pmesh->v_elem[elemID].center;
		Vector correction = CToF * grad;
		FF.SetValue(faceID, phiC + correction);
		this->v_NumerBC[i].con1 = 1.0;
		this->v_NumerBC[i].con2 = correction;
	}
	return;
}

//Calclate field values at boundary in case of general-type boundary conditon
template<>
void BoundaryFaceZoneField<Vector>::CalculateByABC(ElementField<Vector>& EF, FaceField<Vector>& FF)
{
	if (fbtGeneral != this->fbt_BCType)
	{
		FatalError("Cannot calculate boundary values by CalculateByABC(), boundary type is not given as a general ABC form.");
	}
	for (int i = 0; i < (int)this->p_faceZone->v_faceID.size(); i++)
	{
		int faceID = this->p_faceZone->v_faceID[i];
		Face& patchFace = this->p_blockMesh->v_face[faceID];
		int ownerID = patchFace.n_owner;
		Element& thisElement = this->p_blockMesh->v_elem[ownerID];
        std::vector<Vector> v_FaceValue;
        std::vector<Node> v_FaceArea;
		//collecting informations for calculation
		for (int j = 0; j < (int)thisElement.v_faceID.size(); j++)
		{
			int thisFaceID = thisElement.v_faceID[j];
			Face& thisFace = this->p_blockMesh->v_face[thisFaceID];
			if (((thisFace.center - thisElement.center) & thisFace.area) > 0)
			{
				v_FaceArea.push_back(thisFace.area);
			}
			else
			{
				v_FaceArea.push_back(-(thisFace.area));
			}
			v_FaceValue.push_back(FF.GetValue(thisFaceID));
		}
		Vector OToF = patchFace.center - thisElement.center;
		Vector efUnit = OToF.GetNormal();
		Vector nf;
		if ((efUnit&patchFace.area) > 0)
		{
			nf = patchFace.area.GetNormal();
		}
		else
		{
			nf = -patchFace.area.GetNormal();
		}
		Vector ef = efUnit / (efUnit & nf);
		Vector tf = nf - ef;
		Vector temp = (v_FaceArea[0] & tf)*v_FaceValue[0];
		for (int j = 1; j < (int)v_FaceArea.size(); j++)
		{
			temp = temp + (v_FaceArea[j] & tf)*v_FaceValue[j];
		}
		Scalar tempScalar = v_NumerBC[i].b*ef.Mag() / OToF.Mag();
		Scalar c1 = tempScalar / (v_NumerBC[i].a + tempScalar);
		Vector c2 = (v_NumerBC[i].c - v_NumerBC[i].b*temp / thisElement.volume) / (v_NumerBC[i].a + tempScalar);
		v_NumerBC[i].con1 = c1;
		v_NumerBC[i].con2 = c2;
		Vector updatedValue = v_NumerBC[i].con1*EF.GetValue(ownerID) + v_NumerBC[i].con2;
		FF.SetValue(faceID, updatedValue);
	}
}

//calculating face values on the patch accroding to boundary conditions 
template<>
void BoundaryFaceZoneField<Vector>::CalculateValue(ElementField<Vector>& EF, FaceField<Vector>& FF, const std::string& option)
{
	//Decide which calculation method to use from the preset boundary condtion and the given option
	FieldBoundaryType presetType = this->fbt_BCType;
	if ("boundaryLeft" == option)
	{
        presetType = fbtCaclulated;
	}
	else if ("copy" == option)
	{
        presetType = fbtCopy;
	}
	else if ("extrapolate" == option)
	{
        presetType = fbtExtrapolate;
	}
    else if ("boundaryCondition" == option && presetType != fbtGeneral)
	{
        presetType = fbtCopy;
	}
	else if ("default" != option)
	{
		FatalError("unsupported option " + option + " in BoundaryFaceZoneField::CalculateValue()");
	}
	//Select the corresponding calculation method according to the decision.
    if (fbtCaclulated == presetType)
	{
		this->KeepCurrentValue(FF);
	}
    else if (fbtCopy == presetType || fbtNotSet == presetType)
	{
		this->CopyFromElement(EF, FF);
	}
    else if (fbtExtrapolate == presetType)
	{
		this->ExtrapolateFromElement(EF, FF);
	}
    else if (fbtGeneral == presetType)
	{
		this->CalculateByABC(EF, FF);
	}
	return;
}