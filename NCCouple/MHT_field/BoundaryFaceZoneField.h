/*---------------------------------------------------------------------------*\
File Name: BoundaryFaceZoneField.h
Author: Kong Ling
Description:
a utility of boundary conditions for FVM 
\*---------------------------------------------------------------------------*/

#ifndef _BoundaryFaceZoneField_
#define _BoundaryFaceZoneField_

#include <vector>

#include "../MHT_common/Configuration.h"
#include "../MHT_common/Vector.h"
#include "../MHT_field/NumerBC.h"

class Mesh;

class FaceZone;

template<class Type>
class ElementField;

template<class Type>
class FaceField;

//a class handling external face zones, playing as a part of faceField 
template<class Type>
class BoundaryFaceZoneField
{
public:
	// Field status enument 
	enum FieldBoundaryType
	{
		fbtGeneral,
		fbtCopy,
		fbtExtrapolate,
		fbtCaclulated,
		fbtNotSet
	};
	//pointer to dependent mesh 
	Mesh* p_blockMesh;
	//pointer to corresponding face zone
	FaceZone* p_faceZone;
	//Boundary type for a specific field
	FieldBoundaryType fbt_BCType;
	//General boundary conditions on all composing faces, empty for non-general field boundary-condition type
    std::vector<NumerBC<Type> > v_NumerBC;

	//constructor
	BoundaryFaceZoneField(Mesh* UGB, FaceZone* FZ);

	//return true if boundary condition is imposed otherwise return false
    bool BCExist() const;

	//set an uniform numerical boundary condition
	void SetNumerBC(Scalar, Scalar, Type);

	//set a non-uniform numerical boundary condition
	void SetNumerBC(Scalar, Scalar, Type(*udf)(Scalar,Scalar,Scalar));

	//set a non-uniform numerical boundary condition
    void SetNumerBC(std::vector<NumerBC<Type> >&);

	//set special-type boundary condition like copy and extrapolate
	void SetSpecialBC(FieldBoundaryType);

	//Copy boundary condition to a target boundary for a specific field
	void CopyBoundaryConditionTo(BoundaryFaceZoneField& targetBoundary);

	//keep current value and write corresponding coefficients
	void KeepCurrentValue(FaceField<Type>&);

	//copy values from owner elements
	void CopyFromElement(ElementField<Type>& EF, FaceField<Type>& FF);

	//extrapolate values linearly from elements
	void ExtrapolateFromElement(ElementField<Type>& EF, FaceField<Type>& FF);

	//Calclate field values at boundary in case of general-type boundary conditon
	void CalculateByABC(ElementField<Type>& EF, FaceField<Type>& FF);

	//calculating face values on the patch accroding to boundary conditions 
    void CalculateValue(ElementField<Type>& EF, FaceField<Type>& FF, const std::string& option="default");
	
	//End of member functions
};

#endif
