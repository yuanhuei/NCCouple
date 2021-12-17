/*---------------------------------------------------------------------------*\
File Name:
	FaceField.h

Description:
	defining classes playing as composing parts in Field

	Author:		Kong Ling
	Date: 2016-11-28

Revised:
	Description:	1.Convert FaceField&ElementField(bool b_exist;	bool b_assigned;) to enum FieldStatus. And 
					made corresponding modefication. 
					2.Convert FaceField&ElementField(WriteValue) to SetValue.
					3.Made a base for Field (BaseField)
	Revisor:		ShuaiZhang
	Modified Date:	2017-01-04
\*---------------------------------------------------------------------------*/

#ifndef _FaceField_
#define _FaceField_

#include "../MHT_common/Configuration.h"
#include "../MHT_common/Vector.h"
#include "../MHT_common/Tensor.h"
#include "../MHT_field/BaseField.h"

//face field
template<class Type>
class FaceField : public BaseField<Type>
{
public:
	// Constructor
	FaceField(Mesh* UGB);

	// Constructor
	FaceField(Mesh* UGB, Type);

	//Create and Initialize with a zero value
	void Initialize();

	//Initialize with a fixed value
	void Initialize(Type initalValue);

	//Initialize with a list of values
	void Initialize(std::vector<Type> valueList);

	//Initialize with a given function
	void Initialize(Type(udf)(Scalar, Scalar, Scalar));

	// "="overload
    FaceField<Type>& operator = (const FaceField<Type>& rhs);

	//extract values on a specfic boundary
	std::vector<Type> ExtractBoundary(const std::string&);

	//load values on a specific boundary with a ordered value list
	void LoadBoundary(const std::string&, const std::vector<Type>&);

	//Kong Ling supplemented on 2019/12/23
	void AddToBoundary(const std::string&, const std::vector<Type>&);

	//Kong Ling supplemented on 2020/01/17
	void CopyBoundaryFieldsFrom(const FaceField<Type>&);

};

//====================External functions=======================
template<class Type>
FaceField<Type> operator + (const FaceField<Type>&, const FaceField<Type>&);

template<class Type>
FaceField<Type> operator - (const FaceField<Type>&);

template<class Type>
FaceField<Type> operator - (const FaceField<Type>&, const FaceField<Type>&);

template<class Type>
FaceField<Type> operator * (const Scalar, const FaceField<Type>&);

template<class Type>
FaceField<Type> operator *(const FaceField<Scalar>&, const FaceField<Type>&);

template<class Type>
FaceField<Type> operator / (const FaceField<Type>&, const FaceField<Scalar>&);

FaceField<Scalar> operator / (Scalar, const FaceField<Scalar>&);

FaceField<Scalar> operator & (const FaceField<Vector>&, const FaceField<Vector>&);

template<class Type>
std::pair<int, Scalar> Compare(const FaceField<Type>&, const FaceField<Type>&);

#endif
