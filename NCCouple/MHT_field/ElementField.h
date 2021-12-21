/*---------------------------------------------------------------------------*\
File Name:
ElementField.h

Description:
defining classes playing as composing parts in Field
\*---------------------------------------------------------------------------*/

#ifndef _ElementField_
#define _ElementField_

#include "../MHT_common/Configuration.h"
#include "../MHT_common/Vector.h"
#include "../MHT_common/Tensor.h"
#include "../MHT_field/BaseField.h"

//element field
template<class Type>
class ElementField : public BaseField<Type>
{
public:
	//constructor
	ElementField(Mesh* UGB);

	// Constructor
	ElementField(Mesh* UGB, Type);

	//Create and Initialize with a zero value
	void Initialize();

	//Initialize with a fixed value (create first if not created)
	void Initialize(Type initalValue);

	//Initialize with a list of values (create first if not created)
	void Initialize(std::vector<Type> valueList);

	//Initialize with a given function (create first if not created)
	void Initialize(Type(udf)(Scalar, Scalar, Scalar));

	//interpolate the value on an lagrange point from a given element
	Type LagrangeInterpolation(Vector lagPoint,	int elemID);

	Type Norm();

	Type AverageValue();

	// "="overload
    ElementField<Type>& operator = (const ElementField<Type>& rhs);

};

//=======================External functions=========================
template<class Type>
ElementField<Type> operator + (const ElementField<Type>&, const ElementField<Type>&);

template<class Type>
ElementField<Type> operator - (const ElementField<Type>&);

template<class Type>
ElementField<Type> operator - (const ElementField<Type>&, const ElementField<Type>&);

template<class Type>
ElementField<Type> operator * (const Scalar, const ElementField<Type>&);

template<class Type>
ElementField<Type> operator * (const ElementField<Scalar>&, const ElementField<Type>&);

template<class Type>
ElementField<Type> operator / (const ElementField<Type>&, const ElementField<Scalar>&);

ElementField<Scalar> operator / (Scalar, const ElementField<Scalar>&);

ElementField<Vector> operator * (const ElementField<Tensor>&, const ElementField<Vector>&);

ElementField<Scalar> operator & (const ElementField<Vector>&, const ElementField<Vector>&);

template<class Type>
std::pair<int, Scalar> Compare(const ElementField<Type>&, const ElementField<Type>&);

#endif
