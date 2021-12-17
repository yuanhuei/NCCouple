
#ifndef _FieldOperator_
#define _FieldOperator_

#include "../MHT_common/Configuration.h"
#include "../MHT_common/Vector.h"
#include "../MHT_field/Field.h"

//=================basic operators=================
template<class Type>
Field<Type> operator + (const Field<Type>&, const Field<Type>&);

template<class Type>
Field<Type> operator - (const Field<Type>&);

template<class Type>
Field<Type> operator - (const Field<Type>&, const Field<Type>&);

template<class Type>
Field<Type> operator * (const Scalar, const Field<Type>&);

template<class Type>
Field<Type> operator * (const Field<Scalar>&, const Field<Type>&);

template<class Type>
Field<Type> operator / (const Field<Type>&, const Field<Scalar>&);

Field<Scalar> operator / (Scalar, const Field<Scalar>&);

Field<Scalar> operator & (const Field<Vector>&, const Field<Vector>&);

#endif
