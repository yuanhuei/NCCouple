
#ifndef _NodeField_
#define _NodeField_

#include "../MHT_common/Configuration.h"
#include "../MHT_common/Vector.h"
#include "../MHT_common/Tensor.h"

#include "../MHT_field/BaseField.h"

//element field
template<class Type>
class NodeField : public BaseField<Type>
{
public:
	// Constructor
	NodeField(Mesh* UGB);

	// Constructor
	NodeField(Mesh* UGB, Type);

	//Create and initialize with a zero value
	void Initialize();

	//Initialize with a fixed value
	void Initialize(Type initalValue);

	//Initialize with a list of values
	void Initialize(std::vector<Type> valueList);

	//Initialize with a given function
	void Initialize(Type(udf)(Scalar, Scalar, Scalar));

	// "="overload
    NodeField<Type>& operator = (const NodeField<Type>& rhs);
};

//====================External functions=======================
template<class Type>
NodeField<Type> operator + (const NodeField<Type>&, const NodeField<Type>&);

template<class Type>
NodeField<Type> operator - (const NodeField<Type>&);

template<class Type>
NodeField<Type> operator - (const NodeField<Type>&, const NodeField<Type>&);

template<class Type>
NodeField<Type> operator * (const Scalar, const NodeField<Type>&);

template<class Type>
NodeField<Type> operator * (const NodeField<Scalar>&, const NodeField<Type>&);

template<class Type>
NodeField<Type> operator / (const NodeField<Type>&, const NodeField<Scalar>&);

NodeField<Scalar> operator / (Scalar, const NodeField<Scalar>&);

NodeField<Scalar> operator & (const NodeField<Vector>&, const NodeField<Vector>&);

#endif
