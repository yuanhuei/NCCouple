
#ifndef _BaseField_
#define _BaseField_

#include <vector>
#include <string>
#include "../MHT_common/Configuration.h"
#include "../MHT_common/Vector.h"
class Mesh;

// Base class Field
template<class Type>
class BaseField
{
public:

	// Field status enument 
	enum FieldStatus
	{
		fsNotExist,
		fsCreated,
		fsAssigned
	};

	//pointer to dependent block mesh
	Mesh* p_blockMesh;

	//indicating whether the value list status
	FieldStatus fs_status;

	//values of a corresponding face/node/element field
    std::vector<Type> v_value;

	// Constructor
	BaseField(Mesh* UGB);

	//check wether the field is existing;
    bool Existence() const;

	//check wether the field have values assigned;
    bool Assignment() const;

	//clear the value list
	void Destroy();

	//reading the value on a face/node/element field with a given ID
    Type GetValue(const int ElementID) const;

	//writing value on on a face/node/element field a with given ID;
    void SetValue(const int ElementID, const Type& value);

	//set min and max for a Scalar field
	void Bound(Scalar, Scalar);

	void WriteData(std::string);

	void ReadData(std::string);

};

template<>
Vector BaseField<Vector>::GetValue(const int ElementID) const;
template<>
void BaseField<Vector>::SetValue(const int ElementID, const Vector& value);
#endif
