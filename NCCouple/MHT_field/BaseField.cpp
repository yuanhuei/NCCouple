/*---------------------------------------------------------------------------*\
File Name:
	BaseField.cpp

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
	
Revised:
	Description:
	1. Convert "non const" input parameters to const;
	2. Convert "non const" function to const function;
	
	Revisor:		Shuai Zhang
	Modified Date:	2018-11-28
\*---------------------------------------------------------------------------*/

#include <algorithm>
#include "../MHT_field/BaseField.h"
#include "../MHT_common/Vector.h"
#include "../MHT_common/SystemControl.h"

// =======================================Scalar Type BaseField=============================================

// Constructor
template<>
BaseField<Scalar>::BaseField(Mesh* UGB)
	:p_blockMesh(UGB),
	fs_status(fsNotExist)
{}

//check wether the field is existing;
template<>
bool BaseField<Scalar>::Existence() const
{
	if (BaseField<Scalar>::fsNotExist == this->fs_status)
	{
		return false;
	}
	else
	{
		return true;
	}
}

//check wether the field have values assigned;
template<>
bool BaseField<Scalar>::Assignment() const
{
	if (BaseField<Scalar>::fsAssigned != this->fs_status)
	{
		return false;
	}
	else
	{
		return true;
	}
}

//clear the value list
template<>
void BaseField<Scalar>::Destroy()
{
	this->fs_status = BaseField<Scalar>::fsNotExist;
	this->v_value.clear();
}

//reading the value on a face/node/element field with a given ID
template<>
Scalar BaseField<Scalar>::GetValue(const int ElementID) const
{
    if (false == Assignment())
    {
        FatalError("Can get value. values are not assingned");
    }
    return this->v_value[ElementID];
}

//writing value on on a face/node/element field a with given ID;
template<>
void BaseField<Scalar>::SetValue(const int ElementID,const Scalar& value)
{
	if (false == Existence())
	{
		FatalError("In BaseField<Scalar>::SetValue, value list is not created");
	}
	this->v_value[ElementID] = value;
}

template<>
void BaseField<Scalar>::Bound(Scalar minValue, Scalar maxValue)
{
	if (minValue > maxValue)
	{
		FatalError("In BaseField<Scalar>::Bound, minimum value is greater than the maximum one");
	}
	if (false == Assignment())
	{
		FatalError("In BaseField<Scalar>::Bound, value list is not given");
	}
	for (int i = 0; i < (int)this->v_value.size(); i++)
	{
		Scalar phiC = this->GetValue(i);
		Scalar correctedValue = Max(minValue, phiC);
		this->SetValue(i, Min(Max(minValue, phiC), maxValue));
	}
	return;
}

template<>
void BaseField<Scalar>::WriteData(std::string filename)
{
	int dataNumber = this->v_value.size();
	std::ofstream outFile(filename);
	outFile << dataNumber << std::endl;
	for (int i = 0; i < dataNumber; i++)
	{
		outFile << this->GetValue(i) << std::endl;
	}
	outFile.close();
	return;
}

template<>
void BaseField<Scalar>::ReadData(std::string filename)
{
	int dataNumber = 0;
	std::ifstream inFile(filename);
	inFile >> dataNumber;
	if (dataNumber != this->v_value.size())
	{
		FatalError("in BaseField::ReadData(), the numbers of data points are not consitent");
	}
	for (int i = 0; i < dataNumber; i++)
	{
		Scalar value = 0.0;
		inFile >> value;
		this->SetValue(i, value);
	}
	inFile.close();
	return;
}

template<>
void BaseField<Vector>::WriteData(std::string filename)
{
	int dataNumber = this->v_value.size();
	std::ofstream outFile(filename);
	outFile << dataNumber << std::endl;
	for (int i = 0; i < dataNumber; i++)
	{
		Vector temp = this->GetValue(i);
		outFile << temp.x_ << "\t";
		outFile << temp.y_ << "\t";
		outFile << temp.z_ << std::endl;
	}
	outFile.close();
	return;
}

template<>
void BaseField<Vector>::ReadData(std::string filename)
{
	int dataNumber = 0;
	std::ifstream inFile(filename);
	inFile >> dataNumber;
	if (dataNumber != this->v_value.size())
	{
		FatalError("in BaseField::ReadData(), the numbers of data points are not consitent");
	}
	for (int i = 0; i < dataNumber; i++)
	{
		Vector value(0.0, 0.0, 0.0);
		inFile >> value.x_ >> value.y_ >> value.z_;
		this->SetValue(i, value);
	}
	inFile.close();
	return;
}

// =======================================Vector Type BaseField=============================================

// Constructor
template<>
BaseField<Vector>::BaseField(Mesh* UGB)
	:p_blockMesh(UGB),
	fs_status(fsNotExist)
{}

template<>
bool BaseField<Vector>::Existence() const
{
	if (BaseField<Vector>::fsNotExist == this->fs_status)
	{
		return false;
	}
	else
	{
		return true;
	}
}

template<>
bool BaseField<Vector>::Assignment() const
{
	if (BaseField<Vector>::fsAssigned != this->fs_status)
	{
		return false;
	}
	else
	{
		return true;
	}
}

//clear the value list
template<>
void BaseField<Vector>::Destroy()
{
	this->fs_status = BaseField<Vector>::fsNotExist;
	this->v_value.clear();
}

//reading the value on a face/node/element field with a given ID
template<>
Vector BaseField<Vector>::GetValue(const int ElementID) const
{
    if (false == Assignment())
    {
        FatalError("Can get value. values are not assingned");
    }
	return this->v_value[ElementID];
}

//writing value on on a face/node/element field a with given ID;
template<>
void BaseField<Vector>::SetValue(const int ElementID, const Vector& value)
{
    if (false == Existence())
    {
        FatalError("Can not set value. value list is not created");
    }
    this->v_value[ElementID] = value;
}