#include "../MHT_field/ElementField.h"
#include "../MHT_common/InterpolationTools.h"
#include "../MHT_common/SystemControl.h"
#include "../MHT_mesh/Mesh.h"

// =======================================Scalar Type ElementField=============================================
//constructor
template<>
ElementField<Scalar>::ElementField(Mesh* UGB)
	:BaseField<Scalar>(UGB)
{}


//Initialize with a fixed value
template<>
void ElementField<Scalar>::Initialize(Scalar initalValue)
{
    if (false == Existence())
    {
        this->v_value.resize(this->p_blockMesh->v_elem.size());
        this->fs_status = BaseField<Scalar>::fsCreated;
    }
    for (int i = 0; i < (int)this->v_value.size(); i++)
    {
        SetValue(i, initalValue);
    }
    this->fs_status = fsAssigned;
}
//constructor
template<>
ElementField<Scalar>::ElementField(Mesh* UGB, Scalar value)
	: BaseField<Scalar>(UGB)
{
	this->Initialize(value);
}




//Create and Initialize with a zero value
template<>
void ElementField<Scalar>::Initialize()
{
    this->v_value.resize(this->p_blockMesh->v_elem.size());
    this->fs_status = BaseField<Scalar>::fsCreated;
    for (int i = 0; i < (int)this->v_value.size(); i++)
    {
        SetValue(i, 0.0);
    }
    this->fs_status = fsAssigned;
}


//Initialize with a list of values
template<>
void ElementField<Scalar>::Initialize(std::vector<Scalar> valueList)
{
	if (false == Existence())
	{
		this->v_value.resize(this->p_blockMesh->v_elem.size());
		this->fs_status = BaseField<Scalar>::fsCreated;
	}
	if (v_value.size() != valueList.size())
	{
		FatalError("Cannot initialize field, Numbers of field points and given valueList are not consistent");
	}
	for (int i = 0; i < (int)this->v_value.size(); i++)
	{
		SetValue(i, valueList[i]);
	}
	this->fs_status = fsAssigned;
}

//Initialize with a given function (for 3D case)
template<>
void ElementField<Scalar>::Initialize(Scalar(*udf)(Scalar, Scalar, Scalar))
{
	if (false == Existence())
	{
		this->v_value.resize(this->p_blockMesh->v_elem.size());
		this->fs_status = BaseField<Scalar>::fsCreated;
	}
	for (int i = 0; i < (int)this->v_value.size(); i++)
	{
		Node& xyz = this->p_blockMesh->v_elem[i].center;
		SetValue(i, udf(xyz.x_, xyz.y_, xyz.z_));
	}
	this->fs_status = fsAssigned;
}


//interpolate the value on an lagrange point from this element field
template<>
Scalar ElementField<Scalar>::LagrangeInterpolation
(
Vector lagPoint,	//lagrange point
int elemID			//id of the nearest element
)
{
	if (false == Assignment())
	{
		FatalError("Can not interpolate value to a point. element values are not assingned");
	}
    std::vector<int> elemList = this->p_blockMesh->SearchElementNeighbor(elemID);
	elemList.push_back(elemID);
    std::vector<Vector> rePosition;
    std::vector<Scalar> valueList;
	for (int i = 0; i < (int)elemList.size(); i++)
	{
		rePosition.push_back(this->p_blockMesh->v_elem[elemList[i]].center - lagPoint);
		valueList.push_back(this->v_value[elemList[i]]);
	}
	return Interpolation(rePosition, valueList);
}

template<>
Scalar ElementField<Scalar>::Norm()
{
	Scalar temp = 0.0;

	for (int i = 0; i < (int)this->v_value.size(); i++)
	{
		temp += fabs(this->v_value[i])*this->p_blockMesh->v_elem[i].volume;
		
	}
	return temp;
}

template<>
Scalar ElementField<Scalar>::AverageValue()
{
	Scalar numerator = this->v_value[0] * this->p_blockMesh->v_elem[0].volume;
	Scalar totalVolume = this->p_blockMesh->v_elem[0].volume;
	for (int i = 1; i < (int)this->v_value.size(); i++)
	{
		numerator = numerator + this->v_value[i] * this->p_blockMesh->v_elem[i].volume;
		totalVolume = totalVolume + this->p_blockMesh->v_elem[i].volume;
	}
	return numerator / totalVolume;
}

// "="overload
template<>
ElementField<Scalar>& ElementField<Scalar>::operator = (const ElementField<Scalar>& rhs)
{
    if (false == rhs.Assignment())
	{
		FatalError("Can not assign an empty field");
	}
	this->p_blockMesh = rhs.p_blockMesh;
	this->Initialize(rhs.v_value);
	this->fs_status = BaseField<Scalar>::fsAssigned;
    return *this;
}

// =======================================Vector Type ElementField=============================================
//constructor
template<>
ElementField<Vector>::ElementField(Mesh* UGB)
	:BaseField<Vector>(UGB)
{}

template<>
void ElementField<Vector>::Initialize(Vector initalValue)
{
    if (false == Existence())
    {
        this->v_value.resize(this->p_blockMesh->v_elem.size());
        this->fs_status = BaseField<Vector>::fsCreated;
    }
    for (int i = 0; i < (int)this->v_value.size(); i++)
    {
        SetValue(i, initalValue);
    }
    this->fs_status = fsAssigned;
}
//constructor
template<>
ElementField<Vector>::ElementField(Mesh* UGB, Vector value)
	: BaseField<Vector>(UGB)
{
	this->Initialize(value);
}

//Create and Initialize with a zero value
template<>
void ElementField<Vector>::Initialize()
{
	this->v_value.resize(this->p_blockMesh->v_elem.size());
	this->fs_status = BaseField<Vector>::fsCreated;
	for (int i = 0; i < (int)this->v_value.size(); i++)
	{
		SetValue(i, Vector(0.0, 0.0, 0.0));
	}
	this->fs_status = fsAssigned;
}



//Initialize with a list of values
template<>
void ElementField<Vector>::Initialize(std::vector<Vector> valueList)
{
	if (false == Existence())
	{
		this->v_value.resize(this->p_blockMesh->v_elem.size());
		this->fs_status = BaseField<Vector>::fsCreated;
	}
	if (v_value.size() != valueList.size())
	{
		FatalError("Cannot initialize field, Numbers of field points and given valueList are not consistent");
	}
	for (int i = 0; i < (int)this->v_value.size(); i++)
	{
		SetValue(i, valueList[i]);
	}
	this->fs_status = fsAssigned;
}

//Initialize with a given function (for 3D case)
template<>
void ElementField<Vector>::Initialize(Vector(*udf)(Scalar, Scalar, Scalar))
{
	if (false == Existence())
	{
		this->v_value.resize(this->p_blockMesh->v_elem.size());
		this->fs_status = BaseField<Vector>::fsCreated;
	}
	for (int i = 0; i < (int)this->v_value.size(); i++)
	{
		Node& xyz = this->p_blockMesh->v_elem[i].center;
		SetValue(i, udf(xyz.x_, xyz.y_, xyz.z_));
	}
	this->fs_status = fsAssigned;
}

//interpolate the value on an lagrange point from this element field
template<>
Vector ElementField<Vector>::LagrangeInterpolation
(
Vector lagPoint,	//lagrange point
int elemID			//id of the nearest element
)
{
	if (false == Assignment())
	{
		FatalError("Can not interpolate value to a point. element values are not assingned");
	}
    std::vector<int> elemList = this->p_blockMesh->SearchElementNeighbor(elemID);
	elemList.push_back(elemID);
    std::vector<Vector> rePosition;
    std::vector<Vector> valueList;
	for (int i = 0; i < (int)elemList.size(); i++)
	{
		rePosition.push_back(this->p_blockMesh->v_elem[elemList[i]].center - lagPoint);
		valueList.push_back(this->v_value[elemList[i]]);
	}
	return Interpolation(rePosition, valueList);
}

template<>
Vector ElementField<Vector>::Norm()
{
	Vector temp(0.0, 0.0, 0.0);
	for (int i = 0; i < (int)this->v_value.size(); i++)
	{
		temp.x_ += this->v_value[i].x_*this->p_blockMesh->v_elem[i].volume;
		temp.y_ += this->v_value[i].y_*this->p_blockMesh->v_elem[i].volume;
		temp.z_ += this->v_value[i].z_*this->p_blockMesh->v_elem[i].volume;
	}
	return temp;
}

template<>
Vector ElementField<Vector>::AverageValue()
{
	Vector numerator = this->v_value[0] * this->p_blockMesh->v_elem[0].volume;
	Scalar totalVolume = this->p_blockMesh->v_elem[0].volume;
	for (int i = 1; i < (int)this->v_value.size(); i++)
	{
		numerator = numerator + this->v_value[i] * this->p_blockMesh->v_elem[i].volume;
		totalVolume = totalVolume + this->p_blockMesh->v_elem[i].volume;
	}
	return numerator / totalVolume;
}

// "="overload
template<>
ElementField<Vector>& ElementField<Vector>::operator = (const ElementField<Vector>& rhs)
{
    if (false == rhs.Assignment())
	{
		FatalError("Can not assign an empty field");
	}
	this->p_blockMesh = rhs.p_blockMesh;
	this->Initialize(rhs.v_value);
	this->fs_status = BaseField<Vector>::fsAssigned;
    return *this;
}