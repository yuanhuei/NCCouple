/*---------------------------------------------------------------------------*\
File Name:
	ElementFieldOperator.cpp

Description:
	

	Author:		Kong Ling
	Date: 

	
Revised:
	Description:
	1. Convert "non const" input parameters to const;
	2. Convert "non const" function to const function;
	
	Revisor:		Shuai Zhang
	Modified Date:	2018-11-28
\*---------------------------------------------------------------------------*/

#include "../MHT_field/NodeField.h"
#include "../MHT_common/SystemControl.h"
#include "../MHT_mesh/Mesh.h"

// =======================================Scalar Type NodeField=============================================
// Constructor

template<>
NodeField<Scalar>::NodeField(Mesh* UGB): BaseField<Scalar>(UGB)
{}
//Initialize with a fixed value
template<>
void NodeField<Scalar>::Initialize(Scalar initalValue)
{
    if (false == Existence())
    {
        this->v_value.resize(this->p_blockMesh->v_node.size());
        this->fs_status = BaseField<Scalar>::fsCreated;
    }
    for (int i = 0; i < (int)this->v_value.size(); i++)
    {
        SetValue(i, initalValue);
    }
    this->fs_status = fsAssigned;
}
template<>
NodeField<Scalar>::NodeField(Mesh* UGB, Scalar value)
	: BaseField<Scalar>(UGB)
{
	Initialize(value);
}


//Create and Initialize with a zero value
template<>
void NodeField<Scalar>::Initialize()
{
	if (false == Existence())
	{
		this->v_value.resize(this->p_blockMesh->v_node.size());
		this->fs_status = BaseField<Scalar>::fsCreated;
	}
	for (int i = 0; i < (int)this->v_value.size(); i++)
	{
		SetValue(i, 0.0);
	}
	this->fs_status = fsAssigned;
}



//Initialize with a list of values
template<>
void NodeField<Scalar>::Initialize(std::vector<Scalar> valueList)
{
	if (false == Existence())
	{
		this->v_value.resize(this->p_blockMesh->v_node.size());
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
void NodeField<Scalar>::Initialize(Scalar(*udf)(Scalar, Scalar, Scalar))
{
	if (false == Existence())
	{
		this->v_value.resize(this->p_blockMesh->v_node.size());
		this->fs_status = BaseField<Scalar>::fsCreated;
	}
	for (int i = 0; i < (int)this->v_value.size(); i++)
	{
		Node& xyz = this->p_blockMesh->v_node[i];
		SetValue(i, udf(xyz.x_, xyz.y_, xyz.z_));
	}
	this->fs_status = fsAssigned;
}

// "="overload
template<>
NodeField<Scalar>& NodeField<Scalar>::operator = (const NodeField<Scalar>& rhs)
{
    if (false == static_cast<NodeField<Scalar> > (rhs).Assignment())
	{
		FatalError("Can not assign an empty field");
	}
	this->p_blockMesh = rhs.p_blockMesh;
	this->Initialize(rhs.v_value);
	this->fs_status = BaseField<Scalar>::fsAssigned;
    return *this;
}

// =======================================Vector Type NodeField=============================================
template<>
NodeField<Vector>::NodeField(Mesh* UGB): BaseField<Vector>(UGB)
{}
template<>
void NodeField<Vector>::Initialize(Vector initalValue)
{
    if (false == Existence())
    {
        this->v_value.resize(this->p_blockMesh->v_node.size());
        this->fs_status = BaseField<Vector>::fsCreated;
    }
    for (int i = 0; i < (int)this->v_value.size(); i++)
    {
        SetValue(i, initalValue);
    }
    this->fs_status = fsAssigned;
}
template<>
NodeField<Vector>::NodeField(Mesh* UGB, Vector value)
	: BaseField<Vector>(UGB)
{
	Initialize(value);
}

//Create and Initialize with a zero value
template<>
void NodeField<Vector>::Initialize()
{
	if (false == Existence())
	{
		this->v_value.resize(this->p_blockMesh->v_node.size());
		this->fs_status = BaseField<Vector>::fsCreated;
	}
	for (int i = 0; i < (int)this->v_value.size(); i++)
	{
		SetValue(i, Vector(0.0,0.0,0.0));
	}
	this->fs_status = fsAssigned;
}



//Initialize with a list of values
template<>
void NodeField<Vector>::Initialize(std::vector<Vector> valueList)
{
	if (false == Existence())
	{
		this->v_value.resize(this->p_blockMesh->v_node.size());
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
void NodeField<Vector>::Initialize(Vector(*udf)(Scalar, Scalar, Scalar))
{
	if (false == Existence())
	{
		this->v_value.resize(this->p_blockMesh->v_node.size());
		this->fs_status = BaseField<Vector>::fsCreated;
	}
	for (int i = 0; i < (int)this->v_value.size(); i++)
	{
		Node& xyz = this->p_blockMesh->v_node[i];
		SetValue(i, udf(xyz.x_, xyz.y_, xyz.z_));
	}
	this->fs_status = fsAssigned;
}

// "="overload
template<>
NodeField<Vector>& NodeField<Vector>::operator = (const NodeField<Vector>& rhs)
{
    if (false == static_cast<NodeField<Vector> > (rhs).Assignment())
	{
		FatalError("Can not assign an empty field");
	}
	this->p_blockMesh = rhs.p_blockMesh;
	this->Initialize(rhs.v_value);
	this->fs_status = BaseField<Vector>::fsAssigned;
    return *this;
}