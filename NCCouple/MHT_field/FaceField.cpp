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
	
Revised:
	Description:
	1. Convert "non const" input parameters to const;
	2. Convert "non const" function to const function;
	
	Revisor:		Shuai Zhang
	Modified Date:	2018-11-28
\*---------------------------------------------------------------------------*/

#include "../MHT_field/FaceField.h"
#include "../MHT_common/SystemControl.h"
#include "../MHT_mesh/Mesh.h"

// =======================================Scalar Type FaceField=============================================

// Constructor
template<>
FaceField<Scalar>::FaceField(Mesh* UGB)
	:BaseField<Scalar>(UGB)
{}
//Initialize with a fixed value
template<>
void FaceField<Scalar>::Initialize(Scalar initalValue)
{
    if (false == Existence())
    {
        this->v_value.resize(this->p_blockMesh->v_face.size());
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
FaceField<Scalar>::FaceField(Mesh* UGB, Scalar value)
	: BaseField<Scalar>(UGB)
{
	this->Initialize(value);
}

//Create and Initialize with a zero value
template<>
void FaceField<Scalar>::Initialize()
{
	this->v_value.resize(this->p_blockMesh->v_face.size());
	this->fs_status = BaseField<Scalar>::fsCreated;
	for (int i = 0; i < (int)this->v_value.size(); i++)
	{
		SetValue(i, 0.0);
	}
	this->fs_status = fsAssigned;
}


//Initialize with a list of values
template<>
void FaceField<Scalar>::Initialize(std::vector<Scalar> valueList)
{
	if (false == Existence())
	{
		this->v_value.resize(this->p_blockMesh->v_face.size());
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
void FaceField<Scalar>::Initialize(Scalar(*udf)(Scalar, Scalar, Scalar))
{
	if (false == Existence())
	{
		this->v_value.resize(this->p_blockMesh->v_face.size());
		this->fs_status = BaseField<Scalar>::fsCreated;
	}
	for (int i = 0; i < (int)this->v_value.size(); i++)
	{
		Node& xyz = this->p_blockMesh->v_face[i].center;
		SetValue(i, udf(xyz.x_, xyz.y_, xyz.z_));
	}
	this->fs_status = fsAssigned;
}

// "="overload
template<>
FaceField<Scalar>& FaceField<Scalar>::operator = (const FaceField<Scalar>& rhs)
{
    if (false == static_cast<FaceField<Scalar> > (rhs).Assignment())
	{
		FatalError("Can not assign an empty field");
	}
	this->p_blockMesh = rhs.p_blockMesh;
	this->Initialize(rhs.v_value);
	this->fs_status = BaseField<Scalar>::fsAssigned;
    return *this;
}

template<>
std::vector<Scalar> FaceField<Scalar>::ExtractBoundary(const std::string& boundaryName)
{
	Mesh* pmesh = this->p_blockMesh;
	if (false == this->Assignment())
	{
		FatalError("Cannot ExtractBoundary for not passing Assignment check");
	}
	int boundaryID = pmesh->GetBoundaryFaceZoneID(boundaryName);
	if (-1 == boundaryID)
	{
		FatalError("Cannot ExtractBoundary for not finding the boundary you are extract from");
	}
	int faceNum = (int)pmesh->v_boundaryFaceZone[boundaryID].v_faceID.size();
	std::vector<Scalar> valueList;
	valueList.resize(faceNum);

	for (int i = 0; i < faceNum; i++)
	{
		int faceID = pmesh->v_boundaryFaceZone[boundaryID].v_faceID[i];
		valueList[i] = this->GetValue(faceID);
	}
	return valueList;
}

//load values on a specific boundary with a ordered value list
template<>
void FaceField<Scalar>::LoadBoundary(const std::string& boundaryName, const std::vector<Scalar>& valueList)
{
	Mesh* pmesh = this->p_blockMesh;
	if (false == this->Assignment())
	{
		FatalError("Cannot LoadBoundary for not passing Assignment check");
	}
	int boundaryID = pmesh->GetBoundaryFaceZoneID(boundaryName);
	if (-1 == boundaryID)
	{
		FatalError("Cannot LoadBoundary for not finding the boundary you are intending to integrate on");
	}
	int faceNum = (int)pmesh->v_boundaryFaceZone[boundaryID].v_faceID.size();
	if (faceNum != (int)valueList.size())
	{
		FatalError("Cannot LoadBoundary. The number of values is not consistent with that of faces on the boundary");
	}
	for (int i = 0; i < faceNum; i++)
	{
		int faceID = pmesh->v_boundaryFaceZone[boundaryID].v_faceID[i];
		this->SetValue(faceID, valueList[i]);
	}
	return;
}

//load values on a specific boundary with a ordered value list
template<>
void FaceField<Scalar>::AddToBoundary(const std::string& boundaryName, const std::vector<Scalar>& valueList)
{
	Mesh* pmesh = this->p_blockMesh;
	if (false == this->Assignment())
	{
		FatalError("Cannot AddToBoundary for not passing Assignment check");
	}
	int boundaryID = pmesh->GetBoundaryFaceZoneID(boundaryName);
	if (-1 == boundaryID)
	{
		FatalError("Cannot AddToBoundary for not finding the boundary you are intending to integrate on");
	}
	int faceNum = (int)pmesh->v_boundaryFaceZone[boundaryID].v_faceID.size();
	if (faceNum != (int)valueList.size())
	{
		FatalError("Cannot AddToBoundary. The number of values is not consistent with that of faces on the boundary");
	}
	for (int i = 0; i < faceNum; i++)
	{
		int faceID = pmesh->v_boundaryFaceZone[boundaryID].v_faceID[i];
		this->SetValue(faceID, this->GetValue(faceID) +valueList[i]);
	}
	return;
}

template<>
void FaceField<Scalar>::CopyBoundaryFieldsFrom(const FaceField<Scalar>& phif)
{
	if (phif.p_blockMesh != this->p_blockMesh)
	{
		FatalError("Cannot CopyBoundaryFieldsFrom, grids are not consistent!");
	}
	if (false == phif.Assignment())
	{
		FatalError("Cannot CopyBoundaryFieldsFrom for not passing Assignment check");
	}
	if (false == this->Assignment())
	{
		this->Initialize();
	}
	for (int i = 0; i < (int)p_blockMesh->v_boundaryFaceZone.size(); i++)
	{
		for (int j = 0; j < (int)p_blockMesh->v_boundaryFaceZone[i].v_faceID.size(); j++)
		{
			int faceID = p_blockMesh->v_boundaryFaceZone[i].v_faceID[j];
			this->SetValue(faceID, phif.GetValue(faceID));
		}
	}
	return;
}

// =======================================Vector Type FaceField=============================================

// Constructor
template<>
FaceField<Vector>::FaceField(Mesh* UGB)
	:BaseField<Vector>(UGB)
{}

template<>
void FaceField<Vector>::Initialize(Vector initalValue)
{
    if (false == Existence())
    {
        this->v_value.resize(this->p_blockMesh->v_face.size());
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
FaceField<Vector>::FaceField(Mesh* UGB, Vector value)
	: BaseField<Vector>(UGB)
{
	this->Initialize(value);
}

//Create and Initialize with a zero value
template<>
void FaceField<Vector>::Initialize()
{
	this->v_value.resize(this->p_blockMesh->v_face.size());
	this->fs_status = BaseField<Vector>::fsCreated;
	for (int i = 0; i < (int)this->v_value.size(); i++)
	{
		SetValue(i, Vector(0.0,0.0,0.0));
	}
	this->fs_status = fsAssigned;
}

//Initialize with a list of values
template<>
void FaceField<Vector>::Initialize(std::vector<Vector> valueList)
{
	if (false == Existence())
	{
		this->v_value.resize(this->p_blockMesh->v_face.size());
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
void FaceField<Vector>::Initialize(Vector(*udf)(Scalar, Scalar, Scalar))
{
	if (false == Existence())
	{
		this->v_value.resize(this->p_blockMesh->v_face.size());
		this->fs_status = BaseField<Vector>::fsCreated;
	}
	for (int i = 0; i < (int)this->v_value.size(); i++)
	{
		Node& xyz = this->p_blockMesh->v_face[i].center;
		SetValue(i, udf(xyz.x_, xyz.y_, xyz.z_));
	}
	this->fs_status = fsAssigned;
}

// "="overload
template<>
FaceField<Vector>& FaceField<Vector>::operator = (const FaceField<Vector>& rhs)
{
    if (false == static_cast<FaceField<Vector> > (rhs).Assignment())
	{
		FatalError("Can not assign an empty field");
	}
	this->p_blockMesh = rhs.p_blockMesh;
	this->Initialize(rhs.v_value);
	this->fs_status = BaseField<Vector>::fsAssigned;
    return *this;
}

template<>
std::vector<Vector> FaceField<Vector>::ExtractBoundary(const std::string& boundaryName)
{
	Mesh* pmesh = this->p_blockMesh;
	if (false == this->Assignment())
	{
		FatalError("Cannot ExtractBoundary for not passing Assignment check");
	}
	int boundaryID = pmesh->GetBoundaryFaceZoneID(boundaryName);
	if (-1 == boundaryID)
	{
		FatalError("Cannot ExtractBoundary for not finding the boundary you are extract from");
	}
	int faceNum = (int)pmesh->v_boundaryFaceZone[boundaryID].v_faceID.size();
	std::vector<Vector> valueList;
	valueList.resize(faceNum);

	for (int i = 0; i < faceNum; i++)
	{
		int faceID = pmesh->v_boundaryFaceZone[boundaryID].v_faceID[i];
		valueList[i] = this->GetValue(faceID);
	}
	return valueList;
}

//load values on a specific boundary with a ordered value list
template<>
void FaceField<Vector>::LoadBoundary(const std::string& boundaryName, const std::vector<Vector>& valueList)
{
	Mesh* pmesh = this->p_blockMesh;
	if (false == this->Assignment())
	{
		FatalError("Cannot LoadBoundary for not passing Assignment check");
	}
	int boundaryID = pmesh->GetBoundaryFaceZoneID(boundaryName);
	if (-1 == boundaryID)
	{
		FatalError("Cannot LoadBoundary for not finding the boundary you are intending to integrate on");
	}
	int faceNum = (int)pmesh->v_boundaryFaceZone[boundaryID].v_faceID.size();
	if (faceNum != (int)valueList.size())
	{
		FatalError("Cannot LoadBoundary. The number of values is not consistent with that of faces on the boundary");
	}
	for (int i = 0; i < faceNum; i++)
	{
		int faceID = pmesh->v_boundaryFaceZone[boundaryID].v_faceID[i];
		this->SetValue(faceID, valueList[i]);
	}
	return;
}

//load values on a specific boundary with a ordered value list
template<>
void FaceField<Vector>::AddToBoundary(const std::string& boundaryName, const std::vector<Vector>& valueList)
{
	Mesh* pmesh = this->p_blockMesh;
	if (false == this->Assignment())
	{
		FatalError("Cannot AddToBoundary for not passing Assignment check");
	}
	int boundaryID = pmesh->GetBoundaryFaceZoneID(boundaryName);
	if (-1 == boundaryID)
	{
		FatalError("Cannot AddToBoundary for not finding the boundary you are intending to integrate on");
	}
	int faceNum = (int)pmesh->v_boundaryFaceZone[boundaryID].v_faceID.size();
	if (faceNum != (int)valueList.size())
	{
		FatalError("Cannot AddToBoundary. The number of values is not consistent with that of faces on the boundary");
	}
	for (int i = 0; i < faceNum; i++)
	{
		int faceID = pmesh->v_boundaryFaceZone[boundaryID].v_faceID[i];
		this->SetValue(faceID, this->GetValue(faceID) + valueList[i]);
	}
	return;
}

template<>
void FaceField<Vector>::CopyBoundaryFieldsFrom(const FaceField<Vector>& phif)
{
	if (phif.p_blockMesh != this->p_blockMesh)
	{
		FatalError("Cannot CopyBoundaryFieldsFrom, grids are not consistent!");
	}
	if (false == phif.Assignment())
	{
		FatalError("Cannot CopyBoundaryFieldsFrom for not passing Assignment check");
	}
	if (false == this->Assignment())
	{
		this->Initialize();
	}
	for (int i = 0; i < (int)p_blockMesh->v_boundaryFaceZone.size(); i++)
	{
		for (int j = 0; j < (int)p_blockMesh->v_boundaryFaceZone[i].v_faceID.size(); j++)
		{
			int faceID = p_blockMesh->v_boundaryFaceZone[i].v_faceID[j];
			this->SetValue(faceID, phif.GetValue(faceID));
		}
	}
	return;
}

template<>
std::pair<int, Scalar> Compare(const FaceField<Scalar>& phi1, const FaceField<Scalar>& phi2)
{
	std::pair<int, Scalar> result;
	int faceNum = phi1.p_blockMesh->n_faceNum;
	int imax = 0;
	Scalar diffMax = 0.0;
	for (int i = 0; i < faceNum; i++)
	{
		Scalar thisDiff = fabs(phi1.GetValue(i) - phi2.GetValue(i));
		if (thisDiff > diffMax)
		{
			diffMax = thisDiff;
			imax = i;
		}
	}
	result.first = imax;
	result.second = diffMax;
	return result;
}

template<>
std::pair<int, Scalar> Compare(const FaceField<Vector>& phi1, const FaceField<Vector>& phi2)
{
	std::pair<int, Scalar> result;
	int faceNum = phi1.p_blockMesh->n_faceNum;
	int imax = 0;
	Scalar diffMax = 0.0;
	for (int i = 0; i < faceNum; i++)
	{
		Scalar thisDiff = (phi1.GetValue(i) - phi2.GetValue(i)).Mag();
		if (thisDiff > diffMax)
		{
			diffMax = thisDiff;
			imax = i;
		}
	}
	result.first = imax;
	result.second = diffMax;
	return result;
}
