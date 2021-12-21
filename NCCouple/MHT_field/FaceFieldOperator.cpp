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
#include "../MHT_common/Vector.h"
#include "../MHT_common/SystemControl.h"
#include "../MHT_field/FaceField.h"
#include "../MHT_mesh/Mesh.h"

template<>
FaceField<Scalar> operator + (const FaceField<Scalar>& EF1, const FaceField<Scalar>& EF2)
{
	if (EF1.p_blockMesh != EF2.p_blockMesh)
	{
		FatalError("in operation + between two Element Fields, the grids are not consistent");
	}

    if (false == EF1.Assignment() || false == EF2.Assignment())
	{
		FatalError("Can not proceed operation + between two Element Fields, did not pass Assignment Check");
	}
	FaceField<Scalar> result(EF1.p_blockMesh);
    std::vector<Scalar> List;
	for (int i = 0; i < EF1.p_blockMesh->n_faceNum; i++)
	{
        List.push_back(EF1.v_value[i] + EF2.v_value[i]);
	}
	result.Initialize(List);
	return result;
}

template<>
FaceField<Vector> operator + (const FaceField<Vector>& EF1, const FaceField<Vector>& EF2)
{
	if (EF1.p_blockMesh != EF2.p_blockMesh)
	{
		FatalError("in operation + between two Element Fields, the grids are not consistent");
	}

    if (false == EF1.Assignment() || false == EF2.Assignment())
	{
		FatalError("Can not proceed operation + between two Element Fields, did not pass Assignment Check");
	}
    FaceField<Vector> result(EF1.p_blockMesh);
    std::vector<Vector> List;
	for (int i = 0; i < EF1.p_blockMesh->n_faceNum; i++)
    {
        List.push_back(EF1.v_value[i] + EF2.v_value[i]);
	}
	result.Initialize(List);
	return result;
}

template<>
FaceField<Scalar> operator - (const FaceField<Scalar>& FF)
{
	if (false == FF.Assignment())
	{
		FatalError("Can not proceed operation - for a field, did not pass Assignment Check");
	}
	FaceField<Scalar> result(FF.p_blockMesh);
	std::vector<Scalar> List;
	for (int i = 0; i < FF.p_blockMesh->n_faceNum; i++)
	{
		List.push_back(-FF.v_value[i]);
	}
	result.Initialize(List);
	return result;
}

template<>
FaceField<Vector> operator - (const FaceField<Vector>& FF)
{
	if (false == FF.Assignment())
	{
		FatalError("Can not proceed operation - for a Field, did not pass Assignment Check");
	}
	FaceField<Vector> result(FF.p_blockMesh);
	std::vector<Vector> List;
	for (int i = 0; i < FF.p_blockMesh->n_faceNum; i++)
	{
		List.push_back(Vector(0.0, 0.0, 0.0) -FF.v_value[i]);
	}
	result.Initialize(List);
	return result;
}

template<>
FaceField<Scalar> operator - (const FaceField<Scalar>& EF1, const FaceField<Scalar>& EF2)
{
	if (EF1.p_blockMesh != EF2.p_blockMesh)
	{
		FatalError("in operation - between two Element Fields, the grids are not consistent");
	}

    if (false == EF1.Assignment() || false == EF2.Assignment())
	{
		FatalError("Can not proceed operation - between two Element Fields, did not pass Assignment Check");
	}
	FaceField<Scalar> result(EF1.p_blockMesh);
    std::vector<Scalar> List;
	for (int i = 0; i < EF1.p_blockMesh->n_faceNum; i++)
	{
		List.push_back(EF1.v_value[i] - EF2.v_value[i]);
	}
	result.Initialize(List);
	return result;
}

template<>
FaceField<Vector> operator - (const FaceField<Vector>& EF1, const FaceField<Vector>& EF2)
{
	if (EF1.p_blockMesh != EF2.p_blockMesh)
	{
		FatalError("in operation - between two Element Fields, the grids are not consistent");
	}

    if (false == EF1.Assignment() || false == EF2.Assignment())
	{
		FatalError("Can not proceed operation - between two Element Fields, did not pass Assignment Check");
	}
	FaceField<Vector> result(EF1.p_blockMesh);
    std::vector<Vector> List;
	for (int i = 0; i < EF1.p_blockMesh->n_faceNum; i++)
	{
		List.push_back(EF1.v_value[i] - EF2.v_value[i]);
	}
	result.Initialize(List);
	return result;
}

//=======================External Functions===========================
template<>
FaceField<Scalar> operator * (const Scalar coef, const FaceField<Scalar>& phi)
{
    if (false == phi.Assignment())
	{
		FatalError("Can not proceed operation * between a Scalar and a FaceField, did not pass Assignment Check");
    }
	FaceField<Scalar> result(phi);
	for (int i = 0; i < phi.p_blockMesh->n_faceNum; i++)
    {
        result.SetValue(i, coef*phi.GetValue(i));
	}
	return result;
}

template<>
FaceField<Vector> operator * (const Scalar coef, const FaceField<Vector>& phi)
{
    if (false == phi.Assignment())
	{
		FatalError("Can not proceed operation * between a Scalar and a FaceField, did not pass Assignment Check");
	}
	FaceField<Vector> result(phi);
	for (int i = 0; i < phi.p_blockMesh->n_faceNum; i++)
	{
		result.SetValue(i, coef*phi.GetValue(i));
	}
	return result;
}

// "*"overload
template<>
FaceField<Scalar> operator * (const FaceField<Scalar>& FF1, const FaceField<Scalar>& FF2)
{
	if (FF1.p_blockMesh != FF2.p_blockMesh)
	{
		FatalError("in operation * between two Face Fields, the grids are not consistent");
	}

    if (false == FF1.Assignment() || false == FF2.Assignment())
	{
		FatalError("Can not proceed operation * between two Face Fields, did not pass Assignment Check");
	}
	FaceField<Scalar> result(FF1.p_blockMesh);
    std::vector<Scalar> scalarList;
	for (int i = 0; i < FF1.p_blockMesh->n_faceNum; i++)
	{
		scalarList.push_back(FF1.v_value[i] * FF2.v_value[i]);
	}
	result.Initialize(scalarList);
	return result;
}

template<>
FaceField<Vector> operator * (const FaceField<Scalar>&EF1, const FaceField<Vector>&EF2)
{
	if (EF1.p_blockMesh != EF2.p_blockMesh)
	{
		FatalError("in operation * between a Scalar Field and a Vector one, the grids are not consistent");
	}
    if (false == EF1.Assignment() || false == EF2.Assignment())
	{
		FatalError("Can not proceed operation * between a Scalar Field and a Vector one for not pass Assignment Check");
	}
	FaceField<Vector> result(EF1.p_blockMesh);
    std::vector<Vector> vectorList;
	for (int i = 0; i < EF1.p_blockMesh->n_faceNum; i++)
	{
		vectorList.push_back(EF1.v_value[i] * EF2.v_value[i]);
	}
	result.Initialize(vectorList);
	return result;
}

template<>
FaceField<Scalar> operator / (const FaceField<Scalar>& phi1, const FaceField<Scalar>& phi2)
{
	if (phi1.p_blockMesh != phi2.p_blockMesh)
	{
		FatalError("in operation / between two Scalar Fields, the grids are not consistent");
	}
    if (false == phi1.Assignment() || false == phi2.Assignment())
	{
		FatalError("Can not proceed operation / between two Scalar Fields for not pass Assignment Check");
	}
	FaceField<Scalar> result(phi1.p_blockMesh);
	result.Initialize();
	for (int i = 0; i < phi1.p_blockMesh->n_faceNum; i++)
	{
		result.SetValue(i, phi1.v_value[i] / phi2.v_value[i]);
	}
	return result;
}

template<>
FaceField<Vector> operator / (const FaceField<Vector>& phi1, const FaceField<Scalar>& phi2)
{
	if (phi1.p_blockMesh != phi2.p_blockMesh)
	{
		FatalError("in operation / between two Scalar Fields, the grids are not consistent");
    }

    if (false == phi1.Assignment() || false == phi2.Assignment())
	{
		FatalError("Can not proceed operation / between two Scalar Fields for not pass Assignment Check");
	}
	FaceField<Vector> result(phi1.p_blockMesh);
	result.Initialize();
	for (int i = 0; i < phi1.p_blockMesh->n_faceNum; i++)
	{
		result.SetValue(i, phi1.v_value[i] / phi2.v_value[i]);
	}
	return result;
}

FaceField<Scalar> operator / (Scalar alpha, const FaceField<Scalar>& phi)
{
	if (false == phi.Assignment())
	{
		FatalError("Can not proceed operation / between a Scalar and a Scalar FaceField");
	}
	FaceField<Scalar> result(phi.p_blockMesh);
	result.Initialize();
	for (int i = 0; i < phi.p_blockMesh->n_faceNum; i++)
	{
		result.SetValue(i, alpha / phi.GetValue(i));
	}
	return result;
}

FaceField<Scalar> operator & (const FaceField<Vector>& phi1, const FaceField<Vector>& phi2)
{
	if (phi1.p_blockMesh != phi2.p_blockMesh)
	{
		FatalError("in operation & between two Vector fields, the grids are not consistent");
	}
    if (false == phi1.Assignment() || false == phi1.Assignment())
	{
		FatalError("Can not proceed operation & between two Vector fields for not pass Assignment Check");
	}
	FaceField<Scalar> result(phi1.p_blockMesh);
	result.Initialize();
	for (int i = 0; i < phi1.p_blockMesh->n_faceNum; i++)
	{
		result.SetValue(i, phi1.GetValue(i) & phi2.GetValue(i));
	}
	return result;
}
