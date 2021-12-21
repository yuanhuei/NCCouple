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
#include "../MHT_field/ElementField.h"
#include "../MHT_common/SystemControl.h"
#include "../MHT_mesh/Mesh.h"

template<>
ElementField<Scalar> operator + (const ElementField<Scalar>& EF1, const ElementField<Scalar>& EF2)
{
	if (EF1.p_blockMesh != EF2.p_blockMesh)
	{
		FatalError("in operation + between two Element Fields, the grids are not consistent");
	}
	if (false == EF1.Assignment() || false == EF2.Assignment())
	{
		FatalError("Can not proceed operation + between two Element Fields, did not pass Assignment Check");
	}
	ElementField<Scalar> result(EF1.p_blockMesh);
    std::vector<Scalar> List;
    for (int i = 0; i < (int)EF1.p_blockMesh->n_elemNum; i++)
	{
		List.push_back(EF1.v_value[i] + EF2.v_value[i]);
	}
	result.Initialize(List);
	return result;
}

template<>
ElementField<Vector> operator + (const ElementField<Vector>& EF1, const ElementField<Vector>& EF2)
{
	if (EF1.p_blockMesh != EF2.p_blockMesh)
	{
		FatalError("in operation + between two Element Fields, the grids are not consistent");
	}
	if (false == EF1.Assignment() || false == EF2.Assignment())
	{
		FatalError("Can not proceed operation + between two Element Fields, did not pass Assignment Check");
	}
	ElementField<Vector> result(EF1.p_blockMesh);
    std::vector<Vector> List;
    for (int i = 0; i < (int)EF1.p_blockMesh->n_elemNum; i++)
	{
		List.push_back(EF1.v_value[i] + EF2.v_value[i]);
	}
	result.Initialize(List);
	return result;
}

template<>
ElementField<Scalar> operator - (const ElementField<Scalar>& EF)
{
	if (false == EF.Assignment())
	{
		FatalError("Can not proceed operation - between, did not pass Assignment Check");
	}
	ElementField<Scalar> result(EF.p_blockMesh);
	std::vector<Scalar> List;
	for (int i = 0; i < EF.p_blockMesh->n_elemNum; i++)
	{
		List.push_back(-EF.v_value[i]);
	}
	result.Initialize(List);
	return result;
}

template<>
ElementField<Vector> operator - (const ElementField<Vector>& EF)
{
	if (false == EF.Assignment())
	{
		FatalError("Can not proceed operation - between, did not pass Assignment Check");
	}
	ElementField<Vector> result(EF.p_blockMesh);
	std::vector<Vector> List;
	for (int i = 0; i < EF.p_blockMesh->n_elemNum; i++)
	{
		List.push_back(Vector(0.0, 0.0, 0.0) - EF.v_value[i]);
	}
	result.Initialize(List);
	return result;
}

template<>
ElementField<Scalar> operator - (const ElementField<Scalar>& EF1, const ElementField<Scalar>& EF2)
{
	if (EF1.p_blockMesh != EF2.p_blockMesh)
	{
		FatalError("in operation - between two Element Fields, the grids are not consistent");
	}
	if (false == EF1.Assignment() || false == EF2.Assignment())
	{
		FatalError("Can not proceed operation - between two Element Fields, did not pass Assignment Check");
	}
	ElementField<Scalar> result(EF1.p_blockMesh);
    std::vector<Scalar> List;
	for (int i = 0; i < EF1.p_blockMesh->n_elemNum; i++)
	{
		List.push_back(EF1.v_value[i] - EF2.v_value[i]);
	}
	result.Initialize(List);
	return result;
}

template<>
ElementField<Vector> operator - (const ElementField<Vector>& EF1, const ElementField<Vector>& EF2)
{
	if (EF1.p_blockMesh != EF2.p_blockMesh)
	{
		FatalError("in operation - between two Element Fields, the grids are not consistent");
	}
	if (false == EF1.Assignment() || false == EF2.Assignment())
	{
		FatalError("Can not proceed operation - between two Element Fields, did not pass Assignment Check");
	}
	ElementField<Vector> result(EF1.p_blockMesh);
    std::vector<Vector> List;
	for (int i = 0; i < EF1.p_blockMesh->n_elemNum; i++)
	{
		List.push_back(EF1.v_value[i] - EF2.v_value[i]);
	}
	result.Initialize(List);
	return result;
}

template<>
ElementField<Scalar> operator * (const Scalar coef, const ElementField<Scalar>& phi)
{
	if (false == phi.Assignment())
	{
		FatalError("Can not proceed operation * between a Scalar and a ElementField, did not pass Assignment Check");
	}
	ElementField<Scalar> result(phi);
	for (int i = 0; i < phi.p_blockMesh->n_elemNum; i++)
	{
		result.SetValue(i, coef*phi.GetValue(i));
	}
	return result;
}

template<>
ElementField<Vector> operator * (const Scalar coef, const ElementField<Vector>& phi)
{
	if (false == phi.Assignment())
	{
		FatalError("Can not proceed operation * between a Scalar and a ElementField, did not pass Assignment Check");
	}
	ElementField<Vector> result(phi);
	for (int i = 0; i < phi.p_blockMesh->n_elemNum; i++)
	{
		result.SetValue(i, coef*phi.GetValue(i));
	}
	return result;
}

// "*"overload
template<>
ElementField<Scalar> operator * (const ElementField<Scalar>& EF1, const ElementField<Scalar>& EF2)
{
	if (EF1.p_blockMesh != EF2.p_blockMesh)
	{
		FatalError("in operation * between two Element Fields, the grids are not consistent");
	}
	if (false == EF1.Assignment() || false == EF2.Assignment())
	{
		FatalError("Can not proceed operation * between two Element Fields, did not pass Assignment Check");
	}
	ElementField<Scalar> result(EF1.p_blockMesh);
    std::vector<Scalar> scalarList;
	for (int i = 0; i < EF1.p_blockMesh->n_elemNum; i++)
	{
		scalarList.push_back(EF1.v_value[i] * EF2.v_value[i]);
	}
	result.Initialize(scalarList);
	return result;
}

template<>
ElementField<Vector> operator * (const ElementField<Scalar>& EF1, const ElementField<Vector>& EF2)
{
	if (EF1.p_blockMesh != EF2.p_blockMesh)
	{
		FatalError("in operation * between a Scalar Field and a Vector one, the grids are not consistent");
	}
	if (false == EF1.Assignment() || false == EF2.Assignment())
	{
		FatalError("Can not proceed operation * between a Scalar Field and a Vector one for not pass Assignment Check");
	}
	ElementField<Vector> result(EF1.p_blockMesh);
    std::vector<Vector> vectorList;
	for (int i = 0; i < EF1.p_blockMesh->n_elemNum; i++)
	{
		vectorList.push_back(EF1.v_value[i] * EF2.v_value[i]);
	}
	result.Initialize(vectorList);
	return result;
}

template<>
ElementField<Scalar> operator / (const ElementField<Scalar>& phi1, const ElementField<Scalar>& phi2)
{
	if (phi1.p_blockMesh != phi2.p_blockMesh)
	{
		FatalError("in operation / between two Scalar Fields, the grids are not consistent");
	}
	if (false == phi1.Assignment() || false == phi2.Assignment())
	{
		FatalError("Can not proceed operation / between two Scalar Fields for not pass Assignment Check");
	}
	ElementField<Scalar> result(phi1.p_blockMesh);
	result.Initialize();
	for (int i = 0; i < phi1.p_blockMesh->n_elemNum; i++)
	{
		result.SetValue(i, phi1.v_value[i] / phi2.v_value[i]);
	}
	return result;
}

template<>
ElementField<Vector> operator / (const ElementField<Vector>& phi1, const ElementField<Scalar>& phi2)
{
	if (phi1.p_blockMesh != phi2.p_blockMesh)
	{
		FatalError("in operation / between two Scalar Fields, the grids are not consistent");
	}
	if (false == phi1.Assignment() || false == phi2.Assignment())
	{
		FatalError("Can not proceed operation / between two Scalar Fields for not pass Assignment Check");
	}
	ElementField<Vector> result(phi1.p_blockMesh);
	result.Initialize();
	for (int i = 0; i < phi1.p_blockMesh->n_elemNum; i++)
	{
		result.SetValue(i, phi1.v_value[i] / phi2.v_value[i]);
	}
	return result;
}

ElementField<Scalar> operator / (Scalar alpha, const ElementField<Scalar>& phi)
{
	if (false == phi.Assignment())
	{
		FatalError("Can not proceed operation / between a Scalar and a Scalar ElementField");
	}
	ElementField<Scalar> result(phi.p_blockMesh);
	result.Initialize();
	for (int i = 0; i < phi.p_blockMesh->n_elemNum; i++)
	{
		result.SetValue(i, alpha / phi.GetValue(i));
	}
	return result;
}

ElementField<Scalar> operator & (const ElementField<Vector>& phi1, const ElementField<Vector>& phi2)
{
	if (phi1.p_blockMesh != phi2.p_blockMesh)
	{
		FatalError("in operation & between two Vector fields, the grids are not consistent");
	}
	if (false == phi1.Assignment() || false == phi1.Assignment())
	{
		FatalError("Can not proceed operation & between two Vector fields for not pass Assignment Check");
	}
	ElementField<Scalar> result(phi1.p_blockMesh);
	result.Initialize();
	for (int i = 0; i < phi1.p_blockMesh->n_elemNum; i++)
	{
		result.SetValue(i, phi1.GetValue(i) & phi2.GetValue(i));
	}
	return result;
}
