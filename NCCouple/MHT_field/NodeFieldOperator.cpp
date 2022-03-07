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

template<>
NodeField<Scalar> operator + (const NodeField<Scalar>& EF1, const NodeField<Scalar>& EF2)
{
	if (EF1.p_blockMesh != EF2.p_blockMesh)
	{
		FatalError("in operation + between two Element Fields, the grids are not consistent");
	}
	if (false == EF1.Assignment() || false == EF2.Assignment())
	{
		FatalError("Can not proceed operation + between two Element Fields, did not pass Assignment Check");
	}
	NodeField<Scalar> result(EF1.p_blockMesh);
    std::vector<Scalar> List;
	for (int i = 0; i < EF1.p_blockMesh->n_nodeNum; i++)
	{
		List.push_back(EF1.v_value[i] + EF2.v_value[i]);
	}
	result.Initialize(List);
	return result;
}

template<>
NodeField<Vector> operator + (const NodeField<Vector>& EF1, const NodeField<Vector>& EF2)
{
	if (EF1.p_blockMesh != EF2.p_blockMesh)
	{
		FatalError("in operation + between two Element Fields, the grids are not consistent");
	}
	if (false == EF1.Assignment() || false == EF2.Assignment())
	{
		FatalError("Can not proceed operation + between two Element Fields, did not pass Assignment Check");
	}
	NodeField<Vector> result(EF1.p_blockMesh);
    std::vector<Vector> List;
	for (int i = 0; i < EF1.p_blockMesh->n_nodeNum; i++)
	{
		List.push_back(EF1.v_value[i] + EF2.v_value[i]);
	}
	result.Initialize(List);
	return result;
}

template<>
NodeField<Scalar> operator - (const NodeField<Scalar>& NF)
{
	if (false == NF.Assignment())
	{
		FatalError("Can not proceed operation - for a field, did not pass Assignment Check");
	}
	NodeField<Scalar> result(NF.p_blockMesh);
	std::vector<Scalar> List;
	for (int i = 0; i < NF.p_blockMesh->n_nodeNum; i++)
	{
		List.push_back(-NF.v_value[i]);
	}
	result.Initialize(List);
	return result;
}

template<>
NodeField<Vector> operator - (const NodeField<Vector>& NF)
{
	if (false == NF.Assignment())
	{
		FatalError("Can not proceed operation - for a field, did not pass Assignment Check");
	}
	NodeField<Vector> result(NF.p_blockMesh);
	std::vector<Vector> List;
	for (int i = 0; i < NF.p_blockMesh->n_nodeNum; i++)
	{
		List.push_back(Vector(0.0, 0.0, 0.0) - NF.v_value[i]);
	}
	result.Initialize(List);
	return result;
}

template<>
NodeField<Scalar> operator - (const NodeField<Scalar>& EF1, const NodeField<Scalar>& EF2)
{
	if (EF1.p_blockMesh != EF2.p_blockMesh)
	{
		FatalError("in operation - between two Element Fields, the grids are not consistent");
	}
	if (false == EF1.Assignment() || false == EF2.Assignment())
	{
		FatalError("Can not proceed operation - between two Element Fields, did not pass Assignment Check");
	}
	NodeField<Scalar> result(EF1.p_blockMesh);
    std::vector<Scalar> List;
	for (int i = 0; i < EF1.p_blockMesh->n_nodeNum; i++)
	{
		List.push_back(EF1.v_value[i] - EF2.v_value[i]);
	}
	result.Initialize(List);
	return result;
}

template<>
NodeField<Vector> operator - (const NodeField<Vector>& EF1, const NodeField<Vector>& EF2)
{
	if (EF1.p_blockMesh != EF2.p_blockMesh)
	{
		FatalError("in operation - between two Element Fields, the grids are not consistent");
	}
	if (false == EF1.Assignment() || false == EF2.Assignment())
	{
		FatalError("Can not proceed operation - between two Element Fields, did not pass Assignment Check");
	}
	NodeField<Vector> result(EF1.p_blockMesh);
    std::vector<Vector> List;
	for (int i = 0; i < EF1.p_blockMesh->n_nodeNum; i++)
	{
		List.push_back(EF1.v_value[i] - EF2.v_value[i]);
	}
	result.Initialize(List);
	return result;
}

//=======================External Functions===========================
template<>
NodeField<Scalar> operator * (const Scalar coef, const NodeField<Scalar>& phi)
{
	if (false == phi.Assignment())
	{
		FatalError("Can not proceed operation * between a Scalar and a nodeField, did not pass Assignment Check");
	}
	NodeField<Scalar> result(phi);
	for (int i = 0; i < phi.p_blockMesh->n_nodeNum; i++)
	{
		result.SetValue(i, coef*phi.GetValue(i));
	}
	return result;
}

template<>
NodeField<Vector> operator * (const Scalar coef, const NodeField<Vector>& phi)
{
	if (false == phi.Assignment())
	{
		FatalError("Can not proceed operation * between a Scalar and a nodeField, did not pass Assignment Check");
	}
	NodeField<Vector> result(phi);
	for (int i = 0; i < phi.p_blockMesh->n_nodeNum; i++)
	{
		result.SetValue(i, coef*phi.GetValue(i));
	}
	return result;
}

// "*"overload
template<>
NodeField<Scalar> operator * (const NodeField<Scalar>& NF1, const NodeField<Scalar>& NF2)
{
	if (NF1.p_blockMesh != NF2.p_blockMesh)
	{
		FatalError("in operation * between two Node Fields, the grids are not consistent");
	}
	if (false == NF1.Assignment() || false == NF2.Assignment())
	{
		FatalError("Can not proceed operation * between two Node Fields, did not pass Assignment Check");
	}
	NodeField<Scalar> result(NF1.p_blockMesh);
    std::vector<Scalar> scalarList;
	for (int i = 0; i < NF1.p_blockMesh->n_nodeNum; i++)
	{
		scalarList.push_back(NF1.v_value[i] * NF2.v_value[i]);
	}
	result.Initialize(scalarList);
	return result;
}

template<>
NodeField<Vector> operator * (const NodeField<Scalar>&EF1, const NodeField<Vector>&EF2)
{
	if (EF1.p_blockMesh != EF2.p_blockMesh)
	{
		FatalError("in operation * between a Scalar Field and a Vector one, the grids are not consistent");
	}
	if (false == EF1.Assignment() || false == EF2.Assignment())
	{
		FatalError("Can not proceed operation * between a Scalar Field and a Vector one for not pass Assignment Check");
	}
	NodeField<Vector> result(EF1.p_blockMesh);
    std::vector<Vector> vectorList;
	for (int i = 0; i < EF1.p_blockMesh->n_nodeNum; i++)
	{
		vectorList.push_back(EF1.v_value[i] * EF2.v_value[i]);
	}
	result.Initialize(vectorList);
	return result;
}

template<>
NodeField<Scalar> operator / (const NodeField<Scalar>&phi1, const NodeField<Scalar>&phi2)
{
	if (phi1.p_blockMesh != phi2.p_blockMesh)
	{
		FatalError("in operation / between two Scalar Fields, the grids are not consistent");
	}
	if (false == phi1.Assignment() || false == phi2.Assignment())
	{
		FatalError("Can not proceed operation / between two Scalar Fields for not pass Assignment Check");
	}
	NodeField<Scalar> result(phi1.p_blockMesh);
	result.Initialize();
	for (int i = 0; i < phi1.p_blockMesh->n_nodeNum; i++)
	{
		result.SetValue(i, phi1.v_value[i] / phi2.v_value[i]);
	}
	return result;
}

template<>
NodeField<Vector> operator / (const NodeField<Vector>&phi1, const NodeField<Scalar>&phi2)
{
	if (phi1.p_blockMesh != phi2.p_blockMesh)
	{
		FatalError("in operation / between two Scalar Fields, the grids are not consistent");
	}
	if (false == phi1.Assignment() || false == phi2.Assignment())
	{
		FatalError("Can not proceed operation / between two Scalar Fields for not pass Assignment Check");
	}
	NodeField<Vector> result(phi1.p_blockMesh);
	result.Initialize();
	for (int i = 0; i < phi1.p_blockMesh->n_nodeNum; i++)
	{
		result.SetValue(i, phi1.v_value[i] / phi2.v_value[i]);
	}
	return result;
}

NodeField<Scalar> operator / (Scalar alpha, const NodeField<Scalar>& phi)
{
	if (false == phi.Assignment())
	{
		FatalError("Can not proceed operation / between a Scalar and a Scalar NodeField");
	}
	NodeField<Scalar> result(phi.p_blockMesh);
	result.Initialize();
	for (int i = 0; i < phi.p_blockMesh->n_nodeNum; i++)
	{
		result.SetValue(i, alpha / phi.GetValue(i));
	}
	return result;
}

NodeField<Scalar> operator & (const NodeField<Vector>& phi1, const NodeField<Vector>& phi2)
{
	if (phi1.p_blockMesh != phi2.p_blockMesh)
	{
		FatalError("in operation & between two Vector fields, the grids are not consistent");
	}
	if (false == phi1.Assignment() || false == phi1.Assignment())
	{
		FatalError("Can not proceed operation & between two Vector fields for not pass Assignment Check");
	}
	NodeField<Scalar> result(phi1.p_blockMesh);
	result.Initialize();
	for (int i = 0; i < phi1.p_blockMesh->n_nodeNum; i++)
	{
		result.SetValue(i, phi1.GetValue(i) & phi2.GetValue(i));
	}
	return result;
}
