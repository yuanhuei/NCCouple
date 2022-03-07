/*---------------------------------------------------------------------------*\
File Name:
	ElementFieldOperator.cpp

Description:
    1. This file defines operators for fields

	Author:		Kong Ling
    Date:   2018-10-02

	
Revised:
	Description:
	1. Convert "non const" input parameters to const;
	2. Convert "non const" function to const function;
	
	Revisor:		Shuai Zhang
	Modified Date:	2018-11-28
\*---------------------------------------------------------------------------*/

#include "../MHT_field/FieldOperator.h"
#include "../MHT_common/SystemControl.h"

//========================Field operators===========================
//==========Field<Type> operator + (Field<Type>&, Field<Type>&)=============
template<>
Field<Scalar> operator + (const Field<Scalar>& phi1, const Field<Scalar>& phi2)
{
	Field<Scalar> resultField(phi1.p_blockMesh);
	if (phi1.p_blockMesh != phi2.p_blockMesh)
	{
		FatalError("in operator + between two Scalar fields, the grids are not consistent");
	}
	if (true == phi1.elementField.Assignment() && true == phi2.elementField.Assignment())
	{
		resultField.elementField = phi1.elementField + phi2.elementField;
	}
	if (true == phi1.faceField.Assignment() && true == phi2.faceField.Assignment())
	{
		resultField.faceField = phi1.faceField + phi2.faceField;
	}
	if (true == phi1.nodeField.Assignment() && true == phi2.nodeField.Assignment())
	{
		resultField.nodeField = phi1.nodeField + phi2.nodeField;
	}
	return resultField;
}

template<>
Field<Vector> operator + (const Field<Vector>& phi1, const Field<Vector>& phi2)
{
	Field<Vector> resultField(phi1.p_blockMesh);
	if (phi1.p_blockMesh != phi2.p_blockMesh)
	{
		FatalError("in operator + between two Vector fields, the grids are not consistent");
	}
	if (true == phi1.elementField.Assignment() && true == phi2.elementField.Assignment())
	{
		resultField.elementField = phi1.elementField + phi2.elementField;
	}
	if (true == phi1.faceField.Assignment() && true == phi2.faceField.Assignment())
	{
		resultField.faceField = phi1.faceField + phi2.faceField;
	}
	if (true == phi1.nodeField.Assignment() && true == phi2.nodeField.Assignment())
	{
		resultField.nodeField = phi1.nodeField + phi2.nodeField;
	}
	return resultField;
}

//==========Field<Type> operator - (Field<Type>&)=============
template<>
Field<Scalar> operator - (const Field<Scalar>& phi1)
{
	Field<Scalar> resultField(phi1.p_blockMesh);
	if (true == phi1.elementField.Assignment())
	{
		resultField.elementField = -phi1.elementField;
	}
	if (true == phi1.faceField.Assignment())
	{
		resultField.faceField = -phi1.faceField;
	}
	if (true == phi1.nodeField.Assignment())
	{
		resultField.nodeField = -phi1.nodeField;
	}
	return resultField;
}

template<>
Field<Vector> operator - (const Field<Vector>& phi1)
{
	Field<Vector> resultField(phi1.p_blockMesh);
	if (true == phi1.elementField.Assignment())
	{
		resultField.elementField = -phi1.elementField;
	}
	if (true == phi1.faceField.Assignment())
	{
		resultField.faceField = -phi1.faceField;
	}
	if (true == phi1.nodeField.Assignment())
	{
		resultField.nodeField = -phi1.nodeField;
	}
	return resultField;
}

//==========Field<Type> operator - (Field<Type>&, Field<Type>&)=============
template<>
Field<Scalar> operator - (const Field<Scalar>& phi1, const Field<Scalar>& phi2)
{
	Field<Scalar> resultField(phi1.p_blockMesh);
	if (phi1.p_blockMesh != phi2.p_blockMesh)
	{
		FatalError("in operator - between two Scalar fields, the grids are not consistent");
	}
	if (true == phi1.elementField.Assignment() && true == phi2.elementField.Assignment())
	{
		resultField.elementField = phi1.elementField - phi2.elementField;
	}
	if (true == phi1.faceField.Assignment() && true == phi2.faceField.Assignment())
	{
		resultField.faceField = phi1.faceField - phi2.faceField;
	}
	if (true == phi1.nodeField.Assignment() && true == phi2.nodeField.Assignment())
	{
		resultField.nodeField = phi1.nodeField - phi2.nodeField;
	}
	return resultField;
}

template<>
Field<Vector> operator - (const Field<Vector>& phi1, const Field<Vector>& phi2)
{
	Field<Vector> resultField(phi1.p_blockMesh);
	if (phi1.p_blockMesh != phi2.p_blockMesh)
	{
		FatalError("in operator - between two Vector fields, the grids are not consistent");
	}
	if (true == phi1.elementField.Assignment() && true == phi2.elementField.Assignment())
	{
		resultField.elementField = phi1.elementField - phi2.elementField;
	}
	if (true == phi1.faceField.Assignment() && true == phi2.faceField.Assignment())
	{
		resultField.faceField = phi1.faceField - phi2.faceField;
	}
	if (true == phi1.nodeField.Assignment() && true == phi2.nodeField.Assignment())
	{
		resultField.nodeField = phi1.nodeField - phi2.nodeField;
	}
	return resultField;
}

//==========Field<Type> operator * (const Scalar, Field<Type>&)=============
template<>
Field<Scalar> operator * (const Scalar coeff, const Field<Scalar>& phi)
{
	Field<Scalar> resultField(phi);
	if (true == phi.elementField.Assignment())
	{
		resultField.elementField = coeff*resultField.elementField;
	}
	if (true == phi.faceField.Assignment())
	{
		resultField.faceField = coeff*resultField.faceField;
	}
	if (true == phi.nodeField.Assignment())
	{
		resultField.nodeField = coeff*resultField.nodeField;
	}
	return resultField;
}

template<>
Field<Vector> operator * (const Scalar coeff, const Field<Vector>& phi)
{
	Field<Vector> resultField(phi);
	if (true == phi.elementField.Assignment())
	{
		resultField.elementField = coeff*resultField.elementField;
	}
	if (true == phi.faceField.Assignment())
	{
		resultField.faceField = coeff*resultField.faceField;
	}
	if (true == phi.nodeField.Assignment())
	{
		resultField.nodeField = coeff*resultField.nodeField;
	}
	return resultField;
}

//==========Field<Type> operator * (Field<Scalar>&, Field<Type>&)==========
template<>
Field<Scalar> operator * (const Field<Scalar>& phi1, const Field<Scalar>& phi2)
{
	Field<Scalar> resultField(phi1.p_blockMesh);
	if (phi1.p_blockMesh != phi2.p_blockMesh)
	{
		FatalError("in operator * between two Scalar fields, the grids are not consistent");
	}
	if (true == phi1.elementField.Assignment() && true == phi2.elementField.Assignment())
	{
		resultField.elementField = phi1.elementField * phi2.elementField;
	}
	if (true == phi1.faceField.Assignment() && true == phi2.faceField.Assignment())
	{
		resultField.faceField = phi1.faceField * phi2.faceField;
	}
	if (true == phi1.nodeField.Assignment() && true == phi2.nodeField.Assignment())
	{
		resultField.nodeField = phi1.nodeField * phi2.nodeField;
	}
	return resultField;
}

template<>
Field<Vector> operator * (const Field<Scalar>& phi1, const Field<Vector>& phi2)
{
	Field<Vector> resultField(phi1.p_blockMesh);
	if (phi1.p_blockMesh != phi2.p_blockMesh)
	{
		FatalError("in operator * between Scalar and Vector fields, the grids are not consistent");
	}
	if (true == phi1.elementField.Assignment() && true == phi2.elementField.Assignment())
	{
		resultField.elementField = phi1.elementField * phi2.elementField;
	}
	if (true == phi1.faceField.Assignment() && true == phi2.faceField.Assignment())
	{
		resultField.faceField = phi1.faceField * phi2.faceField;
	}
	if (true == phi1.nodeField.Assignment() && true == phi2.nodeField.Assignment())
	{
		resultField.nodeField = phi1.nodeField * phi2.nodeField;
	}
	return resultField;
}

//==============Field<Type> operator / (Field<Type>&, Field<Scalar>&)==========
template<>
Field<Scalar> operator / (const Field<Scalar>& phi1, const Field<Scalar>& phi2)
{
	Field<Scalar> resultField(phi1.p_blockMesh);
	if (phi1.p_blockMesh != phi2.p_blockMesh)
	{
		FatalError("in operator / between two Scalar fields, the grids are not consistent");
	}
	if (true == phi1.elementField.Assignment() && true == phi2.elementField.Assignment())
	{
		resultField.elementField = phi1.elementField / phi2.elementField;
	}
	if (true == phi1.faceField.Assignment() && true == phi2.faceField.Assignment())
	{
		resultField.faceField = phi1.faceField / phi2.faceField;
	}
	if (true == phi1.nodeField.Assignment() && true == phi2.nodeField.Assignment())
	{
		resultField.nodeField = phi1.nodeField / phi2.nodeField;
	}
	return resultField;
}

Field<Scalar> operator / (Scalar alpha, const Field<Scalar>& phi)
{
	Field<Scalar> resultField(phi.p_blockMesh);
	if (true == phi.elementField.Assignment())
	{
		resultField.elementField = alpha / phi.elementField;
	}
	if (true == phi.faceField.Assignment())
	{
		resultField.faceField = alpha / phi.faceField;
	}
	if (true == phi.nodeField.Assignment())
	{
		resultField.nodeField = alpha / phi.nodeField;
	}
	return resultField;
}

template<>
Field<Vector> operator / (const Field<Vector>& phi1, const Field<Scalar>& phi2)
{
	Field<Vector> resultField(phi1.p_blockMesh);
	if (phi1.p_blockMesh != phi2.p_blockMesh)
	{
		FatalError("in operator / between two Scalar fields, the grids are not consistent");
	}
	if (true == phi1.elementField.Assignment() && true == phi2.elementField.Assignment())
	{
		resultField.elementField = phi1.elementField / phi2.elementField;
	}
	if (true == phi1.faceField.Assignment() && true == phi2.faceField.Assignment())
	{
		resultField.faceField = phi1.faceField / phi2.faceField;
	}
	if (true == phi1.nodeField.Assignment() && true == phi2.nodeField.Assignment())
	{
		resultField.nodeField = phi1.nodeField / phi2.nodeField;
	}
	return resultField;
}

//==============Field<Scalar> operator & (Field<Vector>&, Field<Vector>&)==========
Field<Scalar> operator & (const Field<Vector>& phi1, const Field<Vector>& phi2)
{
	Field<Scalar> resultField(phi1.p_blockMesh);
	if (phi1.p_blockMesh != phi2.p_blockMesh)
	{
		FatalError("in operator & between two Vector fields, the grids are not consistent");
	}
	if (true == phi1.elementField.Assignment() && true == phi2.elementField.Assignment())
	{
		resultField.elementField = phi1.elementField & phi2.elementField;
	}
	if (true == phi1.faceField.Assignment() && true == phi2.faceField.Assignment())
	{
		resultField.faceField = phi1.faceField & phi2.faceField;
	}
	if (true == phi1.nodeField.Assignment() && true == phi2.nodeField.Assignment())
	{
		resultField.nodeField = phi1.nodeField & phi2.nodeField;
	}
	return resultField;
}
