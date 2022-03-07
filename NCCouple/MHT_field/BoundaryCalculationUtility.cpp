/*---------------------------------------------------------------------------*\
File Name:
    BoundaryCalculationUtility.cpp

Description:
    This file defines utilities for field boundary calculation

    Author:		Kong Ling
    Date: 2018-10-07
\*---------------------------------------------------------------------------*/


#include "../MHT_field/BoundaryCalculationUtility.h"
#include "../MHT_common/SystemControl.h"
#include "../MHT_field/FaceField.h"
#include "../MHT_field/ElementField.h"

int ElementToFaceParameterCheck(std::string input)
{
	if ("linear" == input)
	{
		return 0;
	}
	else if ("hamonic" == input)
	{
		return 0;
	}
	else if ("default" == input)
	{
		return 1;
	}
	else if ("copy" == input)
	{
		return 1;
	}
	else if ("extrapolate" == input)
	{
		return 1;
	}
	else if ("boundaryCondition" == input)
	{
		return 1;
	}
	else if ("boundaryLeft" == input)
	{
		return 1;
	}
	else
	{
		FatalError("unsupported parameter " + input + "found in ElementToFace");
		return -1;
	}
}

template<>
void ElementToFaceCheck(ElementField<Scalar>& EF, FaceField<Scalar>& FF)
{
	//step 1: check status of element field, which is interpolated from and therefore required to be assigned
	if (false == EF.Assignment())
	{
		FatalError("Cannot proceed ElementToFace for not passing Assignment Check");
	}
	//step 2: create face field if face field is not existing
	if (false == FF.Existence())
	{
		FF.Initialize();
	}
	return;
}

template<>
void ElementToFaceCheck(ElementField<Vector>& EF, FaceField<Vector>& FF)
{
	//step 1: check status of element field, which is interpolated from and therefore required to be assigned
	if (false == EF.Assignment())
	{
		FatalError("Cannot proceed ElementToFace for not passing Assignment Check");
	}
	//step 2: create face field if face field is not existing
	if (false == FF.Existence())
	{
		FF.Initialize();
	}
	return;
}

