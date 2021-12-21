/*---------------------------------------------------------------------------*\
File Name:
    NumerBC.cpp

Description:
    1. A general template for one-side boundary condition using 3 parameters
    a*phi + b* \partial phi/\partial n = c;

    Author:		Kong Ling
    Date:   2018-10-03
\*---------------------------------------------------------------------------*/


#include "../MHT_field/NumerBC.h"
#include "../MHT_common/Vector.h"
#include "../MHT_common/Tensor.h"

// ==========Scalar Type NumerBC=====================
//constructor
template<>
NumerBC<Scalar>::NumerBC()
	:a(0.0), b(0.0), c(0.0), con1(1.0), con2(0.0)
{}

// ==========Vector Type NumerBC====================
//constructor
template<>
NumerBC<Vector>::NumerBC()
	:a(0.0), b(0.0), c(Vector(0.0, 0.0, 0.0)), con1(1.0), con2(Vector(0.0, 0.0, 0.0))
{}

// ==========Scalar Type NumerBC=====================
//constructor
template<>
NumerBC<Scalar>::NumerBC(Scalar A, Scalar B, Scalar C)
	:a(A), b(B), c(C)
{
	this->con1 = 0.0;
	this->con2 = 0.0;
}

// ==========Vector Type NumerBC====================
//constructor
template<>
NumerBC<Vector>::NumerBC(Scalar A, Scalar B, Vector C)
	:a(A), b(B), c(C)
{
	this->con1 = 0.0;
	this->con2 = Vector(0.0, 0.0, 0.0);
}