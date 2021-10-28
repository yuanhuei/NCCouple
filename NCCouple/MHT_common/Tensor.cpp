/*---------------------------------------------------------------------------*\
Class:
    Tensor
File Name:
    Tensor.h
Description:
    A basic class for Tensor.

    Author:		ShuaiZhang, KongLing
    Revisor:	ZhiHaoJia
    Date: Long ago

Revised:
    Description:
    1. Convert "non const" input parameters to const;
    2. Convert "non const" function to const function;

    Revisor:		Shuai Zhang
    Modified Date:	2018-12-03
\*---------------------------------------------------------------------------*/

#include "Tensor.h"

#include <cstdlib>
#include <cmath>
#include <iostream>

Tensor::Tensor()
{
	mat_[0][0] = 0.0;
	mat_[0][1] = 0.0;
	mat_[0][2] = 0.0;
	mat_[1][0] = 0.0;
	mat_[1][1] = 0.0;
	mat_[1][2] = 0.0;
	mat_[2][0] = 0.0;
	mat_[2][1] = 0.0;
	mat_[2][2] = 0.0;
}

Tensor::Tensor(Scalar a11, Scalar a12, Scalar a13, Scalar a21, Scalar a22, Scalar a23, Scalar a31, Scalar a32, Scalar a33)
{
	mat_[0][0] = a11;
	mat_[0][1] = a12;
	mat_[0][2] = a13;
	mat_[1][0] = a21;
	mat_[1][1] = a22;
	mat_[1][2] = a23;
	mat_[2][0] = a31;
	mat_[2][1] = a32;
	mat_[2][2] = a33;
}

Scalar Tensor::norm() const
{
	Scalar a = 0.0;
	Scalar sum;
	for (int i = 0; i<3; i++)
	{
		sum = 0.0;
        for (int j = 0; j<3; j++) sum += fabs(mat_[i][j]);
		if (sum > a) a = sum;
	}
	return (a);
}

//Added by Kong Ling, at 20160204
Vector Tensor::Row(int n) const
{
	if (n < 1 || n>3)
	{
        std::cout << "Invalid row number is detected in Tensor::Row" << std::endl;
		exit(0);
	}
	return Vector(mat_[n - 1][0], mat_[n - 1][1], mat_[n - 1][2]);
}

//Added by Kong Ling, at 20160204
Vector Tensor::Column(int n) const
{
	if (n < 1 || n>3)
	{
        std::cout << "Invalid column number is detected in Tensor::Column." << std::endl;
		exit(0);
	}
	return Vector(mat_[0][n - 1], mat_[1][n - 1], mat_[2][n - 1]);
}

//Added by Kong Ling, at 20160215
Tensor Tensor::Transpose() const
{
	return Tensor(mat_[0][0], mat_[1][0], mat_[2][0], mat_[0][1], mat_[1][1], mat_[2][1], mat_[0][2], mat_[1][2], mat_[2][2]);
}

Scalar Tensor::Det() const
{
	return (
		mat_[0][0] * (mat_[1][1] * mat_[2][2] - mat_[2][1] * mat_[1][2])
		- mat_[0][1] * (mat_[1][0] * mat_[2][2] - mat_[2][0] * mat_[1][2])
		+ mat_[0][2] * (mat_[1][0] * mat_[2][1] - mat_[2][0] * mat_[1][1])
		);
}

Tensor Tensor::Accompany() const
{
	/*
	| a11  a12  a13 |			| M11  M21  M31 |
	matric =	| a21  a22  a23 |	A* =	| M12  M22  M32 |
	| a31  a32  a33 |			| M13  M23  M33 |
	*/

	return Tensor
		(
		(mat_[1][1] * mat_[2][2] - mat_[2][1] * mat_[1][2]),
		-(mat_[0][1] * mat_[2][2] - mat_[2][1] * mat_[0][2]),
		(mat_[0][1] * mat_[1][2] - mat_[1][1] * mat_[0][2]),

		-(mat_[1][0] * mat_[2][2] - mat_[2][0] * mat_[1][2]),
		(mat_[0][0] * mat_[2][2] - mat_[2][0] * mat_[0][2]),
		-(mat_[0][0] * mat_[1][2] - mat_[1][0] * mat_[0][2]),

		(mat_[1][0] * mat_[2][1] - mat_[2][0] * mat_[1][1]),
		-(mat_[0][0] * mat_[2][1] - mat_[2][0] * mat_[0][1]),
		(mat_[0][0] * mat_[1][1] - mat_[1][0] * mat_[0][1])
		);
}

//Added by Siyuan Yang, at 20191018

Scalar Tensor::Trace() const
{	
	return (mat_[0][0] + mat_[1][1] + mat_[2][2]);
}

//Added by Siyuan Yang, at 20191018
Tensor Tensor::Inverse() const
{
	return (*this).Accompany() / (*this).Det();
}

//Added by Siyuan Yang, at 20191018
Tensor Tensor::Symm() const 
{
	return (0.5 * ((*this) + (*this).Transpose()));
	
}
//Added by Siyuan Yang, at 20191018
Tensor Tensor::AntiSymm() const
{
	return (0.5 * ((*this) - (*this).Transpose()));

}
//Added by Siyuan Yang, at 20191018
Tensor Tensor::Dev() const
{
	Tensor I(1.0, 0.0, 0.0,
		0.0, 1.0, 0.0,
		0.0, 0.0, 1.0);

	return ((*this)- 1.0 / 3.0*(*this).Trace()*I);

}

//Added by Siyuan Yang, at 20191026
Tensor Tensor::Dev2() const
{
	Tensor I(1.0, 0.0, 0.0,
		0.0, 1.0, 0.0,
		0.0, 0.0, 1.0);

	return ((*this) - 2.0 / 3.0*(*this).Trace()*I);

}


Tensor Tensor::operator+(const Tensor& t) const
{
	Tensor t1(*this);
	for (int i = 0; i<3; i++)
		for (int j = 0; j<3; j++)
			t1.mat_[i][j] += t.mat_[i][j];
	return t1;
}

Tensor& Tensor::operator+=(const Tensor& t)
{
	for (int i = 0; i<3; i++)
		for (int j = 0; j<3; j++)
			mat_[i][j] += t.mat_[i][j];
	return *this;
}

Tensor Tensor::operator-(const Tensor& t) const
{
	Tensor t1(*this);
	for (int i = 0; i<3; i++)
		for (int j = 0; j<3; j++)
			t1.mat_[i][j] -= t.mat_[i][j];
	return t1;
}

Tensor& Tensor::operator-=(const Tensor& t)
{
	for (int i = 0; i<3; i++)
		for (int j = 0; j<3; j++)
			mat_[i][j] -= t.mat_[i][j];
	return *this;
}

Tensor Tensor::operator*(const Scalar s) const
{
	Tensor t1(*this);
	for (int i = 0; i<3; i++)
		for (int j = 0; j<3; j++)
			t1.mat_[i][j] *= s;
	return t1;
}

Vector Tensor::operator*(const Vector v) const
{
	return Vector
		(
		v.x_*mat_[0][0] + v.y_*mat_[0][1] + v.z_*mat_[0][2],
		v.x_*mat_[1][0] + v.y_*mat_[1][1] + v.z_*mat_[1][2],
		v.x_*mat_[2][0] + v.y_*mat_[2][1] + v.z_*mat_[2][2]
		);
}

Tensor Tensor::operator*(const Tensor& t) const
{
	return Tensor
		(
		mat_[0][0] * t.mat_[0][0] + mat_[0][1] * t.mat_[1][0] + mat_[0][2] * t.mat_[2][0],
		mat_[0][0] * t.mat_[0][1] + mat_[0][1] * t.mat_[1][1] + mat_[0][2] * t.mat_[2][1],
		mat_[0][0] * t.mat_[0][2] + mat_[0][1] * t.mat_[1][2] + mat_[0][2] * t.mat_[2][2],

		mat_[1][0] * t.mat_[0][0] + mat_[1][1] * t.mat_[1][0] + mat_[1][2] * t.mat_[2][0],
		mat_[1][0] * t.mat_[0][1] + mat_[1][1] * t.mat_[1][1] + mat_[1][2] * t.mat_[2][1],
		mat_[1][0] * t.mat_[0][2] + mat_[1][1] * t.mat_[1][2] + mat_[1][2] * t.mat_[2][2],

		mat_[2][0] * t.mat_[0][0] + mat_[2][1] * t.mat_[1][0] + mat_[2][2] * t.mat_[2][0],
		mat_[2][0] * t.mat_[0][1] + mat_[2][1] * t.mat_[1][1] + mat_[2][2] * t.mat_[2][1],
		mat_[2][0] * t.mat_[0][2] + mat_[2][1] * t.mat_[1][2] + mat_[2][2] * t.mat_[2][2]
		);
}

Tensor Tensor::operator/(const Scalar s) const
{
	Tensor t1(*this);
	for (int i = 0; i<3; i++)
		for (int j = 0; j<3; j++)
			t1.mat_[i][j] /= s;
	return t1;
}

std::ostream& operator<<(std::ostream& os, const Tensor& t)
{
	os << "[\t" << t.mat_[0][0] << "\t" << t.mat_[0][1] << "\t" << t.mat_[0][2] << "\n"
		<< t.mat_[1][0] << "\t" << t.mat_[1][1] << "\t" << t.mat_[1][2] << "\n"
		<< t.mat_[2][0] << "\t" << t.mat_[2][1] << "\t" << t.mat_[2][2] << "\t]" << "\n";
	return os;
}

std::istream& operator>>(std::istream& is, Tensor& t)
{
	t = Tensor();
	is >> t.mat_[0][0] >> t.mat_[0][1] >> t.mat_[0][2]
		>> t.mat_[1][0] >> t.mat_[1][1] >> t.mat_[1][2]
		>> t.mat_[2][0] >> t.mat_[2][1] >> t.mat_[2][2];
	return is;
}

//outer product of two vectors
Tensor operator * (const Vector& p1, const Vector& p2)
{
	return Tensor(p1.x_*p2.x_, p1.x_*p2.y_, p1.x_*p2.z_,
		p1.y_*p2.x_, p1.y_*p2.y_, p1.y_*p2.z_,
		p1.z_*p2.x_, p1.z_*p2.y_, p1.z_*p2.z_);
}

Vector operator * (const Vector&v, const Tensor& t)
{
	return Vector
		(
		v.x_*t.mat_[0][0] + v.y_*t.mat_[1][0] + v.z_*t.mat_[2][0],
		v.x_*t.mat_[0][1] + v.y_*t.mat_[1][1] + v.z_*t.mat_[2][1],
		v.x_*t.mat_[0][2] + v.y_*t.mat_[1][2] + v.z_*t.mat_[2][2]
		);
}

Scalar Tensor::operator && (const Tensor& t1) const
{
	return t1.mat_[0][0] * this->mat_[0][0] + t1.mat_[0][1] * this->mat_[0][1] + t1.mat_[0][2] * this->mat_[0][2] +
		t1.mat_[1][0] * this->mat_[1][0] + t1.mat_[1][1] * this->mat_[1][1] + t1.mat_[1][2] * this->mat_[1][2] +
		t1.mat_[2][0] * this->mat_[2][0] + t1.mat_[2][1] * this->mat_[2][1] + t1.mat_[2][2] * this->mat_[2][2];
}
