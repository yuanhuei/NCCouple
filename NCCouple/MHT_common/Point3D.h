/*---------------------------------------------------------------------------*\
Class:
	Point3D
File Name:
	Point3D.h
Description:
	A basic class describe 3D point.
	1.Construction from 3 components;
	2.Operator including :
	"<<"	Output stream
	"*"		Scalar and vector product, or Vector(Point3D) and scalar product
	"^"		Vector(Point3D) cross product operators
	"/"		Vector(Point3D) devided by scalar, or scalar and Vector(Point3D) product
	"&"		Vector(Point3D) inner-product (dot-product)
	"+"		Vector(Point3D) add Operator
	"-"		Vector(Point3D) subtract Operator
	"+="	Vector(Point3D) add Operator 
	"-="	Vector(Point3D) subtract Operator
	"*="	Vector(Point3D) mutiply Operator
	"/="	Vector(Point3D) divide Operator
	"=="	Whether these two Vector(Point3D) is the same

	Author:		ShuaiZhang
	Revisor:	ZhiHaoJia, KongLing
	Date: 2016-01-21

Revised:	
	Description:	Changed void Normalize();  to Point3D& Normalize();	
	Revisor:		ShuaiZhang
	Modified Date:	2016-09-08
Revised:
	Description:	adding Point3D<Type> GetNormal() as a member function;
					inline friend Point3D<Type> operator / (Scalar, Point3D<Type>) deleted;
	Revisor:		Kong Ling
	Modified Date:	2016-11-22
\*---------------------------------------------------------------------------*/

#pragma   once

#ifndef _Point3D_
#define _Point3D_


#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <ctime>
#include <math.h>

#include "MaxMin.h"
#include "AllocateArray.h"
#include "Configuration.h"

template<class Type>
class Point3D
{
public:
	Type	x_, y_, z_;
	Point3D() :x_(0.0), y_(0.0), z_(0.0){}
	Point3D(Type X, Type Y, Type Z) :x_(X), y_(Y), z_(Z){}

	Scalar Mag() const;
	Point3D& Normalize();
	Point3D GetNormal();
	void Projection(Point3D<Type>& facenorm);

	//==================================friend function=======================================
	//"<<"overload
	inline friend std::ostream& operator << (std::ostream& out, const Point3D<Type> & rhs)
	{
		out << "( " << rhs.x_ << " , " << rhs.y_ << " , " << rhs.z_ << " ) ";
		return out;
	}

	//Scalar multiply Vector, Vector multiply Scalar
    inline friend Point3D<Type> operator * (const Scalar d, const Point3D<Type>& p)
	{
		return Point3D(d * p.x_, d * p.y_, d * p.z_);
	}
    inline friend Point3D<Type> operator * (const Point3D<Type>& p, const Scalar d)
	{
		return Point3D(p.x_ * d, p.y_ * d, p.z_ * d);
	}

	//cross product"^"overload
	inline friend Point3D<Type> operator ^ (const Point3D<Type>& v1, const Point3D<Type>& v2)
	{
		return Point3D<Type>
			(
			(v1.y_ * v2.z_ - v1.z_ * v2.y_),
			(v1.z_ * v2.x_ - v1.x_ * v2.z_),
			(v1.x_ * v2.y_ - v1.y_ * v2.x_)
			);
	}

	//"/"overload
	inline friend Point3D<Type> operator / (const Point3D<Type>& p, const Scalar d)
	{
		if (fabs(d) < SMALL) //revised by Kong Ling, 2016/02/02
		{
			std::cout << "Warning:\n" << "Divisor 'd' is smaller than 'SMALL'.\n" << "Please Check The Divisor!"<< std::endl;
		}
		return Point3D<Type>(p.x_ / d, p.y_ / d, p.z_ / d);
	}

	//"&" overload, inner product
    inline friend Scalar operator & (const Point3D<Type>& lhs, const Point3D<Type>& rhs)
	{
		return (Scalar)(lhs.x_ * rhs.x_ + lhs.y_ * rhs.y_ + lhs.z_ * rhs.z_);
	}

	//==================================friend function=======================================

public:

	//"+"overload
	Point3D<Type> operator + (const Point3D<Type>& rhs) const;

	//"-"overload
	Point3D<Type> operator - (const Point3D<Type>& rhs) const;

	//"-"negtive sign overload
	Point3D<Type> operator - ();

	//"+="overload
	void operator += (const Point3D<Type>& rhs);

	//"-="overload
	void operator -= (const Point3D<Type>& rhs);

	//"*="overload
    void operator *= (const Scalar rhs);

	//"/="overload
    void operator /= (const Scalar rhs);

	//"=="overload
	bool operator == (const Point3D<Type>& rhs);
};


#endif
