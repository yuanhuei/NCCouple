
#ifndef _Tensor_
#define _Tensor_

#include "Vector.h"

class Tensor
{
public:
	Scalar	mat_[3][3];

	Tensor();
	Tensor(Scalar, Scalar, Scalar, Scalar, Scalar, Scalar, Scalar, Scalar, Scalar);
    ~Tensor(){}

    Scalar norm() const;
    Vector Row(int) const;
    Vector Column(int) const;
    Tensor Transpose() const;
    Scalar Det() const;
    Tensor Accompany() const;
    Tensor Inverse() const;
	
	
	Scalar Trace() const;
	Tensor Symm() const;
	Tensor AntiSymm() const;
	Tensor Dev() const;
	Tensor Dev2() const;

	//Operator
	Tensor operator + (const Tensor&) const;
	Tensor& operator += (const Tensor&);
	Tensor operator - (const Tensor&) const;
	Tensor& operator -= (const Tensor&);
	Tensor operator * (const Scalar) const;
	Vector operator * (const Vector) const;
	Tensor operator * (const Tensor&) const;
	Tensor operator / (const Scalar) const;
	Scalar operator && (const Tensor&) const;

	inline friend Tensor operator * (const Scalar d, const Tensor& t)
	{
		return t*d;
	}
};

std::ostream& operator << (std::ostream& os, const Tensor& t);
std::istream& operator >> (std::istream& os, Tensor& t);

Vector operator * (const Vector&v, const Tensor& t);
Tensor operator * (const Vector& p1, const Vector& p2);

#endif
