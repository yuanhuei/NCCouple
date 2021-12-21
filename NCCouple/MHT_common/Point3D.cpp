#include "Point3D.h"

// =======================================Scalar Type BaseField=============================================

template<>
Scalar Point3D<Scalar>::Mag() const
{
	return sqrt(x_*x_ + y_*y_ + z_*z_);
}

template<>
Point3D<Scalar>& Point3D<Scalar>::Normalize()
{
	Scalar length = sqrt(x_*x_ + y_*y_ + z_*z_);
	x_ = x_ / length;
	y_ = y_ / length;
	z_ = z_ / length;

	return *this;
}

template<>
Point3D<Scalar> Point3D<Scalar>::GetNormal()
{
	return (*this) / this->Mag();
}

template<>
void Point3D<Scalar>::Projection(Point3D<Scalar>& facenorm)
{
	Point3D<Scalar> norm = facenorm.GetNormal();
	Scalar d = norm.x_*x_ + norm.y_*y_ + norm.z_*z_;
	x_ = x_ - d*norm.x_;
	y_ = y_ - d*norm.y_;
	z_ = z_ - d*norm.z_;
}

//"+"overload
template<>
Point3D<Scalar> Point3D<Scalar>::operator + (const Point3D<Scalar>& rhs) const
{
	return Point3D<Scalar>(x_ + rhs.x_, y_ + rhs.y_, z_ + rhs.z_);
}

//"-"overload
template<>
Point3D<Scalar> Point3D<Scalar>::operator - (const Point3D<Scalar>& rhs) const
{
	return Point3D<Scalar>(x_ - rhs.x_, y_ - rhs.y_, z_ - rhs.z_);
}

//"-"overload
template<>
Point3D<Scalar> Point3D<Scalar>::operator - ()
{
	return Point3D<Scalar>(-x_, -y_, -z_);
}

//"+="overload
//modified by Kong Ling, "rhs==>rhs.x_", at 20160215
template<>
void Point3D<Scalar>::operator += (const Point3D<Scalar>& rhs)
{
	this->x_ += rhs.x_;
	this->y_ += rhs.y_;
	this->z_ += rhs.z_;
}

//"-="overload
//modified by Kong Ling, "rhs==>rhs.x_", at 20160215
template<>
void Point3D<Scalar>::operator -= (const Point3D<Scalar>& rhs)
{
	this->x_ -= rhs.x_;
	this->y_ -= rhs.y_;
	this->z_ -= rhs.z_;
}

//"*="overload
template<>
void Point3D<Scalar>::operator *= (const Scalar rhs)
{
	this->x_ *= rhs;
	this->y_ *= rhs;
	this->z_ *= rhs;
}

//"/="overload
template<>
void Point3D<Scalar>::operator /= (const Scalar rhs)
{
	this->x_ /= rhs;
	this->y_ /= rhs;
	this->z_ /= rhs;
}

//"=="overload
template<>
bool Point3D<Scalar>::operator == (const Point3D& rhs)
{
	return SMALL > (*this - rhs).Mag() ? true : false;
}

