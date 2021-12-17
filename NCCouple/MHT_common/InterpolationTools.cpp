/*---------------------------------------------------------------------------*\
File Name:
	InterpolationTools.h

Description:
	Commonly used methods for interpolations

	Author:		Kong Ling
	Date: 2016-10-08

Revised:
    Description:
    1. Convert "non const" input parameters to const;

    Revisor:		Shuai Zhang
    Modified Date:	2018-12-03
\*---------------------------------------------------------------------------*/

#include <iostream>
#include <cmath>
#include <vector>

#include "../MHT_common/AllocateArray.h"
#include "../MHT_common/InterpolationTools.h"
#include "../MHT_common/SystemControl.h"

template<>
Scalar DistanceInterpolation
(
const std::vector<Scalar>& distance,
const std::vector<Scalar>& value
)
{
	//Check existence and and the number of interpolating points
	if (distance.size() != value.size())
	{
		FatalError("in DistancInterpolation(): Numbers of distances and values are inconsistent");
	}
	else if (0 == distance.size())
	{
		FatalError("in DistancInterpolation(): no data is provided for interpolation!");
	}
	
	Scalar smallest = distance[0];
	int n_smallest = 0;
	for (int i = 0; i < (int)distance.size(); i++)
	{
		if (distance[i] < 0)
		{
            std::cout << "Error in DistancInterpolation(): ";
            std::cout << "minus distance found!" << std::endl;
			exit(1);
		}
		if (distance[i] < smallest)
		{
			smallest = distance[i];
			n_smallest = i;
		}
	}

	Scalar numerator = value[n_smallest];
	Scalar dominator=1.0;
	for (int i = 0; i < (int)distance.size(); i++)
	{
		if (i == n_smallest) continue;
		Scalar added;
		if (distance[i] < SMALL)
		{
			added = 1;
		}
		else
		{
			added = smallest / distance[i];
		}
		dominator += added;
		numerator = numerator + added*value[i];
	}
	return numerator / dominator;
}

template<>
Vector DistanceInterpolation
(
const std::vector<Scalar>& distance,
const std::vector<Vector>& value
)
{
	//Check existence and and the number of interpolating points
	if (distance.size() != value.size())
	{
		FatalError("in DistancInterpolation(): Numbers of distances and values are inconsistent");
	}
	else if (0 == distance.size())
	{
		FatalError("in DistancInterpolation(): no data is provided for interpolation!");
	}

	Scalar smallest = distance[0];
	int n_smallest = 0;
	for (unsigned int i = 0; i < distance.size(); i++)
	{
		if (distance[i] < 0)
		{
            std::cout << "Error in DistancInterpolation(): ";
            std::cout << "minus distance found!" << std::endl;
			exit(1);
		}
		if (distance[i] < smallest)
		{
			smallest = distance[i];
			n_smallest = i;
		}
	}

	Vector numerator = value[n_smallest];
	Scalar dominator = 1.0;
	for (int i = 0; i < (int)distance.size(); i++)
	{
		if (i == n_smallest) continue;
		Scalar added;
		if (distance[i] < SMALL)
		{
			added = 1;
		}
		else
		{
			added = smallest / distance[i];
		}
		dominator += added;
		numerator = numerator + added*value[i];
	}
	return numerator / dominator;
}


template<>
Tensor DistanceInterpolation
(
const std::vector<Scalar>& distance,
const std::vector<Tensor>& value
)
{
	//Check existence and and the number of interpolating points
	if (distance.size() != value.size())
	{
		FatalError("in DistancInterpolation(): Numbers of distances and values are inconsistent");
	}
	else if (0 == distance.size())
	{
		FatalError("in DistancInterpolation(): no data is provided for interpolation!");
	}

	Scalar smallest = distance[0];
	int n_smallest = 0;
	for (int i = 0; i < (int)distance.size(); i++)
	{
		if (distance[i] < 0)
		{
            std::cout << "Error in DistancInterpolation(): ";
            std::cout << "minus distance found!" << std::endl;
			exit(1);
		}
		if (distance[i] < smallest)
		{
			smallest = distance[i];
			n_smallest = i;
		}
	}

	Tensor numerator = value[n_smallest];
	Scalar dominator = 1.0;
	for (int i = 0; i < (int)distance.size(); i++)
	{
		if (i == n_smallest) continue;
		Scalar added;
		if (distance[i] < SMALL)
		{
			added = 1;
		}
		else
		{
			added = smallest / distance[i];
		}
		dominator += added;
		numerator = numerator + added*value[i];
	}
	return numerator / dominator;
}

//Gaussian elimination
template<>
std::vector<Scalar> Gaussin_L(std::vector<std::vector<Scalar> >& matrix, std::vector<Scalar>& b)
{
	int n = (int)b.size();
	//Check sizes of the matrix
	if (n != (int)matrix.size())
	{
		FatalError("can't use Gaussian elimination: error found in the size of the matrix");
	}
	for (int k = 0; k < n; k++)
	{
		if (n != (int)matrix[k].size())
		{
			FatalError("can't use Gaussian elimination: error found in the size of the matrix");
		}
	}
	for (int k = 0; k < n - 1; k++)
	{
		//Find Max value position
		int row = k;
		Scalar ave = 0;
		for (int i = k; i < n; i++)
		{
			if (fabs(matrix[i][k]) > ave)
			{
				ave = fabs(matrix[i][k]);
				row = i;
			}
		}
		if (fabs(matrix[row][k]) < SMALL)
		{
			FatalError("can't use Gaussian elimination: one row in the matrix is found all zero");
		}
		if (k != row)
		{
			for (int i = 0; i < n; i++)
			{
				Scalar temp = matrix[row][i];
				matrix[row][i] = matrix[k][i];
				matrix[k][i] = temp;
			}
			Scalar temp = b[k];
			b[k] = b[row];
			b[row] = temp;
		}
		//Eliminate process
        std::vector<Scalar> c;
		c.resize(n);
		for (int j = k + 1; j < n; j++)
		{
			c[j] = matrix[j][k] / matrix[k][k];
		}
		for (int i = k + 1; i < n; i++)
		{
			for (int j = 1; j < n; j++)
			{
				matrix[i][j] = matrix[i][j] - c[i] * matrix[k][j];
			}
			b[i] = b[i] - c[i] * b[k];
		}
	}
	//Back substitution process
    std::vector<Scalar> x;
	x.resize(n);
	x[n - 1] = b[n - 1] / matrix[n - 1][n - 1];
	for (int i = n - 2; i >= 0; i--)
	{
		Scalar sum = 0;
		for (int j = i + 1; j < n; j++)
		{
			sum += matrix[i][j] * x[j];
		}
		x[i] = (b[i] - sum) / matrix[i][i];
	}
	return x;
}

template<>
std::vector<Vector> Gaussin_L(std::vector<std::vector<Scalar> >& matrix, std::vector<Vector>& b)
{
	int n = (int)b.size();
	//Check sizes of the matrix
	if (n != (int)matrix.size())
	{
		FatalError("can't use Gaussian elimination: error found in the size of the matrix");
	}
	for (int k = 0; k < n; k++)
	{
		if (n != (int)matrix[k].size())
		{
			FatalError("can't use Gaussian elimination: error found in the size of the matrix");
		}
	}
	for (int k = 0; k < n - 1; k++)
	{
		//Find Max value pos
		int row = k;
		Scalar ave = 0;
		for (int i = k; i < n; i++)
		{
			if (fabs(matrix[i][k]) > ave)
			{
				ave = fabs(matrix[i][k]);
				row = i;
			}
		}
		if (fabs(matrix[row][k]) < SMALL)
		{
			FatalError("can't use Gaussian elimination: one row in the matrix is found all zero");
		}
		if (k != row)
		{
			for (int i = 0; i < n; i++)
			{
				Scalar temp = matrix[row][i];
				matrix[row][i] = matrix[k][i];
				matrix[k][i] = temp;
			}
			Vector temp = b[k];
			b[k] = b[row];
			b[row] = temp;
		}
		
		//Eliminate process
        std::vector<Scalar> c;
		c.resize(n);
		for (int j = k + 1; j < n; j++)
		{
			c[j] = matrix[j][k] / matrix[k][k];
		}
		for (int i = k + 1; i < n; i++)
		{
			for (int j = 1; j < n; j++)
			{
				matrix[i][j] = matrix[i][j] - c[i] * matrix[k][j];
			}
			b[i] = b[i] - c[i] * b[k];
		}
	}
	
	//Back substitution process
    std::vector<Vector> x;
	x.resize(n);
	x[n - 1] = b[n - 1] / matrix[n - 1][n - 1];
	for (int i = n - 2; i >= 0; i--)
	{
		Vector sum = Vector(0.0, 0.0, 0.0);
		for (int j = i + 1; j < n; j++)
		{
			sum += matrix[i][j] * x[j];
		}
		x[i] = (b[i] - sum) / matrix[i][i];
	}
	return x;
}

Scalar BasicFunction(const Scalar x, const Scalar sigma)
{
	Scalar temp = x / sigma;
	return exp(-temp*temp);
	//return pow((x / sigma), 1.5);
}

//Radio Basis Function Interpolation
template<>
Scalar RBFInterpolation
(
//list of positions interpolated from
const std::vector<Vector>& rePosition,
//list of corresponding values
const std::vector<Scalar>& value
)
{
	//Check existence and and the number of interpolating points
	if (rePosition.size() != value.size())
	{
		FatalError("in RBFInterpolation(): Numbers of distances and values are inconsistent");
	}
	else if (0 == rePosition.size())
	{
		FatalError("in RBFInterpolation(): no data is provided for interpolation!");
	}
	int size = (int)rePosition.size();
	//calculate average distance and average value
	Scalar aveDis = 0.0;
	Scalar aveValue = 0.0;
	for (int i = 0; i < size; i++)
	{
		aveDis += rePosition[i].Mag();
		aveValue += value[i];
	}
	aveDis = aveDis / (Scalar)size;
	Scalar sigma = 1.5*aveDis;
	aveValue = aveValue / (Scalar)size;
	//calculate relative values
    std::vector<Scalar> reValue;
	reValue.resize(size);
	for (int i = 0; i < size; i++)
	{
		reValue[i] = value[i] - aveValue;
	}
	//Generate Linear Equation
    std::vector<std::vector<Scalar> > matrix;
	AllocateDim(matrix, size, size);
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			matrix[i][j] = matrix[j][i] = BasicFunction((rePosition[i] - rePosition[j]).Mag(), sigma);
		}
	}
	//solving linear equation
    std::vector<Scalar> omega = Gaussin_L(matrix, reValue);
	Scalar fout = 0.0;
	for (int k = 0; k < size; k++)
	{
		fout += omega[k] * BasicFunction(rePosition[k].Mag(), sigma);
	}
	fout += aveValue;

	return fout;
}

//Calculate gradient a set of surrounding points using least square fit
Vector LeastSquareGradient2D
(
const std::vector<Vector>& rePosition,
const std::vector<Scalar>& value
)
{
	int numPoint = (int)rePosition.size();
	if (0 == numPoint)
	{
		FatalError("in LeastSquareGradient, no neighbor point is given");
	}
	if ((int)value.size() != numPoint)
	{
		FatalError("in LeastSquareGradient, the numbers of points and values are not consistent");
	}
    std::vector<Scalar> omega;
	omega.resize(numPoint);
	for (int k = 0; k < numPoint; k++)
	{
		omega[k] = 1.0 / rePosition[k].Mag();
	}
	//create and initialize matrix for least square fit
    std::vector<std::vector<Scalar> > matrix;
    std::vector<Scalar> rhs;
	AllocateDim(matrix, 2, 2);
	AllocateDim(rhs, 2);
	for (int i = 0; i < 2; i++)
	{
		rhs[i] = 0.0;
		for (int j = 0; j < 2; j++)
		{
			matrix[i][j] = 0.0;
		}
	}
	//assign values for the matrix
	for (int k = 0; k < numPoint; k++)
	{
		rhs[0] += omega[k] * rePosition[k].x_*value[k];
		rhs[1] += omega[k] * rePosition[k].y_*value[k];
		matrix[0][0] += omega[k] * rePosition[k].x_*rePosition[k].x_;
		matrix[1][1] += omega[k] * rePosition[k].y_*rePosition[k].y_;
		Scalar tempxy = omega[k] * rePosition[k].x_*rePosition[k].y_;
		matrix[0][1] += tempxy;
		matrix[1][0] += tempxy;
	}
	//solving linear equation
    std::vector<Scalar> component = Gaussin_L(matrix, rhs);
	return Vector(component[0], component[1], 0.0);
}

//Calculate gradient a set of surrounding points using least square fit
Vector LeastSquareGradient3D
(
const std::vector<Vector>& rePosition,
const std::vector<Scalar>& value
)
{
	int numPoint = (int)rePosition.size();
	if (0 == numPoint)
	{
		FatalError("in LeastSquareGradient, no neighbor point is given");
	}
	if ((int)value.size() != numPoint)
	{
		FatalError("in LeastSquareGradient, the numbers of points and values are not consistent");
	}
    std::vector<Scalar> omega;
	omega.resize(numPoint);
	for (int k = 0; k < numPoint; k++)
	{
		omega[k] = 1.0 / rePosition[k].Mag();
	}
	//create and initialize matrix for least square fit
    std::vector<std::vector<Scalar> > matrix;
    std::vector<Scalar> rhs;
	AllocateDim(matrix, 3, 3);
	AllocateDim(rhs, 3);
	for (int i = 0; i < 3; i++)
	{
		rhs[i] = 0.0;
		for (int j = 0; j < 3; j++)
		{
			matrix[i][j] = 0.0;
		}
	}
	//assign values for the matrix
	for (int k = 0; k < numPoint; k++)
	{
		rhs[0] += omega[k] * rePosition[k].x_*value[k];
		rhs[1] += omega[k] * rePosition[k].y_*value[k];
		rhs[2] += omega[k] * rePosition[k].z_*value[k];
		matrix[0][0] += omega[k] * rePosition[k].x_*rePosition[k].x_;
		matrix[1][1] += omega[k] * rePosition[k].y_*rePosition[k].y_;
		matrix[2][2] += omega[k] * rePosition[k].z_*rePosition[k].z_;
		Scalar tempxy, tempyz, tempxz;
		tempxy = omega[k] * rePosition[k].x_*rePosition[k].y_;
		tempyz = omega[k] * rePosition[k].y_*rePosition[k].z_;
		tempxz = omega[k] * rePosition[k].x_*rePosition[k].z_;
		matrix[0][1] += tempxy;
		matrix[1][0] += tempxy;
		matrix[0][2] += tempxz;
		matrix[2][0] += tempxz;
		matrix[1][2] += tempyz;
		matrix[2][1] += tempyz;
	}
	//solving linear equation
    std::vector<Scalar> component = Gaussin_L(matrix, rhs);
	return Vector(component[0], component[1], component[2]);
}

//Calculate gradient a set of surrounding points using least square fit in two dimensions
Tensor LeastSquareGradient2D
(
const std::vector<Vector>& rePosition,
const std::vector<Vector>& value
)
{
	int numPoint = (int)rePosition.size();
	if (0 == numPoint)
	{
		FatalError("in LeastSquareGradient, no neighbor point is given");
	}
	if ((int)value.size() != numPoint)
	{
		FatalError("in LeastSquareGradient, the numbers of points and values are not consistent");
	}
    std::vector<Scalar> omega;
	omega.resize(numPoint);
	for (int k = 0; k < numPoint; k++)
	{
		omega[k] = 1.0 / rePosition[k].Mag();
	}
	//create and initialize matrix for least square fit
    std::vector<std::vector<Scalar> > matrix;
    std::vector<Vector> rhs;
	AllocateDim(matrix, 2, 2);
	AllocateDim(rhs, 2);
	for (int i = 0; i < 2; i++)
	{
		rhs[i] = Vector(0.0, 0.0, 0.0);
		for (int j = 0; j < 2; j++)
		{
			matrix[i][j] = 0.0;
		}
	}
	//assign values for the matrix
	for (int k = 0; k < numPoint; k++)
	{
		matrix[0][0] += omega[k] * rePosition[k].x_*rePosition[k].x_;
		matrix[1][1] += omega[k] * rePosition[k].y_*rePosition[k].y_;
		Scalar tempxy = omega[k] * rePosition[k].x_*rePosition[k].y_;
		matrix[0][1] += tempxy;
		matrix[1][0] += tempxy;
		rhs[0] += omega[k] * rePosition[k].x_*value[k];
		rhs[1] += omega[k] * rePosition[k].y_*value[k];
	}
	//solving linear equation
    std::vector<Vector> component = Gaussin_L(matrix, rhs);
	return Tensor
		(component[0].x_, component[0].y_, 0.0,
		component[1].x_, component[1].y_, 0.0,
		0.0, 0.0, 0.0);
}

//Calculate gradient a set of surrounding points using least square fit in three dimensions
Tensor LeastSquareGradient3D
(
const std::vector<Vector>& rePosition,
const std::vector<Vector>& value
)
{
	int numPoint = (int)rePosition.size();
	if (0 == numPoint)
	{
		FatalError("in LeastSquareGradient, no neighbor point is given");
	}
	if ((int)value.size() != numPoint)
	{
		FatalError("in LeastSquareGradient, the numbers of points and values are not consistent");
	}
    std::vector<Scalar> omega;
	omega.resize(numPoint);
	for (int k = 0; k < numPoint; k++)
	{
		omega[k] = 1.0 / rePosition[k].Mag();
	}
	//create and initialize matrix for least square fit
    std::vector<std::vector<Scalar> > matrix;
    std::vector<Vector> rhs;
	AllocateDim(matrix, 3, 3);
	AllocateDim(rhs, 3);
	for (int i = 0; i < 3; i++)
	{
		rhs[i] = Vector(0.0, 0.0, 0.0);
		for (int j = 0; j < 3; j++)
		{
			matrix[i][j] = 0.0;
		}
	}

	for (int k = 0; k < numPoint; k++)
	{
		//assign values for the matrix
		matrix[0][0] += omega[k] * rePosition[k].x_*rePosition[k].x_;
		matrix[1][1] += omega[k] * rePosition[k].y_*rePosition[k].y_;
		matrix[2][2] += omega[k] * rePosition[k].z_*rePosition[k].z_;
		Scalar tempxy, tempyz, tempxz;
		tempxy = omega[k] * rePosition[k].x_*rePosition[k].y_;
		tempyz = omega[k] * rePosition[k].y_*rePosition[k].z_;
		tempxz = omega[k] * rePosition[k].x_*rePosition[k].z_;
		matrix[0][1] += tempxy;
		matrix[1][0] += tempxy;
		matrix[0][2] += tempxz;
		matrix[2][0] += tempxz;
		matrix[1][2] += tempyz;
		matrix[2][1] += tempyz;
		//assign values for rhs of the matrix
		rhs[0] += omega[k] * rePosition[k].x_*value[k];
		rhs[1] += omega[k] * rePosition[k].y_*value[k];
		rhs[2] += omega[k] * rePosition[k].z_*value[k];
	}

	//solving linear equation
    std::vector<Vector> component = Gaussin_L(matrix, rhs);

	return Tensor
		(component[0].x_, component[0].y_, component[0].z_,
		component[1].x_, component[1].y_, component[1].z_,
		component[2].x_, component[2].y_, component[2].z_);
}

//Standard interpolation
//interpolating from a list of points around as well as corresponding values. 
template<>
Scalar Interpolation
(
    //list of relative positions(the coordinate of the point interpolated for is (0,0,0))
    const std::vector<Vector>& rePosition,
    //list of corresponding values
    const std::vector<Scalar>& value
)
{
	//Check existence and and the number of interpolating points
	if (rePosition.size() != value.size())
	{
		FatalError("in Interpolation(): Numbers of distances and values are inconsistent");
	}
	else if (0 == rePosition.size())
	{
		FatalError("in Interpolation(): no data is provided for interpolation!");
	}

	//temporally use DistanceInterpolation, and of course needs to be extended; 
    std::vector<Scalar> dis;
	for (int i = 0; i < (int)rePosition.size(); i++)
	{
		dis.push_back(rePosition[i].Mag());
	}
	return DistanceInterpolation(dis, value);
}

template<>
Vector Interpolation
(
    //list of relative positions(the coordinate of the point interpolated for is (0,0,0))
    const std::vector<Vector>& rePosition,
    //list of corresponding values
    const std::vector<Vector>& value
)
{
	//Check existence and and the number of interpolating points
	if (rePosition.size() != value.size())
	{
		FatalError("in Interpolation(): Numbers of distances and values are inconsistent");
	}
	else if (0 == rePosition.size())
	{
		FatalError("in Interpolation(): no data is provided for interpolation!");
	}

	//temporally use DistanceInterpolation, and of course needs to be extended; 
    std::vector<Scalar> dis;
	for (int i = 0; i < (int)rePosition.size(); i++)
	{
		dis.push_back(rePosition[i].Mag());
	}
	return DistanceInterpolation(dis, value);
}

template<>
Tensor Interpolation
(
    //list of relative positions(the coordinate of the point interpolated for is (0,0,0))
    const std::vector<Vector>& rePosition,
    //list of corresponding values
    const std::vector<Tensor>& value
)
{
	//Check existence and and the number of interpolating points
	if (rePosition.size() != value.size())
	{
		FatalError("in Interpolation(): Numbers of distances and values are inconsistent");
	}
	else if (0 == rePosition.size())
	{
		FatalError("in Interpolation(): no data is provided for interpolation!");
	}

	//temporally use DistanceInterpolation, and of course needs to be extended; 
    std::vector<Scalar> dis;
	for (int i = 0; i < (int)rePosition.size(); i++)
	{
		dis.push_back(rePosition[i].Mag());
	}
	return DistanceInterpolation(dis, value);
}

//========================Linear Interpolation=======================
template<>
Scalar LinearInterpolation
(
    //distances to the point to be interpolated
    const std::pair<Scalar, Scalar>& distances,
    //values interpolated from
    const std::pair<Scalar, Scalar>& values
)
{
	Scalar numerator = distances.first * values.second + distances.second*values.first;
	Scalar denominator = distances.first + distances.second;
	return numerator / denominator;
}

template<>
Vector LinearInterpolation
(
    //distances to the point to be interpolated
    const std::pair<Scalar, Scalar>& distances,
    //values interpolated from
    const std::pair<Vector, Vector>& values
)
{
	Vector numerator = distances.first * values.second + distances.second*values.first;
	Scalar denominator = distances.first + distances.second;
	return numerator / denominator;
}

template<>
Tensor LinearInterpolation
(
    //distances to the point to be interpolated
    const std::pair<Scalar, Scalar>& distances,
    //values interpolated from
    const std::pair<Tensor, Tensor>& values
)
{
	Tensor numerator = distances.first * values.second + distances.second*values.first;
	Scalar denominator = distances.first + distances.second;
	return numerator / denominator;
}

//======================Hamonic Interpolation===========================
template<>
Scalar HamonicInterpolation
(
    //distances to the point to be interpolated
    const std::pair<Scalar, Scalar>& distances,
    //values interpolated from
    const std::pair<Scalar, Scalar>& values
)
{
	Scalar numerator = (distances.first + distances.second)*values.first*values.second;
	Scalar denominator = distances.first * values.second + distances.second*values.first;
	return numerator / denominator;
}



