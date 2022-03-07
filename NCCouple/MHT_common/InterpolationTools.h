
#ifndef _InterpolationTools_
#define _InterpolationTools_

#include <vector>

#include "../MHT_common/Configuration.h"
#include "../MHT_common/Vector.h"
#include "../MHT_common/Tensor.h"

template<class Type>
Type DistanceInterpolation
(
    const std::vector<Scalar>& distance,
    const std::vector<Type>& value
);

//Radio Basis Function Interpolation
template<class Type>
Type RBFInterpolation
(
    //list of positions interpolated from
    const std::vector<Vector>&,
    //list of corresponding values
    const std::vector<Type>&
);

//Gaussian elimination
template<class Type>
std::vector<Type> Gaussin_L
(
    std::vector<std::vector<Scalar> >&,
    std::vector<Type>&
);

//Calculate gradient a set of surrounding points using least square fit in two dimensions
Vector LeastSquareGradient2D
(
    const std::vector<Vector>& rePosition,
    const std::vector<Scalar>& value
);

//Calculate gradient a set of surrounding points using least square fit in three dimensions
Vector LeastSquareGradient3D
(
    const std::vector<Vector>& rePosition,
    const std::vector<Scalar>& value
);

//Calculate gradient a set of surrounding points using least square fit in two dimensions
Tensor LeastSquareGradient2D
(
    const std::vector<Vector> &rePosition,
    const std::vector<Vector> &value
);

//Calculate gradient a set of surrounding points using least square fit in three dimensions
Tensor LeastSquareGradient3D
(
    const std::vector<Vector>& rePosition,
    const std::vector<Vector>& value
);

//Standard interpolation
//interpolating from a list of points around as well as corresponding values. 
template<class Type>
Type Interpolation
(
    //list of relative positions(the coordinate of the point interpolated for is (0,0,0))
    const std::vector<Vector>& rePosition,
    //list of corresponding values
    const std::vector<Type>& value
);

template<class Type>
Type LinearInterpolation
(
    //distances to the point to be interpolated
    const std::pair<Scalar, Scalar>& distances,
    //values interpolated from
    const std::pair<Type, Type>& values
);

template<class Type>
Type HamonicInterpolation
(
    //distances to the point to be interpolated
    const std::pair<Scalar, Scalar>& distances,
    //values interpolated from
    const std::pair<Type, Type>& values
);

#endif
