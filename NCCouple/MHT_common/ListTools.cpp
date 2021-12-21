/*---------------------------------------------------------------------------*\
File Name:
	ListTools.cpp

Description:
	Operators between map type objects.

	Author:		Kong Ling
	Date: 2018-03-31

Revised:
    Description:
    1. Convert "non const" input parameters to const;
    2. Convert "non const" function to const function;
    3. Overload Point3D<Scalar> operator + (const Point3D<Scalar>& lhs, const Point3D<Scalar>& rhs)
    4. Overload Point3D<Scalar> operator - (const Point3D<Scalar>& lhs, const Point3D<Scalar>& rhs)

    Revisor:		Shuai Zhang
    Modified Date:	2018-11-28

Revised:
    Description:
    1. Convert "std::map<T1, T2> operator + (std::map<T1, T2>& left, std::map<T1, T2>& right);" To
    "std::map<T1, T2> operator + (const std::map<T1, T2>& left, const std::map<T1, T2>& right);"
    Notice:"const_cast" should be used likes "psmall = const_cast<std::map<unsigned int, Scalar> * > (&right);"
    2. Convert "std::map<T1, T2> operator - (std::map<T1, T2>& left, std::map<T1, T2>& right);" To
    "std::map<T1, T2> operator - (const std::map<T1, T2>& left, const std::map<T1, T2>& right);"
    Notice:"const_iterator" should be used likes "std::map<unsigned int, Scalar>::const_iterator it1;"

    Revisor:		Shuai Zhang
    Modified Date:	2018-12-03
\*---------------------------------------------------------------------------*/

#include "../MHT_common/ListTools.h"
#include "../MHT_common/Vector.h"
#include <vector>

template<>
std::vector<Scalar> operator + (const std::vector<Scalar>& left, const std::vector<Scalar>& right)
{
	int num = (int)left.size();
	std::vector<Scalar> result;
	result.resize(num);
	for (int i = 0; i < num; i++)
	{
		result[i] = left[i] + right[i];
	}
	return result;
}

template<>
std::vector<Vector> operator + (const std::vector<Vector>& left, const std::vector<Vector>& right)
{
	int num = (int)left.size();
	std::vector<Vector> result;
	result.resize(num);
	for (int i = 0; i < num; i++)
	{
		result[i] = left[i] + right[i];
	}
	return result;
}

template<>
std::vector<Scalar> operator - (const std::vector<Scalar>& left, const std::vector<Scalar>& right)
{
	int num = (int)left.size();
	std::vector<Scalar> result;
	result.resize(num);
	for (int i = 0; i < num; i++)
	{
		result[i] = left[i] - right[i];
	}
	return result;
}

template<>
std::vector<Vector> operator - (const std::vector<Vector>& left, const std::vector<Vector>& right)
{
	int num = (int)left.size();
	std::vector<Vector> result;
	result.resize(num);
	for (int i = 0; i < num; i++)
	{
		result[i] = left[i] - right[i];
	}
	return result;
}

//Kong Ling supplemented on 2019/12/23
template<>
std::vector<Scalar> operator *(const Scalar& left, const std::vector<Scalar>& right)
{
	int num = (int)right.size();
	std::vector<Scalar> result;
	result.resize(num);
	for (int i = 0; i < num; i++)
	{
		result[i] = left * right[i];
	}
	return result;
}

//Kong Ling supplemented on 2019/12/23
template<>
std::vector<Vector> operator *(const Scalar& left, const std::vector<Vector>& right)
{
	int num = (int)right.size();
	std::vector<Vector> result;
	result.resize(num);
	for (int i = 0; i < num; i++)
	{
		result[i] = left * right[i];
	}
	return result;
}

template<>
std::vector<Scalar> operator * (const std::vector<Scalar>& left, const std::vector<Scalar>& right)
{
	int num = (int)left.size();
	std::vector<Scalar> result;
	result.resize(num);
	for (int i = 0; i < num; i++)
	{
		result[i] = left[i] * right[i];
	}
	return result;
}

template<>
std::vector<Vector> operator * (const std::vector<Scalar>& left, const std::vector<Vector>& right)
{
	int num = (int)left.size();
	std::vector<Vector> result;
	result.resize(num);
	for (int i = 0; i < num; i++)
	{
		result[i] = left[i] * right[i];
	}
	return result;
}

template<>
std::vector<Scalar> operator / (const std::vector<Scalar>& left, const std::vector<Scalar>& right)
{
	int num = (int)left.size();
	std::vector<Scalar> result;
	result.resize(num);
	for (int i = 0; i < num; i++)
	{
		result[i] = left[i] / right[i];
	}
	return result;
}

template<>
std::vector<Vector> operator / (const std::vector<Vector>& left, const std::vector<Scalar>& right)
{
	int num = (int)left.size();
	std::vector<Vector> result;
	result.resize(num);
	for (int i = 0; i < num; i++)
	{
		result[i] = left[i] / right[i];
	}
	return result;
}

//Kong Ling supplemented on 2019/12/23
std::vector<Scalar> operator & (const std::vector<Vector>& left, const std::vector<Vector>& right)
{
	int num = (int)left.size();
	std::vector<Scalar> result;
	result.resize(num);
	for (int i = 0; i < num; i++)
	{
		result[i] = left[i] & right[i];
	}
	return result;
}

std::vector<int> GetNonRepeatedList(const std::vector<int>& originalList)
{
	std::vector<int> v_nonRepeatList;
	for (int i = 0; i < (int)originalList.size(); i++)
	{
		int currentID = originalList[i];
		bool repeated = false;
		for (int j = 0; j < (int)v_nonRepeatList.size(); j++)
		{
			if (currentID == v_nonRepeatList[j])
			{
				repeated = true;
				break;
			}
		}
		if (false == repeated)
		{
			v_nonRepeatList.push_back(currentID);
		}
	}
	return v_nonRepeatList;
}

//Operator +
template<>
std::map<unsigned int, Scalar> operator + (const std::map<unsigned int, Scalar>& left, const std::map<unsigned int, Scalar>& right)
{
	std::map<unsigned int, Scalar> target;
    std::map<unsigned int, Scalar> *psmall;
	if (left.size() > right.size())
	{
		target = left;
        psmall = const_cast<std::map<unsigned int, Scalar> * > (&right);
	}
	else
	{
		target = right;
        psmall = const_cast<std::map<unsigned int, Scalar> * > (&left);
	}
	std::map<unsigned int, Scalar>::iterator it1, it2;
	for (it1 = (*psmall).begin(); it1 != (*psmall).end(); it1++)
	{ 
		it2 = target.find(it1->first);
		if (target.end() == it2)
		{
			target.insert(*it1);
		}
		else
		{
			(*it2).second += (*it1).second;
		}
	}
	return target;
}


//Operator += 
template<>
void operator += (std::map<unsigned int, Scalar>& left, const std::map<unsigned int, Scalar>& right)
{
	//Version 1 (direct approach): adding right to left directly 
    std::map<unsigned int, Scalar>::const_iterator it1;
    std::map<unsigned int, Scalar>::iterator it2;

	for (it1 = right.begin(); it1 != right.end(); it1++)
	{
		it2 = left.find(it1->first);
		if (left.end() == it2)
		{
			left.insert(*it1);
		}
		else
		{
            (*it2).second += (*it1).second;
		}
	}
}

//Operator -
template<>
std::map<unsigned int, Scalar> operator - (const std::map<unsigned int, Scalar>& left, const std::map<unsigned int, Scalar>& right)
{
    std::map<unsigned int, Scalar> target = left;
    std::map<unsigned int, Scalar>::const_iterator it1;
    std::map<unsigned int, Scalar>::iterator it2;
	for (it1 = right.begin(); it1 != right.end(); it1++)
	{
        it2 = target.find(it1->first);
		if (target.end() == it2)
		{
			target.insert(std::pair<unsigned int, Scalar>((*it1).first, -(*it1).second));
		}
		else
		{
			(*it2).second -= (*it1).second;
		}
	}
	return target;
}

//Operator -= 
template<>
void operator -= (std::map<unsigned int, Scalar>& left, const std::map<unsigned int, Scalar>& right)
{
    std::map<unsigned int, Scalar>::const_iterator it1;
    std::map<unsigned int, Scalar>::iterator it2;
	for (it1 = right.begin(); it1 != right.end(); it1++)
	{
		it2 = left.find(it1->first);
		if (left.end() == it2)
		{
			left.insert(std::pair<unsigned int, Scalar>((*it1).first, -(*it1).second));
		}
		else
		{
			(*it2).second -= (*it1).second;
		}
	}
}

//Operator *
template<>
std::map<unsigned int, Scalar> operator * (const std::map<unsigned int, Scalar>& left, const std::map<unsigned int, Scalar>& right)
{
	std::map<unsigned int, Scalar> *psmall, *plarge;
	if (left.size() > right.size())
	{
        psmall = const_cast<std::map<unsigned int, Scalar> * > (&right);
        plarge = const_cast<std::map<unsigned int, Scalar> * > (&left);
	}
	else
	{
        psmall = const_cast<std::map<unsigned int, Scalar> * > (&left);
        plarge = const_cast<std::map<unsigned int, Scalar> * > (&right);
	}
	std::map<unsigned int, Scalar> target;
	std::map<unsigned int, Scalar>::iterator it1, it2;
	for (it1 = (*psmall).begin(); it1 != (*psmall).end(); it1++)
	{
		it2 = (*plarge).find(it1->first);
		if ((*plarge).end() != it2)
		{
			target.insert(std::pair<unsigned int, Scalar>((*it1).first, (*it1).second*(*it2).second));
		}
	}
	return target;
}

//Operator *=
template<>
void operator *= (std::map<unsigned int, Scalar>& left, const std::map<unsigned int, Scalar>& right)
{
	std::map<unsigned int, Scalar> *psmall, *plarge;
	if (left.size() > right.size())
	{
        psmall = const_cast<std::map<unsigned int, Scalar> * > (&right);
        plarge = &left;
	}
	else
	{
        psmall = &left;
        plarge = const_cast<std::map<unsigned int, Scalar> * > (&right);
	}
	std::map<unsigned int, Scalar> target;
	std::map<unsigned int, Scalar>::iterator it1, it2;
	for (it1 = (*psmall).begin(); it1 != (*psmall).end(); it1++)
	{
		it2 = (*plarge).find(it1->first);
		if ((*plarge).end() != it2)
		{
			target.insert(std::pair<unsigned int, Scalar>((*it1).first, (*it1).second*(*it2).second));
		}
	}
	left=target;
}

template<>
void DisplayMaximumNormOf(std::vector<Scalar>& scalarList)
{
	int num = scalarList.size();
	Scalar valueMax = 0.0;
	int idMax = 0;
	for (int i = 0; i < num; i++)
	{
		if (fabs(scalarList[i])>valueMax)
		{
			valueMax = fabs(scalarList[i]);
			idMax = i;
		}
	}
	std::cout << "max value was found to be " << valueMax << " at id " << idMax << std::endl;
	return;
}

template<>
void DisplayMaximumNormOf(std::vector<Vector>& valueList)
{
	int num = valueList.size();
	Scalar valueMax = 0.0;
	int idMax = 0;
	for (int i = 0; i < num; i++)
	{
		if (valueList[i].Mag() > valueMax)
		{
			valueMax = valueList[i].Mag();
			idMax = i;
		}
	}
	std::cout << "max value was found to be " << valueMax << " at id " << idMax << std::endl;
	return;
}


