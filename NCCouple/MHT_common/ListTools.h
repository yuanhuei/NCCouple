#ifndef _ListTools_
#define _ListTools_

#include <vector>
#include <map>

#include "../MHT_common/Configuration.h"
#include "../MHT_common/Vector.h"

std::vector<int> GetNonRepeatedList(const std::vector<int>& originalList);

//Operators for vector<Type>
template<class Type>
std::vector<Type> operator + (const std::vector<Type>&, const std::vector<Type>&);

template<class Type>
std::vector<Type> operator - (const std::vector<Type>&, const std::vector<Type>&);

//Kong Ling supplemented on 2019/12/23
template<class Type>
std::vector<Type> operator *(const Scalar&, const std::vector<Type>&);

template<class Type>
std::vector<Type> operator * (const std::vector<Scalar>&, const std::vector<Type>&);

template<class Type>
std::vector<Type> operator / (const std::vector<Type>&, const std::vector<Scalar>&);

//Kong Ling supplemented on 2019/12/22
std::vector<Scalar> operator & (const std::vector<Vector>&, const std::vector<Vector>&);

//Operator +
template<class T1, class T2>
std::map<T1, T2> operator + (const std::map<T1, T2>& left, const std::map<T1, T2>& right);

//Operator += 
template<class T1, class T2>
void operator += (std::map<T1, T2>& left, const std::map<T1, T2>& right);

//Operator -
template<class T1, class T2>
std::map<T1, T2> operator - (const std::map<T1, T2>& left, const std::map<T1, T2>& right);

//Operator -= 
template<class T1, class T2>
void operator -= (std::map<T1, T2>& left, const std::map<T1, T2>& right);

//Operator *
template<class T1, class T2>
std::map<T1, T2> operator * (const std::map<T1, T2>& left, const std::map<T1, T2>& right);

//Operator *=
template<class T1, class T2>
void operator *= (std::map<T1, T2>& left, const std::map<T1, T2>& right);

template<class Type>
void DisplayMaximumNormOf(std::vector<Type>&);

#endif
