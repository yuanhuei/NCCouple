/*---------------------------------------------------------------------------*\
File Name:
	MaxMin.h
Description:
	A basic function of Max,Min,Sign.

	Author:		ShuaiZhang
	Date: Long ago
\*---------------------------------------------------------------------------*/

#pragma   once

#ifndef _MaxMin_
#define _MaxMin_

template<class Type>
inline Type Min(Type FirstVal, Type SecondVal)
{
	return FirstVal < SecondVal ? FirstVal : SecondVal;
}

//Return bigger value
template<class Type>
inline Type Max(Type FirstVal, Type SecondVal)
{
	return FirstVal > SecondVal ? FirstVal : SecondVal;
}

//Return sign of value
template<class Type>
inline Type Sign(Type value)
{
	return value > 0 ? 1 : (-1);
}

#endif