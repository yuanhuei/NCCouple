/*---------------------------------------------------------------------------*\
File Name:
	AllocateArray.h
Description:
	A basic function of array allocation.

	Author:		ShuaiZhang, KongLing
	Date: Long ago
\*---------------------------------------------------------------------------*/

#pragma once

#ifndef _AllocateArray_
#define _AllocateArray_

#include <vector>

using namespace std;

// Allocate array for 1D,2D and 3D
template<class Type>
void AllocateDim(std::vector<Type>& CrtVector, int M)
{
	CrtVector.resize(M);
}

template<class Type>
void AllocateDim(std::vector<std::vector<Type> >& CrtVector, int M, int N)
{
	CrtVector.resize(M);
	for (unsigned int i = 0; i < CrtVector.size(); i++)
	{
		CrtVector[i].resize(N);
	}
}

template<class Type>
void AllocateDim(std::vector<std::vector<std::vector<Type> > >& CrtVector, int M, int N, int P)
{
	CrtVector.resize(M);
	for (unsigned int i = 0; i < CrtVector.size(); i++)
	{
		CrtVector[i].resize(N);
		for (unsigned int j = 0; j < CrtVector[i].size(); j++)
		{
			CrtVector[i][j].resize(P);
		}
	}
}

#endif
