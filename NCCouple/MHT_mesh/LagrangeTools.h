/*---------------------------------------------------------------------------*\
File Name:
	LagrangeTools.h

Description:
	A set of tools dealing with lagrange points floating in the computational domain

	Author:		Kong Ling
	Date: 2016-12-02
\*---------------------------------------------------------------------------*/

#ifndef _LagrangeTools_
#define _LagrangeTools_

#include "../MHT_mesh/Mesh.h"
#include "../MHT_common/Configuration.h"
#include "../MHT_common/Vector.h"

class Mesh;

//find the nearest element from a given lagrange point;
int NearestElement
(
    //background grid block
    const Mesh* UGB,
    //coordinate of the given lagrange point
    const Vector &lagPoint,
    //ID of element where searching starts from
    int startEID
);

//find the nearest face from a given lagrange point on a specific patch;
int SearchNearestBoundaryFace
(
    //background grid block
    const Mesh *UGB,
    //coordinate of the given lagrange point
    const Vector &lagPoint,
    //ID of the face which the search starts from
    int startID
);

#endif
