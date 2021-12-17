/*---------------------------------------------------------------------------*\
File Name:
	LagrangeTools.cpp

Description:
	A set of tools dealing with lagrange points floating in the computational domain

	Author:		Kong Ling
	Date: 2016-12-02

Revised:
    Description:
    1. Convert "non const" input parameters to const;

    Revisor:		Shuai Zhang
    Modified Date:	2018-12-03
\*---------------------------------------------------------------------------*/

#include "../MHT_mesh/LagrangeTools.h"
#include "../MHT_common/SystemControl.h"

//find the nearest element from a given lagrange point;
int NearestElement
(
    //background grid block
    const Mesh* UGB,
    //coordinate of the given lagrange point
    const Vector& lagPoint,
    //ID of element where searching starts from
    int startEID
)
{
	//check startEID
	if (startEID<0 || startEID>(int)UGB->v_elem.size())
	{
		FatalError("given element ID out of range in NearestElement()");
	}
	int currentID = startEID;
	while (true)
	{
		Scalar currentDis = (UGB->v_elem[currentID].center - lagPoint).Mag();
        std::vector<int> nbEID = UGB->SearchElementNeighbor(currentID);
		int chosenID = currentID;
		for (int i = 0; i < (int)nbEID.size(); i++)
		//check neighbors of current element to see wether closer element can be found
		{
			int newID = nbEID[i];
			Scalar dis = (UGB->v_elem[newID].center - lagPoint).Mag();
			if (dis < currentDis)
			{
				currentDis = dis;
				chosenID = newID;
			}
		}
		//if found, replace the current element; otherwise the current element will be
		//considered as the nearest one, thus break loop and return 
		if (chosenID == currentID)
		{
			break;
		}
		else
		{
			currentID = chosenID;
		}
	}
	return currentID;
}

//find the nearest element from a given lagrange point;
int SearchNearestBoundaryFace
(
    //background grid block
    const Mesh* pmesh,
    //coordinate of the given lagrange point
    const Vector &lagPoint,
    //ID of the face which the search starts from
    int startID
)
{
	int currentID = startID;
	while (true)
	{
		Scalar currentDis = (pmesh->v_face[currentID].center - lagPoint).Mag();
        std::vector<int> v_nbID = pmesh->SearchBounaryFaceNeighbor(currentID);
		int chosenID = currentID;
		for (int i = 0; i < (int)v_nbID.size(); i++)
		//check neighbors of current faces to see wether a closer one can be found
		{
			int newID = v_nbID[i];
			Scalar dis = (pmesh->v_face[newID].center - lagPoint).Mag();
			if (dis < currentDis)
			{
				currentDis = dis;
				chosenID = newID;
			}
		}
		//if found, replace the current face; otherwise the current face will be
		//considered as the nearest one, thus break loop and return 
		if (chosenID == currentID)
		{
			break;
		}
		else
		{
			currentID = chosenID;
		}
	}
	return currentID;
}
