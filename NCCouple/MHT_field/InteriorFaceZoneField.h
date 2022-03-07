
#ifndef _InteriorFaceZoneField_
#define _InteriorFaceZoneField_

#include <string>

#include "../MHT_common/Configuration.h"

class Mesh;

template<class Type>
class ElementField;

template<class Type>
class FaceField;

//a class handling interior face zones, playing as a part of faceField 
template<class Type>
class InteriorFaceZoneField
{
public:
	//pointer to dependent mesh 
	Mesh* p_blockMesh;
	//constructor
	InteriorFaceZoneField(Mesh* UGB);

	//interpolating face values from adjacent elements
    void CalculateValue(ElementField<Type>& EF, FaceField<Type>& FF, const std::string& approach);
};

#endif
