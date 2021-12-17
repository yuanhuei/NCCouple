
#ifndef _BoundaryCalculationUtility_
#define _BoundaryCalculationUtility_

#include <string>
#include "../MHT_common/Configuration.h"
#include "../MHT_field/ElementField.h"
#include "../MHT_field/FaceField.h"

int ElementToFaceParameterCheck(std::string input);

template<class Type>
void ElementToFaceCheck(ElementField<Type>& EF, FaceField<Type>& FF);

#endif
