
#ifndef _NumerBC_
#define _NumerBC_

#include "../MHT_common/Configuration.h"
#include "../MHT_common/Vector.h"

template<class Type>
class NumerBC
{
public:
	//parameter a in a*phi + b* \partial phi/\partial n = c;
	Scalar a;
	//parameter b in a*phi + b* \partial phi/\partial n = c;
	Scalar b;
	//parameter c in a*phi + b* \partial phi/\partial n = c;
	Type c;
	//parameter c1 in discretisized form phi_C = c1*phi_b + c2;
	Scalar con1;
	//parameter c2 in discretisized form phi_C = c1*phi_b + c2;
	Type con2;
	//constructor
	NumerBC();
	//constructor
	NumerBC(Scalar, Scalar, Type);
};

#endif
