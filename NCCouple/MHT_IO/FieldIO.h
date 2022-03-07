
#ifndef _FieldIO_
#define _FieldIO_

#include <vector>
#include <string>
#include "../MHT_common/Configuration.h"
#include "../MHT_common/Vector.h"

class Mesh;

template<class Type>
class Field;

class FieldIO
{
	//pointer to dependent mesh 
	Mesh* p_blockMesh;

	// Element fields pointer for those which should be updated.
	std::vector<Field<Scalar>* > v_scalarField;
	std::vector<Field<Vector>* > v_vectorField;

public:
	//constructor
	FieldIO();

	// Push_back
	void push_backScalarField(Field<Scalar>& ScalarField);
	void push_backVectorField(Field<Vector>& VectorField);

    void WriteTecplotField(const std::string& outMshFileName);
	void WriteTecplotField(const std::string& outMshFileName, Scalar Time);
	void WriteVTKField(const std::string& outFileName);

	void FieldElementToNode();
};

#endif
