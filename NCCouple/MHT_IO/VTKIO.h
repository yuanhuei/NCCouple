#ifndef _VTKIO_H_
#define _VTKIO_H_

#include "../MHT_field/Field.h"
#include "../MHT_IO/FieldIO.h"
#include "../MHT_mesh/Mesh.h"
#include "../MHT_mesh/StandardMeshBlock.h"
#include "../MHT_mesh/UnGridFactory.h"
#include "../MHT_mesh/RegionConnection.h"
#include "../MHT_common/SystemControl.h"
#include <memory>
class Mesh;

template<class Type>
class Field;

class MHTVTKReader
{
public:
	MHTVTKReader();
	MHTVTKReader(std::string MshFileName,std::string vVTKFileName, std::vector<std::string>& vFiedNameList);

	MHTVTKReader(std::string MshFileName, std::vector<std::string>& vFiedNameList);

	MHTVTKReader(std::string MshFileName);
	~MHTVTKReader();

	std::vector<FieldIO> GetFieldIOList() { return v_FieldIO; }
	FieldIO GetFieldIO(int num) { return v_FieldIO[num]; }

	std::vector<StandardMeshBlock> GetMeshList() { return v_stdMesh; }
	std::vector<Mesh*> GetMeshListPtr() { return v_pmesh; }

	const std::vector < std::vector<Field<Scalar>>> GetFieldList() { return vv_scalarFieldList; }
	Field<Scalar> GetField(int regionNum, int fieldNum) { return vv_scalarFieldList[regionNum][fieldNum]; }

	void WriteDataFile(std::string DataFileName);
private:

	std::vector < std::vector<Field<Scalar>>> vv_scalarFieldList;
	std::vector<FieldIO> v_FieldIO;
	std::vector<StandardMeshBlock> v_stdMesh;
	std::vector<Mesh*> v_pmesh;
	void ReadVTKFile(std::string VTKFileName, std::vector<std::string>& vFiedNameList);
	void ReadMSHFile(std::string MeshFileName);
	void ReadDataFile(std::string DataFileName, std::vector<std::string>& vFiedNameList);
	void InitializeEmptyField(std::vector<std::string>& vFiedNameList);
};



#endif