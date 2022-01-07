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
	MHTVTKReader(std::string MshFileName, std::vector<std::string>vVTKFileName, std::vector<std::string>& vFiedNameList);
	MHTVTKReader(std::string MshFileName, std::vector<std::string>& vFiedNameList);
	MHTVTKReader(std::string MshFileName);

	~MHTVTKReader();
	std::vector<FieldIO> GetFieldIOList() { return v_FieldIO; }
	FieldIO GetFieldIO(int num) { return v_FieldIO[num]; }
	std::vector<StandardMeshBlock> GetMeshList() { return v_stdMesh; }
	std::vector<Mesh*> GetMeshListPtr() { return v_pmesh; }
	const std::vector < std::vector<Field<Scalar>>> GetFieldList() { return vv_scalarFieldList; }
	Field<Scalar> GetField(int regionNum, int fieldNum) { return vv_scalarFieldList[regionNum][fieldNum]; }
	std::vector<int> GetMeshID() { return v_meshID; }

	std::vector<int> GetMeshID() { return v_meshID; }

	void WriteDataFile(std::string DataFileName);
	void ReadVTKFile(std::vector<std::string>, std::vector<std::string>& vFiedNameList);		//read total region field
	void ReadVTKFile(std::vector<std::string>, std::vector<int> vMeshID, std::vector<std::string>& vFiedNameList);		//read field by mesh ID
private:

	std::vector < std::vector<Field<Scalar> > > vv_scalarFieldList;			//size same like region in mesh
	std::vector<FieldIO> v_FieldIO;
	std::vector<StandardMeshBlock> v_stdMesh;
	std::vector<Mesh*> v_pmesh;
	std::vector<int> v_meshID;

	void ReadMSHFile(std::string MeshFileName);
	void ReadDataFile(std::string DataFileName, std::vector<std::string>& vFiedNameList);
	void InitializeEmptyField(std::vector<std::string>& vFiedNameList);
};



#endif