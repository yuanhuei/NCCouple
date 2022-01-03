#ifndef _VTKIO_H_
#define _VTKIO_H_

#include "../MHT_field/Field.h"
#include "../MHT_IO/FieldIO.h"
#include "../MHT_mesh/Mesh.h"
#include "../MHT_mesh/StandardMeshBlock.h"
#include "../MHT_mesh/UnGridFactory.h"
#include "../MHT_mesh/RegionConnection.h"
#include "../MHT_common/SystemControl.h"
class Mesh;

template<class Type>
class Field;

class MHTVTKReader
{
public:
	MHTVTKReader();
	MHTVTKReader(std::vector<std::string>& vMshFileName,std::vector<std::string>& vVTKFileName, std::vector<std::string>& vFiedNameList);
	~MHTVTKReader();

	std::vector<FieldIO> GetFieldIOList() { return v_FieldIO; }
	std::vector<StandardMeshBlock> GetMeshList() { return v_stdMesh; }
	std::vector<Field<Scalar>> GetFieldList() { return v_scalarFieldList; }


	FieldIO GetFieldIO(int num) { return v_FieldIO[num]; }
	StandardMeshBlock GetMesh(int num) { return v_stdMesh[num]; }
private:

	std::vector<Field<Scalar>> v_scalarFieldList;
	
	std::vector<FieldIO> v_FieldIO;
	std::vector<StandardMeshBlock> v_stdMesh;

	void ReadVTKFile(std::vector<std::string> vVTKFileName, std::vector<std::string>& vFiedNameList);
	void ReadMSHFile(std::string MeshFileName, int Num);
};



#endif