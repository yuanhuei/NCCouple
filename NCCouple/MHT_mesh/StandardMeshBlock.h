
#ifndef _StandardMeshBlock_
#define _StandardMeshBlock_

#include <string>
#include "../MHT_common/Configuration.h"
#include "../MHT_mesh/Mesh.h"

class StandardMeshBlock :public Mesh
{
public:

	int nFaceNodeNum;

    StandardMeshBlock(const std::string& MshFileName);

    void ReadMesh(std::ifstream& inFile);
    void WriteTecplotMesh(const std::string& outMshFileName);

    void WriteTecplotMeshOneElem(const std::string& outMshFileName, int elemID);

	StandardMeshBlock()
	{}
};

#endif
