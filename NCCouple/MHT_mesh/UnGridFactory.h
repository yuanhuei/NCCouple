
#ifndef _UnGridFactory_
#define _UnGridFactory_

#include <memory>
#include <vector>
#include "../MHT_common/Configuration.h"
#include "../MHT_mesh/StandardMeshBlock.h"
#include "../MHT_mesh/FluentMeshBlock.h"

class Mesh;

class UnGridFactory
{
public:

	// Enum for UnGrid type
    enum UnGridType
	{
		ugtFluent,
		ugtUGrid,
		ugtStdGrid,
		ugtNoType
    };

	// Variable indicated current grid type
	UnGridType ugt_gridType;

	// Auto pointer for base grid block
	std::shared_ptr<Mesh> p_mesh;

	// Construction function
    UnGridFactory(const std::string& MshFileName, UnGridFactory::UnGridType GridType = UnGridFactory::ugtFluent);

	// Construction function
	// Construct and get mesh pointer
    UnGridFactory(const std::string& MshFileName, std::vector<Mesh*>& v_pmesh, UnGridFactory::UnGridType GridType = UnGridFactory::ugtFluent);

	Mesh* GetPtr();
};

#endif
