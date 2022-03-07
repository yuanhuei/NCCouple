/*---------------------------------------------------------------------------*\
Class:
	UnGridFactory
File Name:
	UnGridFactory.h
Description:
	Class describe UnGridFactory.
	1.Factory of UnGrid;
	2.Generate corresponding object for different mesh type;

	Author:		ShuaiZhang
	Date: 2017-01-12

Revised:
	Description:
	1.Made addation in construction function which allow function get current mesh pointer
	Revisor:	ShuaiZhang
	Modified Date: 2017-01-13
\*---------------------------------------------------------------------------*/

#include "../MHT_mesh/UnGridFactory.h"
#include "../MHT_common/SystemControl.h"

// Construction function
UnGridFactory::UnGridFactory(const std::string& MshFileName, UnGridFactory::UnGridType GridType)
	:
	ugt_gridType(GridType)
{
	switch (this->ugt_gridType)
	{
	case UnGridFactory::ugtStdGrid:

		this->p_mesh.reset(new StandardMeshBlock(MshFileName));
		break;

	case UnGridFactory::ugtFluent:

		this->p_mesh.reset(new FluentMeshBlock(MshFileName));
		break;

	case UnGridFactory::ugtUGrid:

		//this->p_mesh.reset(new UGridMeshBlock(MshFileName));
		break;

	case UnGridFactory::ugtNoType:

        FatalError(std::string("Grid type is not supported. Please check input grid type.\n")
            + std::string("Current grid name:\n") + MshFileName);
		break;
	}
}

// Construction function
// Construct and get mesh pointer
UnGridFactory::UnGridFactory(const std::string& MshFileName, std::vector<Mesh*>& v_pmesh, UnGridFactory::UnGridType GridType)
	:
	ugt_gridType(GridType)
{
	switch (this->ugt_gridType)
	{
	case UnGridFactory::ugtStdGrid:

		this->p_mesh.reset(new StandardMeshBlock(MshFileName));
		break;

	case UnGridFactory::ugtFluent:

		this->p_mesh.reset(new FluentMeshBlock(MshFileName));
		break;

	case UnGridFactory::ugtUGrid:

		//this->p_mesh.reset(new UGridMeshBlock(MshFileName));
		break;

	case UnGridFactory::ugtNoType:

        FatalError(std::string("Grid type is not supported. Please check input grid type.\n")
            + std::string("Current grid name:\n") + MshFileName);
		break;
	}

	// Get pointer
	v_pmesh.push_back(p_mesh.get());
}

Mesh* UnGridFactory::GetPtr()
{
	return this->p_mesh.get();
}


