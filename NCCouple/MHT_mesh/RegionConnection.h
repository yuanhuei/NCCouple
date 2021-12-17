
#ifndef _RegionConnection_
#define _RegionConnection_

#include "../MHT_common/Configuration.h"
#include "../MHT_mesh/Mesh.h"

class CoupledBoundary
{
public:
	Mesh* p_yangMesh;
	Mesh* p_yinMesh;
    std::string st_yangName;
    std::string st_yinName;
    std::vector<std::pair<int, int> > v_coupledFaceID;
	CoupledBoundary
		(
		Mesh* pmesh1,
		Mesh* pmesh2,
        std::string yangName,
        std::string yinName
		)
		:p_yangMesh(pmesh1), p_yinMesh(pmesh2)
	{
		st_yangName = yangName;
		st_yinName = yinName;
	}

	std::pair<int, int> GetLocalBoundaryZoneID()
	{
		int yangBCID = this->p_yangMesh->GetBoundaryFaceZoneID(st_yangName);
		int yinBCID = this->p_yinMesh->GetBoundaryFaceZoneID(st_yinName);
        return std::pair<int, int>(yangBCID, yinBCID);
	}
};

class RegionConnection
{
public:
    std::vector<CoupledBoundary> v_coupledBoundary;
	RegionConnection()
	{}
    int GetCoupledBoundary(const std::string& yangName, const std::string& yinName)
	{
		int found = -1;
		for (int i = 0; i < (int)v_coupledBoundary.size(); i++)
		{
			if (yangName == v_coupledBoundary[i].st_yangName && yinName == v_coupledBoundary[i].st_yinName)
			{
				found = i;
				break;
			}
		}
		return found;
	}
};

#endif
