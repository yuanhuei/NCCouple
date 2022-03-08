#ifndef STRUCTURE_HEADER
#define STRUCTURE_HEADER
#include<functional>
#include<vector>
#include<string>
extern int g_iMpiID;
extern int g_iNumProcs;
extern std::string configFile;
class MOCMesh;
class CFDMesh;

struct pair_hash
{
	template<class T1, class T2>
	std::size_t operator() (const std::pair<T1, T2>& p) const
	{
		auto h1 = std::hash<T1>{}(p.first);
		auto h2 = std::hash<T2>{}(p.second);
		return h1 ^ h2;
	}
};
void ConvergeMocMapInfor();
struct SMocIndex {
	int iAssemblyIndex;
	int iCellIndex;
	int iMocIndex;
	SMocIndex()
	{
		iAssemblyIndex = -1;
		iCellIndex = -1;
		iMocIndex = -1;
	};
	SMocIndex(int i, int j, int k) :iAssemblyIndex(i), iCellIndex(j), iMocIndex(k) {};
	bool operator==(const SMocIndex& p) const
	{
		return iAssemblyIndex == p.iAssemblyIndex && iCellIndex == p.iCellIndex && iMocIndex == p.iMocIndex;
	};
};
namespace std
{
	template<>
	struct hash<SMocIndex>
	{
		size_t operator()(const SMocIndex& p) const {
			return hash<int>()(p.iAssemblyIndex) ^ hash<int>()(p.iCellIndex) ^ hash<int>()(p.iMocIndex);
		};

	};
};

void WriteTotecplot( MOCMesh& mocmesh,  CFDMesh& cfdmesh,
	const std::vector<int> &v_iCFD, const std::vector<SMocIndex>& vSMocindex, std::string fileName);

/*
struct hash_name {
	hash_name() {};
	size_t operator()(const SMocIndex& p) const {
		return std::hash<int>()(p.iAssemblyIndex) ^ std::hash<int>()(p.iCellIndex) ^ std::hash<int>()(p.iMocIndex);
	};
};
*/
#endif