#ifndef STRUCTURE_HEADER
#define STRUCTURE_HEADER

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

struct hash_name {
	size_t operator()(const SMocIndex& p) const {
		return hash<int>()(p.iAssemblyIndex) ^ hash<int>()(p.iCellIndex) ^ hash<int>()(p.iMocIndex);
	};
};

#endif