#include "./Solver.h"
#include <string>

//detailed requirements are needed
//by Ling Kong
void CreateSolver
(
	std::string aplFileName,
	std::string mshFileName,
	std::string weightFileName
)
{
	//read msh file and create CFDMesh list
	//read apl file and create MOCMesh
	//create MOCIndex
	//create solver list (search by region names and create solver for each available pair)
	//save weight factors in files (one file per region)
	return;
}

//txt存储功率(W/cm3)
//vtk存储功率(W/m3)
void MOCFieldsToCFD
(
	std::string weightFileName,
	std::string txtFileName,
	std::string vtkFileName
)
{
	//call ReSetUp()
	//read txt file to get heat power in MOCMesh (W/cm3)
	//interpolation of heat power to CFDMesh list
	//transfer heat power to (W/M3)
	//create fields (one field per region)
	//write vtk file for heat power
	return;
}

//vtk中提温度（K），密度（kg/m3）
//inp中写取温度（K），核子密度（1/cm3）,流体区域密度会变，固体区域密度为常数（直接从原始的inp文件中抄）
void CFDFieldsToMOC
(
	std::string weightFileName,
	std::string vtkFileName,
	std::string inpFileNameIn,
	std::string inpFileNameOut
)
{
	//call ReSetUp()
	//read vtk file to get temperature and density in CFDMeshs
	//interpolation of temperature and density to MOCFile
	//read inp file (in) to get density (1/cm3) in solid regions (if any)
	//transfer density from kg/m3 power to 1/cm3
	//write inp file (out) for temperature (K) and density (1/cm3)
	return;
}

//by Yuan Hui
void ReSetUp()
{
	//read msh file and create CFDMeshes
	//read apl file and create MOCMesh
	//read weight-factor in files
	return;
}

//By Li zi jian
void WriteMultipleRegionVtkFile
(
	const std::vector<FieldIO >& fieldIOList,
	std::string vtkFileName
)
{
	return;
}