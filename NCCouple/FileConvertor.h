
#ifndef FILECONVERTOR_HEADER
#define FILECONVERTOR_HEADER
#include "./Solver.h"
#include <string>
enum class Material
{
	H2O,
	Zr4,
	UO2
};
//detailed requirements are needed


void MOCFieldsToCFD
(
	std::string  strInput_aplFileName,
	std::string strInput_inpFileName,
	std::string strInput_txtFileName,
	std::string strInput_meshFileName,
	std::string strOutput_vtkFileName,
	bool bRenew=false//判断是否是第一次做插值，false为第一次
);

//vtk中提温度（K），密度（kg/m3）
//inp中写取温度（K），核子密度（1/cm3）,流体区域密度会变，固体区域密度为常数（直接从原始的inp文件中抄）
void CFDFieldsToMOC
(
	std::string strInput_meshFileName,
	std::string strInput_vtkFileName,
	std::string  strInput_aplFileName,
	std::string strOutput_inpFileName,
	bool bRenew = false//判断是否是第一次做插值，false为第一次
);
void ConservationValidation
(
	const GeneralMesh& sourceMesh,
	const GeneralMesh& targetMesh,
	ValueType vt,
	std::string strZoneName
);
#endif