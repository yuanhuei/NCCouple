
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

void DisplayHelpInfo();

void ConservationValidation
(
	const GeneralMesh& sourceMesh,
	const GeneralMesh& targetMesh,
	ValueType vt,
	std::string strZoneName,
	Solver& solverMapper
);

void RegisterMapper
(
	std::string strInput_aplFileName,
	std::string strInput_inpFileName,
	std::string strInput_meshFileName
);

void CreateMapper();

void MOCFieldsToCFD();

void CFDFieldsToMOC();

void ClearMapFiles();

#endif