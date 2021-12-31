
#ifndef FILECONVERTOR_HEADER
#define FILECONVERTOR_HEADER
#include "./Solver.h"
#include <string>

//detailed requirements are needed
void CreateSolver
(
	std::string aplFileName,
	std::string auxFileName,
	std::string mshFileName
);

void MOCFieldsToCFD
(
	std::string  strInput_aplFileName,
	std::string strInput_inpFileName,
	std::string strInput_meshFileName,
	std::string strOutput_vtkFileName
);

void CFDFieldsToMOC
(
	std::string strInput_meshFileName,
	std::string strInput_vtkFileName,
	std::string  strInput_aplFileName,
	std::string strOutput_inpFileName

);

#endif