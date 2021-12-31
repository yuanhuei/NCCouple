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
	Solver& solver
	std::string inpFileName,
	std::string vtkFileName
	);

void CFDFieldsToMOC
(
	Solver& solver
	std::string inpFileName,
	std::string vtkFileName
);