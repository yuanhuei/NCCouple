
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
	bool bRenew=false//�ж��Ƿ��ǵ�һ������ֵ��falseΪ��һ��
);

void CFDFieldsToMOC
(
	std::string strInput_meshFileName,
	std::string strInput_vtkFileName,
	std::string  strInput_aplFileName,
	std::string strOutput_inpFileName,
	bool bRenew = false//�ж��Ƿ��ǵ�һ������ֵ��falseΪ��һ��
);
void ConservationValidation
(
	const GeneralMesh& sourceMesh,
	const GeneralMesh& targetMesh,
	ValueType vt
);
#endif