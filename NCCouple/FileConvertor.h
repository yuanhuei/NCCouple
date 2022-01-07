
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
	std::string strZoneName
);

void CreateMapper
(
	std::string strInput_aplFileName,
	std::string strInput_inpFileName,
	std::string strInput_meshFileName
);

void MOCFieldsToCFD();

//vtk�����¶ȣ�K�����ܶȣ�kg/m3��
//inp��дȡ�¶ȣ�K���������ܶȣ�1/cm3��,���������ܶȻ�䣬���������ܶ�Ϊ������ֱ�Ӵ�ԭʼ��inp�ļ��г���
void CFDFieldsToMOC();

#endif