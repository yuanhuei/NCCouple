#include "VTKIO.h"

MHTVTKReader::MHTVTKReader()
{
}


MHTVTKReader::MHTVTKReader(std::vector<std::string>& vMshFileName, std::vector<std::string>& vVTKFileName,  std::vector<std::string>& vFiedNameList)
{
	v_stdMesh.resize(vMshFileName.size());

	for (size_t i = 0; i < vMshFileName.size(); i++)
	{
		ReadMSHFile(vMshFileName[i], i);
	}
	
	v_FieldIO.resize(vVTKFileName.size());

	ReadVTKFile(vVTKFileName, vFiedNameList);
	
}

MHTVTKReader::~MHTVTKReader()
{
	std::vector<Field<Scalar>> ().swap(v_scalarFieldList);

	std::vector<FieldIO>().swap(v_FieldIO);
	std::vector<StandardMeshBlock>().swap(v_stdMesh);
}

void MHTVTKReader::ReadVTKFile(std::vector<std::string> vVTKFileName, std::vector<std::string>& vFiedNameList)
{
	std::cout << "start Read Field" << std::endl;
	for (size_t i = 0; i < vVTKFileName.size(); i++)
	{
		vtkObject::GlobalWarningDisplayOff();
		vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
		reader->SetFileName(vVTKFileName[i].c_str());
		reader->SetReadAllColorScalars(true);
		reader->SetReadAllFields(true);
		reader->SetReadAllScalars(true);
		reader->Update();

		vtkSmartPointer<vtkUnstructuredGrid> Grid;
		Grid = reader->GetOutput();
	

		for (size_t j = 0; j < vFiedNameList.size(); j++)
		{
			Field<Scalar> thisField(&v_stdMesh[i], 0.0, vFiedNameList[j]);
			v_scalarFieldList.push_back(std::move(thisField));

			v_scalarFieldList[i * vFiedNameList.size() + j].ReadVTKGridField(Grid, vFiedNameList[j]);	
		}

	}
	for (size_t i = 0; i < vVTKFileName.size(); i++)
	{
		for (size_t j = 0; j < vFiedNameList.size(); j++)
		{
			v_FieldIO[i].push_backScalarField(v_scalarFieldList[i * vFiedNameList.size() + j]);
		}
	}

}

void MHTVTKReader::ReadMSHFile(std::string MeshFileName, int Num)
{
	std::cout<<"start Read Mesh" << std::endl;
	UnGridFactory meshFactoryCon(MeshFileName, UnGridFactory::ugtFluent);
	FluentMeshBlock* FluentPtrCon = dynamic_cast<FluentMeshBlock*>(meshFactoryCon.GetPtr());
	RegionConnection Bridges;
	FluentPtrCon->Decompose(Bridges);
	v_stdMesh[Num] = FluentPtrCon->v_regionGrid[0];
}

