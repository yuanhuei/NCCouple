#include "VTKIO.h"

MHTVTKReader::MHTVTKReader()
{
}


MHTVTKReader::MHTVTKReader(std::string MshFileName, std::string VTKFileName,  std::vector<std::string>& vFiedNameList)
{


	ReadMSHFile(MshFileName);

//	ReadVTKFile(VTKFileName, vFiedNameList);
	
	ReadDataFile(VTKFileName, vFiedNameList);
}

MHTVTKReader::~MHTVTKReader()
{
	std::vector<Field<Scalar>> ().swap(v_scalarFieldList);

	std::vector<FieldIO>().swap(v_FieldIO);
	std::vector<StandardMeshBlock>().swap(v_stdMesh);
}

void MHTVTKReader::ReadVTKFile(std::string vVTKFileName, std::vector<std::string>& vFiedNameList)
{
	std::cout << "start Read Field" << std::endl;
/*
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
	}*/

}

void MHTVTKReader::ReadMSHFile(std::string MeshFileName)
{
	std::cout<<"start Read Mesh" << std::endl;
	UnGridFactory meshFactoryCon(MeshFileName, UnGridFactory::ugtFluent);
	FluentMeshBlock* FluentPtrCon = dynamic_cast<FluentMeshBlock*>(meshFactoryCon.GetPtr());
	RegionConnection Bridges;
	FluentPtrCon->Decompose(Bridges);
	v_stdMesh.resize(FluentPtrCon->v_regionGrid.size());
	for (size_t i = 0; i < FluentPtrCon->v_regionGrid.size(); i++)
	{
		v_stdMesh[i] = std::move(FluentPtrCon->v_regionGrid[i]);
		v_pmesh.push_back(&v_stdMesh[i]);
	}
}

void MHTVTKReader::ReadDataFile(std::string DataFileName, std::vector<std::string>& vFiedNameList)
{
	std::cout << "start Read Data Field" << std::endl;

	v_FieldIO.resize(v_pmesh.size());
	
	int nTotalFieldNum(v_pmesh.size()* vFiedNameList.size());

	std::ifstream inFile(DataFileName);

	std::vector<std::string> vFieldName();

	std::string sMeshName;
	std::string sFieldName;
	int nFieldNum;
	for (size_t i = 0; i < v_pmesh.size(); i++)
	{
		for (size_t j = 0; j < vFiedNameList.size(); j++)
		{	
			inFile >> sMeshName;
			inFile >> sFieldName;
			inFile >> nFieldNum;

			if (vFiedNameList[j]!= sFieldName)
			{
				FatalError("read field name and given field name are not same please check your input or dataFile");
			}

			Field<Scalar> thisField(v_pmesh[i], 0.0, sFieldName);
			
			for (size_t fieldIndex = 0; fieldIndex < nFieldNum; fieldIndex++)
			{
				double dFieldData;
				inFile >> dFieldData;
				thisField.elementField.SetValue(fieldIndex, dFieldData);
			}
			v_scalarFieldList.push_back(std::move(thisField));
		}
	}
	//do this must after v_scalarFieldList.push_back because point will change
	for (size_t i = 0; i < v_pmesh.size(); i++)
	{
		for (size_t j = 0; j < vFiedNameList.size(); j++)
		{
			v_FieldIO[i].push_backScalarField(v_scalarFieldList[i * vFiedNameList.size() + j]);
		}
	}
}

