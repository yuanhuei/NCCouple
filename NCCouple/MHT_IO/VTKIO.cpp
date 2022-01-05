#include "VTKIO.h"

MHTVTKReader::MHTVTKReader()
{
}


MHTVTKReader::MHTVTKReader(std::string MshFileName, std::vector<std::string> vVTKFileName,  std::vector<std::string>& vFiedNameList)
{


	ReadMSHFile(MshFileName);

	this->vv_scalarFieldList.resize(this->v_pmesh.size());

	ReadVTKFile(vVTKFileName, vFiedNameList);
}

MHTVTKReader::MHTVTKReader(std::string MshFileName, std::vector<std::string>& vFiedNameList)
{
	ReadMSHFile(MshFileName);

	this->vv_scalarFieldList.resize(this->v_pmesh.size());

	InitializeEmptyField(vFiedNameList);
}

MHTVTKReader::MHTVTKReader(std::string MshFileName)
{
	ReadMSHFile(MshFileName);

	this->vv_scalarFieldList.resize(this->v_pmesh.size());
}

MHTVTKReader::~MHTVTKReader()
{
	std::vector < std::vector<Field<Scalar>>>().swap(vv_scalarFieldList);
	std::vector<FieldIO>().swap(v_FieldIO);
	std::vector<StandardMeshBlock>().swap(v_stdMesh);
}

void MHTVTKReader::WriteDataFile(std::string DataFileName)
{

	std::ofstream outFile(DataFileName);

	for (size_t i = 0; i < this->vv_scalarFieldList.size(); i++)
	{
		for (size_t j = 0; j < this->vv_scalarFieldList[i].size(); j++)
		{
			outFile << vv_scalarFieldList[i][j].p_blockMesh->st_meshName<<"\t";
			outFile << vv_scalarFieldList[i][j].st_name << "\t";
			outFile << vv_scalarFieldList[i][j].p_blockMesh->n_elemNum << std::endl;

			for (size_t k = 0; k < vv_scalarFieldList[i][j].p_blockMesh->n_elemNum; k++)
			{
				outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << vv_scalarFieldList[i][j].elementField.GetValue(k)<<"\t";
			}
			outFile << std::endl;
		}
	}
}

void MHTVTKReader::ReadVTKFile(std::vector<std::string> vVTKFileName, std::vector<std::string>& vFiedNameList)
{
	std::cout << "start Read Field" << std::endl;

	v_FieldIO.resize(v_pmesh.size());

	if(vVTKFileName.size() != vv_scalarFieldList.size())
	{
		FatalError("mesh region number not same with vtk file number");
	}

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
			Field<Scalar> thisField(v_pmesh[i], 0.0, vFiedNameList[j]);
			vv_scalarFieldList[i].push_back(std::move(thisField));

			vv_scalarFieldList[i][j].ReadVTKGridField(Grid, vFiedNameList[j]);
		}

	}
	for (size_t i = 0; i < vVTKFileName.size(); i++)
	{
		for (size_t j = 0; j < vFiedNameList.size(); j++)
		{
			v_FieldIO[i].push_backScalarField(vv_scalarFieldList[i][j]);
		}
	}

	std::cout<<"read file succeed" << std::endl;
}

void MHTVTKReader::ReadVTKFile(std::vector<std::string>vVTKFileName, std::vector<int> vMeshID, std::vector<std::string>& vFiedNameList)
{
	v_meshID = vMeshID;
	std::cout << "start Read Field" << std::endl;

	v_FieldIO.resize(v_pmesh.size());

	if (vVTKFileName.size() != v_meshID.size())
	{
		FatalError("mesh ID number not same with vtk file number");
	}

	for (size_t i = 0; i < v_meshID.size(); i++)
	{
		if (v_meshID[i]>=v_pmesh.size())
		{
			FatalError("mesh ID out of range");
		}

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
			Field<Scalar> thisField(v_pmesh[v_meshID[i]], 0.0, vFiedNameList[j]);
			vv_scalarFieldList[v_meshID[i]].push_back(std::move(thisField));

			vv_scalarFieldList[v_meshID[i]][j].ReadVTKGridField(Grid, vFiedNameList[j]);
		}

	}
	for (size_t i = 0; i < v_meshID.size(); i++)
	{
		for (size_t j = 0; j < vFiedNameList.size(); j++)
		{
			v_FieldIO[v_meshID[i]].push_backScalarField(vv_scalarFieldList[v_meshID[i]][j]);
		}
	}

	std::cout << "read file succeed" << std::endl;
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

	std::ifstream inFile(DataFileName);

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
			vv_scalarFieldList[i].push_back(std::move(thisField));
		}
	}
	//do this must after v_scalarFieldList.push_back because point will change
	for (size_t i = 0; i < v_pmesh.size(); i++)
	{
		for (size_t j = 0; j < vFiedNameList.size(); j++)
		{
			v_FieldIO[i].push_backScalarField(vv_scalarFieldList[i][j]);
		}
	}
}

void MHTVTKReader::InitializeEmptyField(std::vector<std::string>& vFiedNameList)
{
	std::cout << "InitializeEmptyField Start" << std::endl;

	v_FieldIO.resize(v_pmesh.size());

	int nTotalFieldNum(v_pmesh.size() * vFiedNameList.size());

	for (size_t i = 0; i < v_pmesh.size(); i++)
	{
		for (size_t j = 0; j < vFiedNameList.size(); j++)
		{
			Field<Scalar> thisField(v_pmesh[i], 0.0, vFiedNameList[j]);

			vv_scalarFieldList[i].push_back(std::move(thisField));
		}
	}

	for (size_t i = 0; i < v_pmesh.size(); i++)
	{
		for (size_t j = 0; j < vFiedNameList.size(); j++)
		{
			v_FieldIO[i].push_backScalarField(vv_scalarFieldList[i][j]);
		}
	}
}

