#include "VTKIO.h"

MHTVTKReader::MHTVTKReader()
{
}


MHTVTKReader::MHTVTKReader(std::string MshFileName, std::vector<std::string> vVTKFileName, std::vector<std::string>& vFiedNameList)
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
	std::vector<Mesh*>().swap(v_pmesh);
	std::vector<int>().swap(v_meshID);
}

void MHTVTKReader::WriteDataFile(std::string DataFileName)
{

	std::ofstream outFile(DataFileName);

	for (size_t i = 0; i < this->vv_scalarFieldList.size(); i++)
	{
		for (size_t j = 0; j < this->vv_scalarFieldList[i].size(); j++)
		{
			outFile << vv_scalarFieldList[i][j].p_blockMesh->st_meshName << "\t";
			outFile << vv_scalarFieldList[i][j].st_name << "\t";
			outFile << vv_scalarFieldList[i][j].p_blockMesh->n_elemNum << std::endl;

			for (size_t k = 0; k < vv_scalarFieldList[i][j].p_blockMesh->n_elemNum; k++)
			{
				outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << vv_scalarFieldList[i][j].elementField.GetValue(k) << "\t";
			}
			outFile << std::endl;
		}
	}
}

void MHTVTKReader::ReadVTKFile(std::vector<std::string> vVTKFileName, std::vector<std::string>& vFiedNameList)
{

	std::cout << "start Read Field" << std::endl;

	v_FieldIO.resize(v_pmesh.size());

	if (vVTKFileName.size() != v_pmesh.size())
	{
		FatalError("mesh number not same with vtk file number");
	}

	for (size_t i = 0; i < v_pmesh.size(); i++)
	{
		for (size_t j = 0; j < vFiedNameList.size(); j++)
		{
			std::ifstream inFile(vVTKFileName[i]);
			ReadVTKMeshFormat(inFile);
			Field<Scalar> thisField(v_pmesh[i], 0.0, vFiedNameList[j]);
			vv_scalarFieldList[i].push_back(std::move(thisField));
			vv_scalarFieldList[i][j].ReadVTKGridField(inFile);
		}

	}
	std::cout << "succeed" << std::endl;

	for (size_t i = 0; i < v_pmesh.size(); i++)
	{
		for (size_t j = 0; j < vFiedNameList.size(); j++)
		{
			v_FieldIO[i].push_backScalarField(vv_scalarFieldList[i][j]);
		}
	}


	std::cout<<"read file succeed" << std::endl;
}

void MHTVTKReader::ReadVTKFile(std::vector<std::string>vVTKFileName, std::vector<int> vMeshID, std::vector<std::string>& fieldNameList)
{
	v_meshID = vMeshID;
	std::cout << "Reading Fields from vtk files" << std::endl;

	v_FieldIO.resize(v_pmesh.size());

	if (vVTKFileName.size() != v_meshID.size())
	{
		FatalError("number of mesh IDs is not the same with that of vtk files");
	}

	for (size_t i = 0; i < v_meshID.size(); i++)
	{
		int meshID = v_meshID[i];
		if (meshID >=v_pmesh.size())
		{
			FatalError("mesh ID " + std::to_string(meshID) + " is out of range");
		}

		for (size_t j = 0; j < fieldNameList.size(); j++)
		{
			std::ifstream inFile(vVTKFileName[i]);
			ReadVTKMeshFormat(inFile);
			Field<Scalar> thisField(v_pmesh[meshID], 0.0, fieldNameList[j]);
			vv_scalarFieldList[meshID].push_back(std::move(thisField));
			vv_scalarFieldList[meshID][j].ReadVTKGridField(inFile);
		}
	}
	std::cout<<"succeed" << std::endl;

	for (size_t i = 0; i < v_meshID.size(); i++)
	{
		int meshID = v_meshID[i];
		for (size_t j = 0; j < fieldNameList.size(); j++)
		{
			v_FieldIO[meshID].push_backScalarField(vv_scalarFieldList[meshID][j]);
		}
	}

	std::cout << "read file succeed" << std::endl;
}

void MHTVTKReader::ReadMSHFile(std::string MeshFileName)
{
	std::cout << "start Read Mesh" << std::endl;
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

			if (vFiedNameList[j] != sFieldName)
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

void MHTVTKReader::ReadVTKMeshFormat(std::ifstream& inFile)
{
	std::string comment;
	getline(inFile, comment);


	std::string WriterName;
	getline(inFile, WriterName);


	std::string codingFormat;
	inFile >> codingFormat;


	std::string DataSet;
	inFile >> DataSet;
	inFile >> DataSet;


	while (true)
	{
		std::string dataComment;
		inFile >> dataComment;

		if (dataComment == "POINTS")
		{

			int pointNum;
			std::string  pointNumType;

			inFile >> pointNum >> pointNumType;

			for (size_t i = 0; i < pointNum; i++)
			{
				Scalar point_x, point_y, point_z;
				inFile >> point_x>> point_y>> point_z;
			}
		}
		else if (dataComment == "CELLS")
		{

			int nCellNum,nCellTotalNum;
			inFile >> nCellNum >> nCellTotalNum;

			for (size_t i = 0; i < nCellNum; i++)
			{
				int nCellNodeNum,nNodeID;
				inFile >> nCellNodeNum;
				for (size_t j = 0; j < nCellNodeNum; j++)
				{
					inFile >> nNodeID;
				}
			}
		}
		else if (dataComment == "CELL_TYPES")
		{

			int nCellTypeNum,nCellType;
			inFile >> nCellTypeNum;
			for (size_t i = 0; i < nCellTypeNum; i++)
			{
				inFile >> nCellType;
			}
			return;
		}
		else
		{
			break;
		}
	}
}
