#include "VTKIO.h"

MHTVTKReader::MHTVTKReader()
	:scaleRatio(1.0),translation(Vector(0.0,0.0,0.0))
{
}


MHTVTKReader::MHTVTKReader
(
	std::string MshFileName, 
	std::vector<std::string> vVTKFileName, 
	std::vector<std::string>& vFiedNameList,
	Scalar ratio
)
	:translation(Vector(0.0, 0.0, 0.0))
{
	ReadMSHFile(MshFileName);
	if (ratio < 0)
	{
		FatalError("invalid scale ratio: " + std::to_string(ratio));
	}
	else
	{
		this->scaleRatio = ratio;
	}
	this->GetTranslation();
	this->TransformToMOC();
	this->vv_scalarFieldList.resize(this->v_pmesh.size());
	ReadVTKFile(vVTKFileName, vFiedNameList);
}

MHTVTKReader::MHTVTKReader
(
	std::string MshFileName, 
	std::vector<std::string>& vFiedNameList,
	Scalar ratio
)
	:scaleRatio(1.0), translation(Vector(0.0, 0.0, 0.0))
{
	ReadMSHFile(MshFileName);
	if (ratio < 0)
	{
		FatalError("invalid scale ratio: " + std::to_string(ratio));
	}
	else
	{
		this->scaleRatio = ratio;
	}
	this->GetTranslation();
	this->TransformToMOC();
	this->vv_scalarFieldList.resize(this->v_pmesh.size());
	InitializeEmptyField(vFiedNameList);
}

MHTVTKReader::MHTVTKReader
(
	std::string MshFileName,
	Scalar ratio
)
	:translation(Vector(0.0, 0.0, 0.0))
{
	ReadMSHFile(MshFileName);
	if (ratio < 0)
	{
		FatalError("invalid scale ratio: " + std::to_string(ratio));
	}
	else
	{
		this->scaleRatio = ratio;
	}
	this->GetTranslation();
	this->TransformToMOC();
	this->vv_scalarFieldList.resize(this->v_pmesh.size());
}

MHTVTKReader::MHTVTKReader(std::vector<std::string>& vtkFileNameList, Scalar ratio)
	:translation(Vector(0.0, 0.0, 0.0))
{
	ReadMSHFromVTKFile(vtkFileNameList);
	//because of haven't face info,can't run these code
//	if (ratio < 0)
//	{
//		FatalError("invalid scale ratio: " + std::to_string(ratio));
//	}
//	else
//	{
//		this->scaleRatio = ratio;
//	}
//	this->GetTranslation();
//	this->TransformToMOC();
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

void MHTVTKReader::GetTranslation()
{
	Scalar xmin = 1e10;
	Scalar ymin = 1e10;
	Scalar zmin = 1e10;
	for (size_t i = 0;i < this->v_pmesh.size();i++)
	{
		Mesh* pmesh = v_pmesh[i];
		for (int j = 0;j < pmesh->v_node.size();j++)
		{
			Vector node = pmesh->v_node[j];
			xmin = Min(xmin, node.x_);
			ymin = Min(ymin, node.y_);
			zmin = Min(zmin, node.z_);
		}
	}
	this->translation = -Vector(xmin, ymin, zmin);
	return;
}

void MHTVTKReader::TransformToMOC()
{
	this->GetTranslation();
	for (size_t i = 0;i < this->v_pmesh.size();i++)
	{
		Mesh* pmesh = v_pmesh[i];
		pmesh->Translate(this->translation);
		pmesh->Scale(this->scaleRatio);
	}
	return;
}

void MHTVTKReader::TransformBack()
{
	for (size_t i = 0;i < this->v_pmesh.size();i++)
	{
		Mesh* pmesh = v_pmesh[i];
		pmesh->Scale(1.0 / this->scaleRatio);
		pmesh->Translate(-this->translation);
	}
	return;
}

int MHTVTKReader::GetIDOfRegion(std::string regionName)
{
	int ID = -1;
	for (size_t i = 0;i < this->v_pmesh.size();i++)
	{
		if (v_pmesh[i]->st_meshName == regionName)
		{
			ID = i;
			break;
		}
	}
	if (-1 == ID)
	{
		FatalError("region name " + regionName + " not found in CFD mesh");
	}
	return ID;
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
			if (!inFile.is_open())
			{
				FatalError(vVTKFileName[i] + " is not found");
			}
			ReadVTKMeshFormat(inFile);
			Field<Scalar> thisField(v_pmesh[meshID], 0.0, fieldNameList[j]);
			vv_scalarFieldList[meshID].push_back(std::move(thisField));
			vv_scalarFieldList[meshID][j].ReadVTKGridField(inFile);
			inFile.close();
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
	return;
}

void MHTVTKReader::VTKMeshEstablishVertice(Mesh* pmeh)
{
	int nodeNum = (int)pmeh->v_node.size();
	pmeh->v_vertice.resize(nodeNum);
	for (int i = 0; i < (int)pmeh->v_elem.size(); i++)
	{
		//visiting one element
		for (int j = 0; j < (int)pmeh->v_elem[i].v_nodeID.size(); j++)
		{
			int verticeNum = pmeh->v_elem[i].v_nodeID[j];
			pmeh->v_vertice[verticeNum].v_elemID.push_back(i);
		}
	}	
}

void GetFaceVerticeList
(
	Mesh* pmesh,
	int elemID,
	std::vector<std::vector<int> >& verticeList
)
{
	//假设网格没有存储face

	Element curElement = pmesh->v_elem[elemID];

	if (curElement.est_shapeType == Element::ElemShapeType::estTetrahedral)
	{
		verticeList.resize(4);
		//face 0
		verticeList[0].push_back(curElement.v_nodeID[0]);
		verticeList[0].push_back(curElement.v_nodeID[1]);
		verticeList[0].push_back(curElement.v_nodeID[2]);
		//face 1
		verticeList[1].push_back(curElement.v_nodeID[0]);
		verticeList[1].push_back(curElement.v_nodeID[1]);
		verticeList[1].push_back(curElement.v_nodeID[3]);
		//face 2
		verticeList[2].push_back(curElement.v_nodeID[0]);
		verticeList[2].push_back(curElement.v_nodeID[2]);
		verticeList[2].push_back(curElement.v_nodeID[3]);
		//face 3
		verticeList[3].push_back(curElement.v_nodeID[1]);
		verticeList[3].push_back(curElement.v_nodeID[2]);
		verticeList[3].push_back(curElement.v_nodeID[3]);
	}
	else if (curElement.est_shapeType == Element::ElemShapeType::estHexahedral)
	{
		verticeList.resize(6);
		//face 0
		verticeList[0].push_back(curElement.v_nodeID[0]);
		verticeList[0].push_back(curElement.v_nodeID[1]);
		verticeList[0].push_back(curElement.v_nodeID[2]);
		verticeList[0].push_back(curElement.v_nodeID[3]);
		//face 1
		verticeList[1].push_back(curElement.v_nodeID[0]);
		verticeList[1].push_back(curElement.v_nodeID[1]);
		verticeList[1].push_back(curElement.v_nodeID[5]);
		verticeList[1].push_back(curElement.v_nodeID[4]);

		//face 2
		verticeList[2].push_back(curElement.v_nodeID[0]);
		verticeList[2].push_back(curElement.v_nodeID[3]);
		verticeList[2].push_back(curElement.v_nodeID[7]);
		verticeList[2].push_back(curElement.v_nodeID[4]);

		//face 3
		verticeList[3].push_back(curElement.v_nodeID[1]);
		verticeList[3].push_back(curElement.v_nodeID[2]);
		verticeList[3].push_back(curElement.v_nodeID[6]);
		verticeList[3].push_back(curElement.v_nodeID[5]);

		//face 4
		verticeList[4].push_back(curElement.v_nodeID[2]);
		verticeList[4].push_back(curElement.v_nodeID[3]);
		verticeList[4].push_back(curElement.v_nodeID[7]);
		verticeList[4].push_back(curElement.v_nodeID[6]);

		//face 5
		verticeList[5].push_back(curElement.v_nodeID[4]);
		verticeList[5].push_back(curElement.v_nodeID[5]);
		verticeList[5].push_back(curElement.v_nodeID[6]);
		verticeList[5].push_back(curElement.v_nodeID[7]);

	}
	return;
}

void MHTVTKReader::VTKCreateFaces(Mesh* pmesh)
{
	std::cout << "elemNum = " << pmesh->v_elem.size() << std::endl;
	std::cout << "faceNum = " << pmesh->v_face.size() << std::endl;
	std::cout << "nodeNum = " << pmesh->v_node.size() << std::endl;
	std::cout << "verticeNum = " << pmesh->v_vertice.size() << std::endl;
	/*
	for (int i = 0;i < pmesh->v_elem.size();i++)
	{
		std::cout << "element #" << i << " is composed of faces" << std::endl;
		for (int j = 0;j < pmesh->v_elem[i].v_faceID.size();j++)
		{
			std::cout << pmesh->v_elem[i].v_faceID[j] << "\t";
		}
		std::cout << std::endl;
		std::cout << "and nodes" << std::endl;
		for (int j = 0;j < pmesh->v_elem[i].v_nodeID.size();j++)
		{
			std::cout << pmesh->v_elem[i].v_nodeID[j] << "\t";
		}
		std::cout << std::endl;
		system("pause");
	}
	*/
	for (int i = 0;i < pmesh->v_vertice.size();i++)
	{
		if (pmesh->v_vertice[i].v_elemID.size() <= 4) continue;
		std::cout << "node #" << i << " is composed of faces" << std::endl;
		for (int j = 0;j < pmesh->v_vertice[i].v_faceID.size();j++)
		{
			std::cout << pmesh->v_vertice[i].v_faceID[j] << "\t";
		}
		std::cout << std::endl;
		std::cout << "and elements" << std::endl;
		for (int j = 0;j < pmesh->v_vertice[i].v_elemID.size();j++)
		{
			std::cout << pmesh->v_vertice[i].v_elemID[j] << "\t";
		}
		std::cout << std::endl;
		system("pause");
	}
	return;
}

void MHTVTKReader::ReadMSHFromVTKFile(std::vector<std::string>& vtkFileNameList)
{
	v_stdMesh.resize(vtkFileNameList.size());
	for (size_t i = 0; i < vtkFileNameList.size(); i++)
	{
		v_pmesh.push_back(&v_stdMesh[i]);
		std::ifstream inFile(vtkFileNameList[i]);
		ReadVTKMeshFormat(inFile, v_pmesh[i]);
		this->VTKMeshEstablishVertice(v_pmesh[i]);
		this->VTKCreateFaces(v_pmesh[i]);
	}
	std::cout<<"end of read vtk mesh file" << std::endl;
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

void MHTVTKReader::ReadVTKMeshFormat(std::ifstream& inFile, Mesh* pmesh)
{
	pmesh->md_meshDim = Mesh::md3D;

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
			pmesh->n_nodeNum = pointNum;
			pmesh->v_node.resize(pointNum);
			for (size_t i = 0; i < pointNum; i++)
			{
				Scalar point_x, point_y, point_z;
				inFile >> point_x >> point_y >> point_z;
				
				pmesh->v_node[i].x_ = point_x;
				pmesh->v_node[i].y_ = point_y;
				pmesh->v_node[i].z_ = point_z;
			}
		}
		else if (dataComment == "CELLS")
		{
		
			int nCellNum, nCellTotalNum;
			inFile >> nCellNum >> nCellTotalNum;
			pmesh->n_elemNum = nCellNum;
			pmesh->v_elem.resize(nCellNum);

			for (size_t i = 0; i < nCellNum; i++)
			{
				int nCellNodeNum, nNodeID;
				inFile >> nCellNodeNum;
				pmesh->v_elem[i].v_nodeID.resize(nCellNodeNum);
				for (size_t j = 0; j < nCellNodeNum; j++)
				{
					inFile >> nNodeID;
					pmesh->v_elem[i].v_nodeID[j] = nNodeID;
				}
			}
		}
		else if (dataComment == "CELL_TYPES")
		{

			int nCellTypeNum, nCellType,firstCellType;
			bool ifMix(false);
			inFile >> nCellTypeNum;
		
			for (size_t i = 0; i < nCellTypeNum; i++)
			{
				inFile >> nCellType;

				if (i==0)
				{
					firstCellType = nCellType;
				}
				else
				{
					if (firstCellType != nCellType)
					{
						ifMix = true;
					}
				}

				if (nCellType == 5)
				{
					pmesh->v_elem[i].est_shapeType = Element::ElemShapeType::estTetrahedral;
				}
				else if (nCellType == 9)
				{
					pmesh->v_elem[i].est_shapeType = Element::ElemShapeType::estQuadrilateral;
				}
				else if (nCellType == 10)
				{
					pmesh->v_elem[i].est_shapeType = Element::ElemShapeType::estTetrahedral;
				}
				else if (nCellType == 14)
				{
					pmesh->v_elem[i].est_shapeType = Element::ElemShapeType::estPyramid;
				}
				else if (nCellType == 13)
				{
					pmesh->v_elem[i].est_shapeType = Element::ElemShapeType::estWedge;
				}
				else if (nCellType == 12)
				{
					pmesh->v_elem[i].est_shapeType = Element::ElemShapeType::estHexahedral;
				}
			}
		
			if (ifMix)
			{
				pmesh->est_shapeType = Element::ElemShapeType::estMixed;
			}
			else
			{
				if (nCellType == 5)
				{
					pmesh->est_shapeType = Element::ElemShapeType::estTetrahedral;
				}
				else if (nCellType == 9)
				{
					pmesh->est_shapeType = Element::ElemShapeType::estQuadrilateral;
				}
				else if (nCellType == 10)
				{
					pmesh->est_shapeType = Element::ElemShapeType::estTetrahedral;
				}
				else if (nCellType == 14)
				{
					pmesh->est_shapeType = Element::ElemShapeType::estPyramid;
				}
				else if (nCellType == 13)
				{
					pmesh->est_shapeType = Element::ElemShapeType::estWedge;
				}
				else if (nCellType == 12)
				{
					pmesh->est_shapeType = Element::ElemShapeType::estHexahedral;
				}
			}

			return;
		}
		else
		{
			break;
		}
	}
}
