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

MHTVTKReader::MHTVTKReader(std::vector<std::string>& vtkFileNameList, std::vector<std::string>& vFiedNameList, Scalar ratio)
{
	ReadMSHFromVTKFile(vtkFileNameList);
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
	ReadVTKFile(vtkFileNameList, vFiedNameList);
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

int GetFaceNum
(
	Mesh* pmesh,
	int elemID
)
{
	Element& curElement = pmesh->v_elem[elemID];
	if (curElement.est_shapeType == Element::ElemShapeType::estTetrahedral)
	{
		return 4;
	}
	else if (curElement.est_shapeType == Element::ElemShapeType::estHexahedral)
	{
		return 6;
	}
	else
	{
		FatalError("no more element shape type supported");
		return 0;
	}
}

void GetFaceVerticeList
(
	Mesh* pmesh,
	int elemID,
	std::vector<std::vector<int> >& verticeList
)
{
	Element& curElement = pmesh->v_elem[elemID];
	if (curElement.est_shapeType == Element::ElemShapeType::estTetrahedral)
	{
		verticeList.resize(4);
		for (int faceID = 0;faceID < 4;faceID++)
		{
			verticeList[faceID].clear();
		}
		//face 0
		verticeList[0].push_back(curElement.v_nodeID[0]);
		verticeList[0].push_back(curElement.v_nodeID[2]);
		verticeList[0].push_back(curElement.v_nodeID[1]);
		//face 1
		verticeList[1].push_back(curElement.v_nodeID[0]);
		verticeList[1].push_back(curElement.v_nodeID[1]);
		verticeList[1].push_back(curElement.v_nodeID[3]);
		//face 2
		verticeList[2].push_back(curElement.v_nodeID[0]);
		verticeList[2].push_back(curElement.v_nodeID[3]);
		verticeList[2].push_back(curElement.v_nodeID[2]);
		//face 3
		verticeList[3].push_back(curElement.v_nodeID[1]);
		verticeList[3].push_back(curElement.v_nodeID[2]);
		verticeList[3].push_back(curElement.v_nodeID[3]);
	}
	else if (curElement.est_shapeType == Element::ElemShapeType::estHexahedral)
	{
		verticeList.resize(6);
		for (int faceID = 0;faceID < 6;faceID++)
		{
			verticeList[faceID].clear();
		}
		//face 0
		verticeList[0].push_back(curElement.v_nodeID[0]);
		verticeList[0].push_back(curElement.v_nodeID[3]);
		verticeList[0].push_back(curElement.v_nodeID[2]);
		verticeList[0].push_back(curElement.v_nodeID[1]);
		//face 1
		verticeList[1].push_back(curElement.v_nodeID[0]);
		verticeList[1].push_back(curElement.v_nodeID[1]);
		verticeList[1].push_back(curElement.v_nodeID[5]);
		verticeList[1].push_back(curElement.v_nodeID[4]);

		//face 2
		verticeList[2].push_back(curElement.v_nodeID[0]);
		verticeList[2].push_back(curElement.v_nodeID[4]);
		verticeList[2].push_back(curElement.v_nodeID[7]);
		verticeList[2].push_back(curElement.v_nodeID[3]);

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

bool CompareFaceVertice
(
	std::vector<std::vector<int> >& verticeList,
	int ID1, 
	int ID2
)
{
	std::vector<int>& list1 = verticeList[ID1];
	std::vector<int>& list2 = verticeList[ID2];
	if (0 == list1.size() || 0 == list2.size()) return false;
	bool found = false;
	int startID = list1[0];
	//first, startID must be found in list2: if found, we can preset it true
	int startIDInList2 = -1;
	for (int i = 0;i < list2.size();i++)
	{
		if (startID == list2[i])
		{
			found = true;
			startIDInList2 = i;
			break;
		}
	}
	if (false == found) return false;
	//second, the following sequence must be the same: if any one is found different, it is then false
	for (int i = 1;i < list1.size();i++)
	{
		int verticeID1 = list1[i];
		//find the corresponding ID in list2 by the rule that IDs in list2 are inversedly saved
		int correspondingID = (startIDInList2 - i) % list2.size();
		int verticeID2 = list2[correspondingID];
		if (verticeID1 != verticeID2)
		{
			found = false;
			break;
		}
	}
	return found;
}

void InsertWhenNotFound(std::vector<int>& IDList, int newID)
{
	bool found = false;
	for (int j = 0;j < IDList.size();j++)
	{
		if (newID == IDList[j])
		{
			found = true;
			break;
		}
	}
	if (false == found)
	{
		IDList.push_back(newID);
	}
	return;
}

void MHTVTKReader::VTKCreateFaces(Mesh* pmesh)
{
	std::cout << "elemNum = " << pmesh->v_elem.size() << std::endl;
	std::cout << "faceNum = " << pmesh->v_face.size() << std::endl;
	std::cout << "nodeNum = " << pmesh->v_node.size() << std::endl;
	std::cout << "verticeNum = " << pmesh->v_vertice.size() << std::endl;
	int totalFaceNum = 0;
	for (int i = 0;i < pmesh->v_elem.size();i++)
	{
		totalFaceNum += GetFaceNum(pmesh, i);
	}
	std::cout << "in vtk file we have totally " << totalFaceNum << " faces" << std::endl;
	std::vector<std::vector<int> > faceVerticeList(totalFaceNum);
	std::vector<int> ownerElemIDList(totalFaceNum);
	std::vector<int> matchFaceIDList(totalFaceNum);
	std::vector<int> startIDs(pmesh->v_elem.size() + 1);
	for (int i = 0;i < totalFaceNum;i++)
	{
		matchFaceIDList[i] = -1;
	}
	int count = 0;
	for (int i = 0;i < pmesh->n_elemNum;i++)
	{
		startIDs[i] = count;
		int faceNum = GetFaceNum(pmesh, i);
		count += faceNum;
	}
	startIDs[pmesh->n_elemNum]= count;
	
	std::vector<std::vector<int> > localFaceVerticeID;
	for (int i = 0;i < pmesh->v_elem.size();i++)
	{
		int faceNum = GetFaceNum(pmesh, i);
		GetFaceVerticeList(pmesh, i, localFaceVerticeID);
		for (int j = startIDs[i];j < startIDs[i+1]; j++)
		{
			int verticeNum = localFaceVerticeID[j].size();
			//set corresponding element ID
			ownerElemIDList[j] = i;
			//writting vertice IDs
			faceVerticeList[j]= localFaceVerticeID[j- startIDs[i]];
		}
	}
	//search coincident faces and mark
	for (int i = 0;i < pmesh->v_vertice.size();i++)
	{
		std::vector<int> LocalfaceIDs;
		//collect local face IDs
		for (int j = 0;j < pmesh->v_vertice[i].v_elemID.size();j++)
		{
			int elemID = pmesh->v_vertice[i].v_elemID[j];
			for (int k = startIDs[elemID];k < startIDs[elemID + 1];k++)
			{
				LocalfaceIDs.push_back(k);
			}
		}
		int numToCompare = LocalfaceIDs.size();

		int elemNumInCompare = pmesh->v_vertice[i].v_elemID.size();
		/*
		std::cout << "number of cell in compare = " << elemNumInCompare << std::endl;
		std::cout << "number of faces in compare = " << numToCompare << std::endl;
		for (int j = 0;j < numToCompare;j++)
		{
			int localFaceID = LocalfaceIDs[j];
			for (int k = 0;k < faceVerticeList[localFaceID].size();k++)
			{
				std::cout << faceVerticeList[localFaceID][k] << "\t";
			}
			std::cout << std::endl;
		}
		*/
		//int numOfMatchFound = 0;
		for (int j = 0;j < numToCompare - 1;j++)
		{
			int ID1 = LocalfaceIDs[j];
			if (-1 != matchFaceIDList[ID1]) continue;
			for (int k = j + 1;k < numToCompare;k++)
			{
				int ID2 = LocalfaceIDs[k];
				if (-1 != matchFaceIDList[ID2]) continue;
				if (CompareFaceVertice(faceVerticeList, ID1, ID2))
				{
					matchFaceIDList[ID1] = ID2;
					matchFaceIDList[ID2] = ID1;
					//numOfMatchFound++;
				}
			}
		}

		//std::cout << "number match found = " << numOfMatchFound << std::endl;
		//system("pause");
	}
	//create faces in mesh, face IDs are also written in elements and FaceZones
	pmesh->v_face.clear();
	pmesh->v_boundaryFaceZone.resize(1);
	pmesh->v_boundaryFaceZone[0].name = "BOUNDARY";
	for (int i = 0;i < totalFaceNum;i++)
	{
		int twinID = matchFaceIDList[i];
		if (twinID == -1)
		{
			Face newFace;
			int ownerID = ownerElemIDList[i];
			newFace.v_nodeID = faceVerticeList[i];
			newFace.n_owner = ownerID;
			newFace.n_neighbor = -1;
			if (newFace.v_nodeID.size() == 3)
			{
				newFace.ft_faceType = Face::ftTrangular;
			}
			else if (newFace.v_nodeID.size() == 4)
			{
				newFace.ft_faceType = Face::ftQuadrilateral;
			}

			pmesh->v_face.push_back(newFace);
			int faceID = pmesh->v_face.size() - 1;
			InsertWhenNotFound(pmesh->v_elem[ownerID].v_faceID, faceID);
			pmesh->v_boundaryFaceZone[0].v_faceID.push_back(faceID);
		}
		else if (i < twinID)
		{
			Face newFace;
			int ownerID = ownerElemIDList[i];
			int nbID = ownerElemIDList[twinID];
			newFace.v_nodeID = faceVerticeList[i];
			newFace.n_owner = ownerID;
			newFace.n_neighbor = nbID;
			pmesh->v_face.push_back(newFace);
			int faceID = pmesh->v_face.size() - 1;
			InsertWhenNotFound(pmesh->v_elem[ownerID].v_faceID, faceID);
			InsertWhenNotFound(pmesh->v_elem[nbID].v_faceID, faceID);
			pmesh->fz_interiorFaceZone.v_faceID.push_back(faceID);
		}
	}
	pmesh->n_faceNum = pmesh->v_face.size();

	//write faceIDs in v_vertice
	for (int faceID = 0;faceID < pmesh->n_faceNum;faceID++)
	{
		for (int i = 0;i < pmesh->v_face[faceID].v_nodeID.size();i++)
		{
			int nodeID = pmesh->v_face[faceID].v_nodeID[i];
			InsertWhenNotFound(pmesh->v_vertice[nodeID].v_faceID, faceID);
		}
	}

	int isolatedFace = 0;
	for (int i = 0;i < totalFaceNum;i++)
	{
		int twinID = matchFaceIDList[i];
		if (twinID == -1)
		{
			isolatedFace++;
		}
		else
		{
			if (matchFaceIDList[twinID] != i)
			{
				FatalError("twins are not consistent");
			}
		}
	}
	std::cout << "isolatedFace number = " << isolatedFace << std::endl;
	std::cout << "face number in mesh = " << pmesh->v_face.size() << std::endl;
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

	for (size_t i = 0;i < v_pmesh.size();i++)
	{
		v_pmesh[i]->CalculateCenterAndVolume();
		Scalar totalVolume = 0;
		for (int j = 0;j < v_pmesh[i]->v_elem.size();j++)
		{
			totalVolume += v_pmesh[i]->v_elem[j].volume;
		}
		std::cout << "total volume of mesh #" << i << " = " << totalVolume << std::endl;
	}
	std::cout << "end of read vtk mesh file" << std::endl;
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
