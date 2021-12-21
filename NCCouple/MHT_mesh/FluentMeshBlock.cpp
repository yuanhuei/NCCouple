/*---------------------------------------------------------------------------*\
File Name:
	FluentMeshBlock.h

Description:
	Derivation of Mesh to make it suited to FLUENT mesh file

	Author:		Shuai Zhang, Kong Ling
	Date: once upon a time (I can't remember)

Revised:
	Description:
	1. Display functions were removed;
	2. Styles of information display during mesh reading were entirely checked and revised,
	where logLevels and systemControl were applied.

	Revisor:		Kong Ling, Shuai Zhang
	Modified Date:	2017-01-14
	
Revised:
	Description:
	1. Modification was made in ReadMesh(F_FilePointer);
	2. Convert symbol of "/n" to "/r" and "/n" in "ReadNew GetData" function;

	Revisor:		ZiJian Li
	Modified Date:	2018-11-14
\*---------------------------------------------------------------------------*/

#include "FluentMeshBlock.h"
#include <iomanip>
#include "../MHT_common/LogLevel.h"
#include "../MHT_common/StringTools.h"
#include "../MHT_common/SystemControl.h"
#include "../MHT_mesh/StandardMeshBlock.h"
#include "../MHT_mesh/RegionConnection.h"

FluentMeshFileSection::FluentMeshFileSection()
{
    n_index = -1;
    ss_details.str("");
    ss_details.clear();
}
int FluentMeshFileSection::ReadNew(std::ifstream& inFile)
{
    n_index = -1;
    ss_details.str("");
    ss_details.clear();
    std::stringstream sstream;
    if (1 == GetBySymbol(inFile, sstream,'(',')'))
    {
        return 1;
    }
    sstream >> n_index;
    if (0 == n_index)
    {
        GetBySymbol(sstream, ss_details, '\"', '\"');
    }
    else if (2 == n_index)
    {
        v_parameters.resize(1);
        sstream >> v_parameters[0];
    }
    else if (10 == n_index || 12 == n_index || 13 == n_index)
    {
        std::stringstream sstreamTemp;
        GetBySymbol(sstream, sstreamTemp,'(',')');
        v_parameters.resize(5);
        for (int i = 0; i < 5; i++)
        {
            sstreamTemp >> std::hex >> v_parameters[i];
        }
        if (0 != v_parameters[0])
        {
            GetBySymbol(sstream, ss_details,'(',')');
        }
    }
    else if (39 == n_index || 45 == n_index)
    {
        GetBySymbol(sstream, ss_details,'(',')');
    }
    return 0;
}

FluentMeshBlock::FluentMeshBlock(const std::string& MshFileName)
	:Mesh(MshFileName)
{
    std::ifstream	inFile(MshFileName.c_str());
	FILE *F_FilePointer;
	F_FilePointer = fopen(MshFileName.c_str(),"r");
	if (inFile.fail())
	{
        FatalError(std::string("Can not find Mesh file:\n") + MshFileName);
	}
	//read .msh file
    //this->ReadMesh(inFile);
    this->ReadMesh(F_FilePointer);
	//establish zones in FLUENT mesh
	//Get index from mesh
	this->EstablishZones();
	//establish topology information regarding vertices
    this->EstablishVertice();
	//Calculate numbers of elements, faces and vertices in this mesh
	this->CalculateParameters();
	//Calculate geometrical parameters
	this->CalculateCenterAndVolume();

}
void FluentMeshBlock::ReadMesh(FILE *F_FilePointer)
{
    std::string s_Index;
	char temp=' ';
	s_Index = fgetc(F_FilePointer);

//	int mCharLocation(0), mStringStart(0);

	while (1)
	{
		while (s_Index != "(")							//erase space before "("
		{
			s_Index = fgetc(F_FilePointer);
			temp = s_Index[0];
			if (temp == -1)
			{
                std::cout <<"run end of file" << std::endl;
				break;
			}
		}
		if (temp == -1)
			break;
		switch (GetIndex(F_FilePointer))
		{
		case 0:
		{			
			ReadComment(F_FilePointer);								
			break;
		}
		case 2:
		{
			ReadDimensionNumber(F_FilePointer);
			break;
		}
		case 10:
		{
			ReadNode(F_FilePointer);
			break;
		}
		case 12:
		{
			ReadElement(F_FilePointer);
			break;
		}
		case 13:
		{
			ReadFace(F_FilePointer);
			break;
		}
		case 39:
		{
			ReadZone(F_FilePointer);
			break;
		}
		case 45:
		{
			ReadZone(F_FilePointer);
			break;
		}
		}
		s_Index = "";
	}
	RefineArrayIndex();

	PopulateCellNodes();
}
void FluentMeshBlock::ReadMesh(std::ifstream& inFile)
{
	FluentMeshFileSection fileSection;
	while (0 == fileSection.ReadNew(inFile))
	{
		if (0 == fileSection.n_index)
		{
			ReadComment(fileSection);
		}
		else if (2 == fileSection.n_index)
		{
			ReadDimensionNumber(fileSection);
		}
		else if (10 == fileSection.n_index)
		{
			ReadNode(fileSection);
		}
		else if (12 == fileSection.n_index)
		{
			ReadElement(fileSection);
		}
		else if (13 == fileSection.n_index)
		{
			ReadFace(fileSection);
		}
		else if (39 == fileSection.n_index || 45 == fileSection.n_index)
		{
			ReadZone(fileSection);
		}
	}

	RefineArrayIndex();

	PopulateCellNodes();
}

int FluentMeshBlock::GetIndex(FILE *F_FilePointer)
{
	int n_Index;
	fscanf(F_FilePointer, "%d", &n_Index);
	return n_Index;
}

void FluentMeshBlock::ReadComment(FluentMeshFileSection& inSection)
{
	std::cout << LogLevel(LogLevel::mlInfo, inSection.ss_details.str()) << std::endl;
}
void FluentMeshBlock::ReadComment(FILE *F_FilePointer)
{
	char c_Temp;
	c_Temp = fgetc(F_FilePointer);
	
	while (c_Temp != ')')
	{
        if (c_Temp != '\r'&&c_Temp != '\"'&&c_Temp !='\n')									//erase effect of "Enter"
		{
            std::cout << c_Temp;
		}
		c_Temp = fgetc(F_FilePointer);
	}
    std::cout << std::endl;
	
}
void FluentMeshBlock::ReadDimensionNumber(FILE *F_FilePointer)
{
	char c_temp;
	fscanf(F_FilePointer, "%d%c", &this->md_meshDim, &c_temp);
	if (c_temp != ')')
	{
        FatalError(std::string(") is missing."));
	}
	std::cout << LogLevel(LogLevel::mlInfo, "Mesh demension number is") << this->md_meshDim << std::endl;
}
void FluentMeshBlock::ReadDimensionNumber(FluentMeshFileSection& inSection)
{
	if (0 == inSection.v_parameters.size())
	{
        FatalError(std::string("Dimension number is missing."));
	}

	this->md_meshDim = (Mesh::MeshDim)inSection.v_parameters[0];

	if (this->md_meshDim != this->md2D && this->md_meshDim != this->md3D)
	{
        FatalError(std::string("Mesh dimension number is incorrect."));
	}

	std::cout << LogLevel(LogLevel::mlInfo, "Mesh demension number is") << this->md_meshDim << std::endl << std::endl;
}

Scalar FluentMeshBlock::GetData(FILE *F_FilePointer, std::string &s_InputData, int n_DataType)
{
    std::string s_OutPutData;
	s_OutPutData = "";
	if (n_DataType==1)
	{ 
    while (s_InputData == "\t" || s_InputData == "\r" || s_InputData == " "|| s_InputData == "\n")
	{
		s_InputData = fgetc(F_FilePointer);
	}
    while (s_InputData != "\t" && s_InputData != "\r" && s_InputData != " "&&s_InputData != ")"&&s_InputData != "\n")
	{
		s_OutPutData = s_OutPutData + s_InputData;
		s_InputData = fgetc(F_FilePointer);
	}
	}
    return std::stod(s_OutPutData);
}
int FluentMeshBlock::GetData(FILE *F_FilePointer, std::string &s_InputData)
{
    std::string s_OutPutData;
	s_OutPutData = "";
	int n_Dst(0), n_Temp;

    while (s_InputData == "\t" || s_InputData == "\r" || s_InputData == " "|| s_InputData == "\n")
	{
		s_InputData = fgetc(F_FilePointer);
	}
    while (s_InputData != "\t" && s_InputData != "\r" && s_InputData != " "&&s_InputData != ")"&&s_InputData != "\n")
	{
		s_OutPutData = s_OutPutData + s_InputData;
		s_InputData = fgetc(F_FilePointer);
	}

    for (int i = 0; i < (int)s_OutPutData.length(); i++)
	{
		if (s_OutPutData[i] <= '9')
		{
			n_Temp = s_OutPutData[i] - '0';
		}
		else
			n_Temp = s_OutPutData[i] - 'a' + 10;
		n_Dst = n_Dst*16 + n_Temp;
	}
		return n_Dst;
}

void FluentMeshBlock::ReadNew(FILE *F_FilePointer, int n_Index)
	{
		int n_Dst(0), n_Temp;
        std::string s_Temp,s_TotalChar;
		v_parameters.resize(5);
		if (10 == n_Index )
		{
			for (int i = 0; i < 5; i++)
			{
				fscanf(F_FilePointer, "%x", &v_parameters[i]);
			}
		}
		else if (12 == n_Index || 13 == n_Index)
		{
            std::string v_parameters_Data;
			for (int i = 0; i < 4; i++)
			{
				fscanf(F_FilePointer, "%x", &v_parameters[i]);
			}
			s_Temp = fgetc(F_FilePointer);
            while (s_Temp == "\t" || s_Temp == "\r" || s_Temp == " "||s_Temp == "\n")
			{
				s_Temp = fgetc(F_FilePointer);
			}
			if (")" != s_Temp)
			{
				s_TotalChar = s_TotalChar + s_Temp;
				s_Temp = fgetc(F_FilePointer);
                while (s_Temp != "\t" && s_Temp != "\r" && s_Temp != " "&&s_Temp != ")"&&s_Temp != "\n")
				{
					s_TotalChar = s_TotalChar + s_Temp;
				}
				for (int i = 0; i < s_TotalChar.length(); i++)
				{
					if (s_TotalChar[i] <= '9')
					{
						n_Temp = s_TotalChar[i] - '0';
					}
					else
						n_Temp = s_TotalChar[i] - 'a' + 10;
					n_Dst = n_Dst * 16 + n_Temp;
				}
				v_parameters[4] = n_Dst;
				/*cout << v_parameters [4]<< endl;
				cout << s_TotalChar << endl;*/
			}
			else
				v_parameters[4] = 0;
		}
	}

void FluentMeshBlock::ReadNode(FluentMeshFileSection& inSection)
{
	if (inSection.v_parameters.size() != 5)
	{
        FatalError(std::string("Incorrect parameter number detected in node section."));
	}

	if (0 != inSection.v_parameters[0])
	{
		std::cout << LogLevel(LogLevel::mlOK, "Reading node coordinates") << std::endl;

		int beginNum = inSection.v_parameters[1];
		int endNum = inSection.v_parameters[2];

		this->n_nodeNum += endNum - beginNum + 1;

		if (this->md_meshDim == Mesh::md2D)
		{
			for (int i = beginNum; i <= endNum; i++)
			{
				Node nodeTemp;
				inSection.ss_details >> nodeTemp.x_ >> nodeTemp.y_;
				this->v_node.push_back(nodeTemp);
			}
		}
		else if (this->md_meshDim == Mesh::md3D)
		{
			for (int i = beginNum; i <= endNum; i++)
			{
				Node nodeTemp;
				inSection.ss_details >> nodeTemp.x_ >> nodeTemp.y_ >> nodeTemp.z_;
				this->v_node.push_back(nodeTemp);
			}
		}
	}

}
void FluentMeshBlock::ReadNode(FILE *F_FilePointer)
{
	
    std::string s_leftChar;
	int n_ZoneId,
		n_FirstIndex,
		n_LastIndex, n_Type,
		n_dimensionalityOfGrid;
	s_leftChar = fgetc(F_FilePointer);
	while (s_leftChar != "(")							//erase space before "("
	{
		s_leftChar = fgetc(F_FilePointer);
	}
	ReadNew(F_FilePointer,10);
	if (v_parameters.size() != 5)
	{
        FatalError(std::string("Incorrect parameter number detected in node section."));
	}
	n_ZoneId = v_parameters[0];
	n_FirstIndex = v_parameters[1];
	n_LastIndex = v_parameters[2];
	n_Type = v_parameters[3];
	n_dimensionalityOfGrid = v_parameters[4];
    /*cout << "n_ZoneId: " << n_ZoneId << endl;
	cout << "n_FirstIndex: " << n_FirstIndex << endl;
	cout << "n_LastIndex: " << n_LastIndex << endl;
	cout << "n_Type: " << n_Type << endl;
    cout << "n_dimensionalityOfGrid: " << n_dimensionalityOfGrid << endl;*/
    if (0 != n_ZoneId)
	{
		std::cout << LogLevel(LogLevel::mlOK, "Reading node coordinates") << std::endl;
		s_leftChar = fgetc(F_FilePointer);
		while (s_leftChar != "(")							//erase space before "("
		{
			s_leftChar = fgetc(F_FilePointer);
		}
		s_leftChar = fgetc(F_FilePointer);

		switch (n_dimensionalityOfGrid)
		{
			case Mesh::md2D:
            {
				for (int i = 0; i < (n_LastIndex - n_FirstIndex + 1); i++)
				{
                    Node nodeTemp;
					nodeTemp.x_=GetData(F_FilePointer, s_leftChar, 1);
					nodeTemp.y_=GetData(F_FilePointer, s_leftChar, 1);
					this->v_node.push_back(nodeTemp);
                //	system("pause");
                }

				break;
			}
			case Mesh::md3D:
			{
				for (int i = 0; i < (n_LastIndex - n_FirstIndex + 1); i++)
				{
					Node nodeTemp;
					nodeTemp.x_ = GetData(F_FilePointer, s_leftChar, 1);
					nodeTemp.y_ = GetData(F_FilePointer, s_leftChar, 1);
					nodeTemp.z_ = GetData(F_FilePointer, s_leftChar, 1);
					this->v_node.push_back(nodeTemp);
				//	system("pause");
				}
				break;
			}
		}

	}
}

void FluentMeshBlock::ReadElement(FluentMeshFileSection& inSection)
{
	if (inSection.v_parameters.size() != 5)
	{
        FatalError(std::string("Incorrect parameter number detected in element section."));
	}
	if (0 == inSection.v_parameters[0])
	{
		std::cout << LogLevel(LogLevel::mlOK, "Reading element declaration") << std::endl;
		int beginNum = inSection.v_parameters[1];
		int endNum = inSection.v_parameters[2];
		this->n_elemNum = endNum - beginNum + 1;
	}
	else
	{
		std::cout << LogLevel(LogLevel::mlOK, "Reading shape types of the elements") << std::endl;

		MshElementZone czTemp;
		czTemp.n_zoneID = inSection.v_parameters[0];
		czTemp.n_begID = inSection.v_parameters[1];
		czTemp.n_endID = inSection.v_parameters[2];
		czTemp.n_active = inSection.v_parameters[3];
		czTemp.n_elemShape = inSection.v_parameters[4];
		czTemp.n_elemNum = czTemp.n_endID - czTemp.n_begID + 1;

		czTemp.est_shapeType = (Element::ElemShapeType)czTemp.n_elemShape;

        if (this->est_shapeType == Element::estNoType)
		{
			this->est_shapeType = (Element::ElemShapeType)czTemp.n_elemShape;
		}

		v_ElementZoneInfo.push_back(czTemp);

        if (czTemp.est_shapeType == Element::estMixed)
		{
			for (int i = 0; i < czTemp.n_elemNum; i++)
			{
				Element elemTemp;
				int	nElemType(-1);
				inSection.ss_details >> nElemType;
				elemTemp.est_shapeType = (Element::ElemShapeType)nElemType;
				this->v_elem.push_back(elemTemp);
			}
		}
		else
		{
			for (int i = 0; i < czTemp.n_elemNum; i++)
			{
				Element elemTemp;
				elemTemp.est_shapeType = (Element::ElemShapeType)czTemp.n_elemShape;
				this->v_elem.push_back(elemTemp);
			}
		}
	}
}
void FluentMeshBlock::ReadElement(FILE *F_FilePointer)
{
    std::string s_leftChar;
	while (s_leftChar != "(")							//erase space before "("
	{
		s_leftChar = fgetc(F_FilePointer);
	}
	ReadNew(F_FilePointer, 12);
	if (v_parameters.size() != 5)
	{
        FatalError(std::string("Incorrect parameter number detected in element section."));
	}
	if (0 == v_parameters[0])
	{
		std::cout << LogLevel(LogLevel::mlOK, "Reading element declaration") << std::endl;
		int beginNum = v_parameters[1];
		int endNum = v_parameters[2];
		this->n_elemNum = endNum - beginNum + 1;
	}
	else
	{
		std::cout << LogLevel(LogLevel::mlOK, "Reading shape types of the elements") << std::endl;

		MshElementZone czTemp;
		czTemp.n_zoneID = v_parameters[0];
		czTemp.n_begID = v_parameters[1];
		czTemp.n_endID = v_parameters[2];
		czTemp.n_active = v_parameters[3];
		czTemp.n_elemShape = v_parameters[4];
		czTemp.n_elemNum = czTemp.n_endID - czTemp.n_begID + 1;

		czTemp.est_shapeType = (Element::ElemShapeType)czTemp.n_elemShape;

		if (this->est_shapeType == Element::estNoType)
		{
			this->est_shapeType = (Element::ElemShapeType)czTemp.n_elemShape;
		}

		v_ElementZoneInfo.push_back(czTemp);
		
		if (czTemp.est_shapeType == Element::estMixed)
		{
            std::cout << "run here" << std::endl;
			s_leftChar = fgetc(F_FilePointer);
			while (s_leftChar != "(")							//erase space before "("
			{
				s_leftChar = fgetc(F_FilePointer);
			}
			s_leftChar = fgetc(F_FilePointer);
			for (int i = 0; i < czTemp.n_elemNum; i++)
			{
				Element elemTemp;
				int	nElemType(-1);
				nElemType=GetData(F_FilePointer, s_leftChar);
				elemTemp.est_shapeType = (Element::ElemShapeType)nElemType;
				this->v_elem.push_back(elemTemp);
			}
		}
		else
		{
			for (int i = 0; i < czTemp.n_elemNum; i++)
			{
				Element elemTemp;
				elemTemp.est_shapeType = (Element::ElemShapeType)czTemp.n_elemShape;
				this->v_elem.push_back(elemTemp);
			}
        }
	}
}


void FluentMeshBlock::ReadFace(FluentMeshFileSection& inSection)
{
	if (inSection.v_parameters.size() != 5)
	{
        FatalError(std::string("Incorrect parameter number detected in face section."));
	}
	if (0 == inSection.v_parameters[0])
	{
		std::cout << LogLevel(LogLevel::mlOK, "Reading face declaration") << std::endl;
		int beginNum = inSection.v_parameters[1];
		int endNum = inSection.v_parameters[2];
		this->n_faceNum = endNum - beginNum + 1;
	}
	else
	{
		std::cout << LogLevel(LogLevel::mlOK, "Reading face topologies") << std::endl;
		MshFaceZone fzTemp;

		//(zone-id first-n_index last-n_index type element-type)
		//type : mixed,linear,triangular,quadrilateral
		fzTemp.n_zoneID = inSection.v_parameters[0];
		fzTemp.n_begID = inSection.v_parameters[1];
		fzTemp.n_endID = inSection.v_parameters[2];
		fzTemp.n_BCType = inSection.v_parameters[3];
		fzTemp.n_faceShape = inSection.v_parameters[4];

		fzTemp.n_faceNum = fzTemp.n_endID - fzTemp.n_begID + 1;

		v_FaceZoneInfo.push_back(fzTemp);

		//Read face
		if (0 == fzTemp.n_faceShape)
		{
			for (int i = 0; i < fzTemp.n_faceNum; i++)
			{
				int nFaceType;
				Face faceTemp;

				inSection.ss_details >> std::hex >> nFaceType;
				faceTemp.v_nodeID.resize(nFaceType);
				if (2 == nFaceType)
				{
					faceTemp.ft_faceType = Face::ftLinear;
					inSection.ss_details >> std::hex >> faceTemp.v_nodeID[0] >> faceTemp.v_nodeID[1];
					inSection.ss_details >> std::hex >> faceTemp.n_owner >> faceTemp.n_neighbor;
				}
					
				else if (3 == nFaceType)
				{
                    faceTemp.ft_faceType = Face::ftTrangular;
					inSection.ss_details >> std::hex >> faceTemp.v_nodeID[0] >> faceTemp.v_nodeID[1] >> faceTemp.v_nodeID[2];
					inSection.ss_details >> std::hex >> faceTemp.n_owner >> faceTemp.n_neighbor;
				}
				else if (4 == nFaceType)
				{
                    faceTemp.ft_faceType = Face::ftQuadrilateral;
					inSection.ss_details >> std::hex >> faceTemp.v_nodeID[0] >> faceTemp.v_nodeID[1] >> faceTemp.v_nodeID[2] >> faceTemp.v_nodeID[3];
					inSection.ss_details >> std::hex >> faceTemp.n_owner >> faceTemp.n_neighbor;
				}
				else
				{
					faceTemp.ft_faceType = Face::ftPolygon;
					for (int i = 0; i < nFaceType; i++)
					{
						inSection.ss_details >> std::hex >> faceTemp.v_nodeID[i];
					}
					inSection.ss_details >> std::hex >> faceTemp.n_owner >> faceTemp.n_neighbor;
				}

				if (0 == faceTemp.n_owner)
				{
					faceTemp.n_owner = faceTemp.n_neighbor;
					faceTemp.n_neighbor = 0;
                    std::vector<int> temp = faceTemp.v_nodeID;
					faceTemp.v_nodeID.assign(temp.rbegin(), temp.rend());
				}

				this->v_face.push_back(faceTemp);
			}
		}
		else if (2 == fzTemp.n_faceShape)
		{
			for (int i = 0; i < fzTemp.n_faceNum; i++)
			{
				Face faceTemp;
				faceTemp.v_nodeID.resize(2);

                faceTemp.ft_faceType = Face::ftLinear;
				inSection.ss_details >> std::hex >> faceTemp.v_nodeID[0] >> faceTemp.v_nodeID[1];
				inSection.ss_details >> std::hex >> faceTemp.n_owner >> faceTemp.n_neighbor;

				if (0 == faceTemp.n_owner)
				{
					faceTemp.n_owner = faceTemp.n_neighbor;
					faceTemp.n_neighbor = 0;
                    std::vector<int> temp = faceTemp.v_nodeID;
					faceTemp.v_nodeID.assign(temp.rbegin(), temp.rend());
				}

				this->v_face.push_back(faceTemp);
			}

		}
		else if (3 == fzTemp.n_faceShape)
		{
			for (int i = 0; i < fzTemp.n_faceNum; i++)
			{
				Face faceTemp;
				faceTemp.v_nodeID.resize(3);

                faceTemp.ft_faceType = Face::ftTrangular;
				inSection.ss_details >> std::hex >> faceTemp.v_nodeID[0] >> faceTemp.v_nodeID[1] >> faceTemp.v_nodeID[2];
				inSection.ss_details >> std::hex >> faceTemp.n_owner >> faceTemp.n_neighbor;

				if (0 == faceTemp.n_owner)
				{
					faceTemp.n_owner = faceTemp.n_neighbor;
					faceTemp.n_neighbor = 0;
                    std::vector<int> temp = faceTemp.v_nodeID;
					faceTemp.v_nodeID.assign(temp.rbegin(), temp.rend());
				}

				this->v_face.push_back(faceTemp);
			}
		}
		else if (4 == fzTemp.n_faceShape)
		{
			for (int i = 0; i < fzTemp.n_faceNum; i++)
			{
				Face faceTemp;
				faceTemp.v_nodeID.resize(4);

                faceTemp.ft_faceType = Face::ftQuadrilateral;
				inSection.ss_details >> std::hex >> faceTemp.v_nodeID[0] >> faceTemp.v_nodeID[1] >> faceTemp.v_nodeID[2] >> faceTemp.v_nodeID[3];
				inSection.ss_details >> std::hex >> faceTemp.n_owner >> faceTemp.n_neighbor;

				if (0 == faceTemp.n_owner)
				{
					faceTemp.n_owner = faceTemp.n_neighbor;
					faceTemp.n_neighbor = 0;
                    std::vector<int> temp = faceTemp.v_nodeID;
					faceTemp.v_nodeID.assign(temp.rbegin(), temp.rend());
				}

				this->v_face.push_back(faceTemp);
			}
		}
		else if (5 == fzTemp.n_faceShape)
		{
			for (int i = 0; i < fzTemp.n_faceNum; i++)
			{
				int nFaceType;
				Face faceTemp;

				inSection.ss_details >> std::hex >> nFaceType;
				faceTemp.v_nodeID.resize(nFaceType);

				faceTemp.ft_faceType = Face::ftPolygon;
				for (int i = 0; i < nFaceType; i++)
				{
					inSection.ss_details >> std::hex >> faceTemp.v_nodeID[i];
				}
				inSection.ss_details >> std::hex >> faceTemp.n_owner >> faceTemp.n_neighbor;

				if (0 == faceTemp.n_owner)
				{
					faceTemp.n_owner = faceTemp.n_neighbor;
					faceTemp.n_neighbor = 0;
                    std::vector<int> temp = faceTemp.v_nodeID;
					faceTemp.v_nodeID.assign(temp.rbegin(), temp.rend());
				}

				this->v_face.push_back(faceTemp);
			}
		}
	}

	/*
	//swap owner and neighbor at boundary (For Gambit)
	for (int i = 0; i < (int)this->v_face.size(); i++)
	{
		int nOwner = this->v_face[i].n_owner;
		if (0 == nOwner)
		{
			this->v_face[i].n_owner = this->v_face[i].n_neighbor;
			this->v_face[i].n_neighbor = 0;
			vector<int> temp = this->v_face[i].v_nodeID;
			this->v_face[i].v_nodeID.assign(temp.rbegin(), temp.rend());
		}
	}
	*/
}
void FluentMeshBlock::ReadFace(FILE *F_FilePointer)
{
    std::string s_leftChar;
	while (s_leftChar != "(")							//erase space before "("
	{
		s_leftChar = fgetc(F_FilePointer);
	}
	ReadNew(F_FilePointer, 13);
	if (v_parameters.size() != 5)
	{
        FatalError(std::string("Incorrect parameter number detected in face section."));
	}
	
	if (0 == v_parameters[0])
	{
		std::cout << LogLevel(LogLevel::mlOK, "Reading face declaration") << std::endl;
		int beginNum = v_parameters[1];
		int endNum = v_parameters[2];
		this->n_faceNum = endNum - beginNum + 1;
	}
	else
	{
		s_leftChar = fgetc(F_FilePointer);
		while (s_leftChar != "(")							//erase space before "("
		{
			s_leftChar = fgetc(F_FilePointer);
		}
		s_leftChar = fgetc(F_FilePointer);
		std::cout << LogLevel(LogLevel::mlOK, "Reading face topologies") << std::endl;
		MshFaceZone fzTemp;

		//(zone-id first-n_index last-n_index type element-type)
		//type : mixed,linear,triangular,quadrilateral
		fzTemp.n_zoneID = v_parameters[0];
		fzTemp.n_begID = v_parameters[1];
		fzTemp.n_endID = v_parameters[2];
		fzTemp.n_BCType = v_parameters[3];
        fzTemp.n_faceShape = v_parameters[4];

		fzTemp.n_faceNum = fzTemp.n_endID - fzTemp.n_begID + 1;

		v_FaceZoneInfo.push_back(fzTemp);

		//Read face
		if (0 == fzTemp.n_faceShape)
		{
			for (int i = 0; i < fzTemp.n_faceNum; i++)
			{
				int nFaceType;
				Face faceTemp;

				nFaceType = GetData(F_FilePointer, s_leftChar);
				faceTemp.v_nodeID.resize(nFaceType);
				if (2 == nFaceType)
				{
					faceTemp.ft_faceType = Face::ftLinear;
					faceTemp.v_nodeID[0] = GetData(F_FilePointer, s_leftChar);
					faceTemp.v_nodeID[1] = GetData(F_FilePointer, s_leftChar);
					faceTemp.n_owner = GetData(F_FilePointer, s_leftChar);
					faceTemp.n_neighbor = GetData(F_FilePointer, s_leftChar);
				}

				else if (3 == nFaceType)
				{
					faceTemp.ft_faceType = Face::ftTrangular;
					faceTemp.v_nodeID[0] = GetData(F_FilePointer, s_leftChar);
					faceTemp.v_nodeID[1] = GetData(F_FilePointer, s_leftChar);
					faceTemp.v_nodeID[2] = GetData(F_FilePointer, s_leftChar);
					faceTemp.n_owner = GetData(F_FilePointer, s_leftChar);
					faceTemp.n_neighbor = GetData(F_FilePointer, s_leftChar);
				}
				else if (4 == nFaceType)
				{
					faceTemp.ft_faceType = Face::ftQuadrilateral;
					faceTemp.v_nodeID[0] = GetData(F_FilePointer, s_leftChar);
					faceTemp.v_nodeID[1] = GetData(F_FilePointer, s_leftChar); 
					faceTemp.v_nodeID[2] = GetData(F_FilePointer, s_leftChar); 
					faceTemp.v_nodeID[3] = GetData(F_FilePointer, s_leftChar);
					faceTemp.n_owner = GetData(F_FilePointer, s_leftChar); 
					faceTemp.n_neighbor = GetData(F_FilePointer, s_leftChar);;
				}
				else
				{
					faceTemp.ft_faceType = Face::ftPolygon;
					for (int i = 0; i < nFaceType; i++)
					{
						faceTemp.v_nodeID[i] = GetData(F_FilePointer, s_leftChar);
					}
					faceTemp.n_owner = GetData(F_FilePointer, s_leftChar); 
					faceTemp.n_neighbor = GetData(F_FilePointer, s_leftChar);
				}

				if (0 == faceTemp.n_owner)
				{
					faceTemp.n_owner = faceTemp.n_neighbor;
					faceTemp.n_neighbor = 0;
                    std::vector<int> temp = faceTemp.v_nodeID;
					faceTemp.v_nodeID.assign(temp.rbegin(), temp.rend());
                }

				this->v_face.push_back(faceTemp);
			}
		}
		else if (2 == fzTemp.n_faceShape)
		{
			for (int i = 0; i < fzTemp.n_faceNum; i++)
			{
				Face faceTemp;
				faceTemp.v_nodeID.resize(2);

				faceTemp.ft_faceType = Face::ftLinear;
                faceTemp.v_nodeID[0] = GetData(F_FilePointer, s_leftChar);
				faceTemp.v_nodeID[1] = GetData(F_FilePointer, s_leftChar);
				faceTemp.n_owner = GetData(F_FilePointer, s_leftChar); 
				faceTemp.n_neighbor = GetData(F_FilePointer, s_leftChar);
				if (0 == faceTemp.n_owner)
				{
					faceTemp.n_owner = faceTemp.n_neighbor;
					faceTemp.n_neighbor = 0;
                    std::vector<int> temp = faceTemp.v_nodeID;
					faceTemp.v_nodeID.assign(temp.rbegin(), temp.rend());
				}

				this->v_face.push_back(faceTemp);
			}

		}
		else if (3 == fzTemp.n_faceShape)
		{
			for (int i = 0; i < fzTemp.n_faceNum; i++)
			{
				Face faceTemp;
				faceTemp.v_nodeID.resize(3);

				faceTemp.ft_faceType = Face::ftTrangular;
				faceTemp.v_nodeID[0] = GetData(F_FilePointer, s_leftChar); 
				faceTemp.v_nodeID[1] = GetData(F_FilePointer, s_leftChar); 
				faceTemp.v_nodeID[2] = GetData(F_FilePointer, s_leftChar);
				faceTemp.n_owner = GetData(F_FilePointer, s_leftChar); 
				faceTemp.n_neighbor = GetData(F_FilePointer, s_leftChar);

				if (0 == faceTemp.n_owner)
				{
					faceTemp.n_owner = faceTemp.n_neighbor;
					faceTemp.n_neighbor = 0;
                    std::vector<int> temp = faceTemp.v_nodeID;
					faceTemp.v_nodeID.assign(temp.rbegin(), temp.rend());
				}

				this->v_face.push_back(faceTemp);
			}
		}
		else if (4 == fzTemp.n_faceShape)
		{
			for (int i = 0; i < fzTemp.n_faceNum; i++)
			{
				Face faceTemp;
				faceTemp.v_nodeID.resize(4);

				faceTemp.ft_faceType = Face::ftQuadrilateral;
				faceTemp.v_nodeID[0] = GetData(F_FilePointer, s_leftChar); 
				faceTemp.v_nodeID[1] = GetData(F_FilePointer, s_leftChar); 
				faceTemp.v_nodeID[2] = GetData(F_FilePointer, s_leftChar); 
				faceTemp.v_nodeID[3] = GetData(F_FilePointer, s_leftChar);
				faceTemp.n_owner = GetData(F_FilePointer, s_leftChar); 
				faceTemp.n_neighbor = GetData(F_FilePointer, s_leftChar);

				if (0 == faceTemp.n_owner)
				{
					faceTemp.n_owner = faceTemp.n_neighbor;
					faceTemp.n_neighbor = 0;
                    std::vector<int> temp = faceTemp.v_nodeID;
					faceTemp.v_nodeID.assign(temp.rbegin(), temp.rend());
				}

				this->v_face.push_back(faceTemp);
			}
		}
		else if (5 == fzTemp.n_faceShape)
		{
			for (int i = 0; i < fzTemp.n_faceNum; i++)
			{
				int nFaceType;
				Face faceTemp;

				nFaceType = GetData(F_FilePointer, s_leftChar);
				faceTemp.v_nodeID.resize(nFaceType);

				faceTemp.ft_faceType = Face::ftPolygon;
				for (int i = 0; i < nFaceType; i++)
				{
					faceTemp.v_nodeID[i] = GetData(F_FilePointer, s_leftChar);
				}
				faceTemp.n_owner = GetData(F_FilePointer, s_leftChar); 
				faceTemp.n_neighbor = GetData(F_FilePointer, s_leftChar);

				if (0 == faceTemp.n_owner)
				{
					faceTemp.n_owner = faceTemp.n_neighbor;
					faceTemp.n_neighbor = 0;
                    std::vector<int> temp = faceTemp.v_nodeID;
					faceTemp.v_nodeID.assign(temp.rbegin(), temp.rend());
				}

				this->v_face.push_back(faceTemp);
			}
		}
	}
}

void FluentMeshBlock::ReadZone(FluentMeshFileSection& inSection)
{
	int zoneNumber;
	inSection.ss_details >> std::dec >> zoneNumber;
    std::string zoneName;
	while (getline(inSection.ss_details, zoneName, ' '))
	{
	}

	for (int i = 0; i < (int)v_ElementZoneInfo.size(); i++)
	{
		if (zoneNumber == v_ElementZoneInfo[i].n_zoneID)
		{
			v_ElementZoneInfo[i].st_name = zoneName;
		}
	}

	for (int i = 0; i < (int)v_FaceZoneInfo.size(); i++)
	{
		if (zoneNumber == v_FaceZoneInfo[i].n_zoneID)
		{
			v_FaceZoneInfo[i].st_name = zoneName;
		}
	}
}
void FluentMeshBlock::ReadZone(FILE *F_FilePointer)
{
	int zoneNumber;
    std::string c_Temp;
	c_Temp = fgetc(F_FilePointer);
	while (c_Temp !="(")
	{
		c_Temp = fgetc(F_FilePointer);
	}
    std::string  zoneName_Second;
	fscanf(F_FilePointer, "%d", &zoneNumber);				
	c_Temp = fgetc(F_FilePointer);
	while (c_Temp != ")")
	{
		zoneName_Second = "";
        while (c_Temp == "\t" || c_Temp == "\n" || c_Temp == " ")
		{
			c_Temp = fgetc(F_FilePointer);
		}
        while (c_Temp != "\t" && c_Temp != "\n" && c_Temp != " "&&c_Temp != ")")
		{ 
			zoneName_Second = zoneName_Second + c_Temp;
			c_Temp = fgetc(F_FilePointer);
		}

	}
//	cout << "\t" << zoneName_Second << endl;
//	system("pause");

	for (int i = 0; i < (int)v_ElementZoneInfo.size(); i++)
	{
		if (zoneNumber == v_ElementZoneInfo[i].n_zoneID)
		{
			v_ElementZoneInfo[i].st_name = zoneName_Second;
		}
	}
	for (int i = 0; i < (int)v_FaceZoneInfo.size(); i++)
	{
		if (zoneNumber == v_FaceZoneInfo[i].n_zoneID)
		{
			v_FaceZoneInfo[i].st_name = zoneName_Second;
		}
	}
	if(c_Temp == ")")
	{
		c_Temp = fgetc(F_FilePointer);
		while (c_Temp != "(")
		{
			c_Temp = fgetc(F_FilePointer);
		}
		while (c_Temp != ")")
		{
			c_Temp = fgetc(F_FilePointer);
		}
		c_Temp = fgetc(F_FilePointer);
		while (c_Temp != ")")
		{
			c_Temp = fgetc(F_FilePointer);
		}
	}
	else
	{
		c_Temp = fgetc(F_FilePointer);
		while (c_Temp != ")")
		{
			c_Temp = fgetc(F_FilePointer);
		}
		while (c_Temp != "(")
		{
			c_Temp = fgetc(F_FilePointer);
		}
		while (c_Temp != ")")
		{
			c_Temp = fgetc(F_FilePointer);
		}
		c_Temp = fgetc(F_FilePointer);
		while (c_Temp != ")")
		{
			c_Temp = fgetc(F_FilePointer);
		}
	}
}

void FluentMeshBlock::EstablishZones()
{
	for (int i = 0; i < (int)this->v_ElementZoneInfo.size(); i++)
	{
		ElementZone ezTemp;
		ezTemp.name = v_ElementZoneInfo[i].st_name;
		for (int j = v_ElementZoneInfo[i].n_begID; j <= v_ElementZoneInfo[i].n_endID; j++)
		{
			ezTemp.v_elementID.push_back(j);
			v_elem[j].n_ElementZoneID = i;
		}
		v_elementZone.push_back(ezTemp);
	}

	for (int i = 0; i < (int)this->v_FaceZoneInfo.size(); i++)
	{
		FaceZone fzTemp;
		fzTemp.name = v_FaceZoneInfo[i].st_name;
		fzTemp.bc_Type = (FaceZone::BCType)v_FaceZoneInfo[i].n_BCType;
		for (int j = v_FaceZoneInfo[i].n_begID; j <= v_FaceZoneInfo[i].n_endID; j++)
		{
			fzTemp.v_faceID.push_back(j);
			v_face[j].n_FaceZoneID = i;
		}
		//insert into the facezone list
		v_faceZone.push_back(fzTemp);
	}
	//Determing the face zone type: interior? external? internal? or mixed?
	for (int i = 0; i < (int)this->v_FaceZoneInfo.size(); i++)
	{
		if (FaceZone::bcInterior == v_faceZone[i].bc_Type)
		{
			v_faceZone[i].faceZoneType = FaceZone::fztInterior;
		}
		else
		{
			int FirstFaceID = v_faceZone[i].v_faceID[0];
			if (-1 != v_face[FirstFaceID].n_owner && -1 == v_face[FirstFaceID].n_neighbor)
			{
				v_faceZone[i].faceZoneType = FaceZone::fztExternal;
				for (int j = 1; j < (int)v_faceZone[i].v_faceID.size(); j++)
				{
					int FaceID = v_faceZone[i].v_faceID[j];
					std::pair<int, int> oAndN = SearchOwnerNeighbor(FaceID);
					if (-1 == oAndN.first || -1 != oAndN.second)
					{
						v_faceZone[i].faceZoneType = FaceZone::fztMixed;
                        std::string info = "face zone " + v_faceZone[i].name + " is found to be mixed";
						WarningPause(info);
						break;
					}
				}
			}
			else
			{
				v_faceZone[i].faceZoneType = FaceZone::fztInternal;
				for (int j = 1; j < (int)v_faceZone[i].v_faceID.size(); j++)
				{
					int FaceID = v_faceZone[i].v_faceID[j];
					std::pair<int, int> oAndN = SearchOwnerNeighbor(FaceID);
					if (-1 == oAndN.first || -1 == oAndN.second)
					{
                        std::string info = "face zone " + v_faceZone[i].name + " is found to be mixed";
						WarningPause(info);
						break;
					}
				}
			}
		}
	}
}

void FluentMeshBlock::RefineArrayIndex()
{
	for (int i = 0; i < (int)this->v_FaceZoneInfo.size(); i++)
	{
		this->v_FaceZoneInfo[i].n_begID -= 1;
		this->v_FaceZoneInfo[i].n_endID -= 1;
	}

	for (int i = 0; i < (int)this->v_ElementZoneInfo.size(); i++)
	{
		this->v_ElementZoneInfo[i].n_begID -= 1;
		this->v_ElementZoneInfo[i].n_endID -= 1;
	}

	//Refine Node label in Face
	for (int i = 0; i < (int)this->v_face.size(); i++)
	{
		for (int j = 0; j < (int)this->v_face[i].v_nodeID.size(); j++)
		{
			this->v_face[i].v_nodeID[j] -= 1;
		}

		this->v_face[i].n_owner -= 1;
		this->v_face[i].n_neighbor -= 1;
	}

	//Put all face ID in corresponding Element
	for (int i = 0; i < (int)this->v_face.size(); i++)
	{
		int nOwn = this->v_face[i].n_owner;
		int nNei = this->v_face[i].n_neighbor;

		if (nNei != -1)
		{
			this->v_elem[nNei].v_faceID.push_back(i);
		}

		this->v_elem[nOwn].v_faceID.push_back(i);
	}
}

void FluentMeshBlock::PopulateCellNodes()
{

	for (int i = 0; i < (int)this->v_elem.size(); i++)
	{
        switch (this->v_elem[i].est_shapeType)
		{
        case Element::estTriangular:
			this->PopulateTriangleCell(i);
			break;

        case Element::estTetrahedral:
			this->PopulateTetraCell(i);
			break;

        case Element::estQuadrilateral:
			this->PopulateQuadCell(i);
			break;

        case Element::estHexahedral:
			this->PopulateHexahedronCell(i);
			break;

        case Element::estPyramid:
			this->PopulatePyramidCell(i);
			break;

        case Element::estWedge:
			this->PopulateWedgeCell(i);
			break;

        case Element::estPolyhedron:
			this->PopulatePolyhedronCell(i);
			break;

        case Element::estNoType:
            WarningContinue(std::string("Element ElemShapeType is \"estNoType\"."));
			break;

        default:
            break;

		}
    }
}

void FluentMeshBlock::PopulateTriangleCell(int i)
{
	if (this->v_elem[i].v_faceID.size() != 3)
	{
        std::stringstream info;
		info << "Current Triangle Cell ID is:" << i + 1 << std::endl;
		info << "Current Triangle Cell Face Num is:" << this->v_elem[i].v_faceID.size() << std::endl;
		FatalError(info.str());
	}

	this->v_elem[i].v_nodeID.resize(3);
	int nFstFace = this->v_elem[i].v_faceID[0];
	int nSecFace = this->v_elem[i].v_faceID[1];

	if (this->v_face[nFstFace].n_owner == i)
	{
		this->v_elem[i].v_nodeID[0] = this->v_face[nFstFace].v_nodeID[0];
		this->v_elem[i].v_nodeID[1] = this->v_face[nFstFace].v_nodeID[1];
	}
	else
	{
		this->v_elem[i].v_nodeID[1] = this->v_face[nFstFace].v_nodeID[0];
		this->v_elem[i].v_nodeID[0] = this->v_face[nFstFace].v_nodeID[1];
	}

	if (this->v_face[nSecFace].v_nodeID[0] != this->v_elem[i].v_nodeID[0] &&
		this->v_face[nSecFace].v_nodeID[0] != this->v_elem[i].v_nodeID[1])
	{
		this->v_elem[i].v_nodeID[2] = this->v_face[nSecFace].v_nodeID[0];
	}
	else
	{
		this->v_elem[i].v_nodeID[2] = this->v_face[nSecFace].v_nodeID[1];
	}
}

void FluentMeshBlock::PopulateTetraCell(int i)
{
	if (this->v_elem[i].v_faceID.size() != 4)
	{
        std::stringstream info;
		info << "Current Tetrahedral Cell ID is:" << i + 1 << std::endl;
		info << "Current Tetrahedral Cell Face Num is:" << this->v_elem[i].v_faceID.size() << std::endl;
		FatalError(info.str());
	}

	this->v_elem[i].v_nodeID.resize(4);
	int nFstFace = this->v_elem[i].v_faceID[0];
	int nSecFace = this->v_elem[i].v_faceID[1];

	if (this->v_face[nFstFace].n_owner == i)
	{
		this->v_elem[i].v_nodeID[0] = this->v_face[nFstFace].v_nodeID[0];
		this->v_elem[i].v_nodeID[1] = this->v_face[nFstFace].v_nodeID[1];
		this->v_elem[i].v_nodeID[2] = this->v_face[nFstFace].v_nodeID[2];
	}
	else
	{
		this->v_elem[i].v_nodeID[2] = this->v_face[nFstFace].v_nodeID[0];
		this->v_elem[i].v_nodeID[1] = this->v_face[nFstFace].v_nodeID[1];
		this->v_elem[i].v_nodeID[0] = this->v_face[nFstFace].v_nodeID[2];
	}

	if (this->v_face[nSecFace].v_nodeID[0] != this->v_elem[i].v_nodeID[0] &&
		this->v_face[nSecFace].v_nodeID[0] != this->v_elem[i].v_nodeID[1] &&
		this->v_face[nSecFace].v_nodeID[0] != this->v_elem[i].v_nodeID[2])
	{
		this->v_elem[i].v_nodeID[3] = this->v_face[nSecFace].v_nodeID[0];
	}
	else if (this->v_face[nSecFace].v_nodeID[1] != this->v_elem[i].v_nodeID[0] &&
		this->v_face[nSecFace].v_nodeID[1] != this->v_elem[i].v_nodeID[1] &&
		this->v_face[nSecFace].v_nodeID[1] != this->v_elem[i].v_nodeID[2])
	{
		this->v_elem[i].v_nodeID[3] = this->v_face[nSecFace].v_nodeID[1];
	}
	else
	{
		this->v_elem[i].v_nodeID[3] = this->v_face[nSecFace].v_nodeID[2];
	}
}

void FluentMeshBlock::PopulateQuadCell(int i)
{
	if (this->v_elem[i].v_faceID.size() != 4)
	{
        std::stringstream info;
		info << "Current Quadrilateral Cell ID is:" << i + 1 << std::endl;
		info << "Current Quadrilateral Cell Face Num is:" << this->v_elem[i].v_faceID.size() << std::endl;
		FatalError(info.str());
	}

	this->v_elem[i].v_nodeID.resize(4);
	int nFstFace = this->v_elem[i].v_faceID[0];
	int nSecFace = this->v_elem[i].v_faceID[1];
	int n3rdFace = this->v_elem[i].v_faceID[2];
	int n4thFace = this->v_elem[i].v_faceID[3];

	if (this->v_face[nFstFace].n_owner == i)
	{
		this->v_elem[i].v_nodeID[0] = this->v_face[nFstFace].v_nodeID[0];
		this->v_elem[i].v_nodeID[1] = this->v_face[nFstFace].v_nodeID[1];
	}
	else
	{
		this->v_elem[i].v_nodeID[1] = this->v_face[nFstFace].v_nodeID[0];
		this->v_elem[i].v_nodeID[0] = this->v_face[nFstFace].v_nodeID[1];
	}

	if ((this->v_face[nSecFace].v_nodeID[0] != this->v_elem[i].v_nodeID[0] &&
		this->v_face[nSecFace].v_nodeID[0] != this->v_elem[i].v_nodeID[1]) &&
		(this->v_face[nSecFace].v_nodeID[1] != this->v_elem[i].v_nodeID[0] &&
		this->v_face[nSecFace].v_nodeID[1] != this->v_elem[i].v_nodeID[1]))
	{
		if (this->v_face[nSecFace].n_owner == i)
		{
			this->v_elem[i].v_nodeID[2] = this->v_face[nSecFace].v_nodeID[0];
			this->v_elem[i].v_nodeID[3] = this->v_face[nSecFace].v_nodeID[1];
		}
		else
		{
			this->v_elem[i].v_nodeID[3] = this->v_face[nSecFace].v_nodeID[0];
			this->v_elem[i].v_nodeID[2] = this->v_face[nSecFace].v_nodeID[1];
		}
	}
	else if ((this->v_face[n3rdFace].v_nodeID[0] != this->v_elem[i].v_nodeID[0] &&
		this->v_face[n3rdFace].v_nodeID[0] != this->v_elem[i].v_nodeID[1]) &&
		(this->v_face[n3rdFace].v_nodeID[1] != this->v_elem[i].v_nodeID[0] &&
		this->v_face[n3rdFace].v_nodeID[1] != this->v_elem[i].v_nodeID[1]))
	{
		if (this->v_face[n3rdFace].n_owner == i)
		{
			this->v_elem[i].v_nodeID[2] = this->v_face[n3rdFace].v_nodeID[0];
			this->v_elem[i].v_nodeID[3] = this->v_face[n3rdFace].v_nodeID[1];
		}
		else
		{
			this->v_elem[i].v_nodeID[3] = this->v_face[n3rdFace].v_nodeID[0];
			this->v_elem[i].v_nodeID[2] = this->v_face[n3rdFace].v_nodeID[1];
		}
	}
	else
	{
		if (this->v_face[n4thFace].n_owner == i)
		{
			this->v_elem[i].v_nodeID[2] = this->v_face[n4thFace].v_nodeID[0];
			this->v_elem[i].v_nodeID[3] = this->v_face[n4thFace].v_nodeID[1];
		}
		else
		{
			this->v_elem[i].v_nodeID[3] = this->v_face[n4thFace].v_nodeID[0];
			this->v_elem[i].v_nodeID[2] = this->v_face[n4thFace].v_nodeID[1];
		}
	}
}

void FluentMeshBlock::PopulateHexahedronCell(int i)
{
	if (this->v_elem[i].v_faceID.size() != 6)
	{
        std::stringstream info;
		info << "Current Hexahedral Cell ID is:" << i + 1 << std::endl;
		info << "Current Hexahedral Cell Face Num is:" << this->v_elem[i].v_faceID.size() << std::endl;
		FatalError(info.str());
	}

	this->v_elem[i].v_nodeID.resize(8);
	int nFstFace = this->v_elem[i].v_faceID[0];

	if (this->v_face[nFstFace].n_owner == i)
	{
		for (int j = 0; j < 4; j++)
		{
			this->v_elem[i].v_nodeID[j] = this->v_face[nFstFace].v_nodeID[j];
		}
	}
	else
	{
		for (int j = 3; j >= 0; j--)
		{
			this->v_elem[i].v_nodeID[3 - j] = this->v_face[nFstFace].v_nodeID[j];
		}
	}

	//Look for opposite face of hexahedron
	for (int j = 1; j < 6; j++)
	{
		int flag = 0;
		int nCurFace = this->v_elem[i].v_faceID[j];

		for (int k = 0; k < 4; k++)
		{
			if ((this->v_elem[i].v_nodeID[0] == this->v_face[nCurFace].v_nodeID[k]) ||
				(this->v_elem[i].v_nodeID[1] == this->v_face[nCurFace].v_nodeID[k]) ||
				(this->v_elem[i].v_nodeID[2] == this->v_face[nCurFace].v_nodeID[k]) ||
				(this->v_elem[i].v_nodeID[3] == this->v_face[nCurFace].v_nodeID[k]))
			{
				flag = 1;
			}
		}

		if (flag == 0)
		{
			if (this->v_face[nCurFace].n_owner == i)
			{
				for (int k = 7; k >= 4; k--)
				{
					this->v_elem[i].v_nodeID[k] = this->v_face[nCurFace].v_nodeID[7 - k];
				}
			}
			else
			{
				for (int k = 4; k < 8; k++)
				{
					this->v_elem[i].v_nodeID[k] = this->v_face[nCurFace].v_nodeID[k - 4];
				}
			}
		}
	}

	//  Find the face with points 0 and 1 in them.
	int f01[4] = { -1, -1, -1, -1 };
	for (int j = 1; j < 6; j++)
	{
		int flag0 = 0;
		int flag1 = 0;
		int nCurFace = this->v_elem[i].v_faceID[j];

		for (int k = 0; k < 4; k++)
		{
			if (this->v_elem[i].v_nodeID[0] == this->v_face[nCurFace].v_nodeID[k])
			{
				flag0 = 1;
			}

			if (this->v_elem[i].v_nodeID[1] == this->v_face[nCurFace].v_nodeID[k])
			{
				flag1 = 1;
			}
		}

		if ((flag0 == 1) && (flag1 == 1))
		{
			if (this->v_face[nCurFace].n_owner == i)
			{
				for (int k = 0; k<4; k++)
				{
					f01[k] = this->v_face[nCurFace].v_nodeID[k];
				}
			}
			else
			{
				for (int k = 3; k >= 0; k--)
				{
					f01[k] = this->v_face[nCurFace].v_nodeID[k];
				}
			}
		}
	}

	//  Find the face with points 0 and 3 in them.
	int f03[4] = { -1, -1, -1, -1 };
	for (int j = 1; j < 6; j++)
	{
		int flag0 = 0;
		int flag1 = 0;
		int nCurFace = this->v_elem[i].v_faceID[j];

		for (int k = 0; k < 4; k++)
		{
			if (this->v_elem[i].v_nodeID[0] == this->v_face[nCurFace].v_nodeID[k])
			{
				flag0 = 1;
			}
			if (this->v_elem[i].v_nodeID[3] == this->v_face[nCurFace].v_nodeID[k])
			{
				flag1 = 1;
			}
		}

		if ((flag0 == 1) && (flag1 == 1))
		{
			if (this->v_face[nCurFace].n_owner == i)
			{
				for (int k = 0; k<4; k++)
				{
					f03[k] = this->v_face[nCurFace].v_nodeID[k];
				}
			}
			else
			{
				for (int k = 3; k >= 0; k--)
				{
					f03[k] = this->v_face[nCurFace].v_nodeID[k];
				}
			}
		}
	}

	// What point is in f01 and f03 besides 0 ... this is point 4
	int p4 = 0;
	for (int k = 0; k < 4; k++)
	{
		if (f01[k] != this->v_elem[i].v_nodeID[0])
		{
			for (int n = 0; n < 4; n++)
			{
				if (f01[k] == f03[n])
				{
					p4 = f01[k];
				}
			}
		}
	}

	// Since we know point 4 now we check to see if points
	//  4, 5, 6, and 7 are in the correct positions.
	int t[8];
	t[4] = this->v_elem[i].v_nodeID[4];
	t[5] = this->v_elem[i].v_nodeID[5];
	t[6] = this->v_elem[i].v_nodeID[6];
	t[7] = this->v_elem[i].v_nodeID[7];
	if (p4 == this->v_elem[i].v_nodeID[5])
	{
		this->v_elem[i].v_nodeID[4] = t[5];
		this->v_elem[i].v_nodeID[5] = t[6];
		this->v_elem[i].v_nodeID[6] = t[7];
		this->v_elem[i].v_nodeID[7] = t[4];
	}
	else if (p4 == this->v_elem[i].v_nodeID[6])
	{
		this->v_elem[i].v_nodeID[4] = t[6];
		this->v_elem[i].v_nodeID[5] = t[7];
		this->v_elem[i].v_nodeID[6] = t[4];
		this->v_elem[i].v_nodeID[7] = t[5];
	}
	else if (p4 == this->v_elem[i].v_nodeID[7])
	{
		this->v_elem[i].v_nodeID[4] = t[7];
		this->v_elem[i].v_nodeID[5] = t[4];
		this->v_elem[i].v_nodeID[6] = t[5];
		this->v_elem[i].v_nodeID[7] = t[6];
	}
	// else point 4 was lined up so everything was correct.
}

void FluentMeshBlock::PopulatePyramidCell(int i)
{
	if (this->v_elem[i].v_faceID.size() != 5)
	{
        std::stringstream info;
		info << "Current Pyramid Cell ID is:" << i + 1 << std::endl;
		info << "Current Pyramid Cell Face Num is:" << this->v_elem[i].v_faceID.size() << std::endl;
		FatalError(info.str());
	}

	this->v_elem[i].v_nodeID.resize(5);

	//The quad face will be the base of the pyramid
	for (int j = 0; j < (int)this->v_elem[i].v_faceID.size(); j++)
	{
		int nCurFace = this->v_elem[i].v_faceID[j];

		if (this->v_face[nCurFace].v_nodeID.size() == 4)
		{
			if (this->v_face[nCurFace].n_owner == i)
			{
				for (int k = 0; k < 4; k++)
				{
					this->v_elem[i].v_nodeID[k] = this->v_face[nCurFace].v_nodeID[k];
				}
			}
			else
			{
				for (int k = 0; k < 4; k++)
				{
					this->v_elem[i].v_nodeID[3 - k] = this->v_face[nCurFace].v_nodeID[k];
				}
			}
		}
	}

	//Just need to find point 4
	for (int j = 0; j < (int)this->v_elem[i].v_faceID.size(); j++)
	{
		int nCurFace = this->v_elem[i].v_faceID[j];

		if (this->v_face[nCurFace].v_nodeID.size() == 3)
		{
			for (int k = 0; k < 3; k++)
			{
				if ((this->v_face[nCurFace].v_nodeID[k] != this->v_elem[i].v_nodeID[0]) &&
					(this->v_face[nCurFace].v_nodeID[k] != this->v_elem[i].v_nodeID[1]) &&
					(this->v_face[nCurFace].v_nodeID[k] != this->v_elem[i].v_nodeID[2]) &&
					(this->v_face[nCurFace].v_nodeID[k] != this->v_elem[i].v_nodeID[3]))
				{
					this->v_elem[i].v_nodeID[4] = this->v_face[nCurFace].v_nodeID[k];
				}
			}
		}
	}
}

void FluentMeshBlock::PopulateWedgeCell(int i)
{
	if (this->v_elem[i].v_faceID.size() != 5)
	{
        std::stringstream info;
		info << "Current Wedge Cell ID is:" << i + 1 << std::endl;
		info << "Current Wedge Cell Face Num is:" << this->v_elem[i].v_faceID.size() << std::endl;
		FatalError(info.str());
	}

	this->v_elem[i].v_nodeID.resize(6);

	//  Find the first triangle face and make it the base.
	int base = 0;
	int first = 0;
	for (int j = 0; j < (int)this->v_elem[i].v_faceID.size(); j++)
	{
		int nCurFace = this->v_elem[i].v_faceID[j];

		if ((this->v_face[nCurFace].v_nodeID.size() == 3) && (first == 0))
		{
			base = this->v_elem[i].v_faceID[j];
			first = 1;
		}
	}

	//  Find the second triangle face and make it the top.
	int top = 0;
	int second = 0;
	for (int j = 0; j < (int)this->v_elem[i].v_faceID.size(); j++)
	{
		int nCurFace = this->v_elem[i].v_faceID[j];

		if ((this->v_face[nCurFace].v_nodeID.size() == 3) && (second == 0) && (this->v_elem[i].v_faceID[j] != base))
		{
			top = this->v_elem[i].v_faceID[j];
			second = 1;
		}
	}

	// Load Base nodes into the nodes std::vector
	if (this->v_face[base].n_owner == i)
	{
		for (int j = 0; j < 3; j++)
		{
			this->v_elem[i].v_nodeID[j] = this->v_face[base].v_nodeID[j];
		}
	}
	else
	{
		for (int j = 2; j >= 0; j--)
		{
			this->v_elem[i].v_nodeID[2 - j] = this->v_face[base].v_nodeID[j];
		}
	}

	// Load Top nodes into the nodes std::vector
	if (this->v_face[top].n_neighbor == i)
	{
		for (int j = 3; j < 6; j++)
		{
			this->v_elem[i].v_nodeID[j] = this->v_face[top].v_nodeID[j - 3];
		}
	}
	else
	{
		for (int j = 3; j < 6; j++)
		{
			this->v_elem[i].v_nodeID[j] = this->v_face[top].v_nodeID[5 - j];
		}
	}

	//  Find the quad face with points 0 and 1 in them.
	int w01[4] = { -1, -1, -1, -1 };
	for (int j = 0; j < (int)this->v_elem[i].v_faceID.size(); j++)
	{
		int nCurFace = this->v_elem[i].v_faceID[j];

		if (nCurFace != base && nCurFace != top)
		{
			int wf0 = 0;
			int wf1 = 0;
			for (int k = 0; k < 4; k++)
			{
				if (this->v_elem[i].v_nodeID[0] == this->v_face[nCurFace].v_nodeID[k])
				{
					wf0 = 1;
				}

				if (this->v_elem[i].v_nodeID[1] == this->v_face[nCurFace].v_nodeID[k])
				{
					wf1 = 1;
				}

				if ((wf0 == 1) && (wf1 == 1))
				{
					for (int n = 0; n<4; n++)
					{
						w01[n] = this->v_face[nCurFace].v_nodeID[n];
					}
				}
			}
		}
	}

	//  Find the quad face with points 0 and 2 in them.
	int w02[4] = { -1, -1, -1, -1 };
	for (int j = 0; j < (int)this->v_elem[i].v_faceID.size(); j++)
	{
		int nCurFace = this->v_elem[i].v_faceID[j];

		if (nCurFace != base &&	nCurFace != top)
		{
			int wf0 = 0;
			int wf2 = 0;
			for (int k = 0; k < 4; k++)
			{
				if (this->v_elem[i].v_nodeID[0] == this->v_face[nCurFace].v_nodeID[k])
				{
					wf0 = 1;
				}
				if (this->v_elem[i].v_nodeID[2] == this->v_face[nCurFace].v_nodeID[k])
				{
					wf2 = 1;
				}
				if ((wf0 == 1) && (wf2 == 1))
				{
					for (int n = 0; n<4; n++)
					{
						w02[n] = this->v_face[nCurFace].v_nodeID[n];
					}
				}
			}
		}
	}

	// Point 3 is the point that is in both w01 and w02

	// What point is in f01 and f02 besides 0 ... this is point 3
	int p3 = 0;
	for (int k = 0; k < 4; k++)
	{
		if (w01[k] != this->v_elem[i].v_nodeID[0])
		{
			for (int n = 0; n < 4; n++)
			{
				if (w01[k] == w02[n])
				{
					p3 = w01[k];
				}
			}
		}
	}

	// Since we know point 3 now we check to see if points
	//  3, 4, and 5 are in the correct positions.
	int t[6];
	t[3] = this->v_elem[i].v_nodeID[3];
	t[4] = this->v_elem[i].v_nodeID[4];
	t[5] = this->v_elem[i].v_nodeID[5];
	if (p3 == this->v_elem[i].v_nodeID[4])
	{
		this->v_elem[i].v_nodeID[3] = t[4];
		this->v_elem[i].v_nodeID[4] = t[5];
		this->v_elem[i].v_nodeID[5] = t[3];
	}
	else if (p3 == this->v_elem[i].v_nodeID[5])
	{
		this->v_elem[i].v_nodeID[3] = t[5];
		this->v_elem[i].v_nodeID[4] = t[3];
		this->v_elem[i].v_nodeID[5] = t[4];
	}
	// else point 3 was lined up so everything was correct.

}

void FluentMeshBlock::PopulatePolyhedronCell(int i)
{
    if (this->v_elem[i].v_faceID.size() < 5)
    {
        std::stringstream info;
        info << "Current Polyhedron Cell ID is:" << i + 1 << std::endl;
        info << "Current Polyhedron Cell Face Num is:" << this->v_elem[i].v_faceID.size() << std::endl;
        FatalError(info.str());
    }

	for (int j = 0; j < (int)this->v_elem[i].v_faceID.size(); j++)
	{
		int CurFaceID(this->v_elem[i].v_faceID[j]);
		for (int k = 0; k < (int)this->v_face[CurFaceID].v_nodeID.size(); k++)
		{
			int flag(0);
			// Is the node already in the cell?
			for (int n = 0; n < (int)this->v_elem[i].v_nodeID.size(); n++)
			{
				if (this->v_elem[i].v_nodeID[n] ==
					this->v_face[CurFaceID].v_nodeID[k])
				{
					flag = 1;
				}
			}
			if (flag == 0)
			{
				//No match - insert node into cell.
				this->v_elem[i].v_nodeID.push_back(this->v_face[CurFaceID].v_nodeID[k]);
			}
		}
	}
}

void FluentMeshBlock::WriteTecplotMesh(const std::string& outMshFileName)
{
    std::ofstream	outFile(outMshFileName.c_str());

	outFile << "TITLE =\"" << "UnGrid Output" << "\"" << std::endl;

	if (this->md_meshDim == Mesh::md2D)
	{
		outFile << "VARIABLES = " << "\"x [m]\"," << "\"y [m]\"" << std::endl;

        if (this->est_shapeType == Element::estTriangular)
		{
			outFile << "ZONE T = " << "\"" << "Tri_mesh" << "\"," << " DATAPACKING = POINT, N = " << this->n_nodeNum << ", E = " << this->n_elemNum
				<< ", ZONETYPE = FETRIANGLE" << std::endl;
		}
        else if (this->est_shapeType == Element::estQuadrilateral)
		{
			outFile << "ZONE T = " << "\"" << "Quad_mesh" << "\"," << " DATAPACKING = POINT, N = " << this->n_nodeNum << ", E = " << this->n_elemNum
				<< ", ZONETYPE = FEQUADRILATERAL" << std::endl;
		}
        else if (this->est_shapeType == Element::estMixed)
		{
			outFile << "ZONE T = " << "\"" << "Mix_mesh" << "\"," << " DATAPACKING = POINT, N = " << this->n_nodeNum << ", E = " << this->n_elemNum
				<< ", ZONETYPE = FEQUADRILATERAL" << std::endl;
		}

		for (int i = 0; i < (int)this->n_nodeNum; i++)
		{
			outFile << this->v_node[i].x_ << " " << this->v_node[i].y_ << std::endl;
		}

		for (int i = 0; i < (int)this->n_elemNum; i++)
		{
            if (this->est_shapeType == Element::estTriangular)
			{
				outFile << this->v_elem[i].v_nodeID[0] + 1 << " " << this->v_elem[i].v_nodeID[1] + 1 << " " << this->v_elem[i].v_nodeID[2] + 1 << std::endl;
			}
            else if (this->est_shapeType == Element::estQuadrilateral)
			{
				outFile << this->v_elem[i].v_nodeID[0] + 1 << " " << this->v_elem[i].v_nodeID[1] + 1 << " " << this->v_elem[i].v_nodeID[2] + 1 << " " << this->v_elem[i].v_nodeID[3] + 1 << std::endl;
			}
            else if (this->est_shapeType == Element::estMixed)
			{
                if (this->v_elem[i].est_shapeType == Element::estTriangular)
				{
					outFile << this->v_elem[i].v_nodeID[0] + 1 << " " << this->v_elem[i].v_nodeID[1] + 1 << " " << this->v_elem[i].v_nodeID[2] + 1 << " " << this->v_elem[i].v_nodeID[2] + 1 << std::endl;
				}
                else if (this->v_elem[i].est_shapeType == Element::estQuadrilateral)
				{
					outFile << this->v_elem[i].v_nodeID[0] + 1 << " " << this->v_elem[i].v_nodeID[1] + 1 << " " << this->v_elem[i].v_nodeID[2] + 1 << " " << this->v_elem[i].v_nodeID[3] + 1 << std::endl;
				}
			}
		}
	}
	else if (this->md_meshDim == Mesh::md3D)
	{
		if (this->est_shapeType == Element::estPolyhedron)
		{
			int nFaceNodeNum(0);
			for (int i = 0; i < (int)this->v_face.size(); i++)
			{
				for (int j = 0; j < (int)this->v_face[i].v_nodeID.size(); j++)
				{
					nFaceNodeNum++;
				}
			}

			std::cout << "Nodes=" << this->n_nodeNum << ", Faces=" << this->n_faceNum << ", Elements=" << this->n_faceNum << ", ZONETYPE=FEPolyhedron" << std::endl;

			outFile << "TITLE =\"" << "UnGrid Output" << "\"" << std::endl;
			outFile << "VARIABLES = " << std::endl << "\"X\"" << std::endl << "\"Y\"" << std::endl << "\"Z\"" << std::endl;
			outFile << "DATASETAUXDATA Common.VectorVarsAreVelocity = \"TRUE\"" << std::endl;
			outFile << "ZONE T = \"fluid\"" << std::endl;
			outFile << "STRANDID=0, SOLUTIONTIME=0" << std::endl;
			outFile << "Nodes=" << this->n_nodeNum << ", Faces=" << this->n_faceNum << ", Elements=" << this->n_faceNum << ", ZONETYPE=FEPolyhedron" << std::endl;
			outFile << "DATAPACKING=BLOCK" << std::endl;
			outFile << "TotalNumFaceNodes=" << nFaceNodeNum << ", NumConnectedBoundaryFaces=0, TotalNumBoundaryConnections=0" << std::endl;
			outFile << "DT=(SINGLE SINGLE SINGLE )";

			// Write node X
			for (int i = 0; i < (int) this->v_node.size(); i++)
			{
				if (i % 5 == 0) outFile << std::endl;
                outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << this->v_node[i].x_ << " ";
			}

			// Write node Y
			for (int i = 0; i < (int) this->v_node.size(); i++)
			{
				if (i % 5 == 0) outFile << std::endl;
                outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << this->v_node[i].y_ << " ";
			}

			// Write node Z
			for (int i = 0; i < (int) this->v_node.size(); i++)
			{
				if (i % 5 == 0) outFile << std::endl;
                outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << this->v_node[i].z_ << " ";
			}

			outFile << std::endl << "# node count per face";
			for (int i = 0; i < (int) this->v_face.size(); i++)
			{
				if (i % 10 == 0) outFile << std::endl;
				outFile << this->v_face[i].v_nodeID.size() << " ";
			}

			int nFaceNodeCheck(0);
			outFile << std::endl << "# face nodes";
			for (int i = 0; i < (int) this->v_face.size(); i++)
			{
				for (int j = 0; j < (int)this->v_face[i].v_nodeID.size(); j++)
				{
					if (nFaceNodeCheck % 10 == 0) outFile << std::endl;
					outFile << this->v_face[i].v_nodeID[j] << " ";
					nFaceNodeCheck++;
				}
			}

			outFile << std::endl << "# left elements";
			for (int i = 0; i < (int) this->v_face.size(); i++)
			{
				if (i % 10 == 0) outFile << std::endl;
				outFile << this->v_face[i].n_owner << " ";
			}

			outFile << std::endl << "# node count per face";
			for (int i = 0; i < (int) this->v_face.size(); i++)
			{
				if (i % 10 == 0) outFile << std::endl;
				outFile << this->v_face[i].n_neighbor << " ";
			}
		}
		else
		{
			outFile << "VARIABLES = " << "\"x [m]\"," << "\"y [m]\"," << "\"z [m]\"" << std::endl;

			if (this->est_shapeType == Element::estTetrahedral)
			{
				outFile << "ZONE T = " << "\"" << "Tetrahedral_mesh" << "\"," << " DATAPACKING = POINT, N = " << this->n_nodeNum << ", E = " << this->n_elemNum
					<< ", ZONETYPE = FETETRAHEDRON" << std::endl;
			}
			else if (this->est_shapeType == Element::estHexahedral)
			{
				outFile << "ZONE T = " << "\"" << "Hexahedral_mesh" << "\"," << " DATAPACKING = POINT, N = " << this->n_nodeNum << ", E = " << this->n_elemNum
					<< ", ZONETYPE = FEBRICK" << std::endl;
			}
			else if (this->est_shapeType == Element::estMixed)
			{
				outFile << "ZONE T = " << "\"" << "Tet_mesh" << "\"," << " DATAPACKING = POINT, N = " << this->n_nodeNum << ", E = " << this->n_elemNum
					<< ", ZONETYPE = FEBRICK" << std::endl;
			}
			else if (this->est_shapeType == Element::estWedge)
			{
				outFile << "ZONE T = " << "\"" << "Wedge_mesh" << "\"," << " DATAPACKING = POINT, N = " << this->n_nodeNum << ", E = " << this->n_elemNum
					<< ", ZONETYPE = FEBRICK" << std::endl;
			}

			for (int i = 0; i < this->n_nodeNum; i++)
			{
				outFile << this->v_node[i].x_ << " " << this->v_node[i].y_ << " " << this->v_node[i].z_ << std::endl;
			}

			if (this->est_shapeType == Element::estTetrahedral)
			{
				for (int i = 0; i < (int)this->v_elem.size(); i++)
				{
					outFile << this->v_elem[i].v_nodeID[0] + 1 << " " << this->v_elem[i].v_nodeID[1] + 1 << " " << this->v_elem[i].v_nodeID[2] + 1 << " " << this->v_elem[i].v_nodeID[3] + 1 << std::endl;
				}
			}
			else if (this->est_shapeType == Element::estHexahedral)
			{
				for (int i = 0; i < (int)this->v_elem.size(); i++)
				{
					outFile << this->v_elem[i].v_nodeID[0] + 1 << " " << this->v_elem[i].v_nodeID[1] + 1 << " " << this->v_elem[i].v_nodeID[2] + 1 << " " << this->v_elem[i].v_nodeID[3] + 1 << " "
						<< this->v_elem[i].v_nodeID[4] + 1 << " " << this->v_elem[i].v_nodeID[5] + 1 << " " << this->v_elem[i].v_nodeID[6] + 1 << " " << this->v_elem[i].v_nodeID[7] + 1 << std::endl;
				}
			}
			else if (this->est_shapeType == Element::estMixed)
			{
				for (int i = 0; i < (int)this->v_elem.size(); i++)
				{
					if (this->v_elem[i].est_shapeType == Element::estTetrahedral)
					{
						outFile << this->v_elem[i].v_nodeID[0] + 1 << " " << this->v_elem[i].v_nodeID[1] + 1 << " " << this->v_elem[i].v_nodeID[2] + 1 << " " << this->v_elem[i].v_nodeID[2] + 1 << " "
							<< this->v_elem[i].v_nodeID[3] + 1 << " " << this->v_elem[i].v_nodeID[3] + 1 << " " << this->v_elem[i].v_nodeID[3] + 1 << " " << this->v_elem[i].v_nodeID[3] + 1 << std::endl;
					}
					else if (this->v_elem[i].est_shapeType == Element::estPyramid)
					{
						outFile << this->v_elem[i].v_nodeID[0] + 1 << " " << this->v_elem[i].v_nodeID[1] + 1 << " " << this->v_elem[i].v_nodeID[2] + 1 << " " << this->v_elem[i].v_nodeID[3] + 1 << " "
							<< this->v_elem[i].v_nodeID[4] + 1 << " " << this->v_elem[i].v_nodeID[4] + 1 << " " << this->v_elem[i].v_nodeID[4] + 1 << " " << this->v_elem[i].v_nodeID[4] + 1 << std::endl;
					}
					else if (this->v_elem[i].est_shapeType == Element::estWedge)
					{
						outFile << this->v_elem[i].v_nodeID[0] + 1 << " " << this->v_elem[i].v_nodeID[1] + 1 << " " << this->v_elem[i].v_nodeID[4] + 1 << " " << this->v_elem[i].v_nodeID[3] + 1 << " "
							<< this->v_elem[i].v_nodeID[2] + 1 << " " << this->v_elem[i].v_nodeID[2] + 1 << " " << this->v_elem[i].v_nodeID[5] + 1 << " " << this->v_elem[i].v_nodeID[5] + 1 << std::endl;
					}
					else if (this->v_elem[i].est_shapeType == Element::estHexahedral)
					{
						outFile << this->v_elem[i].v_nodeID[0] + 1 << " " << this->v_elem[i].v_nodeID[1] + 1 << " " << this->v_elem[i].v_nodeID[2] + 1 << " " << this->v_elem[i].v_nodeID[3] + 1 << " "
							<< this->v_elem[i].v_nodeID[4] + 1 << " " << this->v_elem[i].v_nodeID[5] + 1 << " " << this->v_elem[i].v_nodeID[6] + 1 << " " << this->v_elem[i].v_nodeID[7] + 1 << std::endl;
					}
				}
			}
			else if (this->est_shapeType == Element::estWedge)
			{
				for (int i = 0; i < (int)this->v_elem.size(); i++)
				{
					outFile << this->v_elem[i].v_nodeID[0] + 1 << " " << this->v_elem[i].v_nodeID[1] + 1 << " " << this->v_elem[i].v_nodeID[5] + 1 << " " << this->v_elem[i].v_nodeID[3] + 1 << " "
						<< this->v_elem[i].v_nodeID[2] + 1 << " " << this->v_elem[i].v_nodeID[2] + 1 << " " << this->v_elem[i].v_nodeID[4] + 1 << " " << this->v_elem[i].v_nodeID[4] + 1 << std::endl;
				}
			}
		}
	}

	outFile.close();

	//system(outMshFileName.c_str());
}


void FluentMeshBlock::Decompose(RegionConnection& regionConnectionInfo)
{
	// Get number of zone, and check if mesh could be decomposed.
	int zoneNum = (int)this->v_elementZone.size();
	std::cout << LogLevel(LogLevel::mlOK, "Current mesh element zone number is: ") << zoneNum << std::endl;


	if (this->v_elementZone.size() < 1)
	{
		WarningContinue(std::string("Element zone number is:") + std::to_string(zoneNum)
			+ std::string(" less than 1. \nPlease Check current mesh."));
	}

	//1. establish global-region two-way index infomations
	std::vector<std::vector<std::pair<int, int> > > globalNodeIndex;
	std::vector<std::vector<std::pair<int, int> > > globalFaceIndex;
	std::vector<std::pair<int, int> > globalElementIndex;
	globalNodeIndex.resize((int)this->v_node.size());
	globalFaceIndex.resize((int)this->v_face.size());
	globalElementIndex.resize((int)this->v_elem.size());

	this->GetZoneNumber(zoneNum);

	this->DistributeElements(globalElementIndex);

	this->DistributeFaces(globalFaceIndex);

	this->DistributeNodes_Vertices(globalNodeIndex);

	this->RewriteTopologicalInformation(globalNodeIndex, globalFaceIndex, globalElementIndex);

	//---------end of rewritting topological infomation---------

	this->DistributeInteriorFaceZones(globalFaceIndex);

	this->SplitExternalFaceZones(globalFaceIndex);

	this->SplitInternalFaceZones(globalFaceIndex, globalElementIndex, regionConnectionInfo);

	this->SideCheck();

	this->WriteFaceID();

	//---distribute region mesh parameters (by the way)---
	for (int i = 0; i < zoneNum; i++)
	{
		this->v_regionGrid[i].md_meshDim = this->md_meshDim;
		this->v_regionGrid[i].est_shapeType = this->v_ElementZoneInfo[i].est_shapeType;
		this->v_regionGrid[i].n_elemNum = (int)this->v_regionGrid[i].v_elem.size();
		this->v_regionGrid[i].n_faceNum = (int)this->v_regionGrid[i].v_face.size();
		this->v_regionGrid[i].n_nodeNum = (int)this->v_regionGrid[i].v_node.size();
	}

	//Clear Global Fluent mesh
	this->ClearGlobalFluentMesh();
}

//Clear Global Fluent mesh
void FluentMeshBlock::ClearGlobalFluentMesh()
{
	this->n_elemNum = 0;
	this->n_faceNum = 0;
	this->n_nodeNum = 0;

    std::vector<Element>().swap(this->v_elem);
    std::vector<Face>().swap(this->v_face);
    std::vector<Node>().swap(this->v_node);
    std::vector<Vertice>().swap(this->v_vertice);

    std::vector<ElementZone>().swap(this->v_elementZone);
    std::vector<FaceZone>().swap(this->v_faceZone);

    std::vector<MshFaceZone>().swap(this->v_FaceZoneInfo);
    std::vector<MshElementZone>().swap(this->v_ElementZoneInfo);
}

void FluentMeshBlock::GetZoneNumber(int &zoneNum)
{
	this->v_regionGrid.resize(zoneNum);
	//---distribute region names (by the way)---
	for (int i = 0; i < zoneNum; i++)
	{
		this->v_regionGrid[i].st_meshName = this->v_elementZone[i].name;
	}
}
void FluentMeshBlock::DistributeElements(std::vector<std::pair<int, int> > &globalElementIndex)
{
	for (int i = 0; i < (int)this->v_elem.size(); i++)
	{
		int regionID = this->v_elem[i].n_ElementZoneID;
		this->v_regionGrid[regionID].v_elem.push_back(this->v_elem[i]);
		//write two-way index informations
		int IDinRegion = (int)this->v_regionGrid[regionID].v_elem.size() - 1;
		globalElementIndex[i] = std::pair<int, int>(regionID, IDinRegion);
	}
}
void FluentMeshBlock::DistributeFaces(std::vector<std::vector<std::pair<int, int> > > &globalFaceIndex)
{
	for (int i = 0; i < (int)this->v_face.size(); i++)
	{
		//4.1 collect associated regions
		std::vector<int> v_regionID;
		std::pair<int, int> oAndN = SearchOwnerNeighbor(i);
		if (oAndN.first != -1)
		{
			int regionID = this->v_elem[oAndN.first].n_ElementZoneID;
			v_regionID.push_back(regionID);
		}
		if (oAndN.second != -1)
		{
			int regionID = this->v_elem[oAndN.second].n_ElementZoneID;
			if (1 == (int)v_regionID.size())
			{
				if (regionID != v_regionID[0])
				{
					v_regionID.push_back(regionID);
				}
			}
		}

		//4.2 copy this face to each of the associated regions
		for (int j = 0; j < (int)v_regionID.size(); j++)
		{
			int regionID = v_regionID[j];
			this->v_regionGrid[regionID].v_face.push_back(this->v_face[i]);

			//write two-way index informations
			int IDinRegion = (int)this->v_regionGrid[regionID].v_face.size() - 1;
			globalFaceIndex[i].push_back(std::pair<int, int>(regionID, IDinRegion));
		}
		
	}
}

void FluentMeshBlock::DistributeNodes_Vertices(std::vector<std::vector<std::pair<int, int> > > &globalNodeIndex)
{
	for (int i = 0; i < (int)this->v_node.size(); i++)
	{
		//5.1 collect associated regions
		std::vector<int> v_regionID;
		for (int j = 0; j < (int)this->v_vertice[i].v_elemID.size(); j++)
		{
			int thisElemID = this->v_vertice[i].v_elemID[j];
			int locatedRegionID = this->v_elem[thisElemID].n_ElementZoneID;
			bool reapeated = false;
			for (int k = 0; k < (int)v_regionID.size(); k++)
			{
				if (locatedRegionID == v_regionID[k])
				{
					reapeated = true;
					break;
				}
			}
			if (false == reapeated)
			{
				v_regionID.push_back(locatedRegionID);
			}
		}
		//5.2 copy this node and vertice to each of the associated regions
		for (int j = 0; j < (int)v_regionID.size(); j++)
		{
			int regionID = v_regionID[j];
			this->v_regionGrid[regionID].v_node.push_back(this->v_node[i]);
			this->v_regionGrid[regionID].v_vertice.push_back(this->v_vertice[i]);
			//write two-way index informations
			int IDinRegion = (int)this->v_regionGrid[regionID].v_node.size() - 1;
			globalNodeIndex[i].push_back(std::pair<int, int>(regionID, IDinRegion));
		}
	}
}
void FluentMeshBlock::RewriteTopologicalInformation(std::vector<std::vector<std::pair<int, int> > > &globalNodeIndex, std::vector<std::vector<std::pair<int, int> > > &globalFaceIndex, std::vector<std::pair<int, int> > &globalElementIndex)
{
	for (int i = 0; i < (int)this->v_regionGrid.size(); i++)
	{
		int regionID = i;
		StandardMeshBlock& thisRegion = this->v_regionGrid[i];
		//6.1 elements in this region-----------
		for (int j = 0; j < (int)thisRegion.v_elem.size(); j++)
		{
			Element& thisElem = thisRegion.v_elem[j];
			//6.1.1 v_faceID in elements
			for (int k = 0; k < (int)thisElem.v_faceID.size(); k++)
			{
				int globalID = thisElem.v_faceID[k];
				for (int m = 0; m < (int)globalFaceIndex[globalID].size(); m++)
				{
					if (regionID == globalFaceIndex[globalID][m].first)
					{
						thisElem.v_faceID[k] = globalFaceIndex[globalID][m].second;
						break;
					}
				}
			}
			//6.1.2 v_nodeID in elements
			for (int k = 0; k < (int)thisElem.v_nodeID.size(); k++)
			{
				int globalID = thisElem.v_nodeID[k];
				for (int m = 0; m < (int)globalNodeIndex[globalID].size(); m++)
				{
					if (regionID == globalNodeIndex[globalID][m].first)
					{
						thisElem.v_nodeID[k] = globalNodeIndex[globalID][m].second;
						break;
					}
				}
			}
		}

		//6.2 -----------faces in this region-------------
		for (int j = 0; j < (int)thisRegion.v_face.size(); j++)
		{
			Face& thisFace = thisRegion.v_face[j];
			//6.2.1 v_nodeID in this face
			for (int k = 0; k < (int)thisFace.v_nodeID.size(); k++)
			{
				int globalID = thisFace.v_nodeID[k];
				for (int m = 0; m < (int)globalNodeIndex[globalID].size(); m++)
				{
					if (regionID == globalNodeIndex[globalID][m].first)
					{
						thisFace.v_nodeID[k] = globalNodeIndex[globalID][m].second;
						break;
					}
				}
			}
			//6.2.2 n_owner
			int globalID = thisFace.n_owner;
			if (-1 != globalID)
			{
				int ownerRegionID = globalElementIndex[globalID].first;
				if (ownerRegionID == regionID)
				{
					thisFace.n_owner = globalElementIndex[globalID].second;
				}
				else
				{
					thisFace.n_owner = -1;
				}
			}
			//6.2.3 n_neighbor
			globalID = thisFace.n_neighbor;
			if (-1 != globalID)
			{
				int nbRegionID = globalElementIndex[globalID].first;
				if (nbRegionID == regionID)
				{
					thisFace.n_neighbor = globalElementIndex[globalID].second;
				}
				else
				{
					thisFace.n_neighbor = -1;
				}
			}
		}

		//6.3-------vertices in this region--------
		for (int j = 0; j < (int)thisRegion.v_vertice.size(); j++)
		{
			Vertice& thisVertice = thisRegion.v_vertice[j];
			std::vector<int> tempList;
			//6.3.1 element
			//(1)----------remove elemIDs not belonging to this region------------
			for (int k = 0; k < (int)thisVertice.v_elemID.size(); k++)
			{
				int globalID = thisVertice.v_elemID[k];
				if (globalElementIndex[globalID].first == regionID)
				{
					tempList.push_back(globalID);
				}
			}
			thisVertice.v_elemID.clear();
			for (int k = 0; k < (int)tempList.size(); k++)
			{
				thisVertice.v_elemID.push_back(tempList[k]);
			}
			tempList.clear();
			//(2)----------modify elemIDs left in this region------------
			for (int k = 0; k < (int)thisVertice.v_elemID.size(); k++)
			{
				int globalID = thisVertice.v_elemID[k];
				thisVertice.v_elemID[k] = globalElementIndex[globalID].second;
			}
			//6.3.2 face
			//(1)----------remove faceIDs not belonging to this region------------
			for (int k = 0; k < (int)thisVertice.v_faceID.size(); k++)
			{
				int globalID = thisVertice.v_faceID[k];
				for (int m = 0; m < (int)globalFaceIndex[globalID].size(); m++)
				{
					if (globalFaceIndex[globalID][m].first == regionID)
					{
						tempList.push_back(globalID);
						break;
					}
				}
			}
			thisVertice.v_faceID.clear();
			for (int k = 0; k < (int)tempList.size(); k++)
			{
				thisVertice.v_faceID.push_back(tempList[k]);
			}
			tempList.clear();
			//(2)----------modify faceIDs left in this region------------
			for (int k = 0; k < (int)thisVertice.v_faceID.size(); k++)
			{
				int globalID = thisVertice.v_faceID[k];
				for (int m = 0; m < (int)globalFaceIndex[globalID].size(); m++)
				{
					if (globalFaceIndex[globalID][m].first == regionID)
					{
						thisVertice.v_faceID[k] = globalFaceIndex[globalID][m].second;
						break;
					}
				}
			}
		}
	}
}
void FluentMeshBlock::DistributeInteriorFaceZones(std::vector<std::vector<std::pair<int, int> > > &globalFaceIndex)
{
	for (int i = 0; i < (int)this->v_faceZone.size(); i++)
	{
		if (FaceZone::fztInterior == v_faceZone[i].faceZoneType)
		{
			for (int j = 0; j < (int)this->v_faceZone[i].v_faceID.size(); j++)
			{
				int FaceID = v_faceZone[i].v_faceID[j];
				int regionID = globalFaceIndex[FaceID][0].first;
				StandardMeshBlock& thisRegion = this->v_regionGrid[regionID];

				int globalID = this->v_faceZone[i].v_faceID[j];
				int localID = globalFaceIndex[globalID][0].second;
				thisRegion.fz_interiorFaceZone.v_faceID.push_back(localID);
			}
		}
	}
}

void FluentMeshBlock::SplitExternalFaceZones(std::vector<std::vector<std::pair<int, int> > > &globalFaceIndex)
{
	for (int i = 0; i < (int)this->v_faceZone.size(); i++)
	{
		if (FaceZone::fztExternal == v_faceZone[i].faceZoneType)
		{
			std::vector<int> regionIDList;
			std::vector<FaceZone> tempFaceZoneList;
			for (int j = 0; j < (int)this->v_faceZone[i].v_faceID.size(); j++)
			{
				int globalID = this->v_faceZone[i].v_faceID[j];
				int regionID = globalFaceIndex[globalID][0].first;
				int localID = globalFaceIndex[globalID][0].second;
				bool faceZoneCreated = false;
				//check whether the detected zone has the corresponding face zone created;
				for (int k = 0; k < (int)regionIDList.size(); k++)
				{
					if (regionID == regionIDList[k])
					{
						faceZoneCreated = true;
						break;
					}
				}
				//create face zone at region, if not found
				if (false == faceZoneCreated)
				{
					FaceZone tempFZ;
					//continue from here
					tempFZ.name = this->v_faceZone[i].name + "_" + this->v_regionGrid[regionID].st_meshName;
					tempFaceZoneList.push_back(tempFZ);
					regionIDList.push_back(regionID);
				}

				//search id in tempFaceZoneList
				for (int k = 0; k < (int)regionIDList.size(); k++)
				{
					if (regionID == regionIDList[k])
					{
						//distribute the specific faceID (local) to its host regional face zone (created)
						tempFaceZoneList[k].v_faceID.push_back(localID);
						break;
					}
				}
			}
			for (int j = 0; j < (int)regionIDList.size(); j++)
			{
				this->v_regionGrid[regionIDList[j]].v_boundaryFaceZone.push_back(tempFaceZoneList[j]);
			}
		}
	}
}

void FluentMeshBlock::SplitInternalFaceZones(std::vector<std::vector<std::pair<int, int> > > &globalFaceIndex, std::vector<std::pair<int, int> > &globalElementIndex, RegionConnection& regionConnectionInfo)
{
	for (int i = 0; i < (int)this->v_faceZone.size(); i++)
	{
		if (FaceZone::fztInternal == v_faceZone[i].faceZoneType)
		{
			int nZoneNum(0);
			if (this->v_faceZone[i].v_faceID.size())
			{
				int nTemp = this->v_faceZone[i].v_faceID[0];
				if (1 == globalFaceIndex[nTemp].size())
				{
					nZoneNum = 1;
				}
				else if (2 == globalFaceIndex[nTemp].size())
				{
					nZoneNum = 2;
				}
			}
			// Internal face in one region 
			if (1 == nZoneNum)
			{
				int nTemp = this->v_faceZone[i].v_faceID[0];
				int regionID = globalFaceIndex[nTemp][0].first;
				FaceZone headFaceZone;
				headFaceZone.name = this->v_faceZone[i].name + "_HEAD";
				FaceZone tailFaceZone;
				tailFaceZone.name = this->v_faceZone[i].name + "_TAIL";

				for (int j = 0; j < (int)this->v_faceZone[i].v_faceID.size(); j++)
				{
					int globalFaceID = this->v_faceZone[i].v_faceID[j];
					int localFaceID = globalFaceIndex[globalFaceID][0].second;
					Face tempFace = this->v_regionGrid[regionID].v_face[localFaceID];
					this->v_regionGrid[regionID].v_face[localFaceID].n_neighbor = -1;
					tempFace.n_owner = -1;
					this->v_regionGrid[regionID].v_face.push_back(tempFace);
					headFaceZone.v_faceID.push_back(localFaceID);
					int newFaceID = this->v_regionGrid[regionID].v_face.size() - 1;
					globalFaceIndex[globalFaceID].push_back(std::pair<int, int>(regionID, newFaceID));
					tailFaceZone.v_faceID.push_back(newFaceID);
				}
				this->v_regionGrid[regionID].v_boundaryFaceZone.push_back(headFaceZone);
				this->v_regionGrid[regionID].v_boundaryFaceZone.push_back(tailFaceZone);
			}
			else if (2 == nZoneNum)
			{
				std::vector<int> regionIDList;
				std::vector<FaceZone> tempFaceZoneList;
				for (int j = 0; j < (int)this->v_faceZone[i].v_faceID.size(); j++)
				{
					int globalID = this->v_faceZone[i].v_faceID[j];
					for (int k = 0; k < 2; k++)
					{
						int regionID = globalFaceIndex[globalID][k].first;
						int localID = globalFaceIndex[globalID][k].second;
						bool faceZoneCreated = false;
						//check whether the detected zone has the corresponding face zone created;
						for (int m = 0; m < (int)regionIDList.size(); m++)
						{
							if (regionID == regionIDList[m])
							{
								faceZoneCreated = true;
								break;
							}
						}
						//create face zone at region, if not found
						if (false == faceZoneCreated)
						{
							FaceZone tempFZ;
							//continue from here
							tempFZ.name = this->v_faceZone[i].name + "_" + this->v_regionGrid[regionID].st_meshName;
							tempFaceZoneList.push_back(tempFZ);
							regionIDList.push_back(regionID);
						}
						//search id in tempFaceZoneList
						
						for (int m = 0; m < (int)regionIDList.size(); m++)
						{
							
							if (regionID == regionIDList[m]) 
							{ 
								tempFaceZoneList[m].v_faceID.push_back(localID); 
								break;
							}
						}
						//distribute the specific faceID (local) to its host regional face zone (created)
						
					}
				}
				for (int j = 0; j < (int)regionIDList.size(); j++)
				{
					this->v_regionGrid[regionIDList[j]].v_boundaryFaceZone.push_back(tempFaceZoneList[j]);
				}
			}
		}
	}

	//9.2 write connection information into given RegionConnection
	for (int i = 0; i < (int)this->v_faceZone.size(); i++)
	{
		if (FaceZone::fztInternal == v_faceZone[i].faceZoneType)
		{
			FaceZone& thisFaceZone = this->v_faceZone[i];
			//Get region IDs, boundary names of yang and yin faces from a sample face
			int sampleFaceID = this->v_faceZone[i].v_faceID[0];
			std::pair<int, int> yangLocalIndex = globalFaceIndex[sampleFaceID][0];
			std::pair<int, int> yinLocalIndex = globalFaceIndex[sampleFaceID][1];
			if (yangLocalIndex.first>yinLocalIndex.first)
			{
				yangLocalIndex = globalFaceIndex[sampleFaceID][1];
				yinLocalIndex = globalFaceIndex[sampleFaceID][0];
			}
			int yangRegionID = yangLocalIndex.first;
			std::string yangBCName;
			for (int j = 0; j < (int)this->v_regionGrid[yangRegionID].v_boundaryFaceZone.size(); j++)
			{
				int found = -1;
				for (int k = 0; k < (int)this->v_regionGrid[yangRegionID].v_boundaryFaceZone[j].v_faceID.size(); k++)
				{
					if (yangLocalIndex.second == this->v_regionGrid[yangRegionID].v_boundaryFaceZone[j].v_faceID[k])
					{
						found = j;
						break;
					}
				}
				if (-1 != found)
				{
					yangBCName = this->v_regionGrid[yangRegionID].v_boundaryFaceZone[found].name;
					break;
				}
			}
			int yinRegionID = yinLocalIndex.first;
			std::string yinBCName;
			for (int j = 0; j < (int)this->v_regionGrid[yinRegionID].v_boundaryFaceZone.size(); j++)
			{
				int found = -1;
				for (int k = 0; k < (int)this->v_regionGrid[yinRegionID].v_boundaryFaceZone[j].v_faceID.size(); k++)
				{
					if (yinLocalIndex.second == this->v_regionGrid[yinRegionID].v_boundaryFaceZone[j].v_faceID[k])
					{
						found = j;
						break;
					}
				}
				if (-1 != found)
				{
					yinBCName = this->v_regionGrid[yinRegionID].v_boundaryFaceZone[found].name;
					break;
				}
			}
			//Create a new coupled boundary pair
			Mesh* pmesh1 = &(this->v_regionGrid[yangRegionID]);
			Mesh* pmesh2 = &(this->v_regionGrid[yinRegionID]);
			CoupledBoundary tempCB(pmesh1, pmesh2, yangBCName, yinBCName);
			//Write local face IDs in the created boundary pair
			for (size_t j = 0; j < this->v_faceZone[i].v_faceID.size(); j++)
			{
				int faceID = this->v_faceZone[i].v_faceID[j];
				std::pair<int, int> yangLocalIndex = globalFaceIndex[faceID][0];
				std::pair<int, int> yinLocalIndex = globalFaceIndex[faceID][1];
				if (yangLocalIndex.first>yinLocalIndex.first)
				{
					yangLocalIndex = globalFaceIndex[faceID][1];
					yinLocalIndex = globalFaceIndex[faceID][0];
				}
				tempCB.v_coupledFaceID.push_back(std::pair<int, int>(yangLocalIndex.second, yinLocalIndex.second));
			}
			//insert the created boundary pair into given regionConnectionInfo
			regionConnectionInfo.v_coupledBoundary.push_back(tempCB);
		}
	}
	return;
}

void FluentMeshBlock::WriteFaceID()
{
	for (int i = 0; i < (int)this->v_regionGrid.size(); i++)
	{
		//specify face zone ID as -1 for interior faces
		for (int k = 0; k < (int)this->v_regionGrid[i].fz_interiorFaceZone.v_faceID.size(); k++)
		{
			int faceID = this->v_regionGrid[i].fz_interiorFaceZone.v_faceID[k];
			this->v_regionGrid[i].v_face[faceID].n_FaceZoneID = -1;
		}
		//for a boundary face, specify face zone ID as the ID of the boundary which it is located on
		for (int j = 0; j < (int)this->v_regionGrid[i].v_boundaryFaceZone.size(); j++)
		{
			for (int k = 0; k < (int)this->v_regionGrid[i].v_boundaryFaceZone[j].v_faceID.size(); k++)
			{
				int faceID = this->v_regionGrid[i].v_boundaryFaceZone[j].v_faceID[k];
				this->v_regionGrid[i].v_face[faceID].n_FaceZoneID = j;
			}
		}
	}
}

void FluentMeshBlock::SideCheck()
{
	for (int i = 0; i < (int)this->v_regionGrid.size(); i++)
	{
		for (int j = 0; j < (int)this->v_regionGrid[i].v_boundaryFaceZone.size(); j++)
		{
			FaceZone& thisFaceZone = this->v_regionGrid[i].v_boundaryFaceZone[j];
			for (int k = 0; k < (int)thisFaceZone.v_faceID.size(); k++)
			{
				int faceID = thisFaceZone.v_faceID[k];
				this->v_regionGrid[i].v_face[faceID].CorrectSide();
			}
		}
	}

	for (int i = 0; i < (int)this->v_regionGrid.size(); i++)
	{
		this->v_regionGrid[i].CorrectSide();
	}
	return;
}

