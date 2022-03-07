#include "MOCMesh.h"
#include <iostream> 
#include <fstream>
#include <map>
#include <tuple>
#include <regex>
#include <unordered_set>
#include "Logger.h"
#include "index.h"
#include "Solver.h"
#include "ConfigurationFile.h"
#ifdef _WIN32
#include "boost/boost/regex.hpp"
#else
#include "boost/regex.hpp"
#endif
//using namespace std;

//#define PI 3.14159265358979323846
#define NA 6.022e23
#define BARN 1.0e-24



MOCMesh::MOCMesh(const std::vector<std::string>& vMaterailName)
{
		m_firstCreated = false;
		readMapFile(vMaterailName);
		GetAllMocIndex(m_vSMocIndex);
}

void MOCMesh::reWriteAplOutputFile(std::string outAplFileName)
{
	RenameFile(outAplFileName, outAplFileName + "_temp");
	ofstream outFile_new(outAplFileName+"_new");
	ifstream infile_(outAplFileName + "_temp");

	ofstream outFile(outAplFileName);
	ifstream infile(outAplFileName + "_new");// +"_temp");
	
	std::stringstream firstNt_Assembly,second_Nt_Assembly, end_Nt_Assembly,outAssembly;
	//std::vector<std::stringstream> vSecond_Nt_Assembly;
	vector< shared_ptr<stringstream>> vSecond_Nt_Assembly;
	if (!infile_.is_open())
	{
		Logger::LogError("cannot find the MOC mesh file:" + outAplFileName);
		exit(EXIT_FAILURE);
	}
	string line;

	while (getline(infile_, line))  //read mesh data
	{
	loop_Nt_Assembly:
		if (line.find("*Nt_Assembly")!=line.npos && line.find("*Nt_Assembly_") == line.npos)
		{
			second_Nt_Assembly.str("");
			second_Nt_Assembly<< line<<endl;
			while (getline(infile_, line))
			{ 
				if (line.find("*Boundary_condition_number") != line.npos)
				{
					end_Nt_Assembly << line << endl;
					while (getline(infile_, line))
					{
						if (line.find("Nt_Assembly") != line.npos && line.find("*Nt_Assembly_") == line.npos)
						{
							vSecond_Nt_Assembly.push_back(make_shared<stringstream>());
							*vSecond_Nt_Assembly.back()<< second_Nt_Assembly.str();
							goto loop_Nt_Assembly;
						}
						end_Nt_Assembly << line << endl;
					}
				}
				else
					second_Nt_Assembly << line << endl;
			}
			vSecond_Nt_Assembly.push_back(make_shared<stringstream>());
			*vSecond_Nt_Assembly.back() << second_Nt_Assembly.str();
			//vSecond_Nt_Assembly.push_back(make_shared<stringstream>(second_Nt_Assembly));
		}
		else
			firstNt_Assembly << line << endl;
	}
	outFile_new << firstNt_Assembly.str();
	for (int i = 0; i < m_vAssembly.size(); i++)
	{
		outFile_new << vSecond_Nt_Assembly[m_vAssembly[i].iAssemblyType-1]->str() << endl;
	}
	outFile_new << end_Nt_Assembly.str();
	outFile_new.close();
	infile_.close();
	
	while (getline(infile, line))  //read mesh data
	{
		outFile << line << endl;
		//stringstream stringline(line);
		//string token;
		//stringline >> token;
		if (line.find("*Nt_Assembly_Cat_Num")!=std::string::npos)
		{
			getline(infile, line);
			outFile << m_vAssembly.size() << endl;
			int iAssemblyIndex = 1;
			while (getline(infile, line))
			{
				outFile << line << endl;
			Nt_Assembly_line:
				if (line.find("*Nt_Assembly") != std::string::npos) 
				{
					//outAssembly << line << endl;
					getline(infile, line);
					outFile << iAssemblyIndex << endl;
					//outAssembly << iAssemblyIndex << endl;
					while (getline(infile, line))
					{
						outFile << line << endl;
						outAssembly << line << endl;
					Material_name_of_each_mesh_line:
						if (line.find("*Material_name_of_each_mesh") != std::string::npos)
						{
							int ID_order = 1;
							string tokenMaterialType = "";
							int out0 = 1;
							while (out0 != 0)
							{
								std::streampos Material_pos;
								//Material_pos = infile.tellg();
								getline(infile, line);
								stringstream stringlineMaterialType(line);
								while (stringlineMaterialType >> tokenMaterialType)
								{
									if (tokenMaterialType.find("*") != std::string::npos)
									{
										out0 = 0;
										//stringlineMaterialType >> token;
										//infile.seekg(Material_pos);
										break;
									}

									if (tokenMaterialType.find("mMOD_") != tokenMaterialType.npos)
									{

										tokenMaterialType += "_" + std::to_string(iAssemblyIndex);
									}
									outFile << tokenMaterialType << "  ";
									outAssembly << tokenMaterialType << "  ";
								}
								if (out0 != 0)
								{
									outFile << endl;
									outAssembly << endl;
								}
								else
								{
									outFile << line << endl;
									outAssembly << line << endl;
									goto Material_name_of_each_mesh_line;
									
								}
							}

						}
						else if (line.find("*Temperature_name_of_each_mesh") != std::string::npos)
						{
							//int ID_order = 1;
							string tokenTemperatureName = "";
							int out0 = 1;
							while (out0 != 0)
							{
								//assembly_pos = infile.tellg();
								getline(infile, line);
								stringstream stringlineTemperatureName(line);
								while (stringlineTemperatureName >> tokenTemperatureName)
								{
									if (tokenTemperatureName.find("*") != std::string::npos)
									{
										out0 = 0;
										//stringlineTemperatureName >> token;
										break;
									}

									if (tokenTemperatureName.find("tMOD_")!= tokenTemperatureName.npos
										|| tokenTemperatureName.find("tFUEL_")!= tokenTemperatureName.npos)
									{
										tokenTemperatureName += "_" + std::to_string(iAssemblyIndex);
										//outFile << tokenTemperatureName << "_" << ID_order << "  ";
										//ID_order++;
									}

									outAssembly << tokenTemperatureName << "  ";
									outFile << tokenTemperatureName << "  ";
									//meshTemperatureNameTemperary.push_back(tokenTemperatureName);

								}
								if (out0 != 0)
								{
									outFile << endl;
								}
								else
								{
									outFile << line << endl;
									outAssembly << line << endl;
									goto Material_name_of_each_mesh_line;
								}
							}
						}
						else if (line.find("*Nt_Assembly") != std::string::npos)
						{
							iAssemblyIndex++;
							goto Nt_Assembly_line;
						}
					}
				}
			}
		}
	}

	
	infile.close();
	outFile.close();
	//RemoveFile(outAplFileName + "_temp");

}

MOCMesh::MOCMesh(std::string meshFileName, std::string outAplFileName, MeshKernelType kernelType)
{
	Logger::LogInfo("MOCMesh generation begin");

	ofstream outFile(outAplFileName);
	ifstream infile(meshFileName);
	if (!infile.is_open())
	{
		Logger::LogError("cannot find the MOC mesh file:" + meshFileName);
		exit(EXIT_FAILURE);
	}
	string line;
	m_pAssemblyIndex = std::make_shared<AssemblyIndex>(*this);

	int xDirection_Number = 0, yDirection_Number = 0;//assembly numble on x,y direction
	int iNt_Assembly_index = 0;
	std::streampos strpos, assembly_pos;
	bool bNextAssembly = false;
	while (getline(infile, line))  //read mesh data
	{
		outFile << line << endl;
		stringstream stringline(line);
		string token;
		while (stringline >> token)
		{
			if (token == "*")
			{
				stringline >> token;
				if (token == "Mesh")
				{
					getline(infile, line);
					outFile << line << endl;
					setMeshInformation(line);
					break;
				}
				if (token == "AXIAL")
				{
					getline(infile, line);
					outFile << line << endl;
					setAxialInformation(line);

					m_vAssemblyType.resize(1);
					m_vAssemblyType[0].iAssemblyType = 1;
					m_vAssembly.resize(1);
					m_vAssembly[0].vAssembly_LeftDownPoint = Vector(0, 0, 0);
					m_vAssembly[0].iAssemblyType = 1;
					m_pAssemblyIndex->m_assemblyIndex.resize(1);
					m_pAssemblyIndex->m_assemblyIndex[0].resize(1);
					m_pAssemblyIndex->m_assemblyIndex[0][0] = 0;
					//vNumber_of_each_coarse_mesh.resize(1);
					//vNumber_of_each_coarse_mesh[0] = layerMeshNum;
					break;
				}
				if (token == "EDGE")
				{
					//strpos = infile.tellg();
					infile.seekg(strpos);
					goto loop;
				}
			}
			if (token == "*Core_General")
			{
				m_bSingleCell = false;
				getline(infile, line);
				outFile << line << endl;
				stringstream stringline(line);
				string token;
				stringline >> token;
				stringline >> token;
				xDirection_Number = std::stod(token);
				stringline >> token;
				yDirection_Number = std::stod(token);
				m_vAssembly.resize(xDirection_Number * yDirection_Number);
				break;
			}
			if (token == "*Core_Nt_Axial_Mesh")
			{
				getline(infile, line);
				outFile << line << endl;
				setAxialInformation(line);
				break;
			}
			if (token == "*Assembly_Positon_Index")
			{
				//
				break;
			}
			if (token == "*Nt_Assembly_Layout")
			{
				int k = 0;
				m_pAssemblyIndex->m_assemblyIndex.resize(xDirection_Number);
				for (int i = 0; i < xDirection_Number; i++)
					m_pAssemblyIndex->m_assemblyIndex[i].resize(yDirection_Number);
				int kkk = 1;
				for (int i = 0; i < yDirection_Number; i++)
				{
					getline(infile, line);
					//outFile << line << endl;
					stringstream stringnumline(line);
					string tokennum;
					
					for (int j = 0; j < xDirection_Number; j++)
					{
						stringnumline >> tokennum;
						m_vAssembly[k].iAssemblyType = std::stod(tokennum);
						m_pAssemblyIndex->m_assemblyIndex[j][yDirection_Number - i - 1] = k;
						k++;
						outFile << kkk++ << "   ";
					}
					outFile << std::endl;
				}
				break;
			}
			if (token == "*Nt_Assembly_Cat_Num")
			{
				getline(infile, line);
				outFile << line << endl;
				stringstream stringline(line);
				string token;
				stringline >> token;
				m_vAssemblyType.resize(std::stod(token));
				break;
			}
			if (token == "*Nt_Assembly")
			{
			loop:
				bNextAssembly = false;
				vector<string> meshMaterialNameTemperary;
				vector<string> meshTemperatureNameTemperary;
				vector<string> meshFaceTypeTemperary;
				vector<string> meshFaceTemperatureNameTemperary;
				vector<string> fileNameTemperary;
				vector< shared_ptr<stringstream>> vStreamTemperay;
				vector<int> meshIDTemperary;
				int nFineMesh = kernelType == MeshKernelType::CGAL_KERNEL ? 4 : 1;  //fine mesh
				std::vector<Surface>allMeshFaces;   //all face objects
				std::vector<MOCEdge>allEdges;    //all edge objects
				std::vector<int> vNumber_of_each_coarse_mesh;

				if (!m_bSingleCell)
				{
					getline(infile, line);
					outFile << line << endl;
					stringstream stringline(line);
					string token;
					stringline >> token;
					m_vAssemblyType[iNt_Assembly_index].iAssemblyType = std::stod(token);
				}
				else
				{
					vNumber_of_each_coarse_mesh.resize(1);
					vNumber_of_each_coarse_mesh[0] = layerMeshNum;
				}

				while (getline(infile, line))
				{
					outFile << line << endl;
					stringstream stringline(line);
					string token;
					stringline >> token;
					if (token == "*RECTANGULAR")
					{
						stringline >> token;
						stringline >> token;
						m_vAssemblyType[iNt_Assembly_index].xLength = std::stod(token)/10;
						stringline >> token;
						m_vAssemblyType[iNt_Assembly_index].yLength = std::stod(token)/10;

					}
					if (token == "*Mesh_number_of_each_coarse_mesh")//only on multi  cell 
					{
						std::streampos pos;
						while (getline(infile, line))
						{
							outFile << line << endl;
							if (line.find("*") == std::string::npos)
							{
								stringstream stringline(line);
								string token;
								while (stringline >> token)
								{
									vNumber_of_each_coarse_mesh.push_back(stod(token));
								}
							}
							else
							{
								//infile.seekg(pos);
								break;
							}
							//pos = infile.tellg();
						}
					}
					if (line.find("Mesh_no_of_each_coarse_mesh") != std::string::npos)
					{
						string tokenMeshId = "";
						int out0 = 1;
						std::streampos pos;
						while (out0 != 0)
						{
							pos = infile.tellg();
							getline(infile, line);
							outFile << line << endl;
							stringstream stringlineMeshID(line);
							while (stringlineMeshID >> tokenMeshId)
							{
								if (tokenMeshId.find("*") != std::string::npos)
								{
									out0 = 0;
									stringlineMeshID >> token;
									//infile.seekg(pos);
									break;
								}
								meshIDTemperary.push_back(stod(tokenMeshId));
							}
						}
					}
					if (line.find("EDGE")!=std::string::npos && line.find("*EDGE_")==std::string::npos)
					{
						if (m_bSingleCell)
						{
							stringline >> token;
							stringline >> token;
						}
						else
						{
							stringline >> token;
						}
						int edgeIDTemperary = stod(token);
						string lineType;
						string linePosition;
						getline(infile, lineType);
						outFile << lineType << endl;
						getline(infile, line);
						outFile << line << endl;
						getline(infile, linePosition);
						outFile << linePosition << endl;
						setEdgeInformation(lineType, linePosition, edgeIDTemperary, allEdges, nFineMesh);
						continue;
					}
					if (line.find("Material_name_of_each_mesh") != std::string::npos)
					{
						int ID_order = 1;
						string tokenMaterialType = "";
						int out0 = 1;
						while (out0 != 0)
						{
							std::streampos Material_pos;
							//Material_pos = infile.tellg();
							getline(infile, line);
							stringstream stringlineMaterialType(line);
							while (stringlineMaterialType >> tokenMaterialType)
							{
								if (tokenMaterialType.find("*") != std::string::npos)
								{
									out0 = 0;
									stringlineMaterialType >> token;
									//infile.seekg(Material_pos);
									break;
								}
								
								if (tokenMaterialType == "mMOD")
								{

									tokenMaterialType += "_" + std::to_string(ID_order);
								}
								meshMaterialNameTemperary.push_back(tokenMaterialType);
								outFile << tokenMaterialType << "  ";
								ID_order++;
							}
							if (out0 != 0)
							{
								outFile << endl;
							}
							else
							{
								outFile << line << endl;
							}
						}
					}
					if (line.find("Temperature_name_of_each_mesh") != std::string::npos)
					{
						int ID_order = 1;
						string tokenTemperatureName = "";
						int out0 = 1;
						while (out0 != 0)
						{
							assembly_pos = infile.tellg();
							getline(infile, line);
							stringstream stringlineTemperatureName(line);
							while (stringlineTemperatureName >> tokenTemperatureName)
							{
								if (tokenTemperatureName.find("*") != std::string::npos)
								{
									out0 = 0;
									stringlineTemperatureName >> token;
									break;
								}
								
								if (tokenTemperatureName == "tMOD" || tokenTemperatureName == "tFUEL")
								{
									tokenTemperatureName += "_" + std::to_string(ID_order);
									//outFile << tokenTemperatureName << "_" << ID_order << "  ";
									ID_order++;
								}
								
								outFile << tokenTemperatureName << "  ";
								meshTemperatureNameTemperary.push_back(tokenTemperatureName);
								
							}
							if (out0 != 0)
							{
								outFile << endl;
							}
							else
							{
								outFile << line << endl;
							}
						}
					}
					if (line.find("*Nt_Assembly") != std::string::npos)
					{
						//infile.seekg(assembly_pos);
						//outFile << line << endl;
						bNextAssembly = true;
						break;
					}
					assembly_pos = infile.tellg();

				}

				//generate mesh for every assembly type
				
				int iTotalMeshNum = 0;
				for (int i = 0; i < vNumber_of_each_coarse_mesh.size(); i++)
				{
					iTotalMeshNum = iTotalMeshNum + vNumber_of_each_coarse_mesh[i];
				}
				layerMeshNum = iTotalMeshNum;
				
				for (int j = 0; j < axialNum; j++)
				{
					for (int i = 0; i < layerMeshNum; i++)
					{
						int meshIDtemp_ = meshIDTemperary[i] + j * layerMeshNum - 1;
						meshFaceTypeTemperary.push_back(meshMaterialNameTemperary[meshIDtemp_]);
						meshFaceTemperatureNameTemperary.push_back(meshTemperatureNameTemperary[meshIDtemp_]);
					}
				}
				vStreamTemperay.resize(layerMeshNum* axialNum);
				setMeshFaceInformation(meshIDTemperary, meshFaceTypeTemperary, meshFaceTemperatureNameTemperary, allMeshFaces, allEdges);
				ThreeDemMeshOutput(vStreamTemperay, allMeshFaces, meshFaceTypeTemperary, nFineMesh);
				//ThreeDemMeshOutput(vStreamTemperay, allMeshFaces, meshMaterialNameTemperary, nFineMesh);

				m_vAssemblyType[iNt_Assembly_index].v_Cell.resize(vNumber_of_each_coarse_mesh.size());
				//m_meshPointPtrVec.resize(axialNum * layerMeshNum);
				int iMeshID_index = 0;
				for (int k = 0; k < vNumber_of_each_coarse_mesh.size(); k++)
				{
					m_vAssemblyType[iNt_Assembly_index].v_Cell[k].vMeshPointPtrVec.resize(axialNum * vNumber_of_each_coarse_mesh[k]);
					int iMeshpointIndex = 0;
					for (int j = 0; j < axialNum; j++)
					{
						for (int i = 0; i < vNumber_of_each_coarse_mesh[k]; i++)
						{
							int meshIDtemp_ = meshIDTemperary[iMeshID_index + i] + j * layerMeshNum;
							//int iMeshpointIndex = i + j * vNumber_of_each_coarse_mesh[k];
							int index = iMeshID_index + i + j * layerMeshNum;

							if (kernelType == MeshKernelType::MHT_KERNEL) {
								Vector point, norm;
								for (auto& edge : allMeshFaces[iMeshID_index+i].faceEdges) {
									if (edge.edgeType == 3) {
										point = edge.arcCenter;
										norm = edge.arcAxisDir;
										break;
									}
								}
								//std::ifstream ifs(fileNameTemperary[index]);
								std::stringstream& ifs = *vStreamTemperay[index];
								m_vAssemblyType[iNt_Assembly_index].v_Cell[k].vMeshPointPtrVec[iMeshpointIndex++] = std::make_shared<MHTMocMeshPoint>(
									meshIDtemp_, ifs, allMeshFaces[iMeshID_index+i].curveInfo, point, norm,
									meshFaceTypeTemperary[index], meshFaceTemperatureNameTemperary[index]);
							}
							/*
							const char* removeFile = fileNameTemperary[index].data();
							if (remove(removeFile)) //delete file
							{
								cout << "delete file fail" << endl;
							}*/
						}
					}
					iMeshID_index = iMeshID_index + vNumber_of_each_coarse_mesh[k];
				}
				
				iNt_Assembly_index++;
				if (bNextAssembly)
					goto loop;
			}
			break;
		}
		strpos = infile.tellg();
	}
	outFile.close();
	InitAssembly();
	m_pAssemblyIndex->buildIndex();
	GetAllMocIndex(m_vSMocIndex);
	reWriteAplOutputFile(outAplFileName);

	Logger::LogInfo("MOCMesh generation ends");
}

void MOCMesh::GetAllMocIndex(std::vector< SMocIndex>& vSMocIndex)
{
	SMocIndex sTemp;
	if (m_firstCreated)
	{
		for (int i = 0; i < m_vAssembly.size(); i++)
		{
			//m_vAssembly[i].v_field.resize(m_vAssembly[i].pAssembly_type->v_Cell.size());
			for (int j = 0; j < m_vAssembly[i].pAssembly_type->v_Cell.size(); j++)
			{
				//m_vAssembly[i].v_field[j].resize(64);
				for (int k = 0; k < m_vAssembly[i].pAssembly_type->v_Cell[j].vMeshPointPtrVec.size(); k++)
				{
					sTemp.iAssemblyIndex = i;
					sTemp.iCellIndex = j;
					sTemp.iMocIndex = k;
					vSMocIndex.push_back(sTemp);
				}
			}
		}
	}
	else
	{
		for (int i = 0; i < m_vAssemblyField.size(); i++)
		{
			for (int j = 0; j < m_vAssemblyField[i].size(); j++)
			{
				for (int k = 0; k < m_vAssemblyField[i][j].size(); k++)
				{
					if (m_vAssemblyField[i][j][k])
						vSMocIndex.push_back(SMocIndex(i, j, k));
				}
			}
		}
	}
}
void MOCMesh::GetMocIndexByMaterial(std::vector< SMocIndex>& vSMocIndex, std::string strMaterial)
{
	SMocIndex sTemp;
	if (m_firstCreated)
	{
		for (int i = 0; i < m_vAssembly.size(); i++)
		{
			//m_vAssembly[i].v_field.resize(m_vAssembly[i].pAssembly_type->v_Cell.size());
			for (int j = 0; j < m_vAssembly[i].pAssembly_type->v_Cell.size(); j++)
			{
				//m_vAssembly[i].v_field[j].resize(64);
				for (int k = 0; k < m_vAssembly[i].pAssembly_type->v_Cell[j].vMeshPointPtrVec.size(); k++)
				{
					const MOCMeshPoint& mocPoint = dynamic_cast<const MOCMeshPoint&>
						(*m_vAssembly[i].pAssembly_type->v_Cell[j].vMeshPointPtrVec[k]);
					if(mocPoint.GetMaterialName() == strMaterial)
					{
						sTemp.iAssemblyIndex = i;
						sTemp.iCellIndex = j;
						sTemp.iMocIndex = k;
						vSMocIndex.push_back(sTemp);
					}
				}
			}
		}
	}
	else
	{
		for (int i = 0; i < m_vAssemblyField.size(); i++)
		{
			for (int j = 0; j < m_vAssemblyField[i].size(); j++)
			{
				for (int k = 0; k < m_vAssemblyField[i][j].size(); k++)
				{
					if (m_vAssemblyField[i][j][k]->GetMaterialName() == strMaterial)
						vSMocIndex.push_back(SMocIndex(i, j, k));
				}
			}
		}
	}


}

void MOCMesh::InitAssembly()
{
    //caculate the height of assebmly
	std::pair<int, Scalar> mocAxial = GetAxialInformation();
	double mocHeight = mocAxial.first * mocAxial.second;
	//caculate the leftdown and rightup corner of assembly on original coordinate from datafile
	for(int i=0;i<m_vAssemblyType.size();i++)
	{
		double xAssembly_Min = 10000, yAssembly_Min = 10000,xAssembly_Max = -10000, yAssembly_Max = -10000;
		for (int j=0;j < m_vAssemblyType[i].v_Cell.size(); j++)
		{
			double xMin = 10000, yMin = 10000, xMax = -10000, yMax = -10000, zMin=1000, zMax=-1000;
			for (int k = 0; k < m_vAssemblyType[i].v_Cell[j].vMeshPointPtrVec.size(); k++)
			{
				const MeshPoint& mocPoint = *m_vAssemblyType[i].v_Cell[j].vMeshPointPtrVec[k];

				for (int m = 0; m < mocPoint.VerticesNum(); m++)
				{
					xMin = min(xMin, mocPoint.VerticeCoordinate(m).x_);
					yMin = min(yMin, mocPoint.VerticeCoordinate(m).y_);
					xMax = max(xMax, mocPoint.VerticeCoordinate(m).x_);
					yMax = max(yMax, mocPoint.VerticeCoordinate(m).y_);
				}
			}
			// coordinate of leftdown and rightup corner of cell
			m_vAssemblyType[i].v_Cell[j].vCell_LeftDownPoint = Vector(xMin, yMin,0);
			m_vAssemblyType[i].v_Cell[j].vCell_RightUpPoint = Vector(xMax, yMax, mocHeight);

			xAssembly_Min = min(xAssembly_Min, xMin);
			yAssembly_Min = min(yAssembly_Min,yMin);
			xAssembly_Max = max(xAssembly_Max, xMax);
			yAssembly_Max = max(yAssembly_Max, yMax);
			if (m_bSingleCell)
			{
				m_vAssemblyType[i].xLength = xAssembly_Max - xAssembly_Min;
				m_vAssemblyType[i].yLength = yAssembly_Max - yAssembly_Min;
			}
		}
		m_vAssemblyType[i].vAssemblyType_LeftDownPoint = Vector(xAssembly_Min, yAssembly_Min, 0);
		m_vAssemblyType[i].vAssemblyType_RightUpPoint = Vector(xAssembly_Max, yAssembly_Max, mocHeight);
	}
	for (int i = 0; i < m_vAssembly.size(); i++)
	{
		m_vAssembly[i].pAssembly_type = GetAssemblyTypePointer(m_vAssembly[i].iAssemblyType);
	}
	//caculate the leftdown and rightup corner of assembly on system coordinate
	if (m_bSingleCell)
	{
			m_vAssembly[0].vAssembly_LeftDownPoint = Vector(0, 0, 0);
			double xLength = m_vAssembly[0].pAssembly_type->v_Cell[0].vCell_RightUpPoint.x_ -
			m_vAssembly[0].pAssembly_type->v_Cell[0].vCell_LeftDownPoint.x_;
			double yLength = m_vAssembly[0].pAssembly_type->v_Cell[0].vCell_RightUpPoint.y_ -
				m_vAssembly[0].pAssembly_type->v_Cell[0].vCell_LeftDownPoint.y_;
			m_vAssembly[0].vAssembly_RightUpPoint.x_ = xLength;
			m_vAssembly[0].vAssembly_RightUpPoint.y_ = yLength;
			m_vAssembly[0].vAssembly_RightUpPoint.z_ = mocHeight;
	}
	else
	{
		for (int xIndex = 0; xIndex < m_pAssemblyIndex->m_assemblyIndex.size(); xIndex++)
		{
			for (int yIndex = 0; yIndex < m_pAssemblyIndex->m_assemblyIndex[xIndex].size(); yIndex++)
			{
				int iAssemblyIndex = m_pAssemblyIndex->getAssemblyIndex(xIndex, yIndex);
				m_vAssembly[iAssemblyIndex].vAssembly_RightUpPoint.z_ = mocHeight;
				if (xIndex == 0 && yIndex == 0)
				{
					int iAssembly = m_pAssemblyIndex->getAssemblyIndex(0, 0);
					m_vAssembly[iAssembly].vAssembly_LeftDownPoint = Vector(0, 0, 0);
					m_vAssembly[iAssembly].vAssembly_RightUpPoint.x_ =m_vAssembly[iAssembly].pAssembly_type->xLength;
					m_vAssembly[iAssembly].vAssembly_RightUpPoint.y_ = m_vAssembly[iAssembly].pAssembly_type->yLength;
				}
				else if(xIndex == 0 && yIndex != 0)
				{
					int iDownPreIndex = m_pAssemblyIndex->getAssemblyIndex(xIndex, yIndex - 1);

					m_vAssembly[iAssemblyIndex].vAssembly_LeftDownPoint.x_ = 0;
					m_vAssembly[iAssemblyIndex].vAssembly_RightUpPoint.x_ = m_vAssembly[iAssemblyIndex].vAssembly_LeftDownPoint.x_ 
						+ m_vAssembly[iAssemblyIndex].pAssembly_type->xLength;

					m_vAssembly[iAssemblyIndex].vAssembly_LeftDownPoint.y_ = m_vAssembly[iDownPreIndex].vAssembly_LeftDownPoint.y_ +
						m_vAssembly[iDownPreIndex].pAssembly_type->yLength;
					m_vAssembly[iAssemblyIndex].vAssembly_RightUpPoint.y_ = m_vAssembly[iAssemblyIndex].vAssembly_LeftDownPoint.y_ + m_vAssembly[iAssemblyIndex].pAssembly_type->yLength;


				}
				else if(xIndex != 0 && yIndex == 0)
				{
					int iLeftPreIndex = m_pAssemblyIndex->getAssemblyIndex(xIndex - 1, yIndex);

					m_vAssembly[iAssemblyIndex].vAssembly_LeftDownPoint.y_ = 0;
					m_vAssembly[iAssemblyIndex].vAssembly_RightUpPoint.y_ = m_vAssembly[iAssemblyIndex].vAssembly_LeftDownPoint.y_ 
						+ m_vAssembly[iAssemblyIndex].pAssembly_type->yLength;

					m_vAssembly[iAssemblyIndex].vAssembly_LeftDownPoint.x_ = m_vAssembly[iLeftPreIndex].vAssembly_LeftDownPoint.x_ +
						m_vAssembly[iLeftPreIndex].pAssembly_type->xLength;
					m_vAssembly[iAssemblyIndex].vAssembly_RightUpPoint.x_ = m_vAssembly[iAssemblyIndex].vAssembly_LeftDownPoint.x_ + m_vAssembly[iAssemblyIndex].pAssembly_type->xLength;

				}
				else
				{
					int iLeftPreIndex = m_pAssemblyIndex->getAssemblyIndex(xIndex - 1, yIndex);
					int iDownPreIndex = m_pAssemblyIndex->getAssemblyIndex(xIndex, yIndex - 1);

					m_vAssembly[iAssemblyIndex].vAssembly_LeftDownPoint.x_ = m_vAssembly[iLeftPreIndex].vAssembly_LeftDownPoint.x_ +
						m_vAssembly[iLeftPreIndex].pAssembly_type->xLength;
					m_vAssembly[iAssemblyIndex].vAssembly_LeftDownPoint.y_ = m_vAssembly[iDownPreIndex].vAssembly_LeftDownPoint.y_ +
						m_vAssembly[iDownPreIndex].pAssembly_type->yLength;

					m_vAssembly[iAssemblyIndex].vAssembly_RightUpPoint.x_ = m_vAssembly[iAssemblyIndex].vAssembly_LeftDownPoint.x_ + m_vAssembly[iAssemblyIndex].pAssembly_type->xLength;
					m_vAssembly[iAssemblyIndex].vAssembly_RightUpPoint.y_ = m_vAssembly[iAssemblyIndex].vAssembly_LeftDownPoint.y_ + m_vAssembly[iAssemblyIndex].pAssembly_type->yLength;
				}
			}
		}
	}
	//sort meshpoint on pointid from apl file
	
	for (int i = 0;i < m_vAssemblyType.size();i++)
	{
		for (int j=0;j< m_vAssemblyType[i].v_Cell.size(); j++)
		{
			std::sort(m_vAssemblyType[i].v_Cell[j].vMeshPointPtrVec.begin(), m_vAssemblyType[i].v_Cell[j].vMeshPointPtrVec.end(), [](std::shared_ptr<MeshPoint> a, std::shared_ptr<MeshPoint> b) {
				return a->PointID() < b->PointID();
				});
		}
	}
	//init field
	for (int i = 0; i < m_vAssembly.size(); i++)
	{
		m_vAssembly[i].v_field.resize(m_vAssembly[i].pAssembly_type->v_Cell.size());
		for (int j = 0; j < m_vAssembly[i].pAssembly_type->v_Cell.size(); j++)
		{
			m_vAssembly[i].v_field[j].resize(m_vAssembly[i].pAssembly_type->v_Cell[j].vMeshPointPtrVec.size());
		}
	}

	
}

void MOCMesh::setMeshInformation(string line)
{
	stringstream stringline(line);
	string token;
	stringline >> token;
	layerMeshNum = stod(token);  //mesh node number
	stringline >> token;
	EdgeNum = stod(token);  //edge number
	stringline >> token;
	coarseMeshNum = stod(token);  //coarse mesh number
}

void MOCMesh::setAxialInformation(string line)
{
	stringstream stringline(line);
	string token;
	stringline >> token;
	stringline >> token;    //delete 0
	stringline >> token;
	axialNum = stod(token);  //layer number in axial direction
	for (int i = 0; i < axialNum; i++)
	{
		std::pair<int, double> temp0;
		stringline >> token;
		temp0.second  = stod(token);  //mesh height in axial direction  
		stringline >> token;
		temp0.first = stod(token);  // mesh number in axial direction(no use for this parameter )
		axialInformation.push_back(temp0);
	}
}


//set all edge objects
void MOCMesh::setEdgeInformation(string lineType, string linePosition, int edgeIDTemperary, std::vector<MOCEdge>& allEdges, int nFineMesh)
{
	vector<int> meshIDTemperary;
	stringstream stringline(lineType);
	string token;
	stringline >> token;
	int edgeTypeTemperary = stod(token);  //edge type
	int AditionArcEdgeID = 0;
	if (edgeTypeTemperary == 1)  //line
	{
		stringline >> token;
		meshIDTemperary.push_back(stod(token));  //node in the right
		stringline >> token;
		meshIDTemperary.push_back(stod(token));  //node in the left
		stringstream stringline0(linePosition);
		stringline0 >> token;
		double beginPoint_x = stod(token);  //coordinate of first point
		stringline0 >> token;
		double beginPoint_y = stod(token);
		stringline0 >> token;
		double endPoint_x = stod(token) + beginPoint_x;   //coordinate of second point
		stringline0 >> token;
		double endPoint_y = stod(token) + beginPoint_y;
		std::array<double, 3> beginPoint{ beginPoint_x, beginPoint_y, 0.0 };
		std::array<double, 3> endPoint{ endPoint_x, endPoint_y, 0.0 };
		MOCEdge edge0 = MOCEdge(beginPoint, endPoint, meshIDTemperary, edgeIDTemperary, edgeTypeTemperary);  //edge object
		edge0.arcCenter = Vector(0.0, 0.0, 0.0);
		edge0.arcAxisDir = Vector(0.0, 0.0, 0.0);
		allEdges.push_back(edge0);  //conserve edge object
	}

	if (edgeTypeTemperary == 2)//circle 
	{

	}
	if (edgeTypeTemperary == 3)//arc
	{
		stringline >> token;
		meshIDTemperary.push_back(stod(token));  //node in the right side
		stringline >> token;
		meshIDTemperary.push_back(stod(token));  //node in the left side

		double Center_point_X, Center_point_Y, Radius, Angle, Delta_angle;
		stringstream stringline0(linePosition);
		stringline0 >> token;
		Center_point_X = stod(token);//circle center coordinate
		stringline0 >> token;
		Center_point_Y = stod(token);

		stringline0 >> token;
		Radius = stod(token);//radius

		stringline0 >> token;
		Angle = stod(token);//angle of first point

		stringline0 >> token;
		Delta_angle = stod(token);//included angle
		if (Delta_angle < 0)
		{
			Delta_angle = Delta_angle + 360;
		}
		double fine_Delta_angle = Delta_angle / nFineMesh;
		for (int i = 1; i <= nFineMesh; i++)
		{
			double beginPoint_x = Center_point_X + Radius * cos(Angle * PI / 180.0);  //first point coordinate
			double beginPoint_y = Center_point_Y + Radius * sin(Angle * PI / 180.0);
			double endPoint_x = Center_point_X + Radius * cos((Angle + fine_Delta_angle) * PI / 180.0);  //first second coordinate
			double endPoint_y = Center_point_Y + Radius * sin((Angle + fine_Delta_angle) * PI / 180.0);
			std::array<double, 3> beginPoint{ beginPoint_x, beginPoint_y, 0.0 };
			std::array<double, 3> endPoint{ endPoint_x, endPoint_y, 0.0 };
			MOCEdge edge0 = MOCEdge(beginPoint, endPoint, meshIDTemperary, edgeIDTemperary, edgeTypeTemperary);  //edge object
			if (i > 1)
			{
				AditionArcEdgeID++;//为了保存多出来的弧边的ID
				edge0 = MOCEdge(beginPoint, endPoint, meshIDTemperary, edgeIDTemperary + EdgeNum + AditionArcEdgeID, edgeTypeTemperary);  //edge object
			}
			edge0.arcCenter = Vector(Center_point_X, Center_point_Y, 0.0);
			edge0.arcAxisDir = Vector(0.0, 0.0, 1.0);
			allEdges.push_back(edge0);  //conserve edge object
			Angle = Angle + fine_Delta_angle;
		}
	}

}
//set all surface objects
void MOCMesh::setMeshFaceInformation(vector<int> meshIDTransfer, vector<string> meshFaceTypeTransfer, vector<string> meshFaceTemperatureNameTransfer, std::vector<Surface>& allMeshFaces, std::vector<MOCEdge>& allEdges)
{
	for (int i = 0; i < layerMeshNum; i++)
	{
		Surface face0 = Surface(meshIDTransfer[i], meshIDTransfer[i], allEdges, meshFaceTypeTransfer[i], meshFaceTemperatureNameTransfer[i]);
		allMeshFaces.push_back(face0);
	}
}

void MOCMesh::ThreeDemMeshOutput(vector< shared_ptr<stringstream>>& vStreamTemperay, std::vector<Surface>& allMeshFaces, std::vector<std::string>& meshFaceTypeTransfer, int nFineMesh)
{
	string filename = "";
	double H0 = 0, H1 = 0;
	int index0 = 0;
	for (int j = 0; j < axialNum; j++)
	{
		if (j == 0)
		{
			H0 = 0;
			H1 = axialInformation[j].second;
		}
		else
		{
			H0 = H0 + axialInformation[j - 1].second;
			H1 = H1 + axialInformation[j].second;
		}
		for (int i = 0; i < layerMeshNum; i++)
		{
			vStreamTemperay[index0] = make_shared< stringstream>();

			std::stringstream& outFile= *vStreamTemperay[index0];// vStreamTemperay[index0];
			stringstream ssID;
			stringstream ssaxialID;
			string sID;
			string axialID;
			ssID << allMeshFaces[i].faceID;
			ssID >> sID;
			ssaxialID << (j + 1);
			ssaxialID >> axialID;
			filename = nFineMesh + "_poly" + meshFaceTypeTransfer[index0] + "_" + sID + "_" + axialID;
			index0++;

			//filename = filename + ".off";		
			filename = filename  + ".off";
			
			//fileNameTransfer.push_back(filename);
			//fstream outFile(filename);
			int pointNumPerMesh = allMeshFaces[i].facePointPosition.size() - 1;
			int ArcNumber = 0;
			for (int j = 0; j < allMeshFaces[i].faceEdges.size(); j++)
			{
				if (allMeshFaces[i].faceEdges[j].edgeType == 3)
				{
					ArcNumber++;
				}
			}
			outFile << "OFF" << endl;//off
			outFile << pointNumPerMesh * 2 << "\t" << pointNumPerMesh - (ArcNumber / nFineMesh) * (nFineMesh - 1) + 2 << "\t" << 0 << endl;;//point numbers,face number,0
			for (int j = 0; j < pointNumPerMesh; j++)  //coordinates of bottom points；
			{
				outFile << allMeshFaces[i].facePointPosition[j][0] << "\t" << allMeshFaces[i].facePointPosition[j][1] << "\t" << allMeshFaces[i].facePointPosition[j][2] + H0 << endl;
			}
			for (int j = 0; j < pointNumPerMesh; j++)  //coordinates of top points；
			{
				outFile << allMeshFaces[i].facePointPosition[j][0] << "\t" << allMeshFaces[i].facePointPosition[j][1] << "\t" << allMeshFaces[i].facePointPosition[j][2] + H1 << endl;
			}
			outFile << pointNumPerMesh << "\t";//number of bottom points
			for (int j = 0; j < pointNumPerMesh; j++)  //point id of bottom points;
			{
				outFile << pointNumPerMesh - 1 - allMeshFaces[i].facePointID[j] << "\t";//clockwise
			}
			allMeshFaces[i].curveInfo.push_back(0);
			allMeshFaces[i].curveFaceCenter.push_back(Vector(0.0, 0.0, 0.0));
			allMeshFaces[i].curveFaceAxisDir.push_back(Vector(0.0, 0.0, 0.0));
			outFile << endl;
			outFile << pointNumPerMesh << "\t";//number of top points
			for (int j = 0; j < pointNumPerMesh; j++)  //point id of top points;
			{
				outFile << allMeshFaces[i].facePointID[j] + pointNumPerMesh << "\t";//anticlockwise
			}
			allMeshFaces[i].curveInfo.push_back(0);
			allMeshFaces[i].curveFaceCenter.push_back(Vector(0.0, 0.0, 0.0));
			allMeshFaces[i].curveFaceAxisDir.push_back(Vector(0.0, 0.0, 0.0));
			outFile << endl;

			for (int j = 0; j < pointNumPerMesh; j++)  //side faces；
			{
				int presentEdgeTpye = allMeshFaces[i].faceEdges[j].edgeType;
				if (presentEdgeTpye == 1)//直线
				{
					outFile << 4 << "\t";//points number
					outFile << allMeshFaces[i].facePointID[j] << "\t";//anticlockwise
					outFile << allMeshFaces[i].facePointID[j + 1] << "\t";
					outFile << allMeshFaces[i].facePointID[j + 1] + pointNumPerMesh << "\t";
					outFile << allMeshFaces[i].facePointID[j] + pointNumPerMesh << "\t";
					outFile << endl;
					allMeshFaces[i].curveInfo.push_back(0);
					allMeshFaces[i].curveFaceCenter.push_back(Vector(0.0, 0.0, 0.0));
					allMeshFaces[i].curveFaceAxisDir.push_back(Vector(0.0, 0.0, 0.0));
				}
				else if (presentEdgeTpye == 2)//圆
				{

				}
				else       //圆弧
				{
					outFile << nFineMesh * 2 + 2 << "\t";//points number
					for (int k = 0; k <= nFineMesh; k++)
					{
						outFile << allMeshFaces[i].facePointID[j + k] << "\t";//anticlockwise
					}
					for (int k = nFineMesh; k >= 0; k--)
					{
						outFile << allMeshFaces[i].facePointID[j + k] + pointNumPerMesh << "\t";//anticlockwise
					}
					outFile << endl;
					j = j + nFineMesh - 1;//跳过多加进来的弧边
					allMeshFaces[i].curveInfo.push_back(1);
					allMeshFaces[i].curveFaceCenter.push_back(allMeshFaces[i].faceEdges[j].arcCenter);
					allMeshFaces[i].curveFaceAxisDir.push_back(allMeshFaces[i].faceEdges[j].arcAxisDir);
				}
			}
			//outFile.close();  //close iostream
			//vStreamTemperay.push_back(outFile);
		}
	}
}

void MOCMesh::OutputStatus(std::string outputFileName) const {
	std::ofstream ofs(outputFileName);

	ofs << m_preContext.str();
	ofs << std::endl << "\t\t*   material definition" << std::endl;
	std::unordered_set<std::string> materialNameSet;
	SMocIndex sMocIndex;
	for (auto sMocIndex : m_vSMocIndex) {
		//std::shared_ptr<MOCMeshPoint> p_mocMeshPoint = std::dynamic_pointer_cast<MOCMeshPoint>(meshPointPtr);
		std::string materialName =GetMaterialNameAtIndex(sMocIndex);
		const Medium* p_medium = nullptr;
		auto iter1 = m_mediumMap.find(materialName);
		if (iter1 != m_mediumMap.end())
			p_medium = &(iter1->second);
		
		boost::smatch m;
		if (!boost::regex_search(materialName, m, boost::regex(R"(_\d+)")) && materialName == "mMOD") {
			materialName = GetMaterialNameWithIDInField(sMocIndex);
			materialName += "_" + std::to_string(sMocIndex.iAssemblyIndex+1);
			//materialName += "_" + std::to_string(GetPointIDAtIndex(sMocIndex));
			if (!p_medium) {
				auto iter2 = m_mediumMap.find(materialName);
				if (iter2 != m_mediumMap.end())
					p_medium = &(iter2->second);
			}
		}
		
		if (!materialNameSet.count(materialName)) {
			ofs << std::endl;
			ofs << FormatStr("\t\t\t '%s' = MAT (  /", materialName.c_str());
			for (size_t i = 0; i < p_medium->eleFlagVec.size(); i++) {
				int eleFlag = p_medium->eleFlagVec[i];
				double eleDensity = p_medium->eleDensCalcFunVec[i](GetValueAtIndex(sMocIndex,ValueType::DENSITY));
				ofs << eleFlag << "," << FormatStr("%.6e", eleDensity);

				if (i != p_medium->eleFlagVec.size() - 1)
					ofs << ";" << std::endl << "\t\t\t\t\t";
				else
					ofs << ")" << std::endl;
			}
			materialNameSet.insert(materialName);
		}
	}

	ofs << std::endl << "\t\t*   temperature definition" << std::endl;
	for (auto sMocIndex : m_vSMocIndex) {
		ofs << std::endl;
		//std::shared_ptr<MOCMeshPoint> p_mocMeshPoint = std::dynamic_pointer_cast<MOCMeshPoint>(meshPointPtr);
		std::string tempName = GetTemperatureNameAtIndex(sMocIndex);// p_mocMeshPoint->GetTemperatureName();
		boost::smatch m;
		//if (!std::regex_search(tempName, m, std::regex(R"(_\d+)")))
			//tempName += "_" + std::to_string(GetPointIDAtIndex(sMocIndex));
		tempName +="_" + std::to_string(sMocIndex.iAssemblyIndex+1);
		ofs << FormatStr("\t\t\t '%s' = TEMP(%.6lf)", tempName.c_str(), GetValueAtIndex(sMocIndex,ValueType::TEMPERAURE)) << std::endl;
	}

	ofs << std::endl << m_sufContext.str();

	ofs.close();

	return;
}

void MOCMesh::InitMOCFromInputFile(std::string inputFileName) {
	std::ifstream ifs(inputFileName);
	if (!ifs.is_open()) 
	{
		Logger::LogError("cannot find the MOC data file:" + inputFileName);
		exit(EXIT_FAILURE);
	}
	m_mediumMap.clear();
	m_preContext.clear();
	m_sufContext.clear();

	std::unordered_map<std::string, double> materialDensityMap, temperatureMap;
	enum TOKEN {
		MATERIAL_DEFINITION,
		TEMPERATURE_DEFINITION,
		OTHER_DEFINITION
	} token = OTHER_DEFINITION;
	std::string detailDef;
	std::stringstream* p_currentContext = &m_preContext;
	while (!ifs.eof()) {
		std::string line;
		std::getline(ifs, line);
		TOKEN nextToken = token;
		if (boost::regex_search(line, boost::regex(R"(\*\s*material\s+definition)")))
			nextToken = MATERIAL_DEFINITION;
		if (boost::regex_search(line, boost::regex(R"(\*\s*temperature\s+definition)")))
			nextToken = TEMPERATURE_DEFINITION;
		if (boost::regex_search(line, boost::regex(R"(MODULE)")))
			nextToken = OTHER_DEFINITION;

		if (nextToken != token) {
			p_currentContext = &m_sufContext;
			if (token == MATERIAL_DEFINITION) {
				boost::smatch outer_match;
				while (boost::regex_search(detailDef, outer_match, boost::regex(R"('[\w\.]+')"))) {
					std::string materialName = outer_match.str(0).substr(1, outer_match.str(0).length() - 2);
					detailDef = outer_match.suffix().str();

					boost::regex_search(detailDef, outer_match, boost::regex(R"(MAT[^\)]+\))"));
					std::string matInfo = outer_match.str(0);
					boost::smatch innerMatch;
					Medium& medium = m_mediumMap[materialName];
					while (boost::regex_search(matInfo, innerMatch, boost::regex(R"(([\+|-]?\d+(\.{0}|\.\d+))[Ee]?([\+|-]?\d+))"))) {
						int eleFlag = std::stoi(innerMatch.str(0));
						medium.eleFlagVec.push_back(eleFlag);

						matInfo = innerMatch.suffix().str();
						boost::regex_search(matInfo, innerMatch, boost::regex(R"(([\+|-]?\d+(\.{0}|\.\d+))[Ee]?([\+|-]?\d+))"));
						double eleDens = std::stod(innerMatch.str(0));
						medium.eleDensCalcFunVec.push_back([eleDens](double) {return eleDens; });

						matInfo = innerMatch.suffix().str();
					}

					if (materialName.find("mMOD") != materialName.npos) {
						double slackH2ODensity = 1e3;
						for (size_t i = 0; i < medium.eleFlagVec.size(); i++) {
							int eleFlag = medium.eleFlagVec[i];
							if (eleFlag == 1001) {
								double eleDensity = medium.eleDensCalcFunVec[i](0.0);
								slackH2ODensity = eleDensity / (NA * BARN * 2) * 18000.0;
								medium.eleDensCalcFunVec[i] = [](double density) {
									return density / 18000.0 * NA * BARN * 2;
								};
							}
							else if (eleFlag == 8016) {
								medium.eleDensCalcFunVec[i] = [](double density) {
									return density / 18000.0 * NA * BARN;
								};
							}
						}
						for (size_t i = 0; i < medium.eleFlagVec.size(); i++) {
							int eleFlag = medium.eleFlagVec[i];
							if (eleFlag == 5000) {
								double eleDensity = medium.eleDensCalcFunVec[i](0.0);
								medium.eleDensCalcFunVec[i] = [eleDensity, slackH2ODensity](double density) {
									return density / slackH2ODensity * eleDensity;
								};
							}
						}
						materialDensityMap[materialName] = slackH2ODensity;
					}
					else
						materialDensityMap[materialName] = 0.0;

					detailDef = outer_match.suffix().str();
				}
			}
			else if (token == TEMPERATURE_DEFINITION) {
				boost::smatch outer_match;
				while (boost::regex_search(detailDef, outer_match, boost::regex(R"('[\w\.]+')"))) {
					std::string tempName = outer_match.str(0).substr(1, outer_match.str(0).length() - 2);
					detailDef = outer_match.suffix().str();

					boost::regex_search(detailDef, outer_match, boost::regex(R"(TEMP[^\)]+\))"));
					std::string tempInfo = outer_match.str(0);
					boost::smatch innerMatch;
					boost::regex_search(tempInfo, innerMatch, boost::regex(R"(([\+|-]?\d+(\.{0}|\.\d+))[Ee]?([\+|-]?\d+))"));
					double tempValue = std::stod(innerMatch.str(0));

					temperatureMap[tempName] = tempValue;
					detailDef = outer_match.suffix().str();
				}
			}
			
			detailDef.clear();
			token = nextToken;
		}

		if (token == MATERIAL_DEFINITION ||
			token == TEMPERATURE_DEFINITION) {
			detailDef += line;
		}
		else
			*p_currentContext << line << std::endl;
	}

	ifs.close();
	for (auto sMocIndex : m_vSMocIndex)
	{
		//std::shared_ptr<MOCMeshPoint> p_mocMeshPoint = std::dynamic_pointer_cast<MOCMeshPoint>(GetMocMeshPointPtr(sMocIndex));
		{
			boost::smatch m;
			std::string materialName = GetMaterialNameAtIndex(sMocIndex);
			auto iter1 = materialDensityMap.find(materialName);
			if (iter1 != materialDensityMap.end())
				SetValueAtIndex(sMocIndex, iter1->second, ValueType::DENSITY);
			//p_mocMeshPoint->SetValue(iter1->second, ValueType::DENSITY);

			if (boost::regex_search(materialName, m, boost::regex(R"(_\d+)"))) {
				std::string materialMetaName = m.prefix().str();
				auto iter2 = materialDensityMap.find(materialMetaName);
				if (iter2 != materialDensityMap.end()) {
					//p_mocMeshPoint->SetValue(iter2->second, ValueType::DENSITY);
					SetValueAtIndex(sMocIndex, iter2->second, ValueType::DENSITY);
				}
			}
		}
		{
			boost::smatch m;
			std::string temperatureName = GetTemperatureNameAtIndex(sMocIndex);
			auto iter1 = temperatureMap.find(temperatureName);
			if (iter1 != temperatureMap.end())
				//p_mocMeshPoint->SetValue(iter1->second, ValueType::TEMPERAURE);
				SetValueAtIndex(sMocIndex, iter1->second, ValueType::TEMPERAURE);

			if (boost::regex_search(temperatureName, m, boost::regex(R"(_\d+)"))) {
				std::string tempMetaName = m.prefix().str();
				auto iter2 = temperatureMap.find(tempMetaName);
				if (iter2 != temperatureMap.end()) {
					//p_mocMeshPoint->SetValue(iter2->second, ValueType::TEMPERAURE);
					SetValueAtIndex(sMocIndex, iter2->second, ValueType::TEMPERAURE);
				}
			}

		}
	}


	return;
}

void MOCMesh::InitMOCHeatPower(std::string heatPowerFileName) 
{
	std::ifstream ifs(heatPowerFileName);
	if (!ifs.is_open())
	{
		Logger::LogError("cannot find the moc data file:" + heatPowerFileName);
		exit(EXIT_FAILURE);
	}
	std::vector<std::vector<double>> vHeatPowerValue;
	//vHeatPowerValue.resize(m_vAssemblyField.size());
	
	string line;
	while (getline(ifs, line))  //read mesh data
	{
	loop_:
		int iPos = line.find("_");
		if (iPos != std::string::npos)
		{
			int iAssembly_index = stod(line.substr(0, iPos));
			int iNumb_Value = stod(line.substr(iPos + 1, line.length()));
			std::vector<double> vValue;
			//vValue.resize(iNumb_Value);

			string tokenMeshId = "";
			int out0 = 1;
			std::streampos pos;
			while (getline(ifs, line))
			{
				//pos = infile.tellg();
				//getline(ifs, line);
				//outFile << line << endl;
				stringstream stringlineMeshID(line);
				while (stringlineMeshID >> tokenMeshId)
				{
					if (tokenMeshId.find("_") != std::string::npos)
					{
						vHeatPowerValue.push_back(vValue);
						//out0 = 0;
						//stringlineMeshID >> token;
						goto loop_;
						//infile.seekg(pos);
						break;
					}
					vValue.push_back(stod(tokenMeshId));
				}
				//out0 = 0;
			}
			vHeatPowerValue.push_back(vValue);

			/*
			for (int i = 0; i < iNumb_Value;i++)
			{

				getline(ifs, line);
				vValue[i] = stod(line);
			}*/
			
		}
		
	}
	/*
	if (powerInput.size() != m_vSMocIndex.size())
	{
		Logger::LogError("Wrong number in heatpower.txt");
		exit(EXIT_FAILURE);
	}
	*/
	int kk = 0;
	for (int i = 0; i < m_vAssemblyField.size(); i++)
	{
		for (int j = 0; j < m_vAssemblyField[i].size(); j++)
		{
			for (int k = 0; k < m_vAssemblyField[i][j].size(); k++)
			{
				if (m_vAssemblyField[i][j][k])
				{
					double value = vHeatPowerValue[i][m_vAssemblyField[i][j][k]->m_iPointID - 1];
					//SetValueAtIndex(SMocIndex(i, j, k), powerInput[kk++], ValueType::HEATPOWER);
					SetValueAtIndex(SMocIndex(i, j, k), value, ValueType::HEATPOWER);
				}
			}
		}
	}
	ifs.close();	
}

void MOCMesh::WriteTecplotFile
(
	std::string  strFilename,
	std::string strMaterialType
)
{
	std::ofstream ofile(strFilename);
	ofile << "TITLE =\"" << "polyhedron" << "\"" << endl;
	ofile << "VARIABLES = " << "\"x\"," << "\"y\"," << "\"z\"" << endl;
	for (int i = 0; i < m_vAssembly.size(); i++)
	{
		//if (i == 0)
			//continue;
		//if (i > 0)
			//break;

		//平移坐标,
		double x = m_vAssembly[i].pAssembly_type->vAssemblyType_LeftDownPoint.x_-m_vAssembly[i].vAssembly_LeftDownPoint.x_;
		double y= m_vAssembly[i].pAssembly_type->vAssemblyType_LeftDownPoint.y_-m_vAssembly[i].vAssembly_LeftDownPoint.y_;
		//std::cout << "i="<<i << " x=" << x << std::endl;
		//std::cout << "i="<<i<< " y=" << y << std::endl;
		for (int j = 0; j < m_vAssembly[i].pAssembly_type->v_Cell.size(); j++)
		{
			/*
			std::vector<int> bVect(299, -1);
			std::vector<int> bIndex{ 17, 34, 51, 18, 35, 52, 19, 36, 53 };
			for (auto a : bIndex)
			{
				bVect[a] = 1;
			}
			if (bVect[j] == -1)
				continue;
			*/
			for (int k = 0; k < m_vAssembly[i].pAssembly_type->v_Cell[j].vMeshPointPtrVec.size(); k++)
			{
				const MHTMocMeshPoint& mhtPolyhedron = dynamic_cast<const MHTMocMeshPoint&>
					(*m_vAssembly[i].pAssembly_type->v_Cell[j].vMeshPointPtrVec[k]);
				const MOCMeshPoint& mocPoint = dynamic_cast<const MOCMeshPoint&>
					(*m_vAssembly[i].pAssembly_type->v_Cell[j].vMeshPointPtrVec[k]);
				if(strMaterialType!="" && strMaterialType != mhtPolyhedron.GetMaterialName())continue;
				//move to system coordinate
				
				MHTMocMeshPoint meshPoint = mhtPolyhedron;
				Vector vPoint(-x, -y, 0);
				meshPoint.Move(vPoint);
				meshPoint.WriteTecplotZones(ofile);
			}
		}
	}
	
	/*
	for (int i = 0; i < this->m_meshPointPtrVec.size(); i++)
	{
		const MHTMeshPoint& mhtPolyhedron = dynamic_cast<const MHTMeshPoint&>(*m_meshPointPtrVec[i]);
		const MOCMeshPoint& mocPoint = dynamic_cast<const MOCMeshPoint&>(*m_meshPointPtrVec[i]);
		if (mType != mocPoint.GetMaterialName()) continue;
		mhtPolyhedron.WriteTecplotZones(ofile);
	}
	*/
	ofile.close();
	return;
}

void MOCMesh::WriteSurfaceTecplotFile
(
	std::string fileName
)
{
	//step 1: calculate boundaries
	std::vector<MHT::Polygon> polygonList;
	Vector initialCenter = 0.5 * (m_vAssembly[0].vAssembly_LeftDownPoint + m_vAssembly[0].vAssembly_RightUpPoint);
	Scalar xmin = initialCenter.x_;
	Scalar xmax = initialCenter.x_;
	Scalar ymin = initialCenter.y_;
	Scalar ymax = initialCenter.y_;
	Scalar zmin = initialCenter.z_;
	Scalar zmax = initialCenter.z_;
	for (size_t i = 0; i < 1; i++)
	{
		Vector assemblyMin = m_vAssembly[i].vAssembly_LeftDownPoint;
		xmin = std::min(xmin, assemblyMin.x_);
		ymin = std::min(ymin, assemblyMin.y_);
		zmin = std::min(zmin, assemblyMin.z_);
		Vector assemblyMax = m_vAssembly[i].vAssembly_RightUpPoint;
		xmax = std::max(xmax, assemblyMax.x_);
		ymax = std::max(ymax, assemblyMax.y_);
		zmax = std::max(zmax, assemblyMax.z_);
	}
	Vector globalMin(xmin, ymin, zmin);
	Vector globalMax(xmax, ymax, zmax);
	//step 2: collect faces on assembly boundaries
	for (size_t i = 0; i < m_vAssembly.size(); i++)
	{
		Vector assemblyMin = m_vAssembly[i].vAssembly_LeftDownPoint;
		for (int j = 0; j < m_vAssembly[i].pAssembly_type->v_Cell.size(); j++)
		{
			for (int k = 0; k < m_vAssembly[i].pAssembly_type->v_Cell[j].vMeshPointPtrVec.size(); k++)
			{
				const MHTMocMeshPoint& mhtPolyhedron = dynamic_cast<const MHTMocMeshPoint&>
					(*m_vAssembly[i].pAssembly_type->v_Cell[j].vMeshPointPtrVec[k]);
				MHTMocMeshPoint meshPoint = mhtPolyhedron;
				meshPoint.Move(assemblyMin);
				std::vector<MHT::Polygon> subPolygonList = meshPoint.GetFacesOnBoxBoundary(globalMin, globalMax, 1e-6);
				if (0 != subPolygonList.size())
				{
					polygonList.insert(polygonList.end(), subPolygonList.begin(), subPolygonList.end());
				}
			}
		}
	}

	//step 3: plot faces collected in step 2
	int nodeNum = 0;
	int faceNum = 0;
	int elemNum = (int)polygonList.size();

	for (int polyID = 0; polyID < (int)polygonList.size(); polyID++)
	{
		nodeNum += (int)polygonList[polyID].v_point.size();
		faceNum += (int)polygonList[polyID].v_point.size();
	}
	std::ofstream outFile(fileName);
	outFile << "TITLE = \"PolygonList\" " << std::endl;
	outFile << "VARIABLES = " << "\"x \"," << "\"y \"," << "\"z \"" << std::endl;
	outFile << "ZONE T = \"polygons\"" << std::endl;
	outFile << "STRANDID = 0, SOLUTIONTIME = 0" << std::endl;
	outFile << "Nodes = " << nodeNum << ", Faces = " << faceNum << ", Elements = " << elemNum;
	outFile << ", ZONETYPE = FEPolygon" << std::endl;
	outFile << "NumConnectedBoundaryFaces = 0, TotalNumBoundaryConnections = 0" << std::endl;
	outFile << "DT = (SINGLE SINGLE SINGLE)";
	// Write node X
	int nCountI(0);
	for (int polyID = 0; polyID < (int)polygonList.size(); polyID++)
	{
		for (int i = 0; i < (int)polygonList[polyID].v_point.size(); i++, nCountI++)
		{
			if (nCountI % 5 == 0) outFile << std::endl;
			outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << polygonList[polyID].v_point[i].x_ << " ";
		}
	}
	// Write node Y
	int nCountJ(0);
	for (int polyID = 0; polyID < (int)polygonList.size(); polyID++)
	{
		for (int i = 0; i < (int)polygonList[polyID].v_point.size(); i++, nCountJ++)
		{
			if (nCountJ % 5 == 0) outFile << std::endl;
			outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << polygonList[polyID].v_point[i].y_ << " ";
		}
	}
	// Write node Z
	int nCountK(0);
	for (int polyID = 0; polyID < (int)polygonList.size(); polyID++)
	{
		for (int i = 0; i < (int)polygonList[polyID].v_point.size(); i++, nCountK++)
		{
			if (nCountK % 5 == 0) outFile << std::endl;
			outFile << std::setprecision(8) << std::setiosflags(std::ios::scientific) << polygonList[polyID].v_point[i].z_ << " ";
		}
	}
	outFile << std::endl;
	outFile << "# Face nodes" << std::endl;
	int nLineCount(1);
	for (int polyID = 0; polyID < (int)polygonList.size(); polyID++)
	{
		for (int i = 0; i < (int)polygonList[polyID].v_point.size() - 1; i++)
		{
			outFile << nLineCount + i << " " << nLineCount + i + 1 << " ";
		}
		outFile << nLineCount + polygonList[polyID].v_point.size() - 1 << " " << nLineCount << " ";
		outFile << std::endl;
		nLineCount += (int)polygonList[polyID].v_point.size();
	}

	outFile << "# left elements";
	int count = 0;
	int nElemCount(1);
	for (int polyID = 0; polyID < (int)polygonList.size(); polyID++)
	{
		for (int i = 0; i < (int)polygonList[polyID].v_point.size(); i++)
		{
			if (count % 10 == 0) outFile << std::endl;
			count++;
			outFile << nElemCount << " ";
		}
		nElemCount++;
	}
	outFile << std::endl;

	outFile << "# right elements";
	int nLeftElemCount(0);
	for (int polyID = 0; polyID < (int)polygonList.size(); polyID++)
	{
		for (int i = 0; i < (int)polygonList[polyID].v_point.size(); i++, nLeftElemCount++)
		{
			if (nLeftElemCount % 10 == 0) outFile << std::endl;
			outFile << "0" << " ";
		}
	}
	outFile << std::endl;
	outFile.close();

	return;
}

std::pair<int,Scalar> MOCMesh::GetAxialInformation()
{
	Scalar totalHeight = 0.0;
	int cellNum = this->axialInformation.size();
	for (int i = 0;i < cellNum;i++)
	{
		totalHeight += this->axialInformation[i].second;
	}
	return std::pair<int, Scalar>(cellNum, totalHeight / Scalar(cellNum));
}

void MOCMesh::WriteHeatPowerTxtFile()
{
	std::ofstream ofile("heatPower.txt");
	for (int i = 0; i < this->m_meshPointPtrVec.size(); i++)
	{
		const MHTMeshPoint& mhtPolyhedron = dynamic_cast<const MHTMeshPoint&>(*m_meshPointPtrVec[i]);
		const MOCMeshPoint& mocPoint = dynamic_cast<const MOCMeshPoint&>(*m_meshPointPtrVec[i]);
		Vector PointCenter = mocPoint.Center();
		Vector axisCenter(0.63, 0.63, 0.0);
		Vector axisNorm(0.0, 0.0, 1.0);
		Vector OP = PointCenter - axisCenter;
		Vector axicialLocation = (OP & axisNorm) * axisNorm;
		Scalar radius = (OP - axicialLocation).Mag();
		Scalar sigma = 1.0;
		Scalar hp = (1000.0 / sqrt(2.0 * PI) / sigma) * exp(-radius * radius / 2 / pow(sigma, 2));
		Scalar z = PointCenter.z_;
		hp = hp * z * (5.0 - z) / 6.25;
		ofile << hp << std::endl;
	}
	ofile.close();
	return;
}

void MOCMesh::Display()
{
	std::cout << "m_coarseMeshInfo:" << std::endl;
	for (int i = 0; i < m_coarseMeshInfo.size(); i++)
	{
		std::cout << "i = " << i << std::endl;
		for (int j = 0; j < m_coarseMeshInfo[i].size(); i++)
		{
			std::cout << m_coarseMeshInfo[i][j] << std::endl;
		}
	}
	std::cout << "layerMeshNum = " << layerMeshNum << std::endl;
	std::cout << "EdgeNum = " << EdgeNum << std::endl;
	std::cout << "coarseMeshNum = " << coarseMeshNum << std::endl;
	std::cout << "axialNum = " << axialNum << std::endl;
	std::cout << "axialInformation:" << std::endl;
	for (int i = 0; i < axialInformation.size(); i++)
	{
		std::cout << "i = " << i << std::endl;
		std::cout << axialInformation[i].first << "," << axialInformation[i].second << std::endl;
	}
}

Assembly_Type* MOCMesh::GetAssemblyTypePointer(int iAssemblyType)
{
	for (int i = 0; i < m_vAssemblyType.size(); i++)
	{
		if (iAssemblyType == m_vAssemblyType[i].iAssemblyType)
			return &m_vAssemblyType[i];
	}
	return nullptr;
}
SMocIndex MOCMesh::getIndex(Vector vPoint)
{
	return m_pAssemblyIndex->getIndex(vPoint);
}

Surface::Surface()
{
	facePointPosition.clear();
	facePointID.clear();
	faceID = 0;
	faceEdges.clear();
	faceType = "";
}

Surface::Surface(int faceID0, int nodeID, vector<MOCEdge> allEdgesTransfer, string meshFaceTypeTransfer, string meshTemperatureNameTransfer)
{
	int edgeNum0 = 0;
	faceID = faceID0;
	edgeNum0 = allEdgesTransfer.size();
	faceType = meshFaceTypeTransfer;
	face_temperatureName = meshTemperatureNameTransfer;
	for (int i = 0; i < edgeNum0; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			if (allEdgesTransfer[i].sideMeshID[j] == nodeID)
			{
				faceEdges.push_back(allEdgesTransfer[i]);
			}
		}
	}
	faceEdgeOrder(nodeID);
}

void Surface::faceEdgeOrder(int nodeID)  //save the edge in the anticlockwise
{
	int face_edge_num = faceEdges.size();
	MOCEdge mostLeftEdge;
	double minx = 1.0e8;
	double RoundingError = 1.0e-12;  //Rounding error of computer
	vector <MOCEdge> edgeTemperary;
	for (int i = 0; i < face_edge_num; i++)  //the most left edge and the edge type is line
	{
		for (int j = 0; j < 2; j++)
		{
			if (faceEdges[i].edgePoints[j][0] < minx && faceEdges[i].edgeType == 1)
			{
				minx = faceEdges[i].edgePoints[j][0];
				mostLeftEdge = faceEdges[i];
			}
		}
	}
	edgeTemperary.push_back(mostLeftEdge);
	vector<MOCEdge>faceEdgesTemporary;    //temporary variables for edge object
	for (int i = 0; i < faceEdges.size(); i++)
	{
		if (mostLeftEdge.edgeID == faceEdges[i].edgeID)continue;
		faceEdgesTemporary.push_back(faceEdges[i]);
	}
	std::array<double, 3> connectPoint{ 0.0, 0.0, 0.0 };    //connection point for two edge
	int pointOrderIndex = 0;
	//If the target nodes in the left side of the edge, the direction is anticlockwise
	if (mostLeftEdge.sideMeshID[1] == nodeID)
	{
		facePointID.push_back(pointOrderIndex);
		pointOrderIndex++;
		facePointPosition.push_back(mostLeftEdge.edgePoints[0]);
		facePointID.push_back(pointOrderIndex);
		pointOrderIndex++;
		facePointPosition.push_back(mostLeftEdge.edgePoints[1]);
		connectPoint = mostLeftEdge.edgePoints[1];
	}
	else
	{
		facePointID.push_back(pointOrderIndex);
		pointOrderIndex++;
		facePointPosition.push_back(mostLeftEdge.edgePoints[1]);
		facePointID.push_back(pointOrderIndex);
		pointOrderIndex++;
		facePointPosition.push_back(mostLeftEdge.edgePoints[0]);
		connectPoint = mostLeftEdge.edgePoints[0];
	}
	//Sort other edges and points anticlockwise

	double RoundingError0 = 1.0e-6;
	MOCEdge presentEdge = mostLeftEdge;
	for (int i = 1; i < face_edge_num; i++)
	{
		for (int j = 0; j < faceEdgesTemporary.size(); j++)
		{
			if (faceEdgesTemporary[j].edgeID == presentEdge.edgeID)
			{
				continue;
			}
			if (fabs(faceEdgesTemporary[j].edgePoints[0][0] - connectPoint[0]) < RoundingError0 && fabs(faceEdgesTemporary[j].edgePoints[0][1] - connectPoint[1]) < RoundingError0)
			{
				presentEdge = faceEdgesTemporary[j];
				connectPoint = faceEdgesTemporary[j].edgePoints[1];
				if (i == face_edge_num - 1)   //If it's the last edge, connect to the first point
				{
					facePointID.push_back(facePointID[0]);
				}
				else
				{
					facePointID.push_back(pointOrderIndex);
				}
				pointOrderIndex++;
				facePointPosition.push_back(faceEdgesTemporary[j].edgePoints[1]);
				edgeTemperary.push_back(faceEdgesTemporary[j]);
				faceEdgesTemporary.erase(faceEdgesTemporary.begin() + j);
				break;
			}
			if (fabs(faceEdgesTemporary[j].edgePoints[1][0] - connectPoint[0]) < RoundingError0 && fabs(faceEdgesTemporary[j].edgePoints[1][1] - connectPoint[1]) < RoundingError0)
			{
				presentEdge = faceEdgesTemporary[j];
				connectPoint = faceEdgesTemporary[j].edgePoints[0];
				if (i == face_edge_num - 1)
				{
					facePointID.push_back(facePointID[0]);
				}
				else
				{
					facePointID.push_back(pointOrderIndex);
				}
				pointOrderIndex++;
				facePointPosition.push_back(faceEdgesTemporary[j].edgePoints[0]);
				edgeTemperary.push_back(faceEdgesTemporary[j]);
				faceEdgesTemporary.erase(faceEdgesTemporary.begin() + j);
				break;
			}
		}
	}
	faceEdgesTemporary.clear();

	for (int i = 0; i < face_edge_num; i++)   //������ʱ�뷽���ÿ������ı߽�����������
	{
		faceEdges[i] = edgeTemperary[i];
	}
}
MOCEdge::MOCEdge()
{
	edgePoints.clear();
	sideMeshID.clear();
	edgeID = 0;
	edgeType = 0;
}

MOCEdge::MOCEdge(std::array<double, 3> beginPoint, std::array<double, 3> endPoint, vector<int> meshIDTransfer, int edgeIDTransfer, int edgeTypeTransfer)
{
	int sideMeshIDnum = 0;
	edgeID = edgeIDTransfer;
	edgePoints.push_back(beginPoint);
	edgePoints.push_back(endPoint);
	edgeType = edgeTypeTransfer;
	sideMeshIDnum = meshIDTransfer.size();
	for (int i = 0; i < sideMeshIDnum; i++)
	{
		sideMeshID.push_back(meshIDTransfer[i]);
	}
}

void MOCMesh::readMapFile(const std::vector<std::string>& materialList)
{
	struct IndexAndVaule
	{
		int iCFD;
		SMocIndex iMoc;
		double volume;
		std::string strMaterailName;
		std::string strTemptureName;
		int iPointID;
	};
	std::vector<IndexAndVaule> vIndexValue;
	int iMax_iAssembly = 0, iMax_iCell = 0, iMax_iMoc = 0;

	for (int kk = 0; kk < materialList.size(); kk++)
	{
		std::string fileName = "MapFile_" + materialList[kk] + "_CFDtoMOC";
		ifstream infile(fileName);
		if (!infile.is_open())
		{
			Logger::LogError("cannot find the CFD to MOC map file: " + fileName);
			exit(EXIT_FAILURE);
			return;
		}
		Logger::LogInfo("reading CFD to MOC map file in material: " + materialList[kk]);

		std::string line;

		while (getline(infile, line))
		{
			int i, j, k, m,iPointID;
			double value, valueVolume;
			std::string valueTemptureName, strMaterial;
			stringstream stringline(line);
			stringline >> i >> j >> k >> m >> value>>valueVolume>> valueTemptureName >> iPointID>> strMaterial;
			IndexAndVaule sTemp;
			sTemp.iCFD = i;
			sTemp.iMoc.iAssemblyIndex = j;
			sTemp.iMoc.iCellIndex = k;
			sTemp.iMoc.iMocIndex = m;
			sTemp.volume = valueVolume;
			
			sTemp.strTemptureName = valueTemptureName;
			sTemp.iPointID = iPointID;
			sTemp.strMaterailName = strMaterial;
			vIndexValue.push_back(sTemp);;
			//m_CFD_MOC_Map[i][sTemp] = value;
			iMax_iAssembly = max(iMax_iAssembly, j);
			iMax_iCell = max(iMax_iCell, k);
			iMax_iMoc = max(iMax_iMoc, m);
		}
		infile.close();
	}

	m_vAssemblyField.resize(iMax_iAssembly+1);
	for (int i = 0; i <= iMax_iAssembly; i++)
	{
		m_vAssemblyField[i].resize(iMax_iCell+1);
		for (int j = 0; j <= iMax_iCell; j ++ )
		{
			m_vAssemblyField[i][j].resize(iMax_iMoc+1);
		}
	}

	for (int i = 0; i < vIndexValue.size(); i++)
	{
		m_vAssemblyField[vIndexValue[i].iMoc.iAssemblyIndex][vIndexValue[i].
			iMoc.iCellIndex][vIndexValue[i].iMoc.iMocIndex] = std::make_shared<MocMeshField>();
		m_vAssemblyField[vIndexValue[i].iMoc.iAssemblyIndex][vIndexValue[i].
			iMoc.iCellIndex][vIndexValue[i].iMoc.iMocIndex]->SetVolume(vIndexValue[i].volume);
		m_vAssemblyField[vIndexValue[i].iMoc.iAssemblyIndex][vIndexValue[i].
			iMoc.iCellIndex][vIndexValue[i].iMoc.iMocIndex]->SetMaterialName(vIndexValue[i].strMaterailName);
		m_vAssemblyField[vIndexValue[i].iMoc.iAssemblyIndex][vIndexValue[i].
			iMoc.iCellIndex][vIndexValue[i].iMoc.iMocIndex]->m_iPointID = vIndexValue[i].iPointID;
		m_vAssemblyField[vIndexValue[i].iMoc.iAssemblyIndex][vIndexValue[i].
			iMoc.iCellIndex][vIndexValue[i].iMoc.iMocIndex]->m_temperatureName = vIndexValue[i].strTemptureName;


	}
	return;
}





