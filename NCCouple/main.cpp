#include <iostream>
#include <algorithm>
#include "MOCMesh.h"
#include "MOCIndex.h"
#include "CFDMesh.h"
#include "Solver.h"
#include "Logger.h"
#include "MHT_polyhedron/PolyhedronSet.h"
#include<time.h>
#include "./FileConvertor.h"
#include "./ConfigurationFile.h"

#include <sstream>
#include <mpi.h>
#ifdef _WIN32
#include <process.h>
#include <string>
#else
#include <unistd.h>
#include <string>
#include <string.h>
#endif
#include "./MHT_common/SystemControl.h"
#include "./MHT_mesh/UnGridFactory.h"
#include "./MHT_mesh/RegionConnection.h"
#include "./MHT_mesh/Mesh.h"
#include "./MHT_field/Field.h"
#include "./MHT_IO/FieldIO.h"



#include "./MHT_IO/VTKIO.h"
int g_iMpiID = -1;
int g_iNumProcs = 0;

//this example was designed for test of
//(1) rewriting a apl file
//(2) reading and writting of inp files

void WriteTotecplot( MOCMesh& mocmesh, CFDMesh& cfdmesh,
	const std::vector<int> &v_iCFD,const std::vector<SMocIndex>& vSMocindex,std::string fileName)
{
	std::string  mType;
	std::ofstream ofile(fileName);
	ofile << "TITLE =\"" << "polyhedron" << "\"" << endl;
	ofile << "VARIABLES = " << "\"x\"," << "\"y\"," << "\"z\"" << endl;
	for (int i = 0; i < vSMocindex.size(); i++)
	{

		//move coordinate
		/*
		int iAssembly = vSMocindex[i].iAssemblyIndex;
		double x = mocmesh.m_vAssembly[iAssembly].pAssembly_type->vAssemblyType_LeftDownPoint.x_- mocmesh.m_vAssembly[iAssembly].vAssembly_LeftDownPoint.x_;
		double y= mocmesh.m_vAssembly[iAssembly].pAssembly_type->vAssemblyType_LeftDownPoint.y_-mocmesh.m_vAssembly[iAssembly].vAssembly_LeftDownPoint.y_;
		
		const MHTMocMeshPoint& mocPoint = dynamic_cast<const MHTMocMeshPoint&>(*mocmesh.GetMocMeshPointPtr(vSMocindex[i]));
		*/
		shared_ptr<MHTMocMeshPoint> ptrMocMesh;
		mocmesh.MoveMeshToSysCoordinate(ptrMocMesh, vSMocindex[i]);
		MHTMocMeshPoint& mocPoint = *ptrMocMesh;

		//MHTMocMeshPoint meshPoint = mocPoint;
		//Vector vPoint(-x, -y, 0);
		//meshPoint.Move(vPoint);
		mocPoint.WriteTecplotZones(ofile);
	}
	for (int i = 0; i < v_iCFD.size(); i++)
	{
		const MHTMeshPoint& mhtPolyhedron = dynamic_cast<const MHTMeshPoint&>(*cfdmesh.GetMeshPointPtr(v_iCFD[i]));
		mhtPolyhedron.WriteTecplotZones(ofile);
	}

	ofile.close();
	return;
}

void MOC_APL_INP_FileTest()
{
	WarningContinue("MOC_APL_INP_FileTest");
	MOCMesh mocMesh("3_3cells.apl", "3_3cells_out.apl", MeshKernelType::MHT_KERNEL);
	//mocMesh.WriteSurfaceTecplotFile("filename.plt");
	//MOCMesh mocMesh("pin_c1.apl", "pin_c1.inp", MeshKernelType::MHT_KERNEL); 
	mocMesh.WriteTecplotFile("3_3cells_mMOD.plt","mMOD");
	mocMesh.WriteTecplotFile("3_3cells_mFUEL.plt","mFUEL");
	//mocMesh.WriteTecplotFile("3_3cells_mCLAD.plt", "mCLAD");
	//mocMesh.WriteTecplotFile("mUO2", "c5g721_cell12_mUO2.plt");
	//mocMesh.WriteTecplotFile("Zr4", "zr4.plt");
	//mocMesh.WriteTecplotFile("UO2", "u2o.plt");
	//mocMesh.InitMOCValue("pin_c1.inp","pin_c1.txt");
	//mocMesh.OutputStatus("pin_c1_out.inp");
	return;
}
void MapTest()
{
	//RegisterMapper("c5g72l.apl", "c5g72l_out.apl", "PIN9Coarse.msh");
	//MOCFieldsToCFD();

	std::string configFile = "MapFile_FileNames";
	std::vector<std::vector<std::string> > matches = GetMatchList(configFile);
	std::vector<std::string>& materialList = matches[0];
	std::vector<std::string>& regionList = matches[1];
	std::string mocMeshFile = GetFileName(configFile, "inputApl");
	std::string outMocMeshFile = GetFileName(configFile, "outputApl");

	MOCMesh mocMesh(mocMeshFile, outMocMeshFile);
	std::vector<std::string> inputVTKList;
	inputVTKList.push_back("fluid.vtk");
	inputVTKList.push_back("solid.vtk");
	Scalar ratio = GetValueInFile(configFile, "scaleRatio");
	MHTVTKReader reader(inputVTKList, ratio);
	for (size_t i = 0; i < inputVTKList.size(); i++)
	{
		Mesh* pmesh = reader.GetMeshListPtr()[i];
		//read cfd mesh and create solver
		CFDMesh cfdMesh(pmesh, MeshKernelType::MHT_KERNEL, i);
		cfdMesh.WriteTecplotFile("cfd_50_"+std::to_string(i)+".plt");
	}

	/*
	MHTVTKReader reader(cfdMeshFile, ratio);
	for (size_t i = 0; i < regionList.size(); i++)
	{
		//if (i > 0)break;
		int CFDMeshID = reader.GetIDOfRegion(regionList[i]);
		Mesh* pmesh = reader.GetMeshListPtr()[CFDMeshID];
		//read cfd mesh and create solver
		CFDMesh cfdMesh(pmesh, MeshKernelType::MHT_KERNEL, CFDMeshID);
		
		std::vector<SMocIndex> thh;
		for (int j = 0; j < 64; j++)
		{
			thh.push_back(SMocIndex(0, 51, j));
		}
		WriteTotecplot(mocMesh, cfdMesh, std::vector<int>{2559}, thh, "temp.plt");
		
		//std::string outfilename = "mesh_" + std::to_string(i) + ".plt";
		//cfdMesh.WriteTecplotFile(outfilename);
		Solver solverMapper(mocMesh, cfdMesh, materialList[i], true);
		solverMapper.CheckMappingWeights();
	}
	*/

	//InitCFDMeshValue(cfdMesh);
	//solver.CFDtoMOCinterception(ValueType::DENSITY);

	//mocMesh.OutputStatus("pin_c1.inp");

}
Scalar initialT(Scalar x, Scalar y, Scalar z)
{
	Vector PointCenter(x, y, z);
	Vector axisCenter(0.63, 0.63, 0.0);
	Vector axisNorm(0.0, 0.0, 1.0);
	Vector OP = PointCenter - axisCenter;
	Vector axicialLocation = (OP & axisNorm) * axisNorm;
	Scalar radius = (OP - axicialLocation).Mag();
	Scalar sigma = 1.0;
	Scalar temperature = (1000.0 / sqrt(2.0 * PI) / sigma) * exp(-radius * radius / 2 / pow(sigma, 2));
	return temperature;
}

Scalar initialRho(Scalar x, Scalar y, Scalar z)
{
	Scalar add = 100.0 * z * (5.0 - z);
	return 1000.0 + add;
}

Scalar initialHeatPower(Scalar x, Scalar y, Scalar z)
{
	Vector PointCenter(x, y, z);
	Vector axisCenter(0.63, 0.63, 0.0);
	Vector axisNorm(0.0, 0.0, 1.0);
	Vector OP = PointCenter - axisCenter;
	Vector axicialLocation = (OP & axisNorm) * axisNorm;
	Scalar radius = (OP - axicialLocation).Mag();
	Scalar sigma = 1.0;
	Scalar hp = (1000.0 / sqrt(2.0 * PI) / sigma) * exp(-radius * radius / 2 / pow(sigma, 2));
	hp = hp * z * (5.0 - z) / 6.25;
	return hp;
}

#include "./MHT_IO/VTKIO.h"
#include <fstream>

void VTKReaderTest()
{
	MHTVTKReader reader("pinWR.msh",1.0);
	std::vector<std::string> fileName;
	fileName.push_back("H2O.vtk");
	fileName.push_back("Zr4.vtk");
	std::vector<int> IDList;
	IDList.push_back(0);
	IDList.push_back(1);
	std::vector<std::string> fieldName;
	fieldName.push_back("heatpower");
	reader.ReadVTKFile(fileName, IDList, fieldName);
	reader.GetFieldIO(0).WriteTecplotField("heatpower_0.plt");
	reader.GetFieldIO(1).WriteTecplotField("heatpower_1.plt");
	return;
}

void EntranceOfRegister(std::vector<std::string>& fileNames)
{
	if (2 != fileNames.size())
	{
		std::cout << "Please give 2 file names if you are intending to register solver, like this:" << std::endl;
		std::cout << "NCCouple register (MOCMesh) (outMOCMesh)" << std::endl;
		Logger::LogError("inccorrect number of file names");
	}
	RegisterMapper(fileNames[0], fileNames[1]);
	return;
}

void EntranceOfCreateMapper(std::vector<std::string>& fileNames)
{
	if (0 != fileNames.size())
	{
		std::cout << "Please do not give any file name if you are intending to create mapper, like this:" << std::endl;
std::cout << "NCCouple createmapper" << std::endl;
Logger::LogError("inccorrect number of file names");
	}
	CreateMapper();
	return;
}

void EntranceOfCFDToMOC(std::vector<std::string>& fileNames)
{
	if (0 != fileNames.size())
	{
		std::cout << "Please do not give any file name if you are intending to map CFD to MOC solver, like this:" << std::endl;
		std::cout << "NCCouple cfdtomoc" << std::endl;
		Logger::LogError("inccorrect number of file names");
	}
	CFDFieldsToMOC();
	return;
}

void EntranceOfMOCToCFD(std::vector<std::string>& fileNames)
{
	if (0 != fileNames.size())
	{
		std::cout << "Please do not give any file name if you are intending to map MOC to CFD solver, like this:" << std::endl;
		std::cout << "NCCouple moctocfd" << std::endl;
		Logger::LogError("inccorrect number of file names");
	}
	MOCFieldsToCFD();
	return;
}

void RunWithParameters(std::vector<std::string>& parameters)
{
	if (0 == parameters.size())
	{
		Logger::LogError("no parameter given in RunWithParameters()");
	}
	std::vector<std::string> filenames;
	filenames.resize(parameters.size() - 1);
	for (int i = 1; i < parameters.size(); i++)
	{
		filenames[i - 1] = parameters[i];
	}
	if ("--help" == parameters[0])
	{
		if (0 != filenames.size())
		{
			Logger::LogError("no more parameter need to be given in --help");
		}
		DisplayHelpInfo();
	}
	else if ("clear" == parameters[0])
	{
		if (0 != filenames.size())
		{
			Logger::LogError("no more parameter need to be given in clear");
		}
		ClearMapFiles();
	}
	else if ("register" == parameters[0])
	{
		EntranceOfRegister(filenames);
	}
	else if ("createmapper" == parameters[0])
	{
		EntranceOfCreateMapper(filenames);
	}
	else if ("cfdtomoc" == parameters[0])
	{
		EntranceOfCFDToMOC(filenames);
	}
	else if ("moctocfd" == parameters[0])
	{
		EntranceOfMOCToCFD(filenames);
	}
	else
	{
		std::cout << "the following parameters are acceptabel: " << std::endl;
		std::cout << "(1) cfdtomoc" << std::endl;
		std::cout << "(2) clear" << std::endl;
		std::cout << "(3) createmapper" << std::endl;
		std::cout << "(4) moctocfd" << std::endl;
		std::cout << "(5) register" << std::endl;
		std::cout << "(6) --help" << std::endl;
		Logger::LogError("invalid parameter given: " + parameters[0]);
	}
	return;
}
void writeheatpower()
{
	std::ofstream ofs("heatpower_1cm.txt");
	int iNumber_Assembly = 9;
	for (int i = 1; i <= iNumber_Assembly; i++)
	{
		int iNum = 48*50;
		double iValue = 1.9e+10;
		ofs << i << "_" << iNum << std::endl;
		for (int j = 0; j < 240; j++)
		{
			for(int k = 0; k < 10; k++)
			{
				ofs << iValue << "  ";
			}
			ofs << std::endl;
		}
		//for (int k = 0; k < (iNum-40); k++)
		//{
			//ofs << iValue << "  ";
		//}
		ofs << std::endl;
	}
}

void VTKReadMeshTest()
{
	std::vector<std::string> vVTKname;
	vVTKname.push_back("1_PART-FLUID_couple.vtk");
	vVTKname.push_back("2_solid.vtk");

	std::vector<std::string> vFieldName;
	vFieldName.push_back("temperature");
	MHTVTKReader mhtvtkreader(vVTKname, vFieldName, 1.0);
	std::cout<<"mhtvtkreader.GetMeshListPtr()[0]->v_vertice.size() :" << mhtvtkreader.GetMeshListPtr()[0]->v_vertice.size()<<std::endl;
	mhtvtkreader.GetFieldIO(0).WriteTecplotField("temperature_0.plt");
	mhtvtkreader.GetFieldIO(1).WriteTecplotField("temperature_1.plt");
}

void FinalCFDFieldsToMOC()
{
	std::vector<std::vector<std::string> > matches = GetMatchList(configFile);
	std::string mocFieldFile = GetFileName(configFile, "inputInp");
	std::string outMocFieldFile = GetFileName(configFile, "outputInp");

	//checking file names
	std::vector<std::string> vStrMocFileName;
	for (int i = 1; i < g_iNumProcs; i++)
	{
		vStrMocFileName[i] = "Mocfield_" + std::to_string(i);
	}
		
	MOCMesh mocMesh(vStrMocFileName,true);
	mocMesh.InitMOCFromInputFile(mocFieldFile);
	RenameFile(outMocFieldFile, GetFileNameOfPrevious(outMocFieldFile, "inp"));
	mocMesh.OutputStatus(outMocFieldFile);
}

int main(int argc, char** argv)
{
	//SendFieldForTest(argc, argv);
	//return 0;

	int iNumberOfProcs, iMpiProcessID =-1, source;
	MPI_Status status;
	char message[100];
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &iMpiProcessID);
	MPI_Comm_size(MPI_COMM_WORLD, &iNumberOfProcs);
	//myid = 1;
	g_iMpiID = iMpiProcessID;
	g_iNumProcs = iNumberOfProcs;

	std::vector<std::string> parameterList;

	if (argc > 1)
	{
		parameterList.resize(argc - 1);
		for (int i = 1; i < argc; i++)
			parameterList[i - 1] = argv[i];
	}
	else
	{
		DisplayHelpInfo();
		return 0;
	}
	if (parameterList[0] != "createmapper" && parameterList[0] != "cfdtomoc" && parameterList[0] != "moctocfd")
	{
		RunWithParameters(parameterList);
		return 0;
	}
	if (iMpiProcessID > 0) {
		std::cout << "this is child,pid=" << getpid() << std::endl;
		//for debug
		if (iMpiProcessID == 100)
		{
			int num = 10;
			
			while (num == 10)
			{
			#ifdef _WIN32
				Sleep(10);
			#else
				sleep(10);
			#endif
			}
		}
		std::string strMessage;
		if (parameterList[0] == "createmapper"&&argc==2)
		{
			CreateMapper();
			strMessage = "Process " + std::to_string(iMpiProcessID) + " createmapper finished";
			strcpy(message, strMessage.data());// "caclulation finished!");
			MPI_Send(message, strlen(message) + 1, MPI_CHAR, 0, 99,MPI_COMM_WORLD);
		}
		else if (parameterList[0] == "moctocfd " && argc == 2)
		{
			MOCFieldsToCFD();
			strMessage = "Process " + std::to_string(iMpiProcessID) + " moctocfd finished";
			strcpy(message, strMessage.data());// "caclulation finished!");
			MPI_Send(message, strlen(message) + 1, MPI_CHAR, 0, 99,MPI_COMM_WORLD);
		}
		else if (parameterList[0] == "cfdtomoc" && argc == 2)
		{
			CFDFieldsToMOC();
			strMessage = "Process " + std::to_string(iMpiProcessID) + " cfdtomoc finished";
			strcpy(message, strMessage.data());// "caclulation finished!");
			MPI_Send(message, strlen(message) + 1, MPI_CHAR, 0, 99,MPI_COMM_WORLD);
		}
		else
		{
			Logger::LogError("Wrong input parameter");
			return -1;
		}
	}
	else if(iMpiProcessID ==0) 
	{
		std::cout << "this is 0 Process,pid=" << getpid() << std::endl;
		//if(parameterList[0] != "cfdtomoc")
		//{
			for (source = 1; source < iNumberOfProcs; source++) {
				MPI_Recv(message, 100, MPI_CHAR, source, 99,MPI_COMM_WORLD, &status);
				Logger::LogInfo(FormatStr("Main process received message from No.%d process: %s\n", source, message));
			}
		//}
		/*for debug 
		* int num = 10;
		while (num == 10)
		{
		#ifdef _WIN32
			Sleep(10);
		#else
			sleep(10);
		#endif
		}*/
		if(argc==2 && parameterList[0]=="createmapper")
		{
			ConvergeMocMapInfor();
			Logger::LogInfo("All createmapper finished.");
		}
		if (argc == 2 && parameterList[0] == "cfdtomoc")
		{
			FinalCFDFieldsToMOC();
			Logger::LogInfo("CFD to MOC finished.");
			/*
			std::vector<std::string> vStrMocFileName;
			MOCMesh mocMesh(vStrMocFileName,true);
			
			MPI_Datatype mpiMocField_type;
			InitMocFieldToMpiType(mpiMocField_type)

			std::vector<STRMocField> vReciveField;
			for (source = 1; source < iNumberOfProcs; source++) 
			{
				int iSize;
				MPI_Recv(&iSize, 1, MPI_INT, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				std::cout << "receive isize:" << iSize <<"from process:"<<source<< std::endl;
				vReciveField.resize(iSize);
				MPI_Recv(&vReciveField[0], 2, person_type, source, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				//printf("MPI process %d received person:\n\t- iAssembly = %d\n\t- icell = %d\n\t- name = %s\n", my_rank, vReciveField[1].iAssemblyIndex,
					//vReciveField[1].iCellIndex, vReciveField[1].cMaterialName);
				std::cout << "receive data from process: "<<source<< std::endl;
				SetFieldByMpiType(vReciveField);
			}
			MPI_Type_free(&mpiMocField_type);
			mocMesh.InitMOCFromInputFile(mocFieldFile);
			RenameFile(outMocFieldFile, GetFileNameOfPrevious(outMocFieldFile, "inp"));
			mocMesh.OutputStatus(outMocFieldFile);
			
			*/

		}
	}
	MPI_Finalize();
	return 0;
}