
#ifndef _Field_
#define _Field_

#include <string>
#include <vector>

#include "../MHT_common/Configuration.h"
#include "../MHT_mesh/MeshSupport.h"
#include "../MHT_field/FaceField.h"
#include "../MHT_field/ElementField.h"
#include "../MHT_field/NodeField.h"
#include "../MHT_field/InteriorFaceZoneField.h"
#include "../MHT_field/BoundaryFaceZoneField.h"

#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkFloatArray.h>
#include <vtkDataArray.h>

//Field on unstructured grid block
template<class Type>
class Field
{
public:

	enum outType
	{
		otCellNode,
		otCellCenter
	};

	//pointer to dependent mesh 
	Mesh* p_blockMesh;

	//name of this field
    std::string st_name;

	//element field
	ElementField<Type> elementField;

	//element field on previous time step or interation step
	ElementField<Type> elementField0;

	//face field
	FaceField<Type> faceField;

	//node field
	NodeField<Type> nodeField;

	//corresponding to interior face zone list
	InteriorFaceZoneField<Type> interiorFaceField;

	//corresponding to external face zone list
    std::vector<BoundaryFaceZoneField<Type> > v_BoundaryFaceField;

public:
	//constructor
    Field(Mesh*, const std::string &fileName="unknown Field");

	//constructor, initialization with a uniform field;
    Field(Mesh*, Type, const std::string& fileName = "unknown Field");

	//constructor, initialization with a non-uniform field;
    Field(Mesh*, const std::vector<Type>&, const std::string& fileName = "unknown Field");

	//constructor, initialization with a user defined function;
    Field(Mesh*, Type(Scalar, Scalar, Scalar), const std::string& fileName = "unknown Field");

	//creating boundary fields according to mesh information
    void CreateBoundaryField();

	//specify a special boundary condition on a boundary face zone
    void SetBoundaryCondition(const std::string &, const std::string &);

	//specify a uniform general-type (ABC-form) boundary condition on a boundary face zone
    void SetBoundaryCondition(const std::string&,Scalar,Scalar,Type);

	//specify a non-uniform general-type (ABC-form) boundary condition on a boundary face zone
    void SetBoundaryCondition(const std::string&, Scalar, Scalar, Type(*udf)(Scalar,Scalar,Scalar));

	//specify a non-uniform general-type (ABC-form) boundary condition on a boundary face zone
    void SetBoundaryCondition(const std::string&, Scalar, Scalar, std::vector<Type>);

	//check wether all boundaries have been correctly specified
    void CheckBoundaryCondition() const;

	//copy all boundary conditions to another field
	void CopyBoundaryConditionTo(Field<Type>&);

	//copying elementField to elementField0
	void SaveOld();

	//calculate face values on boundary
    void UpdateBoundary(const std::string& option = "default");

	//calculate face values: interior and boundary.
    void ElementToFace(const std::string& method="linear");

	//calculate face values: interior and boundary
    void ElementToFace(const std::string&, const std::string&);
	
	//interpolate from element to node
	void ElementToNode();

	//interpolation from node to element
	void NodeToElement();

	//Kong Ling supplemented on 2019/12/23
	//extract values on a specfic boundary
	std::vector<Type> ExtractBoundary(const std::string&);

	//Kong Ling modified on 2019/12/23
	//load values on a specific boundary with a ordered value list
    void LoadBoundary(const std::string&, const std::vector<Type>&);

	//Defining field operators
    Field<Type>& operator = (const Field<Type>&);

	//Writing field Tecplot file for visualization
    void WriteTecplotField(const std::string& outMshFileName);

    void WriteTecplotField(const std::string& outMshFileName, outType outFileType);

	//Writing BC Tecplot file for visualization
    void WriteTecplotField(const std::string& patchName, const std::string& outMshFileName);

    void WriteVTK_Field(const std::string& outMshFileName);

	void ReadVTK_Field(const std::string& inVTKFileName);

	void ReadVTKGridField(vtkSmartPointer<vtkUnstructuredGrid> uGrid, const std::string ArryName);
};
template<>
void Field<Scalar>::ReadVTKGridField(vtkSmartPointer<vtkUnstructuredGrid> uGrid, const std::string ArryName);

#endif
