#include "../MHT_field/InteriorFaceZoneField.h"

#include "../MHT_common/SystemControl.h"
#include "../MHT_common/InterpolationTools.h"

#include "../MHT_mesh/Mesh.h"
#include "../MHT_mesh/MeshSupport.h"

#include "../MHT_field/ElementField.h"
#include "../MHT_field/FaceField.h"

// =======================================Scalar Type NodeField=============================================
//constructor
template<>
InteriorFaceZoneField<Scalar>::InteriorFaceZoneField(Mesh* UGB)
	:p_blockMesh(UGB)
{}

//interpolating face values from adjacent elements
template<>
void InteriorFaceZoneField<Scalar>::CalculateValue(ElementField<Scalar>& EF, FaceField<Scalar>& FF, const std::string& approach)
{
	if (false == EF.Assignment())
	{
		FatalError("Can not interpolate face values, element field did not pass Assignment check");
	}
	if (false == FF.Existence())
	{
		FF.Initialize();
	}
    Scalar (*interpolate)(const std::pair<Scalar, Scalar>&, const std::pair<Scalar, Scalar>&) = NULL;
	if ("linear" == approach)
	{
		interpolate = &LinearInterpolation;
	}
	else if ("hamonic"==approach)
	{
		interpolate = &HamonicInterpolation;
	}
	else
	{
		FatalError("unknown interpolation " + approach + " in InteriorFaceZoneField()");
	}
	FaceZone* p_faceZone = &(p_blockMesh->fz_interiorFaceZone);
	for (int i = 0; i < (int)p_faceZone->v_faceID.size(); i++)
	{
		//look up the block mesh for volumes and field values on owner&neighbor elements
		int faceID = p_faceZone->v_faceID[i];
		std::pair<int, int> elementIDs = this->p_blockMesh->SearchOwnerNeighbor(faceID);
		Scalar ownerValue = EF.GetValue(elementIDs.first);
		Scalar neighborValue = EF.GetValue(elementIDs.second);

		//calculate face values according to field values on neigboring elements
		int ownerID = elementIDs.first;
		int neighborID = elementIDs.second;
		Vector faceCenter = this->p_blockMesh->v_face[faceID].center;
		Vector d_Cf = faceCenter - this->p_blockMesh->v_elem[ownerID].center;
		Vector d_fF = this->p_blockMesh->v_elem[neighborID].center - faceCenter;
		//Parameters gc and gf are calculated by P160(6.30)
		Vector Sf = this->p_blockMesh->v_face[faceID].area;
		Vector e_f = Sf / Sf.Mag();
		std::pair<Scalar, Scalar> dis(d_Cf & e_f, d_fF & e_f);
		std::pair<Scalar, Scalar> values(ownerValue, neighborValue);
		Scalar obtainedValue = interpolate(dis, values);
		FF.SetValue(faceID, obtainedValue);
	}

	return;
}

// =======================================Vector Type NodeField=============================================
//constructor
template<>
InteriorFaceZoneField<Vector>::InteriorFaceZoneField(Mesh* UGB)
	:p_blockMesh(UGB)
{}

//interpolating face values from adjacent elements
template<>
void InteriorFaceZoneField<Vector>::CalculateValue(ElementField<Vector>& EF, FaceField<Vector>& FF, const std::string& approach)
{
	if (false == EF.Assignment())
	{
		FatalError("Can not interpolate face values, element field did not pass Assignment check");
	}
	if (false == FF.Existence())
	{
		FF.Initialize();
	}
    Vector(*interpolate)(const std::pair<Scalar, Scalar>&, const std::pair<Vector, Vector>&) = NULL;
	if ("linear" == approach)
	{
		interpolate = &LinearInterpolation;
	}
	else if ("hamonic" == approach)
	{
		FatalError("hamonic interpolation for vectors and tensor are in development");
		//interpolate = &HamonicInterpolation;
	}
	else
	{
		FatalError("unknown interpolation " + approach + " in InteriorFaceZoneField()");
	}
	FaceZone* p_faceZone = &(p_blockMesh->fz_interiorFaceZone);
	for (int i = 0; i < (int)p_faceZone->v_faceID.size(); i++)
	{
		//look up the block mesh for volumes and field values on owner&neighbor elements
		int faceID = p_faceZone->v_faceID[i];
		std::pair<int, int> elementIDs = this->p_blockMesh->SearchOwnerNeighbor(faceID);
		Vector ownerValue = EF.GetValue(elementIDs.first);
		Vector neighborValue = EF.GetValue(elementIDs.second);

		//calculate face values according to field values on neigboring elements
		int ownerID = elementIDs.first;
		int neighborID = elementIDs.second;
		Vector faceCenter = this->p_blockMesh->v_face[faceID].center;
		Vector d_Cf = faceCenter - this->p_blockMesh->v_elem[ownerID].center;
		Vector d_fF = this->p_blockMesh->v_elem[neighborID].center - faceCenter;
		//Parameters gc and gf are calculated by P160(6.30)
		Vector Sf = this->p_blockMesh->v_face[faceID].area;
		Vector e_f = Sf / Sf.Mag();
		std::pair<Scalar, Scalar> dis(d_Cf & e_f, d_fF & e_f);
		std::pair<Vector, Vector> values(ownerValue, neighborValue);
		Vector obtainedValue = interpolate(dis, values);
		FF.SetValue(faceID, obtainedValue);
	}
	return;
}
