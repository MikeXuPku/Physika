/*
* @file  volumetric_mesh_interpolation.h
* @brief get the position of a set of points by applying interpolation in a volumetric mesh.
* @author Leo Xu
*
* This file is part of Physika, a versatile physics simulation library.
* Copyright (C) 2013- Physika Group.
*
* This Source Code Form is subject to the terms of the GNU General Public License v2.0.
* If a copy of the GPL was not distributed with this file, you can obtain one at:
* http://www.gnu.org/licenses/gpl-2.0.html
*
*/

#ifndef PHYSIKA_GEOMETRY_VOLUMETRIC_MESHES_VOLUMETRIC_MESH_INTERPOLATION_H_
#define PHYSIKA_GEOMETRY_VOLUMETRIC_MESHES_VOLUMETRIC_MESH_INTERPOLATION_H_

#include <set>
#include <vector>
#include <string>
#include <fstream>
#include "Physika_Core/Vectors/vector_2d.h"
#include "Physika_Core/Vectors/vector_3d.h"
#include "Physika_Geometry/Volumetric_Meshes/volumetric_mesh.h"
#include "Physika_Geometry/Boundary_Meshes/surface_mesh.h"
using std::vector;
using std::string;

namespace Physika{

	/*
	* 
	*/

	template <typename Scalar, int Dim>
	class VolumetricMeshInterpolation
	{
	public:
		// constructor and destructor
		VolumetricMeshInterpolation();
		VolumetricMeshInterpolation(VolumetricMesh<Scalar, Dim> &vmesh);
		~VolumetricMeshInterpolation();

		//functions
		void getPointsWeights(vector<Vector<Scalar, Dim>> &points, vector<unsigned int>&elements, vector<Scalar> &weights);
		void getSurfaceMeshWeights(SurfaceMesh<Scalar> &sMesh);
		void getPointsWeightsElement(vector<Vector<Scalar, Dim>> &points, vector<unsigned int> &elements);
		int getContainingElement(const Vector<Scalar, Dim>& pos);
		int getClosestElement(const Vector<Scalar, Dim>& pos);
		void interpolate(vector<Vector<Scalar, Dim>> &u, vector<Vector<Scalar, Dim>> &targets, vector<unsigned int>&elements, vector<Scalar> &weights);
		void interpolate(vector<Vector<Scalar, Dim>> &u,SurfaceMesh<Scalar> &sMesh);

		void save(const string &filename);
		void load(const string &filename);


	private:
		VolumetricMesh<Scalar, Dim>* volumetric_mesh_;
		vector<unsigned int> elements_;
		vector<Scalar> weights_;
	};

}  //end of namespace Physika

#endif//PHYSIKA_GEOMETRY_VOLUMETRIC_MESHES_VOLUMETRIC_INTERPOLATION_MESH_H_
