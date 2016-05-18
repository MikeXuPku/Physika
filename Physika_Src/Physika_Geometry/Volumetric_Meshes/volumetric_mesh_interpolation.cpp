/*
* @file  volumetric_mesh_interpolation.cpp
* @brief implementation of methods defined in volumetric_mesh_interpolation.h
* @author Fei Zhu
*
* This file is part of Physika, a versatile physics simulation library.
* Copyright (C) 2013- Physika Group.
*
* This Source Code Form is subject to the terms of the GNU General Public License v2.0.
* If a copy of the GPL was not distributed with this file, you can obtain one at:
* http://www.gnu.org/licenses/gpl-2.0.html
*
*/

#include <algorithm>
#include "Physika_Core/Utilities/physika_assert.h"
#include "Physika_Core/Utilities/physika_exception.h"
#include "Physika_Geometry/Volumetric_Meshes/volumetric_mesh_interpolation.h"


namespace Physika{

template<typename Scalar, int Dim>
VolumetricMeshInterpolation<Scalar, Dim>::VolumetricMeshInterpolation(){
	volumetric_mesh_ = NULL;
	elements_.clear();
	weights_.clear();
}

template<typename Scalar, int Dim>
VolumetricMeshInterpolation<Scalar, Dim>::VolumetricMeshInterpolation(VolumetricMesh<Scalar, Dim>& vmesh){
	volumetric_mesh_ = &vmesh;  //here the equal operator will copy the content of vmesh for volumetric_mesh_
	elements_.clear();
	weights_.clear();
}

template<typename Scalar, int Dim>
VolumetricMeshInterpolation<Scalar, Dim>::~VolumetricMeshInterpolation(){

}

template<typename Scalar, int Dim>
void VolumetricMeshInterpolation<Scalar, Dim>::getPointsWeights(vector<Vector<Scalar, Dim>> &points, vector<unsigned int>&elements,vector<Scalar> &weights){
	weights.clear();
	vector<Scalar> temp;
	for (unsigned int i = 0; i < points.size(); ++i){
		volumetric_mesh_->interpolationWeights(elements[i], points[i], temp);
		for (unsigned int j = 0; j < volumetric_mesh_->eleVertNum(); ++j)weights.push_back(temp[j]);
	}
}

template<typename Scalar, int Dim>
void VolumetricMeshInterpolation<Scalar, Dim>::getSurfaceMeshWeights(SurfaceMesh<Scalar> &sMesh){
	//firstly acquire the containing or closest element information
	elements_.clear();
	for (unsigned int i = 0; i < sMesh.numVertices(); ++i){
		int element = getContainingElement(sMesh.vertexPosition(i));
		if (element == -1) element = getClosestElement(sMesh.vertexPosition(i));
		elements_.push_back(element);
	}

	//get interpolation weights
	weights_.clear();
	vector<Scalar> temp;
	for (unsigned int i = 0; i < sMesh.numVertices(); ++i){
		volumetric_mesh_->interpolationWeights(elements_[i], sMesh.vertexPosition(i), temp);
		for (unsigned int j = 0; j < volumetric_mesh_->eleVertNum(); ++j)weights_.push_back(temp[j]);
	}
}

template<typename Scalar, int Dim>
void VolumetricMeshInterpolation<Scalar, Dim>::getPointsWeightsElement(vector<Vector<Scalar, Dim>>&points, vector<unsigned int>&elements){
	elements.clear();
	for (unsigned int i = 0; i < points.size(); ++i){
		int element = getContainingElement(points[i]);
		if (element == -1) element = getClosestElement(points[i]);
		elements.push_back(element);
	}
}

template<typename Scalar, int Dim>
int VolumetricMeshInterpolation<Scalar, Dim>::getContainingElement(const Vector<Scalar, Dim> &pos){
	vector<Scalar> weights;
	for (unsigned int i = 0; i < volumetric_mesh_->eleNum(); ++i){
		if (volumetric_mesh_->containPoint(i, pos))return i;
	}
	return -1;
}

template<typename Scalar, int Dim>
int VolumetricMeshInterpolation<Scalar, Dim>::getClosestElement(const Vector<Scalar, Dim> &pos){
	double closest_distance = DBL_MAX;
	int closest_element = 0;
	unsigned int ele_vert_num = volumetric_mesh_->eleVertNum();
	for (unsigned int i = 0; i < volumetric_mesh_->eleNum(); ++i){
		Vector<Scalar, Dim> ele_center(0);
		for (unsigned int j = 0; j < ele_vert_num; ++j){
			ele_center += volumetric_mesh_->eleVertPos(i, j);
		}
		ele_center /= ele_vert_num;
		double distance = (ele_center - pos).normSquared();
		if (distance < closest_distance){
			closest_distance = distance;
			closest_element = i;
		}
	}
	return closest_element;
}

template<typename Scalar, int Dim>
void VolumetricMeshInterpolation<Scalar, Dim>::interpolate(vector<Vector<Scalar, Dim>> &u, vector<Vector<Scalar, Dim>> &targets, vector<unsigned int>&elements, vector<Scalar> &weights){
	PHYSIKA_ASSERT(weights.size() == volumetric_mesh_->eleVertNum()*elements.size());
	targets.clear();
	for (unsigned int i = 0; i < elements.size(); ++i){
		Vector<Scalar, Dim> target(0);
		for (unsigned int j = 0; j < volumetric_mesh_->eleVertNum(); ++j){
			unsigned int index = volumetric_mesh_->eleVertIndex(elements[i], j);
			target += u[index] * weights[i*volumetric_mesh_->eleVertNum() + j];
		}
		targets.push_back(target);
	}
}

template<typename Scalar, int Dim>
void VolumetricMeshInterpolation<Scalar, Dim>::interpolate(vector<Vector<Scalar, Dim>> &u, SurfaceMesh<Scalar> &sMesh){
	PHYSIKA_ASSERT(elements_.size() == sMesh.numVertices());
	for (unsigned int i = 0; i < sMesh.numVertices(); ++i){
		Vector<Scalar, Dim> target(0);
		for (unsigned int j = 0; j < volumetric_mesh_->eleVertNum(); ++j){
			unsigned int index = volumetric_mesh_->eleVertIndex(elements_[i], j);
			target += u[index] * weights_[i*volumetric_mesh_->eleVertNum() + j];
		}
		sMesh.setVertexPosition(i,target);
	}
}

template<typename Scalar, int Dim>
void VolumetricMeshInterpolation<Scalar, Dim>::save(const string &filename){
	std::fstream fileout(filename.c_str());
	if (!fileout) {
		std::cout << "open file error:" << filename << std::endl;
		return;
	}
	fileout << elements_.size() << ' ';
	for (unsigned int i = 0; i < elements_.size(); ++i)fileout << elements_[i] << ' ';
	fileout << weights_.size() << ' ';
	for (unsigned int i = 0; i < weights_.size(); ++i)fileout << weights_[i] << ' ';
	fileout.close();
}

template<typename Scalar, int Dim>
void VolumetricMeshInterpolation<Scalar, Dim>::load(const string &filename){
	std::fstream filein(filename.c_str());
	if (!filein){
		std::cout << "open file error:" << filename << std::endl;
		return;
	}
	elements_.clear();
	weights_.clear();
	unsigned int numElements, numWeights;
	filein >> numElements;
	for (unsigned int i = 0; i < numElements; ++i){
		unsigned int temp;
		filein >> temp;
		elements_.push_back(temp);
	}
	filein >> numWeights;
	for (unsigned int i = 0; i < numWeights; ++i){
		Scalar temp;
		filein >> temp;
		weights_.push_back(temp);
	}
	filein.close();
}

template class VolumetricMeshInterpolation<double, 3>;
template class VolumetricMeshInterpolation<float, 3>;

}  //end of namespace Physika
