#include<iostream>
#include <string>
#include <vector>
#include <cmath>
#include "Physika_Core/Vectors/vector_2d.h"
#include "Physika_Core/Vectors/vector_3d.h"
#include "Physika_Core/Vectors/vector_Nd.h"
#include "Physika_Core/Matrices/matrix_MxN.h"
#include "Physika_Core/Matrices/matrix_2x2.h"
#include "Physika_Core/Matrices/matrix_3x3.h"
#include "Physika_Core/Matrices/sparse_matrix.h"
#include "Physika_Core/Matrices/square_matrix.h"
#include "Physika_Core/Timer/timer.h"
#include "Physika_Geometry/Volumetric_Meshes/tri_mesh.h"
#include "Physika_Geometry/Volumetric_Meshes/tet_mesh.h"
#include "Physika_Geometry/Volumetric_Meshes/volumetric_mesh_interpolation.h"
#include "Physika_Geometry/Geometry_Intersections/tetrahedron_tetrahedron_intersection.h"
#include "Physika_IO/Volumetric_Mesh_IO/volumetric_mesh_io.h"
#include "Physika_Dynamics/Constitutive_Models/neo_hookean.h"

using Physika::VolumetricMeshInterpolation;
using Physika::VolumetricMeshIO;
using Physika::MatrixMxN;
using Physika::SquareMatrix;
using Physika::Vector;
using Physika::VectorND;
using Physika::Timer;
using std::fstream;
using std::cout;
using std::endl;
using std::vector;

class Homogenization
{
public:
	Homogenization();
	~Homogenization();
	bool loadFineMesh(const std::string &fine_mesh_file_name);
	bool loadCoarseMesh(const std::string &coarse_mesh_file_name);
	bool loadFineMeshDisplacement(const std::string &fine_displacement_file_name, float amplify = 1.0);
	bool loadFineMeshDisplacementWithWeight(const std::string &fine_displacement_file_name, const std::string &weight_file_name, float amplify = 1.0);
	bool loadFineMeshMaterialParameters(const std::string &youngs_file_name, const std::string &poisson_file_name);
	bool saveCoarseMeshMaterialParameters(const std::string &coarse_material_parameters_file_name, unsigned int start_ele, unsigned int end_ele);
	void solveCoarseMeshElasticTensor(unsigned int start_ele, unsigned int end_ele);
protected:
	void clearMeshData();

	void computeCoarseMeshDisplacementViaInterpolation(); //interpolate fine mesh displacement to coarse mesh
	
	std::vector<Physika::Vector<double, 3> > coarseMeshElementDisplacement(unsigned int ele_idx, unsigned int disp_idx) const;

	std::vector<Physika::Vector<double, 3> > fineMeshElementDisplacement(unsigned int ele_dx, unsigned int disp_idx) const;

	std::vector<unsigned int> elementsOnFineMesh(unsigned int coarse_ele_idx) const; //return indices of elements on fine mesh that correspond to a given element on coarse mesh
	//return indices of fine elements along with their intersected volume, the volume is approximated with monte carlo method

	void elementsOnFineMeshWithIntersectedVolume(unsigned int coarse_ele_idx, std::vector<unsigned int> &fine_elements, std::vector<double> &fine_ele_intersected_volume) const;
	
	//the following methods operate on elements in range [start_ele, end_ele)
	void findElementsOnFineMeshAndIntersectedVolumeForCoarseElements(unsigned int start_ele, unsigned int end_ele);

	void solveCoarseMeshMaterialParametersViaEnergyDeviationMinimization(unsigned int start_ele, unsigned int end_ele); //solve by minimizing potential energy deviation

	Physika::MatrixMxN<double> materialParametersFromYoungsAndPoisson(double youngs, double poisson);

	Physika::SquareMatrix<double, 3> tetDeformationGradient(vector<Vector<double,3>> reference, vector<Vector<double,3>>displacement);

	void prepareAnisotropicDirections(int subdivide_intensity);

protected:
	Physika::VolumetricMesh<double, 3> *fine_mesh_, *coarse_mesh_;
	std::vector<Physika::MatrixMxN<double> > fine_mesh_displacement_, coarse_mesh_displacement_; //disp[i] is the set of displacement on vert i, expressed as a Dimxn matrix
	std::vector<double> displacement_weight_; //i-th value is the weight correpond to i-th displacement for optimization
	std::vector<double> fine_mesh_youngs_modulus_, fine_mesh_poisson_ratio_; //material parameters for each fine mesh element
	std::vector<Physika::MatrixMxN<double> > coarse_mesh_material_parameters_; //material parameters for each coarse mesh element
	std::vector<std::vector<unsigned int> > fine_elements_for_coarse_elements_; //for each coarse element, store fine elements that intersect with it
	std::vector<std::vector<double> > fine_elements_intersected_volume_; //for each coarse element, store intersected fine element volumes
	unsigned int displacement_num_;//displacement number for both fine mesh and coarse mesh
	//for those coarse mesh elements that do not have any fine mesh elements for it, we compute its elastic tensor from default parameters
	double default_youngs_modulus_, default_poisson_ratio_;
	std::vector<Vector<double, 3>> anisotropic_directions_;
};


Homogenization::Homogenization()
:fine_mesh_(NULL), coarse_mesh_(NULL), displacement_num_(0),
default_youngs_modulus_(3.0e5), default_poisson_ratio_(0.3)
{
}

Homogenization::~Homogenization()
{
	clearMeshData();
}

bool Homogenization::loadFineMesh(const std::string &fine_mesh_file_name)
{
	if (fine_mesh_)
		delete fine_mesh_;
	fine_mesh_ = VolumetricMeshIO<double, 3>::load(fine_mesh_file_name);
	return fine_mesh_ == NULL ? false : true;
}

bool Homogenization::loadCoarseMesh(const std::string &coarse_mesh_file_name)
{
	if (coarse_mesh_)
		delete coarse_mesh_;
	coarse_mesh_ = VolumetricMeshIO<double, 3>::load(coarse_mesh_file_name);
	return coarse_mesh_ == NULL ? false : true;
}

bool Homogenization::loadFineMeshDisplacement(const std::string &fine_displacement_file_name, float amplify)
{
	std::fstream rfile;
	unsigned int str_num = 0;
	unsigned int line_num = 0;
	rfile.open(fine_displacement_file_name, std::ios::in);
	if (!rfile)
	{
		std::cerr << "Error: failed to open " << fine_displacement_file_name << endl;
		return false;
	}
	fine_mesh_displacement_.clear();

	unsigned int matrix_row_num;
	unsigned int matrix_col_num;
	rfile >> matrix_row_num >> matrix_col_num;
	std::cout << "row=" << matrix_row_num << "\n";
	std::cout << "col=" << matrix_col_num << "\n";
	this->displacement_num_ = matrix_col_num;
	if (matrix_row_num / 3 != this->fine_mesh_->vertNum()){
		cout << "the format of the fine mesh displacement file is wrong! vertex_num * 3 should equal matrix row num !!!!!!" << endl;
		return false;
	}
	//dis_matrix store each line of the displacement file, so there are three-lines for each vertex

	//each vertex contains a Dim x N matrix  
	unsigned int num = 0;
	for (unsigned int i = 0; i<this->fine_mesh_->vertNum(); ++i)
	{
		MatrixMxN<double> vetex_displacement_matrix(3, matrix_col_num);          //
		for (unsigned int k = 0; k<3; ++k)
		{
			for (unsigned int j = 0; j<matrix_col_num; ++j)
			{
				rfile >> vetex_displacement_matrix(k, j);
				vetex_displacement_matrix(k, j) *= amplify;
			}
		}
		this->fine_mesh_displacement_.push_back(vetex_displacement_matrix);
	}
	rfile.close();
	return true;
}

bool Homogenization::loadFineMeshDisplacementWithWeight(const std::string &fine_displacement_file_name, const std::string &weight_file_name, float amplify){
	bool loadDis = loadFineMeshDisplacement(fine_displacement_file_name, amplify);
	if (!loadDis)return false;
	
	std::fstream read_weight_file;
	read_weight_file.open(weight_file_name, std::ios::in);
	if (!read_weight_file)
	{
		std::cerr << "Error: failed to open " << weight_file_name << "\n";
		return false;
	}
	displacement_weight_.clear();
	for (unsigned int i = 0; i < this->displacement_num_; ++i){
		double weight_value;
		read_weight_file >> weight_value;
		displacement_weight_.push_back(weight_value);
	}
	read_weight_file.close();

	return true;
}

bool Homogenization::loadFineMeshMaterialParameters(const std::string &youngs_file_name, const std::string &poisson_file_name){
	std::fstream rfile;
	rfile.open(youngs_file_name, std::ios::in);
	if (!rfile)
	{
		std::cerr << "Error: failed to open " << youngs_file_name << endl;
		return false;
	}
	fine_mesh_youngs_modulus_.clear();
	while (!rfile.eof())
	{
		double youngs_modulus;
		rfile >> youngs_modulus;
		fine_mesh_youngs_modulus_.push_back(youngs_modulus);
	}
	rfile.close();

	std::fstream rfile1;
	rfile1.open(poisson_file_name, std::ios::in);
	if (!rfile1)
	{
		std::cerr << "Error: failed to open " << poisson_file_name << "\n";
		return false;
	}
	fine_mesh_poisson_ratio_.clear();
	while ((!rfile1.eof()) && (rfile1.peek() != std::ifstream::traits_type::eof()))
	{//file not empty and not reached end of the file
		double possion_rate;
		rfile1 >> possion_rate;
		fine_mesh_poisson_ratio_.push_back(possion_rate);
	}
	rfile1.close();
	return true;
}

bool Homogenization::saveCoarseMeshMaterialParameters(const std::string &coarse_material_parameters_file_name, unsigned int start_ele, unsigned int end_ele){
	PHYSIKA_ASSERT(coarse_mesh_);
	std::ofstream wfile(coarse_material_parameters_file_name);
	if (!wfile)
	{
		std::cerr << "Error: failed to create " << coarse_material_parameters_file_name << "\n";
		return false;
	}
	for (unsigned int coarse_ele_idx = start_ele; coarse_ele_idx<end_ele; ++coarse_ele_idx)
	{
		for (unsigned int i = 0; i<6; ++i)
		{
			wfile << coarse_mesh_material_parameters_[coarse_ele_idx](0, i) << " ";      //only write the lambda,upisilon,C,anisotropic direction v, 6 dimension parameter vector to file.
		}
		wfile << std::endl;
	}
	wfile.close();
	return true;
}

void Homogenization::solveCoarseMeshElasticTensor(unsigned int start_ele, unsigned int end_ele){
	PHYSIKA_ASSERT(fine_mesh_);
	PHYSIKA_ASSERT(coarse_mesh_);

	solveCoarseMeshMaterialParametersViaEnergyDeviationMinimization(start_ele, end_ele);
}

void Homogenization::clearMeshData()
{
	if (fine_mesh_)
		delete fine_mesh_;
	if (coarse_mesh_)
		delete coarse_mesh_;
	fine_mesh_ = NULL;
	coarse_mesh_ = NULL;
}

void Homogenization::computeCoarseMeshDisplacementViaInterpolation()
{
	PHYSIKA_ASSERT(coarse_mesh_);
	PHYSIKA_ASSERT(fine_mesh_);
	std::cout << "Interpolating displacements of fine mesh to coarse mesh..."<<endl;
	coarse_mesh_displacement_.clear();
	coarse_mesh_displacement_.resize(coarse_mesh_->vertNum(), MatrixMxN<double>(3, displacement_num_, 0.0));

	VolumetricMeshInterpolation<double, 3> interpoObj(*fine_mesh_);
	vector<Physika::Vector<double, 3>> coarse_mesh_points;
	vector<unsigned int> elements_index;
	vector<double> weights;
	for (unsigned int i = 0; i < coarse_mesh_->vertNum(); ++i)coarse_mesh_points.push_back(coarse_mesh_->vertPos(i));
	interpoObj.getPointsWeightsElement(coarse_mesh_points, elements_index);
	interpoObj.getPointsWeights(coarse_mesh_points, elements_index, weights);
	
	for (unsigned int i = 0; i < coarse_mesh_->vertNum(); ++i){
		unsigned int elementIndex = elements_index[i];
		for (unsigned int j = 0; j < 4; ++j){
			coarse_mesh_displacement_[i] += fine_mesh_displacement_[fine_mesh_->eleVertIndex(elementIndex, j)] * weights[i * 4 + j];
		}
	}
}

std::vector<Physika::Vector<double, 3> > Homogenization::coarseMeshElementDisplacement(unsigned int ele_idx, unsigned int disp_idx) const
{
	PHYSIKA_ASSERT(coarse_mesh_);
	PHYSIKA_ASSERT(ele_idx < coarse_mesh_->eleNum());
	PHYSIKA_ASSERT(disp_idx < displacement_num_);
	std::vector<Vector<double, 3> > result;
	for (unsigned int idx = 0; idx < coarse_mesh_->eleVertNum(ele_idx); ++idx)
	{
		unsigned int vert_idx = coarse_mesh_->eleVertIndex(ele_idx, idx);
		Vector<double, 3> disp;
		for (unsigned int dim = 0; dim < 3; ++dim)
			disp[dim] = coarse_mesh_displacement_[vert_idx](dim, disp_idx);
		result.push_back(disp);
	}
	return result;
}

std::vector<Vector<double, 3> > Homogenization::fineMeshElementDisplacement(unsigned int ele_idx, unsigned int disp_idx) const
{
	PHYSIKA_ASSERT(fine_mesh_);
	PHYSIKA_ASSERT(ele_idx < fine_mesh_->eleNum());
	PHYSIKA_ASSERT(disp_idx < displacement_num_);
	std::vector<Vector<double, 3> > result;
	for (unsigned int idx = 0; idx < fine_mesh_->eleVertNum(ele_idx); ++idx)
	{
		unsigned int vert_idx = fine_mesh_->eleVertIndex(ele_idx, idx);
		Vector<double, 3> disp;
		for (unsigned int dim = 0; dim < 3; ++dim)
			disp[dim] = fine_mesh_displacement_[vert_idx](dim, disp_idx);
		result.push_back(disp);
	}
	return result;
}

std::vector<unsigned int> Homogenization::elementsOnFineMesh(unsigned int coarse_ele_idx) const
{
	PHYSIKA_ASSERT(coarse_ele_idx < coarse_mesh_->eleNum());
	std::vector<unsigned int> fine_elements;
	//find those elements on fine mesh that intersect with this coarse element
	unsigned int ele_vert_num =  4;
	std::vector<Vector<double, 3>> ele_coarse_trait_3d(ele_vert_num);
	Vector<double, 3> vert_pos;
	for (unsigned int idx = 0; idx < coarse_mesh_->eleVertNum(coarse_ele_idx); ++idx)
	{
		unsigned int coarse_vert_idx = coarse_mesh_->eleVertIndex(coarse_ele_idx, idx);
		vert_pos = coarse_mesh_->vertPos(coarse_vert_idx);
		for (unsigned int dim = 0; dim < 3; ++dim)
			ele_coarse_trait_3d[idx][dim] = vert_pos[dim];
	}
	std::vector<Vector<double, 3> > ele_fine_trait_3d(ele_vert_num);
	for (unsigned int fine_ele_idx = 0; fine_ele_idx < fine_mesh_->eleNum(); ++fine_ele_idx)
	{
		for (unsigned int idx = 0; idx < fine_mesh_->eleVertNum(fine_ele_idx); ++idx)
		{
			unsigned int fine_vert_idx = fine_mesh_->eleVertIndex(fine_ele_idx, idx);
			vert_pos = fine_mesh_->vertPos(fine_vert_idx);
			for (unsigned int dim = 0; dim < 3; ++dim)
				ele_fine_trait_3d[idx][dim] = vert_pos[dim];

		}		
		if (Physika::GeometryIntersections::intersectTetrahedra(ele_coarse_trait_3d, ele_fine_trait_3d))
			fine_elements.push_back(fine_ele_idx);		
	}
	return fine_elements;
}

void Homogenization::elementsOnFineMeshWithIntersectedVolume(unsigned int coarse_ele_idx, std::vector<unsigned int> &fine_elements, std::vector<double> &fine_ele_intersected_volume) const
{
	fine_elements = elementsOnFineMesh(coarse_ele_idx);
	//approximate the intersected volume on each fine element with monte carlo method
	fine_ele_intersected_volume.resize(fine_elements.size());
	std::vector<Vector<double, 3> > sample_points;
	//first generate some sample points inside the coarse element, we use voxel-based uniform sampling

	Vector<double, 3> min_corner((std::numeric_limits<double>::max)()), max_corner((std::numeric_limits<double>::min)()); //the bounding box of the coarse element
	for (unsigned int idx = 0; idx < coarse_mesh_->eleVertNum(coarse_ele_idx); ++idx)
	{
		unsigned int coarse_vert_idx = coarse_mesh_->eleVertIndex(coarse_ele_idx, idx);
		Vector<double, 3> coarse_vert_pos = coarse_mesh_->vertPos(coarse_vert_idx);
		for (unsigned int dim = 0; dim < 3; ++dim)
		{
			min_corner[dim] = min_corner[dim] < coarse_vert_pos[dim] ? min_corner[dim] : coarse_vert_pos[dim];
			max_corner[dim] = max_corner[dim] > coarse_vert_pos[dim] ? max_corner[dim] : coarse_vert_pos[dim];
		}
	}

	unsigned int sample_res = 50; //50 sample points in each direction
	Vector<double, 3> sample_stride = (max_corner - min_corner) / (sample_res - 1.0);
	Vector<double, 3> sample;
	
	{
		for (unsigned int idx_x = 0; idx_x < sample_res; ++idx_x)
		for (unsigned int idx_y = 0; idx_y < sample_res; ++idx_y)
		for (unsigned int idx_z = 0; idx_z < sample_res; ++idx_z)
		{
			sample[0] = min_corner[0] + idx_x*sample_stride[0];
			sample[1] = min_corner[1] + idx_y*sample_stride[1];
			sample[2] = min_corner[2] + idx_z*sample_stride[2];
			if (coarse_mesh_->containPoint(coarse_ele_idx, sample))
				sample_points.push_back(sample);
		}
	}

	//M: total sample points; N: the sample points in fine element as well;
	//intersected volume: N/M*V_coarse
	double coarse_ele_vol = coarse_mesh_->eleVolume(coarse_ele_idx);
	for (unsigned int i = 0; i < fine_elements.size(); ++i)
	{
		unsigned int fine_ele_idx = fine_elements[i];
		unsigned int hit_num = 0;
		for (unsigned int j = 0; j < sample_points.size(); ++j)
		if (fine_mesh_->containPoint(fine_ele_idx, sample_points[j]))
			++hit_num;
		fine_ele_intersected_volume[i] = hit_num*1.0 / sample_points.size()*coarse_ele_vol;
	}
}

void Homogenization::findElementsOnFineMeshAndIntersectedVolumeForCoarseElements(unsigned int start_ele, unsigned int end_ele)
{
	PHYSIKA_ASSERT(coarse_mesh_);
	PHYSIKA_ASSERT(fine_mesh_);
	fine_elements_for_coarse_elements_.resize(coarse_mesh_->eleNum());
	fine_elements_intersected_volume_.resize(coarse_mesh_->eleNum());
	std::cout << "Finding fine elements for coarse elements from element " << start_ele << " to " << end_ele << "\n";
	for (unsigned int coarse_ele_idx = start_ele; coarse_ele_idx <end_ele ; ++coarse_ele_idx)
	{
		elementsOnFineMeshWithIntersectedVolume(coarse_ele_idx, fine_elements_for_coarse_elements_[coarse_ele_idx], fine_elements_intersected_volume_[coarse_ele_idx]);
	}
}

void Homogenization::solveCoarseMeshMaterialParametersViaEnergyDeviationMinimization(unsigned int start_ele, unsigned int end_ele)
{
	unsigned int count_ele_optimized_num = 0;
	double total_compute_time = 0.0;
	computeCoarseMeshDisplacementViaInterpolation();
	unsigned int coarse_mesh_ele_num = coarse_mesh_->eleNum();

	//solve elastic tensor of the coarse element via optimization, using the linear regression
	unsigned int dim_size = 3;
	coarse_mesh_material_parameters_.resize(coarse_mesh_ele_num, MatrixMxN<double>(1, 6));

	//precompute corresponding information for each coarse element
	findElementsOnFineMeshAndIntersectedVolumeForCoarseElements(start_ele, end_ele);
	std::cout << "Compute elastic tensor via energy minimization from element " << start_ele << " to " << end_ele << "\n";
	std::vector<unsigned int> isolate_coarse_elements; //coarse elements that do not have corresponding fine elements

	for (unsigned int coarse_ele_idx = start_ele; coarse_ele_idx <end_ele; ++coarse_ele_idx)
	{
		Timer timer;
		timer.startTimer();

		std::cout << "Computing elastic tensor for element " << coarse_ele_idx << " of coarse mesh...\n";
		std::vector<unsigned int> &ele_on_fine_mesh = fine_elements_for_coarse_elements_[coarse_ele_idx];
		if (ele_on_fine_mesh.empty()) //couldn't find any fine mesh elements for this coarse element
		{
			isolate_coarse_elements.push_back(coarse_ele_idx);
			continue;
		}
		count_ele_optimized_num++;

		//solve via energy deviation minimization
		vector<SquareMatrix<double, 3>> coarse_deformation_gradients;
		vector<Vector<double, 3>> reference, displacement;
		for (unsigned int i = 0; i < 4; ++i)reference.push_back(coarse_mesh_->eleVertPos(coarse_ele_idx, i));
		displacement.resize(4);
		for (unsigned int i = 0; i<displacement_num_; ++i){
			for (unsigned int j = 0; j<4; ++j){
				int vertIndex = coarse_mesh_->eleVertIndex(coarse_ele_idx, j);
				displacement.push_back(Vector<double, 3>(coarse_mesh_displacement_[vertIndex](0, i), coarse_mesh_displacement_[vertIndex](1, i), coarse_mesh_displacement_[vertIndex](2, i)));
			}
			coarse_deformation_gradients.push_back(tetDeformationGradient(reference, displacement));
		}//预先计算好粗糙网格的形变梯度

		//构造A矩阵
		MatrixMxN<double> A(displacement_num_, 3);
		//构造b向量
		VectorND<double> b(displacement_num_, 0);
		PHYSIKA_ASSERT(anisotropic_directions_.size()>0);
		for (unsigned int i = 0; i < displacement_num_; ++i){  //A矩阵前两列与anisotrpic direction无关 b向量也与其无关
			double I1 = (coarse_deformation_gradients[i] * coarse_deformation_gradients[i].transpose()).trace();
			double I3 = coarse_deformation_gradients[i].determinant();
			I3 = I3*I3;
			A(i, 0) = (I1 - log(I3) - 3) / 2;
			A(i, 1) = log(I3)*log(I3) / 8;
			A(i, 0) *= displacement_weight_[i];
			A(i, 1) *= displacement_weight_[i];
		}
		for (unsigned int fineEleIndex = 0; fineEleIndex < fine_elements_for_coarse_elements_[coarse_ele_idx].size(); ++fineEleIndex){
			Physika::NeoHookean<double, 3> fineMaterial;
			fineMaterial.setYoungsModulus(fine_mesh_youngs_modulus_[fineEleIndex]);
			fineMaterial.setPoissonRatio(fine_mesh_poisson_ratio_[fineEleIndex]);
			SquareMatrix<double, 3> fineDeformationGradient;
			vector<Vector<double, 3>> reference, displacement;
			for (unsigned int localIdx = 0; localIdx < 4; ++localIdx)reference.push_back(fine_mesh_->eleVertPos(fineEleIndex, localIdx));

			displacement.resize(4);
			for (unsigned int i = 0; i < displacement_num_; ++i){
				for (unsigned int localIdx = 0; localIdx < 4; ++localIdx){
					unsigned int fineVertIndex = fine_mesh_->eleVertIndex(fineEleIndex, localIdx);
					displacement.push_back(Vector<double, 3>(fine_mesh_displacement_[fineVertIndex](0, i), fine_mesh_displacement_[fineVertIndex](1, i), fine_mesh_displacement_[fineVertIndex](2, i)));
				}
				fineDeformationGradient = tetDeformationGradient(reference, displacement);
				b[i] += displacement_weight_[i]*fineMaterial.energyDensity(fineDeformationGradient)*fine_elements_intersected_volume_[coarse_ele_idx][fineEleIndex] / coarse_mesh_->eleVolume(coarse_ele_idx);
			}
		}


		double errorMax = (std::numeric_limits<double>::max)();
		unsigned int anisoDirectionMark = -1;
		for (unsigned int anisoIndex = 0; anisoIndex < anisotropic_directions_.size(); ++anisoIndex){
			for (unsigned int i = 0; i < displacement_num_; ++i){
				A(i, 2) = ((coarse_deformation_gradients[i] * anisotropic_directions_[anisoIndex]).norm() - 1);
				A(i, 2) = A(i, 2)*A(i, 2)*displacement_weight_[i];
			}
			VectorND<double> materialParameters = (A.transpose()*A).inverse()*(A.transpose()*b); //目前还没有防止过拟合的正则项
			double errorSquare = (A*materialParameters - b).normSquared();
			if (errorSquare < errorMax){
				coarse_mesh_material_parameters_[coarse_ele_idx](0, 0) = materialParameters[0];
				coarse_mesh_material_parameters_[coarse_ele_idx](0, 1) = materialParameters[1];
				coarse_mesh_material_parameters_[coarse_ele_idx](0, 2) = materialParameters[2];
				coarse_mesh_material_parameters_[coarse_ele_idx](0, 3) = anisotropic_directions_[anisoIndex][0];
				coarse_mesh_material_parameters_[coarse_ele_idx](0, 4) = anisotropic_directions_[anisoIndex][1];
				coarse_mesh_material_parameters_[coarse_ele_idx](0, 5) = anisotropic_directions_[anisoIndex][2];
			}
		}

		//energy deviation minimization ends
		timer.stopTimer();
		double time_elapse = timer.getElapsedTime();
		total_compute_time += time_elapse;
		std::cout << "Time for compute elastic tensor of element " << coarse_ele_idx << " is " << time_elapse << "s \n";
		//getchar();
	}

	//for those coarse elements that didn't find fine elements, set its ealsticity tensor to one of its neighbors
	unsigned int ele_num_with_default = 0;
	std::vector<unsigned int> isolate_coarse_elements_second_pass;
	for (unsigned int idx = 0; idx < isolate_coarse_elements.size(); ++idx)
	{
		unsigned int coarse_ele_idx = isolate_coarse_elements[idx];
		//find one element with elasticity tensor that shares the first node with this element
		unsigned int global_vert_idx = coarse_mesh_->eleVertIndex(coarse_ele_idx, 0);
		bool found = false;
		for (unsigned int ele_idx = 0; ele_idx < coarse_mesh_->eleNum(); ++ele_idx)
		{
			if (ele_idx == coarse_ele_idx)
				continue;
			for (unsigned int vert_idx = 0; vert_idx < coarse_mesh_->eleVertNum(ele_idx); ++vert_idx)
			if ((coarse_mesh_->eleVertIndex(ele_idx, vert_idx) == global_vert_idx)
				&& (std::find(isolate_coarse_elements.begin(), isolate_coarse_elements.end(), ele_idx) == isolate_coarse_elements.end()))
			{
				coarse_mesh_material_parameters_[coarse_ele_idx] = coarse_mesh_material_parameters_[ele_idx];
				found = true;
				break;
			}
			if (found)
				break;
		}
		if (!found)
			isolate_coarse_elements_second_pass.push_back(coarse_ele_idx);
	}
	//second pass
	for (unsigned int idx = 0; idx < isolate_coarse_elements_second_pass.size(); ++idx)
	{
		unsigned int coarse_ele_idx = isolate_coarse_elements_second_pass[idx];
		//find one element with elasticity tensor that shares the first node with this element
		unsigned int global_vert_idx = coarse_mesh_->eleVertIndex(coarse_ele_idx, 0);
		bool found = false;
		for (unsigned int ele_idx = 0; ele_idx < coarse_mesh_->eleNum(); ++ele_idx)
		{
			if (ele_idx == coarse_ele_idx)
				continue;
			for (unsigned int vert_idx = 0; vert_idx < coarse_mesh_->eleVertNum(ele_idx); ++vert_idx)
			if ((coarse_mesh_->eleVertIndex(ele_idx, vert_idx) == global_vert_idx)
				&& (std::find(isolate_coarse_elements_second_pass.begin(), isolate_coarse_elements_second_pass.end(), ele_idx) == isolate_coarse_elements_second_pass.end()))
			{
				coarse_mesh_material_parameters_[coarse_ele_idx] = coarse_mesh_material_parameters_[ele_idx];
				found = true;
				break;
			}
			if (found)
				break;
		}
		if (!found)
		{
			++ele_num_with_default;
			coarse_mesh_material_parameters_[coarse_ele_idx] = materialParametersFromYoungsAndPoisson(default_youngs_modulus_, default_poisson_ratio_);
		}
	}
	std::cout << isolate_coarse_elements.size() << " elements with special treatment, ";
	std::cout << isolate_coarse_elements_second_pass.size() << " elements with second pass special treatment, " << ele_num_with_default << " elements with default value.\n";
	std::cout << "Total optimized element No. is: " << count_ele_optimized_num << "\n";
	std::cout << "Total compute elastic tensor time: " << total_compute_time << " s; Average: " << total_compute_time / count_ele_optimized_num << " s/element.\n";
}

Physika::MatrixMxN<double> Homogenization::materialParametersFromYoungsAndPoisson(double youngs, double poisson){
	MatrixMxN<double> result(1, 6, 0.0);
	result(0, 0) = youngs / (2 * (1 + poisson));
	result(0, 1) = youngs*poisson / ((1 + poisson)*(1 - 2 * poisson));
	return result;
}

Physika::SquareMatrix<double, 3> Homogenization::tetDeformationGradient(vector<Vector<double, 3>> reference, vector<Vector<double, 3>>displacement){
	PHYSIKA_ASSERT(reference.size() == 4 && displacement.size() == 4);
	SquareMatrix<double, 3> material_matrix(reference[1] - reference[0], reference[2] - reference[0], reference[3] - reference[0]);
	SquareMatrix<double, 3> shape_matrix(displacement[1] - displacement[0], displacement[2] - displacement[0], displacement[3] - displacement[0]);
	shape_matrix += material_matrix;
	SquareMatrix<double, 3> result = shape_matrix.transpose() * (material_matrix.transpose()).inverse();
	return result;
}

void Homogenization::prepareAnisotropicDirections(int subdivide_intensity){
	anisotropic_directions_.clear();
	anisotropic_directions_.push_back(Vector<double, 3>(0, 0, 0));
	anisotropic_directions_.push_back(Vector<double, 3>(0, 0, 1));
	for (int i = 0; i < subdivide_intensity * 2; ++i)anisotropic_directions_.push_back(Vector<double, 3>(cos(M_PI/2*i / subdivide_intensity), sin(M_PI/2*i / subdivide_intensity), 0));
	for (int i = 1; i < subdivide_intensity; ++i){
		for (int j = 0; j < subdivide_intensity * 4; ++j){
			anisotropic_directions_.push_back(Vector<double,3>(cos(M_PI / 2 / subdivide_intensity*i)*cos(M_PI / 2 / subdivide_intensity*j), cos(M_PI / 2 / subdivide_intensity*i)*cos(M_PI / 2 / subdivide_intensity*j), sin(M_PI / 2 / subdivide_intensity*i)));
		}
	}
}