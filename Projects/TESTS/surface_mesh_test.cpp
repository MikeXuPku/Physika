/*
 * @file surface_mesh_test.cpp 
 * @brief Test the surface_mesh/vertex/edge/triangle class.
 * @author Sheng Yang
 * 
 * This file is part of Physika, a versatile physics simulation library.
 * Copyright (C) 2013- Physika Group.
 *
 * This Source Code Form is subject to the terms of the GNU General Public License v2.0. 
 * If a copy of the GPL was not distributed with this file, you can obtain one at:
 * http://www.gnu.org/licenses/gpl-2.0.html
 *
 */

#include <iostream>
#include "Physika_Core/Arrays/array.h"
#include "Physika_Geometry/Boundary_Meshes/vertex.h"
#include "Physika_Geometry/Boundary_Meshes/edge.h"
#include "Physika_Geometry/Boundary_Meshes/surface_mesh.h"
#include "Physika_IO/Surface_Mesh_IO/surface_mesh_io.h"
#include<fstream>

using namespace std;
using namespace Physika;

int main()
{/*
	SurfaceMesh<double> mesh;
	SurfaceMeshIO<double>::load(string(""), &mesh);
	fstream fileout("vertexData.txt");
	double ymin = 10000;
	for (int i = 0; i < mesh.numVertices(); ++i){
		Vector<double, 3> pos = mesh.vertexPosition(i);
		fileout << pos << endl;
		if (pos[1] < ymin)ymin = pos[1];
	}
	ymin = ymin - 10;
	for (int i = 0; i < mesh.numVertices(); ++i){
		Vector<double, 3> pos = mesh.vertexPosition(i);
		pos[1] = ymin;
		fileout << pos << endl;
	}
	fileout.close();
	return 0;
	*/

	SurfaceMesh<double> mesh;
	SurfaceMeshIO<double>::load(string(""), &mesh);
	double ymin = 10000;
	for (unsigned int i = 0; i < mesh.numVertices(); ++i){
		Vector<double, 3>pos = mesh.vertexPosition(i);
		if (pos[1] < ymin)ymin = pos[1];
	}
	ymin = ymin - 10;
	

	return 0;
}