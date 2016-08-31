#ifndef LPCVT_FUNCTIONS_H
#define LPCVT_FUNCTIONS_H

#include "combinatorics/mesh.h"
#include "combinatorics/delaunay.h"
#include "combinatorics/RVD.h"
#include "algebra/F_Lp.h"
#include "common/line_stream.h"
#include "combinatorics/exact/RVD_predicates.h"
#include <fstream>
#include "algebra/F_Lp.h"



namespace Geex {
	void create_pts(std::vector<std::vector<float> > initial_mesh_vertices, std::vector<vec3>& pts);

	int  test_combinatorics(std::vector<std::vector<float> > initial_mesh_vertices, std::vector<std::vector<float> > initial_mesh_triangles);
	
	void write_RDT(RestrictedVoronoiDiagram& RVD, std::vector<std::vector<float> > cvt_vertices, std::vector<std::vector<float> > cvt_triangles);
	
	int countRDTTriangles(RestrictedVoronoiDiagram& RVD);

	std::vector<std::vector<float> > initializeCubeVertices(float xmin=0, float ymin=0, float zmin=0);

     // Used by save_RDT().
    class WritePrimalTriangle{
    public:
        WritePrimalTriangle(std::vector<std::vector<float> >& triangles) : triangles_(&triangles){}
        void operator()(unsigned int i, unsigned int j, unsigned int k) const
        {
        
                std::cerr << "f " << i+1 << " " << j+1 << " " << k+1 << std::endl ;
		std::vector<float> p;
		const int array_size = 3;
		int pp[array_size] = {0, 0, 0};
		p.assign(pp, pp+array_size);
		p[0] = i+1;
		p[1] = j+1;
		p[2] = k+1;
		(*triangles_).push_back(p);
        }

    private:
        std::vector<std::vector<float> >* triangles_;
    } ;

    



}

#endif
