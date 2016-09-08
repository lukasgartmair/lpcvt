#ifndef LPCVT_FUNCTIONS_H
#define LPCVT_FUNCTIONS_H

#include "combinatorics/mesh.h"
#include "combinatorics/delaunay.h"
#include "combinatorics/RVD.h"
#include "combinatorics/exact/RVD_predicates.h"
#include "algebra/F_Lp.h"
#include "common/line_stream.h"
#include <fstream>
#include "algebra/F_Lp.h"



namespace Geex {
	void create_pts(std::vector<std::vector<float> > seeds, std::vector<vec3>& pts);

	int getCombinatorialStructureOfFLp(std::vector<std::vector<float> > seeds, std::vector<std::vector<float> > initial_mesh_vertices, 
		std::vector<std::vector<float> > initial_mesh_triangles); 
	
	void getCombinatorialStructureOfFLpByReference(std::vector<std::vector<float> > seeds, std::vector<std::vector<float> > initial_mesh_vertices, 
		std::vector<std::vector<float> > initial_mesh_triangles, std::vector<std::vector<float> > &rdt_vertices, std::vector<std::vector<float> > &rdt_triangles);
	
	void write_RDT(RestrictedVoronoiDiagram& RVD, std::vector<std::vector<float> > cvt_vertices, std::vector<std::vector<float> > cvt_triangles);
	
	void write_RDTByReference(RestrictedVoronoiDiagram& RVD, std::vector<std::vector<float> > &rdt_vertices, 
	std::vector<std::vector<float> > &rdt_triangles);
	
	int countRDTTriangles(RestrictedVoronoiDiagram& RVD);

	std::vector<std::vector<float> > initializeCubeVertices(float xmin=0, float ymin=0, float zmin=0);

   	float getAlgebraicStructureOfFLpByReference(std::vector<std::vector<float> > seeds, std::vector<std::vector<float> > initial_mesh_vertices, 
   		std::vector<std::vector<float> > initial_mesh_triangles, std::vector<std::vector<float> > &rdt_vertices, std::vector<std::vector<float> > &rdt_triangles);
	
	float compute_F_g(Mesh* m, const std::vector<vec3>& pts, unsigned int p, bool volume);
	
	void get_combinatorics(Mesh* M, const std::vector<vec3>& pts, std::vector<int>& I, std::vector<vec3>& C, std::vector<int>& F, bool volume);
	
	std::vector<std::vector<float> > generateSeedsLyingOnTriangleSurfaces(std::vector<std::vector<float> > initial_mesh_vertices, std::vector<std::vector<float> > 			initial_mesh_triangles);
	

    class WritePrimalTriangle{
    public:
        WritePrimalTriangle(std::vector<std::vector<float> >& triangles) : triangles_(&triangles){}
        void operator()(unsigned int i, unsigned int j, unsigned int k) const
        {
        
                //std::cerr << "f " << i+1 << " " << j+1 << " " << k+1 << std::endl ;
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
    
    
    class CountPrimalTriangles{
    public:
        CountPrimalTriangles(int& tri_counter) : tri_counter_(&tri_counter){}
        void operator()(unsigned int i, unsigned int j, unsigned int k) const
        {
		(*tri_counter_)+=1;

        }

    private:
        int* tri_counter_;
    } ;
    
    /**
     * Used by get_combinatorics() in surface mode
     */
    class MemorizeIndicesAndFacets{
    public:
        MemorizeIndicesAndFacets(
            const RestrictedVoronoiDiagram& RVD_in,
            std::vector<int>& I_in,
            std::vector<vec3>& C_in,
            std::vector<int>& F_in
        ) : RVD(RVD_in), I(I_in), C(C_in), F(F_in) {
            I.resize(0) ;
            C.resize(0) ;
            F.resize(0) ;
        }

        void operator() (
            unsigned int i, 
            const VertexEdge& v1, 
            const VertexEdge& v2, 
            const VertexEdge& v3
        ) const {
            I.push_back(i) ;
            I.push_back(v1.sym[2]) ;
            I.push_back(v1.sym[1]) ;
            I.push_back(v1.sym[0]) ;
            I.push_back(v2.sym[2]) ;
            I.push_back(v2.sym[1]) ;
            I.push_back(v2.sym[0]) ;
            I.push_back(v3.sym[2]) ;
            I.push_back(v3.sym[1]) ;
            I.push_back(v3.sym[0]) ;
            F.push_back(RVD.current_facet()) ;
            C.push_back(v1) ;
            C.push_back(v2) ;
            C.push_back(v3) ;
        }
    private:
        const RestrictedVoronoiDiagram& RVD ;
        std::vector<int>& I ;
        std::vector<vec3>& C ;
        std::vector<int>& F ;
    } ;


}

#endif 
