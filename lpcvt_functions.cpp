        
#include "combinatorics/mesh.h"
#include "combinatorics/delaunay.h"
#include "combinatorics/RVD.h"
#include "algebra/F_Lp.h"
#include "common/line_stream.h"
#include <fstream>
#include "combinatorics/exact/RVD_predicates.h"
#include "algebra/F_Lp.h"
#include "lpcvt_functions.h"
        
// // create shared library of this with sudo g++ -std=c++11 -fPIC -shared lpcvt_functions.cpp -o /usr/lib/libLpCVT_functions.so

// in 3Depict all the references are required like this 
//sudo g++ -fPIC -shared lpcvt_functions.cpp combinatorics/exact/RVD_predicates.cpp combinatorics/mesh.cpp combinatorics/delaunay.cpp combinatorics/delaunay_CGAL.cpp algebra/F_Lp.cpp  -lCGAL -lgmp -o /usr/lib/libLpCVT_functions.so

namespace Geex {

	void create_pts(std::vector<std::vector<float> > initial_mesh_vertices, std::vector<vec3>& pts)
	{
		pts.clear() ;
		for (int i=0; i<initial_mesh_vertices.size(); i++)
		{
			double x = (double) initial_mesh_vertices[i][0];
			double y = (double) initial_mesh_vertices[i][1];
			double z = (double) initial_mesh_vertices[i][2];
			vec3 p;
			p[0] = x;
			p[1] = y;
			p[2] = z;
			pts.push_back(p) ;
		}
	}
		
	void write_RDT(RestrictedVoronoiDiagram& RVD, std::vector<std::vector<float> > rdt_vertices, 
	std::vector<std::vector<float> > rdt_triangles) 
	{	
	
		std::cerr << "Computing and writing RDT" << std::endl ;
		for(int i=0; i<RVD.delaunay()->nb_vertices(); i++) 
		{
		    rdt_vertices[i][0] = RVD.delaunay()->vertex(i)[0];
		    rdt_vertices[i][1] = RVD.delaunay()->vertex(i)[1];
		    rdt_vertices[i][2] = RVD.delaunay()->vertex(i)[2];
		}

		RVD.for_each_primal_triangle(WritePrimalTriangle(rdt_triangles));
	}
	
	void write_RDTByReference(RestrictedVoronoiDiagram& RVD, std::vector<std::vector<float> > &rdt_vertices, 
	std::vector<std::vector<float> > &rdt_triangles) 
	{	
	
		std::cerr << "Computing and writing RDT by reference" << std::endl ;
		for(int i=0; i<RVD.delaunay()->nb_vertices(); i++) 
		{
		    rdt_vertices[i][0] = RVD.delaunay()->vertex(i)[0];
		    rdt_vertices[i][1] = RVD.delaunay()->vertex(i)[1];
		    rdt_vertices[i][2] = RVD.delaunay()->vertex(i)[2];
		}
		
		RVD.for_each_primal_triangle(WritePrimalTriangle(rdt_triangles));

	}
	
	int countRDTTriangles(RestrictedVoronoiDiagram& RVD)
	{
		
		int tri_counter = 0;
		RVD.for_each_primal_triangle(CountPrimalTriangles(tri_counter));
		return tri_counter;
		
	}
	    
	int getCombinatorialStructureOfFLp(std::vector<std::vector<float> > initial_mesh_vertices, std::vector<std::vector<float> > initial_mesh_triangles) 

	{
		Mesh M ;
		unsigned int nb_borders = M.receiveVerticesAndTriangles(initial_mesh_vertices, initial_mesh_triangles);
		std::vector<vec3> pts ;
		create_pts(initial_mesh_vertices, pts);
		Delaunay* delaunay = Delaunay::create("CGAL") ;
		RestrictedVoronoiDiagram RVD(delaunay, &M) ;

		delaunay->set_vertices(pts) ;
	
		// initialize a new vector which holds the new vertices of the delaunay triangulation
		int number_of_rdt_vertices = RVD.delaunay()->nb_vertices();
	
		int xyzs = 3;
		std::vector<std::vector<float> > rdt_vertices(number_of_rdt_vertices, std::vector<float>(xyzs));

		int number_of_triangles = countRDTTriangles(RVD);
		int number_of_vertex_indices_per_triangle = 3;
		std::vector<std::vector<float> > rdt_triangles(number_of_triangles, std::vector<float>(number_of_vertex_indices_per_triangle));

		write_RDT(RVD, rdt_vertices, rdt_triangles) ;
	
		delete delaunay ;
		
		return rdt_triangles.size();
	}
	
	
	void getCombinatorialStructureOfFLpByReference(std::vector<std::vector<float> > initial_mesh_vertices, std::vector<std::vector<float> > initial_mesh_triangles, 
	std::vector<std::vector<float> > &rdt_vertices, std::vector<std::vector<float> > &rdt_triangles)
	{
		Mesh M ;
		unsigned int nb_borders = M.receiveVerticesAndTriangles(initial_mesh_vertices, initial_mesh_triangles);
		std::vector<vec3> pts ;
		create_pts(initial_mesh_vertices, pts);
		Delaunay* delaunay = Delaunay::create("CGAL") ;
		RestrictedVoronoiDiagram RVD(delaunay, &M) ;

		delaunay->set_vertices(pts) ;
		
		// this is only valid for vertices as they are referenced in writeRDT by index
		// the triangles are passed to RCD.each_primal_triangle which only takes an
		// initialized but not presized vector
		//resize the referenced vectors as the actual sizes are not know until here
		int number_of_rdt_vertices = RVD.delaunay()->nb_vertices();
		int xyzs = 3;
		rdt_vertices.resize(number_of_rdt_vertices);
		for (int i = 0; i < number_of_rdt_vertices; ++i)
		{
   			rdt_vertices[i].resize(xyzs);
		}

		write_RDTByReference(RVD, rdt_vertices, rdt_triangles) ;
	
		delete delaunay ;
	}
	
    /**
     * Gets the combinatorics of the integration simplices,
     * i.e. 10 integers per integration simplex.
     * (see Section 3.1 in the paper)
     * Returns also the array of C vertices (three per integration simplex).
     * Since they are easy to get during the combinatorial phase, they are
     * computed here and kept for the algebraic phase.
     *
     * In 2D mode (volume = false), returns also the array F.
     *   F[i] indicates the facet that contains the i-th integration simplex.
     *
     */
    void get_combinatorics(
        Mesh* M, const std::vector<vec3>& pts, 
        std::vector<int>& I, std::vector<vec3>& C, std::vector<int>& F, bool volume
    ) {
        Delaunay* delaunay = Delaunay::create("CGAL") ;
        delaunay->set_vertices(pts) ;
	RestrictedVoronoiDiagram RVD(delaunay,M) ;
	RVD.set_symbolic(true) ;
	RVD.for_each_triangle(MemorizeIndicesAndFacets(RVD,I,C,F)) ;
           
        delete delaunay ;
    }
    
     //Computes F_{L_p} and its gradient.
    float compute_F_g(Mesh* m, const std::vector<vec3>& pts, unsigned int p, bool volume) {
        std::cerr << "nb pts = " << pts.size() << "   nb facets = " << m->nb_facets() << std::endl ;
        std::vector<int> I ;
        std::vector<vec3> C ;
        std::vector<int> F ;
        get_combinatorics(m, pts, I, C, F, volume) ;
        unsigned int nb_integration_simplices = (unsigned int)I.size() / 10 ;
        std::vector<mat3> M(nb_integration_simplices) ;
        for(unsigned int i=0; i<M.size(); i++) {
            M[i].load_identity() ; 
                // or replace with anisotropy field
                //   In 2D: use F[i] to retreive the index of the facet that contains
                //      the current integration simplex (and access an array of per-facet anisotropy).
                //   In 3D: use geometric search from the centroid of the current
                //      integration simplex.
        }
        std::vector<plane3> Q(m->nb_facets()) ;
        for(unsigned int i=0; i<m->nb_facets(); i++) {
            Q[i] = m->facet_plane(i) ;
        }
        std::vector<double> g(pts.size() * 3) ;
        double f = compute_F_Lp(volume, p, m, I, C, pts, Q, M, g) ;
        double gnorm = 0.0 ;
        for(unsigned int i=0; i<g.size(); i++) {
            gnorm += g[i]*g[i] ;
        }
        gnorm = ::sqrt(gnorm) ;
        std::cerr.precision(16);
        std::cerr << (volume ? "volume " : "surface ") 
                  << "F_L" << p << ":" 
                  << "f=" << std::scientific << f << "  g=" << gnorm << std::endl ;
        return f;
    }

    float test_algebra(std::vector<std::vector<float> > initial_mesh_vertices, std::vector<std::vector<float> > initial_mesh_triangles) {
        Mesh M ;
        unsigned int nb_borders = M.receiveVerticesAndTriangles(initial_mesh_vertices, initial_mesh_triangles);
        std::vector<vec3> pts ;
        create_pts(initial_mesh_vertices, pts);
        std::cerr << "          ========== unit test algebraic surface LpCVT test ======" << std::endl ;
        int p_norm = 2;
        float FL_p = compute_F_g(&M, pts, p_norm, false) ;
        return FL_p;
    }
    
std::vector<std::vector<float> > initializeCubeVertices(float xmin, float ymin, float zmin)
{
	int number_of_vertices = 8;
	int xyz = 3;
	std::vector<std::vector<float> > vertices(number_of_vertices, std::vector<float>(xyz));

	vertices[0][0] = xmin+1;
	vertices[0][1] = ymin;
	vertices[0][2] = zmin;

	vertices[1][0] = xmin+1;
	vertices[1][1] = ymin;
	vertices[1][2] = zmin+1;

	vertices[2][0] = xmin;
	vertices[2][1] = ymin;
	vertices[2][2] = zmin+1;

	vertices[3][0] = xmin;
	vertices[3][1] = ymin;
	vertices[3][2] = zmin;

	vertices[4][0] = xmin+1;
	vertices[4][1] = ymin+1;
	vertices[4][2] = zmin;

	vertices[5][0] = xmin+1;
	vertices[5][1] = ymin+1;
	vertices[5][2] = zmin+1;

	vertices[6][0] = xmin;
	vertices[6][1] = ymin+1;
	vertices[6][2] = zmin+1;

	vertices[7][0] = xmin;
	vertices[7][1] = ymin+1;
	vertices[7][2] = zmin;
	
	return vertices;
}

}
