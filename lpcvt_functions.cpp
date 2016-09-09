        
#include "combinatorics/mesh.h"
#include "combinatorics/delaunay.h"
#include "combinatorics/RVD.h"
#include "algebra/F_Lp.h"
#include "common/line_stream.h"
#include <fstream>
#include "combinatorics/exact/RVD_predicates.h"
#include "algebra/F_Lp.h"
#include "lpcvt_functions.h"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random.hpp>
#include <time.h>
#include <numeric>
#include <math.h>       /* floor */
        
        
// // create shared library of this with sudo g++ -std=c++11 -fPIC -shared lpcvt_functions.cpp -o /usr/lib/libLpCVT_functions.so

// in 3Depict all the references are required like this 
//sudo g++ -fPIC -shared lpcvt_functions.cpp combinatorics/exact/RVD_predicates.cpp combinatorics/mesh.cpp combinatorics/delaunay.cpp combinatorics/delaunay_CGAL.cpp algebra/F_Lp.cpp  -lCGAL -lgmp -o /usr/lib/libLpCVT_functions.so

namespace Geex {

	void create_pts(std::vector<std::vector<float> > seeds, std::vector<vec3>& pts)
	{
		pts.clear() ;
		for (int i=0; i<seeds.size(); i++)
		{
			double x = (double) seeds[i][0];
			double y = (double) seeds[i][1];
			double z = (double) seeds[i][2];
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
	    
	int getCombinatorialStructureOfFLp(std::vector<std::vector<float> > seeds, std::vector<std::vector<float> > initial_mesh_vertices, std::vector<std::vector<float> > 			initial_mesh_triangles) 

	{
		Mesh M ;
		unsigned int nb_borders = M.receiveVerticesAndTriangles(initial_mesh_vertices, initial_mesh_triangles);
		std::vector<vec3> pts ;
		create_pts(seeds, pts);
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
	
	
	void getCombinatorialStructureOfFLpByReference(std::vector<std::vector<float> > seeds, std::vector<std::vector<float> > initial_mesh_vertices, 
		std::vector<std::vector<float> > initial_mesh_triangles, std::vector<std::vector<float> > &rdt_vertices, std::vector<std::vector<float> > &rdt_triangles)
	{
		Mesh M ;
		unsigned int nb_borders = M.receiveVerticesAndTriangles(initial_mesh_vertices, initial_mesh_triangles);
		std::vector<vec3> pts ;
		create_pts(seeds, pts);
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

	float getAlgebraicStructureOfFLpByReference(std::vector<std::vector<float> > seeds, std::vector<std::vector<float> > initial_mesh_vertices, 
		std::vector<std::vector<float> > initial_mesh_triangles, std::vector<std::vector<float> > &rdt_vertices, std::vector<std::vector<float> > &rdt_triangles) {

		Mesh M ;
		unsigned int nb_borders = M.receiveVerticesAndTriangles(initial_mesh_vertices, initial_mesh_triangles);
		std::vector<vec3> pts ;
		create_pts(seeds, pts);
		std::cerr << "          ========== unit test algebraic surface LpCVT test ======" << std::endl ;
		// p has to be even!
		//Assertion `(p/2)*2 == p' failed
		int p_norm = 2;
		float FL_p = compute_F_g(&M, pts, p_norm, false) ;

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

		return FL_p;
	}
	
	//http://math.stackexchange.com/questions/128991/how-to-calculate-area-of-3d-triangle
	
	double calculateHeron(double a, double b, double c)
	{
	    double s = (a + b + c) / 2; 
	    double area = pow((s*(s-a) * (s-b)*(s-c)), 0.5);       
	    return area;
	}

	double calculateDistance3d(double x1, double y1, double z1, double x2, double y2, double z2)
	{
	    double a = pow((x1-x2),2.0)+pow((y1-y2),2.0) + pow((z1-z2),2.0);
	    double d = pow(a,0.5);  
	    return d;
	}
	
	double calculateTriangleArea(std::vector<std::vector<float> > initial_mesh_vertices, std::vector<float> triangle)
	{				

		double x1 = initial_mesh_vertices[triangle[0]][0];
		double x2 = initial_mesh_vertices[triangle[0]][1];
		double x3 = initial_mesh_vertices[triangle[0]][2];
		double y1 = initial_mesh_vertices[triangle[1]][0];
		double y2 = initial_mesh_vertices[triangle[1]][1];
		double y3 = initial_mesh_vertices[triangle[1]][2];
		double z1 = initial_mesh_vertices[triangle[2]][0];
		double z2 = initial_mesh_vertices[triangle[2]][1];
		double z3 = initial_mesh_vertices[triangle[2]][2];
		
		double a = calculateDistance3d(x1,y1,z1,x2,y2,z2);
		double b = calculateDistance3d(x2,y2,z2,x3,y3,z3); 
		double c = calculateDistance3d(x3,y3,z3,x1,y1,z1); 
		double area = calculateHeron(a,b,c);
		
		return area;
	}

	
	std::vector<std::vector<float> > generateSeedsLyingOnTriangleSurfaces(std::vector<std::vector<float> > initial_mesh_vertices, std::vector<std::vector<float> > 			initial_mesh_triangles)
	{
		// naive implementation of seed generation!
		// paper isotropic remeshing yan et al,09
		// how much seeds do i need
		// using boost random because c++11 compiler does not work for unknown reasons
		// in the test file three_holes.pty/obj the seed/initial_verts ratio is about 2.5
		const float seed_initial_verts_ratio = 2.5;
		const int number_of_seeds = floor(initial_mesh_vertices.size()*seed_initial_verts_ratio);
		const int xyzs = 3;
		
		std::vector<std::vector<float> > seeds(number_of_seeds, std::vector<float>(xyzs));
		const int min_triangle_index = 0;
		std::vector<float> rnd_seed_from_triangle_surface(xyzs);
		
		// calculate the triangle areas (calculate the sum of all the weights)
		std::vector<float> triangle_areas;
		for (int i=0;i<initial_mesh_triangles.size();i++)
		{
			double triangle_area = calculateTriangleArea(initial_mesh_vertices, initial_mesh_triangles[i]);
			triangle_areas.push_back(triangle_area);
		}
		
		double total_surface_area = std::accumulate(triangle_areas.begin(), triangle_areas.end(), 0.0);
		
		// http://stackoverflow.com/questions/9294316/distribute-points-on-mesh-according-to-density

		// http://stackoverflow.com/questions/4329284/c-boost-random-numeric-generation-problem
		// http://www.bnikolic.co.uk/blog/cpp-boost-uniform01.html
		static boost::mt19937 generator(static_cast<unsigned int>(time(0)));
		// warning closed range
		//Given the parameters 1 and 6, uniform_int_distribution can can produce any of the values 1, 2, 3, 4, 5, or 6. 
		boost::random::uniform_int_distribution<> int_dist(min_triangle_index, initial_mesh_triangles.size()-1);
		boost::uniform_real<> real_dist(0.0, total_surface_area-0.01);
		boost::uniform_01<boost::random::mt19937> real_dist01(generator);
		
		for (int i=0;i<number_of_seeds;i++)
		{
			// http://stackoverflow.com/questions/1761626/weighted-random-numbers
			// pick a random number that is 0 or greater and is less than the sum of the weights
			float rnd_real_number = real_dist(generator);
			
			//go through the items one at a time, subtracting their weight from your random number,
			// until you get the item where the random number is less than that item's weight
			int rnd_triangle_index = 0;
			for(int k=0; k<initial_mesh_triangles.size(); k++) 
			{
				if(rnd_real_number < triangle_areas[k])
				{
					rnd_triangle_index = k;
				}
				rnd_real_number -= triangle_areas[k];
			}
			
			/*
			// randomly choose triangle index
			int rnd_triangle_index = int_dist(generator);
			*/
			
			double a = real_dist01();
			double b = real_dist01();
			
			// define vertices of the chosen triangle
			std::vector<float> current_triangle = initial_mesh_triangles[rnd_triangle_index];
			
			for (int j=0;j<xyzs;j++)
			{
				//http://math.stackexchange.com/questions/538458/triangle-point-picking-in-3d
				seeds[i][j] = initial_mesh_vertices[current_triangle[0]][j] + 
				a*(initial_mesh_vertices[current_triangle[1]][j] - initial_mesh_vertices[current_triangle[0]][j]) +
				b*(initial_mesh_vertices[current_triangle[2]][j] - initial_mesh_vertices[current_triangle[0]][j]);
			}
		}
		return seeds;
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
