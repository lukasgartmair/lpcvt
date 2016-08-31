        
#include "combinatorics/mesh.h"
#include "combinatorics/delaunay.h"
#include "combinatorics/RVD.h"
#include "algebra/F_Lp.h"
#include "common/line_stream.h"
#include <fstream>
#include "combinatorics/exact/RVD_predicates.h"
#include "algebra/F_Lp.h"
#include "lpcvt_functions.h"
        
// Just pass the mesh vertices to instead a point cloud
// this is a critical question how to handle this library

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
		for(unsigned int i=0; i<RVD.delaunay()->nb_vertices(); i++) 
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
	    
	int test_combinatorics(std::vector<std::vector<float> > initial_mesh_vertices, std::vector<std::vector<float> > initial_mesh_triangles) 

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
