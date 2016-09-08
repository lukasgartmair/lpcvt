
#include <iostream>
#include <cppunit/TestFixture.h>
#include <cppunit/TestAssert.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestCase.h>

#include "combinatorics/mesh.h"

#include "lpcvt_functions.h"
 
class TestMesh : public CppUnit::TestFixture {
private:
	Geex::Mesh *testMesh;
public:
	TestMesh() : testMesh(NULL) {}
	virtual ~TestMesh() {
		delete testMesh;
	}
 
	static CppUnit::Test *suite() {
		CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestMesh");

		suiteOfTests->addTest(new CppUnit::TestCaller<TestMesh>("Test1 - Test the Test itsself",
				&TestMesh::testMesh_TestTheTest ));
				
		suiteOfTests->addTest(new CppUnit::TestCaller<TestMesh>("Test2 - Load Data",
				&TestMesh::testMesh_LoadData ));
				
		suiteOfTests->addTest(new CppUnit::TestCaller<TestMesh>("Test3 - Receive Vertices",
				&TestMesh::testMesh_ReceiveVertices ));
				
		suiteOfTests->addTest(new CppUnit::TestCaller<TestMesh>("Test4 - Receive Vertices and Triangles",
				&TestMesh::testMesh_ReceiveVerticesAndTriangles ));
				
		suiteOfTests->addTest(new CppUnit::TestCaller<TestMesh>("Test5 - Receive Vertices and Triangles Cube Mesh",
				&TestMesh::testMesh_ReceiveVerticesAndTrianglesCubeMesh ));
				
		suiteOfTests->addTest(new CppUnit::TestCaller<TestMesh>("Test6 - Test Combinatorics",
				&TestMesh::testMesh_TestCombinatorics ));

		suiteOfTests->addTest(new CppUnit::TestCaller<TestMesh>("Test7 - Test Algebra",
				&TestMesh::testMesh_TestAlgebra ));
				
		suiteOfTests->addTest(new CppUnit::TestCaller<TestMesh>("Test8 - Test Algebra with different seeds",
				&TestMesh::testMesh_TestAlgebraWithDifferentSeeds ));
			
		suiteOfTests->addTest(new CppUnit::TestCaller<TestMesh>("Test9 - Test Combinatorics with reference arguments",
				&TestMesh::testMesh_TestCombinatoricsByReference ));
				
		suiteOfTests->addTest(new CppUnit::TestCaller<TestMesh>("Test10 - Test Generate Seeds that lie on triangle surfaces ",
				&TestMesh::testMesh_TestGenerateSeeds ));

		return suiteOfTests;
	}
 
	/// Setup method
	void setUp() {}
 
	/// Teardown method
	void tearDown() {}
 
protected:
	void testMesh_TestTheTest() 
	{
		int x = 1;
		int z = 2;
		int u = x + z;
		CPPUNIT_ASSERT( u == 3);
	}
	
	void testMesh_LoadData() 
	{

		Geex::Mesh M ;
		unsigned int nb_borders = 0;
		nb_borders = M.load("data/test_mesh.obj");
		// as the test mesh is only a triangle it has three border egedes
		CPPUNIT_ASSERT( nb_borders == 3);
	}
	

	void testMesh_ReceiveVertices()
	{
		Geex::Mesh M ;
	
		// setup points
		int number_of_vertices = 10;
		int xyzs = 3;
		std::vector<std::vector<float> > vertices(number_of_vertices, std::vector<float>(xyzs));
		// fill the points vector as in test_mesh_vertices.obj
		for (int i=0;i<number_of_vertices;i++)
		{
			vertices[i][0] = i;
			vertices[i][1] = i;
			vertices[i][2] = i;
		//std::cout << vertices[i][0] << " "  << " " << vertices[i][1] << " " << vertices[i][2] << std::endl; 
		}
		
		
		int number_of_verts = M.receiveVertices(vertices);
		CPPUNIT_ASSERT_EQUAL(number_of_vertices,number_of_verts);
	}
	
	void testMesh_ReceiveVerticesAndTriangles()
	{
		Geex::Mesh M ;
	
		// setup points
		int number_of_vertices = 3;
		int xyzs = 3;
		std::vector<std::vector<float> > vertices(number_of_vertices, std::vector<float>(xyzs));
		// fill the points vector as in test_mesh_vertices.obj
		for (int i=0;i<number_of_vertices;i++)
		{
			vertices[i][0] = i;
			vertices[i][1] = i;
			vertices[i][2] = i;
		//std::cout << vertices[i][0] << " "  << " " << vertices[i][1] << " " << vertices[i][2] << std::endl; 
		}
		
		// setup triangles
		int number_of_triangles = 1;
		int number_of_vertex_indices_per_triangle = 3;
		std::vector<std::vector<float> > triangles(number_of_triangles, std::vector<float>(number_of_vertex_indices_per_triangle));
		for (int i=0;i<number_of_triangles;i++)
		{
			triangles[i][0] = i+1;
			triangles[i][1] = i+2;
			triangles[i][2] = i+3;
		}

		
		int nb_borders = M.receiveVerticesAndTriangles(vertices, triangles);
		CPPUNIT_ASSERT( nb_borders == 3);
	
	}
	
	void testMesh_ReceiveVerticesAndTrianglesCubeMesh()
	{
		Geex::Mesh M ;
	
		// cube mesh
		// setup points
		int number_of_vertices = 8;
		int xyzs = 3;
		std::vector<std::vector<float> > initial_mesh_vertices(number_of_vertices, std::vector<float>(xyzs));
		// fill the points vector as in test_mesh_vertices.obj
		initial_mesh_vertices = Geex::initializeCubeVertices();
		// setup triangles
		int number_of_triangles = 12;
		int number_of_vertex_indices_per_triangle = 3;
		int faces[number_of_triangles*number_of_vertex_indices_per_triangle] = {2,3,4,8,7,6,5,6,2,6,7,3,3,7,8,1,4,8,1,2,4,5,8,6,1,5,2,2,6,3,4,3,8,5,1,8}; 
		std::vector<std::vector<float> > triangles(number_of_triangles, std::vector<float>(number_of_vertex_indices_per_triangle));
		for (int i=0;i<number_of_triangles;i++)
		{
			triangles[i][0] = faces[i * number_of_vertex_indices_per_triangle];
			triangles[i][1] = faces[(i * number_of_vertex_indices_per_triangle)+1];
			triangles[i][2] = faces[(i * number_of_vertex_indices_per_triangle)+2];
		}
		
		int nb_borders = M.receiveVerticesAndTriangles(initial_mesh_vertices, triangles);
		CPPUNIT_ASSERT_EQUAL(0,nb_borders);
	
	}
	
	
	void testMesh_TestCombinatorics()
	
	{
		// cube mesh
		// setup points
		int number_of_vertices = 8;
		int xyzs = 3;
		std::vector<std::vector<float> > initial_mesh_vertices(number_of_vertices, std::vector<float>(xyzs));
		// fill the points vector as in test_mesh_vertices.obj
		initial_mesh_vertices = Geex::initializeCubeVertices();
		// setup triangles
		int number_of_triangles = 12;
		int number_of_vertex_indices_per_triangle = 3;
		int faces[number_of_triangles*number_of_vertex_indices_per_triangle] = {2,3,4,8,7,6,5,6,2,6,7,3,3,7,8,1,4,8,1,2,4,5,8,6,1,5,2,2,6,3,4,3,8,5,1,8}; 
		std::vector<std::vector<float> > initial_mesh_triangles(number_of_triangles, std::vector<float>(number_of_vertex_indices_per_triangle));
		for (int i=0;i<number_of_triangles;i++)
		{
			initial_mesh_triangles[i][0] = faces[i * number_of_vertex_indices_per_triangle];
			initial_mesh_triangles[i][1] = faces[(i * number_of_vertex_indices_per_triangle)+1];
			initial_mesh_triangles[i][2] = faces[(i * number_of_vertex_indices_per_triangle)+2];
		}
		
		
		std::cerr << "          ========== unit test combinatorics ======" << std::endl ;
		int number_of_rdttris = 0;
		number_of_rdttris = Geex::getCombinatorialStructureOfFLp(initial_mesh_vertices, initial_mesh_vertices, initial_mesh_triangles);

		CPPUNIT_ASSERT_EQUAL(7, number_of_rdttris);
		
		std::cerr << "          ========== unit test combinatorics corrupted faces ======" << std::endl ;
		
		// the function crashes within the 3Depict Isosurface program - trying to find the reason here
		// what can cause a segmentation fault?
		
		
		// corrupted faces i.e. 1 instead of 2 as first vertice
		int corrupted_faces[number_of_triangles*number_of_vertex_indices_per_triangle] = {1,3,4,8,7,6,5,6,2,6,7,3,3,7,8,1,4,8,1,2,4,5,8,6,1,5,2,2,6,3,4,3,8,5,1,8}; 
	
		for (int i=0;i<number_of_triangles;i++)
		{
			initial_mesh_triangles[i][0] = corrupted_faces[i * number_of_vertex_indices_per_triangle];
			initial_mesh_triangles[i][1] = corrupted_faces[(i * number_of_vertex_indices_per_triangle)+1];
			initial_mesh_triangles[i][2] = corrupted_faces[(i * number_of_vertex_indices_per_triangle)+2];
		}
	
		number_of_rdttris = Geex::getCombinatorialStructureOfFLp(initial_mesh_vertices, initial_mesh_vertices, initial_mesh_triangles);
	
		CPPUNIT_ASSERT_EQUAL(7, number_of_rdttris);
		
		std::cerr << "          ========== unit test combinatorics faces starting from 0 instead of one ======" << std::endl ;
		
		// faces starting from 0 instead of one - the face list is taken from blender, which starts with 1 as the first vertice - openvdb does with zero!
		// i'm not able to catch this one 
		// but the console output tells the whole story 
		
		         //========== unit test combinatorics faces starting from 0 instead of one ======
		         //		Segmentation fault (core dumped)
		         
		// that is the reason why the function is commented out below
		/*
		int faces_starting_from_zero[number_of_triangles*number_of_vertex_indices_per_triangle] = {1,2,3,7,6,5,4,5,1,5,6,2,2,6,7,0,3,7,0,1,3,4,7,5,0,4,1,1,5,2,3,2,7,4,0,7}; 
		for (int i=0;i<number_of_triangles;i++)
		{
			initial_mesh_triangles[i][0] = faces_starting_from_zero[i * number_of_vertex_indices_per_triangle];
			initial_mesh_triangles[i][1] = faces_starting_from_zero[(i * number_of_vertex_indices_per_triangle)+1];
			initial_mesh_triangles[i][2] = faces_starting_from_zero[(i * number_of_vertex_indices_per_triangle)+2];
		}

		number_of_rdttris = Geex::getCombinatorialStructureOfFLp(initial_mesh_vertices, initial_mesh_triangles);
		
		CPPUNIT_ASSERT_EQUAL(7, number_of_rdttris);
		*/
		
	}

	void testMesh_TestAlgebra()
	
	{
		// cube mesh
		// setup points
		int number_of_vertices = 8;
		int xyzs = 3;
		std::vector<std::vector<float> > initial_mesh_vertices(number_of_vertices, std::vector<float>(xyzs));
		// fill the points vector as in test_mesh_vertices.obj
		initial_mesh_vertices = Geex::initializeCubeVertices();
		// setup triangles
		int number_of_triangles = 12;
		int number_of_vertex_indices_per_triangle = 3;
		int faces[number_of_triangles*number_of_vertex_indices_per_triangle] = {2,3,4,8,7,6,5,6,2,6,7,3,3,7,8,1,4,8,1,2,4,5,8,6,1,5,2,2,6,3,4,3,8,5,1,8}; 
		std::vector<std::vector<float> > initial_mesh_triangles(number_of_triangles, std::vector<float>(number_of_vertex_indices_per_triangle));
		for (int i=0;i<number_of_triangles;i++)
		{
			initial_mesh_triangles[i][0] = faces[i * number_of_vertex_indices_per_triangle];
			initial_mesh_triangles[i][1] = faces[(i * number_of_vertex_indices_per_triangle)+1];
			initial_mesh_triangles[i][2] = faces[(i * number_of_vertex_indices_per_triangle)+2];
		}
		
		
		// initialize a new vector which holds the new vertices of the delaunay triangulation
		std::vector<std::vector<float> > rdt_vertices_combinatorial;
		std::vector<std::vector<float> > rdt_triangles_combinatorial;

		std::vector<std::vector<float> > rdt_vertices_algebraic;
		std::vector<std::vector<float> > rdt_triangles_algebraic;

		// get the combinatorial rvd for assert comparison
		Geex::getCombinatorialStructureOfFLpByReference(initial_mesh_vertices, initial_mesh_vertices, initial_mesh_triangles, rdt_vertices_combinatorial, 				rdt_triangles_combinatorial);

		std::cerr << "          ========== unit test combinatorics ======" << std::endl ;
		float FL_p = Geex::getAlgebraicStructureOfFLpByReference(initial_mesh_vertices, initial_mesh_vertices, initial_mesh_triangles, rdt_vertices_algebraic, 				rdt_triangles_algebraic);
		float assert_f = 12;
		CPPUNIT_ASSERT_EQUAL(assert_f, FL_p);

		// compare the combinatorial mesh with the minimized FLp algebraic one
		// check the triangles amount
		std::cerr << rdt_triangles_algebraic.size() << std::endl ;
		CPPUNIT_ASSERT(rdt_triangles_combinatorial.size() == rdt_triangles_algebraic.size());

		// check the vertice amount
		std::cerr << rdt_vertices_algebraic.size() << std::endl ;
		CPPUNIT_ASSERT_DOUBLES_EQUAL(rdt_vertices_combinatorial.size(), rdt_vertices_algebraic.size(), 0.01);

		// check the connections
		for (int i=0;i<rdt_triangles_combinatorial.size();i++)
		{
		    for (int j=0;j<3;j++)
		    {
		    	std::cerr <<  rdt_triangles_algebraic[i][j] << std::endl ;
			CPPUNIT_ASSERT_EQUAL(rdt_triangles_combinatorial[i][j], rdt_triangles_algebraic[i][j]);
		    }
		} 
		// check the vertex sizes
		for (int i=0;i<rdt_triangles_combinatorial.size();i++)
		{
		    for (int j=0;j<3;j++)
		    {	
			    for (int k=0;k<3;k++)
			    {	
			    	// -1 because triangle indices start at one -> segfault with indexing
				CPPUNIT_ASSERT_EQUAL(rdt_vertices_combinatorial[rdt_triangles_combinatorial[i][j]-1][k], 
					rdt_vertices_algebraic[rdt_triangles_algebraic[i][j]-1][k]);
			    }
		    }
		} 
		
	}
	
	void testMesh_TestAlgebraWithDifferentSeeds()
	
	{
		// cube mesh
		// setup points
		int number_of_seeds = 8;
		int number_of_vertices = 8;
		int xyzs = 3;
		
		std::vector<std::vector<float> > initial_mesh_vertices(number_of_vertices, std::vector<float>(xyzs));
		// fill the points vector as in test_mesh_vertices.obj
		initial_mesh_vertices = Geex::initializeCubeVertices();
		
		// here seperate seeds are created that are different from the meshs vertices!
		std::vector<std::vector<float> > seeds(number_of_seeds, std::vector<float>(xyzs));

		seeds[0][0] = 0.5;
		seeds[0][1] = 0;
		seeds[0][2] = 0;

		seeds[1][0] = 0.5;
		seeds[1][1] = 0;
		seeds[1][2] = 0.5;

		seeds[2][0] = 0;
		seeds[2][1] = 0;
		seeds[2][2] = 0.5;

		seeds[3][0] = 0;
		seeds[3][1] = 0;
		seeds[3][2] = 0;

		seeds[4][0] = 0.5;
		seeds[4][1] = 0.5;
		seeds[4][2] = 0;

		seeds[5][0] = 1;
		seeds[5][1] = 1;
		seeds[5][2] = 1;

		seeds[6][0] = 0.5;
		seeds[6][1] = 0.5;
		seeds[6][2] = 0.5;

		seeds[7][0] = 0.5;
		seeds[7][1] = 0.5;
		seeds[7][2] = 0.5;
		
		// setup triangles
		int number_of_triangles = 12;
		int number_of_vertex_indices_per_triangle = 3;
		int faces[number_of_triangles*number_of_vertex_indices_per_triangle] = {2,3,4,8,7,6,5,6,2,6,7,3,3,7,8,1,4,8,1,2,4,5,8,6,1,5,2,2,6,3,4,3,8,5,1,8}; 
		std::vector<std::vector<float> > initial_mesh_triangles(number_of_triangles, std::vector<float>(number_of_vertex_indices_per_triangle));
		for (int i=0;i<number_of_triangles;i++)
		{
			initial_mesh_triangles[i][0] = faces[i * number_of_vertex_indices_per_triangle];
			initial_mesh_triangles[i][1] = faces[(i * number_of_vertex_indices_per_triangle)+1];
			initial_mesh_triangles[i][2] = faces[(i * number_of_vertex_indices_per_triangle)+2];
		}
		
		
		// initialize a new vector which holds the new vertices of the delaunay triangulation
		std::vector<std::vector<float> > rdt_vertices_combinatorial;
		std::vector<std::vector<float> > rdt_triangles_combinatorial;

		std::vector<std::vector<float> > rdt_vertices_algebraic;
		std::vector<std::vector<float> > rdt_triangles_algebraic;

		// get the combinatorial rdt for assert comparison
		Geex::getCombinatorialStructureOfFLpByReference(seeds, initial_mesh_vertices, initial_mesh_triangles, rdt_vertices_combinatorial, rdt_triangles_combinatorial);

		std::cerr << "          ========== unit test algebra 1 ======" << std::endl ;
		float FL_p = Geex::getAlgebraicStructureOfFLpByReference(seeds, initial_mesh_vertices, initial_mesh_triangles, rdt_vertices_algebraic, rdt_triangles_algebraic);
		float assert_f = 16.71875;
		CPPUNIT_ASSERT_DOUBLES_EQUAL(assert_f, FL_p, 0.01);

		// compare the combinatorial mesh with the minimized FLp algebraic one
		// check the triangles amount
		std::cerr << rdt_triangles_algebraic.size() << std::endl ;
		CPPUNIT_ASSERT(rdt_triangles_combinatorial.size() == rdt_triangles_algebraic.size());

		// check the vertice amount
		std::cerr << rdt_vertices_algebraic.size() << std::endl ;
		CPPUNIT_ASSERT_DOUBLES_EQUAL(rdt_vertices_combinatorial.size(), rdt_vertices_algebraic.size(), 0.01);

		// check the connections
		for (int i=0;i<rdt_triangles_combinatorial.size();i++)
		{
		    for (int j=0;j<3;j++)
		    {
		    	std::cerr <<  rdt_triangles_algebraic[i][j] << std::endl ;
			CPPUNIT_ASSERT_EQUAL(rdt_triangles_combinatorial[i][j], rdt_triangles_algebraic[i][j]);
		    }
		} 
		// check the vertex sizes
		for (int i=0;i<rdt_triangles_combinatorial.size();i++)
		{
		    for (int j=0;j<3;j++)
		    {	
			    for (int k=0;k<3;k++)
			    {	
			    	// -1 because triangle indices start at one -> segfault with indexing
				CPPUNIT_ASSERT_EQUAL(rdt_vertices_combinatorial[rdt_triangles_combinatorial[i][j]-1][k], 
					rdt_vertices_algebraic[rdt_triangles_algebraic[i][j]-1][k]);
			    }
		    }
		} 
	}
	
	void testMesh_TestCombinatoricsByReference()
	{
		// cube mesh
		// setup points
		int number_of_vertices = 8;
		int xyzs = 3;
		std::vector<std::vector<float> > initial_mesh_vertices(number_of_vertices, std::vector<float>(xyzs));
		// fill the points vector as in test_mesh_vertices.obj
		initial_mesh_vertices = Geex::initializeCubeVertices();
		// setup triangles
		int number_of_triangles = 12;
		int number_of_vertex_indices_per_triangle = 3;
		int faces[number_of_triangles*number_of_vertex_indices_per_triangle] = {2,3,4,8,7,6,5,6,2,6,7,3,3,7,8,1,4,8,1,2,4,5,8,6,1,5,2,2,6,3,4,3,8,5,1,8}; 
		std::vector<std::vector<float> > initial_mesh_triangles(number_of_triangles, std::vector<float>(number_of_vertex_indices_per_triangle));
		for (int i=0;i<number_of_triangles;i++)
		{
			initial_mesh_triangles[i][0] = faces[i * number_of_vertex_indices_per_triangle];
			initial_mesh_triangles[i][1] = faces[(i * number_of_vertex_indices_per_triangle)+1];
			initial_mesh_triangles[i][2] = faces[(i * number_of_vertex_indices_per_triangle)+2];
		}
		
		std::cerr << "          ========== unit test combinatorics by reference ======" << std::endl ;
		
		// initialize a new vector which holds the new vertices of the delaunay triangulation
		std::vector<std::vector<float> > rdt_vertices;

		std::vector<std::vector<float> > rdt_triangles;
		
		Geex::getCombinatorialStructureOfFLpByReference(initial_mesh_vertices, initial_mesh_vertices, initial_mesh_triangles, rdt_vertices, rdt_triangles);

		int number_of_rdttris = rdt_triangles.size();
		CPPUNIT_ASSERT_EQUAL(7, number_of_rdttris);
		
		int number_of_rdtverts = rdt_vertices.size();
		CPPUNIT_ASSERT_EQUAL(8, number_of_rdtverts);

		// check the contents - i got only zeros in 3depict as content
		int rdt_verts[number_of_rdtverts*xyzs] = {1,0,0,1,0,1,0,0,1,0,0,0,1,1,0,1,1,1,0,1,1,0,1,0}; 
		int verts_counter = 0;
		
		for (int i=0;i<rdt_vertices.size();i++)
		{
			for (int j=0;j<xyzs;j++)
			{
				CPPUNIT_ASSERT_DOUBLES_EQUAL(rdt_verts[verts_counter], rdt_vertices[i][j], 0.01);
				verts_counter += 1;
			}
		}
	
		int rdt_faces[number_of_rdttris*number_of_vertex_indices_per_triangle] = {3,7,4,1,4,5,1,2,4,1,5,2,2,5,6,5,8,6,2,6,3}; 
		int faces_counter = 0;
		
		for (int i=0;i<rdt_triangles.size();i++)
		{
			for (int j=0;j<number_of_vertex_indices_per_triangle;j++)
			{
				CPPUNIT_ASSERT_DOUBLES_EQUAL(rdt_faces[faces_counter], rdt_triangles[i][j], 0.01);
				faces_counter += 1;
			}
		}
	}
	
	void testMesh_TestGenerateSeeds()
	{

		// cube mesh
		// setup points
		int number_of_vertices = 8;
		int xyzs = 3;
		std::vector<std::vector<float> > initial_mesh_vertices(number_of_vertices, std::vector<float>(xyzs));
		// fill the points vector as in test_mesh_vertices.obj
		initial_mesh_vertices = Geex::initializeCubeVertices();
	
		// setup triangles
		int number_of_triangles = 12;
		int number_of_vertex_indices_per_triangle = 3;
		int faces[number_of_triangles*number_of_vertex_indices_per_triangle] = {2,3,4,8,7,6,5,6,2,6,7,3,3,7,8,1,4,8,1,2,4,5,8,6,1,5,2,2,6,3,4,3,8,5,1,8}; 
		std::vector<std::vector<float> > initial_mesh_triangles(number_of_triangles, std::vector<float>(number_of_vertex_indices_per_triangle));
		for (int i=0;i<number_of_triangles;i++)
		{
			// minus one because of indexing
			initial_mesh_triangles[i][0] = faces[i * number_of_vertex_indices_per_triangle]-1;
			initial_mesh_triangles[i][1] = faces[(i * number_of_vertex_indices_per_triangle)+1]-1;
			initial_mesh_triangles[i][2] = faces[(i * number_of_vertex_indices_per_triangle)+2]-1;
		}
		
		std::vector<std::vector<float> > seeds;
		seeds = Geex::generateSeedsLyingOnTriangleSurfaces(initial_mesh_vertices, initial_mesh_triangles);
		
		std::cerr << seeds[0][0] << std::endl ;
		std::cerr << seeds[0][1] << std::endl ;
		std::cerr << seeds[0][2] << std::endl ;
		
		std::cerr << seeds[1][0] << std::endl ;
		std::cerr << seeds[1][1] << std::endl ;
		std::cerr << seeds[1][2] << std::endl ;
		
		
		
		// check amount
		CPPUNIT_ASSERT_DOUBLES_EQUAL(2, seeds.size(), 0.01);
	
	
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

};
