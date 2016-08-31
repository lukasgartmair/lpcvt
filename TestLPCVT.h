
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
		std::vector<std::vector<float> > triangles(number_of_triangles, std::vector<float>(number_of_vertex_indices_per_triangle));
		for (int i=0;i<number_of_triangles;i++)
		{
			triangles[i][0] = faces[i * number_of_vertex_indices_per_triangle];
			triangles[i][1] = faces[(i * number_of_vertex_indices_per_triangle)+1];
			triangles[i][2] = faces[(i * number_of_vertex_indices_per_triangle)+2];
		}
		
		
		std::cerr << "          ========== unit test combinatorics ======" << std::endl ;
		int nov = 0;
		nov = Geex::test_combinatorics(initial_mesh_vertices, triangles);
		
		CPPUNIT_ASSERT_EQUAL(8, nov);
	}

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

};
