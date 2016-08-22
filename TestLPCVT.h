
#include <iostream>
#include <cppunit/TestFixture.h>
#include <cppunit/TestAssert.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestCase.h>

#include "combinatorics/mesh.h"
 
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
				
		suiteOfTests->addTest(new CppUnit::TestCaller<TestMesh>("Test4 - Receive Triangles",
				&TestMesh::testMesh_ReceiveTriangles ));

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
	
		int dimension = 3;
		std::vector<std::vector<float> > vertices(dimension, std::vector<float>(dimension));
		// fill the points vector as in test_mesh_vertices.obj
		for (int i=0;i<dimension;i++)
		{
			vertices[i][0] = i;
			vertices[i][1] = i;
			vertices[i][2] = i;
		//std::cout << vertices[i][0] << " "  << " " << vertices[i][1] << " " << vertices[i][2] << std::endl; 
		}
		
		
		int number_of_verts = M.receiveVertices(vertices);
		CPPUNIT_ASSERT_EQUAL(dimension,number_of_verts);
	}
	
	void testMesh_ReceiveTriangles()
	{
		Geex::Mesh M ;
	
		int dimension = 10;
		std::vector<std::vector<float> > triangles(dimension, std::vector<float>(dimension));
		// fill the points vector as in test_mesh_vertices.obj
		for (int i=0;i<dimension;i++)
		{
			triangles[i][0] = 20+i;
			triangles[i][1] = 20+i+1;
			triangles[i][2] = 20+i+2;
		//std::cout << triangles[i][0] << " "  << " " << triangles[i][1] << " " << triangles[i][2] << std::endl; 
		}
		
		M.receiveTriangles(triangles);

	}
	
	

};
