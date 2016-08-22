#include <cppunit/TestSuite.h>
#include <cppunit/ui/text/TestRunner.h>
 
#include "TestLPCVT.h"


// compile: sudo g++ -I/usr/include/openvdb TestRunner.cpp -lcppunit -ltbb -o RunTests
//latest successful: sudo g++ -I/usr/include/openvdb ~/LpCVT/combinatorics/mesh.cpp TestRunner.cpp -lcppunit -ltbb -o RunTests
//sudo g++ ~/LpCVT/combinatorics/mesh.cpp TestRunner.cpp -lcppunit -o RunTests
 
using namespace std;
 
int main() {
	CppUnit::TextUi::TestRunner runner;
 
	cout << "Creating Test Suites:" << endl;
	runner.addTest(TestMesh::suite());
	cout<< "Running the unit tests."<<endl;
	runner.run();
 
	return 0;
}
