#include <cppunit/TestSuite.h>
#include <cppunit/ui/text/TestRunner.h>
 
#include "TestLPCVT.h"


// compile: sudo g++ -I/usr/include/openvdb TestRunner.cpp -lcppunit -ltbb -o RunTests
//latest successful: sudo g++ -I/usr/include/openvdb ~/LpCVT/combinatorics/mesh.cpp TestRunner.cpp -lcppunit -ltbb -o RunTests
//sudo g++ ~/LpCVT/combinatorics/mesh.cpp TestRunner.cpp -lcppunit -o RunTests
//sudo g++ ~/LpCVT/combinatorics/mesh.cpp lpcvt_functions.cpp TestRunner.cpp -lcppunit -o RunTests
// sudo g++ -lCGAL -lgmp -I/~/LpCVT/combinatorics ~/LpCVT/combinatorics/mesh.cpp ~/LpCVT/combinatorics/delaunay.cpp ~/LpCVT/combinatorics/delaunay_CGAL.cpp ~/LpCVT/combinatorics/exact/RVD_predicates.cpp lpcvt_functions.cpp TestRunner.cpp -lgmp -lCGAL -lcppunit -o RunTests

//sudo g++ -std=c++11 -lCGAL -lgmp -I/~/LpCVT/combinatorics ~/LpCVT/combinatorics/mesh.cpp ~/LpCVT/combinatorics/delaunay.cpp ~/LpCVT/combinatorics/delaunay_CGAL.cpp ~/LpCVT/combinatorics/exact/RVD_predicates.cpp lpcvt_functions.cpp TestRunner.cpp -lgmp -lCGAL -lcppunit -o RunTests

//sudo g++ -lCGAL -lgmp -I/~/LpCVT/combinatorics ~/LpCVT/combinatorics/mesh.cpp ~/LpCVT/combinatorics/delaunay.cpp ~/LpCVT/combinatorics/delaunay_CGAL.cpp ~/LpCVT/combinatorics/exact/RVD_predicates.cpp lpcvt_functions.cpp ~/LpCVT/algebra/F_Lp.cpp TestRunner.cpp -lgmp -lCGAL -lcppunit -o RunTests

//sudo g++ -lCGAL -lgmp -I/combinatorics combinatorics/mesh.cpp combinatorics/delaunay.cpp combinatorics/delaunay_CGAL.cpp combinatorics/exact/RVD_predicates.cpp lpcvt_functions.cpp algebra/F_Lp.cpp TestRunner.cpp -lgmp -lCGAL -lcppunit -o RunTests



 
using namespace std;
 
int main() {
	CppUnit::TextUi::TestRunner runner;
 
	cout << "Creating Test Suites:" << endl;
	runner.addTest(TestMesh::suite());
	cout<< "Running the unit tests."<<endl;
	runner.run();
 
	return 0;
}
