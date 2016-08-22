



// Geometry and combinatorics tests
//
//  Computes and output to a file:
//     Restricted Voronoi Diagram (intersection between a mesh and a 3D Voronoi diagram)
//     Restricted Delaunay Triangulation (dual of a Restricted Voronoi Diagram)
//     Clipped Voronoi Diagram (intersection between the interior of a mesh and a 3D Voronoi diagram)
//
//=================================================================================================== 

// Quick reference:
//
// Mesh: stores a Piecewise Linear Complex (+ functions to read it from a file in .obj format)
// Delaunay: abstract interface to Delaunay triangulation (+ Delaunay_CGAL: implementation with CGAL)
// RestrictedVoronoiDiagram: Given a Mesh and a Delaunay, computes the intersection between a Voronoi diagram and a Mesh. 
// ClippedVoronoiDiagram: computes the intersection between a Voronoi diagram and the interior of a closed Mesh
//
// RestrictedVoronoiDiagram and ClippedVoronoiDiagram are accessed through the for_each_triangle(do_it) function,
// that does the actual computation, and calls the user-provided function do_it(i,j,k,l) for each integration
// simplex. Note that do_it can be a function or an object that overloards operator(). 
// This mechanism can be used for both computing F-Lp and displaying the clipped Voronoi cells as in figures 4,5,6
// in the paper.
