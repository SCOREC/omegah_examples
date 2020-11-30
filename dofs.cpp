#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_for.hpp>
//#include <Omega_h_adj.hpp>
//#include "Omega_h_element.hpp"
#include "Omega_h_shape.hpp"

using namespace	std;
using namespace	Omega_h;

// Finding the centroid of each mesh using parallel for loops
Read<Real> findCentroid(Mesh& mesh, const int& dim) {
  LO numElems = mesh.nelems();
  auto elems2verts = mesh.ask_elem_verts();
  auto coords = mesh.coords();
  Write<Real> centroids_w(dim*numElems,"Element Centroids");
  
  // 1D-Edges
  auto lamb1D = OMEGA_H_LAMBDA(LO i) {
    auto vertIndices = gather_verts<2>(elems2verts,i);
    auto vertCoords = gather_vectors<2,1>(coords, vertIndices);
    for (int j = 0; j < 2; j++) {
      centroids_w[dim*i+0] += vertCoords[j][0]/2;  
    }
  };
  // 2D-Triangles
  auto lamb2D = OMEGA_H_LAMBDA(LO i) { 
    auto vertIndices = gather_verts<3>(elems2verts,i);
    auto vertCoords = gather_vectors<3,2>(coords, vertIndices);
    for (int j = 0; j < 3; j++) {
      centroids_w[dim*i+0] += vertCoords[j][0]/3;
      centroids_w[dim*i+1] += vertCoords[j][1]/3;
    } 
  };
  // 3D-Tetrahedrons
  auto lamb3D = OMEGA_H_LAMBDA(LO i) {
    auto vertIndices = gather_verts<4>(elems2verts,i);
    auto vertCoords = gather_vectors<4,3>(coords, vertIndices);
    for (int j = 0; j < 4; j++) {
      centroids_w[dim*i+0] += vertCoords[j][0]/4;
      centroids_w[dim*i+1] += vertCoords[j][1]/4;
      centroids_w[dim*i+2] += vertCoords[j][2]/4;
    } 
  };

  if (dim == 1) {
    assert(numElems == mesh.nedges());
    parallel_for(numElems,lamb1D,"Element Centroid Calc Edges");
  }
  else if (dim == 2) {
    assert(numElems == mesh.nfaces());
    parallel_for(numElems,lamb2D,"Element Centroid Calc Triangles");
  }
  else if (dim == 3) {
    assert(numElems == mesh.nregions());
    parallel_for(numElems,lamb3D,"Element Centroid Calc Tetrahedrons");
  }
  else {
    fprintf(stderr, "Invalid Dimension of Mesh\n");
    exit(0);
  } 
  Read<Real> centroids(centroids_w);
  return centroids;
}

// Finding the Midpoint of each Edge
Read<Real> findEdgeMidpoint(Mesh& mesh,const int& dim) {
  LO numElems = mesh.nelems();
  LO numEdges = mesh.nedges();
  auto coords = mesh.coords();
  auto edges2verts = mesh.ask_verts_of(EDGE);

  Write<Real> midpoints_w(dim*numEdges,"Edge Midpoints");  
  // 2D-Triangles 
  auto lamb2D = OMEGA_H_LAMBDA(LO i) {
    auto vertIndices = gather_verts<2>(edges2verts,i);
    auto vertCoords = gather_vectors<2,2>(coords, vertIndices);
    for (int j = 0; j < 2; j++) {
      midpoints_w[dim*i+0] += vertCoords[j][0]/2;
      midpoints_w[dim*i+1] += vertCoords[j][1]/2;
    }
  };
  // 3D-Tetrahedron
  auto lamb3D = OMEGA_H_LAMBDA(LO i) {
    auto vertIndices = gather_verts<2>(edges2verts,i);
    auto vertCoords = gather_vectors<2,3>(coords, vertIndices);
    for (int j = 0; j < 2; j++) {
      midpoints_w[dim*i+0] += vertCoords[j][0]/2;
      midpoints_w[dim*i+1] += vertCoords[j][1]/2;
      midpoints_w[dim*i+2] += vertCoords[j][2]/2;
    } 
  };
  if (dim == 2) {
    assert(numElems == mesh.nfaces());
    parallel_for(numEdges,lamb2D,"Edge Midpoint Calc Triangles");
  }
  else if (dim == 3) {
    assert(numElems == mesh.nregions());
    parallel_for(numEdges,lamb3D,"Edge Midpoint Calc Tetrahedrons");
  }
  else {
    fprintf(stderr, "Invalid Dimension of Mesh\n");
    exit(0);
  } 
  
  Read<Real> midpoints(midpoints_w);
  return midpoints;
}

// Finding the centroid of each mesh using parallel for loops
Read<Real> findElementEdgeMidpoint(Mesh& mesh, const int& dim, Read<Real>& midpointEdge) {
  LO numElems = mesh.nelems();
  auto elems2edges = mesh.ask_down(dim,1).ab2b;
  Write<Real> avgMidpoint_w(dim*numElems,"Element Average Midpoint");  
  // 2D-Triangles
  auto lamb2D = OMEGA_H_LAMBDA(LO i) {
    auto vertIndices = gather_verts<3>(elems2edges,i);
    auto midpointCoords = gather_vectors<3,2>(midpointEdge, vertIndices);
    for (int j = 0; j < 3; j++) {
      avgMidpoint_w[dim*i+0] += midpointCoords[j][0]/3;
      avgMidpoint_w[dim*i+1] += midpointCoords[j][1]/3;
    }
  };
  // 3D-Tetrahedron
  auto lamb3D = OMEGA_H_LAMBDA(LO i) {
    auto vertIndices = gather_verts<6>(elems2edges,i);
    auto midpointCoords = gather_vectors<6,3>(midpointEdge, vertIndices);
    for (int j = 0; j < 6; j++) {
      avgMidpoint_w[dim*i+0] += midpointCoords[j][0]/6;
      avgMidpoint_w[dim*i+1] += midpointCoords[j][1]/6;
      avgMidpoint_w[dim*i+2] += midpointCoords[j][2]/6;
    }
  };
  if (dim == 2) {
    assert(numElems == mesh.nfaces());
    parallel_for(numElems,lamb2D,"Element Centroid Calc Triangles");
  }
  else if (dim == 3) {
    assert(numElems == mesh.nregions());
    parallel_for(numElems,lamb3D,"Element Centroid Calc Tetrahedrons");
  }
  else {
    fprintf(stderr, "Invalid Dimension of Mesh\n");
    exit(0);
  } 
  
  Read<Real> avgMidpoint(avgMidpoint_w);
  return avgMidpoint;
}


int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  if(argc!=2) {
    fprintf(stderr, "Usage: %s <input mesh>\n", argv[0]);
    return 0;
  }
  const auto rank = lib.world()->rank();
  const auto inmesh = argv[1];
  Mesh mesh(&lib);
  binary::read(inmesh, lib.world(), &mesh);
  const auto dim = mesh.dim();
  
  Read<Real> centroids = findCentroid(mesh,dim);
  int order = 2;
  if (mesh.nelems() != mesh.nfaces()) {
    order = 3;
  } 

  Read<Real> midpointEdge = findEdgeMidpoint(mesh,dim);
  mesh.add_tag<Real>(EDGE,"Edge Midpoints",dim);
  mesh.set_tag<Real>(EDGE,"Edge Midpoints",midpointEdge);
  
  Read<Real> averageMidpointElem = findElementEdgeMidpoint(mesh,dim,midpointEdge);
  mesh.add_tag<Real>(order,"Average Midpoint per Element",dim);
  mesh.set_tag<Real>(order,"Average Midpoint per Element",averageMidpointElem);
  
  binary::write("./tag.osh",&mesh);

  return 0;
}
