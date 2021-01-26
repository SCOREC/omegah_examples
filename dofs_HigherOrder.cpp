#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_for.hpp>
//#include <Omega_h_adj.hpp>
//#include "Omega_h_element.hpp"
#include "Omega_h_shape.hpp"
#include <math.h>
#include "adjacency.cpp"


using namespace	std;
using namespace	Omega_h;

OMEGA_H_INLINE 
void v(const double* coord, const LO dim, double value[], const LO order) {
  if (dim == 2) {
    value[0] = sin(coord[0])/order+cos(coord[1])*order;
    value[1] = sin(coord[1])/order+cos(coord[0])*order;
  }
  else {
    value[0] = sin(coord[0])/order+cos(coord[1])*order-sin(coord[2])*cos(coord[2]);
    value[1] = sin(coord[1])/order+cos(coord[2])*order-sin(coord[0])*cos(coord[0]);
    value[2] = sin(coord[2])/order+cos(coord[0])*order-sin(coord[1])*cos(coord[1]);
  }
}

Read<Real> createVerts(Mesh& mesh, LO dim, LO order) {
  LO numElems = mesh.nelems();
  LO numVerts = mesh.nverts();
  auto coords = mesh.coords();
  LO numNodes = 1;
  Write<Real> vertDOFs_w(numVerts*numNodes*dim,"Vertex Node DOFs");
  auto lamb2D = OMEGA_H_LAMBDA(LO m) {
    Real coord[2];
    coord[0] = coords[dim*m+0];
    coord[1] = coords[dim*m+1];
    Real value[2];
    v(coord,dim,value,order);
    vertDOFs_w[dim*m+0] = value[0];
    vertDOFs_w[dim*m+1] = value[1];
  };
  auto lamb3D = OMEGA_H_LAMBDA(LO m) {
    Real coord[3];
    coord[0] = coords[dim*m+0];
    coord[1] = coords[dim*m+1];
    coord[2] = coords[dim*m+2];
    Real value[3];
    v(coord,dim,value,order);
    vertDOFs_w[dim*m+0] = value[0];
    vertDOFs_w[dim*m+1] = value[1];
    vertDOFs_w[dim*m+2] = value[2];
  };
  if (dim == 2) {
    assert(numElems == mesh.nfaces());
    parallel_for(numVerts,lamb2D,"Equally Spaced 2D Vertex DOFs");
  }
  else if (dim == 3) {
    assert(numElems == mesh.nregions());
    parallel_for(numVerts,lamb3D,"Equally Spaced 3D Vertex DOFs");
  }
  else {
    fprintf(stderr, "Invalid Dimension of Mesh\n");
    exit(0);
  } 
  
  Read<Real> vertDOFs(vertDOFs_w);
  return vertDOFs;
}

Read<Real> createEdges(Mesh& mesh, LO dim, LO order) {
  LO numElems = mesh.nelems();
  LO numEdges = mesh.nedges();
  auto edges2verts = mesh.ask_verts_of(EDGE);
  auto coords = mesh.coords();
  LO numNodes = order-1;
  Write<Real> edgeDOFs_w(numEdges*numNodes*dim,"Edge Node DOFs");
  auto lamb2D = OMEGA_H_LAMBDA(LO m) {
    auto vertIndices = gather_verts<2>(edges2verts,m);
    auto vertCoord = gather_vectors<2,2>(coords,vertIndices);
    LO n = 0;
    for (int i = order-1; i > 0; i--,n++) {
      int j = order-i;
      Real coord[2];
      coord[0] = (vertCoord[0][0]*i+vertCoord[1][0]*j)/order;
      coord[1] = (vertCoord[0][1]*i+vertCoord[1][1]*j)/order;
      Real value[2];
      v(coord,dim,value,order);
      edgeDOFs_w[dim*(m*numNodes+n)+0] = value[0]; 
      edgeDOFs_w[dim*(m*numNodes+n)+1] = value[1];
    }
  };
  auto lamb3D = OMEGA_H_LAMBDA(LO m) {
    auto vertIndices = gather_verts<2>(edges2verts,m);
    auto vertCoord = gather_vectors<2,3>(coords,vertIndices);
    LO n = 0;
    for (int i = order-1; i > 0; i--,n++) {
      int j = order-i;
      Real coord[3];
      coord[0] = (vertCoord[0][0]*i+vertCoord[1][0]*j)/order;
      coord[1] = (vertCoord[0][1]*i+vertCoord[1][1]*j)/order;
      coord[2] = (vertCoord[0][2]*i+vertCoord[1][2]*j)/order;
      Real value[3];
      v(coord,dim,value,order);
      edgeDOFs_w[dim*(m*numNodes+n)+0] = value[0]; 
      edgeDOFs_w[dim*(m*numNodes+n)+1] = value[1];
      edgeDOFs_w[dim*(m*numNodes+n)+2] = value[2];
    }
  };
  if (dim == 2) {
    assert(numElems == mesh.nfaces());
    parallel_for(numEdges,lamb2D,"Equally Spaced 2D Edge DOFs");
  }
  else if (dim == 3) {
    assert(numElems == mesh.nregions());
    parallel_for(numEdges,lamb3D,"Equally Spaced 3D Edge DOFs");
  }
  else {
    fprintf(stderr, "Invalid Dimension of Mesh\n");
    exit(0);
  } 
  
  Read<Real> edgeDOFs(edgeDOFs_w);
  return edgeDOFs;
}

Read<Real> createTriangles(Mesh& mesh, LO dim, LO order) {
  LO numElems = mesh.nelems();
  LO numTrigs = mesh.nfaces();
  auto trigs2verts = mesh.ask_verts_of(2);
  auto coords = mesh.coords();
  LO numNodes = (order-1)*(order-2)/2;
  Write<Real> trigDOFs_w(numTrigs*numNodes*dim,"Trianlge Node DOFs");

  auto lamb2D = OMEGA_H_LAMBDA(LO m) {
    auto vertIndices = gather_verts<3>(trigs2verts,m);
    auto vertCoord = gather_vectors<3,2>(coords,vertIndices);
    LO n = 0;
    for (int i = order-2; i > 0; i--) {
      for (int j = order-i-1; j > 0; j--,n++) {
        int k = order-i-j;
        Real coord[2];
        coord[0] = (vertCoord[0][0]*i+vertCoord[1][0]*j+
                     vertCoord[2][0]*k)/order;
        coord[1] = (vertCoord[0][1]*i+vertCoord[1][1]*j+
                     vertCoord[2][1]*k)/order;
        Real value[2];
        v(coord,dim,value,order);
        trigDOFs_w[dim*(m*numNodes+n)+0] = value[0];
        trigDOFs_w[dim*(m*numNodes+n)+1] = value[1];
      }
    }
  };
  auto lamb3D = OMEGA_H_LAMBDA(LO m) {
    auto vertIndices = gather_verts<3>(trigs2verts,m);
    auto vertCoord = gather_vectors<3,3>(coords,vertIndices);
    LO n = 0;
    for (int i = order-2; i > 0; i--) {
      for (int j = order-i-1; j > 0; j--,n++) {
        int k = order-i-j;
        Real coord[3];
        coord[0] = (vertCoord[0][0]*i+vertCoord[1][0]*j+
                     vertCoord[2][0]*k)/order;
        coord[1] = (vertCoord[0][1]*i+vertCoord[1][1]*j+
                     vertCoord[2][1]*k)/order;
        coord[2] = (vertCoord[0][2]*i+vertCoord[1][2]*j+
                     vertCoord[2][2]*k)/order;        
        Real value[3];
        v(coord,dim,value,order);
        trigDOFs_w[dim*(m*numNodes+n)+0] = value[0];
        trigDOFs_w[dim*(m*numNodes+n)+1] = value[1];
        trigDOFs_w[dim*(m*numNodes+n)+2] = value[2];
      }
    }
  };
  if (dim == 2) {
    assert(numElems == mesh.nfaces());
    parallel_for(numTrigs,lamb2D,"Equally Spaced 2D Triangle DOFs");
  }
  else if (dim == 3) {
    assert(numElems == mesh.nregions());
    parallel_for(numTrigs,lamb3D,"Equally Spaced 3D Triangle DOFs");
  }
  else {
    fprintf(stderr, "Invalid Dimension of Mesh\n");
    exit(0);
  } 
  
  Read<Real> trigDOFs(trigDOFs_w);
  return trigDOFs;
}

Read<Real> createTetrahedrons(Mesh& mesh, LO dim, LO order) {
  LO numElems = mesh.nelems();
  LO numTets = mesh.nregions();
  auto tets2verts = mesh.ask_verts_of(3);
  auto coords = mesh.coords();
  LO numNodes = (order-1)*(order-2)*(order-3)/6;
  Write<Real> tetDOFs_w(numTets*numNodes*dim,"Tetrahedron Node DOFs");

  auto lamb3D = OMEGA_H_LAMBDA(LO m) {
    auto vertIndices = gather_verts<4>(tets2verts,m);
    auto vertCoord = gather_vectors<4,3>(coords,vertIndices);
    LO n = 0;
    for (int i = order-3; i > 0; i--) {
      for (int j = order-i-2; j > 0; j--) {
        for (int k = order-i-j-1; k > 0; k--,n++) {
        int l = order-i-j-k;
        Real coord[3];
        coord[0] = (vertCoord[0][0]*i+vertCoord[1][0]*j+
                     vertCoord[2][0]*k+vertCoord[3][0]*l)/order;
        coord[1] = (vertCoord[0][1]*i+vertCoord[1][1]*j+
                     vertCoord[2][1]*k+vertCoord[3][1]*l)/order;
        coord[2] = (vertCoord[0][2]*i+vertCoord[1][2]*j+
                     vertCoord[2][2]*k+vertCoord[3][2]*l)/order;        
        Real value[3];
        v(coord,dim,value,order);
        tetDOFs_w[dim*(m*numNodes+n)+0] = value[0];
        tetDOFs_w[dim*(m*numNodes+n)+1] = value[1];
        tetDOFs_w[dim*(m*numNodes+n)+2] = value[2];
        }
      }
    }
  }; 
  if (dim == 3) {
    assert(numElems == mesh.nregions());
    parallel_for(numElems,lamb3D,"Equally Spaced 3D Triangle DOFs");
  }
  else {
    fprintf(stderr, "Invalid Dimension of Mesh\n");
    exit(0);
  } 
  
  Read<Real> tetDOFs(tetDOFs_w);
  return tetDOFs;
}
/*
Read<Real> interpolateMiddle(Mesh& mesh, LO dim, LO order, Read<Real>& vertD,
Read<Real>& edgeD, Read<Real>& trigD, Read<Real>& tetD) {
  LO numElems = mesh.nelems();
  auto elem2verts = mesh.ask_verts_of(dim);
  auto coords = mesh.coords();
  auto elem2edges = mesh.ask_adj(dim,2).ab2b;
  auto edge2verts = mesh.ask_verts_of(2);
  Write<Real> elemInter_r(dim*numElems,"Interpolation of Center Value");
  auto lamb2D = OMEGA_H_LAMBDA(LO m) {
    Real vertDOFs[3*dim];
    int n = 0;
    for (int i = 0; i < 3; i++) {
      for (int d = 0; d < 2; d++,n++) {
        vertDOFs[n] = vertD[2*elem2verts[3*m+i]+d];
      }
    }
    Real edgeDOFs[3*dim*(order-1)];
    n = 0;
    LO elemAdj[2];
    LO edgeAdj[2];
    




    for (int i = 0; i < 3; i++) {
      LO elemAdj[2];
      LO edgeAdj[2];
      for (int j = 0; j < 2; j++) {
        edgeAdj[j] = edge2verts[2*elem2edges[3*m+i]+j];
      }
      for (int j = 0; j < order-1; j++) {
        for (int d = 0; d < 2; d++) {
          edgeDOFs[]
        }
      }
    }
    auto vertIndices = gather_verts<3>(elem2verts,m);
    auto vertDOFs = gather_vectors<3,2>(vertD,vertIndices);
    auto edgeIndices = gather_verts<3>(elem2edges,m);
    auto edgeDOFs = gather_vertors<3>
  }

}*/

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  if(argc!=3) {
    fprintf(stderr, "Usage: %s <input mesh> <polynomial mesh order>\n", argv[0]);
    return 0;
  }
  const auto rank = lib.world()->rank();
  const auto inmesh = argv[1];
  const auto order = atoi(argv[2]);
  Mesh mesh(&lib);
  binary::read(inmesh, lib.world(), &mesh);
  const auto dim = mesh.dim();
 
  Read<Real> vertDOFs = createVerts(mesh,dim,order); 
  Read<Real> edgeDOFs = createEdges(mesh,dim,order);
  Read<Real> trigDOFs = createTriangles(mesh,dim,order);
  Read<Real> tetsDOFs;
  mesh.add_tag<Real>(0,"Vertex DOFs",dim);
  mesh.set_tag<Real>(0,"Vertex DOFs",vertDOFs);
  mesh.add_tag<Real>(1,"Edge DOFs",dim*(order-1));
  mesh.set_tag<Real>(1,"Edge DOFs",edgeDOFs);
  mesh.add_tag<Real>(2,"Triangle DOFs",dim*(order-1)*(order-2)/2);
  mesh.set_tag<Real>(2,"Triangle DOFs",trigDOFs);
  if (dim == 3) {
    tetsDOFs = createTetrahedrons(mesh,dim,order);
    mesh.add_tag<Real>(3,"Tetrahedron DOFs",dim*(order-1)*(order-2)*(order-3)/6);
    mesh.set_tag<Real>(3,"Tetrahedron DOFs",tetsDOFs);
  }
  binary::write("./tag.osh",&mesh);
  
  //*// Manual Testing
  HostRead<Real> vertDOFsHost(vertDOFs);
  HostRead<Real> edgeDOFsHost(edgeDOFs);
  HostRead<Real> trigDOFsHost(trigDOFs);
  
  auto coords = mesh.coords();
  HostRead<Real> coordHost(coords);
  for (int i = 0; i < mesh.nverts(); i++) {
    std::cout << i << ": ";
    for (int d = 0; d < dim; d++) {
      std::cout << coordHost[dim*i+d] << " ";
    }
    std::cout << std::endl << "\t";
    for (int d = 0; d < dim; d++) {
      std::cout << vertDOFsHost[2*i+d] << " "; 
    }
    std::cout << std::endl;
  }
  
  auto edges2verts = mesh.ask_verts_of(EDGE);
  LO numEdges = mesh.nedges();
  for (int i = 0; i < numEdges; i++) {
    std::cout << i << ": ";
    for (int d = 0; d < 2; d++) {
      std::cout << edges2verts.get(2*i+d) << " ";
    }
    std::cout << std::endl << "\t";
    for (int j = 0; j < (order-1); j++) {
      for (int d = 0; d < dim; d++) {
        std::cout << edgeDOFsHost[dim*((order-1)*i+j)+d] << " ";
      }
      std::cout << std::endl << "\t";
    }
    std::cout << std::endl;
  }

  auto trigs2verts = mesh.ask_verts_of(2);
  LO numTrigs = mesh.nfaces();
  for (int i = 0; i < numTrigs; i++) {
    std::cout << i << ": ";
    for (int d = 0; d < 3; d++) {
      std::cout << trigs2verts.get(3*i+d) << " ";
    }
    std::cout << std::endl << "\t";
    for (int j = 0; j < (order-1)*(order-2)/2; j++) {
      for (int d = 0; d < dim; d++) {
        std::cout << trigDOFsHost[dim*((order-1)*(order-2)/2*i+j)+d] << " ";
      }
      std::cout << std::endl << "\t";
    }
    std::cout << std::endl;
  }
  if (dim == 3) {
    HostRead<Real> tetsDOFsHost(tetsDOFs);
    auto tets2verts = mesh.ask_verts_of(3);
    LO numTets = mesh.nregions();
    for (int i = 0; i < numTets; i++) {
      std::cout << i << ": ";
      for (int d = 0; d < 4; d++) {
        std::cout << tets2verts.get(4*i+d) << " ";
      }
      std::cout << std::endl << "\t";
      LO numNodes = (order-1)*(order-2)*(order-3)/6;
      for (int j = 0; j < numNodes; j++) {
        for (int d = 0; d < dim; d++) {
          std::cout << tetsDOFsHost[dim*(numNodes*i+j)+d] << " ";
        }
        std::cout << std::endl << "\t";
      }
      std::cout << std::endl;
    }
  }
  // */
  return 0;
}
