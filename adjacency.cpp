#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_for.hpp>
//#include <Omega_h_adj.hpp>
//#include "Omega_h_element.hpp"
#include "Omega_h_shape.hpp"
#include <math.h>

using namespace std;
using namespace Omega_h;

OMEGA_H_INLINE
void adjEdges(Real* edgeNodes, LO* mesh, LO* adj, LO order, LO dim, Real* out) { 
  int numNodes = order-1;
  if (mesh[0] == adj[0]) {
    for (int i = 0; i < numNodes; i++) {
      for (int d = 0; d < dim; d++) {
        out[i] = edgeNodes[dim*i+d];
      }
    }
  }
  else if (mesh[0] == adj[1]) {
    for (int i = 0; i < numNodes; i++) {
      for (int d = 0; d < dim; d++) {
        out[i] = edgeNodes[dim*(numNodes-i-1)+d];
      }
    }
  }
}

OMEGA_H_INLINE
void adjTrigs(Real* trigNodes, LO* mesh, LO* adj, LO order, LO dim, Real* out) {
  LO n = 0;
  if (mesh[0] == adj[0]) {
    if (mesh[1] == adj[1]) {
      for (int i = order-2; i > 0; i--) {
        for (int j = order-i-1; j > 0; j--,n++) {
          int p = order-i-1;
          for (int d = 0; d < dim; d++) {
            out[dim*(p*(p+1)/2-j)+d] = trigNodes[dim*n+d];
          }
        }
      }
    }
    else if(mesh[2] == adj[1]) {
      for (int i = order-2; i > 0; i--) {
        for (int j = order-i-1; j > 0; j--,n++) {
          int k = order-i-j;
          int p = order-i-1;
          for (int d = 0; d < dim; d++) {
            out[dim*(p*(p+1)/2-k)+d] = trigNodes[dim*n+d];
          }
        }
      }  
    }
  }
  else if (mesh[1] == adj[0]) {
    if (mesh[2] == adj[1]) {
      for (int i = order-2; i > 0; i--) {
        for (int j = order-i-1; j > 0; j--,n++) {
          int k = order-i-j;
          int p = order-j-1;
          for (int d = 0; d < dim; d++) {
            out[dim*(p*(p+1)/2-k)+d] = trigNodes[dim*n+d];
          }
        }
      }
    }
    else if(mesh[0] == adj[1]) {
      for (int i = order-2; i > 0; i--) {
        for (int j = order-i-1; j > 0; j--,n++) {
          int p = order-j-1;
          for (int d = 0; d < dim; d++) {
            out[dim*(p*(p+1)/2-i)+d] = trigNodes[dim*n+d];
          }
        }
      }
    }
  }
  else if (mesh[2] == adj[0]) {
    if (mesh[0] == adj[1]) {
      for (int i = order-2; i > 0; i--) {
        for (int j = order-i-1; j > 0; j--,n++) {
          int k = order-i-j;
          int p = order-k-1;
          for (int d = 0; d < dim; d++) {
            out[dim*(p*(p+1)/2-i)+d] = trigNodes[dim*n+d];
          }
        }
      }
    }
    else if (mesh[1] == adj[1]) {
      for (int i = order-2; i > 0; i--) {
        for (int j = order-i-1; j > 0; j--,n++) {
          int k = order-i-j;
          int p = order-k-1;
          for (int d = 0; d < dim; d++) {
            out[dim*(p*(p+1)/2-j)+d] = trigNodes[dim*n+d];
          }
        }
      }
    }
  }
}
