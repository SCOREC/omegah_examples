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
multi(Real* curr, Real* reff, Real* coord, LO para) {
  Real out = 1;
  for (int i = 0; i < para; i++) {
    if (curr[i] != reff[i])
      continue;
    out *= (coord[i]-reff[i])/(curr[i]-reff[i]);
  }
  return out;
}

OMEGA_H_INLINE
Real monoTerm2D(Real* curr, Real* coord, LO order) {
  Real out = 1;
  Real reff[3] = {1,0,0};
  out *= multi(curr,reff,coord,3);
  reff = {0,1,0};
  out *= multi(curr,reff,coord,3);
  reff = {0,0,1};
  out *= multi(curr,reff,coord,3);

  for (int i = order-1; i > 0; i--) {
    int j = order-i;
    reff = {i/order,j/order,0};
    out *= multi(curr,reff,coord,3);
    reff = {0,i/order,j/order};
    out *= multi(curr,reff,coord,3);
    reff = {j/order,0,i/order};
    out *= multi(curr,reff,coord,3);
  }
  
  for (int i = order-2; i > 0; i--) {
    for (int j = order-i-1; j > 0; j--) {
      int k = order-i-j;
      reff = {i/order,j/order,k/order};
      out *= multi(curr,reff,coord,3);
    }
  }
}
/*
OMEGA_H_INLINE
Real monoTerm3D(Real* curr, Real* coord, LO order) {
  Real out = 1;
  Real reff[4] = {1,0,0,0};
  out *= multi(curr,reff,coord,4);
  reff = {0,1,0,0};
  out *= multi(curr,reff,coord,4);
  reff = {0,0,1,0};
  out *= multi(curr,reff,coord,4);
  reff = {0,0,0,1};
  out *= multi(curr,reff,coord,4);
  
  for (int i = order-1; i > 0; i--) {
    int j = order-i;
    reff = {i/order,j/order,0,0};
    out *= multi(curr,reff,coord,4);
    reff = {i/order,0,j/order,0};
    out *= multi(curr,reff,coord,4);
    reff = {i/order,0,0,j/order};
    out *= multi(curr,reff,coord,4);
    reff = {0,i/order,j/order,0};
    out *= multi(curr,reff,coord,4);
    reff = {0,i/order,0,j/order}
    out *= multi(curr,reff,coord,4);
    reff = {0,0,i/order,j/order};
    out *= multi(curr,reff,coord,4);
  }
  
  for (int i = order-2; i > 0; i--) {
    for (int j = order-i-1; j > 0; j--) {
      int k = order-i-j;
      reff = {i/order,j/order,k/order,0};
      out *= multi(curr,reff,coord,4);
      reff = {i/order,j/order,0,k/order};
      out *= multi(curr,reff,coord,4);
      reff = {i/order,0,j/order,k/order};
      out *= multi(curr,reff,coord,4);
      reff = {0,i/order,j/order,k/order};
      out *= multi(curr,reff,coord,4);
    }
  }
  
  for (int i = order-3; i > 0; i--) {
    for (int j = order-i-2; j > 0; j--) {
      for (int k = order-i-j-1; k > 0; k--) {
        int l = order-i-j-k;
        reff = {i/order,j/order,k/order,l/order};
        out *= multi(curr,reff,coord,4);
      }
    }
  }
}
*/



OMEGA_H_INLINE
void shapeFunction2D(Real* coord, Real* out, Real* vertDOFs, Real* edgeDOFs, Real* trigDOFs, LO order) {
  for (int i = 0; i < 2; i++) {
    out = 0;
  }
 
  Real reff[3] = {1,0,0};
  Real l = monoTerm2D(reff,coord,order);
  out[0] += l*vertDOFs[2*0+0];
  out[1] += l*vertDOFs[2*0+1];
  reff = {0,1,0};
  l = monoTerm2D(reff,coord,order);
  out[0] += l*vertDOFs[2*1+0];
  out[1] += l*vertDOFs[2*1+1];
  reff = {0,0,1};
  l = monoTerm2D(reff,coord,order);
  out[0] += l*vertDOFs[2*2+0];
  out[1] += l*vertDOFs[2*2+1];

  int n = 0;
  int numNodes = order-1;
  for (int i = order-1; i > 0; i--,n++) {
    int j = order-i;
    reff = {i/order,j/order,0};
    l = monoTerm2D(reff,coord,order);
    out[0] += l*edgeDOFs[2*(numNodes*0+n)+0];
    out[1] += l*edgeDOFs[2*(numNodes*0+n)+1];
    reff = {0,i/order,j/order};
    l = monoTerm2D(reff,coord,order);
    out[0] += l*edgeDOFs[2*(numNodes*1+n)+0];
    out[1] += l*edgeDOFs[2*(numNodes*1+n)+1];
    reff = {j/order,0,i/order};
    l = monoTerm2D(reff,coord,order);
    out[0] += l*edgeDOFs[2*(numNodes*1+n)+0];
    out[1] += l*edgeDOFs[2*(numNodes*1+n)+1];
  }
  
  n = 0;
  for (int i = order-2; i > 0; i--) {
    for (int j = order-i-1; j > 0; j--) {
      int k = order-i-j;
      reff = {i/order,j/order,k/order};
      out *= multi(curr,reff,coord,3);
      l = monoTerm2D(reff,coord,order);
      out[0] += l*edgeDOFs[2*n+0];
      out[1] += l*edgeDOFs[2*n+1];
    }
  }
}

/*
OMEGA_H_INLINE
void shapFunc3D(Real* coord, Real* vertDOFs, Real* edgeDOFs, Real* trigDOFs, Real* tetsDOFs, LO order, Real* out) {  
  for (int i = 0; i < dim; i++) {

  }
  for (int i = 0; 
}*/
