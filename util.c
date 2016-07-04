/*
  Copyright (C) 2016 Edmond Shehadi
  This program is distributed under the terms of the GNU GPL v3 License.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "header.h"
#include "util.h"


void initialize(void)
{
  int i, j, k;
  M = (double *)calloc(Nx*Ny*Nz, sizeof(double));

  // fill average inner values
  for (k=1; k<Nz-1; k++){
    for (j=1; j<Ny-1; j++){
      for (i=1; i<Nx-1; i++){
        M[Nx*Ny*k+j*Nx+i] = 273.0; // Kelvin
      }
    }
  }

  // boundary: top surface
  k = Nz-1;
  for (j=1; j<Ny-1; j++){
    for (i=1; i<Nx-1; i++){
      M[Nx*Ny*k+j*Nx+i] = 280.0; // Kelvin
    }
  }
  // boundary: bottom surface
  k = 0;
  for (j=1; j<Ny-1; j++){
    for (i=1; i<Nx-1; i++){
      M[Nx*Ny*k+j*Nx+i] = 250.0; // Kelvin
    }
  }
  // boundary: west surface
  i = 0;
  for (k=1; k<Nz-1; k++){
    for (j=1; j<Ny-1; j++){
      M[Nx*Ny*k+j*Nx+i] = 250.0; // Kelvin
    }
  }
  // boundary: east surface
  i = Nx-1;
  for (k=1; k<Nz-1; k++){
    for (j=1; j<Ny-1; j++){
      M[Nx*Ny*k+j*Nx+i] = 290.0; // Kelvin
    }
  }
  // boundary: north surface
  j = 0;
  for (k=1; k<Nz-1; k++){
    for (i=1; i<Nx-1; i++){
      M[Nx*Ny*k+j*Nx+i] = 250.0;  // Kelvin
    }
  }
  // boundary: south surface
  j = Ny-1;
  for (k=1; k<Nz-1; k++){
    for (i=1; i<Nx-1; i++){
      M[Nx*Ny*k+j*Nx+i] = 290.0; // Kelvin
    }
  }

  // interpolate edges
  j = 0; k = 0;
  for (i=1; i<Nx-1; i++){
    M[Nx*Ny*k+j*Nx+i] = (M[Nx*Ny*(k+1)+j*Nx+i] + M[Nx*Ny*k+(j+1)*Nx+i] )/2.0;
  }
  j = Ny-1; k = 0;
  for (i=1; i<Nx-1; i++){
    M[Nx*Ny*k+j*Nx+i] = (M[Nx*Ny*(k+1)+j*Nx+i] + M[Nx*Ny*k+(j-1)*Nx+i] )/2.0;
  }
  j = 0; k = Nz-1;
  for (i=1; i<Nx-1; i++){
    M[Nx*Ny*k+j*Nx+i] = (M[Nx*Ny*(k-1)+j*Nx+i] + M[Nx*Ny*k+(j+1)*Nx+i] )/2.0;
  }
  j = Ny-1; k = Nz-1;
  for (i=1; i<Nx-1; i++){
    M[Nx*Ny*k+j*Nx+i] = (M[Nx*Ny*(k-1)+j*Nx+i] + M[Nx*Ny*k+(j-1)*Nx+i] )/2.0;
  }
  i = 0; k = 0;
  for (j=1; j<Ny-1; j++){
    M[Nx*Ny*k+j*Nx+i] = (M[Nx*Ny*(k+1)+j*Nx+i] + M[Nx*Ny*k+j*Nx+i+1] )/2.0;
  }
  i = Nx-1; k = 0;
  for (j=1; j<Ny-1; j++){
    M[Nx*Ny*k+j*Nx+i] = (M[Nx*Ny*(k+1)+j*Nx+i] + M[Nx*Ny*k+j*Nx+i-1] )/2.0;
  }
  i = 0; k = Nz-1;
  for (j=1; j<Ny-1; j++){
    M[Nx*Ny*k+j*Nx+i] = (M[Nx*Ny*(k-1)+j*Nx+i] + M[Nx*Ny*k+j*Nx+i+1] )/2.0;
  }
  i = Nx-1; k = Nz-1;
  for (j=1; j<Ny-1; j++){
    M[Nx*Ny*k+j*Nx+i] = (M[Nx*Ny*(k-1)+j*Nx+i] + M[Nx*Ny*k+j*Nx+i-1] )/2.0;
  }

  // interpolate corners + remaining 4 z-axis edges
  i = 0; j = 0;
  for (k=0; k<Nz; k++){
    M[Nx*Ny*k+j*Nx+i] = (M[Nx*Ny*k+j*Nx+i+1] + M[Nx*Ny*k+(j+1)*Nx+i])/2.0;
  }
  i = Nx-1; j = Ny-1;
  for (k=0; k<Nz; k++){
    M[Nx*Ny*k+j*Nx+i] = (M[Nx*Ny*k+j*Nx+i-1] + M[Nx*Ny*k+(j-1)*Nx+i])/2.0;
  }
  i = Nx-1; j = 0;
  for (k=0; k<Nz; k++){
    M[Nx*Ny*k+j*Nx+i] = (M[Nx*Ny*k+j*Nx+i-1] + M[Nx*Ny*k+(j+1)*Nx+i])/2.0;
  }
  i = 0; j = Ny-1;
  for (k=0; k<Nz; k++){
    M[Nx*Ny*k+j*Nx+i] = (M[Nx*Ny*k+j*Nx+i+1] + M[Nx*Ny*k+(j-1)*Nx+i])/2.0;
  }

  // check output
  char checkOUT = 0;
  if (checkOUT){
    printf("\n" );
    for (k=0; k<Nz; k++){
      for (j=0; j<Ny; j++){
        for (i=0; i<Nx; i++){
          // array math: array[M*L*k + L*j + i]
          printf("%3.0f ", M[Nx*Ny*k + j*Nx + i]);
        }
        printf("\n" );
      }
      printf("\n" );
    }
    printf("\n" );
  }

}


void boundary(void)
{
  /*
    stencil kernel:

                               (T)   (N)
                                |   /
                                | /
                     (W) ------ O ------ (E)
                              / |
                            /   |
                         (S)   (B)

    T: top     B: bottom
    N: north   S: south
    W: west    E: east

    orientation:
                        ^  (+k)
                        |
                        |
                        o ------> (+i)
                       /
                     /
                  (+j)

  */

  // this means there is insulation
  // in (inner surface only):
  //                         (jk)-plane @ i=0   (WEST)
  //                         (ik)-plane @ j=0   (NORTH)
  //                         (ij)-plane @ k=0   (BOTTOM)
  //
  insulated(0, 0); // plane _| i (0), j (1), k(2)
  insulated(0, 1); // plane _| i (0), j (1), k(2)
  insulated(0, 2); // plane _| i (0), j (1), k(2)

  convective(Nx-1, 0, h, Tsurr); // plane _| i (0), j (1), k(2)
  convective(Ny-1, 1, h, Tsurr); // plane _| i (0), j (1), k(2)
  convective(Nz-1, 2, h, Tsurr); // plane _| i (0), j (1), k(2)

  // @TODO:
  //       fixed boundary as a function
  // Dirichlet(Nx-1, 0, 200.0);
  // Dirichlet(Ny-1, 1, 250.0);

  allCorners(1, h, Tsurr,
             0, 0, 0,
             0, 0, 0,
             1, h, Tsurr,
             1, h, Tsurr,
             0, 0, 0 );


  // @TODO:
  //       corner function
  // @XXX: done !? (double-check)

  // @TODO:
  //       vertices function

  // @IDEA:
  //       radiation function

  // @IDEA:
  //       heat generation function

}


void allVertices(int isConv_S, double hs, double Ts, // convection south -side
                 int isConv_N, double hn, double Tn, // convection north -side
                 int isConv_W, double hw, double Tw, // convection west  -side
                 int isConv_E, double he, double Te, // convection east  -side
                 int isConv_T, double ht, double Tt, // convection top   -side
                 int isConv_B, double hb, double Tb) // convection bottom-side
{
/*
      stencil kernel:

                                 (T)   (N)
                                  |   /
                                  | /
                       (W) ------ O ------ (E)
                                / |
                              /   |
                           (S)   (B)

      T: top     B: bottom
      N: north   S: south
      W: west    E: east

      vertex diagram:

                               G---(12)-----H
                             / |          / |
                          (8) (9)     (7)  (11)
                        /      |     /      |    <== (BACK)
                       C--(4)--|---D        |
                       |       E---|-(10)---F
        (FRONT) ==>    |      /    |      /
                      (1)  (5)    (3)  (6)       <== (MIDDLE)
                       | /         |  /
                       A----(2)----B
    legend:
      corners:
        A: front, left-bottom   B: front, right-bottom
        C: front, left-top      D: front, right-top
        E: back, left-bottom    F: back, right-bottom
        G: back, left-top       H: back, right-top

      vertices:
           FRONT           MIDDLE            BACK
        -----------    ---------------    ------------
        (1): left      (5): lower-left    (09): left
        (2): bottom    (6): lower-right   (10): bottom
        (3): right     (7): upper-right   (11): right
        (4): top       (8): upper-left    (12): top
*/

  int i, j, k, kc, kw, ke, kn, ks, kt, kb;
  double coef1, coef2, coef3, term, denom;

  assert( isConv_S == 0 || isConv_S == 1 );
  assert( isConv_N == 0 || isConv_N == 1 );
  assert( isConv_W == 0 || isConv_W == 1 );
  assert( isConv_E == 0 || isConv_E == 1 );
  assert( isConv_T == 0 || isConv_T == 1 );
  assert( isConv_B == 0 || isConv_B == 1 );

  // vertex 1
  i = 0; j = Ny-1;

  coef1 = kd*dy*dz/(2*dx);
  coef2 = kd*dx*dz/(2*dy);
  coef3 = kd*dx*dy/(4*dz);

  term = 0.0;
  term += (isConv_S)  ? hs*dx*dz*Ts/(2*Ch*rho)  : 0.0;
  term += (isConv_W)  ? hw*dy*dz*Tw/(2*Ch*rho)  : 0.0;

  denom = coef1 + coef2 + 2*coef3;
  denom += (isConv_S) ? hs*dx*dz/(2*Ch*rho) : 0.0;
  denom += (isConv_W) ? hw*dy*dz/(2*Ch*rho) : 0.0;

  for (k=1; k<Nz-1; k++){
    kc = Nx*Ny*k+j*Nx+i;
    kb = kc-Nx*Ny;
    kn = kc-Nx;
    ke = kc+1;

    M[kc] = ( coef1*M[ke] + \
              coef2*M[kn] + \
              coef3*M[kb] + \
              coef3*M[kt] + \
              term )/denom;
  }
  // vertex 2
  j = Ny-1; k = 0;

  coef1 = kd*dy*dz/(4*dx);
  coef2 = kd*dz*dx/(2*dy);
  coef3 = kd*dy*dx/(2*dz);

  term = 0.0;
  term += (isConv_S)  ? hs*dz*dx*Ts/(2*Ch*rho)  : 0.0;
  term += (isConv_B)  ? hb*dy*dx*Tb/(2*Ch*rho)  : 0.0;

  denom = 2*coef1 + coef2 + coef3;
  denom += (isConv_S) ? hs*dz*dx/(2*Ch*rho) : 0.0;
  denom += (isConv_B) ? hb*dy*dx/(2*Ch*rho) : 0.0;

  for (i=1; i<Nx-1; i++){
    kc = Nx*Ny*k+j*Nx+i;
    kt = kc+Nx*Ny;
    kn = kc-Nx;
    ke = kc+1;
    kw = kc-1;

    M[kc] = ( coef1*M[kw] + \
              coef1*M[ke] + \
              coef2*M[kn] + \
              coef3*M[kt] + \
              term )/denom;
  }
  // vertex 3
  i = Nx-1; j = Ny-1;

  coef1 = kd*dy*dz/(2*dx);
  coef2 = kd*dx*dz/(2*dy);
  coef3 = kd*dx*dy/(4*dz);

  term = 0.0;
  term += (isConv_S)  ? hs*dx*dz*Ts/(2*Ch*rho)  : 0.0;
  term += (isConv_E)  ? he*dy*dz*Te/(2*Ch*rho)  : 0.0;

  denom = coef1 + coef2 + 2*coef3;
  denom += (isConv_S) ? hs*dx*dz/(2*Ch*rho) : 0.0;
  denom += (isConv_E) ? he*dy*dz/(2*Ch*rho) : 0.0;

  for (k=1; k<Nz-1; k++){
    kc = Nx*Ny*k+j*Nx+i;
    kt = kc+Nx*Ny;
    kb = kc-Nx*Ny;
    kn = kc-Nx;
    kw = kc-1;

    M[kc] = ( coef1*M[kw] + \
              coef2*M[kn] + \
              coef3*M[kt] + \
              coef3*M[kb] + \
              term )/denom;
  }
  // vertex 4
  j = Ny-1; k = Nz-1;

  coef1 = kd*dy*dz/(4*dx);
  coef2 = kd*dz*dx/(2*dy);
  coef3 = kd*dy*dx/(2*dz);

  term = 0.0;
  term += (isConv_S)  ? hs*dz*dx*Ts/(2*Ch*rho)  : 0.0;
  term += (isConv_T)  ? ht*dy*dx*Tt/(2*Ch*rho)  : 0.0;

  denom = 2*coef1 + coef2 + coef3;
  denom += (isConv_S) ? hs*dz*dx/(2*Ch*rho) : 0.0;
  denom += (isConv_T) ? ht*dy*dx/(2*Ch*rho) : 0.0;

  for (i=1; i<Nx-1; i++){
    kc = Nx*Ny*k+j*Nx+i;
    kb = kc-Nx*Ny;
    kn = kc-Nx;
    kw = kc-1;
    ke = kc+1;

    M[kc] = ( coef1*M[kw] + \
              coef1*M[ke] + \
              coef2*M[kn] + \
              coef3*M[kb] + \
              term )/denom;
  }
  // vertex 5
  i = 0; k = 0;

  coef1 = kd*dz*dy/(2*dx);
  coef2 = kd*dx*dz/(4*dy);
  coef3 = kd*dx*dy/(2*dz);

  term = 0.0;
  term += (isConv_W)  ? hw*dz*dy*Tw/(2*Ch*rho)  : 0.0;
  term += (isConv_B)  ? hb*dx*dy*Tb/(2*Ch*rho)  : 0.0;

  denom = coef1 + 2*coef2 + coef3;
  denom += (isConv_W) ? hw*dz*dy/(2*Ch*rho) : 0.0;
  denom += (isConv_B) ? hb*dx*dy/(2*Ch*rho) : 0.0;

  for (j=1; j<Ny-1; j++){
    kc = Nx*Ny*k+j*Nx+i;
    kt = kc+Nx*Ny;
    ks = kc+Nx;
    kn = kc-Nx;
    ke = kc+1;

    M[kc] = ( coef1*M[ke] + \
              coef2*M[kn] + \
              coef2*M[ks] + \
              coef3*M[kt] + \
              term )/denom;
  }
  // vertex 6
  i = Nx-1; k = 0;

  coef1 = kd*dz*dy/(2*dx);
  coef2 = kd*dx*dz/(4*dy);
  coef3 = kd*dx*dy/(2*dz);

  term = 0.0;
  term += (isConv_E)  ? he*dz*dy*Te/(2*Ch*rho)  : 0.0;
  term += (isConv_B)  ? hb*dx*dy*Tb/(2*Ch*rho)  : 0.0;

  denom = coef1 + 2*coef2 + coef3;
  denom += (isConv_E) ? ht*dz*dy/(2*Ch*rho) : 0.0;
  denom += (isConv_B) ? hb*dx*dy/(2*Ch*rho) : 0.0;

  for (j=1; j<Ny-1; j++){
    kc = Nx*Ny*k+j*Nx+i;
    kt = kc+Nx*Ny;
    ks = kc+Nx;
    kn = kc-Nx;
    kw = kc-1;

    M[kc] = ( coef1*M[kw] + \
              coef2*M[kn] + \
              coef2*M[ks] + \
              coef3*M[kt] + \
              term )/denom;
  }
  // vertex 7
  i = Nx-1; k = Nz-1;

  coef1 = kd*dz*dy/(2*dx);
  coef2 = kd*dx*dz/(4*dy);
  coef3 = kd*dx*dy/(2*dz);

  term = 0.0;
  term += (isConv_T)  ? ht*dx*dy*Tt/(2*Ch*rho)  : 0.0;
  term += (isConv_N)  ? hn*dz*dy*Tn/(2*Ch*rho)  : 0.0;

  denom = coef1 + 2*coef2 + coef3;
  denom += (isConv_T) ? ht*dx*dy/(2*Ch*rho) : 0.0;
  denom += (isConv_N) ? hn*dz*dy/(2*Ch*rho) : 0.0;

  for (j=1; j<Ny-1; j++){
    kc = Nx*Ny*k+j*Nx+i;
    kb = kc-Nx*Ny;
    ks = kc+Nx;
    kn = kc-Nx;
    kw = kc-1;

    M[kc] = ( coef1*M[kw] + \
              coef2*M[kn] + \
              coef2*M[ks] + \
              coef3*M[kb] + \
              term )/denom;
  }
  // vertex 8

  // vertex 9

  // vertex 10

  // vertex 11

  // vertex 12

}


void allCorners(int isConv_S, double hs, double Ts, // convection south -side
                int isConv_N, double hn, double Tn, // convection north -side
                int isConv_W, double hw, double Tw, // convection west  -side
                int isConv_E, double he, double Te, // convection east  -side
                int isConv_T, double ht, double Tt, // convection top   -side
                int isConv_B, double hb, double Tb) // convection bottom-side
{
/*
    stencil kernel:

                               (T)   (N)
                                |   /
                                | /
                     (W) ------ O ------ (E)
                              / |
                            /   |
                         (S)   (B)

    T: top     B: bottom
    N: north   S: south
    W: west    E: east

    corner diagram:

                          G--------H
                        / |      / |
                      /   |    /   |    <== (BACK)
                     C----|---D    |
                     |    E---|----F
      (FRONT) ==>    |   /    |   /
                     | /      | /
                     A--------B
  legend:

    A: front, left-bottom   B: front, right-bottom
    C: front, left-top      D: front, right-top
    E: back, left-bottom    F: back, right-bottom
    G: back, left-top       H: back, right-top

*/

  //  @TODO:
  //        add radiation ?
  int i, j, k, kc, kw, ke, kn, ks, kt, kb;
  double Coef_ij, Coef_ik, Coef_jk, term;

  double denom, cond0 = kd*( dy*dz/dx + dx*dz/dy + dx*dy/dz )/4.0;

  assert( isConv_S == 0 || isConv_S == 1 );
  assert( isConv_N == 0 || isConv_N == 1 );
  assert( isConv_W == 0 || isConv_W == 1 );
  assert( isConv_E == 0 || isConv_E == 1 );
  assert( isConv_T == 0 || isConv_T == 1 );
  assert( isConv_B == 0 || isConv_B == 1 );

  Coef_ij = (dx/2)*(dy/2)/(Ch*rho); // convection constant for ij-plane
  Coef_ik = (dx/2)*(dz/2)/(Ch/rho); // convection constant for ik-plane
  Coef_jk = (dy/2)*(dz/2)/(Ch*rho); // convection constant for jk-plane

  // corner A
  i = 0; j = Ny-1; k = 0;

  kc = Nx*Ny*k+j*Nx+i;
  kt = kc+Nx*Ny;
  kn = kc-Nx;
  ke = kc+1;

  term = 0.0;
  term += (isConv_W)  ? hw*Coef_jk*Tw  : 0.0;
  term += (isConv_B)  ? hb*Coef_ij*Tb  : 0.0;
  term += (isConv_S)  ? hs*Coef_ik*Ts  : 0.0;

  denom = 0.0;
  denom += (isConv_W) ? hw*Coef_jk  : 0.0;
  denom += (isConv_B) ? hb*Coef_ij  : 0.0;
  denom += (isConv_S) ? hs*Coef_ik  : 0.0;
  denom += cond0;

  M[kc] = ( (kd*dy*dz/(4*dx))*M[ke] + \
            (kd*dx*dz/(4*dy))*M[kn] + \
            (kd*dx*dy/(4*dz))*M[kt] + \
            term )/denom;

  // corner B
  i = Nx-1; j = Ny-1; k = 0;

  kc = Nx*Ny*k+j*Nx+i;
  kt = kc+Nx*Ny;
  kn = kc-Nx;
  kw = kc-1;

  term = 0.0;
  term  = (isConv_S)  ? hs*Coef_ik*Ts  : 0.0;
  term += (isConv_E)  ? he*Coef_jk*Te  : 0.0;
  term += (isConv_B)  ? hb*Coef_ij*Tb  : 0.0;

  denom = 0.0;
  denom += (isConv_S) ? hs*Coef_ik  : 0.0;
  denom += (isConv_E) ? he*Coef_jk  : 0.0;
  denom += (isConv_B) ? hb*Coef_ij  : 0.0;
  denom += cond0;

  M[kc] = ( (kd*dy*dz/(4*dx))*M[kw] + \
            (kd*dx*dz/(4*dy))*M[kn] + \
            (kd*dx*dy/(4*dz))*M[kt] + \
            term )/denom;

  // corner C
  i = 0;  j = Ny-1; k = Nz-1;

  kc = Nx*Ny*k+j*Nx+i;
  kb = kc-Nx*Ny;
  kn = kc-Nx;
  ke = kc+1;

  term = 0.0;
  term += (isConv_S)  ? hs*Coef_ik*Ts : 0.0;
  term += (isConv_W)  ? hw*Coef_jk*Tw : 0.0;
  term += (isConv_T)  ? ht*Coef_ij*Tt : 0.0;

  denom = 0.0;
  denom += (isConv_S) ? hs*Coef_ik  : 0.0;
  denom += (isConv_W) ? hw*Coef_jk  : 0.0;
  denom += (isConv_T) ? ht*Coef_ij  : 0.0;
  denom += cond0;

  M[kc] = ( (kd*dy*dz/(4*dx))*M[ke] + \
            (kd*dx*dz/(4*dy))*M[kn] + \
            (kd*dx*dy/(4*dz))*M[kb] + \
            term )/denom;

  // corner D
  i = Nx-1; j = Ny-1; k = Nz-1;

  kc = Nx*Ny*k+j*Nx+i;
  kb = kc-Nx*Ny;
  kn = kc-Nx;
  kw = kc-1;

  term = 0.0;
  term += (isConv_S)  ? hs*Coef_ik*Ts : 0.0;
  term += (isConv_E)  ? he*Coef_jk*Te : 0.0;
  term += (isConv_T)  ? ht*Coef_ij*Tt : 0.0;

  denom = 0.0;
  denom += (isConv_S) ? hs*Coef_ik  : 0.0;
  denom += (isConv_E) ? he*Coef_jk  : 0.0;
  denom += (isConv_T) ? ht*Coef_ij  : 0.0;
  denom += cond0;

  M[kc] = ( (kd*dy*dz/(4*dx))*M[kw] + \
            (kd*dx*dz/(4*dy))*M[kn] + \
            (kd*dx*dy/(4*dz))*M[kb] + \
            term )/denom;

  // corner E
  i = 0; j = 0; k = 0;

  kc = Nx*Ny*k+j*Nx+i;
  kt = kc+Nx*Ny;
  ks = kc+Nx;
  ke = kc+1;

  term = 0.0;
  term += (isConv_W)  ? hw*Coef_jk*Tw : 0.0;
  term += (isConv_N)  ? hn*Coef_ik*Tn : 0.0;
  term += (isConv_B)  ? hb*Coef_ij*Tb : 0.0;

  denom = 0.0;
  denom += (isConv_W) ? hw*Coef_jk  : 0.0;
  denom += (isConv_N) ? hn*Coef_ik  : 0.0;
  denom += (isConv_B) ? hb*Coef_ij  : 0.0;
  denom += cond0;

  M[kc] = ( (kd*dy*dz/(4*dx))*M[ke] + \
            (kd*dx*dz/(4*dy))*M[ks] + \
            (kd*dx*dy/(4*dz))*M[kt] + \
            term )/denom;

  // corner F
  i = Nx-1; j = 0; k = 0;

  kc = Nx*Ny*k+j*Nx+i;
  kt = kc+Nx*Ny;
  ks = kc+Nx;
  kw = kc-1;

  term = 0.0;
  term += (isConv_E)  ? he*Coef_jk*Te : 0.0;
  term += (isConv_N)  ? hn*Coef_ik*Tn : 0.0;
  term += (isConv_B)  ? hb*Coef_ij*Tb : 0.0;

  denom = 0.0;
  denom += (isConv_E) ? he*Coef_jk  : 0.0;
  denom += (isConv_N) ? hn*Coef_ik  : 0.0;
  denom += (isConv_B) ? hb*Coef_ij  : 0.0;
  denom += cond0;

  M[kc] = ( (kd*dy*dz/(4*dx))*M[kw] + \
            (kd*dx*dz/(4*dy))*M[ks] + \
            (kd*dx*dy/(4*dz))*M[kt] + \
            term )/denom;

  // corner G
  i = 0; j = 0; k = Nz-1;

  kc = Nx*Ny*k+j*Nx+i;
  kb = kc-Nx*Ny;
  ks = kc+Nx;
  ke = kc+1;

  term = 0.0;
  term += (isConv_W)  ? hw*Coef_jk*Tw : 0.0;
  term += (isConv_N)  ? hn*Coef_ik*Tn : 0.0;
  term += (isConv_T)  ? ht*Coef_ij*Tt : 0.0;

  denom = 0.0;
  denom += (isConv_W) ? hw*Coef_jk  : 0.0;
  denom += (isConv_N) ? hn*Coef_ik  : 0.0;
  denom += (isConv_T) ? ht*Coef_ij  : 0.0;
  denom += cond0;

  M[kc] = ( (kd*dy*dz/(4/dx))*M[ke] + \
            (kd*dx*dz/(4*dy))*M[ks] + \
            (kd*dx*dy/(4*dz))*M[kb] + \
            term )/denom;

  // corner H
  i = Nx-1; j = 0; k = Nz-1;

  kc = Nx*Ny*k+j*Nx+i;
  kb = kc-Nx*Ny;
  ks = kc+Nx;
  kw = kc-1;

  term = 0.0;
  term += (isConv_E)  ? he*Coef_jk*Te : 0.0;
  term += (isConv_N)  ? hn*Coef_ik*Tn : 0.0;
  term += (isConv_T)  ? ht*Coef_ij*Tt : 0.0;

  denom = 0.0;
  denom += (isConv_E) ? he*Coef_jk  : 0.0;
  denom += (isConv_N) ? hn*Coef_ik  : 0.0;
  denom += (isConv_T) ? ht*Coef_ij  : 0.0;
  denom += cond0;

  M[kc] = ( (kd*dy*dz/(4*dx))*M[kw] + \
            (kd*dx*dz/(4*dy))*M[ks] + \
            (kd*dx*dy/(4*dz))*M[kb] + \
            term )/denom;


}


void Dirichlet(int fix, int plane, double value)
{
  int i, j, k, kc;

#define X 0 //  plane perpendicular to i-axis
#define Y 1 //  plane perpendicular to j-axis
#define Z 2 //  plane perpendicular to k-axis

  switch (plane) {

    case X:
      i = fix;  // at which level to fix
      assert( i==0 || i==(Nx-1) );

      for (k=1; k<Nz-1; k++){
        for (j=1; j<Ny-1; j++){
          kc = Nx*Ny*k+j*Nx+i;
          M[kc] = value;
        }
      }
      break;

    case Y:
      j = fix; // at which level to fix
      assert( j==0 || j==(Ny-1) );

      for (k=1; k<Nz-1; k++){
        for (i=1; i<Nx-1; i++){
          kc = Nx*Ny*k+j*Nx+i;
          M[kc] = value;
        }
      }
      break;

    case Z:
      k = fix; // at which level to fix
      assert( k==0 || k==(Nz-1) );

      for (j=1; j<Ny-1; j++){
        for (i=1; i<Nx-1; i++){
          kc = Nx*Ny*k+j*Nx+i;
          M[kc] = value;
        }
      }
      break;

    default:
      printf("\nWrong input in 'Dirichlet function' \t plane: %d\n\n", plane);
      exit(1);
  }

}


void convective(int fix, int plane,
                double h, double Tsurr)
{
  /*
    This only takes care of convection in the inner surface of a given plane
    i.e. excluding corners and vertices.
  */
  int i, j, k, kw, ke, kt, kb, ks, kn, kc;
  double Ccp, Cxp, Cyp, Czp, Csp;

#define X 0 // plane perpendicular to i-axis
#define Y 1 // plane perpendicular to j-axis
#define Z 2 // plane perpendicular to k-axis

  switch (plane) {

    case X:
      i = fix; // at which level to fix
      assert( i==0 || i==(Nx-1) );

      Ccp = h*dy*dz/(Ch*rho) + kd*(dx*dz/dy + dy*dz/dx + dx*dy/dz);
      Cxp = kd*dy*dz/dx;
      Cyp = kd*dx*dz/(2*dy);
      Czp = kd*dx*dy/(2*dz);
      Csp = h*dy*dz/(Ch*rho);

      if (i==Nx-1){
        for (k=1; k<Nz-1; k++){
          for (j=1; j<Ny-1; j++){
            kc = Nx*Ny*k+j*Nx+i;
            kt = kc+Nx*Ny; kb = kc-Nx*Ny;
            ks = kc+Nx;    kn = kc-Nx;
            kw = kc-1;
            M[kc] = ( Cxp*M[kw] + Cyp*(M[kn]+M[ks]) + Czp*(M[kb]+M[kt]) + Csp*Tsurr )/Ccp;
          }
        }
      }
      else {
        for (k=1; k<Nz-1; k++){
          for (j=1; j<Ny-1; j++){
            kc = Nx*Ny*k+j*Nx+i;
            kt = kc+Nx*Ny; kb = kc-Nx*Ny;
            ks = kc+Nx;    kn = kc-Nx;
            ke = kc+1;
            M[kc] = ( Cxp*M[ke] + Cyp*(M[kn]+M[ks]) + Czp*(M[kb]+M[kt]) + Csp*Tsurr )/Ccp;
          }
        }
      }
      break;

    case Y:
      j = fix; // at which level to fix
      assert( j==0 || j==(Ny-1) );

      Ccp = h*dx*dz/(Ch*rho) + kd*(dx*dz/dy + dy*dz/dx + dx*dy/dz);
      Cxp = kd*dy*dz/(2*dx);
      Cyp = kd*dx*dz/dy;
      Czp = kd*dx*dy/(2*dz);
      Csp = h*dx*dz/(Ch*rho);

      if (j==Ny-1){
        for (k=1; k<Nz-1; k++){
          for (i=1; i<Nx-1; i++){
            kc = Nx*Ny*k+j*Nx+i;
            kt = kc+Nx*Ny; kb = kc-Nx*Ny;
            ke = kc+1;     kw = kc-1;
            kn = kc-Nx;
            M[kc] = ( Cxp*(M[ke]+M[kw]) + Cyp*M[kn] + Czp*(M[kb]+M[kt]) + Csp*Tsurr )/Ccp;
          }
        }
      }
      else {
        for (k=1; k<Nz-1; k++){
          for (i=1; i<Nx-1; i++){
            kc = Nx*Ny*k+j*Nx+i;
            kt = kc+Nx*Ny; kb = kc-Nx*Ny;
            ke = kc+1;     kw = kc-1;
            ks = kc-Nx;
            M[kc] = ( Cxp*(M[ke]+M[kw]) + Cyp*M[ks] + Czp*(M[kb]+M[kt]) + Csp*Tsurr )/Ccp;
          }
        }
      }
      break;

    case Z:
      k = fix; // at which level to fix
      assert( k==0 || k==(Nz-1) );

      Ccp = h*dx*dy/(Ch*rho) + kd*(dy*dz/dx + dx*dz/dy + dx*dy/dz);
      Cxp = kd*dy*dz/(2*dx);
      Cyp = kd*dx*dz/(2*dy);
      Czp = kd*dx*dy/dz;
      Csp = h*dx*dy/(Ch*rho);

      if (k==Nz-1){
        for (j=1; j<Ny-1; j++){
          for (i=1; i<Nx-1; i++){
            kc = Nx*Ny*k+j*Nx+i;
            ks = kc+Nx;    kn = kc-Nx;
            ke = kc+1;     kw = kc-1;
            kb = kc-Nx*Ny;
            M[kc] = ( Cxp*(M[kw]+M[ke]) + Cyp*(M[kn]+M[ks]) + Czp*M[kb] + Csp*Tsurr )/Ccp;
          }
        }
      }
      else {
        for (j=1; j<Ny-1; j++){
          for (i=1; i<Nx-1; i++){
            kc = Nx*Ny*k+j*Nx+i;
            ks = kc+Nx;    kn = kc-Nx;
            ke = kc+1;     kw = kc-1;
            kt = kc+Nx*Ny;
            M[kc] = ( Cxp*(M[kw]+M[ke]) + Cyp*(M[kn]+M[ks]) + Czp*M[kt] + Csp*Tsurr )/Ccp;
          }
        }
      }
      break;

    default:
      printf("\nWrong input in 'convective function' \t plane: %d\n\n", plane);
      exit(1);
  }

}


void insulated(int fix, int plane)
{
  /*
    This only takes care of insulation in the inner surface of a given plane
    i.e. excluding corners and vertices.
  */
  int i, j, k, step;

#define X 0  // plane perpendicular to i-axis
#define Y 1  // plane perpendicular to j-axis
#define Z 2  // plane perpendicular to k-axis

  switch(plane){

    case X:
      i = fix; // at which level to fix
      assert( i==0 || i==(Nx-1) );
      step = (i == 0)   ? 1 : -1;
      for (k=1; k<Nz-1; k++){
        for (j=1; j<Ny-1; j++){
          M[Nx*Ny*k+j*Nx+i] = M[Nx*Ny*k+j*Nx+i+step];
        }
      }
      break;
    case Y:
      j = fix; // at which level to fix
      assert( j==0 || j==(Ny-1) );
      step = (j == 0)   ? 1 : -1;
      for (k=1; k<Nz-1; k++){
        for (i=1; i<Nx-1; i++){
          M[Nx*Ny*k+j*Nx+i] = M[Nx*Ny*k+(j+step)*Nx+i];
        }
      }
      break;
    case Z:
      k = fix; // at which level to fix
      assert( k==0 || k==(Nz-1) );
      step = (k == 0)   ? 1 : -1;
      for (j=1; j<Ny-1; j++){
        for (i=1; i<Nx-1; i++){
          M[Nx*Ny*k+j*Nx+i] = M[Nx*Ny*(k+step)+j*Nx+i];
        }
      }
      break;

    default:
      printf("\nWrong input in 'insulation function' \t plane: %d\n\n", plane);
      exit(1);
  }

}
