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
  int i, j, k; double temp = 273.15+25;
  M = (double *)calloc(Nx*Ny*Nz, sizeof(double));

  // fill average inner values
  // if (isSourceTerm){ // source term
#if isSourceTerm
  // heat source exists
  #if isRegular
  // box-shaped geometry
    for (k=1; k<Nz-1; k++){
      for (j=1; j<Ny-1; j++){
        for (i=1; i<Nx-1; i++){
          M[Nx*Ny*k+j*Nx+i] = temp + dt*sourceGen(i, j, k)/(rho*Ch); // Kelvin
        }
      }
    }

  #else
  // sphere-shaped geometry
  // box-shaped geometry
    for (k=1; k<Nz-1; k++){
      for (j=1; j<Ny-1; j++){
        for (i=1; i<Nx-1; i++){
          M[Nx*Ny*k+j*Nx+i] = temp + dt*sourceGenIrreg(i, j, k)/(rho*Ch); // Kelvin
        }
      }
    }
  #endif

#else
  // no heat source(s)
  for (k=1; k<Nz-1; k++){
    for (j=1; j<Ny-1; j++){
      for (i=1; i<Nx-1; i++){
        M[Nx*Ny*k+j*Nx+i] = temp; // Kelvin
      }
    }
  }

#endif

  // boundary: top surface
  k = Nz-1;
  for (j=1; j<Ny-1; j++){
    for (i=1; i<Nx-1; i++){
      M[Nx*Ny*k+j*Nx+i] = 273.15+25; // Kelvin
    }
  }
  // boundary: bottom surface
  k = 0;
  for (j=1; j<Ny-1; j++){
    for (i=1; i<Nx-1; i++){
      M[Nx*Ny*k+j*Nx+i] = 273.15+25; // Kelvin
    }
  }
  // boundary: west surface
  i = 0;
  for (k=1; k<Nz-1; k++){
    for (j=1; j<Ny-1; j++){
      M[Nx*Ny*k+j*Nx+i] = 273.15+25; // Kelvin
    }
  }
  // boundary: east surface
  i = Nx-1;
  for (k=1; k<Nz-1; k++){
    for (j=1; j<Ny-1; j++){
      M[Nx*Ny*k+j*Nx+i] = 273.15+25; // Kelvin
    }
  }
  // boundary: north surface
  j = 0;
  for (k=1; k<Nz-1; k++){
    for (i=1; i<Nx-1; i++){
      M[Nx*Ny*k+j*Nx+i] = 273.15+25;  // Kelvin
    }
  }
  // boundary: south surface
  j = Ny-1;
  for (k=1; k<Nz-1; k++){
    for (i=1; i<Nx-1; i++){
      M[Nx*Ny*k+j*Nx+i] = 273.15+25; // Kelvin
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

#if isSourceTerm
    heatGeneration();
#endif


#if isRegular
  assert(length == width && length == height);
#endif
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


void boundary(int condition[6], double h[6],
              double Tsurr[6], double temp[6])
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

#define convection 1  // convection+conduction boundary
#define insulation 2  // insulation boundary
#define dirichCond 3  // fixed temperature boundary

  // 0 : T
  // 1 : B
  // 2 : E
  // 3 : W
  // 4 : S
  // 5 : N
  int isConv[6] = {0}; // NO convection -- DEFAULT

  for (int i=0; i<6; i++){
    assert(condition[i] == 1 || condition[i] == 2 || condition[i] == 3 );

    if (condition[i] == convection){
      isConv[i] = 1;
      surface(i, h[i], Tsurr[i]);
    }
    else if (condition[i] == insulation){
      insulated(i);
    }
    else if (condition[i] == dirichCond){
      Dirichlet(i, temp[i]);
    }

  }

  // @FIXME : debug corner and vertex functions ...
  //
  // @NOTE: prob seems in 'allVertices', model equations seems wrong -- redo !!
  //

  allVertices(isConv[4], h[4], Tsurr[4],   // convection south -side
              isConv[5], h[5], Tsurr[5],   // convection north -side
              isConv[3], h[3], Tsurr[3],   // convection west  -side
              isConv[2], h[2], Tsurr[2],   // convection east  -side
              isConv[0], h[0], Tsurr[0],   // convection top   -side
              isConv[1], h[1], Tsurr[1] ); // convection bottom-side

  allCorners(isConv[4], h[4], Tsurr[4],    // convection south -side
             isConv[5], h[5], Tsurr[5],    // convection north -side
             isConv[3], h[3], Tsurr[3],    // convection west  -side
             isConv[2], h[2], Tsurr[2],    // convection east  -side
             isConv[0], h[0], Tsurr[0],    // convection top   -side
             isConv[1], h[1], Tsurr[1] );  // convection bottom-side
// ---------------------------------------------------------------------
  // @IDEA:
  //       radiation function


}


void heatGeneration(void)
{
/*
  Dependant upon:

    geometry : <#>
                0 : cube/box
                1 : sphere/ellipsoidal

    heatGen  : <#>
                # : source's heat generation value -- constant

    posX     : <#>
                0.0-1.0 : relative position from cube origin -- ith sense

    posY     : <#>
                0.0-1.0 : relative position from cube origin -- jth sense

    posZ     : <#>
                0.0-1.0 : relative position from cube origin -- kth sense

    length   : <#>
                # : source's length in the ith sense from it's center

    width    : <#>
                # : source's length in the jth sense from it's center

    height   : <#>
                # : source's length in the kth sense from it's center


    refer to figure below for clarification,



                ^  (+k)
                |
                |
                |
                |                      (height)      (width)
                |
                |                             |   /
                |                             | /
                |            (length)  ------ O ------  (length)
                |                           / |
                |                         /   |
                |
                |                (width)      (height)
                |
                |
                |
                o -----------------------------------------> (+i)
               /
              /
             /
            /
           /
          /
       (+j)

  orientation and source (w.r.t. overall cube)

*/
  int Cx, Cy, Cz;

  assert(posX > 0.0 && posX < 1.0);
  assert(posY > 0.0 && posY < 1.0);
  assert(posZ > 0.0 && posZ < 1.0);

  Cx = round(Nx*posX);  // center of source geometry in ith sense
  Cy = round(Ny*posY);  // center of source geometry in jth sense
  Cz = round(Nz*posZ);  // center of source geometry in kth sense

  x0 = Cx-round(length*Nx); x1 = Cx+round(length*Nx);
  y0 = Cy- round(width*Ny); y1 = Cy+ round(width*Ny);
  z0 = Cz-round(height*Nz); z1 = Cz+round(height*Nz);

  assert( x0 > 0 && x1 < Nx );
  assert( y0 > 0 && y1 < Ny );
  assert( z0 > 0 && z1 < Nz );

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
        (FRONT) ==>    |     /     |      /
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
  double Cx, Cy, Cz, term, denom, vol;

  assert( isConv_S == 0 || isConv_S == 1 );
  assert( isConv_N == 0 || isConv_N == 1 );
  assert( isConv_W == 0 || isConv_W == 1 );
  assert( isConv_E == 0 || isConv_E == 1 );
  assert( isConv_T == 0 || isConv_T == 1 );
  assert( isConv_B == 0 || isConv_B == 1 );

  // infinitesmal volume
  vol = dx*dy*dz/2.0;

  // vertex 1
  i = 0; j = Ny-1;

  Cx = dt*Kd/(vol*dx);
  Cy = dt*Kd/(vol*dy);
  Cz = dt*Kd/(vol*dz);

  term = 0.0;
  term += (isConv_S)  ? hs*dx*dz*Ts/(2*Ch*rho*vol)  : 0.0;
  term += (isConv_W)  ? hw*dy*dz*Tw/(2*Ch*rho*vol)  : 0.0;
  term *= dt;

  denom = coef1 + coef2 + coef3;
  denom += (isConv_S) ? dt*hs*dx*dz/(2*Ch*rho) : 0.0;
  denom += (isConv_W) ? dt*hw*dy*dz/(2*Ch*rho) : 0.0;

  for (k=1; k<Nz-1; k++){
    kc = Nx*Ny*k+j*Nx+i;
    kt = kc+Nx*Ny;
    kb = kc-Nx*Ny;
    kn = kc-Nx;
    ke = kc+1;

    M[kc] = w*( coef1*M[ke] + \
                coef2*M[kn] + \
                coef3*M[kb] + \
                coef3*M[kt] + \
                term + (1.0 - denom)*M[kc] ) + (1.0 - w)*M[kc];
  }
  // vertex 2
  j = Ny-1; k = 0;

  coef1 = dt*kd/(dx*dx);
  coef2 = dt*kd/(dy*dy);
  coef3 = dt*kd/(dz*dz);

  term = 0.0;
  term += (isConv_S)  ? hs*dz*dx*Ts/(2*Ch*rho)  : 0.0;
  term += (isConv_B)  ? hb*dy*dx*Tb/(2*Ch*rho)  : 0.0;
  term *= dt;

  denom = coef1 + coef2 + coef3;
  denom += (isConv_S) ? dt*hs*dz*dx/(2*Ch*rho) : 0.0;
  denom += (isConv_B) ? dt*hb*dy*dx/(2*Ch*rho) : 0.0;

  for (i=1; i<Nx-1; i++){
    kc = Nx*Ny*k+j*Nx+i;
    kt = kc+Nx*Ny;
    kn = kc-Nx;
    ke = kc+1;
    kw = kc-1;

    M[kc] = w*( coef1*M[kw] + \
                coef1*M[ke] + \
                coef2*M[kn] + \
                coef3*M[kt] + \
                term + (1.0 - denom)*M[kc] ) + (1.0 - w)*M[kc];
  }
  // vertex 3
  i = Nx-1; j = Ny-1;

  coef1 = dt*kd/(dx*dx);
  coef2 = dt*kd/(dy*dy);
  coef3 = dt*kd/(dz*dz);

  term = 0.0;
  term += (isConv_S)  ? hs*dx*dz*Ts/(2*Ch*rho)  : 0.0;
  term += (isConv_E)  ? he*dy*dz*Te/(2*Ch*rho)  : 0.0;
  term *= dt;

  denom = coef1 + coef2 + coef3;
  denom += (isConv_S) ? dt*hs*dx*dz/(2*Ch*rho) : 0.0;
  denom += (isConv_E) ? dt*he*dy*dz/(2*Ch*rho) : 0.0;

  for (k=1; k<Nz-1; k++){
    kc = Nx*Ny*k+j*Nx+i;
    kt = kc+Nx*Ny;
    kb = kc-Nx*Ny;
    kn = kc-Nx;
    kw = kc-1;

    M[kc] = w*( coef1*M[kw] + \
                coef2*M[kn] + \
                coef3*M[kt] + \
                coef3*M[kb] + \
                term + (1.0 - denom)*M[kc] ) + (1.0 - w)*M[kc];
  }
  // vertex 4
  j = Ny-1; k = Nz-1;

  coef1 = dt*kd/(dx*dx);
  coef2 = dt*kd/(dy*dy);
  coef3 = dt*kd/(dz*dz);

  term = 0.0;
  term += (isConv_S)  ? hs*dz*dx*Ts/(2*Ch*rho)  : 0.0;
  term += (isConv_T)  ? ht*dy*dx*Tt/(2*Ch*rho)  : 0.0;
  term *= dt;

  denom = coef1 + coef2 + coef3;
  denom += (isConv_S) ? dt*hs*dz*dx/(2*Ch*rho) : 0.0;
  denom += (isConv_T) ? dt*ht*dy*dx/(2*Ch*rho) : 0.0;

  for (i=1; i<Nx-1; i++){
    kc = Nx*Ny*k+j*Nx+i;
    kb = kc-Nx*Ny;
    kn = kc-Nx;
    kw = kc-1;
    ke = kc+1;

    M[kc] = w*( coef1*M[kw] + \
                coef1*M[ke] + \
                coef2*M[kn] + \
                coef3*M[kb] + \
                term + (1.0 - denom)*M[kc] ) + (1.0 - w)*M[kc];
  }
  // vertex 5
  i = 0; k = 0;

  coef1 = dt*kd/(dx*dx);
  coef2 = dt*kd/(dy*dy);
  coef3 = dt*kd/(dz*dz);

  term = 0.0;
  term += (isConv_W)  ? hw*dz*dy*Tw/(2*Ch*rho)  : 0.0;
  term += (isConv_B)  ? hb*dx*dy*Tb/(2*Ch*rho)  : 0.0;
  term *= dt;

  denom = coef1 + coef2 + coef3;
  denom += (isConv_W) ? dt*hw*dz*dy/(2*Ch*rho) : 0.0;
  denom += (isConv_B) ? dt*hb*dx*dy/(2*Ch*rho) : 0.0;

  for (j=1; j<Ny-1; j++){
    kc = Nx*Ny*k+j*Nx+i;
    kt = kc+Nx*Ny;
    ks = kc+Nx;
    kn = kc-Nx;
    ke = kc+1;

    M[kc] = w*( coef1*M[ke] + \
                coef2*M[kn] + \
                coef2*M[ks] + \
                coef3*M[kt] + \
                term + (1.0 - denom)*M[kc] ) + (1.0 - w)*M[kc];
  }
  // vertex 6
  i = Nx-1; k = 0;

  coef1 = dt*kd/(dx*dx);
  coef2 = dt*kd/(dy*dy);
  coef3 = dt*kd/(dz*dz);

  term = 0.0;
  term += (isConv_E)  ? he*dz*dy*Te/(2*Ch*rho)  : 0.0;
  term += (isConv_B)  ? hb*dx*dy*Tb/(2*Ch*rho)  : 0.0;
  term *= dt;

  denom = coef1 + coef2 + coef3;
  denom += (isConv_E) ? dt*ht*dz*dy/(2*Ch*rho) : 0.0;
  denom += (isConv_B) ? dt*hb*dx*dy/(2*Ch*rho) : 0.0;

  for (j=1; j<Ny-1; j++){
    kc = Nx*Ny*k+j*Nx+i;
    kt = kc+Nx*Ny;
    ks = kc+Nx;
    kn = kc-Nx;
    kw = kc-1;

    M[kc] = w*( coef1*M[kw] + \
                coef2*M[kn] + \
                coef2*M[ks] + \
                coef3*M[kt] + \
                term + (1.0 - denom)*M[kc] ) + (1.0 - w)*M[kc];
  }
  // vertex 7
  i = Nx-1; k = Nz-1;

  coef1 = dt*kd/(dx*dx);
  coef2 = dt*kd/(dy*dy);
  coef3 = dt*kd/(dz*dz);

  term = 0.0;
  term += (isConv_T)  ? ht*dx*dy*Tt/(2*Ch*rho)  : 0.0;
  term += (isConv_N)  ? hn*dz*dy*Tn/(2*Ch*rho)  : 0.0;
  term *= dt;

  denom = coef1 + coef2 + coef3;
  denom += (isConv_T) ? dt*ht*dx*dy/(2*Ch*rho) : 0.0;
  denom += (isConv_N) ? dt*hn*dz*dy/(2*Ch*rho) : 0.0;

  for (j=1; j<Ny-1; j++){
    kc = Nx*Ny*k+j*Nx+i;
    kb = kc-Nx*Ny;
    ks = kc+Nx;
    kn = kc-Nx;
    kw = kc-1;

    M[kc] = w*( coef1*M[kw] + \
                coef2*M[kn] + \
                coef2*M[ks] + \
                coef3*M[kb] + \
                term + (1.0 - denom)*M[kc] ) + (1.0 - w)*M[kc];
  }
  // vertex 8
  i = 0; k = Nz-1;

  coef1 = dt*kd/(dx*dx);
  coef2 = dt*kd/(dy*dy);
  coef3 = dt*kd/(dz*dz);

  term = 0.0;
  term += (isConv_W)  ? hw*dz*dy*Tw/(2*Ch*rho)  : 0.0;
  term += (isConv_T)  ? ht*dx*dy*Tt/(2*Ch*rho)  : 0.0;
  term *= dt;

  denom = coef1 + coef2 + coef3;
  denom += (isConv_W) ? dt*hw*dz*dy/(2*Ch*rho) : 0.0;
  denom += (isConv_T) ? dt*ht*dx*dy/(2*Ch*rho) : 0.0;

  for (j=1; j<Ny-1; j++){
    kc = Nx*Ny*k+j*Nx+i;
    kb = kc-Nx*Ny;
    ks = kc+Nx;
    kn = kc-Nx;
    ke = kc+1;

    M[kc] = w*( coef1*M[ke] + \
                coef2*M[kn] + \
                coef2*M[ks] + \
                coef3*M[kb] + \
                term + (1.0 - denom)*M[kc] ) + (1.0 - w)*M[kc];
  }
  // vertex 9
  i = 0; j = 0;

  coef1 = dt*kd/(dx*dx);
  coef2 = dt*kd/(dy*dy);
  coef3 = dt*kd/(dz*dz);

  term = 0.0;
  term += (isConv_N)  ? hn*dx*dz*Tn/(2*Ch*rho)  : 0.0;
  term += (isConv_W)  ? hw*dy*dz*Tw/(2*Ch*rho)  : 0.0;
  term *= dt;

  denom = coef1 + coef2 + coef3;
  denom += (isConv_N) ? dt*hn*dx*dz/(2*Ch*rho) : 0.0;
  denom += (isConv_W) ? dt*hw*dy*dz/(2*Ch*rho) : 0.0;

  for (k=1; k<Nz-1; k++){
    kc = Nx*Ny*k+j*Nx+i;
    kt = kc+Nx*Ny;
    kb = kc-Nx*Ny;
    ks = kc+Nx;
    ke = kc+1;

    M[kc] = w*( coef1*M[ke] + \
                coef2*M[ks] + \
                coef3*M[kb] + \
                coef3*M[kt] + \
                term + (1.0 - denom)*M[kc] ) + (1.0 - w)*M[kc];
  }
  // vertex 10
  j = 0; k = 0;

  coef1 = dt*kd/(dx*dx);
  coef2 = dt*kd/(dy*dy);
  coef3 = dt*kd/(dz*dz);

  term = 0.0;
  term += (isConv_N)  ? hn*dz*dx*Tn/(2*Ch*rho)  : 0.0;
  term += (isConv_B)  ? hb*dy*dx*Tb/(2*Ch*rho)  : 0.0;
  term *= dt;

  denom = coef1 + coef2 + coef3;
  denom += (isConv_N) ? dt*hn*dz*dx/(2*Ch*rho) : 0.0;
  denom += (isConv_B) ? dt*hb*dy*dx/(2*Ch*rho) : 0.0;

  for (i=1; i<Nx-1; i++){
    kc = Nx*Ny*k+j*Nx+i;
    kt = kc+Nx*Ny;
    ks = kc+Nx;
    ke = kc+1;
    kw = kc-1;

    M[kc] = w*( coef1*M[kw] + \
                coef1*M[ke] + \
                coef2*M[ks] + \
                coef3*M[kt] + \
                term + (1.0 - denom)*M[kc] ) + (1.0 - w)*M[kc];
  }
  // vertex 11
  i = Nx-1; j = 0;

  coef1 = dt*kd/(dx*dx);
  coef2 = dt*kd/(dy*dy);
  coef3 = dt*kd/(dz*dz);

  term = 0.0;
  term += (isConv_N)  ? hn*dx*dz*Tn/(2*Ch*rho)  : 0.0;
  term += (isConv_E)  ? he*dy*dz*Te/(2*Ch*rho)  : 0.0;
  term *= dt;

  denom = coef1 + coef2 + coef3;
  denom += (isConv_N) ? dt*hn*dx*dz/(2*Ch*rho) : 0.0;
  denom += (isConv_E) ? dt*he*dy*dz/(2*Ch*rho) : 0.0;

  for (k=1; k<Nz-1; k++){
    kc = Nx*Ny*k+j*Nx+i;
    kt = kc+Nx*Ny;
    kb = kc-Nx*Ny;
    ks = kc+Nx;
    kw = kc-1;

    M[kc] = w*( coef1*M[kw] + \
                coef2*M[ks] + \
                coef3*M[kt] + \
                coef3*M[kb] + \
                term + (1.0 - denom)*M[kc] ) + (1.0 - w)*M[kc];
  }
  // vertex 12
  j = 0; k = Nz-1;

  coef1 = dt*kd/(dx*dx);
  coef2 = dt*kd/(dy*dy);
  coef3 = dt*kd/(dz*dz);

  term = 0.0;
  term += (isConv_N)  ? hn*dz*dx*Tn/(2*Ch*rho)  : 0.0;
  term += (isConv_T)  ? ht*dy*dx*Tt/(2*Ch*rho)  : 0.0;
  term *= dt;

  denom = coef1 + coef2 + coef3;
  denom += (isConv_N) ? dt*hn*dz*dx/(2*Ch*rho) : 0.0;
  denom += (isConv_T) ? dt*ht*dy*dx/(2*Ch*rho) : 0.0;

  for (i=1; i<Nx-1; i++){
    kc = Nx*Ny*k+j*Nx+i;
    kb = kc-Nx*Ny;
    ks = kc+Nx;
    kw = kc-1;
    ke = kc+1;

    M[kc] = w*( coef1*M[kw] + \
                coef1*M[ke] + \
                coef2*M[ks] + \
                coef3*M[kb] + \
                term + (1.0 - denom)*M[kc] ) + (1.0 - w)*M[kc];
  }

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
  double coef1, coef2, coef3;
  double denom, cond0 = kd*( 1.0/(dx*dx) + 1.0/(dy*dy) + 1.0/(dz*dz) );

  assert( isConv_S == 0 || isConv_S == 1 );
  assert( isConv_N == 0 || isConv_N == 1 );
  assert( isConv_W == 0 || isConv_W == 1 );
  assert( isConv_E == 0 || isConv_E == 1 );
  assert( isConv_T == 0 || isConv_T == 1 );
  assert( isConv_B == 0 || isConv_B == 1 );

  Coef_ij = (dx/2.0)*(dy/2.0)/(Ch*rho); // convection constant for ij-plane
  Coef_ik = (dx/2.0)*(dz/2.0)/(Ch*rho); // convection constant for ik-plane
  Coef_jk = (dy/2.0)*(dz/2.0)/(Ch*rho); // convection constant for jk-plane

  coef1 = dt*kd/(dx*dx);
  coef2 = dt*kd/(dy*dy);
  coef3 = dt*kd/(dz*dz);

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
  term *= dt;

  denom = 0.0;
  denom += (isConv_W) ? hw*Coef_jk  : 0.0;
  denom += (isConv_B) ? hb*Coef_ij  : 0.0;
  denom += (isConv_S) ? hs*Coef_ik  : 0.0;
  denom += cond0;
  denom *= dt;

  M[kc] = w*( coef1*M[ke] + \
              coef2*M[kn] + \
              coef3*M[kt] + \
              term + (1.0 - denom)*M[kc] ) + (1.0 - w)*M[kc];

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
  term *= dt;

  denom = 0.0;
  denom += (isConv_S) ? hs*Coef_ik  : 0.0;
  denom += (isConv_E) ? he*Coef_jk  : 0.0;
  denom += (isConv_B) ? hb*Coef_ij  : 0.0;
  denom += cond0;
  denom *= dt;

  M[kc] = w*( coef1*M[kw] + \
              coef2*M[kn] + \
              coef3*M[kt] + \
              term + (1.0 - denom)*M[kc] ) + (1.0 - w)*M[kc];

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
  term *= dt;

  denom = 0.0;
  denom += (isConv_S) ? hs*Coef_ik  : 0.0;
  denom += (isConv_W) ? hw*Coef_jk  : 0.0;
  denom += (isConv_T) ? ht*Coef_ij  : 0.0;
  denom += cond0;
  denom *= dt;

  M[kc] = w*( coef1*M[ke] + \
              coef2*M[kn] + \
              coef3*M[kb] + \
              term + (1.0 - denom)*M[kc] ) + (1.0 - w)*M[kc];

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
  term *= dt;

  denom = 0.0;
  denom += (isConv_S) ? hs*Coef_ik  : 0.0;
  denom += (isConv_E) ? he*Coef_jk  : 0.0;
  denom += (isConv_T) ? ht*Coef_ij  : 0.0;
  denom += cond0;
  denom *= dt;

  M[kc] = w*( coef1*M[kw] + \
              coef2*M[kn] + \
              coef3*M[kb] + \
              term + (1.0 - denom)*M[kc] ) + (1.0 - w)*M[kc];

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
  term *= dt;

  denom = 0.0;
  denom += (isConv_W) ? hw*Coef_jk  : 0.0;
  denom += (isConv_N) ? hn*Coef_ik  : 0.0;
  denom += (isConv_B) ? hb*Coef_ij  : 0.0;
  denom += cond0;
  denom *= dt;

  M[kc] = w*( coef1*M[ke] + \
              coef2*M[ks] + \
              coef3*M[kt] + \
              term + (1.0 - denom)*M[kc] ) + (1.0 - w)*M[kc];

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
  term *= dt;

  denom = 0.0;
  denom += (isConv_E) ? he*Coef_jk  : 0.0;
  denom += (isConv_N) ? hn*Coef_ik  : 0.0;
  denom += (isConv_B) ? hb*Coef_ij  : 0.0;
  denom += cond0;
  denom *= dt;

  M[kc] = w*( coef1*M[kw] + \
              coef2*M[ks] + \
              coef3*M[kt] + \
              term + (1.0 - denom)*M[kc] ) + (1.0 - w)*M[kc];

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
  term *= dt;

  denom = 0.0;
  denom += (isConv_W) ? hw*Coef_jk  : 0.0;
  denom += (isConv_N) ? hn*Coef_ik  : 0.0;
  denom += (isConv_T) ? ht*Coef_ij  : 0.0;
  denom += cond0;
  denom *= dt;

  M[kc] = w*( coef1*M[ke] + \
              coef2*M[ks] + \
              coef3*M[kb] + \
              term + (1.0 - denom)*M[kc] ) + (1.0 - w)*M[kc];

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
  term *= dt;

  denom = 0.0;
  denom += (isConv_E) ? he*Coef_jk  : 0.0;
  denom += (isConv_N) ? hn*Coef_ik  : 0.0;
  denom += (isConv_T) ? ht*Coef_ij  : 0.0;
  denom += cond0;
  denom *= dt;

  M[kc] = w*( coef1*M[kw] + \
              coef2*M[ks] + \
              coef3*M[kb] + \
              term + (1.0 - denom)*M[kc] ) + (1.0 - w)*M[kc];


}


void Dirichlet(int surface, double value)
{
  int i, j, k, kc, fix, plane;
  assert(surface == 0 || surface == 1 || \
         surface == 2 || surface == 3 || \
         surface == 4 || surface == 5 );
#define X 0 //  plane perpendicular to i-axis
#define Y 1 //  plane perpendicular to j-axis
#define Z 2 //  plane perpendicular to k-axis

  switch(surface) {
    case(0): // case, top plane
      plane = Z; fix = Nz-1;
      break;
    case(1): // case, bottom plane
      plane = Z; fix = 0;
      break;
    case(2): // case, east plane
      plane = X; fix = Nx-1;
      break;
    case(3): // case, west plane
      plane = X; fix = 0;
      break;
    case(4): // case, south plane
      plane = Y; fix = Ny-1;
      break;
    case(5): // case, north plane
      plane = Y; fix = 0;
      break;
    default:
      printf("\nWrong input in Dirichlet function %d", surface);
      exit(1);
  }

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


void surface(int surface, double h, double Tsurr)
{
  /*
    This only takes care of convection/free-surface in the inner surface
    of a given plane -- i.e. excluding corners and vertices.
  */
  int i, j, k, kw, ke, kt, kb, ks, kn, kc, fix, plane;
  double Ccp, Cxp, Cyp, Czp, Csp;
  assert(surface == 0 || surface == 1 || \
         surface == 2 || surface == 3 || \
         surface == 4 || surface == 5 );
#define X 0 // plane perpendicular to i-axis
#define Y 1 // plane perpendicular to j-axis
#define Z 2 // plane perpendicular to k-axis

  switch(surface) {
    case(0): // case, top plane
      plane = Z; fix = Nz-1;
      break;
    case(1): // case, bottom plane
      plane = Z; fix = 0;
      break;
    case(2): // case, east plane
      plane = X; fix = Nx-1;
      break;
    case(3): // case, west plane
      plane = X; fix = 0;
      break;
    case(4): // case, south plane
      plane = Y; fix = Ny-1;
      break;
    case(5): // case, north plane
      plane = Y; fix = 0;
      break;
    default:
      printf("\nWrong input in surface function %d", surface);
      exit(1);
  }

  switch (plane) {

    case X:
      i = fix; // at which level to fix
      assert( i==0 || i==(Nx-1) );

      Ccp = dt*( h*dy*dz/(Ch*rho) + kd*( 1.0/(dx*dx) + 1.0/(dy*dy) + 1.0/(dz*dz) ) );
      Cxp = dt*kd/(dx*dx);
      Cyp = dt*kd/(dy*dy);
      Czp = dt*kd/(dz*dz);
      Csp = dt*h*dy*dz/(Ch*rho);

      if (i==Nx-1){
        for (k=1; k<Nz-1; k++){
          for (j=1; j<Ny-1; j++){
            kc = Nx*Ny*k+j*Nx+i;
            kt = kc+Nx*Ny; kb = kc-Nx*Ny;
            ks = kc+Nx;    kn = kc-Nx;
            kw = kc-1;
            M[kc] = w*( Cxp*M[kw] + Cyp*(M[kn]+M[ks]) + \
                        Czp*(M[kb]+M[kt]) + Csp*Tsurr + \
                        (1.0 - Ccp)*M[kc] ) + \
                        (1.0 - w)*M[kc];
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
            M[kc] = w*( Cxp*M[ke] + Cyp*(M[kn]+M[ks]) + \
                        Czp*(M[kb]+M[kt]) + Csp*Tsurr + \
                        (1.0 - Ccp)*M[kc] ) + \
                        (1.0 - w)*M[kc];
          }
        }
      }
      break;

    case Y:
      j = fix; // at which level to fix
      assert( j==0 || j==(Ny-1) );

      Ccp = dt*( h*dx*dz/(Ch*rho) + kd*( 1.0/(dx*dx) + 1.0/(dy*dy) + 1.0/(dz*dz) ) );
      Cxp = dt*kd/(dx*dx);
      Cyp = dt*kd/(dy*dy);
      Czp = dt*kd/(dz*dz);
      Csp = dt*h*dx*dz/(Ch*rho);

      if (j==Ny-1){
        for (k=1; k<Nz-1; k++){
          for (i=1; i<Nx-1; i++){
            kc = Nx*Ny*k+j*Nx+i;
            kt = kc+Nx*Ny; kb = kc-Nx*Ny;
            ke = kc+1;     kw = kc-1;
            kn = kc-Nx;
            M[kc] = w*( Cxp*(M[ke]+M[kw]) + Cyp*M[kn] + \
                        Czp*(M[kb]+M[kt]) + Csp*Tsurr + \
                        (1.0 - Ccp)*M[kc] ) + \
                        (1.0 - w)*M[kc];
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
            M[kc] = w*( Cxp*(M[ke]+M[kw]) + Cyp*M[ks] + \
                        Czp*(M[kb]+M[kt]) + Csp*Tsurr + \
                        (1.0 - Ccp)*M[kc] ) + \
                        (1.0 - w)*M[kc];
          }
        }
      }
      break;

    case Z:
      k = fix; // at which level to fix
      assert( k==0 || k==(Nz-1) );

      Ccp = dt*( h*dx*dy/(Ch*rho) + kd*( 1.0/(dx*dx) + 1.0/(dy*dy) + 1.0/(dz*dz) ) );
      Cxp = dt*kd/(dx*dx);
      Cyp = dt*kd/(dy*dy);
      Czp = dt*kd/(dz*dz);
      Csp = dt*h*dx*dy/(Ch*rho);

      if (k==Nz-1){
        for (j=1; j<Ny-1; j++){
          for (i=1; i<Nx-1; i++){
            kc = Nx*Ny*k+j*Nx+i;
            ks = kc+Nx;    kn = kc-Nx;
            ke = kc+1;     kw = kc-1;
            kb = kc-Nx*Ny;
            M[kc] = w*( Cxp*(M[kw]+M[ke]) + Cyp*(M[kn]+M[ks]) + \
                        Czp*M[kb] + Csp*Tsurr + \
                        (1.0 - Ccp)*M[kc] ) + \
                        (1.0 - w)*M[kc];
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
            M[kc] = w*( Cxp*(M[kw]+M[ke]) + Cyp*(M[kn]+M[ks]) + \
                        Czp*M[kt] + Csp*Tsurr + \
                        (1.0 - Ccp)*M[kc] ) + \
                        (1.0 - w)*M[kc];
          }
        }
      }
      break;

    default:
      printf("\nWrong input in 'convective function' \t plane: %d\n\n", plane);
      exit(1);
  }

}


void insulated(int surface)
{
  /*
    This only takes care of insulation in the inner surface of a given plane
    i.e. excluding corners and vertices.
  */
  int i, j, k, step, plane, fix;
  double temp = 0.0;
  assert(surface == 0 || surface == 1 || \
         surface == 2 || surface == 3 || \
         surface == 4 || surface == 5 );
#define X 0  // plane perpendicular to i-axis
#define Y 1  // plane perpendicular to j-axis
#define Z 2  // plane perpendicular to k-axis

  switch(surface) {
    case(0): // case, top plane
      plane = Z; fix = Nz-1;
      break;
    case(1): // case, bottom plane
      plane = Z; fix = 0;
      break;
    case(2): // case, east plane
      plane = X; fix = Nx-1;
      break;
    case(3): // case, west plane
      plane = X; fix = 0;
      break;
    case(4): // case, south plane
      plane = Y; fix = Ny-1;
      break;
    case(5): // case, north plane
      plane = Y; fix = 0;
      break;
    default:
      printf("\nWrong input in insulated function %d", surface);
      exit(1);
  }

  switch(plane){

    case X:
      i = fix; // at which level to fix
      assert( i==0 || i==(Nx-1) );
      step = (i == 0)   ? 1 : -1;
      for (k=1; k<Nz-1; k++){
        for (j=1; j<Ny-1; j++){
          temp = M[Nx*Ny*k+j*Nx+i+step] - M[Nx*Ny*k+j*Nx+i];
          M[Nx*Ny*k+j*Nx+i] += 2.0*kd*dt*temp/(dx*dx);
        }
      }
      break;
    case Y:
      j = fix; // at which level to fix
      assert( j==0 || j==(Ny-1) );
      step = (j == 0)   ? 1 : -1;
      for (k=1; k<Nz-1; k++){
        for (i=1; i<Nx-1; i++){
          temp = M[Nx*Ny*k+(j+step)*Nx+i] - M[Nx*Ny*k+j*Nx+i];
          M[Nx*Ny*k+j*Nx+i] += 2.0*kd*dt*temp/(dy*dy);
        }
      }
      break;
    case Z:
      k = fix; // at which level to fix
      assert( k==0 || k==(Nz-1) );
      step = (k == 0)   ? 1 : -1;
      for (j=1; j<Ny-1; j++){
        for (i=1; i<Nx-1; i++){
          temp = M[Nx*Ny*(k+step)+j*Nx+i] - M[Nx*Ny*k+j*Nx+i];
          M[Nx*Ny*k+j*Nx+i] += 2.0*kd*dt*temp/(dz*dz);
        }
      }
      break;

    default:
      printf("\nWrong input in 'insulation function' \t plane: %d\n\n", plane);
      exit(1);
  }

}
