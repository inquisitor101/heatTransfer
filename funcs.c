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
#include "funcs.h"
#include "util.h"

static inline double sourceGenIrreg(int x, int y, int z)
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

  double rX, rY, rZ;

  rX = x*dx-posX*Lx;
  rY = y*dy-posY*Ly;
  rZ = z*dz-posZ*Lz;

  double dist = sqrt(rX*rX + rY*rY + rZ*rZ);
  double gen = (dist <= length) ? genHeat : 0.0; // units: W/m3

  return gen;
}


void SOR(void)
{
  double sum, err = 1.0, TOL = 1e-6, w = 1.9;
  int i, j, k, iter = 0, maxIter = simulationTime;
  int kc, kn, ks, kt, kb, ke, kw;
  double Cx = kd*dt/(dx*dx);
  double Cy = kd*dt/(dy*dy);
  double Cz = kd*dt/(dz*dz);
  double Cc = 2*kd*dt*( 1.0/(dx*dx) + 1.0/(dy*dy) + 1.0/(dz*dz) );
  double Cs = dt/(rho*Ch); // this transforms (W/m3) into (Kelvin)

  double *Mold = (double *)calloc(Nx*Ny*Nz, sizeof(double));
  memcpy(Mold, M, Nx*Ny*Nz*sizeof(double));

#if isRegular
  // regular --> box
  while (iter < maxIter && err > TOL){

    for (k=1; k<Nz-1; k++){
      for (j=1; j<Ny-1; j++){
        for (i=1; i<Nx-1; i++){
          kc = Nx*Ny*k+j*Nx+i;
          kt = kc+Nx*Ny; kb = kc-Nx*Ny;
          ks = kc+Nx;    kn = kc-Nx;
          ke = kc+1;     kw = kc-1;
          M[kc] = w*( Cx*( M[ke] + M[kw] ) \
                    + Cy*( M[ks] + M[kn] ) \
                    + Cz*( M[kt] + M[kb] ) \
                     + Cs*sourceGen(i,j,k) \
                      + (1.0 - Cc)*M[kc] ) \
                      + (1.0 - w )*M[kc];
        }
      }
    }

    boundary(boundCond, h_conv, Tsurr, fixedTemp);

    if (iter%50 == 0){
      // calculate norm of error
      sum = 0.0;
      for (k=0; k<Nz; k++){
        for (j=0; j<Ny; j++){
          for (i=0; i<Nx; i++){
            kc = Nx*Ny*k+j*Nx+i;
            sum += fabs(M[kc] - Mold[kc]);
          }
        }
      }
      err = sum/(Nx*Ny*Nz);
      // update Mold
      memcpy(Mold, M, Nx*Ny*Nz*sizeof(double));
      // monitor progress
      printf("\nstep: %d out of %d, with difference: %e", iter, maxIter, err);
    }

    iter++;
  }

#else
  // irregular --> sphere
  while (iter < maxIter && err > TOL){

    for (k=1; k<Nz-1; k++){
      for (j=1; j<Ny-1; j++){
        for (i=1; i<Nx-1; i++){
          kc = Nx*Ny*k+j*Nx+i;
          kt = kc+Nx*Ny; kb = kc-Nx*Ny;
          ks = kc+Nx;    kn = kc-Nx;
          ke = kc+1;     kw = kc-1;
          M[kc] = w*( Cx*( M[ke] + M[kw] ) \
                    + Cy*( M[ks] + M[kn] ) \
                    + Cz*( M[kt] + M[kb] ) \
                     + Cs*sourceGenIrreg(i,j,k) \
                      + (1.0 - Cc)*M[kc] ) \
                      + (1.0 - w )*M[kc];
        }
      }
    }

    boundary(boundCond, h_conv, Tsurr, fixedTemp);

    if (iter%50 == 0){
      // calculate norm of error
      sum = 0.0;
      for (k=0; k<Nz; k++){
        for (j=0; j<Ny; j++){
          for (i=0; i<Nx; i++){
            kc = Nx*Ny*k+j*Nx+i;
            sum += fabs(M[kc] - Mold[kc]);
          }
        }
      }
      err = sum/(Nx*Ny*Nz);
      // update Mold
      memcpy(Mold, M, Nx*Ny*Nz*sizeof(double));
      // monitor progress
      printf("\nstep: %d out of %d, with difference: %e", iter, maxIter, err);
    }

    iter++;
  }
#endif


  free(Mold);
}
