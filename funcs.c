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


void SOR(void)
{
  double sum, err = 1.0, TOL = 1e-6, w = 1.9;
  int i, j, k, iter = 0, maxIter = simulationTime;
  int kc, kn, ks, kt, kb, ke, kw;
  double Cx = kd*dt/(dx*dx);
  double Cy = kd*dt/(dy*dy);
  double Cz = kd*dt/(dz*dz);
  double Cc = 2*kd*dt*( 1.0/(dx*dx) + 1.0/(dy*dy) + 1.0/(dz*dz) );

#if isSourceTerm
  double Cs = dt/(rho*Ch); // this transforms (W/m3) into (Kelvin)
#endif

  double *Mold = (double *)calloc(Nx*Ny*Nz, sizeof(double));
  memcpy(Mold, M, Nx*Ny*Nz*sizeof(double));

#if isSourceTerm
  // check if source term exists
  #if isRegular
    // regular --> box
    while (iter < maxIter){

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
      // check if steady-state is reached ?
      if (err < TOL){
        printf("\nSteady-state like conditions are reached.\nSimulation over.\n");
        break;
      }
    }

  #else
    // irregular --> sphere
    while (iter < maxIter){

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
      // check if steady-state is reached ?
      if (err < TOL){
        printf("\nSteady-state like conditions are reached.\nSimulation over.\n");
        break;
      }
    }
  #endif

#else
  // no source term
  while (iter < maxIter){

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
    // check if steady-state is reached ?
    if (err < TOL){
      printf("\nSteady-state like conditions are reached.\nSimulation over.\n");
      break;
    }
  }

#endif


  free(Mold);
}
