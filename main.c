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
// #define isRegular 1
// #define isSourceTerm 1

#include "header.h"
#include "funcs.h"
#include "util.h"


int main(int argc, char **argv){


  Nx = 5;  Lx = 1.0; // meters
  Ny = 5;  Ly = 1.0; // meters
  Nz = 3;  Lz = 1.0; // meters
  simulationTime = 100; // total time steps

  // SOR coefficient
  w = 1.9;

  // @FIXME: needs fixing, debug : rho & Ch -- or maybe it is correct ??
  //
  // material: concrete
  rho = 1;//2400.0; // units: Kg/m3
  Ch  = 1;//750.0; // units: J/(K.Kg)
  Kcond = 1.0; // thermal conductivity, units: W/(m.K)


  kd  = Kcond/(rho*Ch); // units: m2/sec
  dx  = Lx/(Nx-1); // units: m
  dy  = Ly/(Ny-1); // units: m
  dz  = Lz/(Nz-1); // units: m
  dt  = 0.001; // units: seconds

  // source term and geometry
  posX = 0.5; length = 0.1; // unitless -- relative
  posY = 0.5; width  = 0.1; // unitless -- relative
  posZ = 0.5; height = 0.1; // unitless -- relative
  genHeat = 1000.0; // units: W/m3

  // radiation terms
  sigma = 5.67e-8; // stefan-boltzmann constant, units: W/(m2.K4)
  ems = 0.91; // emissivity constant concrete (rough)

  // @TODO: change boundary conditions into input if/else
  //

  // legend: 0: T, 1: B, 2: E, 3: W, 4: S, 5: N
  // 1: convection/free, 2: insulation, 3: Dirichlet
  boundCond[0] = 3; // T
  boundCond[1] = 2; // B
  boundCond[2] = 2; // E
  boundCond[3] = 2; // W
  boundCond[4] = 2; // S
  boundCond[5] = 2; // N

  // legend: 0: T, 1: B, 2: E, 3: W, 4: S, 5: N
  // 0: no convection OR free, #: convective coefficient
  // units: W/(m2.K)
  h_conv[0] = 20.0; // T
  h_conv[1] = 20.0; // B
  h_conv[2] = 20.0; // E
  h_conv[3] = 20.0; // W
  h_conv[4] = 20.0; // S
  h_conv[5] = 20.0; // N

  // legend: 0: T, 1: B, 2: E, 3: W, 4: S, 5: N
  // temperature in Kelvin (degC + 273.15)
  Tsurr[0] = 273.15+100; // T
  Tsurr[1] = 273.15+100; // B
  Tsurr[2] = 273.15+100; // E
  Tsurr[3] = 273.15+100; // W
  Tsurr[4] = 273.15+100; // S
  Tsurr[5] = 273.15+100; // N

  // legend: 0: T, 1: B, 2: E, 3: W, 4: S, 5: N
  // surface temperature in Kelvin (degC + 273.15)
  fixedTemp[0] = 273.15+100; // T
  fixedTemp[1] = 273.15+100; // B
  fixedTemp[2] = 273.15+100; // E
  fixedTemp[3] = 273.15+100; // W
  fixedTemp[4] = 273.15+100; // S
  fixedTemp[5] = 273.15+100; // N


  initialize();

  // iterate via Successive Over-Relaxation
  SOR();


  FILE *f = fopen("file.txt", "w");
  if (f == NULL){
    printf("\nError opening file ! " );
    exit(1);
  }
  int i,j,k;
  for (k=0; k<Nz; k++){
    for (j=0; j<Ny; j++){
      for (i=0; i<Nx; i++){
        fprintf(f, "%3.2f,", M[Nx*Ny*k + j*Nx + i] );
      }
      fprintf(f, "\n" );
    }
    fprintf(f, "\n" );
  }

  // ----
  printf("\n\n\n" );
  for (k=0; k<Nz; k++){
    for (j=0; j<Ny; j++){
      for (i=0; i<Nx; i++){
        // array math: array[M*L*k + L*j + i]
        printf("%3.0f ", M[Nx*Ny*k + j*Nx + i]-273.15);
      }
      printf("\n" );
    }
    printf("\n" );
  }
  printf("\n" );

  // printf("\n       \t\t lev 0\n  B B  \n  B B  \n       \n" );
  // printf("\n  N N  \t\t lev 1\nW     E\nW     E\n  S S  \n" );
  // printf("\n  N N  \t\t lev 2\nW     E\nW     E\n  S S  \n" );
  // printf("\n       \t\t lev 3\n  T T  \n  T T  \n       \n" );
  // // ----

  free(M);
  printf("\n");
  return 0;
}
