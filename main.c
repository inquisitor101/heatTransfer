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


  Nx = 50;  Lx = 10.0; // meters
  Ny = 50;  Ly = 10.0; // meters
  Nz = 50;  Lz = 10.0; // meters
  simulationTime = 1000; // total time steps

  // @FIXME: needs fixing, debug : rho & Ch
  //
  rho = 1.0;//2400.0; // units: Kg/m3
  Ch  = 1.0;//750.0; // units: J/(K.Kg)
  kd  = 1.0/(rho*Ch);
  dx  = Lx/(Nx-1);
  dy  = Ly/(Ny-1);
  dz  = Lz/(Nz-1);
  dt  = 0.00001; // seconds

  // source term and geometry
  posX = 0.5; length = 0.2; // unitless -- relative
  posY = 0.5; width  = 0.2; // unitless -- relative
  posZ = 0.5; height = 0.2; // unitless -- relative
  genHeat = 1000.0; // units: W/m3

  // radiation terms
  sigma = 5.67e-8; // units: W/(m2.K4)
  ems = 0.91; // concrete (rough)

  // legend: 0: T, 1: B, 2: E, 3: W, 4: S, 5: N
  // 1: convection/free, 2: insulation, 3: Dirichlet
  boundCond[0] = 2; // T
  boundCond[1] = 2; // B
  boundCond[2] = 1; // E
  boundCond[3] = 1; // W
  boundCond[4] = 2; // S
  boundCond[5] = 2; // N

  // legend: 0: T, 1: B, 2: E, 3: W, 4: S, 5: N
  // 0: no convection OR free, #: convective coefficient
  h_conv[0] = 0.00; // T
  h_conv[1] = 0.00; // B
  h_conv[2] = 20.0; // E
  h_conv[3] = 20.0; // W
  h_conv[4] = 0.00; // S
  h_conv[5] = 0.00; // N

  // legend: 0: T, 1: B, 2: E, 3: W, 4: S, 5: N
  // temperature in Kelvin (degC + 273.15)
  Tsurr[0] = 273.15+25; // T
  Tsurr[1] = 273.15+25; // B
  Tsurr[2] = 273.15+00; // E
  Tsurr[3] = 273.15+00; // W
  Tsurr[4] = 273.15+25; // S
  Tsurr[5] = 273.15+25; // N

  // legend: 0: T, 1: B, 2: E, 3: W, 4: S, 5: N
  // surface temperature in Kelvin (degC + 273.15)
  fixedTemp[0] = 273.15+25; // T
  fixedTemp[1] = 373.15+40; // B
  fixedTemp[2] = 273.15+55; // E
  fixedTemp[3] = 273.15+25; // W
  fixedTemp[4] = 273.15+25; // S
  fixedTemp[5] = 273.15+25; // N


  initialize();
  boundary(boundCond, h_conv, Tsurr, fixedTemp);

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

  free(M);
  printf("\n");
  return 0;
}
