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


int main(int argc, char **argv){

  /* input solver choise */
  int solver = atoi(argv[1]);


  Nx = 50;
  Ny = 50;
  Nz = 50;

  rho = 1.0;
  Ch  = 1.0;
  h   = 20.0;
  kd  = 1.0/(rho*Ch);
  dx  = 1.0/(Nx-1);
  dy  = 1.0/(Ny-1);
  dz  = 1.0/(Nz-1);
  dt  = 0.00001;
  Tsurr = 370.0;  // Kelvin



  initialize();
  boundary();

  switch(solver){
    // use Gauss-Seidel
    case 0:
      GS();
      break;

    // use Successive Over-Relaxation
    case 1:
      SOR();
      break;

    // use conjugate gradient
    case 2:
      CG();
      break;

    // use multigrid
    case 3:
      MG();
      break;

    // use FEM
    case 4:
      FEM();
      break;

    // default - gauss elimination
    default:
      DS();
  }

  FILE *f = fopen("file.txt", "w");
  if (f == NULL){
    printf("\nError opening file ! " );
    exit(1);
  }
  int i,j,k;
  for (k=0; k<Nz; k++){
    for (j=0; j<Ny; j++){
      for (i=0; i<Nx; i++){
        fprintf(f, "%2.2f,", M[Nx*Ny*k+j*Nx+i] );
      }
      fprintf(f, "\n" );
    }
    fprintf(f, "\n" );
  }


  free(M);
  printf("\n");
  return 0;
}
