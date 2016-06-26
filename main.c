
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
