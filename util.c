
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
  //       fixed bounary as a function
  //       corner + vertices function

  // @IDEA:
  //       radiation function 
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
