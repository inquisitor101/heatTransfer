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

#ifndef _HEADER_H
#define _HEADER_H

#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#define isRegular 0
#define isSourceTerm 0

// @TODO: check which constants can be macros or 'const'
//
int Nx, Ny, Nz; // discretization points per axis
double dx, dy, dz, dt, kd, Ch, rho, w;
double *M;

int boundCond[6];
double h_conv[6];
double Tsurr[6];
double ems[6];
double fixedTemp[6];
double Lx, Ly, Lz;
double Kcond;

double posX, posY, posZ;
double length, width, height;
int x0, x1, y0, y1, z0, z1;
double genHeat, simulationTime;

double sigma; // stefan-boltzmann constant
#endif
