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
#ifndef _UTIL_H
#define _UTIL_H


void initialize(void);

void insulated(int surface);

void heatGeneration(void);

void Dirichlet(int surface, double value);

void surface(int surface, double h, double Tsurr);

void boundary(int condition[], double h[],
              double Tsurr[], double temp[]);

void allCorners(int isConv_S, double hs, double Ts, // convection south -side
                int isConv_N, double hn, double Tn, // convection north -side
                int isConv_W, double hw, double Tw, // convection west  -side
                int isConv_E, double he, double Te, // convection east  -side
                int isConv_T, double ht, double Tt, // convection top   -side
                int isConv_B, double hb, double Tb);// convection bottom-side


void allVertices(int isConv_S, double hs, double Ts, // convection south -side
                 int isConv_N, double hn, double Tn, // convection north -side
                 int isConv_W, double hw, double Tw, // convection west  -side
                 int isConv_E, double he, double Te, // convection east  -side
                 int isConv_T, double ht, double Tt, // convection top   -side
                 int isConv_B, double hb, double Tb);// convection bottom-side


static inline double sourceGen(int x, int y, int z)
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

  double gen = (x>x0 && x<x1 && y>y0 && y<y1 && z>z0 && z<z1) ? genHeat  : 0.0; // units: W/m3

  return gen;
}



#endif
