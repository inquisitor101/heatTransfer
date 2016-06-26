
#ifndef _UTIL_H
#define _UTIL_H


void initialize(void);

void boundary(void);

void insulated(int fix, int plane);

void convective(int fix, int plane,
                double h, double Tsurr);

#endif
