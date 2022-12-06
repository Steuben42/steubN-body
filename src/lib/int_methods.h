#ifndef INT_METHODS
#define INT_METHODS

#include "nbody.h"

// arguments interpreted identically to rk4.c
void euler_std(double y[], double dydx[], int n, double x, double h,
	       double yout[], void(*derivs)(double, double [], double []) );

// arguments interpreted identically to rk4.c with addition of d2ydx2 due to
// this being a second order ODE solver
void leap_std(double y[], double dydx[], double d2ydx2[], int n, double x,
	      double h, double yout[], void (*derivs2)(double, double [],
						       double []) );
// leap using the particle type
void leap_pt(PARTICLE *pts, int n, double t, double h,
	     void (*derivs_pt)(double, PARTICLE []));

#endif
