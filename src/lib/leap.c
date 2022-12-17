#include <math.h>

#include "nr.h"
#include "nrutil.h"
#include "int_methods.h"

// takes almost the same arguments as NRIC's rk4.c with the addition of d2ydx2
// before the n parameter (since this is a purely second order ODE solver)
void leap_std(double y[], double dydx[], double d2ydx2[], int n, double x,
	      double h, double yout[], void (*derivs2)(double, double [],
						       double [])) {
  int i;
  // will store the kicked derivative
  double *v;
  v = dvector(1, n);

  // evaluate second derivative
  (*derivs2)(x, y, d2ydx2);
  // kick (find derivative at half step)
  for(i = 1; i <= n; i++) {
    // euler approx. with half step
    v[i] = dydx[i] + (h/2)*d2ydx2[i];
  }
  // drift (full euler step with the half step derivative)
  for(i = 1; i <= n; i++) yout[i] = y[i] + h*v[i];
  // sync (another half step derivative to resync at t=t0+n)
  // evaluate second derivative at end point
  (*derivs2)(x+h, yout, d2ydx2);
  for(i = 1; i <= n; i++) {
    // euler approx. with half step to resync
    dydx[i] = v[i] + (h/2)*d2ydx2[i];
  }

  // save memory
  free_dvector(v, 1, n);
}

// leap adapted to take the particle type instead of an array
void leap_pt(PARTICLE *pts, int n, double t, double h,
	     void (*derivs_pt)(double, PARTICLE pt[])) {
  int i, k;
  
  for(i = 0; i < n; i++)     
    for(k = 0; k < 3; k++) 
      pts[i].vel[k] += (h/2.0)*pts[i].acc[k];
  
  for(i = 0; i < n; i++)
    for(k = 0; k < 3; k++)
      pts[i].pos[k] += h*pts[i].vel[k];

  (*derivs_pt)(t+h, pts);
  for(i = 0; i < n; i++)
    for(k = 0; k < 3; k++) 
      pts[i].vel[k] += (h/2.0)*pts[i].acc[k];
}
