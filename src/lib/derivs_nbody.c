#include "nbody.h"
#include <math.h>

// leapfrog derivatives function
void l_derivs(double t, double y[], double d2ydt2[]) {
  // x**_i = -G*m_j*(x_i - x_j)/(|r_i - r_j|**2 + eta**2)**(3/2)
  // y**_i = -G*m_i*(y_i - y_j)/(|r_i - r_j|**2 + eta**2)**(3/2)
  // z...
  int i, j, k, k_i, k_j;
  double num, den, m_i, m_j;
  // outer loop loops through each point (spans 4 elements)
  for(i = 1; i <= n; i++) {
    // inner loop loops through points as well
    for(j = 1; j <= n; j++) {
      // sets all of the deriv to zero to start
      if(i == 1) for(k = 1; k <= 4; k++) d2ydt2[l_i(j,k)] = 0.0;
      // will skip a point from attracting itself, as well as previously
      // calculated points
      if(i>=j) continue;
      m_i = y[l_i(i,4)]; // mass of i
      m_j = y[l_i(j,4)]; // mass of j
      for(k = 1; k <= 3; k++) {
	k_i = l_i(i,k); // x/y/z of i
	k_j = l_i(j,k); // x/y/z of j
	// Newton's equation in two parts
	num = -G*(y[k_i] - y[k_j]);
	den = pow( pow(y[l_i(i,1)] - y[l_i(j,1)], 2.0)
		   + pow(y[l_i(i,2)] - y[l_i(j,2)], 2.0)
		   + pow(y[l_i(i,3)] - y[l_i(j,3)], 2.0)
		   + pow(eta, 2.0), 1.5);
	d2ydt2[k_i] += m_j*num/den;
	d2ydt2[k_j] -= m_i*num/den;
      }
    }
  }
}

// a simple magnitude function taking an array and a counter
double mag(double a[], int c) {
  double val;
  int i;
  for(i = 0; i < c; i++) val += sqr(a[i]);
  val = sqrt(val);
  return val;
}

// finding the magntiude of the difference of two arrays
double mag_diff(double a[], double b[], int c) {
  double val = 0.0;
  int i;
  for(i = 0; i < c; i++)
    val += sqr(a[i] - b[i]);
  val = sqrt(val);
  return val;
}

// a square function for increased peformance over pow()
double sqr(double a) {
  return a*a;
}

// the bh-mode derivs function
void derivs_pt(double t, PARTICLE pts[]) {
  int i;
  // fill force
  NODE *root;
  root = build_tree(pts);
  
  for(i = 0; i < n; i++) {
    add_force(&(pts[i]), root);
  }

  free_tree(root);
}
