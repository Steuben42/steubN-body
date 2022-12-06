#include "nbody.h"
#include <math.h>
#include <unistd.h>

void create2body(double ecc, unsigned int flags, double l[], double dl[], double u[]) {
  if(flags & DEBUG) printf("nbody_it.create2body@start\n");
  // start
  double x, y, z, vx, vy, vz, m;
  double v;

  // massflow
  if(flags & MASS_NORM) {
    // creates the initial conditions
    if(flags & DEBUG) printf("nbody_it.create2body->massflow->MASS_NORM\n");
    //r1 = (-1/2 0 0)
    //r2 = (1/2 0 0)
    //v1 = (0 -v/2 0)
    //v2 = (0 v/2 0)
    m = 1.0;
    x = 0.5;
    y = 0.0;
    z = 0.0;

    v = sqrt(G*2*(1 - ecc));
    vx = 0.0;
    vy = -v/2;
    vz = 0.0;
  }

  // assign result
  if(flags & DEBUG) printf("nbody_it.create2body@assign result\n");
  // leapfrog assigning
  if(l) {
    l[1] = x;
    l[2] = y;
    l[3] = z;
    l[4] = m;
  
    l[5] = -x;
    l[6] = -y;
    l[7] = -z;
    l[8] = m;

    dl[1] = vx;
    dl[2] = vy;
    dl[3] = vz;
    dl[4] = 0.0;

    dl[5] = -vx;
    dl[6] = -vy;
    dl[7] = -vz;
    dl[8] = 0.0;
  }
  // returns l in format of leapfrog, returns u in format of rk4

  // finish
  if(flags & DEBUG) printf("nbody_it.create2body@finish\n");
}

void create_sphere_simple(int n, unsigned int flags, double l[], double dl[],
			  double u[]) {
  double R = 0.5;
  long idum = -1*((long) getpid());
  int i;

  double x, y, z, m, r, t, V, vx, vy, vz;

  // average mass of particle 0.5
  // total mass expect 0.5*pts
  // total mass within sphere M*V_sphere/V_box = 0.5*pts*(4/3)*pi*r**3/(2R)**3
  //                                           = pts*pi r**3/12R**3 = (pts*pi/12)*r^3
  // V_circ = sqrt(GM/R) = sqrt(G*pts*pi/12r^2) = sqrt(G*pts*pi/12)/r
  
  for(i = 1; i<=n; i++) {
    m = ran1(&idum);
    x = ran1(&idum) - 0.5;
    y = ran1(&idum) - 0.5;
    z = ran1(&idum) - 0.5;
    r = sqrt(pow(x, 2.0) + pow(y, 2.0) + pow(z, 2.0));

    if(flags & DEBUG) printf("Round one generation [%d]\n", i);

    if(R<r) {
      while(R<r) {
	x = ran1(&idum) - 0.5;
	y = ran1(&idum) - 0.5;
	z = ran1(&idum) - 0.5;

	r = sqrt(pow(x, 2.0) + pow(y, 2.0) + pow(z, 2.0));
      }
    }

    V = sqrt(G*n*M_PI/12.0)/r;
    vx = sqrt(V)*y/r;
    vy = -sqrt(V)*x/r;
    vz = 0.0;
      
    if(l) {
      if(!dl) {
	printf("create_sphere_simple error: dl[] not provided\n");
	exit(1);
      }
      l[l_i(i,4)] = m;
      l[l_i(i,1)] = x;
      l[l_i(i,2)] = y;
      l[l_i(i,3)] = z;

      dl[l_i(i,4)] = 0.0;
      dl[l_i(i,1)] = vx;
      dl[l_i(i,2)] = vy;
      dl[l_i(i,3)] = vz;
    }
  }
}

void create_solar(unsigned int flags, PARTICLE pts[]) {
  // 0 sun, 1 mer, 2 ven, 3 ear, 4 mar, 5 jup, 6 sat, 7 nep, 8 ura, 9 plu
  double *r, *v, *theta, *phi, *mass; // theta measured from ecliptic (basically inc)
  r = malloc(sizeof(double)*10);
  v = malloc(sizeof(double)*10);
  theta = malloc(sizeof(double)*10);
  phi = malloc(sizeof(double)*10);
  mass = malloc(sizeof(double)*10);

  double r_u, v_u, m_u;
  r_u = 1.496E+2; // convert from 10E6 km
  v_u = 30.0; // convet from km/s
  m_u = 1.989E+6; // convert from 10E24 kg
  
  // Sun
  r[0] = 0.0;
  v[0] = 0.0;
  theta[0] = 0.0;
  phi[0] = 0.0;
  mass[0] = 1.0;
  // Mercury
  r[1] = 57.9/r_u;
  v[1] = 47.4/v_u;
  theta[1] = 7.0*M_PI/180;
  phi[1] = 0.0;
  mass[1] = 0.33/m_u;
  // Venus
  r[2] = 108.2/r_u;
  v[2] = 35.0/v_u;
  theta[2] = 3.4*M_PI/180.0;
  phi[2] = 0.0;
  mass[2] = 4.87/m_u;
  // Earth
  r[3] = 149.6/r_u;
  v[3] = 29.8/v_u;
  theta[3] = 0.0;
  phi[3] = 0.0;
  mass[3] = 5.97/m_u;
  // Mars
  r[4] = 228.0/r_u;
  v[4] = 24.1/v_u;
  theta[4] = 1.8*M_PI/180.0;
  phi[4] = 0.0;
  mass[4] = 0.642/m_u;
  // Jupiter
  r[5] = 778.5/r_u;
  v[5] = 13.1/v_u;
  theta[5] = 1.3*M_PI/180.0;
  phi[5] = 0.0;
  mass[5] = 1898.0/m_u;
  // Saturn
  r[6] = 1432.0/r_u;
  v[6] = 9.7/v_u;
  theta[6] = 2.5*M_PI/180.0;
  phi[6] = 0.0;
  mass[6] = 568.0/m_u;
  // Neptune
  r[7] = 2867.0/r_u;
  v[7] = 6.8/v_u;
  theta[7] = 0.8*M_PI/180.0;
  phi[7] = 0.0;
  mass[7] = 86.8/m_u;
  // Uranus
  r[8] = 4515.0/r_u;
  v[8] = 5.4/v_u;
  theta[8] = 1.8*M_PI/180.0;
  phi[8] = 0.0;
  mass[8] = 102.0/m_u;
  // Pluto
  r[9] = 5906.4/r_u;
  v[9] = 4.7/v_u;
  theta[9] = 17.2*M_PI/180.0;
  phi[9] = 0.0;
  mass[9] = 0.013/m_u;

  int i;
  for(i = 0; i < n; i++) {
    pts[i].mass = mass[i];
    pts[i].pos[0] = r[i]*cos(phi[i])*sin((M_PI/2.0)-theta[i]);// x
    pts[i].pos[1] = r[i]*sin(phi[i])*sin((M_PI/2.0)-theta[i]);
    pts[i].pos[2] = r[i]*cos((M_PI/2.0)-theta[i]);
    
    pts[i].vel[0] = v[i]*sin(phi[i]);
    pts[i].vel[1] = v[i]*cos(phi[i]);
    pts[i].vel[2] = 0.0;
  }
}

void create_ring(unsigned int flags, double r_inner, double width, double e_max,
		 double height, double m_avg, double m_sig, double m_orb,
		 int group, int n_tmp, PARTICLE pts[]) {
  int i;
  long idum = -1*((long) getpid());
  double a, b, e, r, m, theta, *vals, v;
  a = width/2.0;

  for(i=0;i<n_tmp;i++) {
    r = width*ran1(&idum) + r_inner;
    theta = 2.0*M_PI*ran1(&idum);
    // 'b' as in baking powder
    b = height*(2.0*ran1(&idum) - 1.0);
    e = e_max*ran1(&idum);

    // velocity
    v = sqrt((G*m_orb/r)*(1-e));
    m = m_sig*(2.0*ran1(&idum) - 1.0) + m_avg;
    vals = calloc(7, sizeof(double));

    vals[0] = m;
    vals[1] = r*cos(theta);
    vals[2] = r*sin(theta);
    vals[3] = b*sqrt(1 - sqr((r-r_inner-a)/a));
    vals[4] = -v*sin(theta);
    vals[5] = v*cos(theta);
    vals[6] = 0.0;

    pts[i] = create_point(vals, group);
  }
}

void create_satellite(unsigned int flags, double mass, double a_ax, double ecc,
		      double inc, double theta, double p_mass, PARTICLE *pt) {

}
