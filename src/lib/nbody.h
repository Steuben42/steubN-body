#ifndef NBODY
#define NBODY

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nr.h"
#include "nrutil.h"



/*************************
     DEFINE CONSTANTS
*************************/
// utility constants
#define BUFF   127   // general character buffer
#define F_BUFF 255   // read buffer
#define N_MAX  10000 // max points from file

// calculation constants
#define G 1.0 // mass in M_sol, distance in AU, time in yr/2pi, speeds in 30 km/s
#define G_GAL 0.0 // calculate using galaxy units

/**** flags *****/
// runtime modifiers (4)
#define DEBUG      (1 << 0)
#define MODE_BH    (1 << 1) // Barnes-Hutt mode
#define EXIT_RUN   (1 << 2) // prematurely aborts the program
// output options (4)
#define PRINT_DIR  (1 << 4) // implies directory output mode
#define PRINT_CONS (1 << 5) // console printing
#define PRINT_OPTS (1 << 6) // print options to a file
#define VERBOSE    (1 << 7) // extra output
// generator options
// gen opts - mass (6)
#define MASS_NORM  (1 << 8) // unit mass
#define MASS_UNI   (1 << 9) // uniform distribution of mass
#define MASS_RAY   (1 << 10) // rayleigh distro mass
#define GROUP_MASS (MASS_NORM | MASS_UNI | MASS_RAY)
// gen opts - distro (8)
#define GEN_SPHERE (1 << 14) // generating with sphere
#define GEN_DISK   (1 << 15) // generating with disk
#define GEN_2BODY  (1 << 16) // 2 body generation
#define GEN_SOLAR  (1 << 17) // generate solar setup
#define GROUP_MODE (GEN_2BODY | GEN_SPHERE | GEN_DISK | GEN_SOLAR)
// gen opts - scale (4)
#define SCALE_SOL  (1 << 22)
#define SCALE_GAL  (1 << 23)
// input options (2)
#define READ_FILE  (1 << 26) // implies input mode
#define OPTS_FILE  (1 << 27) // options input file
//...
#define FLAG_31    (1 << 31)

// defaults
#define ETA      0.00 // softening
#define FREQ_OUT 10   // printout frequency (steps/print)
#define H        0.05 // stepsize
#define STEPSIZE 0.05
#define STEPS    1000 
// default flags
#define DEF_OP (0 | MASS_NORM | MODE_BH | GEN_2BODY | PRINT_CONS)


/*************************
   EXTERNS AND STRUCT
*************************/
// extern globals
extern int n;      // particle number
extern double eta; // softening parameter

// Particle struct
typedef struct Particle_t PARTICLE;
struct Particle_t {
  int id;
  int group;
  double pos[3];
  double vel[3];
  double mass;
  double acc[3];
};

// Polar struct - same format as a normal particle, separate for sanity
typedef struct Polar_t POLAR;
struct Polar_t {
  int id;
  int group;
  double pos[3]; // r, phi, theta
  double vel[3]; // v, v_phi, v_theta
  double mass;
  double acc[3]; // a_r, a_phi, a_theta
};

// Node struct
typedef struct Node_t NODE;
struct Node_t {
  short order; // 0 is root
  int id; // may not be necessary
  int branch; // 0 leaf, 1 branch
  double size;
  double pos[3]; // center of node
  double mass;
  double com[3];
  struct Node_t *c_node[8]; // 0th: +++, 1st: ++-; 2nd: +-+; 3rd: +--; 4th: -++;
  // 5th: -+-; 6th: --+; 7th: ---;
  struct Node_t *p_node;
  struct Particle_t *c_pt[8];
  short pts;
};

/**********************
     FILE HEADERS
**********************/
// util_nbody.c
  // indexing
int l_i(int i, int pos);
  // conversion
void pt_to_l(PARTICLE *pts, int n, double *l, double *dl);
void l_to_pt(double *l, double *dl, int n, PARTICLE *pts);
void pol_to_pt(POLAR *pols, int m, PARTICLE *pts);
void pt_to_pol(PARTICLE *pts, int m, POLAR *pols);
  // string manipulation
char *rtrim(char *str);
char *ltrim(char *str);
char *trim(char *str);
  // reading functions
int points_from_file(FILE *fp, unsigned int flags);
void state_from_file(FILE *fp, int *n, unsigned int flags, double *d[],
		     double l[], double dl[], double u[]);
void opts_from_file(FILE *fp, unsigned int *flags, double *soft, double *ecc,
		    char **f_in, double *h, int *out_f, char **d_out,
		    char *f_out, int *steps);
PARTICLE *obj_from_file(FILE *fp, unsigned int *flags);
void raise_error(unsigned short type, char *parent, char *param, char *val, int ln);

// print_nbody.c
void print_vector(double *v, int n);
void print_opts(char *file, unsigned int flags, double soft, double ecc, char *f_in,
		double h, int out_f, char *d_out, char *f_out, int steps);
char *pt_to_str(PARTICLE pt, unsigned int flags, char *str);
void print_frame_to_file(FILE *fp, unsigned int flags, PARTICLE *pts, int n);
//void print_obj();

// derivs_nbody.c
void l_derivs(double t, double y[], double d2ydt2[]);
void f_derivs(double t, double y[], double dydt[]);
double mag(double a[], int c);
double mag_diff(double a[], double b[], int c);
double sqr(double a);
void derivs_pt(double t, PARTICLE pt[]);

// generator_nbody.c
void create2body(double ecc, unsigned int flags, double l[], double dl[],
		 double u[]);
void create_sphere_simple(int n, unsigned int flags, double l[], double dl[],
			  double u[]);
void create_solar(unsigned int flags, PARTICLE pts[]);
void create_ring(unsigned int flags, double r_inner, double width, double e_max,
		 double height, double m_avg, double m_sig, double m_orb,
		 int group, int n_tmp, PARTICLE pts[]);

// bhmode_nbody.c
NODE *create_root_node(PARTICLE *pts, int n);
NODE *create_node(NODE *p_node, int oct);
PARTICLE create_point(double *vals, int group);
POLAR create_polar(double *vals, int group);
void put_pt_in_tree(PARTICLE *pt, NODE *node);
void calc_moment(NODE *node);
NODE *build_tree(PARTICLE *pts);
void add_force(PARTICLE *pt, NODE *node);
void free_tree(NODE *node);

#endif
