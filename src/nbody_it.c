#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>

#include "./lib/nr.h"
#include "./lib/nrutil.h"

#include "./lib/int_methods.h"
#include "./lib/nbody.h"

// globals, number of points and softening param
int n = 0;
double eta = ETA;

void solve_leapfrog(double ecc, double t, uint steps, double h, int out_f,
		    FILE *fp_i, FILE *fp_o, int fp_buff, char *fp_path,
		    uint flags) {
  // start
  if(flags & DEBUG) printf("nbody_it.solve_leapfrog@start\n");
  int i, j;
  double *y, *dydt, *d2ydt2;

  // 2body_init
  if(n==2) {
    if(flags & DEBUG) printf("nbody_it.solve_leapfrog->2body_init\n");
    // initializing vectors for two body
    y = dvector(1, 4*n);
    dydt = dvector(1, 4*n);
    d2ydt2 = dvector(1, 4*n);
    // returns filled in vectors for leapfrog
    create2body(ecc, flags, y, dydt, NULL);
    // init_d2ydt2
    if(flags & DEBUG) printf("nbody_it.solve_leapfrog->2body_init@init_d2ydt2\n");
  } else if(flags & READ_FILE) {
    if(flags & DEBUG) printf("nbody_it.solve_leapfrog->read_file\n");
    // file reading
    // set particles to max reading from file
    n = N_MAX;
    int k;
    // declare the data matrix
    double **d;    
    d = dmatrix(1, 7, 1, n);
    if(flags & DEBUG) printf("nbody_it.solve_leapfrog->read_file@vectors_init\n");
    // functions reads file and outputs to d and gives number of points read in n
    state_from_file(fp_i, &n, flags, d, NULL, NULL, NULL);
    if(flags & DEBUG) printf("nbody_it.solve_leapfrog->read_file@file_read\n");
    // now, initialize vectors
    y = dvector(1, 4*n);
    dydt = dvector(1, 4*n);
    d2ydt2 = dvector(1, 4*n);
    // loop for unzipping d into the leapfrog vectors
    for(k = 1; k <= n; k++) {
      y[l_i(k,4)] = d[1][k];
      dydt[l_i(k,4)] = 0.0;
      y[l_i(k,1)] = d[2][k];
      dydt[l_i(k,1)] = d[5][k];
      y[l_i(k,2)] = d[3][k];
      dydt[l_i(k,2)] = d[6][k];
      y[l_i(k,3)] = d[4][k];
      dydt[l_i(k,3)] = d[7][k];
    }
    free_dmatrix(d, 1, 7, 1, N_MAX);
  } else if(flags & GEN_SPHERE) {
    y = dvector(1, 4*n);
    dydt = dvector(1, 4*n);
    d2ydt2 = dvector(1, 4*n);
    create_sphere_simple(n, flags, y, dydt, NULL);
  } else if(flags & GEN_DISK) {
    // unimplemented
  } else {
    printf("nbody_it: nothing to do...\n");
    exit(1);
  }

  if(flags & DEBUG) {
    printf("initial conditions: \n");
    for(j = 1; j <= n; j++) {
      printf("%.6g %.6g %.6g %.6g %.6g %.6g %.6g\n", y[l_i(j,4)], y[l_i(j,1)],
	     y[l_i(j,2)], y[l_i(j,3)], dydt[l_i(j,1)], dydt[l_i(j,2)],
	     dydt[l_i(j,3)]);
    }
  }
  
  // stepsloop
  for(i = 1; i <= steps; i++) {
    if(flags & DEBUG) printf("nbody_it.solve_leapfrog->stepsloop [%d]\n", i);
    leap_std(y, dydt, d2ydt2, 4*n, t, h, y, l_derivs);
    // leaped
    if(flags & DEBUG) printf("nbody_it.solve_leapfrog->stepsloop@leaped\n");
    t += h;

    // outflow, triggers on multiples of the output frequency
    if(!(i%out_f)) {
      if(flags & PRINT_DIR) {
	// prints multiple files to a directory
	if(flags & DEBUG) printf("nbody_it.solve_leapfrog->stepsloop->outflow->\
PRINT_DIR [%d]\n", i);
	char fp_path_temp[BUFF];
	strcpy(fp_path_temp,fp_path);
	// create the name for this step's file
	snprintf(fp_path_temp+fp_buff, BUFF, "data_s_%d.csv", i);
	if(flags & DEBUG) printf("nbody_it.solve_leapfrog->stepsloop->outflow->\
PRINT_DIR@snprintf\n");
	FILE *fp_o2;
	fp_o2 = fopen(fp_path_temp, "w");
	if(!fp_o2) {
	  printf("nbody_it: error in writing to %s", fp_path);
	  exit(1);
	}
	// now print each particle to a line in the file
	if(flags & DEBUG) printf("nbody_it.solve_leapfrog->stepsloop->outflow->\
PRINT_DIR@path_opened\n");
	for(j = 1; j <= n; j++) {
	  fprintf(fp_o2, "%.6g %.6g %.6g %.6g %.6g %.6g %.6g\n", y[l_i(j,4)],
		  y[l_i(j,1)], y[l_i(j,2)], y[l_i(j,3)], dydt[l_i(j,1)],
		  dydt[l_i(j,2)], dydt[l_i(j,3)]);
	}
	if(flags & DEBUG) printf("nbody_it.solve_leapfrog->stepsloop->outflow->\
PRINT_DIR@print_finish\n");
	fclose(fp_o2);
      }
      else if(flags & PRINT_CONS) {
	// output directly to console
	printf("step: %d, time: %.3f\n", i, t);
	for(j = 1; j <= n; j++) {
	  printf("%.6g %.6g %.6g %.6g %.6g %.6g %.6g\n", y[l_i(j,4)],
		 y[l_i(j,1)], y[l_i(j,2)], y[l_i(j,3)], dydt[l_i(j,1)],
		 dydt[l_i(j,2)], dydt[l_i(j,3)]);
	}
	printf("\n");
      }
    }
  }
}

void solve_leapfrog_pt(double ecc, double t, uint steps, double h, int out_f,
		       FILE *fp_i, FILE *fp_o, int fp_buff, char *fp_path,
		       uint flags) {
  // start
  if(flags & DEBUG) printf("nbody_it.solve_leapfrog_pt@start\n");
  int i, j;
  double *y, *dydt, *d2ydt2;
  PARTICLE *pts;

  // 2body_init
  if(flags & GEN_2BODY) {
    if(flags & DEBUG) printf("nbody_it.solve_leapfrog_pt->2body_init\n");
    // initializing vectors for two body
    y = dvector(1, 4*n);
    dydt = dvector(1, 4*n);
    d2ydt2 = dvector(1, 4*n);
    // returns filled in vectors for leapfrog
    create2body(ecc, flags, y, dydt, NULL);
    pts = (PARTICLE *) malloc(sizeof(PARTICLE)*n);
    l_to_pt(y, dydt, n, pts);
    free_dvector(y, 1, 4*n);
    free_dvector(dydt, 1, 4*n);
    free_dvector(d2ydt2, 1, 4*n);
    // init_d2ydt2
    if(flags & DEBUG) printf("nbody_it.solve_leapfrog_pt->2body_init@init_d2ydt2\
\n");
  } else if(flags & READ_FILE) {
    if(flags & DEBUG) printf("nbody_it.solve_leapfrog_pt->read_file\n");
    // file reading
    // set particles to max reading from file
    n = N_MAX;
    int k;
    // declare the data matrix
    double **d;    
    d = dmatrix(1, 7, 1, n);
    if(flags & DEBUG) printf("nbody_it.solve_leapfrog_pt->read_file@vectors_init\
\n");
    // functions reads file and outputs to d and gives number of points read in n
    state_from_file(fp_i, &n, flags, d, NULL, NULL, NULL);
    if(flags & DEBUG) printf("nbody_it.solve_leapfrog_pt->read_file@file_read\n");
    // now, initialize vectors
    y = dvector(1, 4*n);
    dydt = dvector(1, 4*n);
    d2ydt2 = dvector(1, 4*n);
    // loop for unzipping d into the leapfrog vectors
    for(k = 1; k <= n; k++) {
      y[l_i(k,4)] = d[1][k];
      dydt[l_i(k,4)] = 0.0;
      y[l_i(k,1)] = d[2][k];
      dydt[l_i(k,1)] = d[5][k];
      y[l_i(k,2)] = d[3][k];
      dydt[l_i(k,2)] = d[6][k];
      y[l_i(k,3)] = d[4][k];
      dydt[l_i(k,3)] = d[7][k];
    }
    free_dmatrix(d, 1, 7, 1, N_MAX);
    pts = (PARTICLE *) malloc(sizeof(PARTICLE)*n);
    l_to_pt(y, dydt, n, pts);
    free_dvector(y, 1, 4*n);
    free_dvector(dydt, 1, 4*n);
    free_dvector(d2ydt2, 1, 4*n);
  } else if(flags & GEN_SPHERE) {
    y = dvector(1, 4*n);
    dydt = dvector(1, 4*n);
    d2ydt2 = dvector(1, 4*n);
    create_sphere_simple(n, flags, y, dydt, NULL);
    pts = (PARTICLE *) malloc(sizeof(PARTICLE)*n);
    l_to_pt(y, dydt, n, pts);
    free_dvector(y, 1, 4*n);
    free_dvector(dydt, 1, 4*n);
    free_dvector(d2ydt2, 1, 4*n);
  } else if(flags & GEN_SOLAR) {
    n = 7;
    pts = (PARTICLE *) malloc(sizeof(PARTICLE)*n);
    create_solar(flags, pts);
  } else {
    printf("nbody_it: nothing to do...\n");
    exit(1);
  }

  if(flags & DEBUG) {
    printf("initial conditions: \n");
    for(j = 0; j < n; j++) {
      printf("%.6g %.6g %.6g %.6g %.6g %.6g %.6g\n", pts[j].mass, pts[j].pos[0],
	     pts[j].pos[1], pts[j].pos[2], pts[j].vel[0], pts[j].vel[1],
	     pts[j].vel[2]);
    }
  }

  derivs_pt(t, pts);
  // stepsloop
  for(i = 0; i < steps; i++) {
    if(flags & DEBUG) printf("nbody_it.solve_leapfrog_pt->stepsloop [%d]\n", i+1);
    leap_pt(pts, n, t, h, derivs_pt);
    // leaped
    if(flags & DEBUG) printf("nbody_it.solve_leapfrog_pt->stepsloop@leaped\n");
    t += h;

    // outflow, triggers on multiples of the output frequency
    if(!((i+1)%out_f)) {
      if(flags & PRINT_DIR) {
	// prints multiple files to a directory
	if(flags & DEBUG) printf("nbody_it.solve_leapfrog_pt->stepsloop->\
outflow->PRINT_DIR [%d]\n", i+1);
	char fp_path_temp[BUFF];
	strcpy(fp_path_temp,fp_path);
	// create the name for this step's file
	snprintf(fp_path_temp+fp_buff, BUFF, "data_s_%d.csv", i);
	if(flags & DEBUG) printf("nbody_it.solve_leapfrog_pt->stepsloop->\
outflow->PRINT_DIR@snprintf\n");
	FILE *fp_o2;
	fp_o2 = fopen(fp_path_temp, "w");
	if(!fp_o2) {
	  printf("nbody_it: error in writing to %s", fp_path);
	  exit(1);
	}
	// now print each particle to a line in the file
	if(flags & DEBUG) printf("nbody_it.solve_leapfrog_pt->stepsloop->\
outflow->PRINT_DIR@path_opened\n");
	if(!(flags & VERBOSE)) {
	  for(j = 0; j < n; j++) {
	    fprintf(fp_o2, "%.6g %.6g %.6g %.6g %.6g %.6g %.6g\n", pts[j].mass,
		    pts[j].pos[0], pts[j].pos[1], pts[j].pos[2], pts[j].vel[0],
		    pts[j].vel[1], pts[j].vel[2]);
	  }
	} else {
	  for(j = 0; j < n; j++) {
	    fprintf(fp_o2, "%.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g\n",
		    pts[j].mass, pts[j].pos[0], pts[j].pos[1], pts[j].pos[2],
		    pts[j].vel[0], pts[j].vel[1], pts[j].vel[2], pts[j].acc[0],
		    pts[j].acc[1], pts[j].acc[2]);
	  }
	}
	if(flags & DEBUG) printf("nbody_it.solve_leapfrog_pt->stepsloop->\
outflow->PRINT_DIR@print_finish\n");
	fclose(fp_o2);
      }
      if(flags & PRINT_CONS) {
	// output directly to console
	printf("step: %d, time: %.3f\n", i+1, t);
	for(j = 0; j < n; j++) {
	  printf("%.6g %.6g %.6g %.6g %.6g %.6g %.6g\n", pts[j].mass,
		 pts[j].pos[0], pts[j].pos[1], pts[j].pos[2], pts[j].vel[0],
		 pts[j].vel[1], pts[j].vel[2]);
	}
	printf("\n");
      }
    }
  }
}

int main(int argc, char *argv[]) {
  // initializing options variables, defaults are:
  // output frequency: 10
  // steps: 100000
  // stepsize:
  // eccentricity: 0.0
  // generation: 2body
  char opt, *read_name, *print_name, *opts_name, *opts_file;
  uint flags = DEF_OP;
  int out_f = FREQ_OUT;
  uint steps = STEPS;
  double h = H;
  double ecc = 0.0;
  n = 2;
  flags |= MASS_NORM;
  while( (opt = getopt(argc, argv, "a:bde:f:F:g:h:Hm:n:o:O:p:s:vx")) != -1 ) {
    switch(opt) {
    case 'a':
      // softening parameter
      eta = (double) atof(optarg);
      break;
    case 'b':
      // two body mode
      n = 2;
      flags |= MASS_NORM;
      break;
    case 'd':
      // debug flag
      flags |= DEBUG;
      break;
    case 'e':
      // eccentricity for two body
      ecc = (double) atof(optarg);
      if( (ecc < 0.0) || (ecc > 1.0) ) {
	printf("%s: argument error -%c: eccentricity must be between 0.0 and 1.0\
\n", argv[0], opt);
	return 1;
      }
      break;
    case 'f':
      // file input mode
      flags |= READ_FILE;
      flags &= ~GROUP_MODE;
      read_name = optarg;
      break;
    case 'g':
      flags |= GEN_SPHERE;
      n = atoi(optarg);
      break;
    case 'F':
      // options file input
      flags |= OPTS_FILE;
      opts_name = optarg;
      break;
    case 'h':
      // stepsize input
      h = (double) atof(optarg);
      if(h < 0.0) {
	printf("%s: argument warning -%c: stepsize negative\n", argv[0], opt);
      }
      break;
    case 'n':
      // output frequency in steps per print
      out_f = atoi(optarg);
      if(out_f <= 0) {
	printf("%s: argument error -%c: must have output greater than 0\n",
	       argv[0], opt);
	return 1;
      }
      break;
    case 'o':
      // output directory mode
      flags |= PRINT_DIR;
      print_name = optarg;
      if(print_name[strlen(print_name)-1] != '/') {
	printf("%s: argument error -%c: must be a directory ending in '/'\n",
	       argv[0], opt);
	return 0;
      }
      break;
    case 'p':
      // print options
      flags |= PRINT_OPTS;
      opts_file = optarg;
      break;
    case 's':
      // total steps
      steps = atoi(optarg);
      if(steps < 0) {
	printf("%s: argument error -%c: steps cannot be negative\n", argv[0],
	       opt);
	return 1;
      }
      break;
    case 'v':
      // verbosity
      flags |= VERBOSE;
      break;
    case 'x':
      // prematurely abort the program
      flags |= EXIT_RUN;
      break;
    case '?':
      break;
    default:
      break;
    }
  }

  // file fmt: m x y z vx vy vz
  FILE *fp_i, *fp_o, *fp_opts;
  char fp_path[BUFF];
  int fp_buff;
  // read options from a file (not yet finished)
  if(flags & OPTS_FILE) {
    fp_opts = fopen(opts_name, "r");
    if(!fp_opts) {
      printf("%s: error: options input file \"%s\" does not exit\n", argv[0],
	     opts_name);
      return 1;
    }
    opts_from_file(fp_opts, &flags, &eta, &ecc, &read_name, &h, &out_f,
		   &print_name, print_name, &steps);
    fclose(fp_opts);
  }

  // open the input file
  if(flags & READ_FILE) {
    if(flags & DEBUG) printf("Detected read file...\n");
    fp_i = fopen(read_name, "r");
    if(!fp_i) {
      printf("%s: error: input file \"%s\" does not exist\n", argv[0], read_name);
      return 1;
    }
  }
  // prepare output file/directory
  if(flags & PRINT_DIR) fp_buff = snprintf(fp_path, BUFF, "%s", print_name);
  if(flags & PRINT_OPTS) print_opts(opts_file, flags, eta, ecc, read_name,
				    h, out_f, print_name, print_name, steps);
  
  
  double t;
  t = 0.0;

  if(flags & EXIT_RUN) return 0;
  if(flags & DEBUG) printf("Console mode: %d\n", flags & PRINT_CONS);
  
  // mode_exec
  if(flags & DEBUG) printf("nbody_it.main@mode_exec\n"); 
  if(!(flags & MODE_BH)) solve_leapfrog(ecc, t, steps, h, out_f, fp_i, fp_o,
					fp_buff, fp_path, flags);
  if(flags & MODE_BH) solve_leapfrog_pt(ecc, t, steps, h, out_f, fp_i, fp_o,
					fp_buff, fp_path, flags);
  
  if(flags & READ_FILE) fclose(fp_i);
  
  return 0;
}
