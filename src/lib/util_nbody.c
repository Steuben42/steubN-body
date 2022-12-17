#include "nbody.h"
#include <ctype.h>

/*********************
    FILE CONSTANTS
*********************/
// Error values
#define VALUE_ERROR  0
#define NAME_ERROR   1
#define NAME_WARNING 2
#define EMPTY_ERROR  3
#define FMT_ERROR    4
// Context values
#define CTX_GLOBAL 0
#define CTX_UNSPEC 1
#define CTX_SIM    2
#define CTX_IO     3
#define CTX_GEN    4
#define CTX_GENOBJ 5
// Object context values
#define CTX_2BODY      2
#define CTX_SPHERE     3
#define CTX_GROUP      4
#define CTX_POINT      5
#define CTX_RING       6
#define CTX_SATELLITES 7
// Context flags
#define POLAR_FLAG (1 << 0)
#define CART_FLAG  (1 << 1)


/*********************
  INDEXING FUNCTIONS
*********************/

// helper index functions (since particles in y[] are indexed x, y, z, m)
int l_i(int i, int pos) {
  return 4*i + pos - 4;
}


/***********************
  CONVERSION FUNCTIONS
***********************/

// convert a list of particle to the l/dl format
void pt_to_l(PARTICLE *pts, int n, double *l, double *dl) {
  int i, k;
  // pts is off by one index (starts at 0)
  for(i=1; i<=n; i++) {
    for(k=0; k<3; k++)
      l[l_i(i, k+1)] = pts[i-1].pos[k]; // Position conversion
    l[l_i(i, 4)] = pts[i-1].mass;

    for(k=0; k<3; k++)
      dl[l_i(i, k+1)] = pts[i-1].vel[k]; // Velocity conversion
    dl[l_i(i, 4)] = 0.0;
  }
}

// ditto vice versa
void l_to_pt(double *l, double *dl, int n, PARTICLE *pts) {
  int i, k;
  static int id = 0;
  for(i=1; i<=n; i++) {
    for(k=0; k<3; k++)
      pts[i-1].pos[k] = l[l_i(i, k+1)];
    pts[i-1].mass = l[l_i(i, 4)];

    for(k=0; k<3; k++)
      pts[i-1].vel[k] = dl[l_i(i, k+1)];
    pts[i-1].id = id++;
  }
}

// Polars to points function
void pol_to_pt(POLAR *pols, int m, PARTICLE *pts) {
  int i;
  double r, ph, th, v, ph_v, th_v;
  for(i=0; i<m; i++) {
    r = pols[i].pos[0];
    th = pols[i].pos[1];
    ph = pols[i].pos[2];
    v = pols[i].vel[0];
    th_v = pols[i].vel[1];
    ph_v = pols[i].vel[2];

    pts[i].id = pols[i].id;
    pts[i].mass = pols[i].mass;
    
    pts[i].pos[0] = r*cos(th)*sin(ph);
    pts[i].pos[1] = r*sin(th)*sin(ph);
    pts[i].pos[2] = r*cos(ph);

    pts[i].vel[0] = v*cos(th_v)*sin(ph_v);
    pts[i].vel[1] = v*sin(th_v)*sin(ph_v);
    pts[i].vel[2] = v*cos(ph_v);

    pts[i].group = pols[i].group;
  }
}

void pt_to_pol(PARTICLE *pts, int m, POLAR *pols) {
  int i;
  double x, y, z, v_x, v_y, v_z;
  for(i=0; i<m; i++) {
    x = pts[i].pos[0];
    y = pts[i].pos[1];
    z = pts[i].pos[2];
    v_x = pts[i].vel[0];
    v_y = pts[i].vel[1];
    v_z = pts[i].vel[2];

    pols[i].mass = pts[i].mass;
    pols[i].id = pts[i].id;

    pols[i].pos[0] = sqrt(abs(pow(x, 2.0) + pow(y, 2.0) + pow(z, 2.0)));
    pols[i].pos[1] = atan2(y, x);
    pols[i].pos[2] = acos(z/pols[i].pos[0]);
    if(pols[i].pos[0]==0.0) pols[i].pos[2] = 0.0;

    pols[i].vel[0] = sqrt(abs(pow(v_x, 2.0) + pow(v_y, 2.0) + pow(v_z, 2.0)));
    pols[i].vel[1] = atan2(v_y, v_x);
    pols[i].vel[2] = acos(v_z/pols[i].vel[0]);
    if(pols[i].vel[0]==0.0) pols[i].vel[2] = 0.0;

    pols[i].group = pts[i].group;
  }
}

/**********************
  STRING MANIPULATION
**********************/

// right trim function
char *rtrim(char *str) {
  char *back = str + strlen(str);
  while(isspace(*--back));
  *(back+1) = '\0';
  return str;
}

// left trim function
char *ltrim(char *str) {
  while(isspace(*str)) str++;
  return str;
}

// trim function
char *trim(char *str) {
  return rtrim(ltrim(str));
}

/********************
  READING FUNCTIONS
********************/

// reader function for line counting
int points_from_file(FILE *fp, unsigned int flags) {
  if(flags & DEBUG) printf("nbody_util.points_from_file\n");
  int k;
  rewind(fp);
  if(flags & DEBUG) printf("nbody_util.points_from_file@beginning_scan\n");
  // k will equal number of points
  for(k = 0; fscanf(fp, "%*f %*f %*f %*f %*f %*f %*f") != EOF; k++);
  if(!k) {
    printf("points_from_file error: file could not be read properly\n");
    exit(1);
  }
  if(flags & DEBUG) printf("nbody_util.points_from_file@end_scan\n");
  return k;
}

// reader function
void state_from_file(FILE *fp, int *n, unsigned int flags, double *d[], double l[],
		     double dl[], double u[]) {
  int k;
  rewind(fp);
  // reading when leapfrog vectors are provided
  if(l) {
    for(k = 1; k <= (*n); k++) {
      fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf\n", &l[l_i(k,4)], &l[l_i(k,1)],
	     &l[l_i(k,2)], &l[l_i(k,3)], &dl[l_i(k,1)], &dl[l_i(k,2)],
	     &dl[l_i(k,3)]);
      if(flags & DEBUG) printf("nbody_util.state_from_file->l_loop [%d]\n", k);
      if(flags & DEBUG) printf("point read: x y z - %lf %lf %lf\n", l[l_i(k,1)],
			       l[l_i(k,2)], l[l_i(k,3)]);
      dl[l_i(k,4)] = 0.0;
    }
  } else if (d) {
    // reading when a data matrix is provided
    for(k = 1; fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf\n", &d[1][k], &d[2][k], &d[3][k],
		      &d[4][k], &d[5][k], &d[6][k], &d[7][k]) != EOF; k++);
    (*n) = k-1;
  } else {
    printf("state_from_file error: no arrays initialized\n");
  }
}

// options reading file function
void opts_from_file(FILE *fp, unsigned int *flags, double *soft, double *ecc,
		    char **f_in, double *h, int *out_f, char **d_out,
		    char *f_out, int *steps) {
  rewind(fp);
  char pt_l_buff[BUFF], l_buff[BUFF];
  char *tkn;
  int k = 0;
  int ctx = 0; // context; 0 is global, 1 is unspec local, 2 is simulator, 3 is I/O,
  // 4 is generator, 5 is generator object
  const char del[2] = ":";

  if(*flags & DEBUG) printf("Beginning file loop...\n");
  
  while(fgets(pt_l_buff, BUFF, fp)) {
    k++;

    // error if still in a context with an @ group declaration
    if( ctx && (pt_l_buff[0] == '@') ) {
      raise_error(FMT_ERROR, NULL, NULL, "missing closing bracket\n", k);
      continue;
    }
    // close context if an end of group detected
    if( ctx && (pt_l_buff[0]=='}') ) {
      ctx = 0;
      continue;
    }

    if(ctx) strcpy(l_buff, trim(pt_l_buff));
    else strcpy(l_buff, trim(pt_l_buff));
    
    switch(ctx) {
    case CTX_GLOBAL:
      // global context, skip non-@ led lines
      if(l_buff[0] != '@') continue;
      
      // global context, processing a group
      if(*flags & DEBUG) printf("Found group key...\n");
      if(l_buff[strlen(l_buff)-1] == '{') {
	ctx++;
	char parent[32];
	sscanf(l_buff, "@%[^{]", parent);
	if(*flags & DEBUG) printf("Read parent: %s\n", parent);

	if(!strcmp(parent, "SIM")) ctx = CTX_SIM;
	if(!strcmp(parent, "IO")) ctx = CTX_IO;
	if(!strcmp(parent, "GEN")) ctx = CTX_GEN;
	if(!strcmp(parent, "GAL")) ctx = CTX_GENOBJ;
	if(ctx==CTX_UNSPEC) raise_error(NAME_WARNING, parent, NULL, NULL, k);
	continue;
      }
      break;
    case CTX_UNSPEC:
      
      break;
    case CTX_SIM:
      tkn = strtok(l_buff, del);
      if(*flags & DEBUG) printf("tkn: %s\n", tkn);
      // SOFTENING
      if(!strcmp(tkn, "softening")) {
	*soft = (double) atof(strtok(NULL, del));
	if(*flags & DEBUG) printf("Set value of softening to %f\n", *soft);
	continue;
      }
      // STEPSIZE
      if(!strcmp(tkn, "stepsize")) {
	*h = (double) atof(strtok(NULL, del));
	if(*flags & DEBUG) printf("Set value of stepsize to %f\n", *h);
	continue;
      }
      // STEPS
      if(!strcmp(tkn, "steps")) {
	*steps = (int) atoi(strtok(NULL, del));
	if(*flags & DEBUG) printf("Set value of steps to %d\n", *steps);
	continue;
      }
      // BARNES-HUT MODE
      if(!strcmp(tkn, "bh_mode")) {
	if(atoi(strtok(NULL, del))) *flags |= MODE_BH;
	if((*flags & DEBUG) && (*flags & READ_FILE)) printf("Set to Barnes-Hutt mode\n");
	continue;
      }
      raise_error(NAME_WARNING, "SIM", tkn, NULL, k);
      break;
    case CTX_IO:
      tkn = strtok(l_buff, del);
      // FILE INPUT
      if(!strcmp(tkn, "file_in")) {
	char *in_tmp;
	in_tmp = strtok(NULL, del);
        *f_in = malloc(strlen(in_tmp)+1);
	strcpy(*f_in, in_tmp);
	if(strcmp(*f_in, "none")) *flags ^= (READ_FILE | (DEF_OP & GROUP_MODE));
	if((*flags & DEBUG) && (*flags & READ_FILE)) printf("Set value of f_in to '%s'\n",
						f_in);
	continue;
      }
      // DIR OUTPUT
      if(!strcmp(tkn, "dir_out")) {
	char *out_tmp;
	out_tmp = strtok(NULL, del);
	*d_out = malloc(strlen(out_tmp)+1);
	strcpy(*d_out, out_tmp);
	if(strcmp(*d_out, "none")) *flags |= PRINT_DIR;
	if((*flags & DEBUG) && (*flags & PRINT_DIR)) printf("Set value of d_out to '%s'\n",
						*d_out);
	continue;
      }
      // CONSOLE OUTPUT
      if(!strcmp(tkn, "console")) {
	if(!atoi(strtok(NULL, del))) *flags ^= PRINT_CONS;
	if((*flags & DEBUG) && !(*flags & PRINT_CONS)) printf("Set to not print to console\n");
	continue;
      }
      // OUTPUT FREQUENCY
      if(!strcmp(tkn, "out_freq")) {
	*out_f = atoi(strtok(NULL, del));
	if(*flags & DEBUG) printf("Set value of out_f to %d\n", *out_f);
	continue;
      }
      raise_error(NAME_WARNING, "IO", tkn, NULL, k);
      break;
    case CTX_GEN:
      tkn = strtok(l_buff, del);
      // GENERATION MODE
      if(!strcmp(tkn, "mode")) {
	tkn = strtok(NULL, del);
	if(!strcmp(tkn, "2body")) *flags ^= (GEN_2BODY ^ (*flags & GROUP_MODE));
	else if(!strcmp(tkn, "sphere")) *flags ^= (GEN_SPHERE ^ (*flags & GROUP_MODE));
	else if(!strcmp(tkn, "solar")) *flags ^= (GEN_SOLAR ^ (*flags & GROUP_MODE));
	else raise_error(VALUE_ERROR, "GEN", "mode", tkn, k);
	continue;
      }
      // OBJECT NUMBER
      if(!strcmp(tkn, "objs")) {
	continue;
      }
      // MASS DISTRO
      if(!strcmp(tkn, "mass")) {
	tkn = strtok(NULL, del);
	if(!strcmp(tkn, "uniform")) *flags ^= (MASS_NORM ^ (*flags & GROUP_MASS));
	else if(!strcmp(tkn, "rayleigh")) *flags ^= (MASS_RAY ^ (*flags & GROUP_MASS));
	else if(!strcmp(tkn, "uniform random")) *flags ^= (MASS_UNI ^ (*flags & GROUP_MASS));
	else raise_error(VALUE_ERROR, "GEN", "mass", tkn, k); 
	continue;
      }
      raise_error(NAME_WARNING, "GEN", tkn, NULL, k);
      break;
    case CTX_GENOBJ:
      tkn = strtok(l_buff, del);
      
      raise_error(NAME_WARNING, "GENOBJ", tkn, NULL, k);
      break;
    }
  }
}

// objects from file function
PARTICLE *obj_from_file(FILE *fp, unsigned int *flags) {
  rewind(fp);
  char pt_l_buff[BUFF], l_buff[BUFF];
  char *tkn, *val;
  int i = 0, ptidx = 0, group = 0, k = 0, ctx = 0, n_tmp = 1;
  double r_inner=0.0, width=0.0, e_max=0.0, height=0.0, m_avg=0.0, m_sig=0.0,
    m_orb=0.0;
  unsigned int ctx_flag = 0;
  PARTICLE *pts;
  POLAR *pols;
  const char del[2] = ":";
  const char del_list[2] = ",";
  n = 0;

  if(*flags & DEBUG) printf("Beginning objects file loop...\n");

  while(fgets(pt_l_buff, BUFF, fp)) {
    k++;
    ctx_flag = 0;

    // unclosed group
    if( ctx && (pt_l_buff[0] == '@') ) {
      raise_error(FMT_ERROR, NULL, NULL, "missing closing bracket\n", k);
      continue;
    }
    // close context if an end of group detected
    if( ctx && (pt_l_buff[0]=='}') ) {
      ctx = CTX_GLOBAL;
      group++;
      i = ptidx;
      continue;
    }

    if(ctx) strcpy(l_buff, trim(pt_l_buff));
    else strcpy(l_buff, trim(pt_l_buff));
    // comments
    if( (l_buff[0] == '-') && (l_buff[1] == '-') ) continue;

    switch(ctx) {
    case CTX_GLOBAL:
      if(l_buff[0] != '@') continue;

      if(*flags & DEBUG) printf("Found group key...\n");
      if(l_buff[strlen(l_buff)-1] == '{') {
	ctx++;
	char parent[BUFF];
	int scanres = 0;
	scanres = sscanf(l_buff, "@%[^{(](%d)", parent, &n_tmp);
	if(*flags & DEBUG) printf("Read parent: %s\n", parent);
	if(scanres==1)
	  n_tmp = 1;
	n += n_tmp;
	if(*flags & DEBUG) printf("Incremented n [%d after] by %d\n", n, n_tmp);
        
	if(n==n_tmp) pts = calloc(n, sizeof(PARTICLE));
	else pts = realloc(pts, sizeof(PARTICLE)*n);

	if(!strcmp(parent, "2BODY")) ctx = CTX_2BODY;
	if(!strcmp(parent, "SPHERE")) ctx = CTX_SPHERE;
	if(!strcmp(parent, "GROUP")) ctx = CTX_GROUP;
	if(!strcmp(parent, "POINT")) ctx = CTX_POINT;
	if(!strcmp(parent, "RING")) ctx = CTX_RING;
	if(!strcmp(parent, "SATELLITES")) ctx = CTX_SATELLITES;
	if(ctx==CTX_UNSPEC) raise_error(NAME_WARNING, parent, NULL, NULL, k);
	continue;
      }
      break;
    case CTX_UNSPEC:
      break;
    case CTX_2BODY:
      break;
    case CTX_SPHERE:
      break;
    case CTX_GROUP:
      if(*flags & DEBUG) printf("Entering the group context...\n");
      // coord sys
      tkn = strtok(l_buff, del);
      if(!strcmp(tkn, "polar")) ctx_flag |= POLAR_FLAG;
      else if(!strcmp(tkn, "cart")) ctx_flag |= CART_FLAG;
      // initialize arrays
      if(ctx_flag & POLAR_FLAG) pols = calloc(n_tmp, sizeof(POLAR));
      // Only one form of input - m, r, theta, phi, v, t_v, p_v
      // per instance value array
      double *vals = calloc(7, sizeof(double));
      int k;
      for(k = 0; k < 7; k++) {
	val = strtok(NULL, del_list);
	vals[k] = (double) atof(val);
      }
      if(ctx_flag & POLAR_FLAG) {
	pols[i-ptidx] = create_polar(vals, group);
        pol_to_pt(&(pols[i-ptidx]), 1, &(pts[i]));
      } else if(ctx_flag & CART_FLAG) 
	pts[i] = create_point(vals, group);
      i++;
      free(vals);
      if(i==(n_tmp+ptidx)) {
	ptidx = i;
	pols = realloc(pols, 0);
      }
      continue;
      break;
    case CTX_POINT:
      break;
    case CTX_RING:
      if(*flags & DEBUG) printf("Entering the ring context...\n");
      tkn = strtok(l_buff, del);
      if(!strcmp(tkn, "r_inner")) {
	r_inner = (double) atof(strtok(NULL, del));
	continue;
      }
      if(!strcmp(tkn, "width")) {
	width = (double) atof(strtok(NULL, del));
	continue;
      }
      if(!strcmp(tkn, "e_max")) {
	e_max = (double) atof(strtok(NULL, del));
	continue;
      }
      if(!strcmp(tkn, "height")) {
	height = (double) atof(strtok(NULL, del));
	continue;
      }
      if(!strcmp(tkn, "m_avg")) {
	m_avg = (double) atof(strtok(NULL, del));
	continue;
      }
      if(!strcmp(tkn, "m_sig")) {
	m_sig = (double) atof(strtok(NULL, del));
	continue;
      }
      if(!strcmp(tkn, "m_orb")) {
	m_orb = (double) atof(strtok(NULL, del));
	continue;
      }
      if(!strcmp(tkn, "GEN")) {
        if(*flags & DEBUG) printf("Generating ring...\n");
	create_ring(*flags, r_inner, width, e_max, height, m_avg, m_sig, m_orb, 
	            group, n_tmp, pts+ptidx);
        ptidx += n_tmp;
	continue;
      }
      break;
    case CTX_SATELLITES:
      if(*flags & DEBUG) printf("Entering the satellites context...\n");
      double mass, a_ax, ecc, inc, theta_a, theta_i, theta_0, p_mass;
      tkn = strtok(l_buff, del_list);
      mass = (double) atof(tkn);
      k = 0;
      while( (tkn = strtok(NULL, del_list)) != NULL ) {
	switch(k) {
	case 0:
	  a_ax = (double) atof(tkn);
	case 1:
	  ecc = (double) atof(tkn);
	case 2:
	  inc = (double) atof(tkn);
	case 3:
	  theta_a = (double) atof(tkn);
	case 4:
	  theta_i = (double) atof(tkn);
	case 5:
	  theta_0 = (double) atof(tkn);
	case 6:
	  p_mass = (double) atof(tkn);
	}
	k++;
      }
      create_satellite(*flags, mass, a_ax, ecc, inc, theta_a, theta_i, theta_0,
		       p_mass, group, &(pts[i]));
      i++;
      if(i==(ptidx+n_tmp)) ptidx = i;
      break;
    }
  }
  if(*flags & DEBUG) printf("Object read successfully\n");
  free(pols);

  return pts;
}

// opts helpers
void raise_error(unsigned short type, char *parent, char *param, char *val, int ln) {
  switch(type) {
  case VALUE_ERROR:
    printf("Options file error [line: %d]: %s does not allow value %s\n", ln,
	   param, val);
    exit(1);
  case NAME_ERROR:
    printf("Options file error:\n");
    exit(1);
  case NAME_WARNING:
    if(param) printf("Options file warning [line %d]: '%s' of group '%s' is unrecognized\n",
		     ln, param, parent);
    else printf("Options file warning [line: %d]: '%s' is an uncrecognized group\n", ln,
	   parent);
    break;
  case EMPTY_ERROR:
    if(param) printf("Options file error: %s of group %s expects a value\n",
		     param, parent);
    else printf("Options file error: group %s is empty or not found\n", parent);
    exit(1);
  case FMT_ERROR:
    printf("Options file error [line %d]: %s\n", ln, val);
    exit(1);
  }
}
