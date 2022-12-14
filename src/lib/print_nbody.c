#include "nbody.h"

// vector printing function
void print_vector(double *v, int n) {
  printf("( ");
  for(int i=1;i<=n;i++) printf("%f ", v[i]);
  printf(")\n");
}

// printing options to a file
void print_opts(char *file, unsigned int flags, double soft, double ecc, char *f_in,
		double h, int out_f, char *d_out, char *f_out, int steps) {
  if(flags & DEBUG) printf("Printing options to file '%s'\n", file);
  
  FILE *fp;
  fp = fopen(file, "w");
  if(!fp) {
    printf("Error in printing options: file '%s' could not be opened\n", file);
    exit(1);
  }
  if(flags & DEBUG) printf("Opened file\n");

  //Sim
  fprintf(fp, "Simulator parameters.\n@SIM{\n\tsoftening:%f\n\tstepsize:%f", soft,
	  h);
  if(flags & DEBUG) printf("Line 1\n");
  fprintf(fp, "\n\tsteps:%d\n\tbh_mode:%d\n}\n", steps, (flags & MODE_BH)/MODE_BH);
  if(flags & DEBUG) printf("Line 2\n");
  //IO
  char *file_in, *file_out, *dir_out;
  if(flags & READ_FILE) file_in = f_in;
  else file_in = "none";
  if(flags & FLAG_31) file_out = f_out;
  else file_out = "none";
  if(flags & PRINT_DIR) dir_out = d_out;
  else dir_out = "none";
  
  fprintf(fp, "I/O parameters.\n@IO{\n\tfile_in:%s\n\tfile_out:%s\n\t", file_in,
	  file_out);
  if(flags & DEBUG) printf("Line 3\n");
  fprintf(fp, "dir_out:%s\n\tconsole:%d\n\tout_freq:%d\n}\n", dir_out,
	  (flags & PRINT_CONS)/PRINT_CONS, out_f);
  if(flags & DEBUG) printf("Line 4\n");
  
  //GEN
  char *mode, *mass;
  if(flags & MASS_NORM) mass = "uniform";
  else if(flags & MASS_UNI) mass = "uniform random";
  else if(flags & MASS_RAY) mass = "rayleigh";
  if(flags & GEN_2BODY) mode = "2body";
  else if(flags & GEN_SPHERE) mode = "sphere";
  else if(flags & GEN_DISK) mode = "disk";
  fprintf(fp, "Generator parameters.\n@GEN{\n\tmode:%s\n\tmass:%s\n}\n", mode, mass);
  if(flags & DEBUG) printf("Line 5\n");
  if(flags & DEBUG) printf("%d\n", flags & PRINT_CONS);

  fclose(fp);
}

// converting a point to a formatted string
char *pt_to_str(PARTICLE pt, unsigned int flags, char *str) {
  str = malloc(sizeof(char)*BUFF*3);
  if(flags & VERBOSE)
    snprintf(str, BUFF*3, "%.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %d",
	     pt.mass, pt.pos[0], pt.pos[1], pt.pos[2], pt.vel[0], pt.vel[1],
	     pt.vel[2], pt.acc[0], pt.acc[1], pt.acc[2], pt.group);
  else
    snprintf(str, BUFF*3, "%.6g %.6g %.6g %.6g %.6g %.6g %.6g", pt.mass, pt.pos[0],
	   pt.pos[1], pt.pos[2], pt.vel[0], pt.vel[1], pt.vel[2]);
  return str;
}

// printing a frame to a file in the format of pts
void print_frame_to_file(FILE *fp, unsigned int flags, PARTICLE *pts, int n) {
  int i;
  char *pt;
  for(i=0; i<n; i++) 
    fprintf(fp, "%s\n", pt_to_str(pts[i], flags, pt));
}
