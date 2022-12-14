#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "./lib/nbody.h"

int n = 0;
double eta = 0.0;

int main(int argc, char *argv[]) {
  unsigned int flags = 0;
  char opt, *input_name, *out_name;
  while( (opt = getopt(argc, argv, "cdf:o:")) != -1 ) {
    switch(opt) {
    case 'c':
      flags |= PRINT_CONS;
      break;
    case 'd':
      flags |= DEBUG;
      break;
    case 'f':
      input_name = optarg;
      break;
    case 'o':
      flags |= READ_FILE;
      out_name = optarg;
      break;
    default:
      break;
    }
  }

  FILE *fp_in, *fp_out;
  fp_in = fopen(input_name, "r");
  if(!fp_in) {
    printf("%s: error in opening file '%s'\n", argv[0], input_name);
    return(1);
  }

  if(flags & READ_FILE) fp_out = fopen(out_name, "w");

  printf("Beginning to read from file...\n");
  PARTICLE *pts;
  pts = obj_from_file(fp_in, &flags);
  printf("Finished reading object, printing...\n");
  printf("%d\n", flags & PRINT_CONS);
  if(flags & READ_FILE) {
    //print file
    //print_to_file(fp_out, flags, pts, n);
    print_frame_to_file(fp_out, flags, pts, n);
  }
  if(flags & PRINT_CONS) {
    //print to console
    printf("Output object:\n");
    int i;
    char *pt;
    for(i=0; i<n; i++) {
      pt = pt_to_str(pts[i], flags, pt);
      printf("%s\n", pt);
    }
  }
  
  return 0;
}
