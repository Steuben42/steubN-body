#include "nbody.h"
#include <math.h>
#include <stdio.h>

/***********************
  BARNES-HUT CONSTANTS
***********************/
#define THETA_CRIT 0.05 // critical angle l/d before approximation
// should be 0.05

// Will prompt the user if this value has been changed
void confirm_critical() {
  if(THETA_CRIT != 0.05) {
    printf("WARNING: theta critical in bhmode_nbody.c is set to a value higher \
than 0.05. This will improve performance, but may seriously degrade results. \
Continue? (y/n)\n");
    char op;
    op = getchar();
    while((op!='y') && (op!='n')) printf("Please answer (y)es or (n)o.\n");
    if(op=='y') {
    } else if(op=='n') {
      printf("Aborting.\n");
      exit(1);
    }
  }
}

/************************
  STRUCT CREATION FUNCS
************************/
// root node creation
NODE *create_root_node(PARTICLE *pts, int n) {
  NODE *root;
  root = (NODE *) malloc(sizeof(NODE));
  int i, k;
  double max = 0.0;
  double min = 0.0;
  // This loop and the following increment will ensure that the
  // root node is just slightly larger than the bounds of the particles
  for(i = 0; i < n; i++) {
    for(k = 0; k < 3; k++) {
      if(max < pts[i].pos[k]) max = pts[i].pos[k];
      if(min > pts[i].pos[k]) min = pts[i].pos[k];
    }
  }
  max += 0.01;
  min -= 0.01;

  // Ensure that other important values are initialized
  root->order = 0;
  root->size = max - min;
  for(k = 0; k < 3; k++) root->pos[k] = max - root->size/2.0;

  // Make sure that the children are NULL pointers
  for(i = 0; i < 8; i++) {
    root->c_node[i] = NULL;
    root->c_pt[i] = NULL;
  }
  root->pts = 0;
  root->id = 0;

  return root;
}

// unused
NODE *create_root_node_static(double size[], double pos[]) {
  NODE *root;
  return root;
}

// node creation
NODE *create_node(NODE *p_node, int oct) {
  NODE *node;
  static int n_id = 1;
  node = (NODE *) malloc(sizeof(NODE));
  node->order = p_node->order + 1;
  node->size = p_node->size/2.0;
  node->mass = 0.0;
  double sign[3];
  // Find which octant this node is in, and flip the signs accordingly so
  // that its position values will be accurate
  if(oct > 3) sign[0] = -1.0;
  else sign[0] = 1.0;
  if( (oct == 2) || (oct == 3) || (oct == 6) || (oct == 7) ) sign[1] = -1.0;
  else sign[1] = 1.0;
  if(oct%2) sign[2] = -1.0;
  else sign[2] = 1.0;

  // set the position values and default com to the center of the node
  int i;
  for(i = 0; i < 3; i++) {
    node->pos[i] = p_node->pos[i] + sign[i]*node->size/2.0;
    node->com[i] = p_node->pos[i] + sign[i]*node->size/2.0;
  }

  // refer to the parent node
  node->p_node = p_node;

  // initialize NULL children pointers
  for(i = 0; i < 8; i++) {
    node->c_node[i] = NULL;
    node->c_pt[i] = NULL;
  }
  node->pts = 0;

  // ensure that if this is the first non-root node created, the id will be 1
  if(p_node->id==0 & p_node->pts==0) n_id = 1;

  node->id = n_id++;
  return node;
}

// particle creation
PARTICLE create_point(double *vals, int group) {
  static int pt_id = 0;
  PARTICLE pt;
  pt.id = pt_id++;
  pt.group = group;

  // passing the initialization values of the array vals
  int k;
  for(k = 0; k < 7; k++) {
    if(k==0) pt.mass = vals[k];
    else if(k<4) pt.pos[k-1] = vals[k];
    else pt.vel[k-4] = vals[k];
  }
  
  return pt;
}

// polar creation
POLAR create_polar(double *vals, int group) {
  PARTICLE pt;
  POLAR pol;
  pt = create_point(vals, group);
  
  pol.id = pt.id;
  pol.group = group;
  
  // same code as particle (essentially same type, separated to ensure
  // no conflicting types)
  int k;
  for(k = 0; k < 7; k++) {
    if(k==0) pol.mass = vals[k];
    else if(k<4) pol.pos[k-1] = vals[k];
    else pol.vel[k-4] = vals[k];
  }
  
  return pol;
}


void put_pt_in_tree(PARTICLE *pt, NODE *node) {
  // use >= and < to ensure no border cases
  int i, oct;
  double pos_rel[3];

  // determining the octant this point falls in
  for(i = 0; i < 3; i++) pos_rel[i] = pt->pos[i] - node->pos[i];
  if(pos_rel[0] >= 0.0) { // +
    if(pos_rel[1] >= 0.0) { // ++
      if(pos_rel[2] >= 0.0) oct = 0; // +++
      else oct = 1; // ++-
    } else {
      if(pos_rel[2] >= 0.0) oct = 2; // +-+
      else oct = 3; // +--
    }
  } else {
    if(pos_rel[1] >= 0.0) {
      if(pos_rel[2] >= 0.0) oct = 4; // -++
      else oct = 5; // -+-
    } else {
      if(pos_rel[2] >= 0.0) oct = 6; // --+
      else oct = 7; // ---
    }
  }

  node->pts++;

  // first, make a leaf if there is no child node or point
  // second, make a branch if there is a child point
  // third, sort the point into the branch if there is a child node
  if(!node->c_pt[oct] && !node->c_node[oct]) {
    node->c_pt[oct] = pt;
  } else if(!(node->c_node[oct])) {
    node->c_node[oct] = create_node(node, oct);
    put_pt_in_tree(node->c_pt[oct], node->c_node[oct]);
    // there should be no direct children points if this is a branch
    node->c_pt[oct] = NULL;
    put_pt_in_tree(pt, node->c_node[oct]);
  } else put_pt_in_tree(pt, node->c_node[oct]);
}

// centers of mass
// should this be its own function? Could we wrap this into put_pt_in_tree()
// for increase performance?
void calc_moment(NODE *node) {
  int i, k;
  for(i = 0; i < 8; i++) {
    // if it is a leaf, accumulate the mass to this node
    // if it is a branch, recurse the function and then accumulate
    // the branch's mass to this node
    if(node->c_pt[i]) {
      node->mass += node->c_pt[i]->mass;
      for(k = 0; k < 3; k++) node->com[k] += (node->c_pt[i]->pos[k]
					      *node->c_pt[i]->mass);
    } else if(node->c_node[i]) {
      calc_moment(node->c_node[i]);
      node->mass += node->c_node[i]->mass;
      for(k = 0; k < 3; k++) node->com[k] += (node->c_node[i]->com[k]
					      *node->c_node[i]->mass);
    }
  }
  for(k = 0; k < 3; k++) node->com[k] /= node->mass;
}

NODE *build_tree(PARTICLE *pts) {
  NODE *root;
  // simply create a root node and build off of it
  root = create_root_node(pts, n);
  int i;
  for(i = 0; i < n; i++) put_pt_in_tree(&(pts[i]), root);
  calc_moment(root);

  return root;
}

void add_force(PARTICLE *pt, NODE *node) {
  // find the relative angle
  double theta = node->size/mag_diff(pt->pos, node->pos, 3);

  // set this point's accelration to 0.0 as it will accumulate later
  int i, k;
  if(node->id==0)
    for(k = 0; k < 3; k++)
      pt->acc[k] = 0.0;

  // find if we can approximate; if so, then use the node to calculate
  // if we can't approximate, then check each octant
  // if there is a child point that is not this point, then directly
  // calculate
  // if there is a child node, recurse this function
  if(theta < THETA_CRIT) {
    for(k = 0; k < 3; k++) {
      // -G*mass*mass*positiondiff/(mag)**3
      pt->acc[k] += (-G*node->mass*(pt->pos[k] - node->com[k])
		     /pow(sqr(mag_diff(pt->pos, node->com, 3)) + sqr(eta), 1.5));
    }
  } else {
    for(i = 0; i < 8; i++) {
      if(node->c_pt[i]) {
	if(pt->id==node->c_pt[i]->id) continue;
	for(k = 0; k < 3; k++) {
	  pt->acc[k] += (-G*node->c_pt[i]->mass
			 *(pt->pos[k] - node->c_pt[i]->pos[k])
			 /pow(sqr(mag_diff(pt->pos, node->c_pt[i]->pos, 3))
			      + sqr(eta), 1.5));
	}
      } else if(node->c_node[i]) {
	add_force(pt, node->c_node[i]);
      }
    }
  }
}

// ensure that the tree is fully released following its use
void free_tree(NODE *node) {
  int i;
  for(i = 0; i < 8; i++) {
    if(node->c_node[i]) free_tree(node->c_node[i]);
  }
  free(node);
}
