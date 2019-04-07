#ifndef PROP_H
#define PROP_H

#include <gsl/gsl_matrix.h>

struct propiedades_st
{
  int        npart;
  my_real  pcm[3], vcm[3], sig[3], mostbound[3];
  my_real  r200, m200, v200;
  my_real  rvir, mvir, vvir;
  my_real  vmax;
  my_real  L[3]; 
  my_real  Ep, Ec, lambda;
  my_real  aa, bb, cc;
  my_real  aa_vel, bb_vel, cc_vel;
  gsl_matrix *evec;
  gsl_matrix *evec_vel;
};

struct node
{
  int indx;
  struct node *next;
};

struct stuff
{
  struct node *first;
  unsigned int np;
} *groups;


struct propiedades_st propiedades(struct particle_data *P, struct grupos *g, int gid);
void forma(char *flag, struct particle_data *Q, struct propiedades_st Prop);
void write_properties(FILE *pfout, struct propiedades_st Prop);
#endif
