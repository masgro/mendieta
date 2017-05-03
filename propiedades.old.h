#ifndef PROP_H
#define PROP_H

#include <gsl/gsl_matrix.h>

struct propiedades_st
{
  int        npart;
  type_real  pcm[3], vcm[3], sig[3], mostbound[3];
  type_real  r200, m200, v200;
  type_real  rvir, mvir, vvir;
  type_real  vmax;
  type_real  L[3]; 
  type_real  Ep, Ec, lambda;
  type_real  aa, bb, cc;
  type_real  aa_vel, bb_vel, cc_vel;
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


struct propiedades_st propiedades(struct particle_data *P, int grupo);
void forma(char *flag, struct particle_data *Q, struct propiedades_st Prop);
void write_properties(FILE *pfout, struct propiedades_st Prop);
#endif
