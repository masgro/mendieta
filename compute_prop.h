#ifndef _COMPUTE_PROP_H_
#define _COMPUTE_PROP_H_

#include <gsl/gsl_matrix.h>
#include <stdbool.h>

struct subst{
  double pcm[3];
  double r200;
} *sub_hijos;

struct propiedades_st{
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

bool PROP_FOF;
bool PROP_SUB;

void compute_properties(struct grupos *g);
void propiedades(struct particle_data *P, struct grupos *g, int gid, struct propiedades_st *Prop);
void forma(char *flag, struct particle_data *Q, struct propiedades_st *Prop);
void write_properties(FILE *pfout, struct propiedades_st *Prop);
#endif
