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

bool PROP_FOF;
bool PROP_SUB;

void compute_properties(struct grupos *g);
void propiedades(struct particle_data *P, struct grupos *g, int gid, struct propiedades_st *Prop);
void forma(char *flag, struct particles *Q, struct propiedades_st *Prop);
void write_properties(FILE *pfout, struct propiedades_st *Prop);
#endif
