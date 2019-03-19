#ifndef LIMPIEZA_H
#define LIMPIEZA_H

#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector_long.h>
#include <gsl/gsl_sort_vector_long.h>

my_int npmax;

struct temporary {
  my_int nsub;
	my_int nsub_final;  
  my_int massive_one;
  my_int np_massive_one;
	my_int np;          
	my_int *head;
	my_int *tail;
	my_int *ll;
	my_int *npgrup;
	my_int *grup;
	my_int *in;
	my_real **vcm, **pcm;
} Temp;

struct elemento {
  int id;
  struct elemento *next;
} *head_subgroups, *head_free_particles; 

void make_cleaning(unsigned int ngrupos, unsigned int *contador_subgrupo);
void compute_velocity_position(void);
void free_memory(void);
void reasigna(void);
void reasigna_closest(void);
void limpieza_new(my_int ig, my_int destino);
void compute_potential_energy_subgrupo(my_int npart,struct particle_data *Q, my_int *iEpmin);
void compute_potential_energy(void);
void only_one_substructure(int k, int contador);
#endif
