#ifndef LIMPIEZA_H
#define LIMPIEZA_H

#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector_long.h>
#include <gsl/gsl_sort_vector_long.h>

unsigned int npmax;

struct temporary {
  unsigned int nsub;
	unsigned int nsub_final;  
  unsigned int massive_one;
  unsigned int np_massive_one;
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
void compute_potential_energy_subgrupo(my_int npart,struct particles *Q, my_int *iEpmin);
void compute_cinetical_energy_subgrupo(my_int, struct particles *);
void compute_potential_energy(void);
void only_one_substructure(int k, int contador);
#endif
