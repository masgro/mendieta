#ifndef VARIABLES_H
#define VARIABLES_H

#ifndef POSFACTOR
#define POSFACTOR 1.0
#endif

#ifndef VELFACTOR
#define VELFACTOR 1.0
#endif

#define KERNEL_TABLE 10000

#ifndef NPARTMIN
#define NPARTMIN 20
#endif


/* Precision del codigo (reales) */
#ifdef PRECDOUBLE
typedef double type_real;
#else
typedef float type_real;
#endif

/* Precision del codigo (enteros) */
#ifdef LONGIDS
typedef unsigned long long type_int;
#else
typedef unsigned int type_int;
#endif

size_t size_real;
size_t size_int;

/* Posiciones, velocidades y energias de las partículas */
struct particle_data {
  type_real    Pos[3];
  #ifdef STORE_VELOCITIES
  type_real    Vel[3];
  #endif
  #ifdef STORE_IDS
  type_int     id;
  #endif
  type_int     indx;
  int          gr;
  unsigned int fof;
  int          llfof;
  #ifdef IDENSUB
  unsigned int sub;
  int          llsub;
  #endif
  #ifdef ENERGIES
  double       Ep, Ec;
  #endif
} *P;

struct grupos{
  int llirst;
  unsigned int np;
} *fof, *sub;

unsigned int np_in_fof;
unsigned int np_in_sub;
unsigned int n_grupos_fof;
unsigned int n_grupos_sub;
int  nfrac;
char fof_file[200];
char sub_file[200];

void init_variables(int argc, char **argv);

#endif










