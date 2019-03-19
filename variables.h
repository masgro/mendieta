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
#define NPARTMIN 10
#endif

/* Precision del codigo (reales) */
#ifdef PRECDOUBLE
typedef double my_real;
#else
typedef float my_real;
#endif

/* Precision del codigo (enteros) */
#ifdef LONGIDS
typedef unsigned long my_int;
#else
typedef unsigned int my_int;
#endif

size_t size_real;
size_t size_int;

/* Posiciones, velocidades y energias de las partículas */
struct particle_data {
  my_real Pos[3];
  #ifdef STORE_VELOCITIES
  my_real Vel[3];
  #endif
  #ifdef STORE_IDS
  my_int  id;
  #endif
  //my_int  indx;
  #ifdef READIDENFOF
  my_int  fof;
  #endif
  my_int  llfof;
  #ifdef IDENSUB
  my_int  sub;
  my_int  llsub;
  #endif
  #ifdef ENERGIES
  double  Ep, Ec;
  #endif
  my_int  gr;
} *P;

struct grupos{
  my_int llirst;
  my_int np;
} *fof, *sub, *groups;

my_int np_in_fof;
my_int np_in_sub;
my_int n_grupos_fof;
my_int n_grupos_sub;
int  nfrac;
char fof_file[200];
char sub_file[200];

void init_variables(int argc, char **argv);

#endif










