#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <omp.h>

#include "variables.h"
#include "cosmoparam.h"
#include "limpieza.h"
#include "io.h"
#include "leesnap.h"
#include "deltas.h"
#include "timer.h"
#include "iden.h"
#include "compute_prop.h"
#include "colores.h"
#include "peano.h"
#include "grid.h"

void linkedlist_grupos(int ngrupos, int *head, unsigned int *npgrupo, int npart, int *ll, unsigned int *grupo);
void sub_groups(void);
void fof_groups(void);

int main(int argc, char **argv){
  int    i,init_ifrac,ifrac;
  double frac,start,end;

  TIMER(start);
  
  init_variables(argc,argv);
  omp_set_nested(1);

  #ifdef SUBBOXES
  YELLOW("SUBBOXES ACTIVED\n");
  assert(argc >= 6);
  box.cen[0] = (my_real) atof(argv[2]);
  box.cen[1] = (my_real) atof(argv[3]);
  box.cen[2] = (my_real) atof(argv[4]);
  box.lado   = (my_real) atof(argv[5]);
  box.franja = 1000.0;
  box.min[0] = box.franja;
  box.max[0] = box.franja + 2.0*box.lado;
  box.min[1] = box.franja;
  box.max[1] = box.franja + 2.0*box.lado;
  box.min[2] = box.franja;
  box.max[2] = box.franja + 2.0*box.lado;
  #endif

  /*Lee archivos de la simulacion*/
  read_gadget();

  #ifdef READIDENFOF
  read_idenfof(fof_file);
  #endif

  /*Cambia origen de coordenadas*/
  change_positions(cp.npart);

  peano_hilbert();

  #ifdef COMPUTE_EP
  /*Calcula energia potencial de las particulas*/
  compute_potential_energy();
  #endif

  #ifndef IDENSUB
  nfrac = 0;
  #endif

  #ifdef READIDENFOF
  init_ifrac = 1;
  #else
  init_ifrac = 0;
  n_grupos_sub = 0;
  #endif

  #ifdef IDENSUB
  for(i = 0; i < cp.npart; i++) P[i].sub = 0;
  #endif

  for(ifrac = init_ifrac; ifrac <= nfrac; ifrac++){
    fprintf(stdout, "\n Begins Identification : Step %d of %d \n",ifrac,nfrac);
    
    frac  = 1.0f/(float)(nfrac + 1);
    frac *= (float)(nfrac + 1 - ifrac);
 
    iden.r0  = frac;
    iden.r0 *= 0.2;
    iden.r0 *= cbrt(cp.Mpart*1.0E10/cp.omegam/RHOCRIT)*1000.0;

    #ifdef READIDENFOF
    iden.step  = 0;
    #else
    iden.step  = ifrac;
    #endif

    if(iden.r0 <= cp.soft)
    {
      iden.r0 = cp.soft;
      ifrac = nfrac;
    }

    fprintf(stdout,"Linking length = %f \n",iden.r0);

    iden.nobj = cp.npart;

    identification();

    if(ifrac == 0)
      fof_groups();
    #ifdef IDENSUB
    else
      sub_groups();
    #endif

    #ifdef GETPOSITIONS
    get_positions(ifrac);
    #endif

    free(Temp.head);
    free(Temp.npgrup);
    free(Temp.ll);
  }

  /************* TERMINO LA IDENTIFICACION ***************/
  #ifndef READIDENFOF
  write_idenfof(fof_file);
  #endif

  #ifdef IDENSUB
  write_idensub(sub_file);
  #endif

  #ifdef COMPUTE_FOF_PROPERTIES
  PROP_FOF = true;
  compute_properties(fof);
  PROP_FOF = false;
  #endif

  #ifdef IDENSUB
  #ifdef COMPUTE_SUB_PROPERTIES
  PROP_SUB = true;
  compute_properties(sub);
  PROP_SUB = false;
  #endif
  #endif

  free(fof);
  #ifdef IDENSUB
  free(sub);
  #endif
  free(P);

  grid_free();

  TIMER(end);
  printf("Total time %f\n",end-start);

  return(EXIT_SUCCESS);
}

void linkedlist_grupos(int ngrupos, int *head, unsigned int *npgrupo,
                int npart, int *ll, unsigned int *grupo)
                {
  int i;
  unsigned int g;
  
  for(i = 0; i < npart; i++)
  {
    g = grupo[i];
    ll[i] = head[g];
    head[g] = i;
    npgrupo[g]++;
  }
}

void fof_groups(void){
  int i;

  n_grupos_fof = iden.ngrupos;

  /* Allocatea memoria */
  fof = (struct grupos *) malloc(n_grupos_fof*sizeof(struct grupos));
  assert(fof != NULL);
  #ifdef IDENSUB
  n_grupos_sub = iden.ngrupos;
  sub = (struct grupos *) malloc(n_grupos_sub*sizeof(struct grupos));
  assert(sub != NULL);
  #endif

  Temp.nsub = n_grupos_fof;
  #ifdef limpiamelo
  for(i = 1; i < iden.ngrupos; i++)
    limpieza_new(i,0);
  #endif

  for(i = 0; i < iden.ngrupos; i++)
  {
    fof[i].llirst = -1;
    fof[i].np     =  0;
    #ifdef IDENSUB
    sub[i].llirst = -1;
    sub[i].np     =  0;
    #endif
  }

  for(i = 0; i < cp.npart; i++){
    P[i].fof = P[i].gr;
    #ifdef DEBUG
    assert(P[i].fof < n_grupos_fof);
    #endif
    P[i].llfof = fof[P[i].fof].llirst;
    fof[P[i].fof].llirst = i;
    fof[P[i].fof].np++;

    #ifdef IDENSUB
    P[i].sub = P[i].fof;
    #ifdef DEBUG
    assert(P[i].sub < n_grupos_sub);
    #endif
    P[i].llsub = sub[P[i].sub].llirst;
    sub[P[i].sub].llirst = i;
    sub[P[i].sub].np++;
    #endif
  }
}

#ifdef IDENSUB
void sub_groups(void)
{
  int i, mysub;
  unsigned int contador_subgrupo;

  if(iden.ngrupos != 0)
    make_cleaning(iden.ngrupos,&contador_subgrupo);

  /* LE SUMAMOS 1 PARA TENER EN CUENTA EL GRUPO 0*/
  n_grupos_sub = contador_subgrupo + 1;

  free(sub);
  sub = (struct grupos *) malloc(n_grupos_sub*sizeof(struct grupos));
  assert(sub != NULL);

  for(i = 0; i < n_grupos_sub; i++)
  {
    sub[i].llirst = -1;
    sub[i].np = 0;
  }

  for(i = 0; i < cp.npart; i++){
    mysub = P[i].sub;
    #ifdef DEBUG
    assert(mysub < n_grupos_sub);
    #endif
    P[i].llsub = sub[mysub].llirst;
    sub[mysub].llirst = i;
    sub[mysub].np++;
  }
}
#endif
