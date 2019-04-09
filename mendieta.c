#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdbool.h>
#include <omp.h>

#include "variables.h"
#include "cosmoparam.h"
#include "limpieza.h"
#include "io.h"
#include "leesnap.h"
#include "deltas.h"
#include "timer.h"
#include "identest.h"
#include "compute_prop.h"
#include "colores.h"
#include "grid.h"

void linkedlist_grupos(int *head, unsigned int *npgrupo,int npart, int *ll, unsigned int *grupo);
void linkedlist();
void sub_groups(void);
void fof_groups(void);

int main(int argc, char **argv){
  #ifdef IDENSUB
  my_int i;
  #endif
  int    init_ifrac,ifrac;
  double frac,start,end;

  TIMER(start);
  
  init_variables(argc,argv);

  /*Lee archivos de la simulacion*/
  read_gadget();

  #ifdef READIDENFOF
  read_idenfof(fof_file);
  #endif

  /*Cambia origen de coordenadas*/
  change_positions(cp.npart);

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
    iden.r0 *= 0.17071;
    iden.r0 *= cbrt(cp.Mpart*1.0E10/cp.omegam/RHOCRIT)*1000.0;

    //if(ifrac == 0){
    //  iden.r0 = 0.7937;
    //  iden.r0 *= cbrt(cp.Mpart*1.0E10/cp.omegam/RHOCRIT)*1000.0;
    //}

    //if(ifrac == 1){
    //  iden.r0 = 0.17071;
    //  iden.r0 *= cbrt(cp.Mpart*1.0E10/cp.omegam/RHOCRIT)*1000.0;
    //}

    //if(ifrac > 1){
    //  exit(EXIT_SUCCESS);
    //}

    if(iden.r0 <= cp.soft){
      iden.r0 = cp.soft;
      ifrac = nfrac;
    }

    fprintf(stdout,"Linking length = %f \n",iden.r0);

    iden.nobj = cp.npart;

    identification();

    linkedlist();

    if(ifrac == 0)
      fof_groups();
    #ifdef IDENSUB
    else
      sub_groups();
    #endif

    #ifdef GETPOSITIONS
    get_positions(ifrac);
    #endif

    #ifdef GETINDEX
    get_index(ifrac);
    #endif

    //free(Temp.grup);
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

  //#if defined(COMPUTE_PROPERTIES)
  //PROP_FOF = true;
  //compute_properties(fof);
  //PROP_FOF = false;
  //#ifdef IDENSUB
  //PROP_SUB = true;
  //compute_properties(sub);
  //PROP_SUB = false;
  //#endif
  //#endif

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

void linkedlist(){
  my_int i;
  my_int g;

  printf("Inicia LinkedList %u %u\n",iden.nobj,iden.ngrupos);fflush(stdout);

  Temp.head   = (my_int *) malloc(iden.ngrupos*sizeof(my_int));
  Temp.npgrup = (my_int *) calloc(iden.ngrupos,sizeof(my_int));
  Temp.ll     = (my_int *) malloc(iden.nobj*sizeof(my_int));

  for(g = 0; g < iden.ngrupos; g++) Temp.head[g] = GROUND;

  for(i = 0; i < iden.nobj; i++){
    g = P[i].gr;
    #ifdef DEBUG
    assert((g >= 0) && (g < iden.ngrupos));
    #endif
    Temp.ll[i] = Temp.head[g];
    Temp.head[g] = i;
    Temp.npgrup[g]++;
  }

  #ifdef DEBUG
  for(g = 1; g < iden.ngrupos; g++)
    assert(Temp.npgrup[g] >= NPARTMIN);
  #endif

  ///i = 1;
  ///for(g = 1; g < iden.ngrupos; g++){
  ///  if(Temp.npgrup[g] < NPARTMIN){
  ///    Temp.npgrup[g] = 0;
  ///    Temp.head[g] = GROUND;
  ///    continue;
  ///  }
  ///  Temp.head[i] = Temp.head[g];
  ///  Temp.npgrup[i] = Temp.npgrup[g];
  ///  i++;
  ///}
  ///iden.ngrupos = i;

  ///Temp.head   = (my_int *) realloc(Temp.head,iden.ngrupos*sizeof(my_int));
  ///Temp.npgrup = (my_int *) realloc(Temp.npgrup,iden.ngrupos*sizeof(my_int));

  printf("End LinkedList\n");fflush(stdout);
}

void linkedlist_grupos(int *head, unsigned int *npgrupo,
                int npart, int *ll, unsigned int *grupo){
  int i;
  unsigned int g;
  
  for(i = 0; i < npart; i++){
    g = grupo[i];
    ll[i] = head[g];
    head[g] = i;
    npgrupo[g]++;
  }
}

void fof_groups(void){
  my_int i;

  n_grupos_fof = iden.ngrupos;

  /* Allocatea memoria */
  fof = (struct grupos *) malloc(n_grupos_fof*sizeof(struct grupos));
  assert(fof != NULL);

  for(i = 1; i < n_grupos_fof; i++){
    fof[i].llirst = Temp.head[i];
    fof[i].np = Temp.npgrup[i];
  }

  #ifdef limpiamelo
  for(i = 1; i < iden.ngrupos; i++)
    limpieza_new(i,0);
  #endif

  #ifdef IDENSUB
  n_grupos_sub = n_grupos_fof;
  sub = (struct grupos *) malloc(n_grupos_sub*sizeof(struct grupos));
  assert(sub != NULL);
  #endif

  Temp.nsub = n_grupos_fof;
  for(i = 0; i < iden.ngrupos; i++){
    fof[i].llirst = GROUND;
    fof[i].np     =  0;
    #ifdef IDENSUB
    sub[i].llirst = GROUND;
    sub[i].np     =  0;
    #endif
  }

  for(i = 0; i < cp.npart; i++){
    P[i].fof = P[i].gr;
    P[i].llfof = fof[P[i].gr].llirst;
    fof[P[i].gr].llirst = i;
    fof[P[i].gr].np++;
    #ifdef DEBUG
    assert(P[i].gr < n_grupos_fof);
    #endif

    #ifdef IDENSUB
    P[i].sub = P[i].gr;
    P[i].llsub = sub[P[i].sub].llirst;
    sub[P[i].sub].llirst = i;
    sub[P[i].sub].np++;
    #endif
  }
}

#ifdef IDENSUB
void sub_groups(void){
  my_int i, mysub;
  my_int contador_subgrupo;

  if(iden.ngrupos != 0)
    make_cleaning(iden.ngrupos,&contador_subgrupo);

  /* LE SUMAMOS 1 PARA TENER EN CUENTA EL GRUPO 0*/
  n_grupos_sub = contador_subgrupo + 1;

  free(sub);
  sub = (struct grupos *) malloc(n_grupos_sub*sizeof(struct grupos));
  assert(sub != NULL);

  for(i = 0; i < n_grupos_sub; i++){
    sub[i].llirst = GROUND;
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
