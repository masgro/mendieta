#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "variables.h"
#include "cosmoparam.h"
#include "leesnap.h"
#include "colores.h"

void init_variables(int argc, char **argv){
  FILE *pfin;
  char filename[200];
  int ierr;

  RED("Initializing variables...\n");

  size_real = sizeof(my_real);
  fprintf(stdout,"size_real: %zu\n",size_real);

  size_int  = sizeof(my_int);
  fprintf(stdout,"size_int:  %zu\n",size_int);

  sprintf(filename,"%s",argv[1]);
  pfin = fopen(filename,"r");
  assert(pfin != NULL);

  ierr = fscanf(pfin,"%d  \n",&snap.nfiles);
  assert(ierr == 1);
  ierr = fscanf(pfin,"%s  \n",snap.root);
  assert(ierr == 1);
  ierr = fscanf(pfin,"%s  \n",snap.name);
  assert(ierr == 1);
  ierr = fscanf(pfin,"%s  \n",fof_file);
  assert(ierr == 1);
  ierr = fscanf(pfin,"%s  \n",sub_file);
  assert(ierr == 1);
  ierr = fscanf(pfin,"%lf \n",&cp.soft);
  assert(ierr == 1);
  ierr = fscanf(pfin,"%d  \n",&nfrac);
  assert(ierr == 1);
  fclose(pfin);
  
  BLUE("********** Information ***********\n");
  sprintf(message,"Snapshots directory:     %s\n",snap.root);BLUE(message);
  sprintf(message,"Snapname:                %s\n",snap.name);BLUE(message);
  sprintf(message,"# of snapshots:          %d\n",snap.nfiles);BLUE(message);
  sprintf(message,"Out file fof groups:     %s\n",fof_file);BLUE(message);
  sprintf(message,"Out file sub groups:     %s\n",sub_file);BLUE(message);
  sprintf(message,"Softening of simulation: %lf \n",cp.soft);BLUE(message);
  sprintf(message,"Identification steps:    %d\n",nfrac);BLUE(message);
  BLUE("********** Makefile Options ***********\n");
  #ifdef DEBUG
  BLUE("  DEBUG\n");
  #endif
  #ifdef COMPUTE_EP
  BLUE("  COMPUTE_EP\n");
  #endif
  #ifdef GETPOSITIONS
  BLUE("  GETPOSITIONS\n");
  #endif
  #ifdef SUBBOXES
  BLUE("  SUBBOXES\n");
  #endif
  #ifdef PERIODIC
  BLUE("  PERIODIC\n");
  #endif
  #ifdef PRECDOUBLE
  BLUE("  PRECDOUBLE\n");
  #endif
  #ifdef LONGIDS
  BLUE("  LONGIDS\n");
  #endif
  #ifdef POSFACTOR
  sprintf(message,"  POSFACTOR = %f\n",POSFACTOR);BLUE(message);
  #endif
  #ifdef VELFACTOR
  sprintf(message,"  VELFACTOR = %f\n",VELFACTOR);BLUE(message);
  #endif
  #ifdef STORE_VELOCITIES
  BLUE("  STORE_VELOCITIES\n");
  #endif
  #ifdef ENERGIES
  BLUE("  ENERGIES\n");
  #endif
  #ifdef READIDENFOF
  BLUE("  READIDENFOF\n");
  #endif
  #ifdef NHALO
  sprintf(message,"  NHALO = %d\n",NHALO);BLUE(message);
  #endif
  #ifdef limpiamelo
  BLUE("  limpiamelo\n");
  #endif
  #ifdef COMPUTE_FOF_PROPERTIES
  BLUE("  COMPUTE_FOF_PROPERTIES\n");
  #endif
  #ifdef IDENSUB
  BLUE("  IDENSUB\n");
  #endif
  #ifdef ASSIGN_CLOSEST_GROUP
  BLUE("  ASSIGN_CLOSEST_GROUP\n");
  #endif
  #ifdef PROP
  BLUE("  PROP\n");
  #endif

  GREEN("END\n");
}
