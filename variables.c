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

  #ifdef SUBBOXES
  if(argc < 6){
    fprintf(stdout,"Usage: %s input.conf xc yc zc lbox\n",argv[0]);
    exit(EXIT_FAILURE);
  }
  #else
  if(argc < 2){
    fprintf(stdout,"Usage: %s input.conf\n",argv[0]);
    exit(EXIT_FAILURE);
  }
  #endif

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
  
  #ifdef SUBBOXES
  YELLOW("SUBBOXES ACTIVED\n");
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

  BLUE("********** Information ***********\n");
  sprintf(message,"Snapshots directory:     %s\n",snap.root);BLUE(message);
  sprintf(message,"Snapname:                %s\n",snap.name);BLUE(message);
  sprintf(message,"# of snapshots:          %d\n",snap.nfiles);BLUE(message);
  sprintf(message,"Out file fof groups:     %s\n",fof_file);BLUE(message);
  sprintf(message,"Out file sub groups:     %s\n",sub_file);BLUE(message);
  sprintf(message,"Softening of simulation: %lf \n",cp.soft);BLUE(message);
  sprintf(message,"Identification steps:    %d\n",nfrac);BLUE(message);
  BLUE("**********************************\n");
  RED("******** Makefile Options ********\n");
  #ifdef DEBUG
  RED("  DEBUG\n");
  #endif
  #ifdef COMPUTE_EP
  RED("  COMPUTE_EP\n");
  #endif
  #ifdef GETPOSITIONS
  RED("  GETPOSITIONS\n");
  #endif
  #ifdef SUBBOXES
  RED("  SUBBOXES\n");
  #endif
  #ifdef PERIODIC
  RED("  PERIODIC\n");
  #endif
  #ifdef PRECDOUBLE
  RED("  PRECDOUBLE\n");
  #endif
  #ifdef LONGIDS
  RED("  LONGIDS\n");
  #endif
  #ifdef POSFACTOR
  sprintf(message,"  POSFACTOR = %f\n",POSFACTOR);RED(message);
  #endif
  #ifdef VELFACTOR
  sprintf(message,"  VELFACTOR = %f\n",VELFACTOR);RED(message);
  #endif
  #ifdef STORE_VELOCITIES
  RED("  STORE_VELOCITIES\n");
  #endif
  #ifdef ENERGIES
  RED("  ENERGIES\n");
  #endif
  #ifdef READIDENFOF
  RED("  READIDENFOF\n");
  #endif
  #ifdef NHALO
  sprintf(message,"  NHALO = %d\n",NHALO);RED(message);
  #endif
  #ifdef limpiamelo
  RED("  limpiamelo\n");
  #endif
  #ifdef COMPUTE_FOF_PROPERTIES
  RED("  COMPUTE_FOF_PROPERTIES\n");
  #endif
  #ifdef IDENSUB
  RED("  IDENSUB\n");
  #endif
  #ifdef ASSIGN_CLOSEST_GROUP
  RED("  ASSIGN_CLOSEST_GROUP\n");
  #endif
  #ifdef PROP
  RED("  PROP\n");
  #endif
  RED("**********************************\n");
}
