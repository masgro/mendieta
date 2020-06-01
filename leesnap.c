#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "variables.h"
#include "allocate.h"
#include "leesnap.h"
#include "colores.h"

static struct io_header header;

static void leeheader(const char *filename)
{
  FILE *pf;
  type_int d1,d2;

  pf = fopen(filename,"r");
  if(pf == NULL){
    fprintf(stderr,"can't open file `%s`\n",filename);
    exit(EXIT_FAILURE);
  }

  fread(&d1, sizeof(d1), 1, pf);
  fread(&header, sizeof(header), 1, pf);
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);
  fclose(pf);

  // Definicion estructura cosmoparam
  cp.omegam    = header.Omega0;
  cp.omegal    = header.OmegaLambda;
  cp.omegak    = 1.0 - cp.omegam - cp.omegal;
  cp.hparam    = header.HubbleParam;
  cp.lbox      = header.BoxSize;
  cp.npart     = header.npartTotal[1];
  cp.Mpart     = header.mass[1];
  //cp.Mpart    = 3.143E-4*cp.hparam;  /*A-5*/
  //cp.Mpart    = 3.929E-5*cp.hparam;  /*A-4*/
  cp.redshift  = header.redshift;
  cp.aexp      = ( 1.0 / ( 1.0 + cp.redshift ) );
  cp.Hubble_a  = cp.omegam/cp.aexp/cp.aexp/cp.aexp;
  cp.Hubble_a += cp.omegak/cp.aexp/cp.aexp;
  cp.Hubble_a += cp.omegal;
  cp.Hubble_a  = sqrt(cp.Hubble_a);
  cp.Hubble_a *= 100.0*cp.hparam;

  printf("*********************************** \n");
  printf("*   Parametros de la simulacion   * \n");
  printf("*********************************** \n");
  printf("  Numero de particulas = %u \n", cp.npart);
  printf("  Lado del box = %g \n", cp.lbox);
  printf("  Redshift = %g \n", cp.redshift);
  printf("  Omega Materia = %g \n", cp.omegam);
  printf("  Omega Lambda = %g \n", cp.omegal);
  printf("  Parametro de Hubble = %g \n",cp.hparam);
  printf("  Masa por particula = %g \n",cp.Mpart);
  printf("  Softening = %g\n",cp.soft);
  printf("*********************************** \n");
  printf("*********************************** \n");
}

static void lee(const char *filename, type_int * restrict ind)
{
  FILE *pf;
  type_int d1, d2;
  type_int k, pc, n;

  type_real r[3];
  #ifdef STORE_VELOCITIES
  type_real v[3];
  #endif
  #ifdef STORE_IDS
  type_int id;
  #endif

  pf = fopen(filename,"r");
  if(pf == NULL){
    fprintf(stderr,"can't open file `%s`\n",filename);
    exit(EXIT_FAILURE);
  }

  fprintf(stdout,"Reading file: %s \n",filename); fflush(stdout);

  fread(&d1, sizeof(d1), 1, pf);
  fread(&header, sizeof(header), 1, pf);
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);

  fread(&d1, sizeof(d1), 1, pf);
  for(k = 0, pc = *ind; k < N_part_types; k++){
    for(n = 0; n < header.npart[k]; n++){
      fread(&r[0], sizeof(type_real), 3, pf);
      if(k == 1){ /*ONLY KEEP DARK MATTER PARTICLES*/
      #ifdef COLUMN
        P.x[pc] = r[0]*POSFACTOR;
        P.y[pc] = r[1]*POSFACTOR;
        P.z[pc] = r[2]*POSFACTOR;

#ifdef CHANGE_POSITION
        if(P.x[pc] > pmax[0]) pmax[0] = P.x[pc];
        if(P.x[pc] < pmin[0]) pmin[0] = P.x[pc];
        if(P.y[pc] > pmax[1]) pmax[1] = P.y[pc];
        if(P.y[pc] < pmin[1]) pmin[1] = P.y[pc];
        if(P.z[pc] > pmax[2]) pmax[2] = P.z[pc];
        if(P.z[pc] < pmin[2]) pmin[2] = P.z[pc];
#endif
      #else
        P[pc].pos[0] = r[0]*POSFACTOR;
        P[pc].pos[1] = r[1]*POSFACTOR;
        P[pc].pos[2] = r[2]*POSFACTOR;

#ifdef CHANGE_POSITION
        if(P[pc].pos[0] > pmax[0]) pmax[0] = P[pc].pos[0];
        if(P[pc].pos[0] < pmin[0]) pmin[0] = P[pc].pos[0];
        if(P[pc].pos[1] > pmax[1]) pmax[1] = P[pc].pos[1];
        if(P[pc].pos[1] < pmin[1]) pmin[1] = P[pc].pos[1];
        if(P[pc].pos[2] > pmax[2]) pmax[2] = P[pc].pos[2];
        if(P[pc].pos[2] < pmin[2]) pmin[2] = P[pc].pos[2];
#endif
      #endif
        pc++;
      }
    }
  }
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);

  fread(&d1, sizeof(d1), 1, pf);
#ifdef STORE_VELOCITIES
  for(k = 0, pc = *ind; k < N_part_types; k++){
    for(n = 0; n < header.npart[k]; n++){
      fread(&v[0], sizeof(type_real), 3, pf);
      if(k == 1){ /*ONLY KEEP DARK MATTER PARTICLES*/
      #ifdef COLUMN        
        P.vx[pc] = v[0]*VELFACTOR;
        P.vy[pc] = v[1]*VELFACTOR;
        P.vz[pc] = v[2]*VELFACTOR;
      #else
        P[pc].vel[0] = v[0]*VELFACTOR;
        P[pc].vel[1] = v[1]*VELFACTOR;
        P[pc].vel[2] = v[2]*VELFACTOR;
      #endif
        pc++;
      }
    }
  }
#else
  fseek(pf,d1,SEEK_CUR);
#endif
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);

  fread(&d1, sizeof(d1), 1, pf);
#ifdef STORE_IDS
  for(k = 0, pc = *ind; k < N_part_types; k++){
    for(n = 0; n < header.npart[k]; n++){
      fread(&id, sizeof(type_int), 1, pf);
      if(k == 1){ /*ONLY KEEP DARK MATTER PARTICLES*/
      #ifdef COLUMN
        P.id[pc] = id;
      #else
        P[pc].id = id;
      #endif
        pc++;
      }
    }
  }
#else
  fseek(pf,d1,SEEK_CUR);
#endif
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);

  *ind = pc;
  
  fclose(pf);
}

extern void read_gadget(void)
{
  char filename[200];
  type_int ifile, ind;
  size_t total_memory;

#ifdef CHANGE_POSITION
  for(ind = 0; ind < 3; ind++)
  {
    pmin[ind] =  1.E26; 
    pmax[ind] = -1.E26;
  }
#endif

  if(snap.nfiles>1)
    sprintf(filename,"%s%s.0",snap.root,snap.name);
  else
    sprintf(filename,"%s%s",snap.root,snap.name);

  leeheader(filename);

  /****** ALLOCATACION TEMPORAL DE LAS PARTICULAS ****************/
  total_memory = (float)cp.npart*sizeof(struct particle_data)/1024.0/1024.0/1024.0;
  printf("Allocating %.5zu Gb for %u particles\n",total_memory,cp.npart);
  if(!allocate_particles(&P, cp.npart))  exit(1);

  /***** LEE POS Y VEL DE LAS PARTICULAS ***********************/
  for(ifile = 0, ind = 0; ifile < snap.nfiles; ifile++)
  {
    if(snap.nfiles>1)
      sprintf(filename,"%s%s.%d",snap.root,snap.name,ifile);
    else
      sprintf(filename,"%s%s",snap.root,snap.name);

    lee(filename,&ind);
  }

  cp.lbox *= POSFACTOR;

  fprintf(stdout,"cp.lbox %f....\n",cp.lbox);
  fprintf(stdout,"End reading snapshot file(s)...\n"); fflush(stdout);
}

#ifdef CHANGE_POSITION

extern void change_positions(type_int n)
{
  type_int ip;

  printf("xmin %.1f xmax %.1f\n",pmin[0],pmax[0]);
  printf("ymin %.1f ymax %.1f\n",pmin[1],pmax[1]);
  printf("zmin %.1f zmax %.1f\n",pmin[2],pmax[2]);

  for(ip=0;ip<n;ip++)
  {
#ifdef COLUMN
    P.x[ip] -= pmin[0];
    P.y[ip] -= pmin[1];
    P.z[ip] -= pmin[2];
#else
    P[ip].pos[0] -= pmin[0];
    P[ip].pos[1] -= pmin[1];
    P[ip].pos[2] -= pmin[2];
#endif
  }

  cp.lbox = 0.0f;
  for(ip = 0; ip < 3; ip++)
    if(cp.lbox < (pmax[ip] - pmin[ip])) 
      cp.lbox = (pmax[ip] - pmin[ip]);

  cp.lbox *= 1.001;

  fprintf(stdout,"Changing cp.lbox %f....\n",cp.lbox);

}

#endif
