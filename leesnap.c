#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "variables.h"
#include "cosmoparam.h"
#include "leesnap.h"
#include "colores.h"

void read_gadget(void){
  char   filename[200];
  int    ifile;
  my_int ind;
  size_t total_memory;

  if(snap.nfiles>1)
    sprintf(filename,"%s%s.0",snap.root,snap.name);
  else
    sprintf(filename,"%s%s",snap.root,snap.name);

  leeheader(filename);

  /****** ALLOCATACION TEMPORAL DE LAS PARTICULAS ****************/
  total_memory = (float)cp.npart*sizeof(struct particle_data)/1024.0/1024.0/1024.0;
  printf("Allocating %.5zu Gb for %lu particles\n",total_memory,(unsigned long)cp.npart);
  P = (struct particle_data *) malloc(cp.npart*sizeof(struct particle_data));
  assert(P != NULL);

  /***** LEE POS Y VEL DE LAS PARTICULAS ***********************/
  for(ifile = 0, ind = 0; ifile < snap.nfiles; ifile++){
    if(snap.nfiles>1)
      sprintf(filename,"%s%s.%d",snap.root,snap.name,ifile);
    else
      sprintf(filename,"%s%s",snap.root,snap.name);

    lee(filename,P,&ind);
  }

#ifdef SUBBOXES
  select_particles();
#endif

  fprintf(stdout,"End reading snapshot file(s)...\n"); fflush(stdout);
}

void leeheader(char *filename){
  FILE *pf;
  int d1,d2;
  size_t ierr;

  pf = fopen(filename,"r");
  if(pf == NULL){
    fprintf(stderr,"can't open file `%s`\n",filename);
    exit(EXIT_FAILURE);
  }

  ierr = fread(&d1, sizeof(d1), 1, pf);
  assert(ierr == 1);
  ierr = fread(&header, sizeof(header), 1, pf);
  assert(ierr == 1);
  ierr = fread(&d2, sizeof(d2), 1, pf);
  assert(ierr == 1);
  assert(d1==d2);
  fclose(pf);

  /* Definicion estructura cosmoparam */
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
  printf("  Numero de particulas = %lu\n",(unsigned long)cp.npart);
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

void lee(char *filename, struct particle_data *Q, my_int *ind){
  FILE *pf;
  int d1, d2;
  int k, pc, n;
  size_t ierr;

  my_real r[3],v[3];
  #ifdef STORE_IDS
  my_int id;
  #endif
 
  for(k = 0; k < 3; k++){
    pmin[k] = 1.E26; 
    pmax[k] = -1.E26;
  }

  pf = fopen(filename,"r");
  if(pf == NULL){
    fprintf(stderr,"can't open file `%s`\n",filename);
    exit(EXIT_FAILURE);
  }

  fprintf(stdout,"Reading file: %s \n",filename); fflush(stdout);

  ierr = fread(&d1, sizeof(d1), 1, pf);
  assert(ierr == 1);
  ierr = fread(&header, sizeof(header), 1, pf);
  assert(ierr == 1);
  ierr = fread(&d2, sizeof(d2), 1, pf);
  assert(ierr == 1);
  assert(d1==d2);

  ierr = fread(&d1, sizeof(d1), 1, pf);
  assert(ierr == 1);
  for(k = 0, pc = 0; k < 6; k++){
    for(n = 0; n < header.npart[k]; n++){
      ierr = fread(&r[0], size_real, 3, pf);
      assert(ierr == 3);
      if(k == 1){ /*ONLY KEEP DARK MATTER PARTICLES*/
        Q[*ind+pc].Pos[0] = r[0]*POSFACTOR;
        Q[*ind+pc].Pos[1] = r[1]*POSFACTOR;
        Q[*ind+pc].Pos[2] = r[2]*POSFACTOR;
            
        if(Q[*ind+pc].Pos[0] > pmax[0]) pmax[0] = Q[*ind+pc].Pos[0];
        if(Q[*ind+pc].Pos[0] < pmin[0]) pmin[0] = Q[*ind+pc].Pos[0];
        if(Q[*ind+pc].Pos[1] > pmax[1]) pmax[1] = Q[*ind+pc].Pos[1];
        if(Q[*ind+pc].Pos[1] < pmin[1]) pmin[1] = Q[*ind+pc].Pos[1];
        if(Q[*ind+pc].Pos[2] > pmax[2]) pmax[2] = Q[*ind+pc].Pos[2];
        if(Q[*ind+pc].Pos[2] < pmin[2]) pmin[2] = Q[*ind+pc].Pos[2];

        pc++;
      }
    }
  }
  ierr = fread(&d2, sizeof(d2), 1, pf);
  assert(ierr == 1);
  assert(d1==d2);

  ierr = fread(&d1, sizeof(d1), 1, pf);
  assert(ierr == 1);
#ifdef STORE_VELOCITIES
  for(k = 0, pc = 0; k < 6; k++){
    for(n = 0; n < header.npart[k]; n++){
      ierr = fread(&v[0], size_real, 3, pf);
      assert(ierr == 3);
      if(k == 1){ /*ONLY KEEP DARK MATTER PARTICLES*/
        Q[*ind+pc].Vel[0] = v[0]*VELFACTOR;
        Q[*ind+pc].Vel[1] = v[1]*VELFACTOR;
        Q[*ind+pc].Vel[2] = v[2]*VELFACTOR;
        pc++;
      }
    }
  }
#else
  fseek(pf,d1,SEEK_CUR);
#endif
  ierr = fread(&d2, sizeof(d2), 1, pf);
  assert(ierr == 1);
  assert(d1==d2);

  ierr = fread(&d1, sizeof(d1), 1, pf);
  assert(ierr == 1);
#ifdef STORE_IDS
  for(k = 0, pc = 0; k < 6; k++){
    for(n = 0; n < header.npart[k]; n++){
      ierr = fread(&id, size_int, 1, pf);
      assert(ierr == 1);
      if(k == 1){ /*ONLY KEEP DARK MATTER PARTICLES*/
        Q[*ind+pc].id = id;
        pc++;
      }
    }
  }
#else
  fseek(pf,d1,SEEK_CUR);
#endif
  ierr = fread(&d2, sizeof(d2), 1, pf);
  assert(ierr == 1);
  assert(d1==d2);

  *ind += pc;
  
  fclose(pf);
}

void change_positions(my_int n){
  my_int ip;

  RED("Inicio Change Positions\n");

  for(int idim = 0; idim < 3; idim++){
    pmin[idim] = 1.E26; 
    pmax[idim] = -1.E26;
  }

  for(ip = 0; ip < n; ip++){
    if(P[ip].Pos[0] > pmax[0]) pmax[0] = P[ip].Pos[0];
    if(P[ip].Pos[0] < pmin[0]) pmin[0] = P[ip].Pos[0];
    if(P[ip].Pos[1] > pmax[1]) pmax[1] = P[ip].Pos[1];
    if(P[ip].Pos[1] < pmin[1]) pmin[1] = P[ip].Pos[1];
    if(P[ip].Pos[2] > pmax[2]) pmax[2] = P[ip].Pos[2];
    if(P[ip].Pos[2] < pmin[2]) pmin[2] = P[ip].Pos[2];
  }

  printf("xmin %.1f xmax %.1f\n",pmin[0],pmax[0]);
  printf("ymin %.1f ymax %.1f\n",pmin[1],pmax[1]);
  printf("zmin %.1f zmax %.1f\n",pmin[2],pmax[2]);

  for(ip = 0; ip < n; ip++)
    for(int idim = 0; idim < 3; idim++)
      P[ip].Pos[idim] -= pmin[idim];

  cp.lbox = pmax[0] - pmin[0];
  for(int idim = 1; idim < 3; idim++)
    if(cp.lbox < (pmax[idim] - pmin[idim])) cp.lbox = (pmax[idim] - pmin[idim]);

  cp.lbox *= 1.001;
  fprintf(stdout,"Changing cp.lbox %f....\n",cp.lbox);
  GREEN("Fin Change Positions\n");
}

void re_change_positions(my_int n, struct particle_data *Q){
  my_int ip;
  for(ip = 0; ip < n; ip++)
    for(int idim = 0; idim < 3; idim++)
      Q[ip].Pos[idim] += pmin[idim];
}

void select_particles(void){
  my_int i,j;

  fprintf(stdout,"Selecting particles within a subbox\n"); fflush(stdout);
  
  for(j = 0, i = 0; i < cp.npart; i++){
    if((P[i].Pos[0] - box.cen[0]) >  cp.lbox/2.0) P[i].Pos[0] -= (my_real)cp.lbox;
    if((P[i].Pos[0] - box.cen[0]) < -cp.lbox/2.0) P[i].Pos[0] += (my_real)cp.lbox;
    if((P[i].Pos[1] - box.cen[1]) >  cp.lbox/2.0) P[i].Pos[1] -= (my_real)cp.lbox;
    if((P[i].Pos[1] - box.cen[1]) < -cp.lbox/2.0) P[i].Pos[1] += (my_real)cp.lbox;
    if((P[i].Pos[2] - box.cen[2]) >  cp.lbox/2.0) P[i].Pos[2] -= (my_real)cp.lbox;
    if((P[i].Pos[2] - box.cen[2]) < -cp.lbox/2.0) P[i].Pos[2] += (my_real)cp.lbox;

    ///////////////////////////////////////////////////
    if(P[i].Pos[0] < (box.cen[0]-(box.lado+box.franja)) || 
       P[i].Pos[0] > (box.cen[0]+(box.lado+box.franja)))continue;
    if(P[i].Pos[1] < (box.cen[1]-(box.lado+box.franja)) || 
       P[i].Pos[1] > (box.cen[1]+(box.lado+box.franja)))continue;
    if(P[i].Pos[2] < (box.cen[2]-(box.lado+box.franja)) || 
       P[i].Pos[2] > (box.cen[2]+(box.lado+box.franja)))continue;
    ///////////////////////////////////////////////////

    j++;
  }
  cp.npart = j;
  P = (struct particle_data *) realloc(P,cp.npart*sizeof(struct particle_data));
}
