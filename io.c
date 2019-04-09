#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdbool.h>
#include "variables.h"
#include "cosmoparam.h"
#include "leesnap.h"
#include "io.h"
#include "limpieza.h"

#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector_long.h>
#include <gsl/gsl_sort_vector_long.h>

#ifdef READIDENFOF
void read_idenfof(char *file){
  FILE *pfin;
  unsigned int n,nt,_ng,_npg,count,nmax;
  my_int idt,imax,i,j;

  for(i = 0; i < cp.npart; i++) P[i].fof = 0;

  /**** File FOF Groups ****/
  pfin = fopen(file,"r");
  printf("Open FoF groups file: %s\n",file);

  #ifdef NHALO
  printf("Warning!.. running MENDIETA with NHALO option\n");
  printf("MENDIETA will run only over fof halo %d\n",NHALO);
  #endif

  /**** LEE PARTICULAS EN GRUPOS FOF Y CREA LINKED LIST ****/
  n_grupos_fof  = 0;
  np_in_fof = 0;
  nmax  = 0;
  imax  = 0;
  count = 0;

  fread(&_ng,sizeof(unsigned int),1,pfin);
  fread(&_npg,sizeof(unsigned int),1,pfin);
  while(1){
    fread(&n,sizeof(unsigned int),1,pfin);
    if(feof(pfin))break;
    count++;

    /***********************************************/
    /*ESTO ES PARA SELECCIONAR UN SOLO HALO */
    #ifdef NHALO
    if( count != NHALO ){
      for( i = 0 ; i < n ; i++ ){
        fread(&idt,sizeof(my_int),1,pfin);
      }
      fread(&nt,sizeof(unsigned int),1,pfin);
      assert(n == nt);
      continue;
    }
    #endif
    /***********************************************/
    n_grupos_fof++;

    /*Selecciona Grupo FoF mas masivo*/
    if(n > nmax){nmax = n; imax = n_grupos_fof;}

    for(i = 0; i < n; i++){
      fread(&idt,sizeof(my_int),1,pfin);

      P[idt].fof = n_grupos_fof;

      np_in_fof++;
    }
    fread(&nt,sizeof(unsigned int),1,pfin);
    assert(n == nt);
  }

  fclose(pfin);

  #ifndef NHALO
  assert(np_in_fof==_npg);
  #endif

  printf("End reading FoF groups: \n");
  printf("FoF groups = %u \n",n_grupos_fof);
  printf("Particles in biggest FoF group = %d \n",nmax);
  printf("Biggest FoF group ID = %d \n",imax);

  for(i = 0; i < 3; i++){
    pmin[i] =  1.E26; pmax[i] = -1.E26;
  }

  /*Reallocatea para liberar memoria*/
  j = 0;
  for(i = 0; i < cp.npart; i++){

    if(P[i].fof == 0)
      continue;

    P[j].Pos[0] = P[i].Pos[0];
    P[j].Pos[1] = P[i].Pos[1];
    P[j].Pos[2] = P[i].Pos[2];

    #ifdef STORE_VELOCITIES
    P[j].Vel[0] = P[i].Vel[0];
    P[j].Vel[1] = P[i].Vel[1];
    P[j].Vel[2] = P[i].Vel[2];
    #endif

    if(P[j].Pos[0] > pmax[0]) pmax[0] = P[j].Pos[0];
    if(P[j].Pos[0] < pmin[0]) pmin[0] = P[j].Pos[0];
    if(P[j].Pos[1] > pmax[1]) pmax[1] = P[j].Pos[1];
    if(P[j].Pos[1] < pmin[1]) pmin[1] = P[j].Pos[1];
    if(P[j].Pos[2] > pmax[2]) pmax[2] = P[j].Pos[2];
    if(P[j].Pos[2] < pmin[2]) pmin[2] = P[j].Pos[2];

    P[j].fof    = P[i].fof;
    #ifdef IDENSUB
    P[j].sub    = P[i].fof;
    #endif
    //P[j].indx   = i;
    j++;
  }
  assert(j == np_in_fof);
  fprintf(stdout,"Changing cp.npart %d -> %d....\n",cp.npart,j); 
  cp.npart = j;

  P = (struct particle_data *) realloc(P,cp.npart*sizeof(struct particle_data));

  /* Allocatea memoria */
  n_grupos_fof += 1;
  fof = (struct grupos *) calloc(n_grupos_fof,sizeof(struct grupos));
  #ifdef IDENSUB
  n_grupos_sub = n_grupos_fof;
  sub = (struct grupos *) calloc(n_grupos_sub,sizeof(struct grupos));
  #endif

  for(i = 0; i < n_grupos_fof; i++){
    fof[i].llirst = -1;
    #ifdef IDENSUB
    sub[i].llirst = -1;
    #endif
  }

  for(i = 0; i < cp.npart; i++){
    j = P[i].fof;

    P[i].llfof = fof[j].llirst;
    fof[j].llirst = i;
    fof[j].np++;

    #ifdef IDENSUB
    P[i].sub = P[i].fof;
    P[i].llsub = sub[j].llirst;
    sub[j].llirst = i;
    sub[j].np++;
    #endif
  }
}
#endif

void write_idenfof(char *filename){
  FILE *pf;
  my_int i, l, nobj, ng;
  #ifdef SUBBOXES
  double xcm,ycm,zcm;
  #endif

	gsl_permutation *index;
  gsl_vector_long *group;

  pf = fopen(filename,"w");
  if(pf == NULL){
    fprintf(stderr,"No se pudo abrir archivo: %s\n",filename);
    exit(EXIT_FAILURE);
  }
  printf("Writing output file %s\n",filename);

#ifdef SUBBOXES
  for(i = 1; i < n_grupos_fof; i++){
    if(fof[i].np == 0) continue;

    xcm = 0.0;
    ycm = 0.0;
    zcm = 0.0;

    l = fof[ng].llirst;
    do{
      xcm += (double)P[l].Pos[0];
      ycm += (double)P[l].Pos[1];
      zcm += (double)P[l].Pos[2];

      l = P[l].llfof;
    }while(l != fof[ng].llirst);

    xcm /= (double)fof[i].np;
    ycm /= (double)fof[i].np;
    zcm /= (double)fof[i].np;

    if(xcm >= box.max[0]){fof[i].np=0;continue;}
    if(xcm <  box.min[0]){fof[i].np=0;continue;}
    if(ycm >= box.max[1]){fof[i].np=0;continue;}
    if(ycm <  box.min[1]){fof[i].np=0;continue;}
    if(zcm >= box.max[2]){fof[i].np=0;continue;}
    if(zcm <  box.min[2]){fof[i].np=0;continue;}
  }
#endif

  group = gsl_vector_long_alloc(n_grupos_fof);
  index = gsl_permutation_alloc(n_grupos_fof);

  nobj = 0; ng = 0;   
  for(i = 0; i < n_grupos_fof; i++){
    if(fof[i].np == 0 || i == 0){
      gsl_vector_long_set(group,i,0);
    }else{
      gsl_vector_long_set(group,i,fof[i].np);
      nobj += fof[i].np;
      ng++;
    }
  }

  gsl_sort_vector_long_index (index,group);

  fprintf(stdout,"Numero Final de grupos FoF: %d\n",ng);
  fprintf(stdout,"Numero de particulas en grupos FoF: %d\n",nobj);

  fwrite(&ng,sizeof(unsigned int),1,pf);
  fwrite(&nobj,sizeof(unsigned int),1,pf);
  
  for(i = n_grupos_fof - 1; i > 0; i--){
    ng = gsl_permutation_get(index,i); 
    if(fof[ng].np == 0)continue;

    fwrite(&fof[ng].np,sizeof(unsigned int),1,pf);

    #ifdef DEBUG
    int k = 0;
    #endif
    l = fof[ng].llirst;
    while(l != GROUND){
      fwrite(&l,sizeof(my_int),1,pf);

      //reasgina grupos fof a la particula l para que los grupos
      //este ordenados de mayor a menor
      //P[l].fof = n_grupos_fof - i ;

      #ifdef DEBUG
      k++;
      #endif

      l = P[l].llfof;
      //l = Temp.ll[l];
    }
    fwrite(&fof[ng].np,sizeof(unsigned int),1,pf);
    #ifdef DEBUG
    assert(k == fof[ng].np);
    #endif
  }
  fclose(pf);

  gsl_vector_long_free(group);
  gsl_permutation_free(index);
}

#ifdef IDENSUB
void write_idensub(char *sub_file){
  char filename[200];
  FILE *pfsub, *pffof, *pfascii;
  my_int nobj, ng;
  my_int i,l;
  int j;
  size_t it;

	gsl_permutation *index;
  gsl_vector_long *group;

  sprintf(filename,"%s",sub_file);
  pfsub = fopen(filename,"w");
  if(pfsub == NULL){
    fprintf(stderr,"No se pudo abrir archivo: %s\n",filename);
    exit(EXIT_FAILURE);
  }
  printf("Writing output file %s\n",filename);

  sprintf(filename,"%s.fof",sub_file);
  pffof = fopen(filename,"w");
  if(pffof == NULL){
    fprintf(stderr,"No se pudo abrir archivo: %s\n",filename);
    exit(EXIT_FAILURE);
  }
  printf("Writing output file %s\n",filename);

  pfascii = fopen("sub.ascii","w");

  nobj = 0; ng = 0;   
  for(i = 1 ;i < n_grupos_sub ;i++){
    if(sub[i].np == 0)continue;
    nobj += sub[i].np;
    ng++;
  }

  fprintf(stdout,"Numero Final de Subgrupos: %d\n",ng);
  fprintf(stdout,"Numero de particulas en subgrupos: %d\n",nobj);

  fwrite(&ng,sizeof(unsigned int),1,pffof);

  fwrite(&ng,sizeof(unsigned int),1,pfsub);
  fwrite(&nobj,sizeof(unsigned int),1,pfsub);
  
  Temp.head   = (my_int *) malloc(n_grupos_fof*sizeof(my_int));
  Temp.npgrup = (my_int *) calloc(n_grupos_fof,sizeof(my_int));
  Temp.ll     = (my_int *) malloc(n_grupos_sub*sizeof(my_int));
  Temp.grup   = (my_int *) calloc(n_grupos_sub,sizeof(my_int));

  Temp.grup[0] = 0;
  for(i = 1 ; i < n_grupos_sub ; i++){
    if(sub[i].np == 0) continue;
    j = sub[i].llirst;
    Temp.grup[i] = P[j].fof;
  }

  for(i = 0; i < n_grupos_fof; i++) Temp.head[i] = GROUND;

  for(i = 0; i < n_grupos_sub; i++){
    ng = Temp.grup[i];
    Temp.ll[i] = Temp.head[ng];
    Temp.head[ng] = i;
    Temp.npgrup[ng]++;
  }

  for(i = 1; i < n_grupos_fof; i++){
    if(Temp.npgrup[i] == 0) continue;

    Temp.in  = (my_int *) malloc(Temp.npgrup[i]*sizeof(my_int));

    group = gsl_vector_long_alloc(Temp.npgrup[i]);
    index = gsl_permutation_alloc(Temp.npgrup[i]);

    l = Temp.head[i]; j = 0;
    while(l != GROUND){
      #ifdef DEBUG
      assert(j < Temp.npgrup[i]);
      #endif
      gsl_vector_long_set(group,j,sub[l].np);
      Temp.in[j] = l;

      j++;

      l = Temp.ll[l];
    }

    #ifdef DEBUG
    assert(j == Temp.npgrup[i]);
    #endif

    gsl_sort_vector_long_index(index,group);

    for(j = Temp.npgrup[i] - 1; j >= 0; j--){
      it = gsl_permutation_get(index,j);
      ng = Temp.in[it];

      if(sub[ng].np == 0)continue;

      fwrite(&sub[ng].np,sizeof(unsigned int),1,pfsub);
      fwrite(&P[sub[ng].llirst].fof,sizeof(unsigned int),1,pffof);

      fprintf(pfascii,"%d %d %d\n",i,ng,sub[ng].np);

      #ifdef DEBUG
      int k = 0;
      #endif
      l = sub[ng].llirst;
      while(l != GROUND){
        #ifdef READIDENFOF
        fwrite(&P[l].indx,sizeof(int),1,pfsub);
        #else
        fwrite(&l,sizeof(int),1,pfsub);
        #endif
        l = P[l].llsub;
        #ifdef DEBUG
        k++;
        #endif
      }
      fwrite(&sub[ng].np,sizeof(unsigned int),1,pfsub);
      
      #ifdef DEBUG
      assert(k == sub[ng].np);
      #endif
    }

    gsl_permutation_free(index);
    gsl_vector_long_free(group);

    free(Temp.in);
  }
  fclose(pfsub);
  fclose(pffof);
  fclose(pfascii);

  free(Temp.head);
  free(Temp.npgrup);
  free(Temp.ll);
  free(Temp.grup);
}
#endif

#ifdef GETPOSITIONS
void get_positions(int ifrac){
  char filename[200];
  FILE *pfout;
  my_int i, l;
  sprintf(filename,"pos.%02d",ifrac);
  pfout = fopen(filename,"w");
  for(i = 0; i < n_grupos_sub; i++){
    l = sub[i].llirst;
    while(l != GROUND){
      fprintf(pfout,"%f %f %f %d\n",P[l].Pos[0],P[l].Pos[1],P[l].Pos[2],i);
      l = P[l].llsub;
    }
  }
  fclose(pfout);
}
#endif

#ifdef GETINDEX
void get_index(int ifrac){
  char filename[200];
  FILE *pfout;
  my_int i, l;
  sprintf(filename,"index.%02d",ifrac);
  pfout = fopen(filename,"w");
  for(i = 0; i < n_grupos_fof; i++){
    l = fof[i].llirst;
    while(l != GROUND){
      fprintf(pfout,"%d %d\n",l,i);
      l = P[l].llfof;
    }
  }
  fclose(pfout);
}
#endif
