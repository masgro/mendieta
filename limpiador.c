#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "variables.h"
#include "cosmoparam.h"
#include "leesnap.h"
#include "timer.h"
#include "deltas.h"
#include "limpieza.h"
#include "io.h"

void read_idenfof_bis(void);
void fof_groups(void);

int main(int argc, char **argv){
	int    i;
  double start,end;
  char   filename[200];
  char   fout[200];
  FILE   *pfout, *pfin;

	TIMER(start);

  init_variables(argc,argv);

  sprintf(fout,"%s",sub_file);

  /* Lee archivos de la simulacion */
	read_gadget();

  /* Lee archivo de la identificacion FoF */
  read_idenfof_bis();

  /* Si no encuentra grupos, crea archivo vacio */
	if(n_grupos_fof == 0){
  	printf("No se encontraron grupos fof. Creando archivo vacio y saliendo.\n");
  	pfout = fopen(fout,"w");
  	printf("Abre archivo de salida: %s\n",fout);
  	fwrite(&n_grupos_fof,sizeof(int),1,pfout);
		fclose(pfout);
		return(EXIT_SUCCESS);
	}

  fof_groups();

  write_idenfof(fout);

	TIMER(end);
	printf("\nTiempo Total %lf\n",end-start);

	return(EXIT_SUCCESS);
}

void read_idenfof_bis(void){
  FILE     *pfin;
  int      n, nt;
  my_int i,j;
  my_int indx;

  /**** File Grupos FOF ****/
  pfin = fopen(fof_file,"r");
  printf("Open FoF groups file: %s\n",fof_file);

  /**** LEE PARTICULAS EN GRUPOS FOF Y CREA LINKED LIST ****/

  fread(&n_grupos_fof,sizeof(int),1,pfin);
  fread(&np_in_fof,sizeof(int),1,pfin);

	fprintf(stdout,"Numero de grupos: %u\n",n_grupos_fof);
	fprintf(stdout,"Numero de particulas en grupos: %u\n",np_in_fof);

	if(n_grupos_fof == 0) exit(EXIT_SUCCESS);

  n_grupos_fof += 1;
  Temp.head   = (my_int *) malloc(n_grupos_fof*sizeof(my_int));
  Temp.npgrup = (my_int *) calloc(n_grupos_fof,sizeof(my_int));
  Temp.ll     = (my_int *) malloc(cp.npart*sizeof(my_int));

  np_in_fof = 0;
  for(my_int g = 1; g < n_grupos_fof; g++){
    fread(&n,sizeof(my_int),1,pfin);

    Temp.npgrup[g] = n;

    for(j = 0; j < n ; j++){
      fread(&indx,sizeof(my_int),1,pfin);
      P[indx].gr = g;

      Temp.head[g] = indx;
      np_in_fof++;
    }
    fread(&nt,sizeof(my_int),1,pfin);
    assert(n == nt);
  }

  for(i = 0; i < cp.npart; i++){
    my_int g = P[i].gr;
    Temp.ll[i] = Temp.head[g];
    Temp.head[g] = i;
  }

  fprintf(stdout,"End reading FoF groups: \n");
	fprintf(stdout,"Numero de grupos: %u\n",n_grupos_fof);
	fprintf(stdout,"Numero de particulas en grupos: %u\n",np_in_fof);
}

void fof_groups(void){
  my_int i;

  /* Allocatea memoria */
  fof = (struct grupos *) malloc(n_grupos_fof*sizeof(struct grupos));
  assert(fof != NULL);

  for(i = 1; i < n_grupos_fof; i++){
    printf("%d\n",i);
    limpieza_new(i,0);
  }

  for(i = 0; i < n_grupos_fof; i++){
    fof[i].np     =  0;
  }

  for(i = 0; i < cp.npart; i++){
    fof[P[i].gr].llirst = i;
    fof[P[i].gr].np++;
  }

  for(i = 0; i < cp.npart; i++){
    P[i].llfof = fof[P[i].gr].llirst;
    fof[P[i].gr].llirst = i;
  }
}
