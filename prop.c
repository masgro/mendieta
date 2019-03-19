#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "variables.h"
#include "cosmoparam.h"
#include "leesnap.h"
#include "timer.h"
#include "propiedades.old.h"
#include "deltas.h"

void read_idenfof(char *file);
char ffof[200],bfof[200];

struct subst
{
  double pcm[3];
  double r200;
} *sub_hijos;

int main(int argc, char **argv)
{
	int    i;
  double start,end;
  char   filename[200],fout[200];
  FILE   *pfout, *pfin;
  FILE   *asco;
  struct propiedades_st Prop;

	TIMER(start);

  sprintf(fout,"propiedades.bin");

  init_variables(argc,argv);

  /* Lee archivos de la simulacion */
	read_gadget();

  /* Lee archivo de la identificacion FoF */
  read_idenfof(fof_file);

  /* Si no encuentra grupos, crea archivo vacio */
	if(n_grupos_fof == 0)
	{
  	printf("No se encontraron grupos fof. Creando archivo vacio y saliendo.\n");
  	pfout = fopen(fout,"w");
  	printf("Abre archivo de salida: %s\n",fout);
  	fwrite(&n_grupos_fof,sizeof(int),1,pfout);
		fclose(pfout);
		return(EXIT_SUCCESS);
	}

  #ifdef COMPUTE_EP
  /* Calcula energia potencial de las particulas */
  compute_potential_energy();
  #endif

  pfout = fopen(fout,"w");
  fprintf(stdout,"Abre archivo de salida: %s\n",fout);

  #ifdef IDENSUB
  sub_hijos = (struct subst *) calloc((n_grupos_fof + 1),sizeof(struct subst));
  #endif

  asco = fopen("propiedades.ascii","w");

  fwrite(&n_grupos_fof,sizeof(int),1,pfout);
	for(i = 1; i <= n_grupos_fof; i++)
	{
		fprintf(stdout,"Calculando prop. grupo %6d de %6d, nmem %d\n",
                                           i,n_grupos_fof,fof[i].np);
		Prop = propiedades(P,i);

    //write_properties(pfout,Prop);

    fprintf(asco,"%d %f %f %f %f %f %f\n",
		     Prop.npart,Prop.pcm[0],Prop.pcm[1],Prop.pcm[2],Prop.vcm[0],Prop.vcm[1],Prop.vcm[2]);
//    fprintf(asco,"%f %f %f %f %f %f %f ",
//         sig[0],sig[1],sig[2],L[0],L[1],L[2],lambda);
//    fprintf(asco,"%f %f %f %f %f %f %f ",
//         m200,r200,v200,mvir,rvir,vvir,vmax);
//    fprintf(asco,"%f %f %f %f %f %f %f %f\n",
//         Ep,Ec,aa,bb,cc,aa_vel,bb_vel,cc_vel);
    fflush(asco);

    #ifdef IDENSUB
    sub_hijos[i].pcm[0] = mostbound[0];
    sub_hijos[i].pcm[1] = mostbound[1];
    sub_hijos[i].pcm[2] = mostbound[2];
    sub_hijos[i].r200   = rvir;
    #endif
    
  }

	fclose(pfout);
	fclose(asco);

#ifdef IDENSUB
	sprintf(filename,"%s",fhijos);
  asco = fopen(filename,"w");
  int *hijo,j,k,l;
  double dx,dy,dz,d;

  hijo = (int *) malloc((n_grupos_fof+1)*sizeof(int));

  for(i=0;i<n_grupos_fof;i++)
  {
    for(l=0;l<n_grupos_fof;l++) hijo[l] = -1;
    k=0;
    for(j=i+1;j<n_grupos_fof;j++)
    {
      dx = sub_hijos[i].pcm[0] - sub_hijos[j].pcm[0];
      dy = sub_hijos[i].pcm[1] - sub_hijos[j].pcm[1];
      dz = sub_hijos[i].pcm[2] - sub_hijos[j].pcm[2];
      d = sqrt(dx*dx + dy*dy + dz*dz);
      if(d < sub_hijos[i].r200)
      {
        hijo[k] = j;
        k++;
      }
      
    }
    fprintf(asco,"%d %d\n",i,k);
    for(l=0;l<k;l++)
      fprintf(asco,"%d\n",hijo[l]);
  }
  fclose(asco);
#endif
	TIMER(end);
	printf("\nTiempo Total %lf\n",end-start);

	return(EXIT_SUCCESS);
}

void read_idenfof(void)
{
  char         file[200];
  FILE         *pfin;
  int          i, n, nt;
  unsigned int indx;
  unsigned int  _ng, _npg;
  struct node *item;

  /**** File Grupos FOF ****/
  sprintf(file,"%s%s",bfof,ffof);
  pfin = fopen(file,"r");
  printf("Open FoF groups file: %s\n",file);

  /**** LEE PARTICULAS EN GRUPOS FOF Y CREA LINKED LIST ****/

  fread(&n_grupos_fof,sizeof(int),1,pfin);
  fread(&np_in_fof,sizeof(int),1,pfin);

	fprintf(stdout,"Numero de grupos: %u\n",n_grupos_fof);
	fprintf(stdout,"Numero de particulas en grupos: %u\n",np_in_fof);

	if(n_grupos_fof == 0) exit(EXIT_SUCCESS);

  groups = (struct stuff *) malloc((n_grupos_fof + 1)*sizeof(struct stuff));
  for(i = 0; i < (n_grupos_fof + 1); i++) groups[i].first = NULL;

  n_grupos_fof = 0; np_in_fof = 0;
  do
  {
    fread(&n,sizeof(int),1,pfin);

    if(feof(pfin))break;

    n_grupos_fof++;

    groups[n_grupos_fof].np = n;

    for( i = 0 ; i < n ; i++ )
    {
      fread(&indx,sizeof(int),1,pfin);

      item = (struct node *) malloc(sizeof(struct node));

      item->indx = indx;
      item->next = groups[n_grupos_fof].first;
      groups[n_grupos_fof].first = item;

      np_in_fof++;
    }
    fread(&nt,sizeof(int),1,pfin);

    assert(n == nt);

  }while(1);

  fprintf(stdout,"End reading FoF groups: \n");
	fprintf(stdout,"Numero de grupos: %u\n",n_grupos_fof);
	fprintf(stdout,"Numero de particulas en grupos: %u\n",np_in_fof);
}


