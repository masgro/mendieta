#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "variables.h"
#include "propiedades.h"

#define NPARTMIN_LOCAL 20

struct hijost
{
  int nh;
  int *h;
} *hijos;

struct subst
{
  int np;
  my_real vmax;
  my_int *id;
  my_int *indx;
} *lala;

int main(int argc, void **argv)
{
  char filename[200];
  FILE *pfout,*pfin,*pfout1;
  int j, k;
  unsigned np, npt, n;
  unsigned int ng, nobj;
  int id;
  my_int idt;
  my_int zero = 0;
  my_int indx;
  my_int i;
  int cuenta,cuenta1,tot;
	char path[200];
  char filepegado[200],filesubin[200];
  char fileprop[200],filehijos[200];

  init_variables();

	sprintf(path,"/home/marioagustin/trabajo/idensub/prueba/aquarius/");

  sprintf(filesubin,"%s/idensub/idensub.bin",path);
  sprintf(filehijos,"%s/idensub/propiedades/hijos.ascii",path);
  sprintf(filepegado,"%s/idensub/idensub_pegado.bin",path);
  sprintf(fileprop,"%s/idensub/propiedades/propiedades.ascii",path);

  sprintf(filename,"%s/AqA_haloids_4_SR.dat",path);
  pfout = fopen(filename,"w");
  fprintf(pfout,"# Aquarius subhalo comparison project\n");
  fprintf(pfout,"# Level: 4\n");
  fprintf(pfout,"# Author: Sgró, M. A. & Ruiz, A. N.\n");
  fprintf(pfout,"# Halo finder: Mendieta\n");
  fprintf(pfout,"# Date: 05-15-2012 00:00\n");

  pfout1 = fopen(filepegado,"w");

  pfin = fopen(filesubin,"r");
  fread(&ng,sizeof(unsigned int),1,pfin);
  fread(&nobj,sizeof(unsigned int),1,pfin);

  fwrite(&ng,sizeof(unsigned int),1,pfout1);
  fwrite(&nobj,sizeof(unsigned int),1,pfout1);

  printf("ng: %u nobj %u\n",ng,nobj);

  lala = (struct subst *) malloc(ng*sizeof(struct subst));

  for(i=0;i<ng;i++)
  {
    np = 0; npt = 0;
    fread(&np,sizeof(unsigned int),1,pfin);

    lala[i].np = np;

    lala[i].id = (my_int *) malloc(np*sizeof(my_int));
    lala[i].indx = (my_int *) malloc(np*sizeof(my_int));

    for(j=0;j<np;j++)
    {
      //fread(&id,size_int,1,pfin);
      fread(&indx,size_int,1,pfin);
      assert(idt>=zero);
      //lala[i].id[j] = idt;
      lala[i].indx[j] = indx;
    }
    fread(&npt,sizeof(unsigned int),1,pfin);
    assert(npt==np);
  }
  fclose(pfin);
	printf("termino de leer %s\n",filename);

  hijos = (struct hijost *) malloc(ng*sizeof(struct hijost));

  pfin = fopen(filehijos,"r");
  int nh,h;
  for(i=0;i<ng;i++)
  {
    fscanf(pfin,"%d %d\n",&id,&nh);
    assert(i==id);
    hijos[id].nh = nh;
    hijos[id].h = (int *) calloc(nh,sizeof(int));
    for(j=0;j<nh;j++)
    {
      fscanf(pfin,"%d\n",&h);
      hijos[id].h[j] = h;
    }
  }
  fclose(pfin);
  printf("termino lecturas\n");

  //pfin = fopen(fileprop,"r");
  //int idum;
  //my_real fdum;
  //for (i=0; i<ng; i++)
  //{
  //  fscanf(pfin,"%d %lf %lf %lf %lf %lf %lf ",
	//	     &idum,&fdum,&fdum,&fdum,&fdum,&fdum,&fdum);
  //  fscanf(pfin,"%lf %lf %lf ",&fdum,&fdum,&fdum);
  //  fscanf(pfin,"%lf %lf %lf %lf %lf %lf %lf ",
  //       &fdum,&fdum,&fdum,&fdum,&fdum,&fdum,&fdum);
  //  fscanf(pfin,"%lf %lf %lf %lf %lf %lf %lf ",
  //       &fdum,&fdum,&fdum,&fdum,&fdum,&fdum,&sub[i].vmax);
  //  fscanf(pfin,"%lf %lf %lf %lf %lf %lf %lf %lf\n",
  //       &fdum,&fdum,&fdum,&fdum,&fdum,&fdum,&fdum,&fdum);

  //  assert(idum==sub[i].np);
  //}

  //fclose(pfin);

  n = 0;
  tot = 0;
  for(i=0;i<ng;i++)
  {
    if(lala[i].np < NPARTMIN_LOCAL)continue;

    cuenta = lala[i].np;

    for(k=0;k<hijos[i].nh;k++)
    {
      cuenta += lala[hijos[i].h[k]].np;
    }

    fprintf(pfout,"%03u %07d\n",n,cuenta);
    fwrite(&cuenta,sizeof(int),1,pfout1);
    tot += cuenta;

    cuenta1=0;
    for(j=0;j<lala[i].np;j++)
    {
      fprintf(pfout,"%d\n",lala[i].id[j]);
      //fwrite(&lala[i].id[j],size_int,1,pfout1);
      fwrite(&lala[i].indx[j],size_int,1,pfout1);
      cuenta1++;
    }

    assert(cuenta1==lala[i].np);

    for(k=0;k<hijos[i].nh;k++)
    {
      for(j=0;j<lala[hijos[i].h[k]].np;j++)
      {
        fprintf(pfout,"%d\n",lala[hijos[i].h[k]].id[j]);
        //fwrite(&lala[hijos[i].h[k]].id[j],size_int,1,pfout1);
        fwrite(&lala[hijos[i].h[k]].indx[j],size_int,1,pfout1);
        cuenta1++;
      }
    }

		assert(cuenta1==cuenta);
    fprintf(pfout,"# %03u %07d\n",n,cuenta);
    fwrite(&cuenta,sizeof(int),1,pfout1);
    n++;
  }

  rewind(pfout1);
  fwrite(&n,sizeof(unsigned int),1,pfout1);
  fwrite(&tot,sizeof(unsigned int),1,pfout1);

  printf("ng: %u nobj %d\n",n,tot);

  fclose(pfout);
  fclose(pfout1);

  exit(0);

}
