#include <stdlib.h>
#include <stdio.h>
#include "variables.h"
#include "propiedades.h"

#define NPARTMIN_LOCAL 20

int main(int argc, void **argv)
{
  char filename[200];
  char propfile[200];
  FILE *pfout,*pfin;
  //float Mpart = 0.000229439; /*A-5*/
  float Mpart = 2.86817e-05; /*A-4*/
  //float Mpart = 3.58481e-06; /*A-3*/
  int i, id, ng;
	char path[200];

	sprintf(path,"/home/marioagustin/trabajo/idensub/prueba/aquarius/");
  sprintf(propfile,"%s/idensub/propiedades/propiedades.bin",path);

  init_variables();
  printf("size_real %zu\n",size_real);
  printf("size_int %zu\n",size_int);

  sprintf(filename,"%s/AqA_halodata_4_SR_bis.dat",path);
  pfout = fopen(filename,"w");
  if(pfout==NULL)printf("noooo\n");
  fprintf(pfout,"# Aquarius subhalo comparison project\n");
  fprintf(pfout,"# Level: 4\n");
  fprintf(pfout,"# Author: Sgró, M. A. & Ruiz, A. N.\n");
  fprintf(pfout,"# Halo finder: Mendieta\n");
  fprintf(pfout,"# Date: 01-05-2012 11:11\n");

  pfin = fopen(propfile,"r");
  fread(&ng,sizeof(int),1,pfin);
  printf("ng %d\n",ng);
  id = 0;
  for(i=0;i<ng;i++)
  {
    read_properties(pfin);
    if(npart < NPARTMIN_LOCAL)continue;
    pcm[0] = mostbound[0];
    pcm[1] = mostbound[1];
    pcm[2] = mostbound[2];
    m200 = (float)npart*Mpart*1E10;
    r200 = r200;
    fprintf(pfout,"%04d %13.6f %13.6f %13.6f %11.6f %11.6f %11.6f %.8E %11.6f \n",
             id,pcm[0],pcm[1],pcm[2],vcm[0],vcm[1],vcm[2],m200,r200);
    id++;
  }
  fclose(pfin);
  fclose(pfout);
  exit(0);

}
