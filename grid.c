#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "variables.h"
#include "grid.h"
#include "cosmoparam.h"

/*
i = grid.llirst[igrid];
while(i != -1)
{
  //hace lo que quieras

  i = grid.ll[i];
}
*/

void grid_init(void){
  unsigned long nalloc = grid.ngrid*grid.ngrid*grid.ngrid;
  size_t memsize;

  memsize = (grid.nobj+nalloc)*sizeof(int)/1024.0/1024.0/1024.0;
  printf("allocating %.5zu gb\n",memsize);

	grid.llirst	= (my_int *) malloc(nalloc*sizeof(my_int));
  assert(grid.llirst != NULL);

	grid.ll = (my_int *) malloc(grid.nobj*sizeof(my_int));
  assert(grid.ll != NULL);
}

void grid_build(void){
  my_int i;
  unsigned long ix, iy, iz;
  double fac;
	unsigned long ibox;

  fac = (double)grid.ngrid/(double)cp.lbox ;
	printf("Building Grid..... Ngrid = %lu\n",(unsigned long)grid.ngrid);

  for( i = 0 ; i < grid.nobj ; i++ ){
    ix = (unsigned long)((double)P[i].Pos[0]*fac);
    iy = (unsigned long)((double)P[i].Pos[1]*fac);
    iz = (unsigned long)((double)P[i].Pos[2]*fac);

		ibox = (ix * grid.ngrid + iy) * grid.ngrid + iz;

    #ifdef DEBUG
    assert(ibox >= 0L);
    if(ibox >= grid.ngrid*grid.ngrid*grid.ngrid)
      printf("%d %lu %lu %lu %lu %f %f %f\n",i,ix,iy,iz,ibox,P[i].Pos[0],P[i].Pos[1],P[i].Pos[2]);
    assert(ibox < grid.ngrid*grid.ngrid*grid.ngrid);
    #endif

    grid.llirst[ibox] = i;
  }

  for( i = 0 ; i < grid.nobj ; i++ ){
    ix = (unsigned long)((double)P[i].Pos[0]*fac);
    iy = (unsigned long)((double)P[i].Pos[1]*fac);
    iz = (unsigned long)((double)P[i].Pos[2]*fac);

		ibox = (ix * grid.ngrid + iy) * grid.ngrid + iz;

    grid.ll[i] = grid.llirst[ibox];
    grid.llirst[ibox] = i;
  }
	printf("End Building Grid.....\n");

}

void grid_free(void){
  if(grid.ll!=NULL)free(grid.ll);
  if(grid.llirst!=NULL)free(grid.llirst);
}
