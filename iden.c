#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <omp.h>

#include "variables.h"
#include "cosmoparam.h"
#include "grid.h"
#include "iden.h"

//#define FOF_MARIO


struct global_head_st{
  item *head;
  struct global_head_st *global_next;
};

void identification(void){
  int       i,nthreads,tid;
  bool      *test;
  int *Padre, *Size;
  int *Start, *End;
  int ngrid_old = grid.ngrid;

  grid.ngrid = (int)(cp.lbox/iden.r0);
  test = (bool *) malloc(cp.npart*sizeof(bool));
  Padre = (int *) malloc(cp.npart*sizeof(int));
  Size = (int *)  malloc(cp.npart*sizeof(int));

  if(grid.ngrid > NGRIDMAX){
    fprintf(stdout,"Using NGRIDMAX = %d\n",NGRIDMAX);
    grid.ngrid = NGRIDMAX;
  }

  grid.nobj = iden.nobj;

  if(ngrid_old != grid.ngrid)
  {
    grid_free();
    grid_init();
    grid_build();
  }

  iden.r0 = iden.r0*iden.r0;

  #ifdef NTHREADS
  omp_set_dynamic(0);
  omp_set_num_threads(NTHREADS);
  #endif
  
  for(i=0;i<cp.npart;i++)
  { 
    Padre[i] = i;
    Size[i] = 0;
    test[i] = false;
    P[i].gr = 0;
  }

  fprintf(stdout,"Comienza identificacion.....\n");
  iden.ngrupos = 0;

  #pragma omp parallel default(none) private(tid,i) \
  shared(P,iden,nthreads,Start,End,cp,Padre,stdout,test)   
  {

    tid = omp_get_thread_num();

    if(tid == 0)
    {
      nthreads = omp_get_num_threads();
      printf("\nRunning on %d threads\n",nthreads);
      nthreads = 1024;
 
      Start = (int *) malloc(nthreads*sizeof(int));
      End = (int *) malloc(nthreads*sizeof(int));
 
      for(i=0;i<nthreads;i++)
      {
        Start[i] = i==0 ? 0 : i*floor((float)cp.npart/nthreads);
        End[i] = i==nthreads-1 ? cp.npart : (i+1)*floor((float)cp.npart/nthreads);
      }
    }

    #pragma omp barrier

    //#pragma omp for schedule(guided,32)
    #pragma omp for schedule(static)
    for(tid=0;tid<nthreads;tid++)
      for(i=Start[tid];i<End[tid];i++)
      {

        if(test[i]) continue;  // Salta a la siguiente

        if(i%100000==0){printf("%i\r",i);fflush(stdout);}

        test[i] = true;
        busv_rec(i,Padre,test);
      }

  }  /****** TERMINA SECCION PARALELA ****************************/

  for(i=0;i<cp.npart;i++)
  {
    Padre[i] = Raiz(i,Padre);
    Size[Padre[i]]++;
  }

  for(i=0; i<cp.npart; i++)
  {
    if(Padre[i]==i && Size[i]>1)
    {
      iden.ngrupos++;
      P[i].gr = iden.ngrupos;
    }
  }
   
  for(i=0; i<cp.npart; i++)
  {
    P[i].gr = P[Padre[i]].gr; 
  } 
  
  iden.ngrupos+=1;
  free(Padre);
  free(Size);
  free(test);

  fprintf(stdout,"Nro. grupos identificados = %d\n",iden.ngrupos);
  fprintf(stdout,"Termino identificacion\n"); fflush(stdout);

}

int Raiz(int i, int *Padre)
{

 if(i != Padre[i])
   Padre[i] = Raiz(Padre[i],Padre);

 return Padre[i];

}

void Unir(int u, int v, int *Padre)
{
  int z;

  while(Padre[u] != Padre[v])
  { 
      if(Padre[u] < Padre[v])
      {

          if(u == Padre[u]){
              Padre[u] = Padre[v];   
              break;             
          }

          z = Padre[u];   
          Padre[u] = Padre[v];   
          u = z;

      }else{

          if(v == Padre[v]){
              Padre[v] = Padre[u];   
              break;             
          }

          z = Padre[v];   
          Padre[v] = Padre[u];   
          v = z;
      }
  }

  return;
}

void busv_rec(int ic, int *Padre, bool *test)
{

  long ixc, iyc, izc;
  long ixci, iyci, izci;
  long ixcf, iycf, izcf;
  long ix, iy, iz;
  long ixx, iyy, izz;
  long ibox;
  type_real xx, yy, zz;
  type_real dis;
  int i;
  type_real lbox,fac,lbox2;
  long ngrid;
  item *tmp, *clone;
  unsigned int count;


  ngrid = grid.ngrid;
  lbox  = cp.lbox;
  fac   = (type_real)ngrid/lbox;
  lbox2 = lbox/2.0;

  ixc  = (int)(P[ic].Pos[0]*fac);
  ixci = ixc - 1;
  ixcf = ixc + 1;
  iyc  = (int)(P[ic].Pos[1]*fac);
  iyci = iyc - 1;
  iycf = iyc + 1;
  izc  = (int)(P[ic].Pos[2]*fac);
  izci = izc - 1;
  izcf = izc + 1;

  #ifndef PERIODIC
  if( ixci < 0 ) ixci = 0;
  if( iyci < 0 ) iyci = 0;
  if( izci < 0 ) izci = 0;
  if( ixcf >= ngrid ) ixcf = ngrid - 1;
  if( iycf >= ngrid ) iycf = ngrid - 1;
  if( izcf >= ngrid ) izcf = ngrid - 1;
  #endif

  count = 0;
  for(ixx = ixci; ixx <= ixcf; ixx++){
    ix = ixx;
    #ifdef PERIODIC
    if(ix >= ngrid) ix = ix - ngrid;
    if(ix < 0) ix = ix + ngrid;
    #endif
    for( iyy = iyci ; iyy <= iycf ; iyy++){
      iy = iyy;
      #ifdef PERIODIC
      if(iy >= ngrid) iy = iy - ngrid;
      if(iy < 0) iy = iy + ngrid;
      #endif

      for( izz = izci ; izz <= izcf ; izz++){
        iz = izz;
        #ifdef PERIODIC
        if(iz >= ngrid) iz = iz - ngrid;
        if(iz < 0) iz = iz + ngrid;
        #endif

        ibox = (ix * ngrid + iy) * ngrid + iz ;

        i = grid.llirst[ibox];

        while(i != -1)
        {
          #ifdef IDENSUB
          if(P[i].sub != P[ic].sub){
            i = grid.ll[i];
            continue;
          }
          #endif

          if(test[i] == true)
          {
            i = grid.ll[i];
            continue;
          }      

          xx = P[i].Pos[0] - P[ic].Pos[0];
          yy = P[i].Pos[1] - P[ic].Pos[1];
          zz = P[i].Pos[2] - P[ic].Pos[2];

          #ifdef PERIODIC
          if( xx >  lbox2 ) xx = xx - lbox;
          if( yy >  lbox2 ) yy = yy - lbox;
          if( zz >  lbox2 ) zz = zz - lbox;
          if( xx < -lbox2 ) xx = xx + lbox;
          if( yy < -lbox2 ) yy = yy + lbox;
          if( zz < -lbox2 ) zz = zz + lbox;
          #endif

          dis = xx*xx + yy*yy + zz*zz;

          if(dis < iden.r0)
          {
              //#pragma omp critical
              Unir(ic,i,Padre);

              if(i>ic)
              {
               test[i] = true;
               busv_rec(i,Padre,test);
              }
          }

          i = grid.ll[i];

        } /*fin lazo particulas del grid*/
      } /*fin izz*/
    } /*fin iyy*/
  } /*fin ixx*/

  return 0;
}
