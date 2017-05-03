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
#include "bitmask.h"

struct global_head_st{
  item *head;
  struct global_head_st *global_next;
};

void identification(void){
  unsigned int nvec;
  int  nthreads, tid;
  int  i,j,ngrupo;
  item *tmp,*tmp1;
  item *curr,*head;
  int *test;
  type_real zmin,zmax,zcm;
  int nn;
  int *grupos_per_thread;
  void *puntero;
  int ntotal;
  struct global_head_st **global_head;
  struct global_head_st *curr_global_head;

  unsigned long ngrid_old = grid.ngrid;
  grid.ngrid = (int)(cp.lbox/iden.r0);

  if(grid.ngrid > NGRIDMAX){
    fprintf(stdout,"Using NGRIDMAX = %d\n",NGRIDMAX);
    grid.ngrid = NGRIDMAX;
  }

  grid.nobj = iden.nobj;

  if(ngrid_old != grid.ngrid){
    grid_free();
    grid_init();
    grid_build();
  }

  iden.r0 = iden.r0*iden.r0;

  for(i = 0; i < iden.nobj; i++)
    P[i].gr = 0;

  #ifdef NTHREADS
  omp_set_dynamic(0);
  omp_set_num_threads(NTHREADS);
  #endif

  fprintf(stdout,"Comienza identificacion.....\n");
  iden.ngrupos = 0;

  ntotal = 0; ngrupo = 0;
#pragma omp parallel default(none) reduction(+:ntotal) \
  private(head,curr,nvec,tmp,tmp1,test,tid,i,zcm,nn,   \
          curr_global_head,puntero,j,zmin,zmax,nthreads) \
  shared(P,iden,global_head,cp,grupos_per_thread,ngrupo)                     
{
  tid = omp_get_thread_num();

  nthreads = omp_get_num_threads();
  zmin = (type_real)tid*cp.lbox/(type_real)nthreads;
  zmax = (type_real)(tid+1)*cp.lbox/(type_real)nthreads;

  if(tid == 0){
    printf("\nRunning on %d threads\n",nthreads);
    grupos_per_thread = (int *) malloc(nthreads*sizeof(int));
    global_head = (struct global_head_st **) malloc(nthreads*sizeof(struct global_head_st *));
  }
  #pragma omp barrier

  test = (int *) calloc(iden.nobj/32 + 1,sizeof(int));

  global_head[tid] = NULL;
  for(i =0; i < iden.nobj; i++){
    if(TestBit(test,i))continue;
    if(P[i].Pos[2] < zmin || P[i].Pos[2] >= zmax) continue;

    SetBit(test,i);

    head = (item *) malloc(sizeof(item));
    head->indx = i;
    head->next = NULL;

    curr = head;
    while(curr){
      tmp = curr->next;
      tmp1 = curr;

      nvec = 0;
      busv(&curr, &nvec, test);

      curr->next = tmp;
      curr = tmp1->next;
    }

    zcm = 0.0; nn = 0;
    curr = head;
    while(curr){
      zcm += P[curr->indx].Pos[2]; 
      nn++;

      curr = curr->next;
    }
    zcm /= (type_real)nn;

    if(zcm < zmin || zcm >= zmax || nn < NPARTMIN){
      curr = head;
      while(curr){
        puntero = curr->next;
        free(curr);
        curr = puntero;
      }
    }else{
      ntotal++;
      curr_global_head = (struct global_head_st *) malloc(sizeof(struct global_head_st));
      assert(curr_global_head != NULL);
      curr_global_head->head = head;
      curr_global_head->global_next = global_head[tid];
      global_head[tid] = curr_global_head;
    }
  }
  grupos_per_thread[tid] = ntotal;
  printf("tid %d ntotal %d\n",tid,ntotal);
  free(test);
  #pragma omp barrier

  /***** RENOMBRA *****/
  #pragma omp for schedule(dynamic)
  for(tid = 0; tid < nthreads; tid++){
    ntotal = 0;
    curr_global_head = global_head[tid];
    while(curr_global_head){
      ntotal++;
      head = curr_global_head->head;
      if(P[head->indx].gr == 0){
        #pragma omp critical
        {ngrupo++;j = ngrupo;}
        #ifdef DEBUG
        i = 0;
        #endif
        curr = head;
        while(curr){
          #ifdef DEBUG
          i++;
          assert(P[curr->indx].gr == 0);
          #endif
          P[curr->indx].gr = j;
          curr = curr->next;
        }
        #ifdef DEBUG
        assert(i >= NPARTMIN);
        #endif
      }
      curr_global_head = curr_global_head->global_next;
    }
    //assert(ntotal == grupos_per_thread[tid]);

    curr_global_head = global_head[tid];
    while(curr_global_head){
      curr = curr_global_head->head;
      while(curr){
        puntero = curr->next;
        free(curr);
        curr = puntero;
      }
      puntero = curr_global_head->global_next;
      free(curr_global_head);
      curr_global_head = puntero;
    }
  }
}
  /****** TERMINA SECCION PARALELA ****************************/

  /*Le sumamos 1 para contar el grupo 0*/
  iden.ngrupos = ngrupo + 1;
  free(global_head);
  free(grupos_per_thread);
  fprintf(stdout,"Nro. grupos identificados = %d\n",iden.ngrupos);
  fprintf(stdout,"Termino identificacion\n"); fflush(stdout);
}

void busv(item **curr, unsigned int *nvec, int *test){
  long ixc, iyc, izc;
  long ixci, iyci, izci;
  long ixcf, iycf, izcf;
  long ix, iy, iz;
  long ixx, iyy, izz;
  long ibox;
  type_real xx, yy, zz;
  type_real dis;
  int i;
  type_real lbox,fac;
  long ngrid;
  int ic;
  item *tmp, *clone;
  unsigned int count;

  clone = *curr;
  ic = clone->indx;

  ngrid = grid.ngrid;
  lbox  = cp.lbox;
  fac   = (type_real)ngrid/lbox;
  #ifdef PERIODIC
  type_real lbox2 = lbox/2.0;
  #endif

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
        while(i != -1){
          #ifdef IDENSUB
          if(P[i].sub != P[ic].sub){
            i = grid.ll[i];
            continue;
          }
          #endif
          if(TestBit(test,i)){
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

          if(dis < iden.r0){
            SetBit(test,i);

            tmp = (item *) malloc(sizeof(item));
            #ifdef DEBUG
            assert(tmp != NULL);
            #endif

            tmp->indx = i;
            tmp->next = NULL;

            clone->next = tmp;
            clone = tmp;

            count++;
          }

          i = grid.ll[i];

        } /*fin lazo particulas del grid*/
      } /*fin izz*/
    } /*fin iyy*/
  } /*fin ixx*/

  *curr = clone;
  *nvec += count;
}
