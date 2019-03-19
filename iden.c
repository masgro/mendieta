#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <omp.h>
#include <string.h>

#include "variables.h"
#include "cosmoparam.h"
#include "grid.h"
#include "iden.h"
#include "bitmask.h"

#define MY_DIM 0

struct global_head_st{
  item *head;
  struct global_head_st *global_next;
};

void identification(void){
  unsigned int nvec;
  int  nthreads, tid;
  my_int i,j,ngrupo;
  item *tmp,*tmp1;
  item *curr,*head;
  int *test;
  my_real zcm;
  my_int nn;
  my_int *grupos_per_thread,ntotal;
  void *puntero;
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
  my_real *zmin,*zmax;
#pragma omp parallel default(none) reduction(+:ntotal) \
  private(head,curr,nvec,tmp,tmp1,test,tid,i,zcm,nn,   \
          curr_global_head,puntero,j,nthreads) \
  shared(P,iden,global_head,cp,grupos_per_thread,ngrupo,zmin,zmax,stdout)                     
{
  tid = omp_get_thread_num();

  nthreads = omp_get_num_threads();
  //zmin = (my_real)tid*cp.lbox/(my_real)nthreads;
  //zmax = (my_real)(tid+1)*cp.lbox/(my_real)nthreads;

  if(tid == 0){
    my_int nhist = 160;
    my_int *zhist;
    my_real *zhist_lim;
    zhist = (my_int *) calloc(nhist,sizeof(my_int));
    zhist_lim = (my_real *) malloc(nhist*sizeof(my_real));
    zmin = (my_real *) malloc(nthreads*sizeof(my_real));
    zmax = (my_real *) malloc(nthreads*sizeof(my_real));

    printf("\nRunning on %d threads\n",nthreads);
    grupos_per_thread = (my_int *) malloc(nthreads*sizeof(my_int));
    global_head = (struct global_head_st **) malloc(nthreads*sizeof(struct global_head_st *));

//    my_real med,sig;
//    for(j = 0; j < 3; j++){
//      memset(zhist,0,nhist*sizeof(my_real));
//      for(i = 0; i < iden.nobj; i++)
//        zhist[(my_int)(P[i].Pos[j]/cp.lbox*(my_real)nhist)]++;
//
//      med = 0.;
//      for(i = 0; i < nhist; i++){
//        med += (my_real)zhist[i];
//      }
//      med /= (my_real)nhist;
//      sig = 0.;
//      for(i = 0; i < nhist; i++){
//        sig += (my_real)((zhist[i]-med)*(zhist[i]-med));
//      }
//      sig /= (my_real)(nhist - 1);
//      printf("%d %f %f\n",j,sqrt(sig),med);
//    }

    for(i = 0; i < iden.nobj; i++)
      zhist[(my_int)(P[i].Pos[MY_DIM]/cp.lbox*(my_real)nhist)]++;

    for(i = 0; i <= nhist; i++) zhist_lim[i] = (my_real)i*cp.lbox/(my_real)nhist;

    my_int sum;
    sum = 0;
    for(j = 0, i = 0; j < nthreads; j++){
      zmin[j] = zhist_lim[i];
      while(1){
        sum += zhist[i];
        i++;
        if((sum > (my_int)((my_real)(j+1)*(my_real)iden.nobj/(my_real)nthreads)) || (i >= nhist)){
          break;
        }
      }
      zmax[j] = zhist_lim[i];
    }
  }

  #pragma omp barrier
  test = (int *) calloc(iden.nobj/32 + 1,sizeof(int));

  global_head[tid] = NULL;
  for(i = 0; i < iden.nobj; i++){
    if(TestBit(test,i))continue;
    //if(P[i].Pos[2] < zmin || P[i].Pos[2] >= zmax) continue;
    if(P[i].Pos[MY_DIM] < zmin[tid] || P[i].Pos[MY_DIM] >= zmax[tid]) continue;

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
      zcm += P[curr->indx].Pos[MY_DIM]; 
      nn++;

      curr = curr->next;
    }
    zcm /= (my_real)nn;

    //if(zcm < zmin || zcm >= zmax || nn < NPARTMIN){
    if(zcm < zmin[tid] || zcm >= zmax[tid] || nn < NPARTMIN){
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
  long ixci, iyci, izci;
  long ixcf, iycf, izcf;
  long ix, iy, iz;
  long ixx, iyy, izz;
  long ibox,ngrid;
  my_int i,ic,count;
  my_real xx, yy, zz;
  my_real dis;
  my_real lbox,fac;
  item *tmp, *clone;

  clone = *curr;
  ic = clone->indx;

  ngrid = grid.ngrid;
  lbox  = cp.lbox;
  fac   = (my_real)ngrid/lbox;
  #ifdef PERIODIC
  my_real lbox2 = lbox/2.0;
  #endif

  ixcf = (long)(P[ic].Pos[0]*fac);
  ixci = ixcf - 1;
  ixcf = ixcf + 1;
  iycf = (long)(P[ic].Pos[1]*fac);
  iyci = iycf - 1;
  iycf = iycf + 1;
  izcf = (long)(P[ic].Pos[2]*fac);
  izci = izcf - 1;
  izcf = izcf + 1;

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
        do{
          #ifdef IDENSUB
          if(P[i].sub == P[ic].sub){
          #endif
            if(!TestBit(test,i)){
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
            }
          #ifdef IDENSUB 
          }
          #endif
          i = grid.ll[i];
        }while(i != grid.llirst[ibox]); /*fin lazo particulas del grid*/
      } /*fin izz*/
    } /*fin iyy*/
  } /*fin ixx*/

  *curr = clone;
  *nvec += count;
}
