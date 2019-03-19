#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <omp.h>

#include "variables.h"
#include "cosmoparam.h"
#include "grid.h"
#include "identest.h"
#include "bitmask.h"

struct global_head_st{
  item *head;
  struct global_head_st *global_next;
};

my_int *testo;

void identification(void){
  my_int nvec;
  my_int i,j,ngrupo;
  //my_int *head;
  int    *test;
  my_int nn;
  //void   *puntero;
  my_int ntotal;
  //struct global_head_st *global_head;
  //struct global_head_st *curr_global_head;

  int ngrid_old = grid.ngrid;
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

  for(i = 0; i < iden.nobj; i++){
    P[i].gr = 0;
    P[i].llfof = i;
  }

  #ifdef NTHREADS
  omp_set_dynamic(0);
  omp_set_num_threads(NTHREADS);
  #endif

  fprintf(stdout,"Comienza identificacion.....\n");

  test = (int *) calloc(iden.nobj/32 + 1,sizeof(int));
  //head = (my_int *) malloc(iden.nobj*sizeof(my_int));
  testo = (my_int *) calloc(iden.nobj,sizeof(my_int));

  ngrupo = 0;
  for(i = 0; i < iden.nobj; i++){
    //if(TestBit(test,i))continue;
    if(testo[i] == 1)continue;
    
    //SetBit(test,i);
    testo[i] == 1;
    ngrupo++;

    //groups[ngrupo].llirst = i;
    //groups[ngrupo].np = 1;

    my_int curr = i, tmp, tmp1;
    do{
      tmp = P[curr].llfof;
      tmp1 = curr;

      busv(&curr,test);

      P[curr].llfof = tmp;
      curr = P[tmp1].llfof;
    }while(curr != i);
  }
  free(test);
  free(testo);

  //ntotal = 0;
  //ngrupo = 0;
  //curr_global_head = global_head;
  //while(curr_global_head){
  //  ntotal++;
  //  head = curr_global_head->head;
  //  if(P[head->indx].gr == 0){
  //    ngrupo++;
  //    j = ngrupo;
  //    #ifdef DEBUG
  //    i = 0;
  //    #endif
  //    curr = head;
  //    while(curr){
  //      #ifdef DEBUG
  //      i++;
  //      assert(P[curr->indx].gr == 0);
  //      #endif
  //      P[curr->indx].gr = j;
  //      curr = curr->next;
  //    }
  //    #ifdef DEBUG
  //    assert(i >= NPARTMIN);
  //    #endif
  //  }
  //  curr_global_head = curr_global_head->global_next;
  //}

  //curr_global_head = global_head;
  //while(curr_global_head){
  //  curr = curr_global_head->head;
  //  while(curr){
  //    puntero = curr->next;
  //    free(curr);
  //    curr = puntero;
  //  }
  //  puntero = curr_global_head->global_next;
  //  free(curr_global_head);
  //  curr_global_head = puntero;
  //}
  /****** TERMINA SECCION PARALELA ****************************/

  /*Le sumamos 1 para contar el grupo 0*/
  iden.ngrupos = ngrupo + 1;
  fprintf(stdout,"Nro. grupos identificados = %d\n",iden.ngrupos);
  fprintf(stdout,"Termino identificacion\n"); fflush(stdout);
}

void busv(my_int *curr, int *test){
  long ixc, iyc, izc;
  long ixci, iyci, izci;
  long ixcf, iycf, izcf;
  long ix, iy, iz;
  long ixx, iyy, izz;
  long ibox;
  double xx, yy, zz;
  double dis;
  my_int ic;
  my_real lbox,fac,lbox2;
  long ngrid;
  my_int tmp;
  my_int count;

  ic = *curr;

  ngrid = grid.ngrid;
  lbox  = cp.lbox;
  fac   = (my_real)ngrid/lbox;
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

  unsigned int box_list[27];
  unsigned int j = 0;
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
        box_list[j] = ibox;
        j++;
      }
    }
  }

  count = 0;
  my_int head[27],tail[27];;
  my_int cuentas[27];
  #pragma omp parallel default(none)\
  private(ibox,xx,yy,zz,dis,tmp,j)\
  shared(P,ic,lbox,lbox2,test,grid,head,box_list,iden,curr,stdout,testo,count,cuentas,tail)
  {
    my_real Pos[3];
    my_int  i;
    int tid = omp_get_thread_num();

    Pos[0] = P[ic].Pos[0];
    Pos[1] = P[ic].Pos[1];
    Pos[2] = P[ic].Pos[2];

    #pragma omp for reduction(+:count)
    for(j = 0; j < 27; j++){
      //printf("%d\n",j);fflush(stdout);
      cuentas[j] = 0;
      head[j] = ic;
      tail[j] = ic;

      ibox = box_list[j];
      i = grid.llirst[ibox];
      do{
        #ifdef IDENSUB
        if(P[i].sub == P[ic].sub){
        #endif
          //if(!TestBit(test,i)){
          if(testo[i] == 0){


            xx = P[i].Pos[0] - Pos[0];
            yy = P[i].Pos[1] - Pos[1];
            zz = P[i].Pos[2] - Pos[2];

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
              //SetBit(test,i);
              testo[i] = 1;

              if(head[j] == ic) tail[j] = i;

              P[i].llfof = head[j];
              head[j] = i;

              count++;
              cuentas[j]++;
            }
          }
        #ifdef IDENSUB
        }
        #endif
        i = grid.ll[i];
      }while(i != grid.llirst[ibox]);
    }


  }
  for(j = 1; j < 27; j++) cuentas[0] += cuentas[j];
  assert(cuentas[0] == count);

  //printf("%u entra al master %u\n",ic,count);fflush(stdout);
  //my_int lala = 0;
  //i = ic;
  //for(j = 0; j < 27; j++){
  //  tmp = head[j];
  //  while(tmp != ic){

  //    lala++;

  //    P[i].llfof = tmp;
  //    i = tmp;

  //    tmp = P[tmp].llfof;
  //  }
  //}
  my_int i = ic;
  for(j = 0; j < 27; j++){
    if(head[j] == ic)continue;
    P[i].llfof = head[j];
    i = tail[j];
  }
  *curr = i;
  //printf("%u sale al master %u\n",ic,lala);fflush(stdout);
  //assert(lala == count);
  //printf("%d sale al master\n",i);fflush(stdout);
}
