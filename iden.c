#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <omp.h>

#include "variables.h"
#include "cosmoparam.h"
#include "grid.h"
#include "iden.h"
#include "bitmask.h"
#include "limpieza.h"

void identification(void)
{
  int i, nthreads, tid;
  int *test, *Size;
  int ngrid_old = grid.ngrid;
  #ifdef LOCK
  omp_lock_t *lock; lock = (omp_lock_t *) malloc(cp.npart*sizeof(omp_lock_t));
  #endif

  grid.ngrid = (int)(cp.lbox/iden.r0);
  Size = (int *) calloc(cp.npart,sizeof(int));
  test = (int *) calloc(cp.npart/32 + 1,sizeof(int));

  if(grid.ngrid > NGRIDMAX){
    fprintf(stdout,"Using NGRIDMAX = %d\n",NGRIDMAX);
    grid.ngrid = NGRIDMAX;
  }

  grid.nobj = iden.nobj;
  grid.step = iden.step;

  if(ngrid_old != grid.ngrid || iden.step==1)
  {
    grid_free();
    grid_init();
    grid_build();
  }

  iden.r0 *= iden.r0;

  #ifdef NTHREADS
  omp_set_dynamic(0);
  omp_set_num_threads(NTHREADS);
  #endif

  fprintf(stdout,"Comienza identificacion.....\n");

  printf("\nRunning on %d threads\n",NTHREADS);
  nthreads = iden.step == 0 ? (1024<<NTHREADS) : (1024<<NTHREADS)/(2<<iden.step);
  printf("chunks %d for threads\n",nthreads);

  #ifdef LOCK
    #pragma omp parallel default(none) private(tid,i) \
    shared(P,iden,nthreads,cp,test,Size,lock,stdout)   
  #else
    #pragma omp parallel default(none) private(tid,i) \
    shared(P,iden,nthreads,cp,test,Size,stdout)   
  #endif
  {
    tid = omp_get_thread_num(); 
    //nthreads = omp_get_num_threads();
    
    for(i = tid*floor((float)cp.npart/NTHREADS);
    i<(tid==NTHREADS-1 ? cp.npart : (tid+1)*floor((float)cp.npart/NTHREADS));
    i++)
    {
      P[i].gr = i;
      if(iden.step!=0 && P[i].sub == 0)
       SetBit(test,i);     
      #ifdef LOCK
       omp_init_lock(&(lock[i]));
      #endif
    }

    #pragma omp barrier
    
    //#pragma omp for schedule(dynamic)
    //#pragma omp for schedule(static)
    #pragma omp for schedule(guided)
    for(tid=0;tid<nthreads;tid++)
      for(i = tid*floor((float)cp.npart/nthreads);
      i<(tid==nthreads-1 ? cp.npart : (tid+1)*floor((float)cp.npart/nthreads));
      i++)
      {

        if(TestBit(test,i)) continue ;  // Salta a la siguiente
       
        SetBit(test,i);
        
        #ifdef LOCK
        busv_rec(i,test,lock);
        #else
        busv_rec(i,test);
        #endif
      }

  }  /****** TERMINA SECCION PARALELA ****************************/

  fprintf(stdout,"Sale del paralelo\n"); fflush(stdout);

  #ifdef LOCK
  for(i=0;i<cp.npart;i++) omp_destroy_lock(&(lock[i]));
  free(lock);
  #endif

  linkedlist(test);

  free(Size);
  free(test);

  fprintf(stdout,"Termino identificacion\n"); fflush(stdout);

}

int Raiz(int i)
{

 if(i != P[i].gr)
   P[i].gr = Raiz(P[i].gr);

 return P[i].gr;

}

#ifdef LOCK
void Unir(int u, int v, omp_lock_t *lock)
#else
void Unir(int u, int v)
#endif
{
  int z;

  while(P[u].gr != P[v].gr)
  { 
      if(P[u].gr < P[v].gr)
      {
#ifdef LOCK
          if(u == P[u].gr)
          {
            omp_set_lock(&(lock[u]));
            z = 0;
            if(u == P[u].gr)
            {
                P[u].gr = P[v].gr;  
                z = 1;
            }
            omp_unset_lock(&(lock[u]));
            if(z==1) break;             
          }
#else
          if(u == P[u].gr)
          {
            P[u].gr = P[v].gr;  
            break;             
          }
#endif
          
          z = P[u].gr;   
          P[u].gr = P[v].gr;
          u = z;

      }else{
#ifdef LOCK
          if(v == P[v].gr)
          {
            omp_set_lock(&(lock[v]));
            z = 0;
            if(v == P[v].gr)
            {
                P[v].gr = P[u].gr;   
                z = 1;
            }
            omp_unset_lock(&(lock[v]));
            if(z == 1) break;            
          }
#else
          if(v == P[v].gr)
          {
              P[v].gr = P[u].gr;   
              break;             
          }
#endif

          z = P[v].gr;   
          P[v].gr = P[u].gr;   
          v = z;

      }
  }

  return;
}

#ifdef LOCK
void busv_rec(int ic, int *test, omp_lock_t *lock)
#else
void busv_rec(int ic, int *test)
#endif
{

  long ixc, iyc, izc;
  long ixci, iyci, izci;
  long ixcf, iycf, izcf;
  long ix, iy, iz;
  long ixx, iyy, izz;
  long ibox;
  my_real xx, yy, zz;
  my_real dis;
  int i;
  my_real lbox,fac;
  long ngrid;


  ngrid = grid.ngrid;
  lbox  = cp.lbox;
  fac   = (my_real)ngrid/lbox;
#ifdef PERIODIC
  my_real lbox2 = lbox/2.0;
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
 
          if(TestBit(test,i))
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
#ifdef LOCK
             Unir(ic,i,lock);
#else
             Unir(ic,i);
#endif
             if(i>ic)
             {
              SetBit(test,i);
#ifdef LOCK
              busv_rec(i,test,lock);
#else
              busv_rec(i,test);
#endif
             }
          }

          i = grid.ll[i];

        } /*fin lazo particulas del grid*/
      } /*fin izz*/
    } /*fin iyy*/
  } /*fin ixx*/

}

void linkedlist(int *test)
{
  int i,g;

  Temp.ll = (int *) calloc(cp.npart,sizeof(int));

  iden.ngrupos = 0;
  for(i=0;i<cp.npart;i++)
  {
    P[i].gr = Raiz(i);
    if(TestBit(test,P[i].gr))
    {
      Temp.ll[P[i].gr]++;
      if(Temp.ll[P[i].gr]>=NPARTMIN)
      { 
        iden.ngrupos++;
        Temp.ll[P[i].gr] = iden.ngrupos;
        ClearBit(test,P[i].gr);
      }
    }
  }

  iden.ngrupos++;  // SUMA UNO;

  Temp.head   = (int *) malloc(iden.ngrupos*sizeof(int));
  Temp.npgrup = (unsigned int *) malloc(iden.ngrupos*sizeof(unsigned int));

  for(i=0;i<cp.npart;i++)
  {
    //if(TestBit(test,P[i].gr))
    //{ 
    //  Temp.ll[P[i].gr] = 0;
    //  ClearBit(test,P[i].gr);
    //}

    //P[i].gr = Temp.ll[P[i].gr];
    if(TestBit(test,P[i].gr))
      P[i].gr = 0;
    else
      P[i].gr = Temp.ll[P[i].gr];

    if(i<iden.ngrupos)
    {
      Temp.head[i] = -1;
      Temp.npgrup[i] = 0;
    }
  }

  for(i=0;i<cp.npart;i++)
  {
    g = P[i].gr;

    #ifdef DEBUG
    assert((g >= 0) && (g < iden.ngrupos));
    #endif
    Temp.ll[i] = Temp.head[g];
    Temp.head[g] = i;
    Temp.npgrup[g]++;
  }
}
