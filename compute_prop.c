#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <omp.h>

#include "variables.h"
#include "colores.h"
#include "cosmoparam.h"
#include "leesnap.h"
#include "octree.h"
#include "deltas.h"
#include "timer.h"
#include "compute_prop.h"
#include "deltas.h"

#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector_float.h>
#include <gsl/gsl_sort_vector_float.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_math.h>

/* ============= Unidades =============== 
   Energias = 10^10 Msol (km/seg)^2 
   Mom. Angular = 10^10 Msol kpc km/seg
   M(200,vir) = 10^10 Msol/h 
   R(200,vir) = kpc/h 
   V(200,vir) = km/seg
   ====================================== */

void compute_properties(struct grupos *g){
  int    i;
  char   filename[200];
  char   fout[200];
  FILE   *pfout;
  FILE   *asco;
  struct propiedades_st Prop;
  unsigned int n_grupos;

  Prop.evec = gsl_matrix_alloc(3,3);
  Prop.evec_vel = gsl_matrix_alloc(3,3);

  if(PROP_FOF == true){
    n_grupos = n_grupos_fof;
    sprintf(fout,"properties.fof");
  }

  #ifdef IDENSUB
  if(PROP_SUB == true){
    n_grupos = n_grupos_sub;
    sprintf(fout,"properties.sub");
  }
  #endif

  /* Si no encuentra grupos, crea archivo vacio */
  if(n_grupos == 0){
    printf("Creando archivo vacio y saliendo.\n");
    pfout = fopen(fout,"w");
    printf("Abre archivo de salida: %s\n",fout);
    fwrite(&n_grupos,sizeof(int),1,pfout);
    fclose(pfout);
    exit(EXIT_SUCCESS);
  }

  pfout = fopen(fout,"w");
  fprintf(stdout,"Abre archivo de salida: %s\n",fout);

  #ifdef IDENSUB
  sub_hijos = (struct subst *) calloc((n_grupos_fof + 1),sizeof(struct subst));
  #endif

  #ifdef ASCII_PROPERTIES_FILE
  asco = fopen("propiedades.ascii","w");
  #endif

  fwrite(&n_grupos,sizeof(int),1,pfout);
  for(i = 1; i <= n_grupos; i++){
    if(g[i].np == 0)continue;
    propiedades(P,g,i,&Prop);

    write_properties(pfout,&Prop);

    #ifdef ASCII_PROPERTIES_FILE
    fprintf(asco,"%d %f %f %f %f %f %f ",
         Prop.npart,Prop.pcm[0],Prop.pcm[1],Prop.pcm[2],Prop.vcm[0],Prop.vcm[1],Prop.vcm[2]);
    fprintf(asco,"%f %f %f %f %f %f %f ",
         Prop.sig[0],Prop.sig[1],Prop.sig[2],Prop.L[0],Prop.L[1],Prop.L[2],Prop.lambda);
    fprintf(asco,"%f %f %f %f %f %f %f ",
         Prop.m200,Prop.r200,Prop.v200,Prop.mvir,Prop.rvir,Prop.vvir,Prop.vmax);
    fprintf(asco,"%f %f %f %f %f %f %f %f\n",
         Prop.Ep,Prop.Ec,Prop.aa,Prop.bb,Prop.cc,Prop.aa_vel,Prop.bb_vel,Prop.cc_vel);
    fflush(asco);
    #endif

    #ifdef IDENSUB
    sub_hijos[i].pcm[0] = Prop.mostbound[0];
    sub_hijos[i].pcm[1] = Prop.mostbound[1];
    sub_hijos[i].pcm[2] = Prop.mostbound[2];
    sub_hijos[i].r200   = Prop.rvir;
    #endif
    
  }
  fclose(pfout);

  #ifdef ASCII_PROPERTIES_FILE
  fclose(asco);
  #endif

  #ifdef IDENSUB
  //sprintf(filename,"%s",fhijos);
  sprintf(filename,"hijos.ascii");
  asco = fopen(filename,"w");
  int *hijo,j,k,l;
  double dx,dy,dz,d;

  hijo = (int *) malloc((n_grupos_fof+1)*sizeof(int));

  for(i = 0; i < n_grupos_fof; i++){
    for(l = 0; l < n_grupos_fof; l++) hijo[l] = -1;
    k=0;
    for(j = i + 1; j < n_grupos_fof; j++){
      dx = sub_hijos[i].pcm[0] - sub_hijos[j].pcm[0];
      dy = sub_hijos[i].pcm[1] - sub_hijos[j].pcm[1];
      dz = sub_hijos[i].pcm[2] - sub_hijos[j].pcm[2];
      d = sqrt(dx*dx + dy*dy + dz*dz);
      if(d < sub_hijos[i].r200){
        hijo[k] = j;
        k++;
      }
    }
    fprintf(asco,"%d %d\n",i,k);
    for(l = 0; l < k; l++)
      fprintf(asco,"%d\n",hijo[l]);
  }
  fclose(asco);
  #endif
}

void propiedades(struct particle_data *P, struct grupos *g, int gid, struct propiedades_st *PropIN){
  float  dis,_vmax;
  float  rho;
  float  _t1, _t2;
  float  Theta = 0.45;
  int    dim,l,i,j;
  int    indx,iepmin;
  int    i200, ivir;
  float  drt;
  type_real  dx[3],dv[3];
  float a,lbox2,Epmin;
  gsl_vector_float *_dis;
  gsl_permutation  *_ind;
  struct node *current;
  struct particle_data *Q;
  struct propiedades_st Prop = *PropIN;

  lbox2 = (float)cp.lbox*0.5f;
  a = (1.0 / (1.0 + cp.redshift));

  Prop.npart = g[gid].np;

  for(dim = 0; dim < 3; dim++){
    Prop.pcm[dim] = 0.;
    Prop.vcm[dim] = 0.;
    Prop.L[dim]   = 0.;
    Prop.sig[dim] = 0.;
  }

  Q  = (struct particle_data *) malloc(Prop.npart*sizeof(struct particle_data));

  i = 0;
  l = g[gid].llirst;
  while(l != -1){
    Q[i].Pos[0] = P[l].Pos[0];
    Q[i].Pos[1] = P[l].Pos[1];
    Q[i].Pos[2] = P[l].Pos[2];
    Q[i].Vel[0] = P[l].Vel[0];
    Q[i].Vel[1] = P[l].Vel[1];
    Q[i].Vel[2] = P[l].Vel[2];

    for(dim = 0; dim < 3; dim++){
      drt = Q[i].Pos[dim] - Q[0].Pos[dim];
      Prop.pcm[dim] += Q[i].Pos[dim];
      Prop.vcm[dim] += Q[i].Vel[dim];
    }

    i++;

    if(PROP_FOF == true)
      l = P[l].llfof;

    #ifdef IDENSUB
    if(PROP_SUB == true)
      l = P[l].llsub;
    #endif
   }

  #ifdef DEBUG
  assert(i == Prop.npart);
  #endif

  for(dim = 0; dim < 3; dim++){
    Prop.pcm[dim] /= (type_real)Prop.npart;
    Prop.vcm[dim] /= (type_real)Prop.npart;
  }

  /* Si el halo tiene menos de 1000 particulas calcula la energia
  de forma directa (N^2). Si tiene mas de 1000 particulas usa 
  un octree. */
  if(Prop.npart < 1000){
    #pragma omp parallel for default(none) schedule(dynamic) \
    shared(Prop,Q,cp,a) \
    private(i,j,dim,dis,dx,dv)
    for(i = 0; i < Prop.npart; i++){
      /* Calcula energia cinetica (velocidades en [km / seg]) */
      for(dim = 0; dim < 3; dim++){
        dv[dim]  = cp.Hubble_a*a*(Q[i].Pos[dim]-Prop.pcm[dim])/1000.0f; /* Vel de Hubble */
         dv[dim] += (float)sqrt(a)*(Q[i].Vel[dim]-Prop.vcm[dim]);
      }

      Q[i].Ec  = dv[0]*dv[0];
      Q[i].Ec += dv[1]*dv[1];
      Q[i].Ec += dv[2]*dv[2];
      Q[i].Ec *= 0.5;

      #ifndef COMPUTE_EP
      /* Calcula energia potencial */
      Q[i].Ep  = 0.;
      for(j = 0; j < Prop.npart; j++){
        if(j == i)continue;

        for(dim = 0; dim < 3; dim++)
          dx[dim] = Q[i].Pos[dim] - Q[j].Pos[dim];

        dis  = dx[0]*dx[0];
        dis += dx[1]*dx[1];
        dis += dx[2]*dx[2];
        dis  = (float)sqrt(dis);
        dis += (float)cp.soft;

        Q[i].Ep += 1./dis;

      }
      Q[i].Ep += (1./cp.soft);            /* Autoenergia */
      Q[i].Ep *= (GCONS*cp.Mpart*Msol*1.E10/Kpc/a);
      Q[i].Ep *= (-1.);                   /* Cambio de signo para que Ep sea negativa */
      #endif
    }
  }else{
    #ifndef COMPUTE_EP
    force_treeallocate(2 * Prop.npart + 200);
    force_treebuild(Prop.npart,Q,Theta);
    #endif

    #pragma omp parallel for default(none) schedule(dynamic) \
    shared(Prop,Q,cp,a) \
    private(i,dim,dv)
    for(i = 0; i < Prop.npart; i++){
      #ifndef COMPUTE_EP
      /* Calcula la energia potencial */
      force_treeevaluate_potential(&Q[i].Pos[0], &Q[i].Ep);
      #ifdef DEBUG
      if(isnan(Q[i].Ep))exit(EXIT_FAILURE);
      #endif
      Q[i].Ep -= cp.Mpart/cp.soft;    /* Autoenergia */
      Q[i].Ep *= GCONS/a*Msol*1.E10/Kpc;
      #endif

      /* Calcula la energia cinetica (velocidades en [km / seg]) */
      for(dim = 0; dim < 3; dim++){
        dv[dim]  = cp.Hubble_a*a*(Q[i].Pos[dim]-Prop.pcm[dim])/1000.0f;  /* Vel de Hubble */
        dv[dim] += (float)sqrt(a)*(Q[i].Vel[dim]-Prop.vcm[dim]);
      }

      Q[i].Ec  = dv[0]*dv[0];
      Q[i].Ec += dv[1]*dv[1];
      Q[i].Ec += dv[2]*dv[2];
      Q[i].Ec *= 0.5;
    }

    #ifndef COMPUTE_EP
    force_treefree();
    #endif
  }

  iepmin = -1; Epmin = 0.;
  for(i = 0; i < Prop.npart; i++){
    if(Q[i].Ep < Epmin){
      Epmin  = (float)Q[i].Ep;
      iepmin = i;
    }
  }

  #ifdef DEBUG
  assert(iepmin != -1);
  #endif

  for(dim = 0; dim < 3; dim++)
    Prop.mostbound[dim] = Q[iepmin].Pos[dim];

  /*Ordena por distancia al centro de potencial*/
   _dis = gsl_vector_float_alloc(Prop.npart);
  _ind = gsl_permutation_alloc(Prop.npart);
  Prop.Ep = 0.;
  Prop.Ec = 0.;
  for(i = 0; i < Prop.npart; i++){
    /* Energia cinetica y potencial del halo */
    Prop.Ep += (float)Q[i].Ep;
    Prop.Ec += (float)Q[i].Ec;

    for(dim = 0; dim < 3; dim++){
      dx[dim] = Q[iepmin].Pos[dim] - Q[i].Pos[dim];
      dv[dim] = Q[i].Vel[dim] - Prop.vcm[dim];
      Prop.sig[dim] += Q[i].Vel[dim]*Q[i].Vel[dim];
    }

    Prop.L[0] += dx[1]*dv[2] - dx[2]*dv[1];
    Prop.L[1] += dx[2]*dv[0] - dx[0]*dv[2];
    Prop.L[2] += dx[0]*dv[1] - dx[1]*dv[0];  

    dis  = dx[0]*dx[0];
    dis += dx[1]*dx[1];
    dis += dx[2]*dx[2];
    dis  = (float)sqrt(fabs(dis));
    gsl_vector_float_set(_dis,i,dis);
  }

#ifdef DEBUG
  assert((0 <= iepmin) && (iepmin < Prop.npart));
#endif

  Prop.Ep *= 0.5f*(float)cp.Mpart;
  Prop.Ec *= (float)cp.Mpart;

  /* Sigma y pasaje a coord fisicas del L, pibe! */
  for(dim = 0; dim < 3; dim++){
    Prop.sig[dim]  = (Prop.sig[dim] - Prop.vcm[dim]*Prop.vcm[dim]*(float)Prop.npart);
    Prop.sig[dim] /= (float)(Prop.npart - 1);
    Prop.sig[dim]  = (float)sqrt(Prop.sig[dim]) ;
    Prop.L[dim]   *= a*(float)(sqrt(a)*cp.Mpart);
  }

  gsl_sort_vector_float_index(_ind,_dis);

  /* Parametro adimensional de Spin */
  Prop.lambda  = (float)sqrt(Prop.L[0]*Prop.L[0]+Prop.L[1]*Prop.L[1]+Prop.L[2]*Prop.L[2]);
  Prop.lambda *= (float)sqrt(fabs(Prop.Ec + Prop.Ep));
  Prop.lambda /= (float)(GCONS/Kpc*Msol*1.E10);
  Prop.lambda /= (float)pow((double)Prop.npart*cp.Mpart,2.5);

  i200 = -1;
  ivir = -1;
  Prop.vmax = 0.;
  _t1  = (float)pow(cp.lbox/(float)cp.npart,3);
  _t1 *= 0.75;
  _t1 /= M_PI;
  i = 4;
  do{
    #ifdef DEBUG
    assert(i<Prop.npart);
    #endif
    indx = gsl_permutation_get(_ind,i);
    #ifdef DEBUG
    assert((0<=indx) && (indx<Prop.npart));
    #endif
    dis = gsl_vector_float_get(_dis,indx);
    _t2  = (float)(i+1);
    rho  = _t1;
    rho *= _t2;
    rho /= (float)(dis*dis*dis);

    if (rho < 200./cp.omegam && i200 == -1) 
      i200 = i;

    if (rho < deltavir(cp.redshift) && ivir == -1) 
      ivir = i;

    _vmax = (float)(GCONS/Kpc*Msol*(double)(i+1)*cp.Mpart)/dis;
    if(_vmax > Prop.vmax) Prop.vmax = _vmax;

    i++;
  }while(i < Prop.npart);

  /* R200 y M200 */
  if(i200 == -1) i200 = Prop.npart-1;
  indx = gsl_permutation_get(_ind,i200);
  dis  = gsl_vector_float_get(_dis,indx);
  Prop.r200 = dis;
  Prop.m200 = (float)((double)(i200+1)*cp.Mpart);
  Prop.v200 = (float)(GCONS/Kpc*Msol*Prop.m200/Prop.r200);
  Prop.v200 = (float)(sqrt(Prop.v200)*1.E5);

  /* Rvir y Mvir */
  if(ivir == -1) ivir = Prop.npart-1;
  indx = gsl_permutation_get(_ind,ivir);
  dis  = gsl_vector_float_get(_dis,indx);
  Prop.rvir = dis;
  Prop.mvir = (float)((double)(ivir+1)*cp.Mpart);
  Prop.vvir = (float) (GCONS/Kpc*Msol*Prop.mvir/Prop.rvir);
  Prop.vvir = (float) (sqrt(Prop.vvir)*1.E5);

  /* Velocidad circular maxima */
  Prop.vmax = (float) (sqrt(Prop.vmax)*1.E5);

  gsl_vector_float_free(_dis);
  gsl_permutation_free(_ind);

  /* Autovalores del tensor de forma */
  forma("pos",Q,&Prop);

  /* Autovalores del tensor de velocidad */
  forma("vel",Q,&Prop);
  
  free(Q);
}

void forma(char *flag,struct particle_data *Q, struct propiedades_st *PropIN){
  float I[3][3];
  int   i,j,k,kk;
  double data[9];
  gsl_matrix *evec;
  struct propiedades_st Prop = *PropIN;

  evec = gsl_matrix_alloc(3, 3);

  if( flag == "pos"){
    Prop.aa = 0.0;
    Prop.cc = 0.0;
    Prop.bb = 0.0;
  }
  if( flag == "vel"){
    Prop.aa_vel = 0.0;
    Prop.cc_vel = 0.0;
    Prop.bb_vel = 0.0;
  }

  /* Calc. 'algo asi como' el Tensor de Inercia */
  kk = 0;
  for(i = 0; i < 3; i++)
    for(j = 0; j < 3; j++){
      I[i][j] = 0. ;
      for(k = 0; k < Prop.npart ; k++){
        if( flag == "pos")
          I[i][j] += ((Q[k].Pos[i]-Prop.pcm[i]) * (Q[k].Pos[j]-Prop.pcm[j]));
        if( flag == "vel")
          I[i][j] += ((Q[k].Vel[i]-Prop.vcm[i]) * (Q[k].Vel[j]-Prop.vcm[j])) ;
      }
      data[kk]=I[i][j]/(double)Prop.npart;
      kk++;
    }

  gsl_matrix_view m = gsl_matrix_view_array(data, 3, 3);
  gsl_vector *eval = gsl_vector_alloc (3);
  gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (3);
  gsl_eigen_symmv (&m.matrix, eval, evec, w);
  gsl_eigen_symmv_free (w);
  gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);

  if( flag == "pos"){
    Prop.cc = (float) gsl_vector_get (eval, 0);
    Prop.bb = (float) gsl_vector_get (eval, 1);
    Prop.aa = (float) gsl_vector_get (eval, 2);
    Prop.evec = evec;
  }
  if( flag == "vel"){
    Prop.cc_vel = (float) gsl_vector_get (eval, 0);
    Prop.bb_vel = (float) gsl_vector_get (eval, 1);
    Prop.aa_vel = (float) gsl_vector_get (eval, 2);
    Prop.evec_vel = evec;
  }

  gsl_vector_free (eval);
}

void write_properties(FILE *pfout, struct propiedades_st *PropIN){
  int   jj;
  float tmp[3];
  struct propiedades_st Prop = *PropIN;

  fwrite(&Prop.npart,sizeof(int),1,pfout);
  fwrite(&Prop.pcm,sizeof(float),3,pfout);
  fwrite(&Prop.vcm,sizeof(float),3,pfout);
  fwrite(&Prop.mostbound,sizeof(float),3,pfout);
  fwrite(&Prop.sig,sizeof(float),3,pfout);
  fwrite(&Prop.L,sizeof(float),3,pfout);
  fwrite(&Prop.lambda,sizeof(float),1,pfout);
  fwrite(&Prop.m200,sizeof(float),1,pfout);
  fwrite(&Prop.r200,sizeof(float),1,pfout);
  fwrite(&Prop.v200,sizeof(float),1,pfout);
  fwrite(&Prop.mvir,sizeof(float),1,pfout);
  fwrite(&Prop.rvir,sizeof(float),1,pfout);
  fwrite(&Prop.vvir,sizeof(float),1,pfout);
  fwrite(&Prop.vmax,sizeof(float),1,pfout);
  fwrite(&Prop.Ep,sizeof(float),1,pfout);
  fwrite(&Prop.Ec,sizeof(float),1,pfout);
  fwrite(&Prop.aa,sizeof(float),1,pfout);
  fwrite(&Prop.bb,sizeof(float),1,pfout);
  fwrite(&Prop.cc,sizeof(float),1,pfout);
  //for(jj=0;jj<3;jj++){
  //  tmp[0] = (float)gsl_matrix_get(Prop.evec,0,jj);
  //  tmp[1] = (float)gsl_matrix_get(Prop.evec,1,jj);
  //  tmp[2] = (float)gsl_matrix_get(Prop.evec,2,jj);
  //  fwrite(&tmp,sizeof(float),3,pfout);
  //}
  fwrite(&Prop.aa_vel,sizeof(float),1,pfout);
  fwrite(&Prop.bb_vel,sizeof(float),1,pfout);
  fwrite(&Prop.cc_vel,sizeof(float),1,pfout);
  //for(jj=0;jj<3;jj++){
  //  tmp[0] = (float)gsl_matrix_get(Prop.evec_vel,0,jj);
  //  tmp[1] = (float)gsl_matrix_get(Prop.evec_vel,1,jj);
  //  tmp[2] = (float)gsl_matrix_get(Prop.evec_vel,2,jj);
  //  fwrite(&tmp,sizeof(float),3,pfout);
  //}
}

void read_properties(FILE *pfin,struct propiedades_st Prop){
  int   jj;
  double tmp[3];

  size_t sizei = sizeof(Prop.npart);
  size_t sizer = sizeof(float);
 
  fread(&Prop.npart,sizei,1,pfin);
  fread(&Prop.pcm,size_real,3,pfin);
  fread(&Prop.vcm,size_real,3,pfin);
  fread(&Prop.mostbound,size_real,3,pfin);
  fread(&Prop.sig,size_real,3,pfin);
  fread(&Prop.L,size_real,3,pfin);
  fread(&Prop.lambda,size_real,1,pfin);
  fread(&Prop.m200,size_real,1,pfin);
  fread(&Prop.r200,size_real,1,pfin);
  fread(&Prop.v200,size_real,1,pfin);
  fread(&Prop.mvir,size_real,1,pfin);
  fread(&Prop.rvir,size_real,1,pfin);
  fread(&Prop.vvir,size_real,1,pfin);
  fread(&Prop.vmax,size_real,1,pfin);
  fread(&Prop.Ep,size_real,1,pfin);
  fread(&Prop.Ec,size_real,1,pfin);
  fread(&Prop.aa,sizer,1,pfin);
  fread(&Prop.bb,sizer,1,pfin);
  fread(&Prop.cc,sizer,1,pfin);
  for(jj = 0; jj < 3; jj++){
    fread(&tmp,sizer,3,pfin);
  }
  fread(&Prop.aa_vel,sizer,1,pfin);
  fread(&Prop.bb_vel,sizer,1,pfin);
  fread(&Prop.cc_vel,sizer,1,pfin);
  for(jj = 0; jj < 3; jj++){
    fread(&tmp,sizer,3,pfin);
  }
}
