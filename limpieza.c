#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <stdbool.h>
#include <string.h>
#include "variables.h"
#include "cosmoparam.h"
#include "limpieza.h"
#include "timer.h"
#include "octree.h"
#include "bitmask.h"

#define EMIN 0.0

int n_free_particles;

#ifdef IDENSUB
void make_cleaning(unsigned int ngrupos, unsigned int *contador_subgrupo){
  unsigned int i, k, l, grupo;
  unsigned int contador;
  struct elemento *curr;
  bool *test;

  #ifdef DEBUG
  int j, l1, l2, l3;
  #endif

  fprintf(stdout,"Clean up begins.... n_grupos_sub %d\n",n_grupos_sub);
  test = (bool *) calloc(ngrupos,sizeof(bool));

  for(i = 0; i < cp.npart; i++)
    P[i].sub = 0;

  contador = 0; 
  /*Recorre cada uno de los grupos identificados en el paso ifrac-1
    limpiando por energia cada uno de sus subgrupos*/
  for(k = 1; k < n_grupos_sub; k++){
    if(sub[k].np == 0 || sub[k].llirst == -1) continue;

    /*Para cada grupos identificado en el paso ifrac-1 crea una linked list
     *de particulas sin grupo en el paso ifrac, y una linked list de subgrupos
     *en el paso ifrac*/
    head_subgroups = NULL;
    head_free_particles = NULL;

    n_free_particles = 0;    //Contador de particulas libres
    Temp.nsub = 0;           //Contador de subgrupos
    Temp.massive_one = 0;    //Subgrupo mas masivo
    Temp.np_massive_one = 0; //Cantidad de particulas en subgrupo mas masivo

    #ifdef DEBUG
    j= 0; l1 = 0;
    #endif
    l = sub[k].llirst; 
    while(l != -1){

      #ifdef DEBUG
      assert(P[l].sub == 0);
      #endif

      grupo = P[l].gr;
      if(grupo == 0) { //Checkea si la particula es una particula libre 
        n_free_particles++;
        curr = (struct elemento *) malloc(sizeof(struct elemento));
        curr->id = l;
        curr->next = head_free_particles;
        head_free_particles = curr;
      }
      else if(!TestBit(test,grupo)) { //Checkear que este subgrupo no este en la lista

        #ifdef DEBUG
        assert(Temp.npgrup[grupo] >= NPARTMIN);
        #endif
        SetBit(test,grupo);

        //Agrega este subgrupo a la lista de subgrupos
        curr = (struct elemento *) malloc(sizeof(struct elemento));
        curr->id = grupo;
        curr->next = head_subgroups;
        head_subgroups = curr;

        Temp.nsub++; //Contador de subgrupos

        #ifdef DEBUG
        l3 = 0;
        l2 = Temp.head[grupo];
        while(l2 != -1)
        {
          l1++;
          l2 = Temp.ll[l2];
          l3++;
        }
        assert(l3 == Temp.npgrup[grupo]);
        #endif

        //Checkea si es el grupo mas masivo
        if(Temp.npgrup[grupo] > Temp.np_massive_one) {
          Temp.np_massive_one = Temp.npgrup[grupo];
          Temp.massive_one = grupo;
        }
      }

      #ifdef DEBUG
      j++;
      #endif

      l = P[l].llsub;
    }

    #ifdef DEBUG
    assert(j == sub[k].np);
    assert(Temp.massive_one != 0 || Temp.nsub == 0);
    assert((l1+n_free_particles) == j);
    #endif

    /*** Una sola subestructura ***/
    if(Temp.nsub <= 1){
      contador++;
      only_one_substructure(k,contador);
      free_memory();
      continue;
    }
    /******************************/

    #ifdef ASSIGN_CLOSEST_GROUP
    /* Calcula velocidades y posicion de los subhalos con mas de NPARTMIN part */
    Temp.vcm    = (my_real **) malloc(3*sizeof(my_real*));
    Temp.vcm[0] = (my_real *)  malloc(ngrupos*sizeof(my_real));
    Temp.vcm[1] = (my_real *)  malloc(ngrupos*sizeof(my_real));
    Temp.vcm[2] = (my_real *)  malloc(ngrupos*sizeof(my_real));
    Temp.pcm    = (my_real **) malloc(3*sizeof(my_real*));
    Temp.pcm[0] = (my_real *)  malloc(ngrupos*sizeof(my_real));
    Temp.pcm[1] = (my_real *)  malloc(ngrupos*sizeof(my_real));
    Temp.pcm[2] = (my_real *)  malloc(ngrupos*sizeof(my_real));

    compute_velocity_position();
    #endif

    #ifdef REASIGNA
    #ifdef ASSIGN_CLOSEST_GROUP
    reasigna_closest();
    #else
    /* Reasigna particulas sin grupo */
    reasigna();
    #endif
    #endif

    #ifdef DEBUG
    //assert(Temp.npgrup[Temp.massive_one] == Temp.np_massive_one + n_free_particles);
    #endif

    Temp.np_massive_one = Temp.npgrup[Temp.massive_one];

    /* Limpia subgrupos excepto el mas masivo */
    curr = head_subgroups;
    while(curr != NULL){
      i = curr->id;

      if(i != Temp.massive_one)
        limpieza_new(i,Temp.massive_one);

      curr = curr->next;
    }
    /*Termino de limpiar los subgrupos excepto el mas masivo*/

    /*** Despues de la limpieza queda solo una subestructura ***/
    if(Temp.nsub <= 1){
      contador++;
      only_one_substructure(k,contador);
      free_memory();
      continue;
    }
    /***********************************************************/

    limpieza_new(Temp.massive_one,0);

    curr = head_subgroups;
    while(curr){
      i = curr->id;
      if(Temp.head[i] == -1){
        assert(Temp.npgrup[i] == 0);
        curr = curr->next;
        continue;
      }

      contador++;

      #ifdef DEBUG
      j = 0;
      #endif
      l = Temp.head[i];
      while(l != -1){
        P[l].sub = contador;
        l = Temp.ll[l];
        #ifdef DEBUG
        j++;
        #endif
      }
      #ifdef DEBUG
      assert(j == Temp.npgrup[i]);
      #endif

      curr = curr->next;
    }

    free_memory();
  }

  *contador_subgrupo = contador;

  #ifdef ASSIGN_CLOSEST_GROUP
  free(Temp.pcm[0]); free(Temp.pcm[1]); free(Temp.pcm[2]); free(Temp.pcm);
  free(Temp.vcm[0]); free(Temp.vcm[1]); free(Temp.vcm[2]); free(Temp.vcm);
  #endif

  free(test);
  printf("End cleaning. Subgroups found = %u \n",contador);
}
#endif

#ifdef ASSIGN_CLOSEST_GROUP
void compute_velocity_position(){
  int ig,ip,idim;

  struct elemento *curr;
 
  curr = head_subgroups;
  while(curr != NULL){
    ig = curr->id;

    for(idim = 0; idim < 3; idim++){
      Temp.vcm[idim][ig] = 0.0;
      Temp.pcm[idim][ig] = 0.0;
    }

    ip = Temp.head[ig];
    while(ip != -1){
      for(idim = 0; idim < 3; idim++){
        Temp.vcm[idim][ig] += P[ip].Vel[idim];
        Temp.pcm[idim][ig] += P[ip].Pos[idim];
      }
      ip = Temp.ll[ip];
    }

    for(idim = 0; idim < 3; idim++){
      Temp.vcm[idim][ig] /= (my_real)Temp.npgrup[ig];
      Temp.pcm[idim][ig] /= (my_real)Temp.npgrup[ig];
    }

    curr = curr->next;
  }
}
#endif

#ifdef IDENSUB
void free_memory(void){
  struct elemento *curr, *puntero;

  curr = head_free_particles;
  while(curr){
    puntero = curr->next;
    free(curr);
    curr = puntero;
  }

  curr = head_subgroups;
  while(curr){
    puntero = curr->next;
    free(curr);
    curr = puntero;
  }
}
#endif

#ifdef ASSIGN_CLOSEST_GROUP
void reasigna_closest(void){
  int ig,ip,idim;
  int destino;
  struct elemento *curr, *curr_gr;
  my_real dp[3], dv[3], dis, E, Ep, Ec;

  curr = head_free_particles;
  while(curr != NULL){
    ip = curr->id;

    E = 1.E26;
    destino = Temp.massive_one;
  
    curr_gr = head_subgroups;
    while(curr_gr != NULL){
      ig = curr_gr->id;

      dis = 0.0;
      Ec = 0.0;
      for(idim = 0 ; idim < 3; idim++){
        dp[idim] = Temp.pcm[idim][ig] - P[ip].Pos[idim];
        dis += dp[idim]*dp[idim];

        dv[idim]  = cp.Hubble_a;
        dv[idim] *= cp.aexp;
        dv[idim] *= (double)dp[idim]*0.001;  /* Velocidad de Hubble */
        dv[idim] += sqrt(cp.aexp)*(double)(Temp.vcm[idim][ig] - P[ip].Vel[idim]);
        Ec += dv[idim]*dv[idim];
      }
	    dis  = sqrt(dis);
	  	dis += (my_real)cp.soft;

      Ec *= 0.5;

	  	Ep = (my_real)Temp.npgrup[ig]/dis;
	    Ep += (1./cp.soft);   /* Autoenergia */
	    Ep *= (GCONS*cp.Mpart*Msol*1.E10/Kpc/cp.aexp);
	  	Ep *= (-1.);          /* Cambio de signo para que Ep sea negativa */

      //Ep = dis;

      //if((Ep+Ec) < E)
      //{
      //  E = (Ep+Ec);
      //  destino = ig;
      //}

      if(Ep < E){
        E = Ep;
        destino = ig;
      }

      curr_gr = curr_gr->next;
    }

    #ifdef DEBUG
    assert(destino != 0);
    assert(P[ip].gr == 0);
    #endif
  
    P[ip].gr = destino;
    Temp.ll[ip] = Temp.head[destino];
    Temp.head[destino] = ip;
    Temp.npgrup[destino]++;

    curr = curr->next;
  }
}
#endif

#ifndef  ASSIGN_CLOSEST_GROUP
/*reasigna particulas libres al grupo mas masivo*/
void reasigna(void){
  int k;
  int destino;
  struct elemento *curr;

  destino = Temp.massive_one;

  curr = head_free_particles;
  while(curr != NULL){
    k = curr->id;

    P[k].gr = destino;
    Temp.ll[k] = Temp.head[destino];
    Temp.head[destino] = k;
    Temp.npgrup[destino]++;

    curr = curr->next;
  }
}
#endif

void limpieza_new(my_int ig, my_int destino){
  my_int i,j,k,n_unbound,temp,iEpmin;
  my_int *lista;
  int    dim;
  double pcm[3],vcm[3];
  double dv[3];
  double E;
  struct particle_data *Q;
  my_int npart = Temp.npgrup[ig];
  //int destino = Temp.massive_one;
  my_int *indices;

  my_real xmin[3], xmax[3];
  xmin[0] = 1.E26; xmax[0] = -1.E26;
  xmin[1] = 1.E26; xmax[1] = -1.E26;
  xmin[2] = 1.E26; xmax[2] = -1.E26;

  printf("%u %u %u\n",ig,destino,npart);

  Q	= (struct particle_data *) malloc(npart*sizeof(struct particle_data));
  assert(Q != NULL);
  indices = (my_int *) malloc(npart*sizeof(my_int));
  assert(indices != NULL);

  j = 0;
  k = Temp.head[ig];
  do{
    Q[j].Pos[0] = P[k].Pos[0];
		Q[j].Pos[1] = P[k].Pos[1];
		Q[j].Pos[2] = P[k].Pos[2];
    Q[j].Vel[0] = P[k].Vel[0];
		Q[j].Vel[1] = P[k].Vel[1];
		Q[j].Vel[2] = P[k].Vel[2];
    Q[j].gr     = P[k].gr;
    indices[j]  = k;
    Q[j].Ep     = 0.0;
    Q[j].Ec     = 0.0;

    if(Q[j].Pos[0] < xmin[0]) xmin[0] = Q[j].Pos[0];
    if(Q[j].Pos[0] > xmax[0]) xmax[0] = Q[j].Pos[0];
    if(Q[j].Pos[1] < xmin[1]) xmin[1] = Q[j].Pos[1];
    if(Q[j].Pos[1] > xmax[1]) xmax[1] = Q[j].Pos[1];
    if(Q[j].Pos[2] < xmin[2]) xmin[2] = Q[j].Pos[2];
    if(Q[j].Pos[2] > xmax[2]) xmax[2] = Q[j].Pos[2];

    j++;
    k = Temp.ll[k];
  }while(k != Temp.head[ig]);

  #ifdef DEBUG
  assert(j == npart);
  #endif
	for(i = 0; i < npart; i++){
    Q[i].Pos[0] -= xmin[0];
    Q[i].Pos[1] -= xmin[1];
    Q[i].Pos[2] -= xmin[2];
  }

  pcm[0] = 0.0;
  pcm[1] = 0.0;
  pcm[2] = 0.0;
  vcm[0] = 0.0;
  vcm[1] = 0.0;
  vcm[2] = 0.0;

	for(i = 0; i < npart; i++)
    for(dim = 0; dim < 3; dim++){
      pcm[dim] += (double)Q[i].Pos[dim];
      vcm[dim] += (double)Q[i].Vel[dim];
    }

  for(dim = 0; dim < 3; dim++){
    pcm[dim] /= (double)npart;
    vcm[dim] /= (double)npart;
  }

  compute_potential_energy_subgrupo(npart,Q,&iEpmin);

  pcm[0] = Q[iEpmin].Pos[0];
  pcm[1] = Q[iEpmin].Pos[1];
  pcm[2] = Q[iEpmin].Pos[2];

	for(i = 0; i < npart; i++)
    for(dim = 0; dim < 3; dim++){
      Q[i].Pos[dim] -= (my_real)pcm[dim];
      Q[i].Vel[dim] -= (my_real)vcm[dim];
    }

  lista = (my_int *) calloc(npart,sizeof(my_int));
  assert(lista != NULL);

  n_unbound = 0;
	for(i = 0; i < npart; i++){
    Q[i].Ec = 0.0;
    /* Calcula la energia cinetica (velocidades en [km / seg]) */
    for(dim = 0 ; dim < 3; dim++){
      dv[dim]  = cp.Hubble_a;
      dv[dim] *= cp.aexp;
      dv[dim] *= (double)Q[i].Pos[dim]*0.001;  /* Velocidad de Hubble */
      dv[dim] += sqrt(cp.aexp)*(double)Q[i].Vel[dim];
      Q[i].Ec += dv[dim]*dv[dim];
    }
    
    Q[i].Ec *= 0.5;
    /***********************************************************/

    /* Calcula Energia Total */
    E = Q[i].Ec + Q[i].Ep;
    /*************************/

    if(E > EMIN){
      lista[n_unbound] = i;
      n_unbound++;
    }
  }

  lista = (my_int *) realloc(lista,n_unbound*sizeof(my_int));

  /** Si luego de limpiar el grupo se queda con menos de NPARTMIN
      lo disolvemos totalmente, y todas su particulas se las damos
      al grupo mas masivo **/
  if((npart - n_unbound) < NPARTMIN){
    Temp.nsub--;

    #ifdef DEBUG
    j = 0;
    #endif
    k = Temp.head[ig];
    do{
      temp = Temp.ll[k];
      #ifdef DEBUG
      j++;
      #endif
      P[k].gr = destino;
//      Temp.ll[k] = Temp.head[destino];
//      Temp.head[destino] = k;
//      Temp.npgrup[destino]++;

      k = temp;
    }while(k != Temp.head[ig]);

    #ifdef DEBUG
    assert(j == npart);
    #endif

    Temp.npgrup[ig] = 0;
    Temp.head[ig]   = 0;
  }else if(n_unbound > 0)
  /** Si se queda con mas de NPARTMIN entonces las no-ligadas
      se las damos al grupo mas masivo, y luego reconstruimos
      la linked list para este grupo **/
  {
    for(i = 0; i < n_unbound; i++){
      k = lista[i];
      k = indices[k];
      P[k].gr = destino;
      //Temp.ll[k] = Temp.head[destino];
      //Temp.head[destino] = k;
      //Temp.npgrup[destino]++;
    }

    /***Regenera la linked list para el subgrupo ig******/
    //Temp.head[ig]  = -1;
    Temp.npgrup[ig] = 0;

    j = 0;
    for(i = 0; i < npart; i++){
      if(indices[i] == lista[j]){
        j++;
        continue;
      }

      k = indices[i];
      #ifdef DEBUG
      assert(P[k].gr == ig);
      #endif
      Temp.head[ig] = k;
      Temp.npgrup[ig]++;
    }

    j = 0;
    for(i = 0; i < npart; i++){
      if(indices[i] == lista[j]){
        j++;
        continue;
      }

      k = indices[i];
      #ifdef DEBUG
      assert(P[k].gr == ig);
      #endif

      Temp.ll[k] = Temp.head[ig];
      Temp.head[ig] = k;
      //Temp.npgrup[ig]++;
    }

    #ifdef DEBUG
    assert(Temp.npgrup[ig] >= NPARTMIN);
    #endif
    /****************************************************/
  }

  free(Q);
  free(lista);
}

void compute_potential_energy_subgrupo(my_int npart, struct particle_data *Q, my_int *iEpmin){
  my_int i,j,dim,iEpmin_local = -1;
  float  Theta = 0.45;
  double dx[3],dis,Epmin;

	/* Si el halo tiene menos de 1000 particulas calcula la energia
	de forma directa (N^2). Si tiene mas de 1000 particulas usa 
	un octree. */
  if(npart < 1000){
    #pragma omp parallel for default(none) \
                shared(npart,Q,cp)         \
                private(i,j,dim,dx,dis)    \
                schedule(dynamic)
	  for(i = 0; i < npart; i++){
      /* Calcula energia potencial */
	  	Q[i].Ep  = 0.;

	    for(j = 0; j < npart; j++){
	  	  if(j == i)continue;

	  		for(dim = 0; dim < 3; dim++)
	  			dx[dim] = (double)Q[i].Pos[dim] - (double)Q[j].Pos[dim];

	  	  dis  = dx[0]*dx[0];
	  	  dis += dx[1]*dx[1];
	  	  dis += dx[2]*dx[2];
	  		dis  = sqrt(dis);
	  		dis += (double)cp.soft;

	  	  Q[i].Ep += 1.0/dis;
	    }

	    Q[i].Ep += (1./cp.soft);   /* Autoenergia */
	    Q[i].Ep *= (GCONS*cp.Mpart*Msol*1.E10/Kpc/cp.aexp);
	  	Q[i].Ep *= (-1.);          /* Cambio de signo para que Ep sea negativa */
	  }
  }else{
    force_treeallocate(2 * npart + 100);
    force_treebuild(npart,Q,Theta);

    #pragma omp parallel for default(none) \
                shared(npart,Q,cp)         \
                private(i)                 \
                schedule(dynamic)
    for(i = 0; i < npart; i++){
      Q[i].Ep = 0.0;
      force_treeevaluate_potential(&Q[i].Pos[0], &Q[i].Ep);
       
      Q[i].Ep -= cp.Mpart/cp.soft;    /* Autoenergia */
      Q[i].Ep *= GCONS/cp.aexp*Msol/Kpc;
      Q[i].Ep *= 1.0e10;              /* Para que tenga unidades de Ecin */
    }

    force_treefree();
  }

  Epmin = 1.E26;
  for(i = 0; i < npart; i++){
    if(Q[i].Ep < Epmin){
      Epmin = Q[i].Ep;
      iEpmin_local = i;
    }
  }

  *iEpmin = iEpmin_local;
}

#ifdef COMPUTE_EP
void compute_potential_energy(void){
  int i;
  float Theta = 0.45;

  fprintf(stdout,"Computting potential energy...\n");
  fflush(stdout);

  force_treeallocate(2 * np_in_fof + 100);
  force_treebuild(np_in_fof,P,Theta);

  for(i = 0; i < np_in_fof; i++){
    /* Calcula la energia potencial */
    P[i].Ep = 0.0;
    force_treeevaluate_potential(&P[i].Pos[0], &P[i].Ep);
     
    P[i].Ep -= cp.Mpart/cp.soft;    /* Autoenergia */
    P[i].Ep *= GCONS/cp.aexp*Msol/Kpc;
    P[i].Ep *= 1.0e10;              /* Para que tenga unidades de Ecin */
  }

  force_treefree();
}
#endif

#ifdef IDENSUB
void only_one_substructure(int k, int contador)
{
  int l;
  #ifdef DEBUG
  int j = 0;
  #endif
  l = sub[k].llirst;
  while(l != -1)
  {
    #ifdef DEBUG
    j++;
    #endif
    P[l].sub = contador;
    l = P[l].llsub;
  }
  #ifdef DEBUG
  assert(j == sub[k].np);
  #endif
}
#endif
