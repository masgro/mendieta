#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <stdbool.h>
#include <string.h>
#include "cosmoparam.h"
#include "timer.h"

#define KERNEL_LENGTH   10000
#define KERNEL_TABLE 10000

static struct NODE *nodes;

static struct NODE *last;

static int    numnodestotal;        /* total number of nodes */
static int    MaxNodes;

static float  xmin[3],xmax[3];

static  int    N;

static float   knlrad  [KERNEL_LENGTH+1],
               knlpot  [KERNEL_LENGTH+1];


/* Posiciones, velocidades y energias de las partículas */
struct particle_data 
{
  float  Pos[3];
  int   id;
  double     Ep;
} *P;

struct SnapST
{
  int nfiles;
  char root[200], name[200]; 
} snap;

struct io_header
{
  int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
} header;

struct NODE 
{ 
  float s[3];                         /* center of mass */
  int    partind;
  float center[3],len;                /* center and sidelength of treecubes */
  float mass;                          /* mass*/
  float oc;                           /* variable for opening criter*/
  struct NODE *next,*sibling,*father,*suns[8];
};

void force_treeallocate(int maxnodes);
void force_treefree(void);

void add_particle_props_to_node(struct NODE *no, int p);

int  force_treebuild(int np, float thetamax);
void force_setupnonrecursive(struct NODE *no);
void force_treeevaluate_potential(float *pos, double *pot);

void force_setkernel(void);

void compute_potential_energy(void);
void read_gadget();
void leeheader(char *filename);
void lee(char *filename, struct particle_data *Q, int *ind);

int main(int argc, char **argv)
{
  int    i,l;
  int    init_ifrac, ifrac, nfrac;
  double frac, start, end;
  char   filename[200];
  FILE   *pfin, *pfout;
  char   fof_file[200], sub_file[200];

  TIMER(start);

  sprintf(filename,"%s",argv[1]);
  pfin = fopen(filename,"r");
  fscanf(pfin,"%d  \n",&snap.nfiles);
  fscanf(pfin,"%s  \n",snap.root);
  fscanf(pfin,"%s  \n",snap.name);
  fscanf(pfin,"%s  \n",fof_file);
  fclose(pfin);

  /* Lee archivos de la simulacion */
  read_gadget();

  compute_potential_energy();

  free(P);

  TIMER(end);
  printf("Total time %f\n",end-start);

  return(EXIT_SUCCESS);
}

void compute_potential_energy(void)
{
  int i;
  float  Theta = 0.45;

  fprintf(stdout,"Computting potential energy...\n");
  fflush(stdout);

  force_treeallocate(2 * cp.npart + 100);
  force_treebuild(cp.npart, Theta);

  for( i = 0 ; i < cp.npart ; i++ )
  {
    /* Calcula la energia potencial */
    P[i].Ep = 0.0;
    force_treeevaluate_potential(&P[i].Pos[0], &P[i].Ep);
     
    P[i].Ep -= cp.Mpart/cp.soft;    /* Autoenergia */
    P[i].Ep *= GCONS/cp.aexp*Msol/Kpc;
    P[i].Ep *= 1.0e10;              /* Para que tenga unidades de Ecin */
  }
}

void read_gadget()
{
  char filename[200];
  int  ifile,ind;

  if(snap.nfiles>1)
    sprintf(filename,"%s%s.0",snap.root,snap.name);
  else
    sprintf(filename,"%s%s",snap.root,snap.name);

  leeheader(filename);

  /****** ALLOCATACION TEMPORAL DE LAS PARTICULAS ****************/
  printf("Allocating %zu Gb for %d particles\n",
              cp.npart*sizeof(struct particle_data)/1024/1024/1024,cp.npart);
  P = (struct particle_data *) malloc(cp.npart*sizeof(struct particle_data));
  assert(P != NULL);
  /***** LEE POS Y VEL DE LAS PARTICULAS ***********************/
  for(ifile = 0, ind = 0; ifile < snap.nfiles; ifile++)
  {
    if(snap.nfiles>1)
      sprintf(filename,"%s%s.%d",snap.root,snap.name,ifile);
    else
      sprintf(filename,"%s%s",snap.root,snap.name);

    lee(filename,P,&ind);
  }

  fprintf(stdout,"Termino de leer\n"); fflush(stdout);
}

void leeheader(char *filename)
{
  FILE *pf;
  int d1,d2;

  pf = fopen(filename,"r");
  if(pf == NULL)
  {
    fprintf(stderr,"can't open file `%s`\n",filename);
    exit(EXIT_FAILURE);
  }

  fread(&d1, sizeof(d1), 1, pf);
  fread(&header, sizeof(header), 1, pf);
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);
  fclose(pf);

  /* Definicion estructura cosmoparam */
  cp.omegam    = header.Omega0;
  cp.omegal    = header.OmegaLambda;
  cp.omegak    = 1.0 - cp.omegam - cp.omegal;
  cp.hparam    = header.HubbleParam;
  cp.lbox      = header.BoxSize;
  cp.npart     = header.npartTotal[1];
  cp.Mpart     = header.mass[1];
  cp.redshift  = header.redshift;
  cp.aexp      = ( 1.0 / ( 1.0 + cp.redshift ) );
  cp.Hubble_a  = cp.omegam/cp.aexp/cp.aexp/cp.aexp;
  cp.Hubble_a += cp.omegak/cp.aexp/cp.aexp;
  cp.Hubble_a += cp.omegal;
  cp.Hubble_a  = sqrt(cp.Hubble_a);
  cp.Hubble_a *= 100.0*cp.hparam;

  printf("*********************************** \n");
  printf("*   Parametros de la simulacion   * \n");
  printf("*********************************** \n");
  printf("  Numero de particulas = %d \n", cp.npart);
  printf("  Lado del box = %g \n", cp.lbox);
  printf("  Redshift = %g \n", cp.redshift);
  printf("  Omega Materia = %g \n", cp.omegam);
  printf("  Omega Lambda = %g \n", cp.omegal);
  printf("  Parametro de Hubble = %g \n",cp.hparam);
  printf("  Masa por particula = %g \n",cp.Mpart);
  printf("  Softening = %g\n",cp.soft);
  printf("*********************************** \n");
  printf("*********************************** \n");

}

void lee(char *filename, struct particle_data *Q, int *ind)
{
  FILE *pf;
  int d1, d2;
  int k, pc, n;

  float r[3],v[3];
  int id;
 
  pf = fopen(filename,"r");
  if(pf == NULL)
  {
    fprintf(stderr,"can't open file `%s`\n",filename);
    exit(EXIT_FAILURE);
  }

  fprintf(stdout,"Reading file: %s \n",filename); fflush(stdout);

  fread(&d1, sizeof(d1), 1, pf);
  fread(&header, sizeof(header), 1, pf);
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);

  fread(&d1, sizeof(d1), 1, pf);
  for(k = 0, pc = 0; k < 6; k++)
  {
    for(n = 0; n < header.npart[k]; n++)
    {
      fread(&r[0], sizeof(float), 3, pf);
      if(k == 1)
      {
        Q[*ind+pc].Pos[0] = r[0];
        Q[*ind+pc].Pos[1] = r[1];
        Q[*ind+pc].Pos[2] = r[2];
        pc++;
      }
    }
  }
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);

  fread(&d1, sizeof(d1), 1, pf);
#if defined STORE_VELOCITIES
  for(k = 0, pc = 0; k < 6; k++)
  {
    for(n = 0; n < header.npart[k]; n++)
    {
      fread(&v[0], sizeof(float), 3, pf);
      if(k == 1)
      {
        Q[*ind+pc].Vel[0] = v[0]*VELFACTOR;
        Q[*ind+pc].Vel[1] = v[1]*VELFACTOR;
        Q[*ind+pc].Vel[2] = v[2]*VELFACTOR;
        pc++;
      }
    }
  }
#else
  fseek(pf,d1,SEEK_CUR);
#endif
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);

  fread(&d1, sizeof(d1), 1, pf);
#if defined STORE_IDS
  for(k = 0, pc = 0; k < 6; k++)
  {
    for(n = 0; n < header.npart[k]; n++)
    {
      fread(&id, size_int, 1, pf);
      if(k == 1)
      {
        Q[*ind+pc].id = id;
        pc++;
      }
    }
  }
#else
  fseek(pf,d1,SEEK_CUR);
#endif
  fread(&d2, sizeof(d2), 1, pf);
  assert(d1==d2);

  *ind += pc;
  
  fclose(pf);
}


/*************************************************************************/
/***ALLOCATACION Y FREE***/
void force_treeallocate(int maxnodes) 
{
  MaxNodes = maxnodes;
  nodes = (struct NODE *) malloc(MaxNodes*sizeof(struct NODE));
  assert(nodes!=NULL);
  force_setkernel();
}

void force_treefree(void)
{
  free(nodes);
}
/************************************************************************/

void add_particle_props_to_node(struct NODE *no, int p)
{
  int i;
  for( i = 0 ; i < 3 ; i++)
    no->s[i] += (float)cp.Mpart * (P[p].Pos[i] - no->center[i]);

  no->mass += (float)cp.Mpart;
}

/* packs the particles of group 'gr' into into BH-trees */
int force_treebuild(int np, float thetamax)
{
  int    i,j;
  int    subp,subi,p,subnode,fak;
  float  length;
  float  dx,dy,dz;
  struct NODE *nfree,*th,*nn,*ff;

  N = np;

  nfree = nodes;
  numnodestotal = 0;

  for(j = 0 ; j < 3 ; j++)                        /* find enclosing rectangle */
    xmin[j] = xmax[j] = P[0].Pos[j];

  for(i = 0 ; i < N ; i++)
    for(j = 0 ; j < 3 ; j++)
    {
      if(P[i].Pos[j] > xmax[j]) 
        xmax[j] = P[i].Pos[j];
      if(P[i].Pos[j] < xmin[j]) 
        xmin[j] = P[i].Pos[j];
    }
  
  for(j = 1 , length = xmax[0]-xmin[0] ; j < 3 ; j++)  /* determine maxmimum extension */
    if((xmax[j]-xmin[j]) > length)
      length = xmax[j]-xmin[j];
  
  length *= 1.001;

  /* insert first particle in root node */
  for(j = 0 ; j < 3 ; j++)
    nfree->center[j] = (xmax[j]+xmin[j])/2;
  nfree->len = length;
  
  /*inicializa variables*/
  nfree->father = 0;
  for(i = 0 ; i < 8 ; i++)
    nfree->suns[i] = 0;
  nfree->partind = 0;

  nfree->mass = (float)cp.Mpart;

  for(i = 0 ; i < 3 ; i++)
    nfree->s[i] = (float)cp.Mpart*(P[0].Pos[i] - nfree->center[i]);
  
  /*inicializa la variable que apunta al hermano*/
  nfree->sibling = 0;

  numnodestotal++; nfree++;
  
  if(numnodestotal >= MaxNodes)
  {
    fprintf(stderr,"Maximum number %d of tree-nodes reached. file: %s line: %d\n",
            numnodestotal,__FILE__,__LINE__);
    exit(EXIT_FAILURE);
  }

  /* insert all other particles */
  for(i = 1 ; i < N ; i++)
  {
    th = nodes;
  
    while(1)
    {
      add_particle_props_to_node(th, i);

      if(th->partind >= 0)
        break;
    
      for(j = 0 , subnode = 0 , fak = 1 ; j < 3 ; j++ , fak <<= 1)
        if(P[i].Pos[j] > th->center[j])
          subnode += fak;

      nn = th->suns[subnode];
      if(nn != NULL)
        th = nn;
      else
        break;
    }
      
    if(th->partind >= 0)  /* cell is occcupied with one particle */
    {
      while(1)
      {
        p = th->partind;

        for( j = 0 , subp = 0 , fak = 1 ; j < 3 ; j++ , fak <<= 1)
          if(P[p].Pos[j] > th->center[j])
            subp += fak;

        nfree->father = th;
        
        for( j = 0 ; j < 8 ; j++)
          nfree->suns[j] = 0;

        nfree->sibling = 0;
        
        nfree->len = th->len/2;
    
          for(j = 0 ; j < 3 ; j++)
          nfree->center[j] = th->center[j];

        for(j = 0 ; j < 3 ; j++)
          if(P[p].Pos[j] > nfree->center[j])
            nfree->center[j] += nfree->len/2;
          else
            nfree->center[j] -= nfree->len/2;

        nfree->partind = p;

        nfree->mass = (float)cp.Mpart;

        for(j = 0 ; j < 3 ; j++)
          nfree->s[j] = (float)cp.Mpart*(P[p].Pos[j] - nfree->center[j]);
        
        th->partind = -1;
        th->suns[subp] = nfree;
      
        numnodestotal++; nfree++;
        if(numnodestotal >= MaxNodes)
        {
          fprintf(stderr,"Maximum number %d of tree-nodes reached. file: %s line: %d\n",
                  numnodestotal,__FILE__,__LINE__);
          fprintf(stderr,"i=%d\n", i);
          exit(EXIT_FAILURE);
        }

        for(j = 0 , subi = 0 , fak = 1 ; j < 3 ; j++ , fak <<= 1)
          if(P[i].Pos[j] > th->center[j])
            subi += fak;

        if(subi == subp)   /* the new particle lies in the same sub-cube */
        {
          th = nfree-1;
          add_particle_props_to_node(th,i);      
        }
        else
          break;
      }
    }
      
    for(j = 0 , subi = 0 , fak = 1 ; j < 3 ; j++ , fak <<= 1)
      if(P[i].Pos[j] > th->center[j])
        subi += fak;
      
    nfree->father = th;
      
    for(j = 0 ; j < 8 ; j++)
      nfree->suns[j] = 0;

    nfree->sibling = 0;

    nfree->len = th->len/2;

    for(j = 0 ; j < 3 ; j++)
      nfree->center[j] = th->center[j];

    for(j = 0 ; j < 3 ; j++)
      if(P[i].Pos[j] > nfree->center[j])
        nfree->center[j] += nfree->len/2;
      else
        nfree->center[j] -= nfree->len/2;

    nfree->mass = (float)cp.Mpart;

    for(j = 0 ; j < 3 ; j++)
      nfree->s[j] = (float)cp.Mpart*(P[i].Pos[j] - nfree->center[j]);

    nfree->partind = i;
    th->suns[subi] = nfree;
      
    numnodestotal++; nfree++;

    if(numnodestotal >= MaxNodes)
    {
      fprintf(stderr,"Maximum number %d of tree-nodes reached. file: %s line: %d\n",
              numnodestotal,__FILE__,__LINE__);
      fprintf(stderr,"i=%d\n", i);
      exit(EXIT_FAILURE);
    }
  }

  /* now finish-up center-of-mass and quadrupol computation */
  
  for(i = 0, th = nodes; i < numnodestotal; i++, th++)
  {
    for(j = 0; j < 3; j++)
      th->s[j] /= th->mass;
      
    if(th->partind < 0)   /* cell contains more than one particle */
    {
      dx = th->s[0];
      dy = th->s[1];
      dz = th->s[2];
    
      th->oc  = (float)sqrt(dx*dx + dy*dy + dz*dz);
      th->oc += th->len/(thetamax); 
      th->oc *= th->oc;     /* used in cell-opening criterion */
    }

    th->s[0] += th->center[0];
    th->s[1] += th->center[1];
    th->s[2] += th->center[2];
 
    for(j = 7, nn = 0; j >= 0; j--)    /* preparations for non-recursive walk */
    {
      if(th->suns[j])
      {
        th->suns[j]->sibling = nn;
        nn = th->suns[j];
      }
    }
  }

  last = 0;
  force_setupnonrecursive(nodes);    /* set up non-recursive walk */
  last->next = 0;
  
  for(i = 0, th = nodes; i < numnodestotal; i++, th++)
    if(!(th->sibling))
    {
      ff = th;
      nn = ff->sibling;

      while(!nn)
      {
        ff = ff->father;
        if(!ff)
        break;
        nn = ff->sibling;
      }
  
      th->sibling = nn;
    }
  return numnodestotal;
}


void force_setupnonrecursive(struct NODE *no)
{
  int i;
  struct NODE *nn;
  
  if(last)
    last->next = no;

  last = no;
  
  for(i = 0; i < 8; i++)
  {
    nn = no->suns[i];
    if(nn != NULL)
      force_setupnonrecursive(nn);
  }
}


/****************************************************************/
void force_treeevaluate_potential(float *pos, double *pot)
{
  struct NODE *no;
  float r2,dx,dy,dz,r,u,h,ff;
  float wp;
  float h_inv;
  float local_pot;
  int ii;

  h = 2.8*(float)cp.soft;
  h_inv = 1.0 / h;

  local_pot = 0.;

  no = nodes;
  while(no)
  {
    dx = no->s[0];     
    dy = no->s[1];     
    dz = no->s[2];
    dx -= pos[0];     
    dy -= pos[1];     
    dz -= pos[2];

    r2 = dx*dx + dy*dy + dz*dz;

    if(no->partind >= 0)   /* single particle */
    {
      r = (float)sqrt(r2);  
     
      u = r * h_inv;

      if(u >= 1)
      {
        local_pot -= no->mass/r;
      }
      else
      {
        ff = u*KERNEL_LENGTH;
        ii = (int)ff;
        ff -= knlrad[ii]*KERNEL_LENGTH;
        wp = knlpot[ii] + (knlpot[ii+1] - knlpot[ii])*ff;
        
        local_pot += no->mass*h_inv*wp;
      }
      no = no->sibling;
    }
    else
    {
      if(r2 < no->oc)
      {
        no = no->next;  /* open cell */
      }
      else
      {
        r = (float)sqrt(r2);  
        u = r*h_inv;
    
        if(u >= 1)  /* ordinary quadrupol moment */
        {
          local_pot += -no->mass/r;
        }
        else    /* softened monopole moment */
        {
          ff = u*KERNEL_LENGTH;
          ii = (int)ff;
          ff -= knlrad[ii]*KERNEL_LENGTH;
          wp = knlpot[ii] + (knlpot[ii+1]-knlpot[ii])*ff;

          local_pot += no->mass*h_inv*wp;
        }
        no = no->sibling;
      }
    }
  }

  *pot = local_pot;
}


/****************************************/
void force_setkernel(void) 
{
  int i;
  double u;
  double u2;
  double u3;
  double u4;
  double u5;

  for(i = 0; i <= KERNEL_LENGTH; i++)
  {
    u = (double)i/(double)KERNEL_LENGTH;

    knlrad[i] = u;

    if(u <= 0.5)
    {
      u2 = u*u;
      u4 = u2*u2;
      u5 = u4*u;
      knlpot[i]  = 16.0/3.0*u2;
      knlpot[i] -= 48.0/5.0*u4;
      knlpot[i] += 32.0/5.0*u5;
      knlpot[i] -= 14.0/5.0;
    }
    else
    {
      u2 = u*u;
      u3 = u2*u;
      u4 = u2*u2;
      u5 = u4*u;
      knlpot[i]  = 1.0/15.0/u;
      knlpot[i] += 32.0/3.0*u2;
      knlpot[i] -= 16.0*u3;
      knlpot[i] += 48.0/5.0*u4;
      knlpot[i] -= 32.0/15.0*u5;
      knlpot[i] -= 16.0/5.0;
    }
  }
}
