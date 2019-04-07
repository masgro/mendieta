#ifndef OCTREE_H
#define OCTREE_H

//struct Qst{
//  my_real Pos[3];
//  my_real Vel[3];
//  double    Ep, Ec;
//  int       gr;
//  int       indx;
//} *Q;

struct NODE{ 
  my_real s[3];                         /* center of mass */
  int    partind;
  my_real center[3],len;                /* center and sidelength of treecubes */
  float mass;                          /* mass*/
  float oc;                           /* variable for opening criter*/
  struct NODE *next,*sibling,*father,*suns[8];
};

void force_treeallocate(int maxnodes);
void force_treefree(void);

void add_particle_props_to_node(struct NODE *no, struct particle_data *Q, int p);

int  force_treebuild(int np, struct particle_data *Q, float thetamax);
void force_setupnonrecursive(struct NODE *no);
void force_treeevaluate_potential(my_real *pos, double *pot);

void force_setkernel(void);

#endif
