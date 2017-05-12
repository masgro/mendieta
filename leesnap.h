#ifndef LEESNAP_H
#define LEESNAP_H

void leeheader(char *nombrefile);
void lee(char *filename, struct particle_data *Q, int *ind);
void change_positions(int n);
void re_change_positions(int n, struct particle_data *Q);
void read_gadget();
void select_particles(void);


static type_real pmin[3], pmax[3];

/*Input and output files*/
struct SnapST{
  int nfiles;
  char root[200], name[200]; 
} snap;

struct io_header{
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

struct box_st{
  type_real cen[3];
  type_real max[3];
  type_real min[3];
  type_real franja;
  type_real lado;
} box;
#endif
