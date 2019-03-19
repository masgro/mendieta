#ifndef LEESNAP_H
#define LEESNAP_H

void leeheader(char *nombrefile);
void lee(char *filename, struct particle_data *Q, my_int *ind);
void change_positions(my_int n);
void re_change_positions(my_int n, struct particle_data *Q);
void read_gadget();
void select_particles(void);

my_real pmin[3], pmax[3];

/*Input and output files*/
struct SnapST{
  int nfiles;
  char root[200], name[200]; 
} snap;

struct io_header{
  unsigned int npart[6];
  double       mass[6];
  double       time;
  double       redshift;
  int          flag_sfr;
  int          flag_feedback;
  unsigned int npartTotal[6];
  int          flag_cooling;
  int          num_files;
  double       BoxSize;
  double       Omega0;
  double       OmegaLambda;
  double       HubbleParam; 
  int          flag_stellarage;
  int          flag_metals;
  unsigned int npartTotal_HW[6];
  char         fill[256-6*4-6*8-2*8-2*4-6*4-2*4-4*8-2*4-6*4];  /* fills to 256 Bytes */
} header;

struct box_st{
  my_real cen[3];
  my_real max[3];
  my_real min[3];
  my_real franja;
  my_real lado;
} box;
#endif
