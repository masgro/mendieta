#ifndef LEESNAP_H
#define LEESNAP_H

#ifndef VARIABLES_H
  #include "variables.h"
#endif

struct io_header
{
  type_int npart[N_part_types];
  double   mass[N_part_types];
  double   time;
  double   redshift;
  type_int flag_sfr;
  type_int flag_feedback;
  type_int npartTotal[N_part_types];
  type_int flag_cooling;
  type_int num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  char     fill[256- N_part_types*4- N_part_types*8- 2*8- 2*4- N_part_types*4- 2*4 - 4*8];  /* fills to 256 Bytes */
};

extern void read_gadget(void);
#ifdef CHANGE_POSITION
  extern void change_positions(type_int n);
#endif

#endif
