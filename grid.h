#ifndef GRID_H
#define GRID_H

#ifndef NGRIDMAX
#define NGRIDMAX 512
#endif

struct gridst{
	my_int ngrid;
	my_int nobj;
	my_int *llirst;
	my_int *ll;
} grid;


void grid_init(void);
void grid_build(void);
void grid_free(void);

#endif
