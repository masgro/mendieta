#ifndef IDENTIFICADOR_H
#define IDENTIFICADOR_H

struct iden_st
{
	double       r0;        /*Linking Length para el FoF*/
	unsigned int nobj;
	int          iper;
  unsigned int ngrupos;
} iden;

struct listvec
{
  int indx;
  struct listvec * next;
};

typedef struct listvec item;

void busv(item **curr, unsigned int *nvec, int *test);
//void busv(int *ic, unsigned int *nvec, bool *test, int *ll);
void identification(void);

#endif
