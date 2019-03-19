#ifndef IDENTIFICADOR_H
#define IDENTIFICADOR_H

struct iden_st{
	double r0;        /*Linking Length para el FoF*/
	my_int nobj;
	int    iper;
  my_int ngrupos;
} iden;

struct listvec{
  my_int indx;
  struct listvec * next;
};

typedef struct listvec item;

void busv(my_int *curr, int *test);
//void busv(int *ic, unsigned int *nvec, bool *test, int *ll);
void identification(void);

#endif
