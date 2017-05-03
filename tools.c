#include <stdio.h>
#include "tools.h"

void print_compiler_options(void)
{
  printf("------------ Compiler Options ----------\n");
#ifdef DEBUG
  printf("- DEBUG \n");
#endif
#ifdef PRECDOUBLE
  printf("- PRECDOUBLE \n");
#endif
#ifdef LONGIDS
  printf("- LONGIDS \n");
#endif
#ifdef IPER
  printf("- IPER = %1d\n",IPER);
#endif
#ifdef ncores
  printf("- NCORES = %1d\n",ncores);
#endif
#ifdef NHALO
  printf("- NHALO = %1d\n",NHALO);
#endif
#ifdef limpiamelo
  printf("- limpiamelo\n");
#endif
  printf("----------------------------------------\n");
}
