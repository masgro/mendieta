#ifndef TOOLS_H
#define TOOLS_H

void print_compiler_options(void);

#define my_fread(ptr,size,nmemb,stream) assert(fread(ptr,size,nmemb,stream) == nmemb);

#endif
