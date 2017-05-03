#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <stdarg.h>
#include "linkedlist.h"

void initll(int *head, int *tail, int ng, int *ll, int np){
	int j;

	for(j=0;j<ng;j++){
		head[j]=-1;
		tail[j]=-1;
	}

	for(j=0;j<np;j++){
		ll[j]=-1;
	}
}

void agrega(int *head, int *tail, int i, int *ll){
	if(*head==-1){
		*head=i;
		*tail=i;
	}else{
		ll[*tail]=i;
		*tail=i;
	}
	ll[*tail]=-1;
	return;
}

void borra(int *head, int *tail, int i, int *ll){
	int temp,prev;

	temp=i;
	prev=*head;
	if(temp==prev){
		*head=ll[*head];
		if(*tail==temp){
			*tail=ll[*tail];
		}
	}else{
		while(ll[prev]!=temp){
			prev=ll[prev];
		}
		ll[prev]=ll[temp];
		if(*tail==temp)
			*tail=prev;
	}
	ll[temp]=-1;
}

void checkll(int ll,...){
	char flag[200];
	va_list pa;
	char * t;

	va_start (pa, ll);
	if(t = va_arg(pa, char *)){
		strcpy(flag, t);
  }
  va_end (pa);
	if(ll==-1){
		printf("final de lista encontrado. flag: %s\n",flag);
		exit(0);
	}
}

void borra_ll(int *head, int *tail, int *ll){
	int l;
	int temp;

	if(*head==-1)return;

	l=*head;
	while(l!=-1){
		temp=ll[l];
		ll[l]=-1;
		l=temp;
	}
	*head=-1;
	*tail=-1;
}
