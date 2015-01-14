#ifndef FSTACK_H
#define FSTACK_H

/*
 *
 *	fstack.h
 *  flexible stack declarations
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 07/07/08 08:36:50 CEST  
 *
 */

#include <stdlib.h>
#include <string.h>

#define FSTACKINC 1000

typedef struct {
  size_t size;
  int    alloc;
  int    top;
  void   *elems;
} fstack;


int initfstack(fstack *stack, size_t size, int alloc);
unsigned char fstackisempty(fstack *stack);
int fstackpush(fstack *stack, void* elem);
void* fstackpop(fstack *stack);
void* fstacktop(fstack *stack);
void destructstack(fstack *stack);

#endif

