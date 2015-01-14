
/*
 *  fstack.c
 *  flexible stack implementaion that takes care
 *  of all the allocation work
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 07/07/08 08:35:38 CEST
 *  
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "memory.h"
#include "fstack.h"

int
initfstack(fstack* stack, size_t size, int alloc) {
  
  stack->elems=realloc(NULL, size*alloc);
  if (stack->elems == NULL) return -1;
 
  stack->size = size;
  stack->alloc = alloc;
  stack->top = -1;

  return 0;
}

unsigned char fstackisempty(fstack *stack){
    return ((unsigned char)(stack->top < 0));
}

int
fstackpush(fstack *stack, void* elem) {
    char *ptr;

    if(stack->top >= stack->alloc-1) {
        stack->elems = realloc(stack->elems, stack->size*(stack->alloc+stack->inc)); 
        if (stack->elems == NULL) exit(-1);
        stack->alloc +=  stack->inc;        
    }

    stack->top++;
    ptr = (char*) stack->elems;
    memmove((ptr+(stack->top*stack->size)), elem, stack->size);
    
    return 0;
}

void*
fstackpop(fstack *stack) {
  void *elem;
  char *ptr;

  if(fstackisempty(stack)) return NULL;
  
  /*cleanup*/
  if(stack->top < stack->alloc-stack->inc) {
    ptr = realloc(stack->elems, stack->size*(stack->top+1)); 
    if(ptr == NULL) return NULL;
    stack->elems = ptr;
    stack->alloc = stack->top+1;
  }
  
  
  ptr = (char*) stack->elems;
  elem= (ptr+(stack->size*stack->top));
 
  stack->top--;
  return elem;
}

void*
fstacktop(fstack *stack) {
  void *elem;
  char *ptr;

  if(fstackisempty(stack)) return NULL;
  ptr = (char*) stack->elems;
  elem= (ptr+(stack->size*stack->top));

  return elem;
}

void
destructfstack(void *space, fstack *stack){
    free(stack->elems);
}

