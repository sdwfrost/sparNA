
/*
 *  fqueue.c
 *  implementation for a flexible circular queue that
 *  takes care of memory management
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 07/07/08 09:31:37 CEST
 *  
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "memory.h"
#include "fqueue.h"

unsigned char 
fqueueisempty (fqueue *queue) {
    return((unsigned char) (queue->noofelems==0));
}

int
initfqueue(fqueue *queue, size_t size, int alloc) {
    
  queue->elems=realloc(NULL, size*alloc);
  if (queue->elems == NULL) return -1;
  
  queue->size = size;
  queue->alloc = alloc;
  queue->noofelems = 0;
  queue->enqueueidx = 0;
  queue->dequeueidx = 0;

  return 0;
}


int 
resizefqueue(fqueue *queue) {
  char *q, *src, *dst, *ptr;

  queue->elems = realloc(queue->elems, queue->size*queue->alloc*2);
  if(queue->dequeueidx >= queue->enqueueidx) {
    q=(char*) queue->elems;
    src = q+(queue->dequeueidx*queue->size);
    dst = q+((queue->dequeueidx+queue->alloc)*queue->size);
    ptr = memmove(dst, src, (queue->alloc-queue->dequeueidx)*queue->size);
    if (ptr != dst) return -1;
    queue->dequeueidx += queue->alloc;
  }

  queue->alloc = queue->alloc*2;
  return 0;
}


int
fenqueue(fqueue *queue, void *elem) {
  int r;
  char *ptr;

  if(queue->noofelems == queue->alloc) {
    r = resizefqueue(queue);
    if (r != 0) return -1;  
  }

  ptr = (char*)queue->elems;
  ptr = ptr+(queue->enqueueidx*queue->size);
  ptr = memmove(ptr, elem, queue->size);
  
  if (ptr == NULL) return -2;    
  queue->noofelems++;
  if(queue->enqueueidx == queue->alloc-1) {
    queue->enqueueidx = 0;
  } else {
    queue->enqueueidx++;
  }

  return 0;
}


int
fqueuejump(fqueue *queue, void *elem) {
  int r;
  char *ptr;

  if(queue->noofelems == queue->alloc) {
    r = resizefqueue(queue);
    if (r != 0) return -1;  
  }

  if(queue->dequeueidx == 0) { 
    queue->dequeueidx = queue->alloc-1;
  } else {
    queue->dequeueidx--;
  }

  assert(queue->dequeueidx != queue->enqueueidx);

  ptr = (char*)queue->elems;
  ptr = ptr+(queue->dequeueidx*queue->size);    
  ptr = memmove(ptr, elem, queue->size);
  
  if (ptr == NULL) return -2;    
  queue->noofelems++;
 
  return 0;
}



void*
fdequeue (fqueue *queue) {
    char *elem;
    
    if(fqueueisempty(queue)) {
        return NULL;
    } 

    elem = (char*) queue->elems;
    elem = elem+(queue->dequeueidx*queue->size);
    queue->noofelems--;

    if(queue->dequeueidx == queue->alloc-1) {
        queue->dequeueidx = 0;
    } else {
        queue->dequeueidx++;
    }
    return elem;
}


void*
fqueuefront(fqueue *queue) {
  char *elem;
    
    if(fqueueisempty(queue)) {
        return NULL;
    } 

    elem = (char*) queue->elems;
    elem = elem+(queue->dequeueidx*queue->size);

    return elem;
}

void*
fqueueget(fqueue *queue, unsigned int k) {
 char *elem;
 unsigned int i,j;
    
    if(fqueueisempty(queue) || k >= queue->noofelems) {
        return NULL;
    } 

    i = queue->dequeueidx; 
    for(j=0; j < k; j++) {

      if(i == queue->alloc-1) {
        i = 0;
      } else {
        i++;
      }
    }

    elem = (char*) queue->elems;
    elem = elem+(i*queue->size);

    return elem;
}

void
wrapfqueue(fqueue *queue) {
    free(queue->elems);
}





