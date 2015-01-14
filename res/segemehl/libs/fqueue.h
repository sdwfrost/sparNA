#ifndef FQUEUE_H
#define FQUEUE_H

/*
 *
 *	fqueue.h
 *  declaration for flexible queues
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 07/07/08 09:30:29 CEST  
 *
 *
 */

typedef struct {

    size_t size;
    int alloc;
    int enqueueidx;
    int dequeueidx;
    int noofelems;
    void *elems;
} fqueue;

void wrapfqueue(fqueue *queue);
void* fdequeue (fqueue *queue);
int fenqueue(fqueue *queue, void *elem);
int resizefqueue(fqueue *queue);
int initfqueue(fqueue *queue, size_t size, int alloc);
unsigned char fqueueisempty (fqueue *queue);
void *fqueuefront(fqueue *queue);
void* fqueueget(fqueue *queue, unsigned int i);
int fqueuejump(fqueue *queue, void *elem);
#endif
