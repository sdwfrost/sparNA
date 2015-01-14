#ifndef ALURUSORT_H
#define ALURUSORT_H

/*
 *
 *	aluruSort.h
 *  declarations for aluruSort.c
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 07/11/2007 07:06:04 PM CEST  
 *
 *  SVN
 *  Revision of last commit: $Rev: 31 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-05-16 14:11:09 +0200 (Fri, 16 May 2008) $
 *
 *  Id: $Id: aluruSort.h 31 2008-05-16 12:11:09Z steve $
 *  Url: $URL: file:///homes/bierdepot/steve/svn/segemehl/trunk/libs/aluruSort.h $
 */

#include "mathematics.h"
#include "bitArray.h"

#define BUCKETRET(X, N, I, B)   for(B=0; B < (N); B++) { \
                                    if(X[B].id == (I)) break; \
                                }
                              
#define BUCKETINIT(X,I)       (X).id = (I);\
                              (X).elems=NULL; \
                              (X).noofelems=0;\
                              (X).allocsize=0
#define BUCKETSWAP(X,A,B)     {\
                                Uint resc = (X)->elems[(A)]; \
                                (X)->elems[(A)] = (X)->elems[(B)];\
                                (X)->elems[(B)] = resc;\
                              }     
#define BUCKETINC             1000
#define BUCKETADD(S,X,E)     {\
                                if((X).allocsize <= (X).noofelems) {\
                                    (X).elems = ALLOCMEMORY((S), (X).elems, Uint, ((X).allocsize+=BUCKETINC)); \
                                }\
                                (X).elems[(X).noofelems++] = (E);\
                              }
#define BUCKETFRONT(X)       (X)->front++
#define SUFA(X, N, I)        {\
                                Uint ITER;\
                                Uint CUMSUM=0;\
                                for(ITER=0; ITER < N; ITER++) {\
                                  if(CUMSUM+X[ITER].noofelems > I) {\
                                      break;\
                                    }\
                                    CUMSUM+=X[ITER].noofelems;\
                                }\
                             } X[ITER].elems[I-CUMSUM]
typedef struct {
    Uint id;
    Uint *elems;
    Uint noofelems;
    Uint front;
    Uint allocsize;

} Alurubucket;

Alurubucket* getAluruBuckets(void *, char *, Uint, Uint *, Uint **);
Uint* Qdist(void *, bitarray, Uint, unsigned char);
bitarray classify(void *, char* s, Uint, Uint*, Uint*);
Uint* getAluruArray(void *space, char *s, Uint len, char delim);
vector_t** getQdistList(void *space, Alurubucket*, Uint, Uint *, Uint);
void showAluruBuckets(Alurubucket *bckts, Uint*, Uint n);
void sortAluruSubstrings (void *space, 
    vector_t** qdlist, 
                     Uint nooflists, 
                     Alurubucket *bckts, 
                     Uint noofbckts,
                     Uint *R,
                     char *cl,
                     Uint len,
                     char Q);

Uint* alurusortint(void *space, Uint *s, Uint *l);
Uint* alurusort(void *space, char *s, Uint *l);
extern void getinterval(Uint *s, Uint len, Uint *min, Uint *max);
#endif
