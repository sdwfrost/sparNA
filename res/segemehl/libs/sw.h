#ifndef DPALIGN_H
#define DPALIGN_H

/*
 *
 *	sw.h
 *  
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 06.02.2010 14:14:53 CET  
 *
 */

 


#include "basic-types.h"
#include "alignment.h"

typedef char symtype;

Uint edist(void *, symtype *, Uint, symtype *, Uint, Uint, Uint *, Uint);
int constscr (symtype, symtype, void *);
int constscr_Nmatch (symtype, symtype, void *);
int* swmatrix (void *, symtype*, Uint, symtype*, Uint, int,
                   Sint (*sub)(symtype, symtype, void *), void *);

int* splicedmatrix (void *, symtype*, symtype*, Uint, symtype*, Uint, symtype*, Uint,
		    Uint, Uint, int, Sint (*sub)(symtype, symtype, void *), void *);

int* swgapless (void *, symtype*, Uint, symtype*, Uint, 
                   Sint (*sub)(symtype, symtype, void *), void *);

int* swgaplesstraceback (void *, int *,  
                 symtype *, Uint, symtype *, Uint, 
                             Sint (*sub)(symtype, symtype, void *), void *, int*);

void 
swtraceback (void *space, int *M,  
    symtype *a, Uint m, symtype *b, Uint n, int indel,
    Sint (*sub)(symtype, symtype, void *), void *nfo, Alignment *al);

void 
splicedtraceback (void *space, int *M, symtype *a1, symtype *a2, Uint m,
		  symtype *b1, Uint n1, symtype *b2, Uint n2, Uint strand1, Uint strand2, int indel,
		  Sint (*sub)(symtype, symtype, void *), void *nfo, Alignment *al1, Alignment *al2);

  
void
swtracebacksplit (void *space, int *L,  
    symtype *a, symtype *b, Uint m, symtype *s1, symtype *s2, Uint n, int indel, 
    unsigned char rc,
    Sint (*sub)(symtype, symtype, void *), void *nfo, Alignment *al1, Alignment *al2, FILE *dev);


int*
swsplitalign (void *space, symtype *a, symtype *b, 
    Uint m, symtype *s1, symtype *s2, Uint n, int indel,
    unsigned char rc, Sint (*sub)(symtype, symtype, void *), void *nfo);


void
localsplicedtraceback (void *space, int *M, symtype *a1, symtype *a2, Uint m, symtype *b1, Uint n1,
    symtype *b2, Uint n2, Uint strand1, Uint strand2, int indel,
    Sint (*sub)(symtype, symtype, void *), void *nfo,
    Alignment *al1, Alignment *al2, int* lmv, int *lmr);


int*
localsplicedmatrix (void *space, symtype *a1, symtype *a2, Uint m,
    symtype *b1, Uint n1, symtype *b2, Uint n2, 
    Uint strand1, Uint strand2, int indel, 
    Sint (*sub)(symtype, symtype, void *), void *nfo, int **lv, int **lr);

void
localsplicedtraceback_test (void *space, int *M, symtype *a1, symtype *a2, Uint m, symtype *b1, Uint n1,
    symtype *b2, Uint n2, Uint strand1, Uint strand2, int indel,
    Sint (*sub)(symtype, symtype, void *), void *nfo,
    Alignment *al1, Alignment *al2, int* lmv, int *lmr);


int*
localsplicedmatrix_test (void *space, symtype *a1, symtype *a2, Uint m,
    symtype *b1, Uint n1, symtype *b2, Uint n2, 
    Uint strand1, Uint strand2, int indel,
    Sint (*sub)(symtype, symtype, void *), void *nfo, int **lv, int **lr);

int*
localmultisplicedmatrix (void *space, symtype *a1, symtype *a2, Uint m,
    symtype **b, Uint *n, Uint *strand, Uint noofseqs, int indel, int trans,
    Sint (*sub)(symtype, symtype, void *), void *nfo, int ***lv, int ***lr, int ***lc);

int***
localmultisplicedmatrixopt (void *space, symtype *a1, symtype *a2, Uint qrylen, Uint *m,
    symtype **b, Uint *n, Uint *strand, Uint *qstart, Uint *qend, Uint *tstart, Uint *tend, Uint noofseqs, int indel, int trans,
    Sint (*sub)(symtype, symtype, void *), void *nfo, int ***lv, int ***lr, int ***lc, PairUint **scr, char ****, PairUint *diag);

void
localmultisplicedtraceback (void *space, int *M, symtype *a1, symtype *a2, Uint m, 
    symtype **b, Uint* n, Uint *strand, Uint noofseqs, int indel, int trans,
    Sint (*sub)(symtype, symtype, void *), void *nfo,
    Alignment **al, int **lmv, int **lmr, int **lmc);

char***  
localmultisplicedtracebackopt (void *space, int ***M, symtype *a1, symtype *a2, Uint qrylen, Uint *m, 
    symtype **b, Uint* n, Uint *strand, Uint *qstart, Uint *qend, Uint *tstart, Uint *tend, Uint noofseqs, int indel, int trans,
    Sint (*sub)(symtype, symtype, void *), void *nfo,
    Alignment **al, int **lmv, int **lmr, int **lmc, PairUint *);
#endif

