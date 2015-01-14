
/*
 *
 *	nw.h
 *  
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 11/15/2010 12:52:56 AM CET  
 *
 */
   
#include "basic-types.h"
#include "alignment.h"
#include "sw.h"

  int*
nwmatrix (void *space, symtype *a, Uint m, symtype *b, Uint n, int indel,
    Sint (*sub)(symtype, symtype, void *), void *nfo);

  
  void
nwtraceback (void *space, int *M,  
    symtype *a, Uint m, symtype *b, Uint n, int indel,
    Sint (*sub)(symtype, symtype, void *), void *nfo, Alignment *al);

  int*
sgmatrix (void *space, symtype *a, Uint m, symtype *b, Uint n, int indel,
    Sint (*sub)(symtype, symtype, void *), void *nfo);


  void
sgtraceback (void *space, int *M,  
    symtype *a, Uint m, symtype *b, Uint n, int indel,
    Sint (*sub)(symtype, symtype, void *), void *nfo, Alignment *al);

