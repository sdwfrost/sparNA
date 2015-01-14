
/*
 *  nw.c
 *  needleman wunsch
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 11/15/2010 12:12:27 AM CET
 *  
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "basic-types.h"
#include "memory.h"
#include "mathematics.h"
#include "alignment.h"
#include "sw.h"


/*--------------------------------- nwalign ----------------------------------
 *      
 *  needleman-wunsch global similarity alignment
 *  returns a matrix of size (m+1)*(n+1) where m is length of given sequence a
 *  and n the length of sequence b, respectively. Function expects
 *  a function to calculate a substitution score
 *     
 */

  int*
nwmatrix (void *space, symtype *a, Uint m, symtype *b, Uint n, int indel,
    Sint (*sub)(symtype, symtype, void *), void *nfo)
{
  int i, j, cols, rows, size;
  int *L;

  rows = m+1;
  cols = n+1;

  size = rows*cols;
  L = ALLOCMEMORY(space, NULL, int, size);
  L = memset(L, 0, sizeof(int)*size);

  for(i=1; i < m+1; i++) {
    MATRIX2D(L, cols, i, 0) = i*indel;
    for(j=1; j < n+1; j++) {
      MATRIX2D(L, cols, 0, j) = j*indel;

      MATRIX2D(L, cols, i, j) = 
        MAX3(
            MATRIX2D(L, cols, (i-1), j) + indel ,    
            MATRIX2D(L, cols, i, (j-1)) + indel , 
            MATRIX2D(L, cols, (i-1), (j-1)) + sub(a[i-1], b[j-1], nfo)
            );
    }
  }

  return L;
}


/*------------------------------- nwtraceback --------------------------------
 *      
 *  traceback to find optimal global alignment path
 *   
 */

  void
nwtraceback (void *space, int *M,  
    symtype *a, Uint m, symtype *b, Uint n, int indel,
    Sint (*sub)(symtype, symtype, void *), void *nfo, Alignment *al)
{
  Uint i, j, ncol, cur; 

  ncol = n+1;
  i = m;
  j = n;

  al->uoff = 0;
  al->voff = 0;

  while(i > 0 && j > 0) {

    cur = MATRIX2D(M, ncol, i, j);
 //   fprintf(stderr, "enter while cur = %d, %c-%c, %d\n", cur, a[i-1], b[j-1],
   //   sub(a[i-1], b[j-1], nfo) );
    if (MATRIX2D(M, ncol, i-1, j) + indel == cur){
      insertEop(al, Insertion);
      i--;   

  //        fprintf(stderr, "insertion\n");
    } else {
      if (MATRIX2D(M, ncol, i, j-1) + indel == cur) {
        insertEop(al, Deletion);
        j--;

    //      fprintf(stderr, "deletion\n");
      } else {
        if (MATRIX2D(M, ncol, i-1, j-1)+sub(a[i-1], b[j-1], nfo) 
            == cur){
          insertEop(al, Replacement);
          i--; j--;
      //    fprintf(stderr, "replacement\n");
        }
        else {
          assert(cur == 0);
        //  fprintf(stderr, "asserting.\n");
          al->uoff = i;
          al->voff = j;

          revMeops(al);
          return;
        }
      }
    }
  }

  al->uoff = i;
  al->voff = j;
  revMeops(al);

  return;
}

/*--------------------------------- sgalign ----------------------------------
 *      
 *  semi-global similarity alignment
 *  returns a matrix of size (m+1)*(n+1) where m is length of given sequence a
 *  and n the length of sequence b, respectively. Function expects
 *  a function to calculate a substitution score
 *     
 */

  int*
sgmatrix (void *space, symtype *a, Uint m, symtype *b, Uint n, int indel,
    Sint (*sub)(symtype, symtype, void *), void *nfo)
{
  int i, j, cols, rows, size;
  int *L;

  rows = m+1;
  cols = n+1;

  size = rows*cols;
  L = ALLOCMEMORY(space, NULL, int, size);
  L = memset(L, 0, sizeof(int)*size);

  for(i=1; i < m+1; i++) {
    MATRIX2D(L, cols, i, 0) = i*indel;
    for(j=1; j < n+1; j++) {
      MATRIX2D(L, cols, 0, j) = 0;

      MATRIX2D(L, cols, i, j) = 
        MAX3(
            MATRIX2D(L, cols, (i-1), j) + indel ,    
            MATRIX2D(L, cols, i, (j-1)) + indel , 
            MATRIX2D(L, cols, (i-1), (j-1)) + sub(a[i-1], b[j-1], nfo)
            );
    }
  }

  return L;
}

/*------------------------------- sgtraceback --------------------------------
 *      
 *  traceback to find optimal semi global alignment path
 *   
 */

  void
sgtraceback (void *space, int *M,  
    symtype *a, Uint m, symtype *b, Uint n, int indel,
    Sint (*sub)(symtype, symtype, void *), void *nfo, Alignment *al)
{
  Uint i, j, ncol, cur; 
  Uint maxcol;

  ncol = n+1;
  i = m;
  maxcol = n;

  for(j=0; j < n; j++) {
    if(MATRIX2D(M, ncol, i, maxcol) < MATRIX2D(M, ncol, i, j))
      maxcol = j;
  }

  j = maxcol;

  al->uoff = 0;
  al->voff = 0;

  while(i > 0 && j > 0) {

    cur = MATRIX2D(M, ncol, i, j);
 //   fprintf(stderr, "enter while cur = %d, %c-%c, %d\n", cur, a[i-1], b[j-1],
   //   sub(a[i-1], b[j-1], nfo) );
    if (MATRIX2D(M, ncol, i-1, j) + indel == cur){
      insertEop(al, Insertion);
      i--;   

  //        fprintf(stderr, "insertion\n");
    } else {
      if (MATRIX2D(M, ncol, i, j-1) + indel == cur) {
        insertEop(al, Deletion);
        j--;

    //      fprintf(stderr, "deletion\n");
      } else {
        if (MATRIX2D(M, ncol, i-1, j-1)+sub(a[i-1], b[j-1], nfo) 
            == cur){
          insertEop(al, Replacement);
          i--; j--;
      //    fprintf(stderr, "replacement\n");
        }
        else {
          assert(cur == 0);
        //  fprintf(stderr, "asserting.\n");
          al->uoff = 0;
          al->voff = 0;

          revMeops(al);
          return;
        }
      }
    }
  }

  al->uoff = 0;
  al->voff = 0;
  revMeops(al);

  return;
}

