
/*
 *  sw.c
 *  
 *  local alignments
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 06.02.2010 13:47:30 CET
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
#include "kdchain.h"

/*---------------------------------- edist -----------------------------------
 *    
 * evaluation of the edit distance in O(n) space
 * the function accepts two sequences of symtype and
 * an m-dimensional substitution matrix 
 * 
 */

Uint edist(void *space, symtype *sa, Uint lena, symtype *sb, Uint lenb, 
    Uint indel, Uint *sub, Uint m) {

  Uint minlen, maxlen;
  Uint i, j, r1=0, r2=0;
  Uint pen;
  Uint* col;
  symtype *min;
  symtype *max;

  minlen = (MIN(lena, lenb)) + 1;
  maxlen = (MAX(lena, lenb)) + 1;
  col = ALLOCMEMORY(space, NULL, Uint, minlen);

  if(minlen-1 == lena) {
    min = sa;
    max = sb;
  } else {
    min = sb;
    max = sa;
  }

  for(i=0; i < maxlen; i++) {
    for(j=0; j < minlen; j++) {
      if (i==0) {
        col[j]=j; 
      } else {
        if (j==0) {
          r1 = col[j]++; 
        } else {
          r2 = col[j];
          if (sub == NULL) {
            pen = (min[j-1]==max[i-1]) ? 0 : 1;
          } else {
            pen = (min[j-1]==max[i-1]) ? 0 : 
              MATRIX2D(sub, m, min[j-1], max[i-1]);
          }
          col[j] = MIN ((col[j-1]+indel),  
              (MIN ((r1 + pen), (col[j]+indel))));
          r1 = r2;    
        }
      }
    }
  }

  r1 = col[minlen-1];
  FREEMEMORY(space, col);
  return r1;
}

/*------------------------------- constscr_Nmatch ------------------------------
 *      
 *  a function that assigns constant scores for matches and mismatches
 *  given in info[0] and info[1], respectively.
 *    
 */

  int
constscr_Nmatch (symtype a, symtype b, void *info)
{   
  int* scores;

  scores = (int*) info;
  if(a == b || a == 'N' || b == 'N') 
    return scores[0];

  return scores[1];
}


/*--------------------------------- constscr ----------------------------------
 *      
 *  a function that assigns constant scores for matches and mismatches
 *  given in info[0] and info[1], respectively.
 *    
 */

  int
constscr (symtype a, symtype b, void *info)
{   
  int* scores;

  scores = (int*) info;
  if(a == b) return scores[0];

  return scores[1];
}


/*--------------------------------- swgapless ----------------------------------
 *      
 *  smith-waterman local similarity alignment w/o gaps
 *  returns a matrix of size (m+1)*(n+1) where m is length of given sequence a
 *  and n the length of sequence b, respectively. Function expects
 *  a function to calculate a substitution score
 *     
 */

  int*
swgapless (void *space, symtype *a, Uint m, symtype *b, Uint n,
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
    for(j=1; j < n+1; j++) {

      MATRIX2D(L, cols, i, j) = MAX(0,
          MATRIX2D(L, cols, (i-1), (j-1)) + sub(a[i-1], b[j-1], nfo));
    }
  }

  return L;
}



/*--------------------------------- swalign ----------------------------------
 *      
 *  smith-waterman local similarity alignment
 *  returns a matrix of size (m+1)*(n+1) where m is length of given sequence a
 *  and n the length of sequence b, respectively. Function expects
 *  a function to calculate a substitution score
 *     
 */

  int*
swmatrix (void *space, symtype *a, Uint m, symtype *b, Uint n, int indel,
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
    for(j=1; j < n+1; j++) {

      MATRIX2D(L, cols, i, j) = 
        MAX4(0,
            MATRIX2D(L, cols, (i-1), j) + indel ,    
            MATRIX2D(L, cols, i, (j-1)) + indel , 
            MATRIX2D(L, cols, (i-1), (j-1)) + sub(a[i-1], b[j-1], nfo)
            );
    }
  }

  return L;
}


/*------------------------------- swtraceback --------------------------------
 *      
 *  traceback to find optimal local alignment path
 *   
 */

  void
swtraceback (void *space, int *M,  
    symtype *a, Uint m, symtype *b, Uint n, int indel,
    Sint (*sub)(symtype, symtype, void *), void *nfo, Alignment *al)
{
  Uint i, j, ncol, cur, start; 

  ncol = (n+1);
  start = arraymax(M, (m+1)*ncol);
  i = start / ncol;
  j = start % ncol;

  al->uoff = 0;
  al->voff = 0;

  while(i > 0 && j > 0) {

    cur = MATRIX2D(M, ncol, i, j);
    if (MATRIX2D(M, ncol, i-1, j) + indel == cur){
      insertEop(al, Insertion);
      i--;   
    } else {
      if (MATRIX2D(M, ncol, i, j-1) + indel == cur) {
        insertEop(al, Deletion);
        j--;
      } else {
        if (MATRIX2D(M, ncol, i-1, j-1)+sub(a[i-1], b[j-1], nfo) 
            == cur){
          insertEop(al, Replacement);
          i--; j--;
        }
        else {
          assert(cur == 0);

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

/*------------------------------- splicescore --------------------------------
 *    
 * @brief consensus splice site score
 * @author Steve Hoffmann 
 *   
 */
 
int
splicescore (char a, char b, Uint strand, Uint type)
{

  if((type == 0 && strand == 0) ||
     (type == 1 && strand == 1)) {
    if((a == 'G' && b == 'T') ||
       (a == 'C' && b == 'T')) return 1;
  } else {
    if((a == 'A' && b == 'G') ||
       (a == 'A' && b == 'C')) return 1;
  }
  return 0;
}


#define LOCALMULTISPLICECALCOFF

/*------------------------- localmulitsplicedmatrixopt --------------------------
 *    
 * @brief aligning a read a1 and its reverse complement a2 of length m to 
 * noofseqs loci
 * @author Steve Hoffmann 
 *   
 */

int***
localmultisplicedmatrixopt (void *space, symtype *a1, symtype *a2, Uint qrylen, Uint *m,
    symtype **b, Uint *n, Uint *strand, Uint *qstart, Uint *qend, Uint *tstart, Uint *tend, Uint noofseqs, int indel, int trans,
    Sint (*sub)(symtype, symtype, void *), void *nfo, int ***lv, int ***lr, int ***lc, PairUint **bestscr, char ****KBAND, PairUint *diag){

  int i, j, k, q, cols=0, rows=0, abs, relq, tstartq, tendq, lk, rk;
  unsigned int ovlrange;
  PairUint *scr;
  Uint off,start,margin=50;
  int ***L, **lmr=NULL, **lmv=NULL, **lmc=NULL, tmp;
  symtype cura, curb, r1, r2, l1, l2;

#ifdef DEBUGMULTISPLICEOPT      
  int maxk, maxi, maxj, lastpick, maxpick;
#endif

#ifdef DEBUGKBAND
  char ***K;
  K = ALLOCMEMORY(space, NULL, char**, noofseqs);
#endif

  lmr = ALLOCMEMORY(space, NULL, int*, noofseqs);
  lmv = ALLOCMEMORY(space, NULL, int*, noofseqs);
  lmc = ALLOCMEMORY(space, NULL, int*, noofseqs);
  L = ALLOCMEMORY(space, NULL, int**, noofseqs);
  scr = ALLOCMEMORY(space, NULL, PairUint, noofseqs);
  memset(scr, 0, sizeof(PairUint)*noofseqs);


  for(k=0; k < noofseqs; k++) {
#ifdef DEBUGKBAND
    K[k] = ALLOCMEMORY(space, NULL, char*, m[k]+1);
#endif
    cols += n[k] + 1;
    rows += m[k] + 1;
    lmr[k] = ALLOCMEMORY(space, NULL, int, m[k]+1);
    lmv[k] = ALLOCMEMORY(space, NULL, int, m[k]+1);
    lmc[k] = ALLOCMEMORY(space, NULL, int, m[k]+1);
    L[k] = ALLOCMEMORY(space, NULL, int*, m[k]+1);
    for(i=0; i < m[k]+1; i++) {
#ifdef DEBUGKBAND
      K[k][i] = ALLOCMEMORY(space, NULL, char, n[k]+1);
      memset(K[k][i], 0, sizeof(n[k]+1));
#endif
      L[k][i]=ALLOCMEMORY(space, NULL, int, n[k]+1);
      memset(L[k][i], 0, sizeof(int)*(n[k]+1));
    }
    memset(lmv[k], 0, sizeof(int)*m[k]+1);
    memset(lmr[k], 0, sizeof(int)*m[k]+1);
    memset(lmc[k], 0, sizeof(int)*m[k]+1);
  }

  for(k=0; k < noofseqs; k++) {
#ifdef DEBUGMULTISPLICEOPT      
    maxk = 0;
    maxi = 0;
    maxj = 0;
    lastpick = 0;
    maxpick = 0;
#endif

    ovlrange =0;

    for(q=0; q < k; q++) {
#ifdef LOCALMULTISPLICECALCOFF
      if(strand[q] == 0) {
        tstartq = qstart[q];
        tendq = qend[q];
      } else {
        tstartq = qrylen - (qstart[q] + m[q]);
        tendq = tstartq + m[q] - 1;
      }
      assert(tstartq == tstart[q] && tendq == tend[q]);
#else
      tstartq = tstart[q];
      tendq = tend[q];
#endif

      if(tstartq <= tstart[k] && tstart[k] <= tendq && ovlrange < tendq) { 
        ovlrange = tendq;
      }
    }

    for (i=1; i < m[k]+1; i++) {
    
      lmv[k][i] = lmv[k][i-1];    
      lmr[k][i] = lmr[k][i-1];
      lmc[k][i] = lmc[k][i-1];

#ifdef LOCALMULTISPLICECALCOFF
      if(strand[k] == 0)
        abs = qstart[k] + i;
      else
        abs = qrylen - qstart[k] - m[k] + i;

      assert(tstart[k] +i == abs);
#else
      abs = tstart[k] + i;
#endif

      start = (diag[k].b > diag[k].a) ? diag[k].b - diag[k].a : 0;
      lk = (start + i > margin) ? start + i - margin : 1;
      rk = (start + i + margin < n[k]+1) ? start + i + margin : n[k]+1;

#ifdef NOKBAND
      for (j=1; j < n[k]+1; j++){
#else
      for (j=lk; j < rk; j++){
#endif

#ifdef DEBUGKBAND
        if(j >= lk && j <= rk) K[k][i][j] = 1;
#endif
        l1 = 0;
        l2 = 0;
        r1 = 0;
        r2 = 0;

        if (strand[k] == 0){
          off = qstart[k]; //off supports 0!
          cura = a1[off+i-1]; 
          curb = b[k][j-1];
          if(j+1 < n[k]) {
            r1 = b[k][j];
            r2 = b[k][j+1];
          } 
          if(j > 2) {
            l1 = b[k][j-3];
            l2 = b[k][j-2];
          }
        } else {
          off = qend[k]+1; //off should b m -> qend + 1
          cura = a2[off-i];
          curb = b[k][n[k]-j];
          if(j > 2) {
            l1 = b[k][n[k]-j+1];
            l2 = b[k][n[k]-j+2];
          }
          if(n[k]-j > 1) {
            r1 = b[k][n[k]-j-2];
            r2 = b[k][n[k]-j-1];
          }
        }	

        L[k][i][j] = MAX4(0, L[k][i-1][j]+indel, 
            L[k][i][j-1]+indel, 
            L[k][i-1][j-1]+sub(cura, curb, nfo));

       if(abs <= ovlrange) { 
        for(q=0; q < k; q++) {
#ifdef LOCALMULTISPLICECALCOFF
          if(strand[q] == 0) {
            tstartq = qstart[q];
            tendq = qend[q];
          } else {
            tstartq = qrylen - (qstart[q] + m[q]);
            tendq = tstartq + m[q] - 1;
          }
          assert(tstartq == tstart[q] && tendq == tend[q]);
#else
          tstartq = tstart[q];
          tendq = tend[q];
#endif


          if(tstartq < abs && abs < tendq) { 

            relq = abs - tstartq - 1;
            assert(relq < m[q]);

            tmp = MAX(L[k][i][j],
                lmv[q][relq] + sub(cura, curb, nfo) + trans + splicescore(l1, l2, strand[k], 1)) ;

#ifdef DEBUGMULTISPLICEOPT      
            if(tmp > L[k][i][j] && tmp >= maxpick) {
              maxpick = tmp;
              lastpick = q;
            }
#endif
            L[k][i][j] = tmp;
          }
        }
       }

       if(L[k][i][j] > L[k][scr[k].a][scr[k].b]) {
         scr[k].a = i;
         scr[k].b = j;
       }

        //if strand -  acc1, acc2
        //if strand +  don1, don2
        if (L[k][i][j] 
            + splicescore(r1, r2, strand[k], 0) 
            > lmv[k][i-1]) {
          lmv[k][i] = L[k][i][j] + splicescore(r1, r2, strand[k], 0);
          lmr[k][i] = i;
          lmc[k][i] = j;
        }
#ifdef DEBUGMULTISPLICEOPT      
        if(maxk < L[k][i][j] ) {
          maxk = L[k][i][j];
          maxi = abs;
          maxj = j;
        }
#endif
     
    
      }
    }
#ifdef DEBUGMULTISPLICEOPT      
    fprintf(stdout, "maximum calculated score on matrix  L[%d][%d][%d]=%d, maxtransistion from %d:%d\n", k, maxi, maxj, maxk, lastpick, maxpick);
#endif
  }

  *lv = lmv;
  *lr = lmr; 
  *lc = lmc;
  *bestscr = scr;

#ifdef DEBUGKBAND
  *KBAND=K;
#endif

  return L;
}


/*----------------------- localmultisplicedtracebackopt -------------------------
 *    
 * @brief tracing localmultisplicedmatrix back producing noofseqs split aligns
 * @author Steve Hoffmann 
 *   
 */

char***
localmultisplicedtracebackopt (void *space, int ***M, symtype *a1, symtype *a2, Uint qrylen, Uint *m, 
    symtype **b, Uint* n, Uint *strand, Uint *qstart, Uint *qend, Uint *tstart, Uint *tend, Uint noofseqs, int indel, int trans,
    Sint (*sub)(symtype, symtype, void *), void *nfo,
    Alignment **al, int **lmv, int **lmr, int **lmc, PairUint *scr){

  Uint u=1, v=1, cur, cols=0, rows=0, off=0, q, k, xmax=0, maxk, maxu, maxv, abs, relq, tstartq, tendq;
  int x, y, maxval;
  char breakiter=0;
  symtype cura, curb, l1, l2;//, r1, r2;
  char ***backtracetable = NULL;
#if DEBUGMULTISPLICEOPT
  Uint i,j;
#endif

#if DEBUGKBAND  
  backtracetable = ALLOCMEMORY(space, NULL, char**, noofseqs);
#endif

  maxval = 0;
  maxu = 0;
  maxv = 0;
  maxk = 0;

  for(k=0; k < noofseqs; k++) {

#if DEBUGKBAND
    backtracetable[k] = ALLOCMEMORY(space, NULL, char*, m[k]+1);
    for(i=0; i < m[k]+1; i++) {
      backtracetable[k][i] = ALLOCMEMORY(space, NULL, char, n[k]+1);
      memset(backtracetable[k][i], 0, sizeof(n[k]+1));
    }
#endif
    al[k]->uoff = 0;
    al[k]->voff = 0;
    cols += n[k] + 1;
    rows += m[k] + 1;
    
    if(maxval < M[k][scr[k].a][scr[k].b]) {
      maxk = k;
      maxu = scr[k].a;
      maxv = scr[k].b;
      maxval =  M[k][scr[k].a][scr[k].b];
    }


#if DEBUGMULTISPLICEOPT
    for(i=0; i < m[k]+1; i++) {
      for(j=0; j < n[k]+1; j++) {
        if(maxval < M[k][i][j]) {
            maxk = k;
            maxu = i;
            maxv = j;
            maxval = M[k][i][j];
        }
      }
    }
#endif    
  }


  k=maxk;  
  u = maxu;
  v = maxv;

  //fprintf(stdout, "\nstarting on seq: %d,%d,%d -> %d \n", k, u, v, M[k][u][v]);

  while(u > 0 && v > 0){
    cur = M[k][u][v];

#ifdef DEBUGKBAND
    backtracetable[k][u][v] = 1;
    fprintf(stdout, "continue with matrix[%d,%d,%d]:%d strand:%d [u'=%d]\n", k, u, v, M[k][u][v], strand[k], qstart[k]+u);
#endif

 
    if (M[k][u-1][v] + indel == cur){
      insertEop(al[k], Insertion);
      assert(u > 0);
      u--;
    } else {
      if (M[k][u][v-1] + indel == cur){
        insertEop(al[k], Deletion);
        assert(v > 0);
        v--;
      } else {
        if (M[k][u][v]) {
          l1 = 0;
          l2 = 0;

          if (strand[k] == 0){
            off = qstart[k];
            assert(off+u-1 <= qend[k]);
            cura = a1[off+u-1];
            curb = b[k][v-1];
            
             if(v > 2) {
              l1 = b[k][v-3];
              l2 = b[k][v-2];
            }
          } else {
            off = qend[k]+1;
            assert(off >= u);
            cura = a2[off-u];
            curb = b[k][n[k]-v];
            if(v > 2) {
              l1 = b[k][n[k]-v+1];
              l2 = b[k][n[k]-v+2];
            }
          }

//         fprintf(stdout, "sub[%c,%c]=%d\n", cura, curb, sub(cura,curb,nfo));
          if (M[k][u-1][v-1] + sub(cura, curb, nfo) == cur){
            insertEop(al[k], Replacement);
            assert(u > 0 && v > 0);
            u--; v--;
          } else {

            if(k==0) {
              if (strand[k] == 0){
                //TODO: do we need to reset offsets?
                off = qstart[k];
                al[k]->uoff = off+u;
                al[k]->voff = v;
                revMeops(al[k]);
              } else {  
                off = qend[k]+1; //off equals m in fulllocalmultsplice -> qend[k]+1
                assert(off >= u+getUalignlen(al[k]));  
                al[k]->uoff = off-u-getUalignlen(al[k]);

                assert(al[k]->uoff == qstart[k]+(m[k]-u-getUalignlen(al[k]))); //-1
                assert(n[k] >= v+getValignlen(al[k]));              
                al[k]->voff = n[k]-v-getValignlen(al[k]);
              }

#ifdef DEBUGMULTISPLICEOPT      
              fprintf(stdout, "exit 1 at k=%d\n", k);
#endif
              breakiter=1;
              break;
            }

            x = -1;
            y = -1;
#ifdef LOCALMULTISPLICECALCOFF      
            if(strand[k] == 0)
              abs = qstart[k] + u;
            else
              abs = qrylen - qstart[k] - m[k] + u; 
#else
          abs = tstart[k] + u;
#endif
            for(q=k; q > 0; q--) {   

#ifdef LOCALMULTISPLICECALCOFF
          if(strand[q-1] == 0) {
            tstartq = qstart[q-1];
            tendq = qend[q-1];
          } else {
            tstartq = qrylen - (qstart[q-1] + m[q-1]);
            tendq = tstartq + m[q-1]-1;
          }
 
          assert(tstartq == tstart[q-1] && tend[q-1] == tendq);
#else
          tstartq = tstart[q-1];
          tendq = tend[q-1];
#endif


          if(tstartq < abs && abs < tendq) { 

            relq = abs - tstartq -1;
            assert(relq < m[q-1]);
  
                if (lmv[q-1][relq] + splicescore(l1, l2, strand[k], 1) +
                    sub(cura, curb, nfo) + trans == cur) {
                  x=q-1;
                  y=relq;
                  xmax = lmc[q-1][relq];
                }
              }
            }            
//            fprintf(stderr, "matrix:U[%d,%d+%d=%d]:%d, k:%d, q:%d, x:%d -> [%d,%d]\n", u, cum[k], v, cum[k]+v, cur, k, q, x, lmr[x][u-1], xmax);
            assert(x > -1);
            insertEop(al[k], Replacement);
            assert(u > 0 && v > 0);
            u--; v--;

            if (strand[k] == 0){ 
              off = qstart[k];
              al[k]->uoff = off+u;
              al[k]->voff = v;
              revMeops(al[k]);
            } else {
              off = qend[k]+1; //off=m in fulllocalmultisplice -> qend[k]+1
              assert(off >= u+getUalignlen(al[k]));
              assert(n[k] >= v+getValignlen(al[k]));
              al[k]->uoff= off-u-getUalignlen(al[k]);
              assert(al[k]->uoff == qstart[k]+(m[k]-u-getUalignlen(al[k]))); //-1
              al[k]->voff = n[k]-v-getValignlen(al[k]);
            }
#ifdef DEBUGMULTISPLICEOPT          
            fprintf(stdout, "transition from k:%d - %d:x : m:%d, u:%d, v:%d, ulen:%d, vlen:%d -> u': %d v':%d \n", 
                k, x, m[k], u, v, getUalignlen(al[k]), getValignlen(al[k]), lmr[x][y], xmax);
            
            showAlign(al[k], stdout);
            fprintf(stdout, "\n");
#endif
            k = x;	    
            u = lmr[x][y];
            v = xmax;

            assert(u >= 0 && v >= 0); 
          }
        } else {
          //fprintf(stderr, "reversing end\n"); was k>0
          if (strand[k] == 0){
            off = qstart[k];
            al[k]->uoff = off+u;
            al[k]->voff = v;
            revMeops(al[k]);
          } else { 
            off = qend[k]+1; // off=m;
            assert(off >= u+getUalignlen(al[k]));
            al[k]->uoff = off-u-getUalignlen(al[k]);
            assert(al[k]->uoff == qstart[k]+(m[k]-u-getUalignlen(al[k]))); //-1
            assert(n[k] >= v+getValignlen(al[k]));
            al[k]->voff = n[k]-v-getValignlen(al[k]);
          }
#ifdef DEBUGMULTISPLICEOPT 
          fprintf(stdout, "exit 2 at k=%d\n", k);
#endif
          breakiter=1;
          break;
        }
      }
    }
  }    


//  if(k==0){ 
  if (!breakiter){
    if (strand[k] == 0){
      off=qstart[k];
      al[k]->uoff = off+u;
      al[k]->voff = v;
      revMeops(al[k]);
    } else {
  /*    fprintf(stdout, "finalize alignment at u:%d, v:%d, ulen:%d, vlen:%d\n", u, v, getUalignlen(al[k]), getValignlen(al[k]));*/
      off=qend[k]+1; //off = m;
      assert(off >= u+getUalignlen(al[k]));
      al[k]->uoff = off-u-getUalignlen(al[k]);
      assert(al[k]->uoff == qstart[k]+(m[k]-u-getUalignlen(al[k]))); //-1
      assert(n[k] >= v+getValignlen(al[k]));
      al[k]->voff = n[k]-v-getValignlen(al[k]);
    }
  }

#ifdef DEBUGMULTISPLICEOPT 
   fprintf(stdout, "starting on seq: %d (%d), u:%d, v:%d, uoff:%d, voff:%d, ulen:%d, vlen:%d\nq:%s\n", k, strand[k], u, v, al[k]->uoff, al[k]->voff, getUalignlen(al[k]), getValignlen(al[k]), al[k]->u);            
   showAlign(al[k], stdout);
#endif


  return backtracetable;
}



/*------------------------- localmulitsplicedmatrix --------------------------
 *    
 * @brief aligning a read a1 and its reverse complement a2 of length m to 
 * noofseqs loci
 * @author Steve Hoffmann 
 *   
 */

int*
localmultisplicedmatrix (void *space, symtype *a1, symtype *a2, Uint m,
    symtype **b, Uint *n, Uint *strand, Uint noofseqs, int indel, int trans,
    Sint (*sub)(symtype, symtype, void *), void *nfo, int ***lv, int ***lr, int ***lc){

  int i, j, k, q, cols=0, rows, size;
  Uint *cum=NULL;
//  unsigned char transition;
  int *L, **lmr=NULL, **lmv=NULL, **lmc=NULL;
  symtype cura, curb, r1, r2, l1, l2;

  rows = m + 1;

  cum = ALLOCMEMORY(space, NULL, Uint, noofseqs);
  lmr = ALLOCMEMORY(space, NULL, int*, noofseqs);
  lmv = ALLOCMEMORY(space, NULL, int*, noofseqs);
  lmc = ALLOCMEMORY(space, NULL, int*, noofseqs);

  for(k=0; k < noofseqs; k++) {
    cum[k] = cols;
    cols += n[k] + 1;
    lmr[k] = ALLOCMEMORY(space, NULL, int, rows);
    lmv[k] = ALLOCMEMORY(space, NULL, int, rows);
    lmc[k] = ALLOCMEMORY(space, NULL, int, rows);
    memset(lmv[k], 0, sizeof(int)*rows);
    memset(lmr[k], 0, sizeof(int)*rows);
    memset(lmc[k], 0, sizeof(int)*rows);
  }

  size = rows * cols;
  L = ALLOCMEMORY(space, NULL, int, size);
  memset(L, 0, sizeof(int)*size);

  for(k=0; k < noofseqs; k++) {
    for (i=1; i < rows; i++){
      lmv[k][i] = lmv[k][i-1];    
      lmr[k][i] = lmr[k][i-1];
      lmc[k][i] = lmc[k][i-1];

      for (j=1; j < n[k]+1; j++){
        l1 = 0;
        l2 = 0;
        r1 = 0;
        r2 = 0;

        if (strand[k] == 0){
          cura = a1[i-1];
          curb = b[k][j-1];
          if(j+1 < n[k]) {
            r1 = b[k][j];
            r2 = b[k][j+1];
          } 
          if(j > 2) {
            l1 = b[k][j-3];
            l2 = b[k][j-2];
          }
        } else {
          cura = a2[m-i];
          curb = b[k][n[k]-j];
          if(j > 2) {
            l1 = b[k][n[k]-j+1];
            l2 = b[k][n[k]-j+2];
          }
          if(n[k]-j > 1) {
            r1 = b[k][n[k]-j-2];
            r2 = b[k][n[k]-j-1];
          }
        }	

        MATRIX2D(L, cols, i, (cum[k] + j)) =
          MAX4(0,
              MATRIX2D(L, cols, (i-1), (cum[k] + j)) + indel,
              MATRIX2D(L, cols, i, (cum[k] + j-1)) + indel,
              MATRIX2D(L, cols, (i-1), (cum[k] + j-1)) + sub(cura, curb, nfo));

//        transition = 0;
        for(q=0; q < k; q++) {
//          Uint max=0;
//          max = arraymax(&MATRIX2D(L, cols, lmr[q][i-1], cum[q]), n[q]+1);

//          if(lmv[q][i-1] != MATRIX2D(L, cols, lmr[q][i-1], cum[q]+max)) {
//            fprintf(stderr, "k:%d, n[%d]+1:%d, lmv[%d][%d]:%d != matrix[lmr[%d][%d],%d+%d]:%d, lmc[q:%d][i-1:%d]:%d, i:%d, j:%d\n", 
//                k, q, n[q]+1, q, i-1, lmv[q][i-1], q, i-1, cum[q], max,
//                MATRIX2D(L, cols, lmr[q][i-1], cum[q]+max), q, i-1, lmc[q][i-1], i, j);
//            return NULL;
//          }
        
          //if strand[q] != strand[k]  splicescore = -lms[q][i-1]
          MATRIX2D(L, cols, i, (cum[k] + j)) =
            MAX(
                MATRIX2D(L, cols, i, (cum[k] + j)), 
                lmv[q][i-1] + sub(cura, curb, nfo) + trans + splicescore(l1, l2, strand[k], 1)) ;
        }


        //if strand -  acc1, acc2
        //if strand +  don1, don2
        if (MATRIX2D(L, cols, i, cum[k]+j) + splicescore(r1, r2, strand[k], 0) > lmv[k][i-1]) {
          lmv[k][i] = MATRIX2D(L, cols, i, cum[k]+j) + splicescore(r1, r2, strand[k], 0);
          lmr[k][i] = i;
          lmc[k][i] = j;
        }
      }
    }
  }

  FREEMEMORY(space, cum);
  *lv = lmv;
  *lr = lmr; 
  *lc = lmc;
  return L;
}


/*----------------------- localmultisplicedtraceback -------------------------
 *    
 * @brief tracing localmultisplicedmatrix back producing noofseqs split aligns
 * @author Steve Hoffmann 
 *   
 */

void
localmultisplicedtraceback (void *space, int *M, symtype *a1, symtype *a2, Uint m, 
    symtype **b, Uint* n, Uint *strand, Uint noofseqs, int indel, int trans,
    Sint (*sub)(symtype, symtype, void *), void *nfo,
    Alignment **al, int **lmv, int **lmr, int **lmc){

  Uint u=1, v=1, cur, cols=0, start, p, q, k,  
    *cum, xmax=0;
  int x;
  symtype cura, curb, l1, l2;//, r1, r2;

  cum = ALLOCMEMORY(space, NULL, Uint, noofseqs);
  for(k=0; k < noofseqs; k++) {
    al[k]->uoff = 0;
    al[k]->voff = 0;
    cum[k] = cols;
    cols += n[k] + 1;
  }

  start = arraymax(&MATRIX2D(M, cols, 0, 0), (cols*(m+1)));

  p = start / (cols);
  q = start % (cols);

  k=0;
  while(k+1 < noofseqs && q > cum[k+1]) k++;
  
  u = p;
  v = q-cum[k];

//  fprintf(stderr, "\nstarting on seq: %d (%d), u:%d, v:%d\n", k, strand[k], u, v);
  //if(u < m) {
//    for(i=u; i < m; i++) insertEop(al[k], Insertion);
  //}
  
  while(u > 0 && v > 0){
    cur = MATRIX2D(M, cols, u, (cum[k]+v));
//    fprintf(stderr, "continue on %d with matrix[%d,%d+%d=%d]:%d strand:%d\n", k, u, cum[k], v, cum[k]+v, MATRIX2D(M, cols, u, (cum[k]+v)), strand[k]);
 
    if (MATRIX2D(M, cols, (u-1), (cum[k]+v)) + indel == cur){
      insertEop(al[k], Insertion);
      assert(u > 0);
      u--;
    } else {
      if (MATRIX2D(M, cols, u, (cum[k]+v-1)) + indel == cur){
        insertEop(al[k], Deletion);
        assert(v > 0);
        v--;
      } else {
        if (MATRIX2D(M, cols, u, (cum[k]+v))) {
          l1 = 0;
          l2 = 0;
          /*
          r1 = 0;
          r2 = 0;*/

          if (strand[k] == 0){
            cura = a1[u-1];
            curb = b[k][v-1];
            
            /* not used
            if(v+1 < n[k]) {
              r1 = b[k][v];
              r2 = b[k][v+1];
              } */
            if(v > 2) {
              l1 = b[k][v-3];
              l2 = b[k][v-2];
            }
          } else {
            cura = a2[m-u];
            curb = b[k][n[k]-v];
            if(v > 2) {
              l1 = b[k][n[k]-v+1];
              l2 = b[k][n[k]-v+2];
            }
            /* not used
            if(n[k]-v > 1) {
              r1 = b[k][n[k]-v-2];
              r2 = b[k][n[k]-v-1];
              }*/
          }

//        fprintf(stderr, "sub[%c,%c]=%d\n", cura, curb, sub(cura,curb,nfo));
          if (MATRIX2D(M, cols, (u-1), (cum[k]+v-1)) + 
              sub(cura, curb, nfo) == cur){
            insertEop(al[k], Replacement);
            assert(u > 0 && v > 0);
            u--; v--;
          } else {

            if(k==0) {
              if (strand[k] == 0){
                al[k]->uoff = u;
                al[k]->voff = v;
                revMeops(al[k]);
              } else {  
                assert(m >= u+getUalignlen(al[k]));  
                al[k]->uoff = m-u-getUalignlen(al[k]);
                assert(n[k] >= v+getValignlen(al[k]));              
                al[k]->voff = n[k]-v-getValignlen(al[k]);
              }
              break;
            }

            x = -1;
            for(q=k; q > 0; q--) {            
/*            Uint max;
              max = arraymax(&MATRIX2D(M, cols, lmr[q-1][u-1], cum[q-1]), n[q-1]+1);
              assert(MATRIX2D(M, cols, lmr[q-1][u-1],cum[q-1]+max) 
                  == lmv[q-1][u-1]);
               fprintf(stderr, "transition: matrix[lmr[%d][%d],%d]=%d, sub[%c,%c]=%d\n", q-1, u-1, max, lmv[q-1][u-1], cura, curb, sub(cura,curb,nfo));
*/
              if (
                  //MATRIX2D(M, cols, lmr[q-1][u-1], cum[q-1]+lmc[q-1][u-1]) + 
                  lmv[q-1][u-1] +  splicescore(l1, l2, strand[k], 1) +
                  sub(cura, curb, nfo) + trans == cur) {
                x=q-1;
                xmax = lmc[q-1][u-1];
              }
            }
         
//            fprintf(stderr, "matrix:U[%d,%d+%d=%d]:%d, k:%d, q:%d, x:%d -> [%d,%d]\n", u, cum[k], v, cum[k]+v, cur, k, q, x, lmr[x][u-1], xmax);
            
            assert(x > -1);
            insertEop(al[k], Replacement);
            assert(u > 0 && v > 0);
            u--; v--;

            if (strand[k] == 0){
              al[k]->uoff = u;
              al[k]->voff = v;
              revMeops(al[k]);
            } else {
              assert(m >= u+getUalignlen(al[k]));
              assert(n[k] >= v+getValignlen(al[k]));
              al[k]->uoff= m-u-getUalignlen(al[k]);
              al[k]->voff = n[k]-v-getValignlen(al[k]);
            }
/*            
            fprintf(stderr, "transition from: m:%d, u:%d, v:%d, ulen:%d, vlen:%d -> u': %d v':%d \n", 
                m, u, v, getUalignlen(al[k]), getValignlen(al[k]), lmr[x][u], xmax);
            
            showAlign(al[k], stdout);
            fprintf(stdout, "\n");
*/
            u = lmr[x][u];
            v = xmax;

            assert(u >= 0 && v >= 0);
            k = x;	    
          }
        } else {
          //fprintf(stderr, "reversing end\n");
          if (strand[k] == 0 && k > 0){
            al[k]->uoff = u;
            al[k]->voff = v;
            revMeops(al[k]);
          } else { 
            assert(m >= u+getUalignlen(al[k]));
            al[k]->uoff = m-u-getUalignlen(al[k]);
            assert(n[k] >= v+getValignlen(al[k]));
            al[k]->voff = n[k]-v-getValignlen(al[k]);
          }
          break;
        }
      }
    }
  }    

//  fprintf(stderr, "starting on seq: %d (%d), u:%d, v:%d\n", k, strand[k], u, v);

  if (k==0){
    if (strand[k] == 0){
      al[k]->uoff = u;
      al[k]->voff = v;
      revMeops(al[k]);
    } else {
  /*    fprintf(stdout, "finalize alignment at u:%d, v:%d, ulen:%d, vlen:%d\n", u, v, getUalignlen(al[k]), getValignlen(al[k]));*/
      assert(m >= u+getUalignlen(al[k]));
      al[k]->uoff = m-u-getUalignlen(al[k]);
      assert(n[k] >= v+getValignlen(al[k]));
      al[k]->voff = n[k]-v-getValignlen(al[k]);
    }
  }

  FREEMEMORY(space, cum);
  return;
}


/*---------------------------- localsplicedmatrix ----------------------------
 *    
 * @brief generate local spliced alignment
 * @author Steve Hoffmann 
 *   
 */
 

int*
localsplicedmatrix (void *space, symtype *a1, symtype *a2, Uint m,
    symtype *b1, Uint n1, symtype *b2, Uint n2, 
    Uint strand1, Uint strand2, int indel,
    Sint (*sub)(symtype, symtype, void *), void *nfo, int **lv, int **lr){
  int i, j, k, max, cols1, cols2, cols, rows, size, maxval=0;
  int *L, *lmr=NULL, *lmv=NULL;
  symtype cura, curb;

  rows = m + 1;
  cols1 = n1 + 1;
  cols2 = n2 + 1;
  cols = cols1 + cols2;

  size = rows * cols;
  L = ALLOCMEMORY(space, NULL, int, size);
  memset(L, 0, sizeof(int)*size);
  
  lmr = ALLOCMEMORY(space, NULL, int, rows);
  lmv = ALLOCMEMORY(space, NULL, int, rows);
  memset(lmv, 0, sizeof(int)*rows);
  memset(lmr, 0, sizeof(int)*rows);

  for (i = 1; i < rows; i++){
    lmv[i] = lmv[i-1];    
    lmr[i] = lmr[i-1];

    for (j = 1; j < cols1; j++){
      if (strand1 == 0 || strand1 == strand2){
        cura = a1[i-1];
        curb = b1[j-1];
      }
      else {
        cura = a1[m-i];
        curb = b1[n1-j];
      }	
      
      MATRIX2D(L, cols, i, j) = 
        MAX4(0,
            MATRIX2D(L, cols, (i-1), j) + indel,
            MATRIX2D(L, cols, i, (j-1)) + indel,
            MATRIX2D(L, cols, (i-1), (j-1)) + sub(cura, curb, nfo));	
      
      if (MATRIX2D(L, cols, i, j) > lmv[i-1]) {
        lmv[i] = MATRIX2D(L, cols, i, j);
        lmr[i] = i;
      } 
    }
      
    max = arraymax(&MATRIX2D(L, cols, lmr[i-1], 0), cols1);
    assert(lmv[i-1] == MATRIX2D(L, cols, lmr[i-1],max));

    for (k = 1; k < cols2; k++){
      if (strand2 == 0 || strand1 == strand2){
        cura = a2[i-1];
        curb = b2[k-1];
      }
      else {
        cura = a2[m-i];
        curb = b2[n2-k];
      }

      MATRIX2D(L, cols, i, (cols1 + k)) =
        MAX5(0,
            MATRIX2D(L, cols, (i-1), (cols1+k)) + indel,
            MATRIX2D(L, cols, i, (cols1+k-1)) + indel,
            MATRIX2D(L, cols, (i-1), (cols1+k-1)) + sub(cura, curb, nfo),
            MATRIX2D(L, cols, lmr[i-1], max) + sub(cura, curb, nfo));

      if(MATRIX2D(L, cols, i, (cols1+k)) > maxval) {
        maxval = MATRIX2D(L, cols, i, (cols1+k));
      }
    }
  }

  *lv = lmv;
  *lr = lmr; 
  return L;
}


/*--------------------------- localsplicetraceback ---------------------------
 *    
 * @brief tracing back a local spliced alignment
 * @author Steve Hoffmann 
 *   
 */
 
void
localsplicedtraceback (void *space, int *M, symtype *a1, symtype *a2, Uint m, 
    symtype *b1, Uint n1,
    symtype *b2, Uint n2, Uint strand1, Uint strand2, int indel,
    Sint (*sub)(symtype, symtype, void *), void *nfo,
    Alignment *al1, Alignment *al2, int* lmv, int *lmr){

  Uint i=1, j=1, u=1, v=1, cur, cols1, cols2, cols, start, p, q;
  int max;
  BOOL seq1;
  symtype cura, curb;

  al1->uoff = 0;
  al1->voff = 0;
  al2->uoff = 0;
  al2->voff = 0;

  cols1 = n1 + 1;
  cols2 = n2 + 1;
  cols = cols1 + cols2;

  start = arraymax(&MATRIX2D(M, cols, 0, 0), (cols*(m+1)));
  
  p = start / (cols);
  q = start % (cols);
  
  if (q > cols1) {
    u = p;
    v = q-cols1;
    seq1 = 0;
  } else {
    i = p;
    j = q;
    seq1 = 1;
  }

  while(i > 0 && j > 0 && u > 0 && v > 0){

    if (seq1){
      cur = MATRIX2D(M, cols, i, j);

      if (MATRIX2D(M, cols, (i-1), j) + indel == cur){
        insertEop(al1, Deletion);
        i--;
      } else {
        if (MATRIX2D(M, cols, i, (j-1)) + indel == cur){
          insertEop(al1, Insertion);
          j--;
        } else {
          if (MATRIX2D(M, cols, (i), (j))) {
            if (strand1 == 0 || strand1 == strand2){
              cura = a1[i-1];
              curb = b1[j-1];
            } else {
              cura = a1[m-i];
              curb = b1[n1-j];
            }
            assert(MATRIX2D(M, cols, (i-1), (j-1)) 
                + sub(cura, curb, nfo) == cur);
            insertEop(al1, Replacement);
            i--; j--; 
          } else {
            break;
          }
        }
      }
    } else {

      cur = MATRIX2D(M, cols, u, (cols1+v));
      if (MATRIX2D(M, cols, (u-1), (cols1+v)) + indel == cur){
        insertEop(al2, Deletion);
        u--;
      } else {
        if (MATRIX2D(M, cols, u, (cols1+v-1)) + indel == cur){
          insertEop(al2, Insertion);
          v--;
        } else {
          if (MATRIX2D(M, cols, u, (cols1+v))) {
            if (strand2 == 0 || strand1 == strand2){
              cura = a2[u-1];
              curb = b2[v-1];
            } else {
              cura = a2[m-u];
              curb = b2[n2-v];
            }
            if (MATRIX2D(M, cols, (u-1), (cols1+v-1)) + 
                sub(cura, curb, nfo) == cur){
              insertEop(al2, Replacement);
              u--; v--;
            } else {
              max = arraymax(&MATRIX2D(M, cols, lmr[u-1], 0), cols1);
              assert(MATRIX2D(M, cols, lmr[u-1], max) + 
                  sub(cura, curb, nfo) == cur);
              
              insertEop(al2, Replacement);
              u--; v--;
              al2->uoff = u;
              al2->voff = v;
 
              if (strand2 == 0 || strand1 == strand2){
                revMeops(al2);
              }
              i = lmr[u];
              j = max;
              seq1 = 1;	    
            }
          } else {

            al2->uoff = u;
            al2->voff = v;
            if (strand2 == 0 || strand1 == strand2){
              revMeops(al2);
            }
            break;
          }
        }
      }
    }    
  }
  if (seq1){
    al1->uoff = i;
    al1->voff = j;
    if (strand1 == 0 || strand1 == strand2){
      revMeops(al1);
    }
  } 

  return;
}


/*------------------------------ splicedmatrix -------------------------------
 *    
 * @brief calculate semi-global alignment matrix
 * @author Christian Otto
 *   
 */
 
int*
splicedmatrix (void *space, symtype *a1, symtype *a2, Uint m,
    symtype *b1, Uint n1, symtype *b2, Uint n2, 
    Uint strand1, Uint strand2, int indel,
    Sint (*sub)(symtype, symtype, void *), void *nfo){
  int i, j, k, max, cols1, cols2, cols, rows, size;
  int *L;
  symtype cura, curb;

  rows = m + 1;
  cols1 = n1 + 1;
  cols2 = n2 + 1;
  cols = cols1 + cols2;

  size = rows * cols;
  L = ALLOCMEMORY(space, NULL, int, size);
  L = memset(L, 0, sizeof(int)*size);
  for (i = 1; i < rows; i++){
    MATRIX2D(L, cols, i, 0) = i * indel;
    MATRIX2D(L, cols, i, cols1) = i * indel;
  }

  for (i = 1; i < rows; i++){
    for (j = 1; j < cols1; j++){
      if (strand1 == 0 || strand1 == strand2){
        cura = a1[i-1];
        curb = b1[j-1];
      }
      else {
        cura = a1[m-i];
        curb = b1[n1-j];
      }	
      MATRIX2D(L, cols, i, j) = 
        MAX3(MATRIX2D(L, cols, (i-1), j) + indel,
            MATRIX2D(L, cols, i, (j-1)) + indel,
            MATRIX2D(L, cols, (i-1), (j-1)) + sub(cura, curb, nfo));	
    }
    max = arraymax(&MATRIX2D(L, cols, (i-1), 0), cols1);
    for (k = 1; k < cols2; k++){
      if (strand2 == 0 || strand1 == strand2){
        cura = a2[i-1];
        curb = b2[k-1];
      }
      else {
        cura = a2[m-i];
        curb = b2[n2-k];
      }
      MATRIX2D(L, cols, i, (cols1 + k)) =
        MAX4(MATRIX2D(L, cols, (i-1), (cols1+k)) + indel,
            MATRIX2D(L, cols, i, (cols1+k-1)) + indel,
            MATRIX2D(L, cols, (i-1), (cols1+k-1)) + sub(cura, curb, nfo),
            MATRIX2D(L, cols, (i-1), max) + sub(cura, curb, nfo));
    }
  }
  return L;
}

/*----------------------------- splicedtraceback -----------------------------
 *    
 * @brief tracing back semi-global spliced alignment
 * @author Christian Otto
 *   
 */
 
void
splicedtraceback (void *space, int *M, 
    symtype *a1, symtype *a2, Uint m, 
    symtype *b1, Uint n1,
    symtype *b2, Uint n2, 
    Uint strand1, Uint strand2, 
    int indel, Sint (*sub)(symtype, symtype, void *), void *nfo,
    Alignment *al1, Alignment *al2){

  Uint i, j, k, cur, cols1, cols2, cols;
  int max;
  BOOL seq1;
  symtype cura, curb;

  al1->uoff = 0;
  al1->voff = 0;
  al2->uoff = 0;
  al2->voff = 0;

  cols1 = n1 + 1;
  cols2 = n2 + 1;
  cols = cols1 + cols2;
  i = m;
  j = arraymax(&MATRIX2D(M, cols, i, 0), cols1);
  k = arraymax(&MATRIX2D(M, cols, i, cols1), cols2);

  if (MATRIX2D(M, cols, i, j) > MATRIX2D(M, cols, i, (cols1+k))){
    seq1 = 1;
  }
  else {
    seq1 = 0;
  }

  while(i > 0 && j > 0 && k > 0){

    if (seq1){
      cur = MATRIX2D(M, cols, i, j);
      if (MATRIX2D(M, cols, (i-1), j) + indel == cur){
        insertEop(al1, Deletion);
        i--;
      }
      else {
        if (MATRIX2D(M, cols, i, (j-1)) + indel == cur){
          insertEop(al1, Insertion);
          j--;
        }
        else {
          if (strand1 == 0 || strand1 == strand2){
            cura = a1[i-1];
            curb = b1[j-1];
          }
          else {
            cura = a1[m-i];
            curb = b1[n1-j];
          }
          assert(MATRIX2D(M, cols, (i-1), (j-1)) 
              + sub(cura, curb, nfo) == cur);
          
          insertEop(al1, Replacement);
          i--; j--; 
        }
      }
    } else {
      cur = MATRIX2D(M, cols, i, (cols1+k));
      if (MATRIX2D(M, cols, (i-1), (cols1+k)) + indel == cur){
        insertEop(al2, Deletion);
        i--;
      } else {
        if (MATRIX2D(M, cols, i, (cols1+k-1)) + indel == cur){
          insertEop(al2, Insertion);
          k--;
        } else {
          if (strand2 == 0 || strand1 == strand2){
            cura = a2[i-1];
            curb = b2[k-1];
          } else {
            cura = a2[m-i];
            curb = b2[n2-k];
          }
          if (MATRIX2D(M, cols, (i-1), (cols1+k-1)) + sub(cura, curb, nfo) == cur){
            insertEop(al2, Replacement);
            i--; k--;
          } else {
            max = arraymax(&MATRIX2D(M, cols, (i-1), 0), cols1);
            assert(MATRIX2D(M, cols, (i-1), max) 
                + sub(cura, curb, nfo) == cur);

            insertEop(al2, Replacement);
            i--; k--;
            
            al2->uoff = i;
            al2->voff = k;
            if (strand2 == 0 || strand1 == strand2){
              revMeops(al2);
            }
            j = max;
            seq1 = 1;	    
          }
        }
      }
    }    
  }
  if (seq1){
    al1->uoff = i;
    al1->voff = j;
    if (strand1 == 0 || strand1 == strand2){
      revMeops(al1);
    }
  }

  
  return;
}

/*------------------------------- swgaplesstraceback -----------------------------
 *     
 * traceback to find optimal local alignment path
 *    
 */

  int*
swgaplesstraceback (void *space, int *M,  
    symtype *a, Uint m, symtype *b, Uint n, 
    Sint (*sub)(symtype, symtype, void *), void *nfo, int* alignsize)
{
  Uint i, j, ncol, cur, start;
  int *align = NULL;

  *alignsize = 0;
  ncol = (n+1);
  start = arraymax(M, (m+1)*ncol);
  i = start / ncol;
  j = start % ncol;

  while(i > 0 && j > 0) {

    cur = MATRIX2D(M, ncol, i, j);
    if (cur==0)
      return align;

    align = ALLOCMEMORY(space, align, int, *alignsize+2);
    align[*alignsize]=i;
    align[*alignsize+1]=j;
    *alignsize +=2;

    i--; j--;

  }

  return align;
}




