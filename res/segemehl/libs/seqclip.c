
/*
 *  seqclip.c
 *  
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 24.04.2010 21:24:34 CEST
 *  
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include "biofiles.h"
#include "basic-types.h"
#include "mathematics.h"
#include "memory.h"
#include "alignment.h"
#include "time.h"
#include "sw.h"
#include "seqclip.h"
#include "bitvectoralg.h"

static int clpswscr[]={1,-2};   
static int edstscr[]={1,0};
static char *polyA = 
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
/*---------------------------- bl_seqclipSoft3Prime -----------------------------
 *    
 * @brief clipping adapter from 3'end
 * @return length of clipped part
 * @author Steve Hoffmann 
 *   
 */

Uint
bl_seqclipSoft3Prime(void *space, char *s, Uint len, 
    char *C, Uint clen, Uint minclipacc, Uint polyAlen) 
{
  int *M;
  int mrgn;
  Uint rlen=0, rmargin = 0, rend, valignlen, adapteracc=100;
  int polyAScore=0, adapterScore=0, adapterMatchLen=0;
  Alignment *al;


//  fprintf(stderr, "softclipprime");
  mrgn = MIN((((float)clen)*.2) + clen, len);

//  fprintf(stderr, "mrgn: %d, len:%d, len-mrgn:%d\n", mrgn, len, len-mrgn);
  M = swmatrix(space, C, clen, 
      &s[len-mrgn], mrgn, -2, constscr_Nmatch, clpswscr);
  
  al = ALLOCMEMORY(space, NULL, Alignment, 1);
  initAlignment(al, C, clen, 0, &s[len-mrgn], mrgn, 0); 
  
  swtraceback(space, M, C, clen, 
      &s[len-mrgn], mrgn, -2, constscr_Nmatch, clpswscr, al);

  valignlen = getValignlen(al);

  rmargin = len > 4 + (int)(((double)valignlen)*0.2) ? 
    len - 4 - (int)(((double)valignlen)*0.2) : 0;  

  getSoftClipScores(al, polyAlen, edstscr, 0, 
      &polyAScore, &adapterScore, &adapterMatchLen);
  
  if(adapterMatchLen > 3) { 
    adapteracc = (int)(((double)adapterScore/(double)adapterMatchLen)*100);
  } else {
    adapteracc = 100;
  }

  if(len - mrgn + al->voff + valignlen >= rmargin && 
      (polyAScore >= 5 || adapterScore >= 8) &&
      adapteracc >= minclipacc) 
  {  
        rend = (mrgn > al->voff+ valignlen) ? mrgn - al->voff - valignlen : 0;
        rlen = getValignlen(al) + rend;
  }

  
 //  fprintf(stderr, "c:%s\n", s);
//  fprintf(stderr, 
//  "adapteracc:%d, minclipacc: %d, polyAScore %d, valginlen: %d, mrgn: %d, len: %d, off:%d -> rend: %d, rlen:%d\n", 
//      adapteracc, minclipacc, polyAScore, valignlen, mrgn, len, al->voff, rend, rlen);
//   showAlign(al, stderr);  


  
  wrapAlignment(al);
  FREEMEMORY(space, M);
  FREEMEMORY(space, al);
  return rlen;
}


/*---------------------------- bl_seqclipSoft5Prime -----------------------------
 *    
 * @brief clipping adapter from 5' end
 * @author Steve Hoffmann 
 *   
 */

Uint
bl_seqclipSoft5Prime(void *space, char *s, Uint len, 
    char *C, Uint clen, Uint minclipscr) 
{
  int *M;
  int mrgn, allen=0; 
  Alignment *al;


  mrgn = MIN((((float)clen)*.2) + clen, len);
  
  M = swmatrix(space, C, clen, s, mrgn, -2, constscr_Nmatch, clpswscr);
  
  al = ALLOCMEMORY(space, NULL, Alignment, 1);
  initAlignment(al, C, clen, 0, s, mrgn, 0); 
  
  swtraceback(space, M, C, clen, 
      s, mrgn, -2, constscr_Nmatch, clpswscr, al);
  
  if(getAlignScore(al, edstscr, 0) >= minclipscr) {
    allen = getValignlen(al);
  }

 // showAlign(al, stdout);
  
  wrapAlignment(al);
  FREEMEMORY(space, M);
  FREEMEMORY(space, al);

  return allen;
}

/*----------------------------- bl_seqclipPolyA ------------------------------
 *    
 * @brief clipping polyA tails 
 * @author Steve Hoffmann 
 *   
 */

Uint
bl_seqclipPolyA(void *space, char *s, Uint len, char *clp, Uint clen) 
{
  int *M;
  int mrgn;
  Uint polyAlen;
  char *polyAseq;
  Uint rlen = len, rmargin = 0;
  Alignment *al;

  if (len < 10) return len;

  if(clp) {
    polyAlen = strlen(polyA) + clen;
    polyAseq = ALLOCMEMORY(space, NULL, char, polyAlen+1);
    memset(polyAseq, 0, polyAlen+1);
    memmove(polyAseq, polyA, strlen(polyA));
    memmove(&polyAseq[strlen(polyA)], clp, clen);
  } else {
    polyAseq = polyA;
    polyAlen = strlen(polyA);
  }

  al = ALLOCMEMORY(space, NULL, Alignment, 1);
  mrgn =((float)len*.5);

  initAlignment(al, polyAseq, polyAlen, 0, &s[len-mrgn], mrgn, 0); 

  M = swmatrix(space, polyAseq, polyAlen, &s[len-mrgn], 
      mrgn, -3, constscr_Nmatch, clpswscr);

  swtraceback(space, M, polyAseq, polyAlen, &s[len-mrgn], 
      mrgn, -3, constscr_Nmatch, clpswscr, al);
  
  rmargin = len > 4 + (int)((double)getValignlen(al)*0.2) ? 
    len - 4 - (int)(((double)getValignlen(al))*0.2) : 0;  

  if(len - mrgn + al->voff + getValignlen(al) >= rmargin) {
    rlen = len - mrgn + al->voff;
  }
 
  //showAlign(al, stdout);

  FREEMEMORY(space, M);
  wrapAlignment(al);

  return rlen;
}

/*--------------------------- bl_seqclipHard3Prime ---------------------------
 *    
 * @brief hard clipping of 3'prime end
 * @author Steve Hoffmann 
 *   
 */

Uint
bl_seqclipHard3Prime(Uint len, Uint clen) {
  Uint mrgn;

  mrgn = MIN(len, clen);
  return len-mrgn;
}


/*--------------------------- bl_seqclipHard5Prime ---------------------------
 *    
 * @brief hard clipping of 5'prime end
 * @author Steve Hoffmann 
 *   
 */

char*
bl_seqclipHard5Prime(char *s, Uint len, Uint clen) {
  Uint mrgn;

  mrgn = (clen >= len) ? 0 : clen;
  return &s[mrgn];
}


/*----------------------------- bl_seqClipDecode -----------------------------
 *    
 * @brief find DNA sequence of length len for a code to the base of 5
 * @author Steve Hoffmann 
 *   
 */

  char*
bl_seqclipDecode (Uint code, Uint len)
{
  Uint i=0, n, r;
  char *seq, ch;


  seq = ALLOCMEMORY(space, NULL, char, len+1);
  memset(seq, 'A', sizeof(char)*len);
  seq[len] = 0;
  n = code;

  while(n > 0) {
    r = n % 5;
    n = n / 5;

    switch(r) {
      case 0:
        ch = 'A';
        break;
      case 1:
        ch= 'C';
        break;
      case 2:
        ch = 'G';
        break;
      case 3:
        ch = 'T';
        break;
      default:
        ch = 'N';
    }
    seq[i++] = ch; 
  }

  return seq;
}


/*---------------------------- bl_seqclipGetCode -----------------------------
 *    
 * @brief find code for a DNA sequence of length len
 * @author Steve Hoffmann 
 *   
 */

  Uint
bl_seqclipGetCode (char *seq, Uint len)
{
  Uint i, z, sum=0;

  for(i=0; i < len; i++) {
    switch(seq[i]) {
      case 'A':
        z=0;
        break;
      case 'C':
        z=1;
        break;
      case 'G':
        z=2;
        break;
      case 'T':
        z=3;
        break;
      default:
        z=4;
        break;
    }
    sum += z *pow(5, i);
  }
  return sum;
}



/*---------------------------------- bl_lcs ----------------------------------
 *    
 * @brief get 
 * @author Steve Hoffmann 
 *   
 */
 
Uint
bl_lcsub (void *space, char *s1, Uint l1, char *s2, Uint l2, Uint *l, Uint *r)
{
  Uint i, j, maxlen=0, maxi=0, maxj=0;
  Uint *M;

  M = ALLOCMEMORY(space, NULL, Uint, (l1+1)*(l2+1));
  memset(M, 0, sizeof(Uint)*((l1+1)*(l2+1)));

  for(i=0; i < l1; i++) {
    for(j=0; j < l2;  j++) {
      if(s1[i] == s2[j]) {
        if(i==0 || j==0) {
          MATRIX2D(M, l2+1, i, j) = 1;
        } else {
          MATRIX2D(M, l2+1, i, j) = 
            MATRIX2D(M, l2+1, i-1, j-1) + 1;
        }
        if(MATRIX2D(M, l2+1, i, j) > maxlen) {
          maxlen = MATRIX2D(M, l2+1, i, j);
          maxi = i;
          maxj = j;
        }
      }
    }
  }

  *l = maxi;
  *r = maxj;

  FREEMEMORY(space, M);
  return maxlen;
}



/*-------------------------------- bl_greedy ---------------------------------
 *    
 * @brief assemble sequences greedily
 * @author Steve Hoffmann 
 *   
 */

char*
bl_Find3PrimeGreedy (void *space, PairUint *list, Uint len, Uint ws) {

  Uint i=0, j=0, l, r, k, llen, rlen, start, adapterlen=0, 
            lcs, minovl, last;
  Uint *Lovl, *Rovl, *Lchain, *Rchain, *use;
  char *seq1, *seq2;
  char *adapter = NULL;

  Lovl = ALLOCMEMORY(space, NULL, Uint, len*len);
  Rovl = ALLOCMEMORY(space, NULL, Uint, len*len);
  use = ALLOCMEMORY(space, NULL, Uint, len);
  Lchain = ALLOCMEMORY(space, NULL, Uint, len);
  Rchain = ALLOCMEMORY(space, NULL, Uint, len);

  memset(Lovl, 0,sizeof(Uint)*(len*len));
  memset(Rovl, 0,sizeof(Uint)*(len*len));
  minovl = (ws/3) * 2;


  for(i=0; i < len; i++) {
    seq1 = bl_seqclipDecode(list[i].b, ws);
    for(j=0; j < len; j++) {
      if(i != j) {
        seq2 = bl_seqclipDecode(list[j].b, ws);
        lcs = bl_lcsub(space, seq1, ws, seq2, ws, &l, &r);
        //lcs is overlap
        if((l == ws-1 && r == lcs-1) || (r == ws-1 && l == lcs-1)) {
          if(l > r) {
            MATRIX2D(Rovl, len, i, j) = lcs;
          } else {
            MATRIX2D(Lovl, len, i, j) = lcs;
          }
        }
        FREEMEMORY(space, seq2);
      }
    }
    FREEMEMORY(space, seq1);
  }

  for(i=0; i < 5; i++) {
    memset(use, 0, sizeof(Uint)*(len));
    use[i] = 1;
    last = i;
    rlen = 0;
    Rchain[rlen] = -1;
    
    while(1) {
      for(l=0, j=0; j < len; j++) {
        if (MATRIX2D(Rovl, len, last, j) > l && 
            MATRIX2D(Rovl, len, last, j) >= minovl &&
            !use[j]) {
          if(Rchain[rlen] != -1) {
            use[Rchain[rlen]] = 0;
          }
          Rchain[rlen] = j;
          l = MATRIX2D(Rovl, len, last, j);
          use[j] = 1;
        }
      }
      if(Rchain[rlen] == -1) break;
      last = Rchain[rlen++];
      Rchain[rlen] = -1;
    }

    last = i;
    llen = 0;
    Lchain[llen] = -1;
    
    while(1) {
      for(l=0, j=0; j < len; j++) {
        if (MATRIX2D(Lovl, len, last, j) > l && 
            MATRIX2D(Lovl, len, last, j) >= minovl &&
            !use[j]) {
          if(Lchain[llen] != -1) {
            use[Lchain[llen]] = 0;
          }
          Lchain[llen] = j;
          l = MATRIX2D(Lovl, len, last, j);
          use[j] = 1;
        }
      }
      if(Lchain[llen] == -1) break;
      last = Lchain[llen++];
      if (llen >= len) break;
      Lchain[llen] = -1;
    }
 

    seq2 = ALLOCMEMORY(space, NULL, char, len*ws);
    memset(seq2, 0, len*ws);
    
    start = 0;
    for(j=llen; j >= 1; j--) {
      seq1 = bl_seqclipDecode(list[Lchain[j-1]].b, ws);
    
      if(j==llen) { 
        k = ws; 
      } else {
        k = MATRIX2D(Lovl, len, Lchain[j-1], Lchain[j]);
      }

      memmove(&seq2[start+(ws-k)], seq1, ws);
      start += ws-k;
      FREEMEMORY(space, seq1);
    }
 
    if(Lchain[0] != -1) {
    seq1 = bl_seqclipDecode(list[i].b, ws);
    k = MATRIX2D(Lovl, len, i, Lchain[0]);
    memmove(&seq2[start+(ws-k)], seq1, ws);
    start += ws-k;
    FREEMEMORY(space, seq1);
    }

    if(Rchain[0] != -1) {
    seq1 = bl_seqclipDecode(list[Rchain[0]].b, ws);
    k = MATRIX2D(Rovl, len, i, Rchain[0]);
    memmove(&seq2[start+(ws-k)], seq1, ws);
    start += ws-k;
    FREEMEMORY(space, seq1);
    }

    for(j=0; j+1 < rlen; j++) {
      seq1 = bl_seqclipDecode(list[Rchain[j+1]].b, ws);
      k = MATRIX2D(Rovl, len, Rchain[j], Rchain[j+1]);
      memmove(&seq2[start+(ws-k)], seq1, ws);
      start += ws-k;  
      FREEMEMORY(space, seq1);
    }

//    fprintf(stderr, "seq: %s\n", seq2);
    if(adapterlen < strlen(seq2)) {
      if(adapter) 
        FREEMEMORY(space, adapter);
      adapter = seq2;
      adapterlen = strlen(seq2);
    } else {
      FREEMEMORY(space, seq2);
    }

    seq2 = NULL;
  }

  FREEMEMORY(space, Rovl);
  FREEMEMORY(space, Lovl);
  FREEMEMORY(space, use);
  FREEMEMORY(space, Lchain);
  FREEMEMORY(space, Rchain);

  return adapter;
}



/*---------------------- bl_seqclipFind3PrimeUpdateBest ----------------------
 *    
 * @brief helper function to update the list of most frequent motives
 * @author Steve Hoffmann 
 *   
 */
PairUint *
bl_seqclipFind3PrimeUpdateBest(PairUint *list, Uint len, 
    Uint code, Uint value) {

  Uint i;
  unsigned char shift=0;

  for(i=0; i < len; i++) {

    //eliminate old instance
    if(shift && list[i].b == code) {
      memmove(&list[i], &list[i+1], sizeof(PairUint)*(len-i));
      break;
    }
    
    //no higher ranking list elem < value -> update and cancel
    if(!shift && list[i].a < value && list[i].b == code) {
      list[i].a = value;
      break;
    }

    //a higher ranking list elem < value -> push list down if > 0
    //create a new instance and raise shift flag
    if(!shift && list[i].a < value && list[i].b != code) {

      if(list[i].a > 0) {
        list = ALLOCMEMORY(space, list, PairUint, len+1);
        memmove(&list[i+1], &list[i], sizeof(PairUint)*(len-i));
        shift = 1;
      }

      list[i].a = value;
      list[i].b = code;

      if(shift == 0)
        break;
    }
  }
  
  list = ALLOCMEMORY(space, list, PairUint, len);
  return list;
}

/*--------------------------- bl_seqclipFind3Prime ---------------------------
 *    
 * @brief find the 3 prime adapter
 * @brief find most common sequence of size ws in frame of fs
 * @author Steve Hoffmann 
 *   
 */

char*
bl_seqclipFind3Prime (void *space, fasta_t *set, Uint samplesize, Uint fs, int ws)
{

  char *curseq, *curframe, *curwin;
  uint16_t *C;
  double curent = 0;
  Uint i, elem, size, curlen, *tab, curcode=0, bestvalue=0;
  int curfs, j;
  //Uint bestcode=0;
  PairUint *B;

  //WARN: samplesize > noofseqs
  if(samplesize) {
    size = MIN(samplesize, set->noofseqs);
    srand(time(NULL));
  } else {
    size = set->noofseqs;
  }

  assert(ws <= fs);
  C = ALLOCMEMORY(space, NULL, uint16_t, pow(5, ws));
  B = ALLOCMEMORY(space, NULL, PairUint, 100+1);
  tab = ALLOCMEMORY(space, NULL, Uint, 256);
  memset(C, 0, sizeof(uint16_t)*pow(5,ws));
  memset(B, 0, sizeof(PairUint)*101);
  memset(tab, 0, sizeof(Uint)*256);
  
  tab['A'] = 1;
  tab['C'] = 2;
  tab['G'] = 3;
  tab['T'] = 4;


  for(i=0; i < size; i++) {
    
/*  too costly! new scheme: sample chunks of larger size  
    if(samplesize) {
      elem = rand() % set->noofseqs;
    } else {
      elem = i;
    }
*/
    elem = i;

    curseq = bl_fastaGetSequence(set, elem);
    curlen = bl_fastaGetSequenceLength(set, elem);
 
    // catch case where fs <= curlen is never true
    if(fs >= curlen) {
      curfs = curlen-1;
    } else {
      curfs = fs;
    }

    //WARN: fs > readsize   
    curframe = &curseq[curlen-curfs-1];

    for(j=0; j < curfs-ws+2; j++) {
      curwin = &curframe[j];
      curent = shannonentropy(space, curwin, ws, 6, tab);

      if (curent < 1.6) continue;

      curcode = bl_seqclipGetCode(curwin, ws);
      C[curcode]++;

      if(B[100-1].a < C[curcode]) {
        B = bl_seqclipFind3PrimeUpdateBest(B, 100, curcode, C[curcode]);
      }

      if (bestvalue < C[curcode]) {
        bestvalue = C[curcode];
 //       bestcode = curcode;
      }
      
      if(C[curcode] == 65500) 
        break;
    }
    if(C[curcode] == 65500) 
      break;
  }
/*
   {
   Uint k=0, l=0, r=0;
   char *nextseq;
   for(k=0; k < 99; k++) {
    curseq = bl_seqclipDecode(B[k].b, ws);
    nxtseq = bl_seqclipDecode(B[k+1].b, ws);
    bl_lcsub(space, curseq, 10, nxtseq, 10, &l, &r);
    fprintf(stderr, "best sequence[%d]: %s, val:%d\n", k, curseq, B[k].a);
    fprintf(stderr, "shannon: %f\n", shannonentropy(space, curseq,ws,5,tab));
    FREEMEMORY(space, curseq);
    FREEMEMORY(space, nxtseq);
  }

  curseq = bl_seqclipDecode(bestcode, ws);
  fprintf(stderr, "best sequence: %s\n", curseq);
  }
*/
  
  curseq = bl_Find3PrimeGreedy(space, B, 100, ws);
  

  FREEMEMORY(space, C);
  FREEMEMORY(space, B);
  FREEMEMORY(space, tab);


  return curseq;
}


