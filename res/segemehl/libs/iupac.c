/**
 * iupac.c
 * declarations for IUPAC nucleotide code
 *
 * @author Christian Otto
 * @email christian@bioinf.uni-leipzig.de
 * @company Bioinformatics, University of Leipzig
 * @date Fri Jul 23 15:03:08 CEST 2010
 */

/*
 *  SVN
 *  Revision of last commit: $Rev: 149 $
 *  Author: $Author: steve $
 *  Date: $Date: 2010-09-14 05:45:04 -0400 (Tue, 14 Sep 2010) $
 *  Id: $Id: iupac.c 149 2010-09-14 09:45:04Z steve $
 *  Url: $URL: http://www2.bioinf.uni-leipzig.de/svn5/segemehl/libs/iupac.c $
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "basic-types.h"
#include "debug.h"
#include "iupac.h"

/* defined maximal allowed symbol ambiguity on query sequence */
static Uint maxqryamb;

/* defined maximal allowed symbol ambiguity on subject sequence */
static Uint maxseqamb;

/* indicates whether iupac matching is enabled */
static BOOL iupac;

/* define iupac on ascii sized tabular */
#define IUPACTABSIZE 255

/* define maximal iupac bit */
#define IUPACMAXBIT 8

/*
 * iupac symbol as bit vectors, non-iupac as zeros,
 * mapping of different symbols by AND
 */
static Uint IUPACTAB[IUPACTABSIZE];

/*
 * ambiguity of iupac symbol as number of ones
 * in bitvector
 */
static Uint IUPACAMB[IUPACTABSIZE];

/*--------------------------------- getAmb -------------------------------------
 *
 * @brief       get degree of ambiguity of any char with given Uint value
 * @author      Christian Otto
 *
 */
Uint getAmb(Uint num){
  Uint count = 0;
  while(num != 0){
    count = (num & 1)?count + 1:count;
    num >>= 1;
  }
  return count;
}

/*-------------------------------- initIUPAC -----------------------------------
 *
 * @brief       initializes IUPAC table using one bit for each nucleotide,
 *              note that IUPACTAB is a constant since it does not depend on
 *              any given parameter (only qryamb and seqamb are parameters)
 * @author      Christian Otto
 *
 */
void initIUPAC(Uint qryamb, Uint seqamb){
  Uint i, A, C, G, T;
  memset(IUPACTAB, 0, IUPACTABSIZE * sizeof(Uint));

  maxqryamb = qryamb;
  maxseqamb = seqamb;
  iupac = (maxseqamb > 1 || maxqryamb > 1);
  
  /* define nucleotides */
  A = (1 << 0);
  C = (1 << 1);
  G = (1 << 2);
  T = (1 << 3);

  /* init nucleotides */
  IUPACTAB[(Uint)'A'] = A;
  IUPACTAB[(Uint)'C'] = C;
  IUPACTAB[(Uint)'G'] = G;
  IUPACTAB[(Uint)'T'] = T;
  IUPACTAB[(Uint)'U'] = T;

  /* define symbols of ambiguity 2 */
  IUPACTAB[(Uint)'R'] = (A | G);
  IUPACTAB[(Uint)'Y'] = (C | T);
  IUPACTAB[(Uint)'S'] = (G | C);
  IUPACTAB[(Uint)'W'] = (A | T);
  IUPACTAB[(Uint)'K'] = (G | T);
  IUPACTAB[(Uint)'M'] = (A | C);

  /* define symbols of ambiguity 3 */
  IUPACTAB[(Uint)'B'] = (C | G | T);
  IUPACTAB[(Uint)'D'] = (A | G | T);
  IUPACTAB[(Uint)'H'] = (A | C | T);
  IUPACTAB[(Uint)'V'] = (A | C | G);

  /* define symbol of ambiguity 4 */
  IUPACTAB[(Uint)'N'] = (A | C | G | T);

  /* define lower case chars */
  IUPACTAB[(Uint)'a'] = IUPACTAB[(Uint)'A'] << 4;
  IUPACTAB[(Uint)'c'] = IUPACTAB[(Uint)'C'] << 4;
  IUPACTAB[(Uint)'g'] = IUPACTAB[(Uint)'G'] << 4;
  IUPACTAB[(Uint)'t'] = IUPACTAB[(Uint)'T'] << 4;
  IUPACTAB[(Uint)'u'] = IUPACTAB[(Uint)'U'] << 4;
  IUPACTAB[(Uint)'r'] = IUPACTAB[(Uint)'R'] << 4;
  IUPACTAB[(Uint)'y'] = IUPACTAB[(Uint)'Y'] << 4;
  IUPACTAB[(Uint)'s'] = IUPACTAB[(Uint)'S'] << 4;
  IUPACTAB[(Uint)'w'] = IUPACTAB[(Uint)'W'] << 4;
  IUPACTAB[(Uint)'k'] = IUPACTAB[(Uint)'K'] << 4;
  IUPACTAB[(Uint)'m'] = IUPACTAB[(Uint)'M'] << 4;
  IUPACTAB[(Uint)'b'] = IUPACTAB[(Uint)'B'] << 4;
  IUPACTAB[(Uint)'d'] = IUPACTAB[(Uint)'D'] << 4;
  IUPACTAB[(Uint)'h'] = IUPACTAB[(Uint)'H'] << 4;
  IUPACTAB[(Uint)'v'] = IUPACTAB[(Uint)'V'] << 4;
  IUPACTAB[(Uint)'n'] = IUPACTAB[(Uint)'N'] << 4;

  for (i = 0; i < IUPACTABSIZE; i++){
    IUPACAMB[i] = getAmb(IUPACTAB[i]);
  }
  //DBG("qryamb:%u, seqamb:%u, isallowedIUPAC:%u\n", qryamb, seqamb, isallowedIUPAC());
}

BOOL couldMatchIUPAC(char qrych){
  if (maxseqamb == 1 &&
      (IUPACAMB[(Uint) qrych] == 1 ||
       IUPACAMB[(Uint) qrych] > maxqryamb)){
    return 0;
  }
  else {
    return 1;
  }
}

Uint countAmbChars(char *seq, Uint len){
  Uint i, cur, amb=0;
  for (i = 0; i < len; i++){
    cur = IUPACAMB[(Uint) seq[i]];
    if (cur > 1 && cur <= maxqryamb){
      amb++;
    }
  }
  return amb;
}

Uint countNonMatchingChars(char *seq, Uint len){
  Uint i, cur, cnt=len;
  for (i = 0; i < len; i++){
    cur = IUPACAMB[(Uint)seq[i]];
    if (cur > 0 && cur <= maxqryamb){
      cnt--;
    }
  }
  return cnt;
}

/*-------------------------------- matchIUPAC ----------------------------------
 *
 * @brief       indicates whether a query character matches the subject sequence
 *              character under initialized maximal ambiguity parameters
 * @author      Christian Otto
 *
 */
BOOL matchIUPAC(char qrych, char seqch){
  if (IUPACAMB[(Uint) seqch] <= maxseqamb && IUPACAMB[(Uint) qrych] <= maxqryamb){
    return ((IUPACTAB[(Uint) seqch] & IUPACTAB[(Uint) qrych]) > 0);
  }
  return 0; 
}

/*------------------------------ isallowedIUPAC --------------------------------
 *
 * @brief       check whether any ambigious IUPAC symbol is allowed in matching
 * @author      Christian Otto
 *
 */
BOOL isallowedIUPAC(){
  return iupac;
}

/*--------------------------- iupacshannonentropy ------------------------------
 *
 * @brief       minimal zero order sequence entropy for strings containing
 *              symbols of the IUPAC nucleotide code by maximizing nucleotide
 *              counts using ambigious characters
 * @author      Christian Otto
 *
 */
double minshannonentropy(char *seq, Uint len) {
  Uint i, j, k, *bitcnt, *chcnt, max, sum, isamb=0;
  double *p, H=0;
 
  /* init nucleotide counts (currently lower case differs from upper case!!!) */
  p = malloc(IUPACMAXBIT * sizeof(double));
  memset(p, 0, sizeof(double)*IUPACMAXBIT);

  /* set IUPAC symbol counts (and set isamb if ambigious symbols occur) */
  chcnt = malloc(IUPACTABSIZE * sizeof(Uint));
  memset(chcnt, 0, IUPACTABSIZE * sizeof(Uint));
  for (i = 0; i < len; i++){
    chcnt[(Uint)seq[i]]++;
    if (!isamb && IUPACAMB[(Uint)seq[i]] > 1){
      isamb = 1;
    }
  }
 
  /* 
   * maximize counts for nucleotides if
   * ambigious symbols are occuring -> minimize entropy
   */
  if (isamb){
    /* init bit counts at positions */
    bitcnt = malloc(IUPACMAXBIT * sizeof(Uint));
    
    for (i = 0; i < IUPACMAXBIT; i++){
      /* count bits at each position */
      memset(bitcnt, 0, IUPACMAXBIT * sizeof(Uint));
      for (j = 0; j < IUPACTABSIZE; j++){
	if (IUPACTAB[j] > 0 && chcnt[j] > 0){
	  for (k = 0; k < IUPACMAXBIT; k++){
	    bitcnt[k] += chcnt[j] * (1 & (IUPACTAB[j] >> k));
	  }
	}
      }
      /* get max and sum */
      max = 0; sum = 0;
      for (k = 0; k < IUPACMAXBIT; k++){
	if (bitcnt[k] > bitcnt[max]){
	  max = k;
	}
	sum += bitcnt[k];
      }
      if (sum == 0) break;
      
      /* set symbol count to zero if max-th bit is set */
      for (j = 0; j < IUPACTABSIZE; j++){
	if (1 & (IUPACTAB[j] >> max)){
	  p[max] += chcnt[j];
	  chcnt[j] = 0;
	}
      }
    }
    /* abort if sum > 0 */
    assert(sum == 0);
    free(bitcnt);
  }
  /*
   * otherwise simply count characters
   */
  else {
    sum = 0;
    k = 0;
    for (j = 0; j < IUPACTABSIZE; j++){
      if (chcnt[j] > 0){
	p[k++] = chcnt[j];
	sum += chcnt[j];
      }
      assert(k < IUPACMAXBIT);
    }
    assert(sum == len);
  }
  
  /* normalization and calculation of entropy */
  for (i = 0; i < IUPACMAXBIT; i++){
    if (p[i] > 0){
      //DBG("%u\t%g\t%u\n", i, p[i], len);
      H += (p[i]/len) * log2(p[i]/len);
    }
  }
  
  /* cleanup */
  free(chcnt);
  free(p);
  return -1 * H;
}
