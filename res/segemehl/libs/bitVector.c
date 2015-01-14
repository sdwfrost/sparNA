
/*
 *  bitVector.c
 *  implementations
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 07/14/2007 04:15:14 PM CEST
 *  
 *  SVN
 *  Revision of last commit: $Rev: 93 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-12-07 16:58:47 +0100 (Sun, 07 Dec 2008) $
 *
 *  Id: $Id: bitArray.c 19 2008-05-14 13:43:29Z steve $
 *  Url: $URL: file:///homes/bierdepot/steve/svn/segemehl/trunk/libs/bitArray.c $
 *  
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "basic-types.h"
#include "memory.h"
#include "bitVector.h"

inline bitvector
initbitvector(void *space, Uint len) {
  Uint n;

  bitvector a;

  n = len/BITVECTOR_WORDSIZE;
  n += (len % (BITVECTOR_WORDSIZE) > 0) ? 1 : 0;
  a = calloc(n, sizeof(bitvector_t));

  return a;
}

bitvector
resizebitvector(void *space, bitvector a, Uint len) {
  Uint n;

  n = len/BITVECTOR_WORDSIZE;
  n += (len % BITVECTOR_WORDSIZE > 0) ? 1 : 0;
  a = ALLOCMEMORY(space, a, bitvector_t, n);

  return a;
}

inline void
setbitvector(bitvector a, Uint len, unsigned char val) {
  Uint n;

  n = len/BITVECTOR_WORDSIZE;
  n += (len % BITVECTOR_WORDSIZE > 0) ? 1 : 0;
  memset(a,((val) ? 255 :0), n*(sizeof(bitvector_t)));
}

unsigned char
valbitvector(bitvector a, Uint len, unsigned char val) {
    Uint i;
    bitvector array;

    array = a;

    for(i=0; i < (len/BITVECTOR_WORDSIZE); i++) {  
        if (array[i] != (int) 255) 
          return 0;
    }
    for(i=0; i < (len%BITVECTOR_WORDSIZE); i++){
        if (bitvector_getbit(a, len-i-1)!= val) 
          return 0;
    }
    
    return 1;
}

void
dumpbitvector(bitvector a, Uint len) {
  Uint i;

  for(i=0; i < len; i++) {
    printf("%d ", bitvector_getbit(a, i));    
  }
  printf("\n");
}

void
bitvectorAND(bitvector dest, bitvector a, bitvector b, Uint len) {
    int i;
    int n = len/BITVECTOR_WORDSIZE;

    for(i=0; i <n ; i++) {
        dest[i] = a[i] & b[i];
    }
}

void
bitvectorOR(bitvector dest, bitvector a, bitvector b, Uint len) {
    int i;
    int n;

    n = len/BITVECTOR_WORDSIZE;

    for(i=0; i < n; i++) {
        dest[i] = a[i] | b[i];
    }
}

void
bitvectorNOT(bitvector dest, bitvector a, Uint len) {
    int i;
    int n = len/BITVECTOR_WORDSIZE;

    for(i=0; i < n; i++) {
        dest[i] = ~a[i];
    }
}

void
bitvectorXOR(bitvector dest, bitvector a, bitvector b, Uint len) {
    int i;
    int n = len/BITVECTOR_WORDSIZE;

    for(i=0; i < n; i++) {
        dest[i] = a[i] ^ b[i];
    }
}

unsigned char
bitvectorADD(bitvector dest, bitvector a, bitvector b, Uint len){
    int i,n = len/BITVECTOR_WORDSIZE;
    unsigned char carry = 0;

    for(i=0; i < n; i++) {
        dest[i] = a[i] + b[i] + carry;
        if(carry) {
          carry = ((dest[i] <= a[i]) || (dest[i] <= b[i]));
        } else { 
          carry = ((dest[i] < a[i]) || (dest[i] < b[i]));
        }
    }
    return carry;
}


void
bitvectorLSHIFT(bitvector dest, bitvector a, Uint len, Uint shift) {
    int i;

    int wordshift = shift/BITVECTOR_WORDSIZE;
    int offset = shift % BITVECTOR_WORDSIZE;
    int n = len/BITVECTOR_WORDSIZE;
    int suboffset = BITVECTOR_WORDSIZE - offset;

    if (offset == 0) {
        for(i=n-1; i >= wordshift; --i) {
            dest[i] = a[i-wordshift];
        }
    } else {
        for(i=n-1; i > wordshift; --i) {
            dest[i] = (a[i-wordshift] << offset) | (a[i-wordshift-1] >> suboffset);
        }
        dest[wordshift] = a[0] << offset; 
    }
}

void
bitvectorRSHIFT(bitvector dest, bitvector a, Uint len, Uint shift) {
    int i;

    int wordshift = shift/BITVECTOR_WORDSIZE;
    int offset = shift % BITVECTOR_WORDSIZE;
    int n = len/BITVECTOR_WORDSIZE;
    int limit = n - wordshift -1;
    int suboffset = BITVECTOR_WORDSIZE - offset;

    if (offset == 0) {
        for(i=0; i <= limit; ++i) {
            dest[i] = a[i+wordshift];
        }
    } else {
        for(i=0; i < limit; ++i) {
            dest[i] = (a[i+wordshift] >> offset) | (a[i+wordshift+1] << suboffset);
        }
        dest[limit] = a[n-1] >> offset; 
    }
}



