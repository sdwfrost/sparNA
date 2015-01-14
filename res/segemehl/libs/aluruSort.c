
/*
 *  aluruSort.c
 *  implementation of
 *  space efficient linear time construction
 *  of suffix arrays
 *
 *  Ko, Pang and Aluru, Srinivas
 *  Iowa State University
 *
 *  kind support provided by RSR:
 *  Pravda et Tanger-Glasgow sur couleur3.ch
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 07/11/2007 05:15:28 PM CEST
 *  
 *  SVN
 *  Revision of last commit: $Rev: 55 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-09-11 13:39:04 +0200 (Thu, 11 Sep 2008) $
 *
 *  Id: $Id: aluruSort.c 55 2008-09-11 11:39:04Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/trunk/libs/aluruSort.c $
 *  
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "basic-types.h"
#include "memory.h"
#include "mathematics.h"
#include "aluruSort.h"
#include "sort.h"
#include "bitArray.h"
#include "info.h"
#include "debug.h"

#define INTSIZE 64


inline void getinterval(Uint *s, Uint len, Uint *min, Uint *max) {
  Uint i;
  Lint Max,
       Min,
       resc;

  Max = (Lint)s[0]; 
  Min = (Lint)s[0];


  for(i=1; i < len; i++) {
    resc = Max - (Lint)s[i];
    resc = resc >> (INTSIZE-1);
    Max = Max + (((Lint)s[i] -Max) & resc);

    resc = s[i] - Min;
    resc = resc >> (INTSIZE-1);
    Min = Min - ((Min -(Lint)s[i]) & resc);
  }

  
  *max = Max;
  *min = Min;

}

  Uint*
getAluruArray(void *space, char *s, Uint len, char delim) 
{
  Uint i,r,k=0;
  Uint *ptr;

  ptr = ALLOCMEMORY(space, NULL, Uint*, len);

  for(i=0; i < len; i++) {
    if((r=s[i]) != delim){
      ptr[k]= i;	  
      k++;
    }
  }	
  return ptr;
}




/*-------------------------------- distCount ---------------------------------
 *    
 * @brief helper function to accumulate Qdistances
 * @author Steve Hoffmann 
 *   
 */

  Uint*
distCount (void *space, Uint *Qdist, Uint len, Uint maxdist)
{
  Uint i,
       j=0,
       *distAcc,
       temp,
       prev;

  distAcc = ALLOCMEMORY(space, NULL, Uint, maxdist+1);
  memset(distAcc, 0, sizeof(Uint)*(maxdist+1));
  
  while(Qdist[j] == 0) j++;

  for(i=j; i < len; i++) {
    distAcc[Qdist[i]-1]++;
  }

  prev = distAcc[0];
  distAcc[0] = 0;

  for(i=1; i <= maxdist; i++) {
    temp = distAcc[i];
    distAcc[i] = prev + distAcc[i-1];
    prev = temp;
  }

  return distAcc;
}


/*---------------------------------- Qdist -----------------------------------
 *    
 * @brief getting the Qdistance of each sequence member. 
 * @author Steve Hoffmann 
 *   
 */

Uint* 
Qdist(void *space, bitarray cl, Uint len, unsigned char Q) {
  Uint *sdist,
       i;
  Lint  l=-1;
  

  sdist = ALLOCMEMORY(space, NULL, Uint, len);
  for (i=0; i < len; i++) {
    sdist[i] = (l < 0) ? 0 : i-l;
    if (getbit(cl, (Uint) i) == Q)
      l = i;
  }
  
  return sdist;
}




/*--------------------------------- Qmaxdist ---------------------------------
 *    
 * @brief getting the max(d(s,l)) \forall s,l \in S
 * @author Steve Hoffmann 
 *   
 */

Uint
Qmaxdist(void *space, bitarray cl, Uint len, unsigned char Q) {
  Lint i;
  Lint  tmp=0,
       dist=0,
       pre=1;

  for (i=len-2; i >= 0; i--) {
    tmp = dist - pre;
    tmp = tmp >> (INTSIZE-1);
    dist += ((pre-dist) & tmp);

    tmp = (Q) ? (0 - (Lint)(!getbit(cl, (Uint) i))) : (((Lint)(!getbit(cl, (Uint) i))) - 1);  
    pre = (pre & tmp) + 1;
  }

  return (Uint)dist;
}


bitarray 
classifyint(void *space, Uint* s, Uint len, Uint *noL, Uint *noS) {
  Uint i, j=0, k;
  bitarray a; 

  a = initbitarray(space, len);
  setbitarray(a, len, 0);

  *noL = 0;
  *noS = 0;

  for (i=0; i < len-1; i++) {
    if (s[i] > s[i+1]){
      for(k=i-j; k <= i; k++) {
      setbit(a,k,1);
      *noL += 1;
      }
      j=0;
    }
    else if (s[i] < s[i+1]){
      
      for(k=i-j; k <= i; k++) {
      setbit(a,k,0);
      *noS+=1;
      }
      j=0;
    }
    else {
      j++;
    }
  } 
  
  if(*noS < *noL) {
    setbit(a, len-1, 0);
    *noS+=1;
  }
  else
  {
    setbit(a, len-1, 1);
    *noL+=1;
  }

  return a;
}



bitarray 
classify(void *space, char* s, Uint len, Uint *noL, Uint *noS) {
  Uint i, j=0, k;
  bitarray a; 

  a = initbitarray(space, len);
  setbitarray(a, len, 0);

  *noL = 0;
  *noS = 0;

  for (i=0; i < len-1; i++) {
    if (s[i] > s[i+1]){
      for(k=i-j; k <= i; k++) {
      setbit(a,k,1);
      *noL += 1;
      }
      j=0;
    }
    else if (s[i] < s[i+1]){
      
      for(k=i-j; k <= i; k++) {
      setbit(a,k,0);
      *noS+=1;
      }
      j=0;
    }
    else {
      j++;
    }
  } 
  
  if(*noS < *noL) {
    setbit(a, len-1, 0);
    *noS+=1;
  }
  else
  {
    setbit(a, len-1, 1);
    *noL+=1;
  }

  return a;
}


  Uint*
countingsort(void *space, 
    char *s, 
    Uint len,  
    bitarray bckts,
    Uint alphasize,
    Uint alphaoffset) 
{
  Uint  i, 
  resc,
  offset,
  *buffer,
  *A;

  buffer = ALLOCMEMORY(space, NULL, Uint, len);
  A = ALLOCMEMORY(space, NULL, Uint , len);

  /*use buffer to count chars first*/
  memset(buffer, 0, sizeof(Uint)*alphasize); 
  for(i=0; i < len; i++) {
    resc = (Uint) s[i] - alphaoffset;
    buffer[resc]++;
  }

  offset = buffer[0];
  buffer[0]=0;

  for(i=0; i < alphasize; i++) {
    resc = buffer[i];
    buffer[i] = offset + buffer[i-1];
    offset = resc;
  }

  for(i=0; i < len; i++) {
    resc = (Uint) s[i] - alphaoffset;
    A[buffer[resc]] = i;
    buffer[resc]++;
  }
  /*the bucket borders*/
  setbitarray(bckts, len, 0);
  for(i=0; i < alphasize; i++) {
    setbit(bckts, buffer[i]-1, 1);
  }

  FREEMEMORY(space, buffer);
  return A;
}


  Uint*
getlistsL(void *space, 
    Uint *A, 
    Uint len, 
    Uint *dist, 
    Uint *accdist, 
    Uint maxdist, 
    bitarray bckts,
    bitarray list, 
    Uint listlen) 
{         
  Lint i=len-1,
  j,
  pos,
  tmp,
  start,
  end;
  unsigned char firstelem;

  NFO("getlistsL: memsetting list of %u elements\n", listlen);
  setbitarray(list, listlen, 0);
  
  NFO("getlistsL: iter from %lld down to 0\n", i);
  while(i >= 0) {
    end = i;
    if(i > 0) {
      if(!getbit(bckts, (Uint) i-1))
        firstelem = 0;
      else
        firstelem = 1;
    } else 
      firstelem =1;

    while (!firstelem) {

      tmp = dist[A[i]];
/*      printf("tmp %d\n", tmp);*/
      if (tmp > 0) {
        pos = accdist[tmp-1];
        dist[A[i]]=pos;
        setbit(list, pos, 1);
        accdist[tmp-1] += 1;
      } else {
        dist[A[i]]= -1;
      }

      i--;

      if(i > 0) {
        if(!getbit(bckts, (Uint) i-1))
          firstelem = 0;
        else
          firstelem = 1;
      } else {
        firstelem = 1;
      }
    }

    tmp = dist[A[i]];
    if(tmp != 0) {
      pos = accdist[tmp-1];
      dist[A[i]] = pos;
      setbit(list, pos, 1);
      accdist[tmp-1] += 1;

    } else {
      dist[A[i]] = (Uint)-1;
    }

    start=i;
    for(j=end; j >= start; j--) {
      pos = dist[A[j]];
      if (pos != (Uint)-1 && pos != listlen -1) {
        if(getbit(list, (Uint) pos+1))
          setbit(list, (Uint) pos, 0);
      }
    }
    i--;
  }

  NFO("scanning A (%u elems)\n", len);
  for(i=0; i < len; i++) {
    if(dist[i] != (Uint)-1) {
      A[dist[i]] = i;    
    }
  }

  NFO("scanning accdist (%u elems) (1)\n", maxdist);
  for(i=0; i < maxdist; i++) {
    if(accdist[i] > 0) {
        setbit(list, accdist[i]-1, 1);
    } else {
        setbit(list, accdist[i], 1);
    }
  }

  NFO("scanning accdist (%u elems) (2)\n", maxdist);
  for(i=0; i < maxdist; i++) {
    if (i==0) {
      j=0;
    } else {
      j = accdist[i-1];
    }
    while(j < accdist[i]) {
      A[j] = A[j] -i -1;
      j++;
    }
  }

  MSG("getlistsL: exit\n");
  return A;
}


  Uint*
getlistsS(void *space, 
    Uint *A, 
    Uint len, 
    Uint *dist, 
    Uint *accdist, 
    Uint maxdist, 
    bitarray bckts,
    bitarray list, 
    Uint listlen) 
{
 Uint i=0,
  j,
  pos,
  tmp,
  start,
  end;

  NFO("getlistsS: memsetting list of %u elements\n", listlen);
  setbitarray(list, listlen, 0);

  NFO("getlistsS: iter up to %u\n", len);
  while(i < len) {
    start = i;
    
    while(getbit(bckts, (Uint) i) != 1 && i < len) {
      tmp = dist[A[i]];
      if(tmp != (Uint)-1 && tmp > 0) {
        pos = accdist[tmp-1];
        dist[A[i]] = pos;
        setbit(list, pos, 1);
        accdist[tmp-1] += 1;
      }  else {
        dist[A[i]] = (Uint)-1;
      }
      i++;
    }

    tmp = dist[A[i]];
    if (tmp != 0) {
      pos = accdist[tmp-1];
      dist[A[i]] = pos;
      setbit(list, pos, 1);
      accdist[tmp-1] += 1;
    } else {
      dist[A[i]] = (Uint) -1;
    }
    end = i;

    for(j= start; j < end; j++) {
      pos = dist[A[j]];
      if (pos != (Uint)-1 && pos != listlen -1) {
        if (getbit(list, (Uint) pos+1)) {
          setbit(list,pos, 0);
        }
      }
    }
    i++;
  }

  MSG("getlistsS: scan A\n");
  for(i=0; i < len; i++) {
    if(dist[i] != (Uint)-1)
      A[dist[i]]=i;
  }

  MSG("getlistsS: set accidst\n");
  for(i=0; i < maxdist; i++) {
    if (accdist[i] == 0) NFO("getlistsS: i=%u accdist=0!!\n", i);
    else
    setbit(list, accdist[i]-1, 1);
  }

  for(i=0; i < maxdist; i++) {
    j = (i==0) ? 0 : accdist[i-1];
    while(j < accdist[i]) {
      A[j] = A[j] - i -1;
      j++;
    }
  }

  MSG("getlistsS: exiting\n");
  return A;
}

  void
sortlistS(void *space,
    Uint *B,
    Uint lenB,
    Uint len,
    bitarray bckts,
    Uint *list,
    bitarray listb,
    Uint listlen)
{
  Lint *rev,
      *left,
       bcktno,
        i,
        j,
        new, 
        bcktright;

  MSG("sortlistS: allocating stuff\n");
  rev = ALLOCMEMORY(space,  NULL,  Lint, len);
  left = ALLOCMEMORY(space, NULL, Lint, lenB);

  memset(rev, -1, len*sizeof(Lint));
  memset(left, -1, lenB*sizeof(Lint));

  bcktright= lenB -1;

  NFO("sortlistS: iterating %u elems", lenB);
  for(i=lenB-1; i > 0; i--) {
    rev[B[i]] = bcktright;
    if(getbit(bckts, (Uint) i-1) == 1) {
      left[bcktright] = i;
      bcktright = i-1;
    }
  }

  rev[B[0]] = bcktright;
  left[bcktright] = 0;

  NFO("sortlistS: looping %u elems", listlen);
  i=0;
  while (i < listlen) {
    j=i;
    while(!getbit(listb, (Uint)j)) {
      left[rev[list[j]]] += 1;
      j++;
    }

    left[rev[list[j]]] += 1;

    j=i;
    while(!getbit(listb, (Uint) j)) {
      new = left[rev[list[j]]] - 1;
      rev[list[j]] = new;
      j++;
    }

    new = left[rev[list[j]]] - 1;
    rev[list[j]] = new;

    /*correct the values*/

    j=i;
    while (!getbit(listb, (Uint) j)) {
      new = rev[list[j]];
      if (left[new] == -1) {
        left[new] = new;
      } else {
        left[new] -= 1;
      } 
      setbit(bckts, new, 1);
      j++;
    }

    /*last elem*/
    new = rev[list[j]];
    if(left[new] == -1) {
      left[new] = new;
    } else {
      left[new] -=1;
    }
    
    setbit(bckts, new, 1);
    i=j+1;
  }


  NFO("sortlistS: iterating %u elems", len);
  for(i=0; i < len; i++) {
    bcktno = rev[i];
    if(bcktno > -1) {
        B[left[bcktno]] = i;
        left[bcktno] += 1;
    }
  }

  MSG("sortlistsS: exiting happily!\n");
  FREEMEMORY(space, rev);
  FREEMEMORY(space, left);
}

  void
sortlistL(void *space,
    Uint *B,
    Uint lenB,
    Uint len,
    bitarray bckts,
    Uint *list,
    bitarray listb,
    Uint listlen)
{
  Lint *rev,
      *right,
       bcktno,
        i,
        j,
        new, 
        bcktleft;

  MSG("sortlistL: allocating stuff\n");
  rev =   ALLOCMEMORY(space, NULL, Lint, len);
  right = ALLOCMEMORY(space, NULL, Lint, lenB);

  memset(rev,   -1, len*sizeof(Lint));
  memset(right, -1, lenB*sizeof(Lint));

  NFO("sortlistL: iterating %u elems", lenB);
  bcktleft=0;
  for(i=0; i < lenB; i++) {
    
    rev[B[i]] = bcktleft;
    if(getbit(bckts, (Uint) i) == 1) {
      right[bcktleft] = i;
      bcktleft = i+1;
    }
  }

  NFO("sortlistL: looping %u elems", listlen);
  i=0;
  while (i < listlen) {
    
    j=i;
    while(getbit(listb, (Uint) j) == 0) {
      right[rev[list[j]]] -= 1;
      j++;
    }

    right[rev[list[j]]] -= 1;

    j=i;
    while(getbit(listb, (Uint) j) == 0) {
      new = right[rev[list[j]]] + 1;
      rev[list[j]] = new;
      j++;
    }

    new = right[rev[list[j]]] + 1;
    rev[list[j]] = new;

    /*correct elems*/

    j=i;
    while (getbit(listb, (Uint) j) == 0) {
      new = rev[list[j]];
      if (right[new] == -1) {
        right[new] = new;
      } else {
        right[new] += 1;
      }
      if(new > 0) {
        setbit(bckts, new-1, 1);
      }
      j++;
    }

    /*correct last elem*/

    new = rev[list[j]];
    if(right[new] == -1) {
      right[new] = new;
    } else {
      right[new] +=1;
    }

    if(new > 0) {
      setbit(bckts, new-1, 1);
    }

    i=j+1;
  }

  NFO("sortlistL: iterating %u elems", len);
  for(i=0; i < len; i++) {
    bcktno = rev[i];
    if(bcktno > -1) {
        B[right[bcktno]] = i;
        right[bcktno] -= 1;
    }
  }

  MSG("sortlistsL: exiting happily!\n");
  FREEMEMORY(space, rev);
  FREEMEMORY(space, right);
}


  Uint*
countingsortint(void *space, 
    Uint *s, 
    Uint len,  
    bitarray bckts) 
{
  Uint  i, 
  resc,
  offset,
  *buffer,
  *A,
  min,
  max,
  sigma;

  getinterval(s, len, &min, &max);
  sigma = max - min +1;

  MSG("countingsortint: init buffers and A\n");
  buffer = ALLOCMEMORY(space, NULL, Uint, len);
  A = ALLOCMEMORY(space, NULL, Uint , len);

  MSG( "setting buffer to zero\n");
  /*use buffer to count chars first*/
  memset(buffer, 0, sizeof(Uint)*sigma); 
  
  MSG("countsortint: scanning buffer (1 of 3)\n");
  for(i=0; i < len; i++) {
    resc = (Uint) (s[i] - min);
    buffer[resc]++;
  }

  offset = buffer[0];
  buffer[0]=0;

  MSG("countsortint: scanning buffer (2 of 3)\n");
  for(i=1; i < sigma; i++) {
    resc = buffer[i];
    buffer[i] = offset + buffer[i-1];
    offset = resc;
  }

  MSG("countsortint: scanning buffer (3 of 3)\n");
  for(i=0; i < len; i++) {
    resc = (Uint) (s[i] - min);
    A[buffer[resc]] = i;
    buffer[resc]++;
  }



  MSG("countsortint: scanning buffer (to set borders)\n");
  /*the bucket borders*/
  setbitarray(bckts, len, 0);
  for(i=0; i < sigma; i++) {
    setbit(bckts, buffer[i]-1, 1);
  }

 
  MSG("countsortint: exiting\n");
  FREEMEMORY(space, buffer);
  return A;
}


  Uint*
substringsort(void *space,
    char *s,
    Uint *A,
    bitarray cl,
    Uint len, 
    bitarray bckts,
    Uint bucketno,
    Uint dist,
    Uint Q) 
{
  Uint bufferlen = 255 * 2,
  i,
  j=0,
  idx,
  offset;
  Lint  start,
       end,
       resc,
       prevCount,
       tempBucketTest=0,
       *tmp,
       *skip,
       *buffer;
  Lint  type = (Q) ? 0 : 1;

  MSG("setting bit array to zero\n");
  setbitarray(bckts, bucketno, 0);

  MSG("allocating space for buckets and buffers\n");
  buffer = ALLOCMEMORY(space, NULL, Lint, bufferlen+1);
  skip = ALLOCMEMORY(space, NULL, Lint, bucketno+1);
  tmp = ALLOCMEMORY(space, NULL, Lint, bucketno+1);

  MSG("memsetting\n");

  memset(skip, 0, sizeof(Lint)*(bucketno+1));
  skip[0] = bucketno;

   for(i=0; i <= dist; i++) { /*offset*/
    start  =0;
    offset =0; /*prevPos*/
    while(start < bucketno) {
      offset = start;
      while (skip[start] < 0 && start < bucketno) {
        start = (Lint) -skip[start];
      }
      end = skip[start] - 1;
      skip[offset] = -start;
      
      memset(buffer, 0, sizeof(Lint)*bufferlen);

      if(start < bucketno) {
        
        for(j=start; j <=end; j++) { /*i*/
          tempBucketTest++;        
          tmp[j] = A[j];

          idx = A[j] + i;
          resc =  ((Lint) s[idx]) << 1;
          resc += (!getbit(cl, (Uint) idx)) ? 1 : 0;
          buffer[resc] += 1;
        }

        prevCount = buffer[0];
        buffer[0] = start;

        for(j=1; j < bufferlen; j++) {
          resc = buffer[j];
          buffer[j] = buffer[j-1]+prevCount;
          prevCount = resc;
        }

        for(j=start; j <= end; j++) {
          tempBucketTest++;
          idx = tmp[j] + i;
          resc = ((Lint) s[idx]) << 1;
          resc += (!getbit(cl, (Uint)idx)) ? 1 : 0;
          A[buffer[resc]] = tmp[j];
          buffer[resc]++;
        }

        /*bucket boundaries*/ 

        j=1;
        if(i > 0) {

          if(buffer[type] > start) {
            setbit(bckts, buffer[type]-1, 1);
            skip[start] = -buffer[0];
          }

          for(j=1; j < bufferlen;  j++) {
            
            if(buffer[j] == buffer[j-1]+1) {
              setbit(bckts, buffer[j]-1, 1);
              skip[buffer[j-1]] = -buffer[j];
            
            } else 
              if (buffer[j] > buffer[j-1] + 1) {
                setbit(bckts, buffer[j]-1, 1);
                resc = (Q) ? -((j & 1)^1): -(j & 1);
                resc = (buffer[j] ^ resc) - resc;
                skip[buffer[j-1]] = resc;
              }
          }  
        } 
        /*if first bucket greater start not empty*/ 
        else {
          if(buffer[type] > start) {
            setbit(bckts, buffer[type]-1, 1);
            skip[start] = -buffer[0];
          }    

          for(j=1; j < bufferlen; j++) {
            if(buffer[j] == buffer[j-1] +1) {
              setbit(bckts, buffer[j]-1, 1);
              skip[buffer[j-1]] = -buffer[j];
            } else if(buffer[j] > buffer[j-1] +1) {
              setbit(bckts, buffer[j]-1, 1);
              skip[buffer[j-1]] = buffer[j];
            }
          }
        }
        if(type && start == end) {
          skip[start] = -(end - 1);
          setbit(bckts, start, 1);
        }

        start = end + 1;
      }
    }
  }
  
  MSG("substring sort ... ok\n"); 
  FREEMEMORY(space, skip);
  FREEMEMORY(space, tmp);
  FREEMEMORY(space, buffer);

  return A;
}


  Uint*
arrayB (void *space, 
    Uint *A, 
    Uint lenA, 
    Uint lenB,
    bitarray bcktsA, 
    bitarray bcktsB,
    bitarray cl,
    unsigned char Q) 
{
  Uint i, 
  *B;
  Lint j=0;
  unsigned char type = Q ? 1 : 0;

  NFO("arrayB: allocating B with %u elements\n", lenB);
  B = ALLOCMEMORY(space, NULL, Uint, lenB);
  memset(B, 0, sizeof(Uint)*lenB);
  setbitarray(bcktsB, lenB, 0);
  
  NFO("arrayB: iterating to lenA=%u\n", lenA);
  for(i=0; i < lenA; i++) {
    if (getbit(cl, (Uint) A[i]) == type) { 
      B[j] = A[i];
      j++;
    }
    if (j > lenB) DBG("arrayB: j=%lld in B out of bounds!\n", j);
    /*copy bckt-borders*/
    if(getbit(bcktsA, (Uint)i) ==1 && j-1 >= 0)
      setbit(bcktsB, j-1, 1);
  }

  MSG("arrayB: exiting\n");
  return B;
}


  Uint*
Tprime(void *space, 
    Uint len, 
    Uint *B, 
    Uint lenB, 
    bitarray bcktsB, 
    bitarray cl,
    unsigned char Q) 
{
  Uint i,
  j=0, 
  *tprime;
  Lint *buffer;
  Lint  cur=0,
       inv;

  MSG("tprime: init arrays\n");
  buffer = ALLOCMEMORY(space, NULL, Lint, len);
  tprime = ALLOCMEMORY(space, NULL, Uint, lenB);
  
  memset(buffer, 0, sizeof(Lint)*len);

  MSG( "tprime: scan B\n");
  for(i=0; i < lenB; i++) {
    buffer[B[i]] = cur;
    cur += getbit(bcktsB,(Uint) i) ? 1 : 0;
  }

  NFO("tprime: iterating i=%u elements with lenB=%u\n", len, lenB);
  for(i=0; i < len; i++) {
    cur = (Q) ? (((Lint)!getbit(cl,(Uint)i)) - 1) : -((Lint)!getbit(cl, (Uint)i));
    inv = ~cur;
    tprime[j] = (Uint) ((((Lint)tprime[j]) & inv) | (buffer[i] & cur));
    if (j >= lenB) DBG( "j=%u out of bounds\n", j);
    j += (1 & cur) ? 1 : 0; 
  }

  MSG( "tprime: exit\n");
  FREEMEMORY(space, buffer);
  return tprime;
}


  void
reconstruct(void *space,
    Uint len,
    Uint *B,
    Uint lenB,
    bitarray cl,
    unsigned int Q) 
{
  Lint *conv;
  Lint i,
      j=0,
      cur,
      inv;

  MSG("reconstruct: init\n");
  conv = ALLOCMEMORY(space, NULL, Lint, lenB);
  
  NFO("reconstruct: iteration over %u elems\n", len);
  for(i=0; i < len; i++) {
    cur = (getbit(cl, (Uint)i)) ? 0 : 1;
    cur = cur << (INTSIZE-1);
    cur = cur >> (INTSIZE-1);
    inv = ~cur;
   if(Q) {     
    conv[j] = ((i & inv) | (((Lint)conv[j]) & cur));
    j = j + (1 & inv);
   } else {
    conv[j] = ((i & cur) | (((Lint)conv[j]) & inv));
    j = j + (1 & cur);
   }
  }

  NFO("reconstruct: scan B (size: %u)\n", lenB);
  for(i=0; i < lenB; i++) {
    cur = B[i];
    B[i] = (Uint) conv[cur];
  }

  MSG("reconstruct: exit");
  FREEMEMORY(space, conv);
  return;
}


  Uint*
aluruSuffixArrayS(void *space,
    char* T,
    Uint len,
    Uint *B,
    Uint lenB,
    bitarray cl) 
{
  Lint *count;
  Uint *sarray;
  Lint  i,
       j; 
  Lint  tmp;
  Lint  offset;
  bitarray b;

  count = ALLOCMEMORY(space, NULL, Lint, 256);
  sarray = ALLOCMEMORY(space, NULL, Uint, len);
  b = initbitarray(space, len);

  memset(count, 0, 256*sizeof(Lint));
  setbitarray(b, len, 0);

  for(i=0; i < len; i++) {
    tmp = (Lint) T[i];
    count[tmp]++;
  }

  offset = count[0];
  count[0] = 0;

  for(i=1; i < 255; i++) {
    tmp = count[i];
    count[i] = count[i-1] + offset;
    offset = tmp;
  }

  j=0;

  for(i=0; i < len; i++) {
    if (!getbit(b, (Uint) i)) {
      sarray[i] = B[j];
      setbit(b, i, 1);
      j++;
      offset = (Lint) sarray[i] -1;
      if (offset >= 0) {
        if(getbit(cl, (Uint) offset)) {
          tmp = (Lint) T[offset];
          if (count[tmp] > i) {
            sarray[count[tmp]] = offset;
            setbit(b, count[tmp], 1);
            count[tmp] += 1;
          }
        }
      }
    } else {
      offset = (Lint) sarray[i] -1;
      if(offset >= 0) {
        if(getbit(cl, (Uint) offset)) {
          tmp = (Lint) T[offset];
          if (count[tmp] > i) {
            sarray[count[tmp]] = offset;
            setbit(b, count[tmp], 1);
            count[tmp] += 1;
          }
        }
      }
    }
  }
  FREEMEMORY(space, count);
  FREEMEMORY(space, b);
  return sarray;
}

  Uint*
aluruSuffixArrayL(void *space,
    char* T,
    Uint len,
    Uint *B,
    Uint lenB,
    bitarray cl) 
{
  Lint *count;
  Uint *sarray;
  Lint  i,
       j;
  Lint tmp;
  Lint  offset;
  bitarray b;

  MSG( "aluruSuffixArrayL: initalizning arrays\n");
  count = ALLOCMEMORY(space, NULL, Lint, 255);
  sarray = ALLOCMEMORY(space, NULL, Uint, len);
  b = initbitarray(space, len);

  MSG("aluruSuffixArrayL: memsetting count\n");
  memset(count, 0, 255 *sizeof(Lint));

  MSG("aluruSuffixArrayL: setting b\n");
  setbitarray(b, len, 0);

  for(i=0; i < len; i++) {
    tmp = (Lint) T[i];
    count[tmp]++;
  }

  count[0] = count[0] -1;

  for(i=1; i < 255; i++) {
    count[i] = count[i-1] + count[i];
  }

  j=lenB -1;

  MSG( "aluruSuffixArrayL: iteration\n");
  for(i=len-1; i >= 0; i--) {
    if (!getbit(b, (Uint) i)) {
      sarray[i] = B[j];
      setbit(b, (Uint) i, 1);
      j--;
      offset = (Lint) sarray[i] -1;
      if (offset >= 0) {
        if(!getbit(cl, (Uint) offset)) {
          tmp = (Lint) T[offset];
          if (count[tmp] < i) {
            sarray[count[tmp]] = offset;
            setbit(b, (Uint) count[tmp], 1);
            count[tmp] -= 1;
          }
        }
      }
    } else {
      offset = (Lint) sarray[i] -1;
      if(offset >= 0) {
        if(!getbit(cl, (Uint) offset)) {
          tmp = (Lint) T[offset];
          if (count[tmp] < i) {
            sarray[count[tmp]] = offset;
            setbit(b, (Uint) count[tmp], 1); 
            count[tmp] -= 1;
          }
        }
      }
    }
  }
  
  FREEMEMORY(space, count);
  FREEMEMORY(space, b);
  MSG("aluruSuffixArrayL: exit ok\n");
  return sarray;
}


  Uint*
aluruSuffixArraySint(void *space,
    Uint* T,
    Uint len,
    Uint *B,
    Uint lenB,
    bitarray cl) 
{
  Uint *count;
  Uint *sarray;
  Uint i,
       j, 
       min,
       max,
       sigma;
  Uint tmp;
  Lint  offset;
  bitarray b;

  getinterval(T, len, &min, &max);
  sigma = max - min +1;

  count = ALLOCMEMORY(space, NULL, Uint, sigma);
  sarray = ALLOCMEMORY(space, NULL, Uint, len);
  memset(count, 0, sigma*sizeof(Uint));

  b = initbitarray(space, len);
  setbitarray(b, len, 0);

  for(i=0; i < len; i++) {
    tmp = (Uint) T[i] - min;
    count[tmp]++;
  }

  offset = count[0];
  count[0] = 0;

  for(i=1; i < sigma; i++) {
    tmp = count[i];
    count[i] = count[i-1] + offset;
    offset = tmp;
  }

  j=0;
  for(i=0; i < len; i++) {
    if (!getbit(b, (Uint) i)) {
      sarray[i] = B[j];
      setbit(b, i, 1);
      j++;
      offset = (Lint) sarray[i] -1;
      if (offset >= 0) {
        if(getbit(cl, (Uint) offset)) {
          tmp = (Uint) T[offset];
          if (count[tmp] > i) {
            sarray[count[tmp]] = offset;
            setbit(b, count[tmp], 1);
            count[tmp] += 1;
          }
        }
      }
    } else {
      offset = (Lint) sarray[i] -1;
      if(offset >= 0) {
        if(getbit(cl, (Uint) offset)) {
          tmp = (Lint) T[offset];
          if (count[tmp] > i) {
            sarray[count[tmp]] = offset;
            setbit(b, count[tmp], 1);
            count[tmp] += 1;
          }
        }
      }
    }
  }
  FREEMEMORY(space, count);
  FREEMEMORY(space, b);

  return sarray;
}


  Uint*
aluruSuffixArrayLint(void *space,
    Uint* T,
    Uint len,
    Uint *B,
    Uint lenB,
    bitarray cl) 
{
  Uint *count;
  Uint *sarray;
  Uint min,
       max,
       sigma;
  Uint tmp;
  Lint  offset;
  Lint  i;
  Lint  j;
  bitarray b;
  
  getinterval(T, len, &min, &max);
  sigma = max - min +1;

  count = ALLOCMEMORY(space, NULL, Uint, sigma);
  sarray = ALLOCMEMORY(space, NULL, Uint, len+1);
  memset(sarray, 0, (len+1)*sizeof(Uint));
  memset(count, 0, sigma*sizeof(Uint));
  b = initbitarray(space, len);

  memset(count, 0, sigma*sizeof(Uint));
  setbitarray(b, len, 0);

  for(i=0; i < len; i++) {
    tmp = (Lint) T[i] - min;
    count[tmp]++;
  }

  count[0] = count[0] -1;

  for(i=1; i < sigma; i++) {
    count[i] = count[i-1] + count[i];
  }

  j=lenB -1;

  for(i=len-1; i >= 0; i--) {
    if (!getbit(b, (Uint) i)) {
      sarray[i] = B[j];
      setbit(b, i, 1);
      j--;
      offset = (Lint) sarray[i] -1;
      if (offset >= 0) {
        if(!getbit(cl, (Uint) offset)) {
          tmp = (Lint) T[offset];
          if (count[tmp] < (i)) {
            sarray[count[tmp]] = offset;
            setbit(b, count[tmp], 1);
            count[tmp] -= 1;
          }
        }
      }
    } else {
      offset = (Lint) sarray[i] -1;
      if(offset >= 0) {
        if(!getbit(cl, (Uint) offset)) {
          tmp = (Lint) T[offset];
          if (count[tmp] < i) {
            sarray[count[tmp]] = offset;
            setbit(b, count[tmp], 1);
            count[tmp] -= 1;
          }
        }
      }
    }
  }
  FREEMEMORY(space, count);
  FREEMEMORY(space, b);
  return sarray;
}


Uint*
alurusortint(void *space, Uint *s, Uint *l) {
  bitarray cl=NULL,
           bcktsA=NULL,
           bcktsB=NULL,
           bcktslist=NULL;
  Uint *B,
       *A,
       *tprime; 
  Uint noL;
  Uint noS;
  Uint maxdist,
       *dist,
       *accDist,
       *list,
       listlen;
  Uint len = *l;

  MSG("alurusortint: classify int\n");
  cl = classifyint(space, s, len, &noL, &noS);
  
  MSG("alurusortint: getting bit\n");
  if(!getbit(cl, (Uint) len-1) && noS ==1) {

    MSG("alurusortint: aluruSuffixArraySint\n");
    /*printf("fewintS\n");*/
    B = ALLOCMEMORY(space, NULL, Uint, 1);
    B[0] = len-1;
    A = aluruSuffixArraySint(space, s, len, B, noS, cl);
    FREEMEMORY(space, B);
    FREEMEMORY(space, cl);

    *l = len;
    return A;
  }

  MSG("alurusortint: init bcktsA\n");
  bcktsA = initbitarray(space, len);
  MSG("alurusortint: countingsort\n");
  A=countingsortint(space, s, len, bcktsA);

  /*sort type S suffixes*/
  if(!getbit(cl, (Uint) len-1)) {
   
    MSG("alurusortint: Sorting type S suffixes. Init bcktsB.\n");
    NFO("%d\t%d\t%d\n\n", noS, noL, len);
  
    bcktsB = initbitarray(bcktsB, noS);
    B = arrayB(space, A, len, noS, bcktsA, bcktsB, cl, 0);
    MSG("alurusortint: enter Qmaxdist\n");
    maxdist = Qmaxdist(space, cl, len, 0);

    MSG("alurusortint: enter Qdist\n");
    dist = Qdist(space, cl, len, 0);


    MSG("alurusortint: enter distCount\n");
    accDist = distCount(space, dist, len, maxdist);

    listlen = accDist[maxdist];
    bcktslist = initbitarray(space, listlen);


    MSG("alurusortint: enter get listsS\n");
    /*list points to modified A*/
    list = getlistsS(space, A, len, dist, accDist, maxdist, bcktsA, 
                        bcktslist, listlen);
   
    MSG("alurusortint: freeing stuff\n");
    FREEMEMORY(space, bcktsA);
    FREEMEMORY(space, accDist);
    FREEMEMORY(space, dist);

    MSG("alurusortint: enter sortlistsS\n");
    sortlistS(space, B, noS, len, bcktsB, list, bcktslist, listlen);
    FREEMEMORY(space, list);
    FREEMEMORY(space, bcktslist); 
    
    if(valbitarray(bcktsB, noS, 1)) {

      FREEMEMORY(space, bcktsB);
      /*  printf("valbitarraysortedS.\n");*/

      MSG( "alurusortint: valbitarraysortedS.\n");
      A = aluruSuffixArraySint(space, s, len, B, noS, cl);

      FREEMEMORY(space, B);
      FREEMEMORY(space, cl);

      *l = len;
      return A;
    }

    MSG("alurusortint: enter tprime\n");
    tprime = Tprime(space, len, B, noS, bcktsB, cl, 0);
    FREEMEMORY(space, B);
    FREEMEMORY(space, bcktsB);

    MSG("alurusortint: enter alurusortint\n");
    B = alurusortint(space, tprime, &noS);
    FREEMEMORY(space, tprime);
    
    MSG("reconstructintS\n");
    reconstruct(space, len, B, noS, cl, 0);

    A = aluruSuffixArraySint(space, s, len, B, noS, cl);
    FREEMEMORY(space, B);
    FREEMEMORY(space, cl);

    *l = len;
    return A; 

  } else {
    
    /*type L suffixes*/

    /*printf("Sorting type L suffixes\n"); 
    printf("%d\t%d\t%d\n\n", noS, noL, len);*/

    bcktsB = initbitarray(bcktsB, noL);
    B = arrayB(space, A, len, noL, bcktsA, bcktsB, cl, 1);
 
    maxdist = Qmaxdist(space, cl, len, 1);
    dist = Qdist(space, cl, len, 1);
    accDist = distCount(space, dist, len, maxdist);
    
    listlen = accDist[maxdist];
    bcktslist = initbitarray(bcktslist, listlen);

    /*list points to modified A*/

    MSG("alurusortint: enter get listsL\n");
    list = getlistsL(space, A, len, dist, accDist, maxdist, bcktsA, 
                        bcktslist, listlen);
 
    FREEMEMORY(space, bcktsA);
    FREEMEMORY(space, accDist);
    FREEMEMORY(space, dist);
 

    MSG("alurusortint: sort listsL\n");
    sortlistL(space, B, noL, len, bcktsB, list, bcktslist, listlen);

    FREEMEMORY(space, bcktslist);
    FREEMEMORY(space, list);
    
    if(valbitarray(bcktsB, noL, 1)) {
    
      FREEMEMORY(space, bcktsB);
      MSG("alurusortint: valbitarraysortedL.\n");
      A = aluruSuffixArrayLint(space, s, len, B, noL, cl);

      FREEMEMORY(space, B);
      FREEMEMORY(space, cl);
      *l = len;
      return A;
    }

    MSG("alurusortint: enter tprime\n");
    tprime = Tprime(space, len, B, noL, bcktsB, cl, 1);
    FREEMEMORY(space, bcktsB);
    FREEMEMORY(space, B);


    MSG("alurusortint: enter alurusortint\n");
    B = alurusortint(space, tprime, &noL);
    FREEMEMORY(space, tprime);
 
    MSG("reconstructintL\n");
    reconstruct(space, len, B, noL, cl, 1);
    
    A = aluruSuffixArrayLint(space, s, len, B, noL, cl);
    FREEMEMORY(space, B);
    FREEMEMORY(space, cl);
    
    *l = len;
    return A;
  }
}


Uint*
alurusort(void *space, char *s, Uint *l) {
  bitarray cl=NULL,
           bckts=NULL;
  Uint *B,
       *A,
       *tprime,
       i,
       j;
  Uint noL;
  Uint noS;
  Uint dist;
  Uint len = *l;
  
  space=NULL;

  MSG("alurusort: classify\n");
  cl = classify(space, s, len, &noL, &noS);
  
  MSG( "alurusort: getting bit\n");
  if(!getbit(cl, (Uint) len-1) && noS ==1) {
    MSG("alurusort: fewcharS\n");
    B = ALLOCMEMORY(space, NULL, Uint, 1);
    B[0] = len-1;
    A = aluruSuffixArrayS(space, s, len, B, noS, cl);
    
    FREEMEMORY(space, B);
    FREEMEMORY(space, cl);
    *l = len;
    return A;
  }

  if(!getbit(cl, (Uint) len-1)) {

    NFO("not bit alurusort: alloc B of size %u\n", noS);
    B = ALLOCMEMORY(space, NULL, Uint, noS);

    NFO("alurusort: initbitarray of size %u\n", noS);
    bckts = initbitarray(bckts, noS);

    NFO("alurusort: Qmaxdist in cl of size %u\n", len);
    dist = Qmaxdist(space, cl, len, 0);

    MSG("alurusort: scan B\n");
    
    j=0;
    for(i=0; i < len; i++) {
      B[j] = i;
      if (j > len) DBG("%u > %u\n", j, len);
      j +=  (!getbit(cl, (Uint) i)) ? 1 : 0;
    }
    
    MSG("alurusort: substringsort\n");
    B = substringsort(space, s, B, cl, len, bckts, noS, dist, 0);

    MSG("checking valbitarray\n");
    if(valbitarray(bckts, noS, 1)) {

      FREEMEMORY(space, bckts);
      /*printf("valbitarraysortedcharS.\n");*/

      MSG("aluruSuffixArrayS start (if cond 1)\n");
      A = aluruSuffixArrayS(space, s, len, B, noS, cl);

      FREEMEMORY(space, B);
      FREEMEMORY(space, cl);

      *l = len;
      return A;
    }

    MSG("enter Tprime calculation");
    tprime = Tprime(space, len, B, noS, bckts, cl, 0);
    FREEMEMORY(space, B);
    FREEMEMORY(space, bckts);
    
    MSG("enter alursortint \n");
    B = alurusortint(space, tprime, &noS);
    FREEMEMORY(space, tprime);

    MSG("reconstructcharS\n");
    reconstruct(space, len, B, noS, cl, 0);

    MSG( "enter aluruSuffixArrayS start\n");   
    A = aluruSuffixArrayS(space, s, len, B, noS, cl);
    FREEMEMORY(space, B);
    FREEMEMORY(space, cl);

    *l = len;
    return A; 

  } else {

    NFO("bit alurusort: alloc B of size %u\n", noL);
    B = ALLOCMEMORY(space, NULL, Uint, noL);
    
    NFO( "alurusort: initbitarray of size %u\n", noL);
    bckts = initbitarray(bckts, noL);
    
    NFO( "alurusort: Qmaxdist in cl of size %u\n", len);
    dist = Qmaxdist(space, cl, len, 1);

    j=0;
    for(i=0; i < len; i++) {
      B[j] = i;
      j -= ((Lint) (!getbit(cl, (Uint) i)) -1);
    }

    MSG("enter alurusort: substringsort\n");
    B = substringsort(space, s, B, cl, len, bckts, noL, dist, 1);
    MSG( "checking valbitarray\n");
    if(valbitarray(bckts, noL, 1)) {
      FREEMEMORY(space, bckts);
      /*printf("valbitarraysortedcharL.\n");*/
      MSG( "aluruSuffixArrayL start (if cond 1)\n");
      A = aluruSuffixArrayL(space, s, len, B, noL, cl);

      FREEMEMORY(space, B);
      FREEMEMORY(space, cl);

      *l = len;
      return A;
    }

    MSG("enter Tprime calculation");
    tprime = Tprime(space, len, B, noL, bckts, cl, 1);
    FREEMEMORY(space, B);
    FREEMEMORY(space, bckts);

    MSG("enter alursortint \n");
    B = alurusortint(space, tprime, &noL);
    FREEMEMORY(space, tprime);

    /*printf("reconstructcharL\n");*/

    MSG("enter reconstruction\n");
    reconstruct(space, len, B, noL, cl, 1);
    MSG("enter aluruSuffixArrayL start\n");
    A = aluruSuffixArrayL(space, s, len, B, noL, cl);

    FREEMEMORY(space, B);
    FREEMEMORY(space, cl);
    *l = len;

    return A;
  }
}

void
showQDlist(vector_t **qdlist, Uint n) {
  Uint i;

  for(i=0; i < n; i++) {
    printf("list %d\n", i);
    dumpVector(qdlist[i]);
    printf("\n");
  }
}


void
showAluruBuckets(Alurubucket *bckts, Uint *R, Uint n) {
  Uint i,j,k=0;

  for(i=0; i < n; i++) {
    printf("bucket %d\n", i);
    for(j=0; j < bckts[i].noofelems; j++) {
      printf("A[%d]=%d, R[%d]=%d", k++, bckts[i].elems[j], bckts[i].elems[j], R[bckts[i].elems[j]]);
    }
    printf("\n");
  }

}

int
bcktcmpANSI(const void *a, const void *b) {
  Alurubucket *first = (Alurubucket *)a,
              *second = (Alurubucket *)b; 

  if (first->id < second->id) return -1;
  if (first->id > second->id) return 1;

  return 0;
}


Uint
bcktcmp(Uint a, Uint b, void *arr, void *info) {
  Alurubucket *bckts;

  bckts = (Alurubucket*) arr;

  if (bckts[a].id < bckts[b].id) return 2;
  if (bckts[a].id > bckts[b].id) return 1;

  return 0;
}


  void
sortAluruSubstrings (void *space, 
    vector_t** qdlist, 
    Uint nooflists, 
    Alurubucket *bckts, 
    Uint noofbckts,
    Uint *R,
    char *cl,
    Uint len,
    char Q)
{
  Uint i,j;

  for(i=0; i< nooflists; i++) {
    for(j=0; j < LENGTHVEC(qdlist[i]); j++) {
      if(cl[VECTOR(qdlist[i],j)] == Q) {
        printf("sorting suffix %d at pos %d\n", VECTOR(qdlist[i],j), 
            R[VECTOR(qdlist[i],j)]);    
      }
    }
  }

}

Alurubucket*
getAluruBuckets(void *space, 
    char *s, 
    Uint len, 
    Uint *noofbuckets, 
    Uint **inv) {
  Uint i, j, no=0, b, k=0;
  Uint *srtidx;

  Alurubucket *bckts = NULL;

  for(i=0; i< len; i++) {      
    BUCKETRET(bckts, no, (Uint) s[i], b);
    if (b == no) {
      bckts = ALLOCMEMORY(space, bckts, Alurubucket, ++no);
      BUCKETINIT(bckts[b], (Uint) s[i]);
    } 
    BUCKETADD(space, bckts[b], i);               
  }

  srtidx = quickSort(space, bckts, no, bcktcmp, NULL);
  (*inv) = ALLOCMEMORY(space, NULL, Uint, len);

  for(i=0; i< no; i++) {
    for(j=0; j< bckts[srtidx[i]].noofelems; j++) {   
      printf("R[%d]=%d\n", bckts[srtidx[i]].elems[j], k );
      (*inv)[bckts[srtidx[i]].elems[j]] = k++;    
    }
  }

  qsort(bckts, no, sizeof(Alurubucket), bcktcmpANSI);
  *noofbuckets = no;
  return bckts;
}

vector_t**
getQdistList(void *space, 
    Alurubucket* bckts, 
    Uint bcktno, 
    Uint *d, 
    Uint len) {
  Uint max,
  i,
  j,
  l;
  vector_t **list;

  max = uarraymax (d, len); 
  list = ALLOCMEMORY(space, NULL, vector_t*, d[max]+1);

  for (i=0; i < d[max]+1; i++) {
    list[i] = ALLOCMEMORY(space, NULL, vector_t, 1);
    INITVECTOR(list[i]);
  }

  for(i=0; i < bcktno; i++) {
    for (j=0; j < bckts[i].elems[j]; j++) {
      l = d[bckts[i].elems[j]];
      printf("A[%d]=%d l=%d\n",i, bckts[i].elems[j], l); 
      appendvector(space, list[l], bckts[i].elems[j]);
    }
  }
  return list;
}

