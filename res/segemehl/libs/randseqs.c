
/*
 *  randseqs.c
 *  randomize sequences
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 07/15/2010 10:33:13 AM CEST
 *  
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "zlib.h"
#include "stringutils.h"
#include "basic-types.h"
#include "mathematics.h"
#include "biofiles.h"
#include "fileio.h"
#include "randseqs.h"
#include "charsequence.h"
#include "assert.h"
#include "zran.h"
#include "info.h"
#include "vtprogressbar.h"

double fiveacc = .95;
double threeacc = .9;

/*----------------------------- bl_fastxScramble -----------------------------
 *    
 * @brief scramble a fasta sequence
 * @author Steve Hoffmann 
 *   
 */

Uint
bl_fastxScramble (char *buffer, char *quality, 
    char *template, Uint len, 
    double acc, double Pmis, double Pins, 
    Uint uoff, Uint voff, 
    char *alphabet, Uint alphabetsize, 
    Uint minqual, Uint maxqual,
    char *editstring, Uint *editstringlen, Uint *readerrcnt)
{
  Uint u = uoff, k = (*editstringlen);
  int j;
  char errchr, chr;
  double errtype;

  for (j=0; j < len; j++) {
    chr=template[voff+j]; 
    if (RANDUNIT > acc) { 
      errtype = RANDUNIT;
      (*readerrcnt) += 1;
      
      if (errtype <= Pmis) {
        errchr = chr;
        if (editstring) k += sprintf(&editstring[k], "%d:S;", j); 
        while(errchr == chr)
          errchr = alphabet[RANDINT(alphabetsize-1)];
        quality[u] = minqual + RANDINT(maxqual - minqual);
        buffer[u++] = errchr; 
      
      } else if (errtype <= Pmis + Pins) {
        errchr = alphabet[RANDINT(alphabetsize-1)];
        if (editstring) k += sprintf(&editstring[k], "%d:I;", j); 
        quality[u] = minqual + RANDINT(maxqual - minqual);
        buffer[u++] = errchr;
        if(j < len-1) { 
          quality[u] = minqual + RANDINT(maxqual - minqual);
          buffer[u++] = chr;
          j++;
        }
      
      } else {
        if (editstring) k += sprintf(&editstring[k], "%d:D;", j); 
        j--;
      }
    
    } else {
      quality[u] = minqual + RANDINT(maxqual - minqual);
      buffer[u++] = chr;
    }
  
  }

  (*editstringlen) = k;
  return u;
}

/*-------------------------- bl_fastxPrintMatePairs --------------------------
 *    
 * @brief print mate pairs to dev and matedev
 * @author Steve Hoffmann 
 *   
 */

void
bl_fastxPrintRandomMatePairs(FILE *dev, FILE *matedev, 
    char *sequence, Uint seqlen, Uint n, 
    Uint minlen, Uint maxlen, Uint mindist, Uint maxdist,  
    char *alphabet, Uint alphabetsize,
    double acc,
    double Pmis, double Pins, double Pdel,
    unsigned char fastq,
    Uint minqual, Uint maxqual, 
    char *five, Uint fivelen,
    char *three, Uint threelen, Uint polyAlen)
{
  char *buffer, *quality, *matebuffer, *matequality, *polyAseq, *editstring, *rc,
  startchr='>';
  Uint i, u, start, matestart, len, matelen, readerrcnt, mateerrcnt, 
       fiveseqlen=fivelen, threeseqlen=threelen, polyAseqlen=polyAlen, editstringlen;

  assert(maxlen >= minlen);
  assert(maxdist >= mindist);
  assert(seqlen > 100 && seqlen > 5*(maxdist+2*maxlen+1)); 

  srand((unsigned int)time(NULL));

  if (fastq) {
    startchr = '@';
  }

  buffer = ALLOCMEMORY(space, NULL, char, 2*(maxlen+polyAlen)+threelen+fivelen+1);
  matebuffer = ALLOCMEMORY(space, NULL, char,  2*(maxlen+polyAlen)+threelen+fivelen+1);
  quality =  ALLOCMEMORY(space, NULL, char,  2*(maxlen+polyAlen)+threelen+fivelen+1);
  matequality = ALLOCMEMORY(space, NULL, char,  2*(maxlen+polyAlen)+threelen+fivelen+1);
  editstring = ALLOCMEMORY(space, NULL, char,  40*(maxlen+polyAlen+threelen+fivelen)+1);
  polyAseq = ALLOCMEMORY(space, NULL, char, 2*polyAlen);
  memset(polyAseq, 'A', 2*polyAlen);
  
  for (i=0; i < n; i++) {
    start=seqlen;

    memset(buffer, 0,  2*(maxlen+polyAlen)+threelen+fivelen+1);
    memset(matebuffer, 0,  2*(maxlen+polyAlen)+threelen+fivelen+1);
    memset(quality, 0,  2*(maxlen+polyAlen)+threelen+fivelen+1);
    memset(matequality, 0,  2*(maxlen+polyAlen)+threelen+fivelen+1);
    memset(editstring, 0,  40*(maxlen+polyAlen+threelen+fivelen)+1); 
    editstringlen = 0;
    polyAseqlen = polyAlen;

    while (start > seqlen-(maxdist+2*maxlen+1)) start = RANDINT(seqlen);
    len = minlen + RANDINT(maxlen - minlen); 

    if(RANDUNIT > fiveacc) {
      if(RANDUNIT > .5) {
        fiveseqlen = fivelen - RANDINT(fivelen);
      }
    }

    if(RANDUNIT > threeacc) {
      if(RANDUNIT > .5) {
        threeseqlen = threelen - RANDINT(threelen);
      } 
    }
    
    if(RANDUNIT > acc) {
      if(RANDUNIT > .5) {
        polyAseqlen = polyAlen - RANDINT(polyAlen);
      } else {
        polyAseqlen = polyAlen + RANDINT(polyAlen);
      }
    }
 

    matelen = minlen + RANDINT(maxlen - minlen);
    matestart = start + len - 1 + RANDINT(maxdist - mindist) + matelen - 1;

    fprintf(dev, "%c%d %d (len: %d) [", startchr, i, start+1, len);   
    fprintf(dev, "mate %d (len: %d) on rc] ", matestart+1, matelen); 
    fprintf(dev, "3'prime: %s (l:%d), 5'prime %s (l:%d), polyA: %d /",
        five, fiveseqlen, three, threeseqlen, polyAseqlen);

    fprintf(matedev, "%c%d %d (len: %d) [", startchr, i, start+1, len);   
    fprintf(matedev, "mate %d (len: %d) on rc] ", matestart+1, matelen); 
    fprintf(matedev, "3'prime: %s (l:%d), 5'prime %s (l:%d), polyA: %d /",
        five, fiveseqlen, three, threeseqlen, polyAseqlen);

    readerrcnt = 0;
    u = 0;

    u = bl_fastxScramble (buffer, quality, five, fiveseqlen, fiveacc, Pmis, Pins, 
        u, 0, alphabet, alphabetsize, minqual, maxqual, 
        editstring, &editstringlen, &readerrcnt);

    editstring[editstringlen++] = '|'; 
    u = bl_fastxScramble (buffer, quality, sequence, len, acc, Pmis, Pins, 
        u, start, alphabet, alphabetsize, minqual, maxqual, 
        editstring, &editstringlen, &readerrcnt);

    editstring[editstringlen++] = '|'; 
    u = bl_fastxScramble (buffer, quality, polyAseq, polyAseqlen, acc, Pmis, Pins, 
        u, 0, alphabet, alphabetsize, minqual, maxqual, 
        editstring, &editstringlen, &readerrcnt);

    editstring[editstringlen++] = '|'; 
    u = bl_fastxScramble (buffer, quality, three, threeseqlen, threeacc, Pmis, Pins, 
        u, 0, alphabet, alphabetsize, minqual, maxqual, 
        editstring, &editstringlen, &readerrcnt);
    
    fprintf(dev, " %s", editstring);
    
    memset(editstring, 0,  40*(maxlen+polyAlen+threelen+fivelen)+1); 
    editstringlen = 0;

    mateerrcnt = 0;
    rc = charDNAcomplement(NULL, &sequence[matestart-matelen+1], matelen);
    u = 0;

    u = bl_fastxScramble (matebuffer, matequality, five, fiveseqlen, fiveacc, Pmis, Pins, 
        u, 0, alphabet, alphabetsize, minqual, maxqual, 
        editstring, &editstringlen, &mateerrcnt);

    editstring[editstringlen++] = '|'; 
    u = bl_fastxScramble (matebuffer, matequality, rc, matelen, acc, Pmis, Pins, 
        u, 0, alphabet, alphabetsize, minqual, maxqual, 
        editstring, &editstringlen, &mateerrcnt);

    editstring[editstringlen++] = '|'; 
    u = bl_fastxScramble (matebuffer, matequality, polyAseq, polyAseqlen, acc, Pmis, Pins, 
        u, 0, alphabet, alphabetsize, minqual, maxqual, 
        editstring, &editstringlen, &mateerrcnt);

    editstring[editstringlen++] = '|'; 
    u = bl_fastxScramble (matebuffer, matequality, three, threeseqlen, threeacc, Pmis, Pins, 
        u, 0, alphabet, alphabetsize, minqual, maxqual, 
        editstring, &editstringlen, &mateerrcnt);


    FREEMEMORY(space, rc);


    fprintf(dev, " (r: %f) [r: %f]\n", (double) readerrcnt/(len+1.0),
        (double) mateerrcnt/(matelen+1.0));

    fprintf(matedev, " %s", editstring);
    fprintf(matedev, " (r: %f) [r: %f]\n", (double) readerrcnt/(len+1.0),
        (double) mateerrcnt/(matelen+1.0));
    
    fprintf(dev, "%s\n", buffer);
    fprintf(matedev, "%s\n", matebuffer);
    
    if(fastq) {
      fprintf(dev, "+ \n");
      fprintf(dev, "%s\n\n", quality);
      fprintf(matedev, "+ \n");
      fprintf(matedev, "%s\n\n", matequality);
    }
  }

  FREEMEMORY(space, buffer);
  FREEMEMORY(space, quality);
  FREEMEMORY(space, matebuffer);
  FREEMEMORY(space, matequality);
  FREEMEMORY(space, polyAseq);
  FREEMEMORY(space, editstring);
  return ;
}


/*--------------------------- bl_generateGeneModel ---------------------------
 *    
 * @brief generate a random gene model
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_generateGeneModel (void *space, Uint n, Uint minexoncnt, Uint maxexoncnt, 
    Uint minexonlen, Uint maxexonlen, Uint minexondist, Uint maxexondist)
{
  Uint i, j, exoncnt, exonlen, exondist;


  for(i=0; i < n; i++) { 
    exoncnt = minexoncnt + RANDINT(maxexoncnt-minexoncnt);
    for(j=0; j < exoncnt; j++) { 
      exonlen = minexonlen + RANDINT(maxexonlen - minexonlen);
      if(j>0) {
        exondist = minexondist + RANDINT(maxexondist - minexondist);
        // dummy info to avoid unused variables
        if (0){
          fprintf(stderr, "exoncnt=%u\texonlen=%u\texondist%u\n",
                  exoncnt, exonlen, exondist);
        }
      }
    }
  }
  return ;
}



/*----------------------- bl_fastxSimulateSpliceSites ------------------------
 *    
 * @brief simulate splicesites and split reads accross splice sites
 * @author Steve Hoffmann 
 *   
 */

  void
bl_fastxSimulateSpliceSites (void *space, char *sequence,
    Uint seqlen, Uint n, Uint maxchildren, 
    char *alphabet, Uint alphabetsize, double acc, 
    double Pmis, double Pins, double Pdel,
    Uint minqual, Uint maxqual,
    double Pcis, Uint mincisdist, Uint maxcisdist, 
    double Pstrandswitch, Uint readlen)
{

  Uint parent;
  Uint child=0;
  Uint maxcoverage = 200;
  Uint mincoverage = 1;
  Uint c;
  Uint randrange;
  Uint i, j, k, l, e, m, effsplit;
  int u, v;
  Uint parentlen, lastparentlen = 0, childlen;
  Uint readerrcnt;
  Uint editstringlen;
  double errtype;
  char *buffer, *quality, *editstring, chr, errchr;
  unsigned parentstrand, childrc=0, complement = 0;

  assert(seqlen > 2*readlen);
  randrange = seqlen - (2*readlen);

  buffer = ALLOCMEMORY(space, NULL, char, 2*readlen+1);
  quality =  ALLOCMEMORY(space, NULL, char, 2*readlen+1);
  editstring = ALLOCMEMORY(space, NULL, char, 40*(2*readlen)+1);

  /*generate parent splice site*/
  for(i=0; i < n; i++) {
    parent = readlen + RANDINT(randrange);
    child = 0;
    parentstrand = (RANDUNIT > 0.5) ? 1 : 0;
    m = RANDINT(maxchildren);
    m = MAX(1, m);

    /*generate the children*/
    for(j=0; j < m; j++) {
      childrc = 0;

      if(RANDUNIT < Pcis) {
        child = RANDINT(maxcisdist-mincisdist) + mincisdist;
        if(RANDUNIT > 0.5 && parent + child + readlen < seqlen) {
          child = parent + child; 
        } else if(parent > child + readlen) { 
          child = parent - child;
        } 
      } else {
        if (RANDUNIT > Pstrandswitch) {
          child = RANDINT(maxcisdist-mincisdist) + mincisdist;
          if(RANDUNIT > 0.5 && parent + child + readlen < seqlen) {
            child = parent + child; 
          } else if(parent > child + readlen) { 
            child = parent - child;
          } 
          childrc = 1;
        } else {
          child = readlen + RANDINT(randrange);
        }
      }

      parentstrand = 0;
      childrc=0;
      /*do the sequence work*/
      if(parent && child) { 
        c = RANDINT(maxcoverage); 
        c = MAX(mincoverage, c);
        fprintf(stdout, "subject\t%d\t%d\tsplicetype%d-%d\t%d\n", parent+1, child+1, parentstrand, childrc, c);
        /*cover the splice sites*/
        for(k=0; k < c; k++) {

          memset(buffer, 0, 2*readlen+1);
          memset(quality, 0, 2*readlen+1);
          memset(editstring, 0, 40*(2*readlen)+1);

          do { 
            parentlen = RANDINT(readlen-22);
            parentlen = MAX(22, parentlen); 
          } while (parentlen == lastparentlen);

          lastparentlen = parentlen;

          childlen = readlen-parentlen;
          readerrcnt = 0;
          editstringlen = 0;
          effsplit = 0;
          e = 0;

          if(!parentstrand)
            editstringlen = sprintf(editstring, "splice:%d", parent-parentlen);
          else
            editstringlen = sprintf(editstring, "splice:%d", parent + 1);

 
          if((parentstrand && childrc) || (!parentstrand && !childrc)) {
            editstringlen += sprintf(&editstring[editstringlen], "-%d;%d;", child-childlen, k);
          } else {
            editstringlen += sprintf(&editstring[editstringlen], "-%d;%d;", child+1, k);
          }


          /*construct the sequence*/
          for(l=0; l < parentlen+childlen; l++) { 
 
            if (l == parentlen) effsplit = (e > 0) ? e-1: 0;

            if(l < parentlen) { 
              if(!parentstrand) { 
                v = parent - parentlen + 1;
                u = l;
                complement = 0;
              } else { 
                v = parent + parentlen - 1;
                u = -1*l;  
                complement = 1;
              }
            } else {
              if((parentstrand && childrc) || (!parentstrand && !childrc)) {
                v = child;
                u = l-parentlen;
                complement = 0;
              } else {
                v = child;
                u = -1*(l-parentlen);
                complement = 1;
              }
            }

            chr=sequence[v+u]; 
            if(complement) chr = charComplementChar(chr);

            if (RANDUNIT > acc) { 
              errtype = RANDUNIT;
              readerrcnt++;
              if (errtype <= Pmis) {
                errchr = chr;
                editstringlen += sprintf(&editstring[editstringlen], "%d:S;",u);
                while(errchr == chr) errchr = alphabet[RANDINT(alphabetsize)];
                quality[e] = minqual + RANDINT(maxqual - minqual);
                buffer[e++] = errchr; 
              } else if (errtype <= Pmis + Pins) {
                errchr = alphabet[RANDINT(alphabetsize)];
                editstringlen += sprintf(&editstring[editstringlen], "%d:I;",u);
                quality[e] = minqual + RANDINT(maxqual - minqual);
                buffer[e++] = errchr;
                quality[e] = minqual + RANDINT(maxqual - minqual);
                buffer[e++] = chr;
              } else {
                editstringlen += sprintf(&editstring[editstringlen], "%d:D;",u);
              }
            } else {
              quality[e] = minqual + RANDINT(maxqual - minqual);
              buffer[e++] = chr;
            }

          } /*construct the sequence*/


          if(!parentstrand) {
            editstringlen += sprintf(&editstring[editstringlen], "(+:");
          } else { 
            editstringlen += sprintf(&editstring[editstringlen], "(-:");
          }
          
          if((parentstrand && childrc) || (!parentstrand && !childrc)) {
            editstringlen += sprintf(&editstring[editstringlen], "+)");
          } else {
            editstringlen += sprintf(&editstring[editstringlen], "-)");
          }

         
          if(e-1-effsplit >= 0);
          assert(effsplit <= strlen(buffer));
          fprintf(stderr, "@%s\n%s\n+%s\n%s\n", editstring, buffer, editstring, quality);
        
        } /*end of coverage loop*/
      } /*end of sequence work*/

    }
    fprintf(stderr, "\n");
  }


  return ;
}

/*------------------------ bl_fastxGenerateSplitReads ------------------------
 *    
 * @brief generate a set of spliced reads
 * @author Steve Hoffmann 
 *   
 */
 
  void
bl_fastxPrintRandomSplitReads (FILE *dev, char *sequence, 
    Uint seqlen, Uint n, 
    Uint minspltlen, Uint maxspltlen,  
    char *alphabet, Uint alphabetsize,
    double acc,
    double Pmis, double Pins, double Pdel,
    unsigned char fastq,
    Uint minqual, Uint maxqual,
    char *five, Uint fivelen,
    char *three, Uint threelen, Uint polyAlen)
{
  char *buffer, *rc, *quality, *fivebuffer=NULL, *fivequality=NULL,
  *threebuffer=NULL, *threequality=NULL, 
  *polyAbuffer=NULL, *polyAquality=NULL, *polyAseq=NULL, 
  *editstring, startchr='>', chr, errchr=0;
  unsigned char hasReverse, isPartial, inDownstream;
  double errtype;
  int i, j, u, q, k, start, start2, spltlen, spltlen2,effsplit=0;
  Uint threeseqlen=threelen, fiveseqlen=fivelen, polyAseqlen=polyAlen, 
       readerrcnt=0, editstringlen=0;

  assert(maxspltlen >= minspltlen);
  assert(seqlen > 100 && seqlen > 5*maxspltlen); 

  srand((unsigned int)time(NULL));
  if (five) {
    fivebuffer = ALLOCMEMORY(space, NULL, char, 2*fivelen+3);
    fivequality = ALLOCMEMORY(space, NULL, char, 2*fivelen+3);
  }

  buffer = ALLOCMEMORY(space, NULL, char, 2*maxspltlen+3);
  quality =  ALLOCMEMORY(space, NULL, char, 2*maxspltlen+3);
  editstring = ALLOCMEMORY(space, NULL, char, 40*(2*maxspltlen+polyAlen+threelen+fivelen));

  if(polyAlen) {
    polyAbuffer = ALLOCMEMORY(space, NULL, char, 2*polyAlen+3);
    polyAquality = ALLOCMEMORY(space, NULL, char, 2*polyAlen+3);
  }
  
  if (three) {
    threebuffer = ALLOCMEMORY(space, NULL, char, 2*threelen+3);
    threequality = ALLOCMEMORY(space, NULL, char, 2*threelen+3);
  }

  polyAseq = ALLOCMEMORY(space, NULL, char, 2*polyAlen);
  memset(polyAseq, 'A', 2*polyAlen);

  for (i=0; i < n; i++) {
    start=seqlen;
    memset(buffer, 0, 2*maxspltlen+3);
    memset(quality, 0, 2*maxspltlen+3);
    memset(editstring, 0, 40*(2*maxspltlen+polyAlen+threelen+fivelen));
    editstringlen = 0;
    readerrcnt = 0; 
    polyAseqlen = polyAlen;

    hasReverse = isPartial = inDownstream = 0;
    if(RANDINT(10) > 5) hasReverse = 1;
    if(RANDINT(10) > 5) isPartial = 1;
    if(RANDINT(10) > 5) inDownstream = 1;

    while (start > seqlen-(2*maxspltlen)) start = RANDINT(seqlen);
    spltlen = minspltlen + RANDINT(maxspltlen - minspltlen); 
    start2 = start + spltlen + RANDINT(seqlen - start - spltlen - maxspltlen);
    spltlen2 = minspltlen + RANDINT(maxspltlen - minspltlen);

   if(RANDUNIT > fiveacc) {
      if(RANDUNIT > .5) {
        fiveseqlen = fivelen - RANDINT(fivelen);
      }     
   }

    if(RANDUNIT > threeacc) {
      if(RANDUNIT > .5) {
        threeseqlen = threelen - RANDINT(threelen);
      }  
    }
    
    if(RANDUNIT > acc) {
      if(RANDUNIT > .5) {
        polyAseqlen = polyAlen - RANDINT(polyAlen);
      } else {
        polyAseqlen = polyAlen + RANDINT(polyAlen);
      }
    }

    fprintf(dev, "%c%d %d(len: %d)-", startchr, i, start+1, spltlen);   
    fprintf(dev, "%d(len: %d) ", start2+1, spltlen2);   
    u = 0;

    for (k=0; k < spltlen+spltlen2; k++) {
      j = (k < spltlen) ? k : k-spltlen;
      q = (k < spltlen) ? start : start2;
      if (k == spltlen) effsplit = (u > 0) ? u-1: 0;

      chr=sequence[q+j]; 

      if (RANDUNIT > acc) { 
        errtype = RANDUNIT;
        readerrcnt++;
        if (errtype <= Pmis) {
          errchr = chr;
          editstringlen += sprintf(editstring, "%d:S;",j);
          while(errchr == chr)
            errchr = alphabet[RANDINT(alphabetsize)];
          quality[u] = minqual + RANDINT(maxqual - minqual);
          buffer[u++] = errchr; 
        } else if (errtype <= Pmis + Pins) {
          errchr = alphabet[RANDINT(alphabetsize)];
          editstringlen += sprintf(editstring, "%d:I;",j);
          quality[u] = minqual + RANDINT(maxqual - minqual);
          buffer[u++] = errchr;
          quality[u] = minqual + RANDINT(maxqual - minqual);
          buffer[u++] = chr;
        } else {
          editstringlen += sprintf(editstring, "%d:D;",j);
        }
      } else {
        quality[u] = minqual + RANDINT(maxqual - minqual);
        buffer[u++] = chr;
      }
    }
           
    assert(u-1-effsplit >= 0);
    assert(effsplit <= strlen(buffer));

    if (hasReverse) {
      if(isPartial) {
        if(inDownstream) {
           rc = charDNAcomplement(NULL, &buffer[effsplit], strlen(buffer)-effsplit);
           memmove(&buffer[effsplit], rc, strlen(buffer)-effsplit);
           free(rc);
           editstringlen += sprintf(editstring, " (+/-)");
        } else {
           rc = charDNAcomplement(NULL, buffer, effsplit);
           memmove(buffer, rc, effsplit);
           free(rc);
           editstringlen = sprintf(editstring, " (-/+))");
        }
      } else {
        rc = charDNAcomplement(NULL, buffer, strlen(buffer));
        memmove(buffer, rc, strlen(buffer));
        editstringlen = sprintf(editstring, " (-/-)");
        FREEMEMORY(space, rc);
      }
    } else {
      editstringlen += sprintf(editstring, " (+/+)");
    }

    editstring[editstringlen++] = '|'; 
    if (five) {  
    memset(fivebuffer, 0, 2*fivelen+3);
    memset(fivequality, 0, 2*fivelen+3);
    u = bl_fastxScramble (fivebuffer, fivequality, five, fiveseqlen, 
        fiveacc, Pmis, Pins, 0, 0, alphabet, alphabetsize, minqual, maxqual, 
        editstring, &editstringlen, &readerrcnt);
    }
    
    editstring[editstringlen++] = '|'; 
    if(polyAlen) {
      memset(polyAbuffer, 0, 2*polyAlen+3);  
      memset(polyAquality, 0, 2*polyAlen+3);
      u = bl_fastxScramble (polyAbuffer, polyAquality, polyAseq, polyAseqlen, 
          acc, Pmis, Pins, 0, 0, alphabet, alphabetsize, minqual, maxqual, 
          editstring, &editstringlen, &readerrcnt);
    }

    editstring[editstringlen++] = '|'; 
    if(three) {
    memset(threebuffer, 0, 2*threelen+3);
    memset(threequality, 0, 2*threelen+3);
    u = bl_fastxScramble (threebuffer, threequality, three, threeseqlen, 
        threeacc, Pmis, Pins, 0, 0, alphabet, alphabetsize, minqual, maxqual, 
        editstring, &editstringlen, &readerrcnt);
    }
   
    fprintf(dev, " (r: %f)\n", (double) readerrcnt/(spltlen+spltlen2+1.0));
    if (five) fprintf(dev, "%s", fivebuffer);
    fprintf(dev, "%s", buffer);
    if (polyAlen) fprintf(dev, "%s", polyAbuffer);
    if (three) fprintf(dev, "%s", threebuffer);
    fprintf(dev, "\n");
    
    if(fastq) {
      fprintf(dev, "+ \n");
      if(five) fprintf(dev,"%s", fivequality);
      fprintf(dev, "%s", quality);
      if(three) fprintf(dev, "%s", threequality);
      fprintf(dev, "\n\n");
    }
  }

  if(five) {
    FREEMEMORY(space, fivebuffer);
    FREEMEMORY(space, fivequality);
  }
  
  FREEMEMORY(space, buffer);
  FREEMEMORY(space, quality);
  FREEMEMORY(space, editstring);

  if(three) {
    FREEMEMORY(space, threebuffer);
    FREEMEMORY(space, threequality);
  }
  if (polyAlen) {
    FREEMEMORY(space, polyAbuffer);
    FREEMEMORY(space, polyAquality);

  }
  FREEMEMORY(space, polyAseq);
  return ;
}


/*-------------------------- bl_fastxPrintRefReadSet -------------------------
 *    
 * @brief print a reference and a set of simulated reads in fastq format
 *        to device dev
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_fastxPrintRandomReads(
    FILE *reads, 
    char *sequence,
    Uint seqlen,
    Uint n,
    Uint minlen, 
    Uint maxlen, 
    char *alphabet, Uint alphabetsize,
    double acc, 
    double Pmis, double Pins, double Pdel, 
    unsigned char fastq,
    Uint minqual, Uint maxqual,
    char *five, Uint fivelen,
    char *three, Uint threelen, Uint polyAlen)
{

  Uint i, k, start, len, readerrcnt, fiveseqlen=fivelen, threeseqlen=threelen, 
  polyAseqlen=0, editstringlen=0;
  char *quality, *buffer, *editstring, *polyAseq, startchr = '>';
  
  if (fastq) {
    startchr = '@';
    assert(maxqual <= 126 && minqual >=33);
  }

  assert(seqlen > maxlen);
  assert(minlen <= maxlen);
  assert(alphabetsize > 1);
  assert(Pmis + Pins + Pdel == 1);

  srand((unsigned int)time(NULL));

  quality = ALLOCMEMORY(space, NULL, char, 2*(maxlen+fivelen+threelen+polyAlen));
  buffer = ALLOCMEMORY(space, NULL, char, 2*(maxlen+fivelen+threelen+polyAlen));
  editstring = ALLOCMEMORY(space, NULL, char, 40*(maxlen+fivelen+threelen+polyAlen));
  polyAseq = ALLOCMEMORY(space, NULL, char, 2*polyAlen);
  memset(polyAseq, 'A', 2*polyAlen);

  for(i=0; i < n; i++){
    memset(buffer, 0, sizeof(char)*(2*(maxlen+fivelen+threelen+polyAlen)));
    memset(quality, 0, sizeof(char)*(2*(maxlen+fivelen+threelen+polyAlen)));
    memset(editstring, 0, sizeof(char)*(40*(maxlen+fivelen+threelen+polyAlen)));
    editstringlen = 0;
    polyAseqlen = polyAlen;

    len = minlen + RANDINT(maxlen - minlen);
    start = RANDINT(seqlen-len);
   
    if(RANDUNIT > fiveacc) {
      if(RANDUNIT > .5) {
        fiveseqlen = fivelen - RANDINT(fivelen);
      }
    }

    if(RANDUNIT > threeacc) {
      if(RANDUNIT > .5) {
        threeseqlen = threelen - RANDINT(threelen);
      }
    }
    
    if(RANDUNIT > acc) {
      if(RANDUNIT > .5) {
        polyAseqlen = polyAlen - RANDINT(polyAlen);
      } else {
        polyAseqlen = polyAlen + RANDINT(polyAlen);
      }
    }
    fprintf(reads, "%c%d %d (len: %d) ", startchr, i, start+1, len);

    readerrcnt=0;
    k=0; 
     
    k = bl_fastxScramble (buffer, quality, five, fiveseqlen, fiveacc, Pmis, Pins, 
        k, 0, alphabet, alphabetsize, minqual, maxqual, 
        editstring, &editstringlen, &readerrcnt);
    
    editstring[editstringlen++] = '|'; 

    k = bl_fastxScramble (buffer, quality, sequence, len, acc, Pmis, Pins, 
        k, start, alphabet, alphabetsize, minqual, maxqual, 
        editstring, &editstringlen, &readerrcnt);

    editstring[editstringlen++] = '|'; 

    k = bl_fastxScramble (buffer, quality, polyAseq, polyAseqlen, acc, Pmis, Pins, 
        k, 0, alphabet, alphabetsize, minqual, maxqual, 
        editstring, &editstringlen, &readerrcnt);

    editstring[editstringlen++] = '|'; 
    
    k = bl_fastxScramble (buffer, quality, three, threeseqlen, threeacc, Pmis, Pins, 
        k, 0, alphabet, alphabetsize, minqual, maxqual, 
        editstring, &editstringlen, &readerrcnt);
    
    fprintf(reads, "%s", editstring);
    fprintf(reads, " (r: %f)\n", (double) readerrcnt/(len+1.0));
    fprintf(reads, "%s\n", buffer);
    if (fastq) {
      fprintf(reads, "+\n");
      fprintf(reads, "%s\n", quality);
    }
  }

  FREEMEMORY(space, buffer);
  FREEMEMORY(space, quality);
  FREEMEMORY(space, polyAseq);
  FREEMEMORY(space, editstring);
  
  return ;
}


/*---------------------------- bl_fastaPrintRandom ---------------------------
 *    
 * @brief print a random fasta sequence to device dev
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_fastaPrintRandom(FILE *dev, Uint noofseqs, Uint minlen, Uint maxlen)
{
  Uint i, j, len, chr;
  unsigned int iseed = (unsigned int)time(NULL);
  
  assert(minlen <= maxlen);
  srand(iseed);
  for(i=0; i < noofseqs; i++){
    len = minlen + ((maxlen - minlen) * (rand() / (RAND_MAX + 1.0)));
    fprintf(dev, ">random sequence %d (len: %d)\n", i, len);
    for(j=0; j < len; j++) {
      chr  = 65 + (62 * (rand() / (RAND_MAX + 1.0)));
      fprintf(dev, "%c", chr);
    }
    fprintf(dev,"\n");
  }
  return ;
}




/*------------------------- bl_getTrackFromGeneModel -------------------------
 *    
 * @brief get a annotation track from  a gene model
 * @author Steve Hoffmann 
 *   
 */
 
annotationtrack_t *
bl_getTrackFromGeneModel (void *space, geneset_t *set)
{
  Uint i, j;
  annotationtrack_t *track;

  track = ALLOCMEMORY(space, NULL, annotationtrack_t, 1);
  track->items = ALLOCMEMORY(space, NULL, annotationitem_t, set->noofgenes);

  for(i=0; i < set->noofgenes; i++) {
    
    track->items[i].start = set->genes[i].exons[0].start;
    track->items[i].end = set->genes[i].exons[set->genes[i].noofexons-1].end;
    track->items[i].strand = set->genes[i].direction;
    track->items[i].thickStart = set->genes[i].startcodon;
    track->items[i].thickEnd = set->genes[i].stopcodon;
    track->items[i].blockCount = set->genes[i].noofexons;
    track->items[i].score = 0;
    
    track->items[i].blockSizes = 
      ALLOCMEMORY(space, NULL, Uint, set->genes[i].noofexons);
    track->items[i].blockStarts =
      ALLOCMEMORY(space, NULL, Uint, set->genes[i].noofexons);
    track->items[i].blockRefseqs =
      ALLOCMEMORY(space, NULL, char*, set->genes[i].noofexons);
    track->items[i].blockStrands =
      ALLOCMEMORY(space, NULL, char, set->genes[i].noofexons);
    track->items[i].itemRgb = 
      ALLOCMEMORY(space, NULL, Uint, 3);
    memset(track->items[i].itemRgb, 0, 3*sizeof(Uint));
    track->items[i].chromname = 
      ALLOCMEMORY(space, NULL, char, strlen(set->genes[i].exons[0].refchr)+1);
    memmove(track->items[i].chromname, set->genes[i].exons[0].refchr,
        strlen(set->genes[i].exons[0].refchr)+1);
    track->items[i].chromname[strlen(set->genes[i].exons[0].refchr)]=0;
    
    track->items[i].name = 
      ALLOCMEMORY(space, NULL, char, strlen(set->genes[i].id));
    memmove(track->items[i].name, set->genes[i].id, strlen(set->genes[i].id));
    track->items[i].name[strlen(set->genes[i].id)] = 0;

    for(j=0; j < set->genes[i].noofexons; j++) {
      track->items[i].blockSizes[j] = 
        set->genes[i].exons[j].end - set->genes[i].exons[j].start + 1;
      track->items[i].blockStrands[j] = 
        set->genes[i].exons[j].strand; 
      track->items[i].itemRgb[j] = 0;
      if(strcmp(set->genes[i].exons[j].refchr, track->items[i].chromname) || 
          track->items[i].blockStrands[j] != track->items[i].strand) { 
      track->items[i].blockRefseqs[j] = 
        ALLOCMEMORY(space, NULL, char, strlen(set->genes[i].exons[j].refchr));
      memmove(track->items[i].blockRefseqs[j], set->genes[i].exons[j].refchr, 
          strlen(set->genes[i].exons[j].refchr));
      track->items[i].blockRefseqs[j][strlen(set->genes[i].exons[j].refchr)] = 0;
      
      track->items[i].blockStarts[j] = 
        set->genes[i].exons[j].start;
      } else {
      track->items[i].blockRefseqs[j] = NULL; 
      
      track->items[i].blockStarts[j] = 
        set->genes[i].exons[j].start - track->items[i].start;
      }
    }
  }
 
  track->noofitems = set->noofgenes;
  return track;
}


/*------------------------ bl_getGeneModelFromBEDTrack -------------------------
 *    
 * @brief get the gene model from an annotation track
 * @author Steve Hoffmann 
 *   
 */
 
geneset_t *
bl_getGeneModelFromBEDtrack (void *space, annotationtrack_t *track)
{
  Uint i, j, start, blockCount;//, end;
  Uint *blockSizes;
  Uint *blockStarts;
  gene_t *genes;
  geneset_t *geneset;

  geneset = ALLOCMEMORY(space, NULL, geneset_t, 1);
  genes = ALLOCMEMORY(space, NULL, gene_t, track->noofitems);

  for(i=0; i < track->noofitems; i++) {
    
    start = track->items[i].start;
    // not used: end = track->items[i].end;
    blockCount = track->items[i].blockCount;
    blockSizes = track->items[i].blockSizes;
    blockStarts = track->items[i].blockStarts;

    genes[i].id = track->items[i].name;
    genes[i].startcodon = track->items[i].thickStart;
    genes[i].stopcodon = track->items[i].thickEnd;
    genes[i].noofexons = blockCount;
    genes[i].direction = track->items[i].strand;
    genes[i].exons = ALLOCMEMORY(space, NULL, exon_t, blockCount);

    for(j=0; j < blockCount; j++) {
      
      if (track->items[i].blockRefseqs[j]) { 
        genes[i].exons[j].refchr = track->items[i].blockRefseqs[j];
        genes[i].exons[j].start = blockStarts[j];
      } else {
        genes[i].exons[j].refchr = track->items[i].chromname;
        genes[i].exons[j].start = start + blockStarts[j];
      }

      genes[i].exons[j].strand = track->items[i].blockStrands[j];

      genes[i].exons[j].end = genes[i].exons[j].start + blockSizes[j] -1;
      genes[i].exons[j].noofcds = 0;
      genes[i].exons[j].cds = NULL;
    }  
  }
	 
  geneset->noofgenes = track->noofitems;
  geneset->genes = genes;

  return geneset;
}


/*------------------------------- bl_copyGene --------------------------------
 *    
 * @brief copy the gene
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_copyGene (void *space, gene_t *to, gene_t *from)
{
  Uint i;

  to->id = from->id;
  to->startcodon = from->startcodon;
  to->stopcodon = from->stopcodon;
  to->noofexons = from->noofexons;
  to->direction = from->direction;
  to->exons = ALLOCMEMORY(space, NULL, exon_t, from->noofexons);

  for(i=0; i < from->noofexons; i++) {
    to->exons[i].refchr = from->exons[i].refchr;
    to->exons[i].strand = from->exons[i].strand;
    to->exons[i].start = from->exons[i].start;
    to->exons[i].end = from->exons[i].end;
    to->exons[i].noofcds = 0;
    to->exons[i].cds = NULL;
  }

  return ;
}

/*------------------------- bl_simulateTransSplicing -------------------------
 *    
 * @brief simulate random  transsplicing n events from a given gene model
 * @author Steve Hoffmann 
 *   
 */

  geneset_t*
bl_simulateTransSplicing (void *space, geneset_t *set, char type, Uint n)
{
  geneset_t *transset;
  gene_t* genes;
  Uint i, u, v, w;

  transset = ALLOCMEMORY(space, NULL, geneset_t, 1);
  genes = ALLOCMEMORY(space, NULL, gene_t, n);

  for(i=0; i < n; i++) {
    u = RANDINT(set->noofgenes-1);
    switch(type) {
      
      case 'S':
        bl_copyGene(space, &genes[i], &set->genes[u]);
        w = RANDINT(genes[i].noofexons-1);
         
        if(genes[i].exons[w].strand == '+') { 
          genes[i].exons[w].strand = '-';
          if(genes[i].noofexons == 1) {
            genes[i].direction= '-';
          }
        } else { 
          genes[i].exons[w].strand = '+';
          if(genes[i].noofexons == 1) {
            genes[i].direction = '+';
          }
        }

        break;
      
      case 'D':
        bl_copyGene(space, &genes[i], &set->genes[u]);
        v = RANDINT(set->noofgenes-1);
        w = RANDINT(set->genes[v].noofexons-1);
        genes[i].exons = 
          ALLOCMEMORY(space, genes[i].exons, exon_t, genes[i].noofexons+1);
        genes[i].exons[genes[i].noofexons].refchr = 
          set->genes[v].exons[w].refchr;
        genes[i].exons[genes[i].noofexons].strand =
          set->genes[v].exons[w].strand;
        genes[i].exons[genes[i].noofexons].start =
          set->genes[v].exons[w].start;
        genes[i].exons[genes[i].noofexons].end =
          set->genes[v].exons[w].end;
        genes[i].exons[genes[i].noofexons].noofcds = 0;
        genes[i].exons[genes[i].noofexons].cds = NULL;

        genes[i].noofexons++;
        break;
      
      default:
        break;
    }
  }

  transset->noofgenes = n;
  transset->genes = genes;
  return transset;
}


/*-------------------------- bl_printSplicingEdges ---------------------------
 *    
 * @brief print the splicing edges 
 * @author Steve Hoffmann 
 *   
 */

  void
bl_printSplicingEdges (void *space, FILE *dev, geneset_t *set)
{
  Uint i, j, k, direction, next;
  Uint spliceBpos, spliceApos;
  char *spliceAchr, *spliceBchr;

  for(i=0; i < set->noofgenes; i++) {
    direction = set->genes[i].direction;
    for(j=0; j < set->genes[i].noofexons; j++) {
      next = -1;

      if(direction == '+') {
        k = j;
        if(j < set->genes[i].noofexons-1) { 
          next = j+1;
        }
      } else {
        k = set->genes[i].noofexons - j - 1;
        if(j < set->genes[i].noofexons-1) {
          next = set->genes[i].noofexons - j - 2;
        }
      }

      if(next != -1) { 
        spliceAchr = set->genes[i].exons[k].refchr;
        if(set->genes[i].exons[k].strand == '-') { 
          spliceApos = set->genes[i].exons[k].start;
        } else { 
          spliceApos = set->genes[i].exons[k].end;
        }


        spliceBchr = set->genes[i].exons[next].refchr;
        if(set->genes[i].exons[next].strand == '-') { 
          spliceBpos = set->genes[i].exons[next].end;
        } else { 
          spliceBpos = set->genes[i].exons[next].start;
        }

        fprintf(dev, "%s\t%d\t%c\t%s\t%d\t%c\n", spliceAchr, spliceApos, 
            set->genes[i].exons[k].strand, spliceBchr, spliceBpos, 
            set->genes[i].exons[next].strand);
      }
    }
  }
  return ;
}

/*--------------------------- bl_getGeneSequences ----------------------------
 *    
 * @brief get the sequences from a gene model
 * @author Steve Hoffmann 
 *   
 */
 
char*
bl_getGeneSequence(void *space, fasta_t *reference, gene_t *gene)
{
  Uint i, start, end, seqlen = 0, id=0, k = 0;
  char *chr;
  char *sequence = NULL, *buffer, *test;

  for(i=0; i < gene->noofexons; i++) { 
    
    if(gene->direction == '-') {
      k = gene->noofexons - 1 - i;
    } else { 
      k = i;
    }
 
    start = gene->exons[k].start;
    end = gene->exons[k].end;

    sequence = ALLOCMEMORY(space, sequence, char, seqlen+(end-start)+2);
    id = bl_fastxFindIDIdx (gene->exons[k].refchr, reference);
    chr = bl_fastaGetSequence(reference, id);
    
    if(gene->exons[k].strand == '+') { 
       test = ALLOCMEMORY(space, NULL, char, (end-start)+2);
      memmove(test, &chr[start], (end-start)+1);
      test[(end-start)+1] = 0;
//      fprintf(stdout,"[%d,%d]\ns:%s\n\n", start, end, test);
      memmove(&sequence[seqlen], &chr[start], (end-start)+1);

    } else { 
      test = ALLOCMEMORY(space, NULL, char, (end-start)+2);
      memmove(test, &chr[start], (end-start)+1);
      test[(end-start)+1] = 0;
      buffer = charDNAcomplement(space, &chr[start], (end-start)+1);
//      fprintf(stdout,"[%d,%d]\ns:%s\nb:%s\n\n", start, end, test, buffer);
      memmove(&sequence[seqlen], buffer, (end-start)+1);
      FREEMEMORY(space, buffer);
    }
    seqlen += (end-start)+1;
    sequence[seqlen] = 0;
  }

  return sequence;
}


/*------------------------ bl_simulateGeneSequencing -------------------------
 *    
 * @brief generate a set of simulated reads for a given gene
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_simulateGeneSequencing (void *space, FILE *dev, fasta_t *reference, gene_t *gene, Uint readlen, 
    Uint cov, char *alphabet, Uint alphabetsize, Uint minqual, Uint maxqual,
    double acc, double Pmis, double Pins)
{

  char *quality, *buffer, *template, *editstring;
  char *sequence;
  Uint i, j, seqlen, noofreads, editstringlen, readerrcnt=0;//, buflen;
  int start=0, off=0;
  double lambda;
  double x=0;

  srand((unsigned int)time(NULL));
  sequence =  bl_getGeneSequence(space, reference, gene);
  seqlen = strlen(sequence);

  noofreads = ((Uint)(((double)seqlen+readlen)/(double)readlen))*cov;
  lambda = ((double)1.0/((double)readlen))*((double)cov);

  template = ALLOCMEMORY(space, NULL, char, readlen+1);
   
  buffer = ALLOCMEMORY(space, NULL, char, 2*readlen+1);
  quality = ALLOCMEMORY(space, NULL, char, 2*readlen+1);
  editstring = ALLOCMEMORY(space, NULL, char, 40*readlen+1);
 

  for(i=0; i < noofreads; i++) {
    //start=RANDINT(seqlen-1);

    x += 1.0/lambda;
    off = trunc(x);
    start = MAX(0, ((int)seqlen-1)-off);
    
    memset(template, 'A', sizeof(char)*readlen);
    template[readlen] = 0;
    
    for(j=0; j < readlen && start+j < seqlen; j++) {
      template[j] = sequence[start+j];
    }
   
    memset(buffer, 0, sizeof(char)*(2*readlen)+1);
    memset(quality, 0, sizeof(char)*(2*readlen)+1);
    memset(editstring, 0, sizeof(char)*(readlen)*40+1);
    editstringlen =0;
 
    //not used buflen = bl_fastxScramble(... 
    bl_fastxScramble(buffer, quality, template, readlen, acc, Pmis, Pins, 
        0, 0, 
        alphabet, alphabetsize, 
        minqual, maxqual, 
        editstring, &editstringlen, &readerrcnt);
    
    assert(strlen(buffer) == readlen);
    assert(strlen(quality) == readlen);

    fprintf(dev, "@%s:%d:%d;%s\n", gene->id, start, start+j, editstring);
    fprintf(dev, "%s\n+\n", buffer);
    fprintf(dev, "%s\n", quality);

  }

 
  return ;
}


/*----------------------- bl_simulateGeneSetSequencing -----------------------
 *    
 * @brief generate simulated reads for a set of genes
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_simulateGeneSetSequencing (void *space, FILE *dev, fasta_t *reference, geneset_t *genes, 
    Uint readlen, Uint cov, char *alphabet, Uint alphabetsize, Uint minqual, Uint maxqual,
    double acc, double Pmis, double Pins)
{
  Uint i;
 
  initProgressBarVT();

  for(i=0; i < genes->noofgenes; i++) { 
    progressBarVT("genes simulated.", genes->noofgenes, i, 25);

    bl_simulateGeneSequencing(space, dev, reference, &genes->genes[i], readlen, cov, 
        alphabet, alphabetsize, minqual, maxqual, acc, Pmis, Pins);
  }

  return ;
}



/*---------------------- bl_fastxPrintRandomBisulfiteReads ---------------------
 *    
 * @brief print a reference and a set of simulated bisulfite reads in fastq format
 *        to device dev with given conversion rates at cytosines
 * @author Christian Otto
 *   
 */
 
void
bl_fastxPrintRandomBisulfiteReads(
    FILE *reads, 
    char *sequence,
    char *sequencerates,
    Uint seqlen,
    Uint n,
    Uint minlen, 
    Uint maxlen, 
    char *alphabet, Uint alphabetsize,
    double acc, 
    double Pmis, double Pins, double Pdel, 
    unsigned char fastq,
    Uint minqual, Uint maxqual,
    char *five, Uint fivelen,
    char *three, Uint threelen, Uint polyAlen, char *prefix)
{

  Uint i, k, start, len, readerrcnt, fiveseqlen=fivelen, threeseqlen=threelen, 
  polyAseqlen=0, editstringlen=0;
  char *quality, *buffer, *editstring, *polyAseq, startchr = '>';
  
  if (fastq) {
    startchr = '@';
    assert(maxqual <= 126 && minqual >=33);
  }

  assert(seqlen > maxlen);
  assert(minlen <= maxlen);
  assert(alphabetsize > 1);
  assert(Pmis + Pins + Pdel == 1);

  srand((unsigned int)time(NULL));

  quality = ALLOCMEMORY(space, NULL, char, 2*(maxlen+fivelen+threelen+polyAlen));
  buffer = ALLOCMEMORY(space, NULL, char, 2*(maxlen+fivelen+threelen+polyAlen));
  editstring = ALLOCMEMORY(space, NULL, char, 40*(maxlen+fivelen+threelen+polyAlen));
  polyAseq = ALLOCMEMORY(space, NULL, char, 2*polyAlen);
  memset(polyAseq, 'A', 2*polyAlen);

  for(i=0; i < n; i++){
    memset(buffer, 0, sizeof(char)*(2*(maxlen+fivelen+threelen+polyAlen)));
    memset(quality, 0, sizeof(char)*(2*(maxlen+fivelen+threelen+polyAlen)));
    memset(editstring, 0, sizeof(char)*(40*(maxlen+fivelen+threelen+polyAlen)));
    editstringlen = 0;
    polyAseqlen = polyAlen;

    len = minlen + RANDINT(maxlen - minlen);
    start = RANDINT(seqlen-len);
   
    if(RANDUNIT > fiveacc) {
      if(RANDUNIT > .5) {
        fiveseqlen = fivelen - RANDINT(fivelen);
      }
    }

    if(RANDUNIT > threeacc) {
      if(RANDUNIT > .5) {
        threeseqlen = threelen - RANDINT(threelen);
      }
    }
    
    if(RANDUNIT > acc) {
      if(RANDUNIT > .5) {
        polyAseqlen = polyAlen - RANDINT(polyAlen);
      } else {
        polyAseqlen = polyAlen + RANDINT(polyAlen);
      }
    }
    if (prefix == NULL){
      fprintf(reads, "%c%d %d (len: %d) ", startchr, i, start+1, len);
    }
    else {
      fprintf(reads, "%c%d_%s %d_%s %d (len: %d) ", startchr, i, prefix, i, prefix, start+1, len);
    }
    readerrcnt=0;
    k=0; 
     
    k = bl_fastxScramble (buffer, quality, five, fiveseqlen, fiveacc, Pmis, Pins, 
        k, 0, alphabet, alphabetsize, minqual, maxqual, 
        editstring, &editstringlen, &readerrcnt);
    
    editstring[editstringlen++] = '|'; 

    k = bl_fastxBisulfiteScramble (buffer, quality, sequence, sequencerates, len, acc, Pmis, Pins, 
        k, start, alphabet, alphabetsize, minqual, maxqual, 
        editstring, &editstringlen, &readerrcnt);

    editstring[editstringlen++] = '|'; 

    k = bl_fastxScramble (buffer, quality, polyAseq, polyAseqlen, acc, Pmis, Pins, 
        k, 0, alphabet, alphabetsize, minqual, maxqual, 
        editstring, &editstringlen, &readerrcnt);

    editstring[editstringlen++] = '|'; 
    
    k = bl_fastxScramble (buffer, quality, three, threeseqlen, threeacc, Pmis, Pins, 
        k, 0, alphabet, alphabetsize, minqual, maxqual, 
        editstring, &editstringlen, &readerrcnt);
    
    fprintf(reads, "%s", editstring);
    fprintf(reads, " (r: %f)\n", (double) readerrcnt/(len+1.0));
    fprintf(reads, "%s\n", buffer);
    if (fastq) {
      fprintf(reads, "+\n");
      fprintf(reads, "%s\n", quality);
    }
  }

  FREEMEMORY(space, buffer);
  FREEMEMORY(space, quality);
  FREEMEMORY(space, polyAseq);
  FREEMEMORY(space, editstring);
  
  return ;
}

/*--------------------------- bl_fastxBisulfiteScramble ------------------------
 *    
 * @brief scramble a bisulfite-treated fasta sequence
 * @author Christian Otto
 *   
 */

Uint
bl_fastxBisulfiteScramble (char *buffer, char *quality, 
                           char *template, char *rates, Uint len, 
                           double acc, double Pmis, double Pins, 
                           Uint uoff, Uint voff, 
                           char *alphabet, Uint alphabetsize, 
                           Uint minqual, Uint maxqual,
                           char *editstring, Uint *editstringlen, Uint *readerrcnt)
{
  Uint u = uoff, k = (*editstringlen);
  int j;
  char errchr, chr;
  double errtype, rate;

  for (j=0; j < len; j++) {
    chr=template[voff+j]; 
    rate=(double)((int)rates[voff+j]-1)/100;
    /* 
     * RANDUNIT ==> [0,1)
     * rate == 0 ==> RANDUNIT always >= rate
     * rate == 1 ==> RANDUNIT always < rate
     */
    if (RANDUNIT >= rate){
      assert(chr == 'C' || chr == 'G');
      if (chr == 'C'){
        chr = 'T';
      }
      /* only +RC/-RC reads */
      if (chr == 'G'){
        chr = 'A';
      }
    }

    if (RANDUNIT > acc) { 
      errtype = RANDUNIT;
      (*readerrcnt) += 1;
      
      if (errtype <= Pmis) {
        errchr = chr;
        if (editstring) k += sprintf(&editstring[k], "%d:S;", j); 
        while(errchr == chr)
          errchr = alphabet[RANDINT(alphabetsize-1)];
        quality[u] = minqual + RANDINT(maxqual - minqual);
        buffer[u++] = errchr; 
      
      } else if (errtype <= Pmis + Pins) {
        errchr = alphabet[RANDINT(alphabetsize-1)];
        if (editstring) k += sprintf(&editstring[k], "%d:I;", j); 
        quality[u] = minqual + RANDINT(maxqual - minqual);
        buffer[u++] = errchr;
        if(j < len-1) { 
          quality[u] = minqual + RANDINT(maxqual - minqual);
          buffer[u++] = chr;
          j++;
        }
      
      } else {
        if (editstring) k += sprintf(&editstring[k], "%d:D;", j); 
        j--;
      }
    
    } else {
      quality[u] = minqual + RANDINT(maxqual - minqual);
      buffer[u++] = chr;
    }
  
  }

  (*editstringlen) = k;
  return u;
}
