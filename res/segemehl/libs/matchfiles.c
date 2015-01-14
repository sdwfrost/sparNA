
/*
 *  readmatchfiles.c
 *  read alignment files 
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 06.10.2010 01:02:08 CEST
 *  
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "alignment.h"
#include "debug.h"
#include "stringutils.h"
#include "basic-types.h"
#include "mathematics.h"
#include "sort.h"
#include "matfile.h"
#include "bitVector.h"
#include "info.h"
#include "zran.h"
#include "nw.h"
#include "matchfiles.h"
#include "evalmatchfiles.h"
#include "manout.h"
#include "matchfilesfields.h"
#include "matepairs.h"

/*-------------------------- bl_matchfileInitSplit ---------------------------
 *    
 * @brief initialize a split structure
 * @author Steve Hoffmann 
 *   
 */

void
bl_matchfileAddSplit (void *space, matchfileCross_t *cs, Uint v, 
    char strand, char edgetype, Uint edgechridx, Uint edge, Uint adjoint,
    Uint xstart, Uint xend, Uint row, Uint xno, Uint flags)
{

  matchfileSplit_t *split;
  cs[v].splits = ALLOCMEMORY(space, cs[v].splits, matchfileSplit_t, 
      cs[v].noofsplits+1);
 

  split = &cs[v].splits[cs[v].noofsplits];

  split->strand = strand;
  split->edgetype = edgetype;
  split->xno = xno;
  split->xstart = xstart;
  split->xend = xend;
  split->trans = 0;


  if(edgetype == 'A' && (flags & SPLIT_PREV_PLUS) && strand == '-') {
    split->trans = 1;
  }

  if(edgetype == 'A' && !(flags & SPLIT_PREV_PLUS) && strand == '+') { 
    split->trans = 1;
  }

  if(edgetype == 'D' && (flags & SPLIT_NEXT_PLUS) && strand == '-') {
    split->trans = 1;
  }

  if(edgetype == 'D' && !(flags & SPLIT_NEXT_PLUS) && strand == '+') {
    split->trans = 1;
  }
  
  split->edgechridx = edgechridx;
  split->edge = edge;
  split->adjoint = adjoint;
                
  cs[v].noofsplits++;
  
  return ;
}

/*---------------------- bl_matchfileCmpDeletionLength -----------------------
 *    
 * @brief a compare function for qsort to sort the deletions by length
 * a helper function for matchfilegapalign
 * @author Steve Hoffmann 
 *   
 */
 

Uint 
bl_matchfileCmpDeletionLength(Uint a, Uint b, void *tosort, void *info) {
  matchfileDeletion_t *arr;

  arr = (matchfileDeletion_t*) tosort;
  if(arr[a].len > arr[b].len) return 2;
  if(arr[b].len > arr[a].len) return 1;

  return 0;
}


/*-------------------------- bl_matchfileGapAdjust ---------------------------
 *    
 * @brief a helper function for gap align
 * @author Steve Hoffmann 
 *   
 */
 
char*
bl_matchfileGapAdjust(Alignment *al, Uint *len, unsigned char template) {
  Uint i = 0, j = 0, p = 0, q = 0, r = 0;
  char *string;

  string = malloc(sizeof(char)*(al->vlen+al->ulen+1));
  
  for(i=0; i < al->uoff; i++) {
    string[i] = al->u[i];
  }
 
  for(; i < al->voff; i++) {
    if(template) string[i] = al->v[i];
    else string[i] = '^';
  }

  
  r = i;
  for (i=0; i < al->numofmeops; i++) {
    //if Replacement occured
    if (al->meops[i].eop == Replacement) {
      for (j=0; j < al->meops[i].steps; j++) {
        string[j+r] = al->u[j+p+al->uoff];
      }
      p+=j;
      q+=j;
      r+=j;
    }

    if (al->meops[i].eop == Deletion) {
      //iter over all steps
      for (j=0; j < al->meops[i].steps; j++) {
        string[j+r] = '^';
      }
      //set ptrs
      r+=j;
      q+=j;
    }

    if (al->meops[i].eop == Insertion) {
      for (j=0; j < al->meops[i].steps; j++) {
        string[j+r] = al->u[j+p+al->uoff];
      }
      r+=j;
      p+=j;
    }
  }

  string[r]=0;

  *len = r;
  return string;
}


/*--------------------------- bl_matchfileGapAlign ---------------------------
 *    
 * @brief alignment of gaps
 * @author Steve Hoffmann 
 *   
 */
 

void
bl_matchfileGapAlign(matchfileDeletion_t *dels, Uint noofdels) {
  Uint i, *sidx, templatelen=0;
  int scores[] = {1,-1};
  int  *matrix;
  char *template, *temp;
  Alignment al;


  sidx = quickSort(NULL, dels, noofdels, bl_matchfileCmpDeletionLength, NULL);
  templatelen = dels[sidx[0]].len;
  template = ALLOCMEMORY(NULL, NULL, char, templatelen+1);
  memmove(template, dels[sidx[0]].string, templatelen);
  template[templatelen] = 0;

  /*progressive multiple alignment w/o guide tree. Generate template*/
  for(i=1; i < noofdels; i++) {
    
    matrix = nwmatrix(NULL, template, templatelen, 
        dels[sidx[i]].string, dels[sidx[i]].len, -1, constscr, scores);
    
    initAlignment(&al, template, templatelen, 0, 
        dels[sidx[i]].string, dels[sidx[i]].len, 0);

    nwtraceback(NULL, matrix, template, templatelen, dels[sidx[i]].string, 
        dels[sidx[i]].len, -1, constscr, scores, &al);
     
    temp = bl_matchfileGapAdjust(&al, &templatelen, 1);
    
    FREEMEMORY(space, template);
    template = temp;
    FREEMEMORY(space, matrix);
    wrapAlignment(&al);
  }

  /*progressive multiple alignment w/o guide tree. fit to template*/
  for(i=0; i < noofdels; i++) {
    matrix = nwmatrix(NULL, dels[sidx[i]].string, dels[sidx[i]].len,
        template, templatelen, -1, constscr, scores);
   
    initAlignment(&al, dels[sidx[i]].string, dels[sidx[i]].len, 0,
        template, templatelen, 0);

    nwtraceback(NULL, matrix, dels[sidx[i]].string, dels[sidx[i]].len, 
        template, templatelen, -1, constscr, scores, &al);
    
    temp = bl_matchfileGapAdjust(&al, &(dels[sidx[i]].len), 0);
    
    FREEMEMORY(space, dels[sidx[i]].string);
    dels[sidx[i]].string = temp;
    FREEMEMORY(space, matrix);
    wrapAlignment(&al);
  }

  FREEMEMORY(space, template);
  FREEMEMORY(space, sidx);
}

/*-------------------------- bl_matchfileInitIndex ---------------------------
 *    
 * @brief init the index
 * @author Steve Hoffmann 
 *   
 */

matchfileindex_t* bl_matchfileInitIndex(void *space) {
  
  matchfileindex_t *index;
  
  index = ALLOCMEMORY(space, NULL, matchfileindex_t, 1);
  index->exp = 15;
  index->noofbins = NULL;
  index->maxreadlen = 0;
  index->noofchroms = 0;
  index->bins = NULL;
  index->gzindex = NULL;
  index->chromnames = NULL;
  index->matchstart = NULL;
  index->matchend = NULL;
  index->stats = NULL;
  index->md5 = 0;
  index->submatrix = NULL;
  index->mean_coverage = 0;
  index->mean_qual = 0;
  index->Q_ERR = NULL;
  index->Q_N = NULL;
  index->P_ERR = NULL;

  return index;
}


/*------------------------ bl_matchfileDestructIndex -------------------------
 *    
 * @brief remove the index
 * @author Steve Hoffmann 
 *   
 */
 

void
bl_matchfileDestructIndex(void *space, matchfileindex_t *index) {
  Uint i;

  if(index->gzindex) { 
   if(index->gzindex->list) {
     FREEMEMORY(space, index->gzindex->list);
     index->gzindex->list = NULL;
   }
   FREEMEMORY(space, index->gzindex);
   index->gzindex = NULL;
  }

  for(i=0; i < index->noofchroms; i++) {
    FREEMEMORY(space, index->bins[i]);
    FREEMEMORY(space, index->chromnames[i]);
  }
  
  
  if(index->chromnames) FREEMEMORY(space, index->chromnames);
  if(index->noofbins) FREEMEMORY(space, index->noofbins);
  if(index->bins) FREEMEMORY(space, index->bins);
  
  if(index->Q_ERR) FREEMEMORY(space, index->Q_ERR);
  if(index->P_ERR) FREEMEMORY(space, index->P_ERR);
  if(index->Q_N) FREEMEMORY(space, index->Q_N);
  if(index->matchstart) FREEMEMORY(space, index->matchstart); 
  if(index->matchend) FREEMEMORY(space, index->matchend);
  if(index->submatrix) FREEMEMORY(space, index->submatrix);

  if(index->stats) {
    bl_matchfileDestructSampleStats(space, index->stats);
    FREEMEMORY(space, index->stats);
  }
}


/*--------------------- bl_matchfileGetChromIndexNumber ----------------------
 *    
 * @brief get the number of the chromosome chrname in the file index
 * @author Steve Hoffmann 
 *   
 */
 
Uint
bl_matchfileGetChromIndexNumber(matchfileindex_t *index, char *chromname) {
  Uint k, j, desclen;
  char *desc;

  for(k=0; k < index->noofchroms; k++) {
  
    desc = index->chromnames[k];
    desclen = strlen(index->chromnames[k]);

    if(strcmp(chromname, desc) == 0) {
      break;
    }

    for(j=0; j < desclen; j++) {
      if (isspace(desc[j])) break;
    }

    if(strlen(chromname) <= j && strncmp(chromname, desc, j) == 0) 
      break;
  }

  return k;
}


/*-------------------- bl_matchfileGetChromIndexNumberDB ---------------------
 *    
 * @brief get index number of chromosome in DB this is necessary
 *        because order in DB is not necessarily the same in index
 * @author Steve Hoffmann 
 *   
 */

  Uint
bl_matchfileGetChromIndexNumberDB (fasta_t *set, matchfileindex_t *index, 
    char *chr)
{
  Uint idxchrid=-1;
  Uint dbchrid=-1;

  idxchrid =  bl_matchfileGetChromIndexNumber(index, chr);
  assert(idxchrid != -1);
  if(set) { 
    for(dbchrid=0; dbchrid < set->noofseqs; dbchrid++) { 
      if(strcmp(bl_fastaGetDescription(set,dbchrid),
            index->chromnames[idxchrid])==0) {
        break;
      }
    }
    assert(dbchrid < set->noofseqs);
    return dbchrid;
  }

  return idxchrid;
}




/*--------------------------- bl_matchfileInitBin ----------------------------
 *    
 * @brief initialize a bin within the index for the interval [start,end] ond 
 * chrom no k
 * @author Steve Hoffmann 
 *   
 */

void 
bl_matchfileInitBin(void *space, matchfileindex_t *index, Uint k, Uint bin, 
    Uint start, Uint end) {
 
  Uint i;

  if(bin >= index->noofbins[k]) {

    index->bins[k] = ALLOCMEMORY(space, index->bins[k], 
        matchfileBin_t, bin+1);
    
    for(i=index->noofbins[k]; i <= bin; i++) {
      index->bins[k][i].end = end;
      index->bins[k][i].start = start;
      index->bins[k][i].offset = 0;
      index->bins[k][i].matches = 0;
      index->bins[k][i].endoff = 0;
    }
    index->noofbins[k] = bin + 1;
  }

  return;
}

/*------------------------ bl_matchfileIndexAddChrom -------------------------
 *    
 * @brief add a chromosome name to the index and initialize the bins
 * @author Steve Hoffmann 
 *   
 */
 

Uint
bl_matchfileIndexAddChrom(matchfileindex_t *index, char *chromname) {

  Uint k;

  k = bl_matchfileGetChromIndexNumber(index, chromname);
 
  if(k == index->noofchroms) {

    index->chromnames = 
      ALLOCMEMORY(space, index->chromnames, char*, index->noofchroms+1);
    index->chromnames[k] = ALLOCMEMORY(space, NULL, char, strlen(chromname)+1);
    memmove(index->chromnames[k], chromname, strlen(chromname));
    index->chromnames[k][strlen(chromname)] = 0;

    index->matchstart = 
      ALLOCMEMORY(space, index->matchstart, Uint, index->noofchroms+1);
    index->matchend = 
      ALLOCMEMORY(space, index->matchend, Uint, index->noofchroms+1);
    index->matchstart[k] = 0;
    index->matchend[k] = 0;
    
    index->noofbins =
      ALLOCMEMORY(space, index->noofbins, Uint, index->noofchroms+1);
    index->noofbins[k] = 0;

    index->bins =
      ALLOCMEMORY(space, index->bins, matchfileBin_t*, index->noofchroms+1);
    index->bins[k] = NULL;
    index->noofchroms++;

  }

  return k;
}

/*--------------------- bl_matchfileIndexAddAccessPoint ----------------------
 *    
 * @brief add an access point to the matchfileindex_t structure for 
 * given chromname, bin, start and end
 * @author Steve Hoffmann 
 *   
 */
 
void        
bl_matchfileIndexAddAccessPoint(void *space, matchfileindex_t *index, 
    char *chromname, Uint bin, Uint start, Uint end, off_t off, Uint matches, 
    Uint lastchrom, Uint lastbin, Uint lastend) {

  Uint k,i;

  k = bl_matchfileIndexAddChrom(index, chromname);

  index->matchstart[k] = (index->matchstart[k] > start) ? 
    start : index->matchstart[k];
  index->matchend[k] = (index->matchend[k] < end) ? 
    end : index->matchend[k];
  

  bl_matchfileInitBin(space, index, k, bin, start, end);

//  fprintf(stderr, "adding access point for %s:%d, %d-%d (bin:%d): %llu\n", 
//      chromname, k, start, end, bin, off);

  if(!index->bins[k][bin].offset) {
    index->bins[k][bin].offset = off;
  }

  if(index->bins[k][bin].start > start) {
    index->bins[k][bin].start = start;
  }

  if(index->bins[k][bin].end < end) {
    index->bins[k][bin].end = end;
  }


 // if(bin > 0) {
    //should never be 0
    if(off > 0){ 
      index->bins[lastchrom][lastbin].endoff = off-1;
    }
    index->bins[lastchrom][lastbin].matches = matches;
    index->bins[lastchrom][lastbin].end = lastend;
//    fprintf(stderr, "resetting lastchrom:%d, original matchend:%d, new matchend:%d\n", lastchrom, index->matchend[lastchrom], lastend);
    index->matchend[lastchrom] = (index->matchend[lastchrom] < lastend) ?
      lastend : index->matchend[lastchrom];
//  }

  if(lastchrom != k) { 
    lastbin = 0;
  } else { 
    lastbin +=1;
  }

  for(i=lastbin; i < bin; i++) {
    index->bins[k][i].offset = off;
  }
}


/*------------------- bl_matchfileIndexAccessPointOverlap --------------------
 *    
 * @brief update most distant access point and bin current read overlaps with 
 * @author Steve Hoffmann 
 *   
 */
 

void
bl_matchfileIndexAccessPointOverlap (void *space, 
    matchfileindex_t *index, char *chromname, Uint startbin, Uint endbin, 
    Uint start, Uint end, off_t ovlstart) { 
  
  Uint i, k;

  k = bl_matchfileIndexAddChrom(index, chromname);
  
  index->matchstart[k] = (index->matchstart[k] > start) ? 
    start : index->matchstart[k];
  index->matchend[k] = (index->matchend[k] < end) ? 
    end : index->matchend[k];
  
  bl_matchfileInitBin(space, index, k, endbin, start, end);

  for(i=startbin+1; i <= endbin; i++) {
    if (!index->bins[k][i].offset) {
      index->bins[k][i].offset = ovlstart;
    }
    if(index->bins[k][i].end < end) {
      index->bins[k][i].end = end;
    }
    if(index->bins[k][i].start > start) {
      index->bins[k][i].start = start;
    }
  }
}



/*---------------------------- bl_matchfileIndex -----------------------------
 *    
 * @brief build an index for matchfile_t
 * @author Steve Hoffmann 
 *   
 */


void bl_matchfileIndex (void *space, matchfile_t *file, fasta_t *set) {

  FILE *fp = NULL;
  stringset_t *fields = NULL;
  char *buffer = NULL, ch, first=1;
  char *curchrom = NULL, *filename, *diff, *read, *qual, *aln, *ref, strand=0, 
       *acceptorchr, *donorchr;
  Uint buffersize = 1024, len = 0, curstart = 0, 
       curend = 0, matches = 0, curbin = 0, lastchromidx = 0,
       startbin = 0, endbin = 0, lastbin = 0, k = 0,
    lastend = 0, p, u, s, q, id, allen=0,
       readlen=0, noofreads=0;
  unsigned char minqual = 255, maxqual = 0;

  matchfileindex_t *index;
  unsigned char mdcigarcheck = 1, header = 1;
  char *code;
  int gzlen;
  unsigned char gzip, fmt, getsubmatrix=1;
  struct gzidxfile *gzf = NULL;
  off_t off=0;
  Uint *Q_ERR;
  Uint *P_ERR;
  Uint *Q_N;

  code = getNTcodekey(space);
  filename = file->filename;
  gzip = file->gzip;
  fmt = file->fmt;
  index = bl_matchfileInitIndex(space);

  if (gzip) {
    index->gzindex = bl_zranGetIndex(filename, &gzlen);
    fp = fopen(filename, "rb");
    gzf = bl_initgzidxfile(fp, index->gzindex, 0, CHUNK);
  } else {
    fp = fopen(filename, "r");
  }

  if(fp == NULL) {
    DBGEXIT("Couldn't open file %s. Exit forced!\n", filename);
  }

  buffer = ALLOCMEMORY(space, NULL, char, buffersize);

  if(getsubmatrix) {
    index->submatrix = ALLOCMEMORY(space, NULL, double, 
        6*6*QRNGE*MAXREADLENGTH+1);
    memset(index->submatrix, 0, sizeof(double)*(6*6*QRNGE*MAXREADLENGTH+1));

    P_ERR = ALLOCMEMORY(space, NULL, Uint, MAXREADLENGTH);
    memset(P_ERR, 0, sizeof(Uint)*MAXREADLENGTH);
    Q_ERR = ALLOCMEMORY(space, NULL, Uint, QRNGE);
    memset(Q_ERR, 0, sizeof(Uint)*QRNGE);
    Q_N = ALLOCMEMORY(space, NULL, Uint, QRNGE);
    memset(Q_N, 0, sizeof(Uint)*QRNGE);
  }


  for(u=0; u < set->noofseqs; u++) { 
    bl_matchfileIndexAddChrom(index, bl_fastaGetDescription(set, u));
  }

  while((ch = (gzip) ? bl_getgzidxc(gzf) : getc(fp)) != EOF) {

    if(len == buffersize-1) {
      buffersize = 2*buffersize+1;
      buffer = ALLOCMEMORY(space, buffer, char, buffersize);
    }

    if(ch == '\n' && len > 0) {

      buffer = ALLOCMEMORY(space, buffer, char, len+1);  
      buffer[len] = '\0';
      header = (header) ? bl_matchfileIsHeader(buffer, len, fmt) : header;

      if (!header) { 
        fields = tokensToStringset(space, "\t", buffer, len);
        curstart = bl_matchfileGetStartPos(fields, fmt);
        curend = bl_matchfileGetEndPos(fields, fmt);
        
        if(curstart != curend && strand != -1) { 
          
          curchrom = bl_matchfileGetChrom(fields, fmt);
          strand = bl_matchfileGetStrand(fields, fmt);
          read = bl_matchfileGetRead(fields, fmt);
          qual = bl_matchfileGetQual(fields, fmt);
          diff = bl_matchfileGetDiffString(fields, fmt);
          aln  = bl_matchfileGetAln(fields, fmt);  
          acceptorchr = bl_matchfileGetNextChr(fields, fmt);
          donorchr = bl_matchfileGetPrevChr(fields, fmt);

          if(read) {
            readlen = strlen(read);
            if (readlen > index->maxreadlen) {
              index->maxreadlen = readlen;
            }
            noofreads++; 
          }

          if(curend+1 > 0) { 
            if(aln && curchrom) {

              if(acceptorchr) {
                id = bl_fastxFindIDIdx(acceptorchr, set);
                if (id == -1) {
                  DBGEXIT("reference sequence '%s' not found\n", acceptorchr);
                }
                acceptorchr = bl_fastaGetDescription(set, id);
                bl_matchfileIndexAddChrom(index, acceptorchr);
              }

              if(donorchr) {
                id = bl_fastxFindIDIdx(donorchr, set);
                if (id == -1) {
                  DBGEXIT("reference sequence '%s' not found\n", donorchr);
                }
                donorchr = bl_fastaGetDescription(set, id);
                bl_matchfileIndexAddChrom(index, donorchr);
              }
        
              id = bl_fastxFindIDIdx(curchrom, set);

              if (id == -1) {
                DBGEXIT("reference sequence '%s' not found\n", curchrom);
              }

              curchrom = bl_fastaGetDescription(set, id);
              ref = &bl_fastaGetSequence(set, id)[curstart-1];
              allen =strlen(aln);

              for(u=0, q=0, p=0; u < allen; u++) {
                
                //skip soft clipped part here               
                if(aln[u] == 'S') {
                  q++;
                  continue;
                }

                if(qual[q] < minqual) minqual = qual[q];
                if(qual[q] > maxqual) maxqual = qual[q];


                Q_N[(int)qual[q]-MINQUAL]++;

                if(aln[u] == 'M') {
                  MATRIX4D(index->submatrix, 6, QRNGE, MAXREADLENGTH,
                      (int)code[(int)read[q]], 
                      (int)code[(int)ref[p]], 
                      ((int)qual[q])-MINQUAL, q)+=1; 

                  if(read[q] != ref[p]) {            
                    Q_ERR[(int)qual[q]-MINQUAL]++;
                    if(strand == '+')
                      P_ERR[q]++;
                    else
                      P_ERR[readlen-q-1]++;


                    if(mdcigarcheck && ref[p] != diff[p]) {
                      NFO("MD doesnt match CIGAR in '%s'\n", buffer);
                      NFO("Further MD errors ignored in '%s' !\n", filename);
                      mdcigarcheck = 0;
                    }
                  } 
                } else {
                  Q_ERR[(int)qual[q]-MINQUAL]++;
                  if(strand == '-')
                    P_ERR[q]++;
                  else
                    P_ERR[readlen-q-1]++;
                }

                if(aln[u] == 'I') {
                  MATRIX4D(index->submatrix, 6, QRNGE, MAXREADLENGTH,
                      (int)code[(int)read[q]], 
                      (int)code[(int)'-'], 
                      ((int)qual[q])-MINQUAL, q)+=1; 
                }

                if(aln[u] == 'D') {
                  if(q > 0) s = 1; else s = 0;
                  MATRIX4D(index->submatrix, 6, QRNGE, MAXREADLENGTH,
                      (int)code[(int)'-'], 
                      (int)code[(int)ref[p]], 
                      MIN(((int)qual[q]),((int)qual[q-s]))-MINQUAL, q)+=1; 
                }

                if(aln[u] != 'I' ) {
                  p++;
                }
                if(aln[u] != 'D') {
                  q++;
                }
              }
              FREEMEMORY(space, aln);
            }

            if(diff) {    
              FREEMEMORY(space, diff);
            }

            if(curchrom) {

              k = bl_matchfileGetChromIndexNumber(index, curchrom);
              
              if(first) {
                lastchromidx = k;
                first = 0;
              }

              startbin = (curstart >> index->exp);
              endbin = (curend >> index->exp);

              if(k == index->noofchroms) {
                curbin = 0;
                startbin = 0;
                lastbin = 0;
              }
  
              if (k == index->noofchroms || !index->noofbins[k] || startbin > curbin || lastchromidx != k) {

#ifdef DBGIDX
                DBG("AP add: curchrom:%s curstart:%d, startbin:%d (exp:%d), off:%llu, lastchromidx:%d, lastbin:%d, lastend:%d\n", 
                    curchrom, curstart, startbin, index->exp, off, lastchromidx, lastbin, lastend);
#endif
                bl_matchfileIndexAddAccessPoint (space, index, curchrom, 
                    startbin, curstart, curend, off, matches, lastchromidx, lastbin, lastend);
                
                curbin = startbin;
                lastbin = startbin;
                lastchromidx = k;
                lastend = curend;
                matches = 0;
              }

              if(curend > lastend) {
                lastend = curend;
              }

              if (endbin > curbin) {

                bl_matchfileIndexAccessPointOverlap (space, index, curchrom, 
                    startbin, endbin, curstart, curend, off);
              }

              if(gzip) {
                off = bl_ftellgzidx(gzf);
              } else {
                off = ftello(fp);
                if (off == -1) {
                  DBGEXIT("ftello for '%s' failed. Exit forced.\n", filename);
                }
              }

            }
            matches++;
          }
        }
        destructStringset(space, fields);
      }

      buffer = ALLOCMEMORY(space, buffer, char, buffersize);
      len = 0;

    } else {
      if(ch != '\n') buffer[len++] = ch;
    }
  }


  if(matches && index->noofbins[k] == 0) { 

    bl_matchfileIndexAddAccessPoint (space, index, curchrom, 
        startbin, curstart, curend, off, matches, lastchromidx, 
        lastbin, lastend);
    
    lastend = curend;
  }


  if(gzip) {
    off = bl_ftellgzidx(gzf);
  } else {
    off = ftello(fp);
    if (off == -1) {
      MSG("ftello failed. Exit forced.\n");
      exit(-1);
    }
  }

  if(matches) {
    index->bins[k][curbin].endoff = off;
    index->bins[k][curbin].matches = matches;
    index->bins[k][curbin].end = curend;
    index->matchend[k] = (curend > index->matchend[k]) ? 
      curend : index->matchend[k];
#ifdef DBGIDX
    DBG("setting matchend[%d]=%d (off:%llu) after %d matches, curend %d\n", 
        k, index->matchend[k], off, matches, curend);
#endif
  }

  if(gzip) {
    bl_destructgzidxfile(gzf);
    FREEMEMORY(space, gzf);
  }

  file->index = index;
  FREEMEMORY(space, buffer);
  FREEMEMORY(space, code);
  fclose(fp);

  file->index->Q_ERR = Q_ERR;
  file->index->P_ERR = P_ERR;
  file->index->Q_N = Q_N;
  file->index->noofreads = noofreads;
  file->index->minqual = minqual;
  file->index->maxqual = maxqual;

  return;
}



/*-------------------------- bl_matchfileNextRowGet --------------------------
 *    
 * @brief assign a row to the read for visual representation on screen
 * @author Steve Hoffmann 
 *   
 */

Uint
bl_matchfileNextRowGet(void *space, matchfileCross_t *cs, Uint v, Uint len) {
  Uint i, size1=0, size2=0, max1=0, max2=0, max;
  bitvector a;

  if(v-1 < len) {
    size1 = cs[v-1].len;
    max1 = cs[v-1].maxrow;
  }

  size2 = cs[v].len;
  max2 = cs[v].maxrow;
  max = MAX(max1, max2);

  a = initbitvector(space, max+2);
  for(i=0; i < MAX(size1, size2); i++) {

    if (i < size1) {
      bitvector_setbit(a, cs[v-1].row[i], 1);
    }
    if (i < size2) {  
      bitvector_setbit(a, cs[v].row[i], 1);
    }
  }

  for(i=0; i < max+2; i++) {
    if (bitvector_getbit(a,i) == 0) {
      break;
    }
  }

  FREEMEMORY(space, a);
  return i;
}


/*----------------------------- bl_matchfileRead -----------------------------
 *    
 * @brief read all matches from start to end on chromname
 * @author Steve Hoffmann 
 *   
 */

matchfileCross_t*
bl_matchfileRead(void *space, matchfile_t *file, char *chromname, 
    Uint start, Uint end, Uint maxcover, fasta_t *set, unsigned char fields, 
    matchfileCross_t *input) {

  stringset_t *token;

  Uint buffersize=1024, row=0, startbin, //endbin, 
       startsite, endsite, len=0, i=0, u=0, v=0, k=0,  
//       curcnt=0,
//       bisulfite=0, 
       curalnlen, 
       dellen = 0, 
//       curstart=0, curend=0, 
       startcover=0, 
       curcover=0, 
//       acceptorpos=0, donorpos=0, xstart=0, xend = 0, xno = 0, 
       acceptorchridx, donorchridx, adjoint, 
//       acceptorflg = 0, donorflg = 0, pnext = 0, 
       rnextidx=0, 
       curchromidx; 
  char *buffer = NULL, *delstring=NULL, *delquals=NULL, ch, 
       //*curseq=NULL, *curqual=NULL, *curaln, 
       *filename//, strand
//       ,*acceptorchr = NULL, *donorchr = NULL, *rnext, *curchrom
       ;
  unsigned char header = 1;
  int ret=0, // edist=0, 
      readlen=0;
  matchfileCross_t *cs = NULL;
  matchfileindex_t *index;
  matchfileRec_t *r;
  unsigned char startsiteset=0, curedist=0, gzip, fmt;
  char eop;
  FILE *fp = NULL;
  struct gzidxfile *gzf = NULL;

#ifdef DBGIDX
  Uint acceptorcnt = 0, donorcnt = 0;
#endif

  gzip = file->gzip;
  fmt = file->fmt;
  index = file->index;
  filename = file->filename;
   
  r = ALLOCMEMORY(space, NULL, matchfileRec_t, 1);
  memset(r, 0, sizeof(matchfileRec_t));

  

  startbin = (start >> index->exp);
  k = bl_matchfileGetChromIndexNumber(index, chromname); 

  if(!input) { 
    cs = ALLOCMEMORY(space, NULL, matchfileCross_t, (end-start+1));
    memset(cs, 0, sizeof(matchfileCross_t)*(end-start+1)); 
  } else {
    cs = input;
  }


  if(k >= index->noofchroms || startbin >= index->noofbins[k]) { 
    return cs;
  }

  if (gzip) {  
    fp = fopen(filename, "rb");
    gzf = bl_initgzidxfile(fp, index->gzindex, index->bins[k][startbin].offset, MEDIUMCHUNK);
  } else {
    fp = fopen(filename, "r");
    ret = fseeko (fp, index->bins[k][startbin].offset, SEEK_SET);
 
    if(ret == -1) {
      DBGEXIT("fseeko failed for '%s' Exit forced!\n", filename);
    }
  }

  if(fp == NULL) {
    DBGEXIT("Couldn't open file %s. Exit forced!\n", filename);
  }

  curchromidx = k;
  r->curend = start;
  
  buffer = ALLOCMEMORY(space, NULL, char, buffersize);

  while((ch = (gzip) ? bl_getgzidxc(gzf) : getc(fp)) != EOF && 
       r->curstart <= end && curchromidx == k) {

 
    if(len == buffersize-1) {
      buffersize = 2*buffersize+1;
      buffer = ALLOCMEMORY(space, buffer, char, buffersize);
    }

    if(ch == '\n' && len > 0) {

      buffer = ALLOCMEMORY(space, buffer, char, len+1); 
      buffer[len] = '\0';

#ifdef DBGIDX
      DBG("buffer: %s\n", buffer);
#endif

      if(header) header = bl_matchfileIsHeader(buffer, len, fmt);

      if(!header) { 
        token = tokensToStringset(space, "\t", buffer, len);

#ifdef NEWFUNCT        
        curchrom = bl_matchfileGetChrom(token, fmt);
        curstart = bl_matchfileGetStartPos(token, fmt);
        curend   = bl_matchfileGetEndPos(token, fmt);
#endif

        r = bl_matchfileGetMatchFileRec(r, fields, token, fmt);
        curchromidx =  bl_matchfileGetChromIndexNumber(index, r->curchrom); 


        //fprintf(stderr, "reading %s: %d - %d\n", curchrom, curstart, curend);
        /*last condition to avoid inclusion of 0-alignments in BAM files*/
        if (r->curstart != r->curend && 
            r->curend >= start && r->curstart <= end && 
            r->curend+1 > 0 && curchromidx == k) {
        
          readlen  = strlen(r->curseq);
#ifdef NEWFUNCT              
          curseq   = bl_matchfileGetRead(token, fmt);


          if(fields & MFREAD_QUAL) 
            curqual  = bl_matchfileGetQual(token, fmt);
          
          if(fields & MFREAD_MCNT)
            curcnt   = bl_matchfileGetMatchCnt(token, fmt);

	  // if (fields & MFREAD_BISULFITE)
	  bisulfite = bl_matchfileGetBisulfite(token, fmt);

          curaln   = bl_matchfileGetAln(token, fmt);
          edist    = bl_matchfileGetEdist(token, fmt);
          strand   = bl_matchfileGetStrand(token, fmt);
          rnext    = bl_matchfileGetRNext(token, fmt);
          pnext    = bl_matchfileGetPNext(token, fmt);
   
          if(fields & MFREAD_SPLITS) { 
            acceptorpos = bl_matchfileGetNextPos(token, fmt);
            acceptorchr = bl_matchfileGetNextChr(token, fmt);
            acceptorflg = bl_matchfileGetNextFlag(token, fmt);
            donorpos = bl_matchfileGetPrevPos(token, fmt);
            donorchr = bl_matchfileGetPrevChr(token, fmt);
            donorflg = bl_matchfileGetPrevFlag(token, fmt);
            xstart = bl_matchfileGetSplitStart(token, fmt);
            xend = bl_matchfileGetSplitEnd(token, fmt);
            xno = bl_matchfileGetSplitNumber(token, fmt);
          }
#endif
          if(r->edist > 255) {
            curedist = 255; 
          } else curedist = r->edist;



          curalnlen = strlen(r->curaln);
          v = r->curstart - start;
          curcover = (v < end-start+1) ? cs[v].len : startcover;

          if(curcover < maxcover) { 
            startsite = 0;
            endsite = curalnlen-1;

            if (r->curaln[startsite] != 'M' || r->curaln[endsite] != 'M') {
              startsiteset = 0;

              for(i=0; i < curalnlen; i++) {
                eop = r->curaln[i];
                if(!startsiteset && eop == 'M') {
                  startsite = i;
                  startsiteset = 1;
                } 
                if(eop == 'M') {
                  endsite = i;
                }
              }
            }

            if(fields & (MFREAD_ROW | MFREAD_SPLITS)) { 
              if(v < end-start+1) {
                row = bl_matchfileNextRowGet(space, cs, v, end-start+1);
              } else {
                startcover++;
                row = bl_matchfileNextRowGet(space, cs, 0, end-start+1);
              }
            }

            donorchridx = -1;
            acceptorchridx = -1;

            if(r->donorchr) {
#ifdef DBGIDX
              DBG("setting donor chr:%d\n", r->donorchr);
              if(i==0) donorcnt++;
#endif
              donorchridx=bl_matchfileGetChromIndexNumberDB(set, index, r->donorchr);
            }

            if(r->acceptorchr) { 
#ifdef DBGIDX
              DBG("setting acceptor %d\n", r->acceptorchr);
              if(i==0) acceptorcnt++;
#endif
              acceptorchridx=bl_matchfileGetChromIndexNumberDB(set, index, r->acceptorchr);
            }

            if(r->rnext && r->rnext[0] == '=') { 
              rnextidx = k; 
            } else if(r->rnext && r->rnext[0] != '*') { 
              rnextidx = bl_matchfileGetChromIndexNumberDB(set, index, r->rnext);
            }

            
            for(i=0, u=0; i < curalnlen; i++) {

              if(r->curaln[i] == 'S') {
                u++;
                continue;
              }

              if(fields & MFREAD_SPLITS) { 
                if(r->donorchr && v < end-start+1 && 
                    ((i == 0 && r->strand == '+') || 
                      (i == curalnlen-1 && r->strand == '-'))) {
                  adjoint = -1;

                  if (r->acceptorchr) {
                    if (r->strand == '+')
                      adjoint = r->curend;
                    else
                      adjoint = r->curstart;
                  }

                  bl_matchfileAddSplit(space, cs, v, r->strand, 'A',
                      donorchridx, r->donorpos, adjoint, 
                      r->xstart, r->xend, row, r->xno, r->donorflg);
                }

                if(r->acceptorchr && v < end-start+1 && 
                    (( i == 0 && r->strand == '-') || 
                      (i == curalnlen-1 && r->strand == '+'))) {
                  adjoint = -1;

                  if (r->donorchr) {
                    if(r->strand == '+')
                      adjoint = r->curstart;
                    else
                      adjoint = r->curend;
                  }

                  bl_matchfileAddSplit(space, cs, v, r->strand, 'D',
                      acceptorchridx, r->acceptorpos, adjoint, 
                      r->xstart, r->xend, row, r->xno, r->acceptorflg);

                }
              }
 
              eop = r->curaln[i];


              if(startsite == i && v < end-start+1) {
                cs[v].starts++;        
                if(r->rnext && r->rnext[0] != '*') {
                  bl_matchfileAddMate (space, &cs[v], rnextidx, r->pnext);
                }
              }

              if(endsite == i && v < end-start+1) {
                cs[v].ends++;        
              }

              if (eop != 'I') {

                if (delstring && fields & MFREAD_DELS) {
                  if (v < end-start+1) {
                    cs[v].dels = ALLOCMEMORY(space, cs[v].dels, 
                        matchfileDeletion_t, cs[v].noofdels+1);
                    cs[v].dels[cs[v].noofdels].len = dellen;
                    cs[v].dels[cs[v].noofdels].string = delstring;

                    if(fields & MFREAD_QUAL)
                    cs[v].dels[cs[v].noofdels].quals = delquals;
                    
                    if(fields & MFREAD_RPOS)
                    cs[v].dels[cs[v].noofdels].readpos = u - dellen;

                    if (fields & MFREAD_MCNT)
                    cs[v].dels[cs[v].noofdels].matchcnt = r->curcnt;

		    //if (fields & MFREAD_BISULFITE)
		    cs[v].dels[cs[v].noofdels].bisulfite = r->bisulfite;
                    
                    if (fields & MFREAD_ROW)
                    cs[v].dels[cs[v].noofdels].row = row;

                    cs[v].dels[cs[v].noofdels].edist = curedist;
                    cs[v].noofdels++;
                  } else {
                    
                    FREEMEMORY(space, delstring);
                    if(fields & MFREAD_QUAL)
                    FREEMEMORY(space, delquals);
                  }
                  delstring = NULL;
                  delquals = NULL;
                  dellen = 0;
                } 

                if(v < end-start+1) {
                  cs[v].len++;
                  
                  cs[v].chars = ALLOCMEMORY(space, cs[v].chars, 
                      char, cs[v].len+1);
                  
                  if(fields & MFREAD_FEAT) { 
                    cs[v].feat = ALLOCMEMORY(space, cs[v].feat, 
                        char, cs[v].len+1);
                    if(i==startsite) {
                      cs[v].feat[cs[v].len-1] = '*';
                    } else if(i==endsite) { 
                      cs[v].feat[cs[v].len-1] = '$';
                    } else { 
                      cs[v].feat[cs[v].len-1] = 0;
                    }
                  }
                
                  if(fields & MFREAD_QUAL) 
                  cs[v].quals = ALLOCMEMORY(space, cs[v].quals, 
                      char, cs[v].len+1);
                  
                  cs[v].strands = ALLOCMEMORY(space, cs[v].strands, 
                      char, cs[v].len+1);
                  
                  if(fields & MFREAD_RPOS) 
                  cs[v].readpos = ALLOCMEMORY(space, cs[v].readpos, 
                      uint32_t, cs[v].len+1);

                  if(fields & MFREAD_RLEN)
                  cs[v].readlen = ALLOCMEMORY(space, cs[v].readlen, 
                      uint32_t, cs[v].len+1);
                  
                  if(fields & MFREAD_MCNT)
                  cs[v].matchcnt = ALLOCMEMORY(space, cs[v].matchcnt, 
                      uint32_t, cs[v].len+1);

		          //if (fields & MFREAD_BISULFITE)
		         cs[v].bisulfite = ALLOCMEMORY(space, cs[v].bisulfite,
		             uint32_t, cs[v].len+1);
						
                  
                  if(fields & MFREAD_ROW)
                  cs[v].row = ALLOCMEMORY(space, cs[v].row, 
                      uint32_t, cs[v].len+1);
                  
                  cs[v].edist = ALLOCMEMORY(space, cs[v].edist, 
                      unsigned char, cs[v].len+1);
                }

                if (eop != 'D') {
                  if(v < end-start+1) {

                    cs[v].chars[cs[v].len-1] = r->curseq[u]; 
                    cs[v].chars[cs[v].len] = 0 ;
                    
                    
                    if(fields & MFREAD_QUAL) { 
                        cs[v].quals[cs[v].len-1] = r->curqual[u]; 
                        cs[v].quals[cs[v].len] = 0;
                    }
                    
                    if(fields & MFREAD_ROW) { 
                        cs[v].row[cs[v].len-1] = row;
                        if(cs[v].maxrow < row) {
                            cs[v].maxrow = row;
                        }
                    }
                  
                    if(fields & MFREAD_RLEN)
                    cs[v].readlen[cs[v].len-1] = readlen;
                    
                    cs[v].edist[cs[v].len-1] = curedist;
                    cs[v].strands[cs[v].len-1] = r->strand;

                  }
                  u++;
                } else { 
                  if(v < end-start+1) {
                    cs[v].chars[cs[v].len-1] = '-';
                    cs[v].chars[cs[v].len] = 0 ;

                    if(fields & MFREAD_QUAL) {  
                      cs[v].quals[cs[v].len-1] = r->curqual[u];
                      cs[v].quals[cs[v].len] = 0;
                    }

                    if(fields & MFREAD_ROW) {
                      cs[v].row[cs[v].len-1] = row;
                      if(cs[v].maxrow < row) {
                        cs[v].maxrow = row;
                      }
                    }


                    if(fields & MFREAD_RLEN)
                    cs[v].readlen[cs[v].len-1] = readlen;

                    cs[v].edist[cs[v].len-1] = curedist;
                    cs[v].strands[cs[v].len-1] = r->strand;

                  }
                }

                if(v < end-start+1) {

                  if(fields & MFREAD_RPOS)
                  cs[v].readpos[cs[v].len-1] = u-1;
                  
                  if(fields & MFREAD_RLEN)
                  cs[v].readlen[cs[v].len-1] = readlen;

                  if(fields & MFREAD_MCNT)
                  cs[v].matchcnt[cs[v].len-1] = r->curcnt;

		  //if (fields & MFREAD_BISULFITE)
		  cs[v].bisulfite[cs[v].len-1] = r->bisulfite;

                }
                v++;
              } else {
                if(fields & MFREAD_DELS) { 
                  delstring = ALLOCMEMORY(space, delstring, char, dellen+2);
                  if(fields & MFREAD_QUAL) { 
                    delquals = ALLOCMEMORY(space, delquals, char, dellen+2);
                    delquals[dellen] = r->curqual[u];
                    delquals[dellen+1] = 0;
                  }
                  delstring[dellen] = r->curseq[u];
                  delstring[dellen+1] = 0;

                  dellen++;
                  u++;
                }
              }
            }

            if(delstring) {
              FREEMEMORY(space, delstring);
              if(fields & MFREAD_QUAL)
              FREEMEMORY(space, delquals);
              dellen = 0;
              delstring = NULL;
              delquals = NULL;
            }
          }          
        } 
        FREEMEMORY(space, r->curaln);
        FREEMEMORY(space, r->diff);

        destructStringset(space, token);
      }
      buffer = ALLOCMEMORY(space, buffer, char, buffersize);
      len = 0;
    } else {
      if(ch != '\n') buffer[len++] = ch;
    }
  }
  
  fclose(fp);

  if(gzip) {
    bl_destructgzidxfile(gzf);
    FREEMEMORY(space, gzf);
  }
 
  FREEMEMORY(space, r);
  FREEMEMORY(space, buffer);
  return cs;
}

/*-------------------------- bl_matchfileIndexShow ---------------------------
 *    
 * @brief dump the index for debugging
 * @author Steve Hoffmann 
 *   
 */


void
bl_matchfileIndexShow(void *space, matchfileindex_t *index, char *filename) {
  Uint i, k;
  stringset_t *fields;
  char *line;
  FILE *fp;
  Uint curend; //, curstart;

  line = ALLOCMEMORY(space, NULL, char, 1000);
  fp = fopen(filename, "r");

  for(k=0; k < index->noofchroms; k++) {

    for(i=0; i < index->noofbins[k]; i++) {
      fseeko (fp, index->bins[k][i].offset, SEEK_SET);
      line = fgets(line, 1000, fp);


      if(index->bins[k][i].offset) {
        fields = tokensToStringset(space, "\t", line, strlen(line));
        //not used: curstart = bl_matchfileGetStartPos(fields, 0);
        curend   = bl_matchfileGetEndPos(fields, 0);

        assert(curend >= i *(1 << index->exp));

        destructStringset(space, fields);
      } else {

        DBG( "binstart:%d start:%d, end:%d, matches:%d, offset:%lu\n",
            i*(1 << index->exp),
            index->bins[k][i].start, index->bins[k][i].end, 
            index->bins[k][i].matches, index->bins[k][i].offset);

      }
    }
  }
}




/*-------------------------- bl_matchfileCopyCross ---------------------------
 *    
 * @brief copy the x-section
 * @author Steve Hoffmann 
 *   
 */
 

matchfileCross_t* 
bl_matchfileCopyCross(void *space, matchfileCross_t *xs, matchfileCross_t *cs) {
  
  uint32_t len;
  Uint j;
 
  memmove(xs, cs, sizeof(matchfileCross_t));
  len = cs->len;
  
  xs->noofsplits = cs->noofsplits;
  xs->noofdels = cs->noofdels;
  xs->chars = ALLOCMEMORY(space, NULL, char, len+1);
  memset(xs->chars, 0, sizeof(char)*(len+1));
  xs->feat = ALLOCMEMORY(space, NULL, char, len+1);
  memset(xs->feat, 0, sizeof(char)*(len+1));
  xs->quals = ALLOCMEMORY(space, NULL, char, len+1);
  memset(xs->quals, 0, sizeof(char)*(len+1));
  xs->strands = ALLOCMEMORY(space, NULL, char, len+1);
  memset(xs->strands, 0, sizeof(char)*(len+1));
 
  xs->bisulfite = ALLOCMEMORY(space, NULL, uint32_t, len);
  xs->readpos = ALLOCMEMORY(space, NULL, uint32_t, len);
  xs->readlen = ALLOCMEMORY(space, NULL, uint32_t, len);
  xs->matchcnt = ALLOCMEMORY(space, NULL, uint32_t, len);
  xs->row = ALLOCMEMORY(space, NULL, uint32_t, len);
  xs->edist = ALLOCMEMORY(space, NULL, unsigned char, len);


  memmove(xs->chars, cs->chars, sizeof(char)*len);
  memmove(xs->feat, cs->feat, sizeof(char)*len);
  memmove(xs->quals, cs->quals, sizeof(char)*len);
  memmove(xs->strands, cs->strands, sizeof(char)*len);
  memmove(xs->bisulfite, cs->bisulfite, sizeof(uint32_t)*len);
  memmove(xs->readpos, cs->readpos, sizeof(uint32_t)*len); 
  memmove(xs->readlen, cs->readlen, sizeof(uint32_t)*len);
  memmove(xs->matchcnt, cs->matchcnt, sizeof(uint32_t)*len);
  memmove(xs->row, cs->row, sizeof(uint32_t)*len);
  memmove(xs->edist, cs->edist, sizeof(unsigned char)*len);

  if(cs->dels) {
    xs->dels = ALLOCMEMORY(space, NULL, matchfileDeletion_t, cs->noofdels);
    memmove(xs->dels, cs->dels, sizeof(matchfileDeletion_t)*cs->noofdels);
    
    for(j=0; j < cs->noofdels; j++) {
      xs->dels[j].string = 
        ALLOCMEMORY(space, NULL, char, strlen(cs->dels[j].string)+1);
      memset(xs->dels[j].string, 0, strlen(cs->dels[j].string)+1);
      xs->dels[j].quals = 
        ALLOCMEMORY(space, NULL, char, strlen(cs->dels[j].quals)+1);
      memset(xs->dels[j].quals, 0, strlen(cs->dels[j].quals)+1);
      memmove(xs->dels[j].string, cs->dels[j].string,
          strlen(cs->dels[j].string));
      memmove(xs->dels[j].quals, cs->dels[j].string,
          strlen(cs->dels[j].quals));
    }
  }

  if(cs->matelinks) {
    xs->matelinks = ALLOCMEMORY(space, NULL, matelink_t, cs->noofmatelinks);
    memmove(xs->matelinks, cs->matelinks, sizeof(matelink_t)*cs->noofmatelinks);
    xs->noofmatelinks = cs->noofmatelinks;
  }

  if(cs->splits) {
    xs->splits = ALLOCMEMORY(space, NULL, matchfileSplit_t, cs->noofsplits);
    memmove(xs->splits, cs->splits, sizeof(matchfileSplit_t)*cs->noofsplits);
  }

  return xs;
}



/*------------------------ bl_matchfileDestructCross -------------------------
 *    
 * @brief remove the x-section from the heap
 * @author Steve Hoffmann 
 *   
 */

void
bl_matchfileDestructCross(void *space, matchfileCross_t *cs, Uint len) {
  Uint i,j;

  for(i=0; i < len; i++) {

    if(cs[i].chars) FREEMEMORY(space, cs[i].chars);
    if(cs[i].feat) FREEMEMORY(space, cs[i].feat);
    if(cs[i].quals) FREEMEMORY(space, cs[i].quals);
    if(cs[i].strands) FREEMEMORY(space, cs[i].strands);
    if(cs[i].readpos) FREEMEMORY(space, cs[i].readpos);
    if(cs[i].matchcnt) FREEMEMORY(space, cs[i].matchcnt);
    if(cs[i].bisulfite) FREEMEMORY(space, cs[i].bisulfite);
    if(cs[i].readlen) FREEMEMORY(space, cs[i].readlen);
    if(cs[i].row) FREEMEMORY(space, cs[i].row);
    if(cs[i].edist) FREEMEMORY(space, cs[i].edist);
    if(cs[i].matelinks) FREEMEMORY(space, cs[i].matelinks);


    if(cs[i].dels) {
      for(j=0; j < cs[i].noofdels; j++) {
        FREEMEMORY(space, cs[i].dels[j].string);
        FREEMEMORY(space, cs[i].dels[j].quals);
      }
      FREEMEMORY(space, cs[i].dels);
    }

    if(cs[i].splits) {
      FREEMEMORY(space, cs[i].splits);
    }


 //   memset(&cs[i], 0, sizeof(matchfileCross_t));
  }
}


/*-------------------------- bl_matchfileReadIndex ---------------------------
 *    
 * @brief read the index (including stats if avail) from disk
 * @author Steve Hoffmann 
 *   
 */
 
matchfileindex_t * bl_matchfileReadIndex (void *space, char *filename)
{
 
  FILE *fp;
  unsigned char flags = 0;
  struct point *list;
  Uint i, j, len;
  matchfileindex_t *idx;

  idx = bl_matchfileInitIndex(space);
  
  fp = fopen(filename, "r");
  if (fp == NULL) {
    DBG("Couldn't open file '%s'. Exit forced.\n", filename);
    exit(-1);
  }

  fread(&flags, sizeof(unsigned char), 1, fp);

  if(flags & GZ_IDX_STORED) {
    idx->gzindex = ALLOCMEMORY(space, NULL, struct access, 1);
    fread(&idx->gzindex->have, sizeof(int), 1, fp);
    fread(&idx->gzindex->size, sizeof(int), 1, fp);
    
    list = ALLOCMEMORY(space, NULL, struct point, idx->gzindex->size);

    for(i=0; i < idx->gzindex->size; i++) {
      fread(&list[i].out, sizeof(off_t), 1, fp);
      fread(&list[i].in, sizeof(off_t), 1, fp);
      fread(&list[i].bits, sizeof(int), 1, fp);
      fread(&list[i].window, sizeof(char), WINSIZE, fp);
    }

    idx->gzindex->list = list;
  }

  fread(&idx->exp, sizeof(Uint), 1, fp);
  fread(&idx->noofchroms, sizeof(Uint), 1, fp);
  fread(&idx->noofreads, sizeof(Uint), 1, fp);


  idx->chromnames = ALLOCMEMORY(space, NULL, char*, idx->noofchroms);
  
  for(i=0; i < idx->noofchroms; i++) {
    fread(&len, sizeof(Uint), 1, fp);
    idx->chromnames[i] = ALLOCMEMORY(space, NULL, char, len+1);
    memset(idx->chromnames[i], 0, len+1);
    fread(idx->chromnames[i], sizeof(char), len, fp); 
  }

  idx->matchstart = ALLOCMEMORY(space, NULL, Uint, idx->noofchroms);
  fread(idx->matchstart, sizeof(Uint), idx->noofchroms, fp);
  idx->matchend = ALLOCMEMORY(space, NULL, Uint, idx->noofchroms);
  fread(idx->matchend, sizeof(Uint), idx->noofchroms, fp);
  idx->noofbins= ALLOCMEMORY(space, NULL, Uint, idx->noofchroms);
  fread(idx->noofbins, sizeof(Uint), idx->noofchroms, fp);

  idx->bins = ALLOCMEMORY(space, NULL, matchfileBin_t*, idx->noofchroms);

  for(i=0; i < idx->noofchroms; i++) {
    idx->bins[i] = NULL;
    if(idx->noofbins[i]) { 
      idx->bins[i]=ALLOCMEMORY(space, NULL, matchfileBin_t, idx->noofbins[i]);
      for(j=0; j < idx->noofbins[i]; j++) {
        fread(&idx->bins[i][j].start, sizeof(Uint), 1, fp);
        fread(&idx->bins[i][j].end, sizeof(Uint), 1, fp);
        fread(&idx->bins[i][j].matches, sizeof(Uint), 1, fp);
        fread(&idx->bins[i][j].offset, sizeof(off_t), 1, fp);
        fread(&idx->bins[i][j].endoff, sizeof(off_t), 1, fp);
      }
    }
  }

  fread(&idx->maxreadlen, sizeof(Uint), 1, fp);
  idx->submatrix = ALLOCMEMORY(space, NULL, double, 6*6*QRNGE*MAXREADLENGTH);
  fread(idx->submatrix, sizeof(double), 6*6*QRNGE*MAXREADLENGTH, fp);
  fread(&idx->mean_coverage, sizeof(double), 1, fp);
  fread(&idx->mean_qual, sizeof(double), 1, fp);
  
  idx->P_ERR = ALLOCMEMORY(space, NULL, Uint, MAXREADLENGTH);
  idx->Q_ERR = ALLOCMEMORY(space, NULL, Uint, QRNGE);
  idx->Q_N = ALLOCMEMORY(space, NULL, Uint, QRNGE);

  fread(idx->P_ERR, sizeof(Uint), MAXREADLENGTH, fp);
  fread(idx->Q_ERR, sizeof(Uint), QRNGE, fp);
  fread(idx->Q_N, sizeof(Uint), QRNGE, fp);

  if(flags & STATS_TAB_STORED) {
     
    idx->stats = ALLOCMEMORY(space, NULL, matchfileSampleStats_t, 1);
    memset(idx->stats, 0, sizeof(matchfileSampleStats_t));

    fread(&idx->stats->n, sizeof(Uint), 1, fp);
    fread(&idx->stats->px, sizeof(double), 1, fp);
    if(flags & PXX_STORED) {
        fread(&idx->stats->pxx, sizeof(double), 1, fp);
    }
    fread(&idx->stats->maxcover, sizeof(Uint), 1, fp);
    fread(&idx->stats->mincover, sizeof(Uint), 1, fp);
    
    /*errordenstiy*/
    fread(&idx->stats->e_N, sizeof(Uint), 1, fp);
    NFO("reading %d e-samples from index.\n", idx->stats->e_N);

    fread(&idx->stats->entropydensitystep, sizeof(double), 1, fp);
    fread(&idx->stats->entropydensitylen, sizeof(Uint), 1, fp);

    idx->stats->entropydensity = 
      ALLOCMEMORY(space, NULL, double, idx->stats->entropydensitylen);
    fread(idx->stats->entropydensity, sizeof(double), idx->stats->entropydensitylen, fp);
    
    idx->stats->entropy = ALLOCMEMORY(space, NULL, double, idx->stats->e_N);
    fread(idx->stats->entropy, sizeof(double), idx->stats->e_N, fp);
     
    if(flags & STRAND_STORED) {
      idx->stats->b = ALLOCMEMORY(space, NULL, double, idx->stats->e_N);
      fread(idx->stats->b, sizeof(double), idx->stats->e_N, fp);
      fread(&idx->stats->b_mu, sizeof(double), 1, fp);
      fread(&idx->stats->b_sd, sizeof(double), 1, fp);
      fread(&idx->stats->b_ll, sizeof(double), 1, fp);
    }
   
    idx->stats->eraw = ALLOCMEMORY(space, NULL, double, idx->stats->e_N);
    fread(idx->stats->eraw, sizeof(double), idx->stats->e_N, fp);
    
    idx->stats->e = ALLOCMEMORY(space, NULL, double, idx->stats->e_N);
    fread(idx->stats->e, sizeof(double), idx->stats->e_N, fp);

    idx->stats->e_mu = ALLOCMEMORY(space, NULL, double, 2);
    idx->stats->e_sd = ALLOCMEMORY(space, NULL, double, 2);

    fread(&idx->stats->e_mu[0], sizeof(double), 1, fp);
    fread(&idx->stats->e_mu[1], sizeof(double), 1, fp);
    fread(&idx->stats->e_sd[0], sizeof(double), 1, fp);
    fread(&idx->stats->e_sd[1], sizeof(double), 1, fp);
    fread(&idx->stats->e_ll, sizeof(double), 1, fp);

    idx->stats->e_mu = ALLOCMEMORY(space, NULL, double, 2);
    idx->stats->e_sd = ALLOCMEMORY(space, NULL, double, 2);
    idx->stats->gev_mu = ALLOCMEMORY(space, NULL, double, 2);
    idx->stats->gev_si = ALLOCMEMORY(space, NULL, double, 2); 
    idx->stats->gev_xi = ALLOCMEMORY(space, NULL, double, 2);
    idx->stats->gev_ll = ALLOCMEMORY(space, NULL, double, 2);

    
    if(flags & STATS_GEV_STORED) { 
      
      NFO("reading gev %d.\n", flags);
      fread(&idx->stats->gev_mu[0], sizeof(double), 1, fp); 
      fread(&idx->stats->gev_mu[1], sizeof(double), 1, fp);
      fread(&idx->stats->gev_si[0], sizeof(double), 1, fp); 
      fread(&idx->stats->gev_si[1], sizeof(double), 1, fp);
      fread(&idx->stats->gev_xi[0], sizeof(double), 1, fp); 
      fread(&idx->stats->gev_xi[1], sizeof(double), 1, fp);
      fread(&idx->stats->gev_ll[0], sizeof(double), 1, fp);
      fread(&idx->stats->gev_ll[1], sizeof(double), 1, fp);
    } else {
      NFO("not reading gev %d.\n", flags);
     idx->stats->gev_mu[0] = 0.044763;
     idx->stats->gev_mu[1] = 0.020171;
     idx->stats->gev_si[0] = 0.022864;
     idx->stats->gev_si[1] = 0.031077;
     idx->stats->gev_xi[0] = 0.212219;
     idx->stats->gev_xi[1] = -0.041355; 
     idx->stats->gev_ll[0] = 6291.208397;
     idx->stats->gev_ll[1] = 5908.074411;
    }
    
    fread(&idx->stats->P, sizeof(Uint), 1, fp);
    fread(&idx->stats->X, sizeof(Uint), 1, fp);
    fread(&idx->stats->N, sizeof(Uint), 1, fp);
    
    fread(&idx->stats->RR_N, sizeof(Uint), 1, fp); 
    idx->stats->RR = ALLOCMEMORY(space, NULL, Uint, 11);
    fread(idx->stats->RR, sizeof(Uint), 11, fp);

    fread(&idx->stats->MM_N, sizeof(Uint), 1, fp);
    idx->stats->MM = ALLOCMEMORY(space, NULL, Uint, 51);
    fread(idx->stats->MM, sizeof(Uint), 51, fp);

    fread(&idx->stats->areasize, sizeof(Uint), 1, fp);
    fread(&idx->stats->maxareae, sizeof(double), 1, fp);

    /*substition*/
    idx->stats->S_N = ALLOCMEMORY(space, NULL, Uint, 6);
    fread(idx->stats->S_N, sizeof(Uint), 6, fp);
    idx->stats->S = ALLOCMEMORY(space, NULL, double, 6*6);
    fread(idx->stats->S, sizeof(double), 6*6, fp);

    idx->stats->Sx_N = ALLOCMEMORY(space, NULL, Uint, 6);
    fread(idx->stats->Sx_N, sizeof(Uint), 6, fp);
    idx->stats->Sx = ALLOCMEMORY(space, NULL, double, 6*6);
    fread(idx->stats->Sx, sizeof(double), 6*6, fp);
    
    /*noise*/
    idx->stats->R_N = ALLOCMEMORY(space, NULL, Uint, 100*255);
    fread(idx->stats->R_N, sizeof(Uint), 100*255, fp);
    idx->stats->R = ALLOCMEMORY(space, NULL, Uint, 100*255);
    fread(idx->stats->R, sizeof(Uint), 100*255, fp);
    
    idx->stats->RP_N = ALLOCMEMORY(space, NULL, Uint, 100);
    fread(idx->stats->RP_N, sizeof(Uint), 100, fp);
    idx->stats->RP = ALLOCMEMORY(space, NULL, Uint, 100);
    fread(idx->stats->RP, sizeof(Uint), 100, fp);

    idx->stats->RQ_N = ALLOCMEMORY(space, NULL, Uint, 255);
    fread(idx->stats->RQ_N, sizeof(Uint), 255, fp);
    idx->stats->RQ = ALLOCMEMORY(space, NULL, Uint, 255);
    fread(idx->stats->RQ, sizeof(Uint), 255, fp);
  
    /*read variance*/
    fread(&idx->stats->V_N, sizeof(Uint), 1, fp);
    idx->stats->V = ALLOCMEMORY(space, NULL, double, idx->stats->V_N);
    fread(idx->stats->V, sizeof(double), idx->stats->V_N, fp);
    fread(&idx->stats->V_mu, sizeof(double), 1, fp);
    fread(&idx->stats->V_sd, sizeof(double), 1, fp);
    fread(&idx->stats->V_ll, sizeof(double), 1, fp);

    fread(&idx->stats->Vx_N, sizeof(Uint), 1, fp);
    idx->stats->Vx = ALLOCMEMORY(space, NULL, double, idx->stats->Vx_N);
    fread(idx->stats->Vx, sizeof(double), idx->stats->Vx_N, fp);
    fread(&idx->stats->Vx_mu, sizeof(double), 1, fp);
    fread(&idx->stats->Vx_sd, sizeof(double), 1, fp);
    fread(&idx->stats->Vx_ll, sizeof(double), 1, fp);
 
    if(flags & MOTIF_STORED) { 
      idx->stats->MO_N = ALLOCMEMORY(space, NULL, Uint, 1024);
      fread(idx->stats->MO_N, sizeof(Uint), 1024, fp);
      idx->stats->MO = ALLOCMEMORY(space, NULL, Uint, 1024);
      fread(idx->stats->MO, sizeof(Uint), 1024, fp);
    }
 
    fread(&idx->stats->s_N, sizeof(Uint), 1, fp);
    NFO("reading %d scores from index.\n", idx->stats->s_N);
    idx->stats->s = ALLOCMEMORY(space, NULL, double, idx->stats->s_N);
    fread(idx->stats->s, sizeof(double), idx->stats->s_N, fp);

  }


  fclose(fp);
  return idx ;
}

/*-------------------------- bl_matchfileWriteIndex --------------------------
 *    
 * @brief write the index (including statistics) to disk
 * @author Steve Hoffmann 
 *   
 */
 
void bl_matchfileWriteIndex(matchfileindex_t *idx, char *filename) {
  
  FILE *fp;
  unsigned char flags = 0;
  Uint i, j, nchrom, nreads, exp, len;
  matchfileBin_t *bin;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    DBG("Couldn't open file %s. Exit forced.\n", filename);
    exit(-1);
  }

  if(idx->gzindex) {
    flags |= GZ_IDX_STORED;
  }

  if(idx->stats) {
    flags |= STATS_TAB_STORED;
    flags |= PXX_STORED;
    if(idx->stats->gev_mu) {
      flags |= STATS_GEV_STORED;
    }
    if(idx->stats->MO) {
      flags |= MOTIF_STORED;
    }
    if(idx->stats->b) {
      flags |= STRAND_STORED;
    }
  }

  if(idx->md5) {
    flags |= IDXMD5_STORED;
  }

  fwrite(&flags, sizeof(unsigned char), 1, fp);

  if(idx->gzindex) { 
    fwrite(&idx->gzindex->have, sizeof(int), 1, fp);
    fwrite(&idx->gzindex->size, sizeof(int), 1, fp);

   struct point *list = idx->gzindex->list;
   for(i=0; i < idx->gzindex->size; i++) {

      fwrite(&list[i].out, sizeof(off_t), 1, fp);
      fwrite(&list[i].in, sizeof(off_t), 1, fp);
      fwrite(&list[i].bits, sizeof(int), 1, fp);
      fwrite(list[i].window, sizeof(unsigned char), WINSIZE, fp);
    }
  }

  nchrom = idx->noofchroms;
  nreads = idx->noofreads;
  exp = idx->exp;

  fwrite(&exp, sizeof(Uint), 1, fp);
  fwrite(&nchrom, sizeof(Uint), 1, fp);
  fwrite(&nreads, sizeof(Uint), 1, fp);

  for(i=0; i < nchrom; i++) {
    len = strlen(idx->chromnames[i]);
    fwrite(&len, sizeof(Uint), 1, fp);
    fwrite(idx->chromnames[i], sizeof(char), len, fp);
  }

  fwrite(idx->matchstart, sizeof(Uint), nchrom, fp);
  fwrite(idx->matchend, sizeof(Uint), nchrom, fp);  
  fwrite(idx->noofbins, sizeof(Uint), nchrom, fp);

  for(i=0; i < nchrom; i++) {

    //fwrite(idx->bins[i], sizeof(matchfileBin_t), idx->noofbins[i], fp);
    for(j=0; j < idx->noofbins[i]; j++) {
      bin = idx->bins[i];
      fwrite(&bin[j].start, sizeof(Uint), 1, fp);
      fwrite(&bin[j].end, sizeof(Uint), 1, fp);
      fwrite(&bin[j].matches, sizeof(Uint), 1, fp);
      fwrite(&bin[j].offset, sizeof(off_t), 1, fp);
      fwrite(&bin[j].endoff, sizeof(off_t), 1, fp);
    }
  }

  fwrite(&idx->maxreadlen, sizeof(Uint), 1, fp);
  fwrite(idx->submatrix, sizeof(double), (6*6*QRNGE*MAXREADLENGTH), fp);
  fwrite(&idx->mean_coverage, sizeof(double), 1, fp);
  fwrite(&idx->mean_qual, sizeof(double), 1, fp);

  fwrite(idx->P_ERR, sizeof(Uint), MAXREADLENGTH, fp);
  fwrite(idx->Q_ERR, sizeof(Uint), QRNGE, fp);
  fwrite(idx->Q_N, sizeof(Uint), QRNGE, fp);

  if(idx->stats) {
  
    fwrite(&idx->stats->n, sizeof(Uint), 1, fp);
    fwrite(&idx->stats->px, sizeof(double), 1, fp);
    fwrite(&idx->stats->pxx, sizeof(double), 1, fp);
    fwrite(&idx->stats->maxcover, sizeof(Uint), 1, fp);
    fwrite(&idx->stats->mincover, sizeof(Uint), 1, fp);
    
    /*errordensity*/
    NFO("writing %d e-samples to index.\n", idx->stats->e_N);
    fwrite(&idx->stats->e_N, sizeof(Uint), 1, fp);
    
    fwrite(&idx->stats->entropydensitystep, sizeof(double), 1, fp);
    fwrite(&idx->stats->entropydensitylen, sizeof(Uint), 1, fp); 
    
    fwrite(idx->stats->entropydensity, sizeof(double), idx->stats->entropydensitylen, fp); 
    fwrite(idx->stats->entropy, sizeof(double), idx->stats->e_N, fp); 
    
    if(idx->stats->b) {
     fwrite(idx->stats->b, sizeof(double), idx->stats->e_N, fp); 
     fwrite(&idx->stats->b_mu, sizeof(double), 1, fp); 
     fwrite(&idx->stats->b_sd, sizeof(double), 1, fp);
     fwrite(&idx->stats->b_ll, sizeof(double), 1, fp);
    }
    
    fwrite(idx->stats->eraw, sizeof(double), idx->stats->e_N, fp); 
    fwrite(idx->stats->e, sizeof(double), idx->stats->e_N, fp); 
    fwrite(&idx->stats->e_mu[0], sizeof(double), 1, fp); 
    fwrite(&idx->stats->e_mu[1], sizeof(double), 1, fp);
    fwrite(&idx->stats->e_sd[0], sizeof(double), 1, fp); 
    fwrite(&idx->stats->e_sd[1], sizeof(double), 1, fp);
    fwrite(&idx->stats->e_ll, sizeof(double), 1, fp);

    if(idx->stats->gev_mu) {
      
      fwrite(&idx->stats->gev_mu[0], sizeof(double), 1, fp); 
      fwrite(&idx->stats->gev_mu[1], sizeof(double), 1, fp);
      fwrite(&idx->stats->gev_si[0], sizeof(double), 1, fp); 
      fwrite(&idx->stats->gev_si[1], sizeof(double), 1, fp);
      fwrite(&idx->stats->gev_xi[0], sizeof(double), 1, fp); 
      fwrite(&idx->stats->gev_xi[1], sizeof(double), 1, fp);
      fwrite(&idx->stats->gev_ll[0], sizeof(double), 1, fp);
      fwrite(&idx->stats->gev_ll[1], sizeof(double), 1, fp);
      
    }
    
    fwrite(&idx->stats->P, sizeof(Uint), 1, fp);
    fwrite(&idx->stats->X, sizeof(Uint), 1, fp);
    fwrite(&idx->stats->N, sizeof(Uint), 1, fp);
    
    fwrite(&idx->stats->RR_N, sizeof(Uint), 1, fp);
    fwrite(idx->stats->RR, sizeof(Uint), 11, fp);
    fwrite(&idx->stats->MM_N, sizeof(Uint), 1, fp);
    fwrite(idx->stats->MM, sizeof(Uint), 51, fp);

    fwrite(&idx->stats->areasize, sizeof(Uint), 1, fp);
    fwrite(&idx->stats->maxareae, sizeof(double), 1, fp);

    /*substitution*/
    fwrite(idx->stats->S_N, sizeof(Uint), 6, fp);
    fwrite(idx->stats->S, sizeof(double), 6*6, fp);
    fwrite(idx->stats->Sx_N, sizeof(Uint), 6, fp);
    fwrite(idx->stats->Sx, sizeof(double), 6*6, fp);

    /*noise*/
    fwrite(idx->stats->R_N, sizeof(Uint), 100*255, fp);
    fwrite(idx->stats->R, sizeof(Uint), 100*255, fp);
    fwrite(idx->stats->RP_N, sizeof(Uint), 100, fp);
    fwrite(idx->stats->RP, sizeof(Uint), 100, fp);
    fwrite(idx->stats->RQ_N, sizeof(Uint), 255, fp);
    fwrite(idx->stats->RQ, sizeof(Uint), 255, fp);
    
    /*readvariance*/
    fwrite(&idx->stats->V_N, sizeof(Uint), 1, fp);
    fwrite(idx->stats->V, sizeof(double), idx->stats->V_N, fp); 
    fwrite(&idx->stats->V_mu, sizeof(double), 1, fp); 
    fwrite(&idx->stats->V_sd, sizeof(double), 1, fp); 
    fwrite(&idx->stats->V_ll, sizeof(double), 1, fp); 
  
    fwrite(&idx->stats->Vx_N, sizeof(Uint), 1, fp);
    fwrite(idx->stats->Vx, sizeof(double), idx->stats->Vx_N, fp); 
    fwrite(&idx->stats->Vx_mu, sizeof(double), 1, fp); 
    fwrite(&idx->stats->Vx_sd, sizeof(double), 1, fp); 
    fwrite(&idx->stats->Vx_ll, sizeof(double), 1, fp); 
 
    /*motif*/
    if(idx->stats->MO) { 
        fwrite(idx->stats->MO_N, sizeof(Uint), 1024, fp);
        fwrite(idx->stats->MO, sizeof(Uint), 1024, fp);
    }
  
    NFO("writing %d scores to index.\n", idx->stats->s_N);
    fwrite(&idx->stats->s_N, sizeof(Uint), 1, fp);
    fwrite(idx->stats->s, sizeof(double), idx->stats->s_N, fp); 

  }

  fclose(fp);
  return;
}

Uint
bl_compareIndices(matchfileindex_t *i1, matchfileindex_t *i2){

  Uint i, j;

  if(i1->gzindex == NULL && i2->gzindex == NULL) return 0;

  if(i1->gzindex->have != i2->gzindex->have) return 1;
  if(i1->gzindex->size != i2->gzindex->size) return 2;

  for(i=0; i < i1->gzindex->size; i++) {
    if(i1->gzindex->list[i].out != i1->gzindex->list[i].out) return 3;
    if(i1->gzindex->list[i].in != i1->gzindex->list[i].in) return 4;
    if(i1->gzindex->list[i].bits != i1->gzindex->list[i].bits) return 5;
    
    for(j=0; j < WINSIZE; j++) {
      if(i1->gzindex->list[i].window[j] != 
          i1->gzindex->list[i].window[j]) return 6;
    }
  }

  if(i1->exp != i2->exp) return 7;
  if(i1->noofchroms != i2->noofchroms) return 8;
  if(i1->noofreads != i2->noofreads) return 9;

  for(i=0; i < i1->noofchroms; i++) {
    if(strcmp(i1->chromnames[i], i2->chromnames[i]) != 0) return 10;
    if(i1->matchstart[i] != i2->matchstart[i]) return 11;
    if(i1->matchend[i] != i2->matchend[i]) return 12;
    if(i1->noofbins[i] != i2->noofbins[i]) return 13;

    for(j=0; j < i1->noofbins[i]; j++) {
      if(i1->bins[i][j].start != i2->bins[i][j].start) return 14;
      if(i1->bins[i][j].end != i2->bins[i][j].end) return 15;
      if(i1->bins[i][j].matches != i2->bins[i][j].matches) return 16;
      if(i1->bins[i][j].offset != i2->bins[i][j].offset) return 17;
      if(i1->bins[i][j].endoff != i2->bins[i][j].endoff) return 18;
    }
  }


  if(i1->maxreadlen != i2->maxreadlen) return 19;
  if(i1->mean_coverage != i2->mean_coverage) return 20;
  if(i1->mean_qual != i2->mean_qual) {
    return 21;
  }

  for(i=0; i < MAXREADLENGTH; i++) { 
    if(i1->P_ERR[i] != i2->P_ERR[i]) return 22;
  }

  for(i=0; i < QRNGE; i++) { 
    if(i1->Q_ERR[i] != i2->Q_ERR[i]) return 23;
    if(i1->Q_N[i] != i2->Q_N[i]) return 24;
  }

  for(i=0; i < QRNGE*MAXREADLENGTH*6*6; i++) { 
    if(i1->submatrix[i] != i2->submatrix[i]) return 24;
  }


  if(i1->stats && i2->stats) {
    
    if(i1->stats->n != i2->stats->n) return 25;
    if(i1->stats->px != i2->stats->px) return 26;
    if(i1->stats->maxcover != i2->stats->maxcover) return 27;
    if(i1->stats->mincover != i2->stats->mincover) return 28;
    if(i1->stats->e_N != i2->stats->e_N) return 29;

    for(i=0; i < i1->stats->e_N; i++) {
      if(i1->stats->e[i] != i2->stats->e[i]) return 30;
    }

    if(i1->stats->e_mu[0] != i2->stats->e_mu[0]) return 31;
    if(i1->stats->e_mu[1] != i2->stats->e_mu[1]) return 32;
    if(i1->stats->e_sd[0] != i2->stats->e_sd[0]) return 33;
    if(i1->stats->e_sd[1] != i2->stats->e_sd[1]) return 34;

    if(i1->stats->e_ll != i2->stats->e_ll) return 35;
    if(i1->stats->areasize != i2->stats->areasize) return 36;
    if(i1->stats->maxareae != i2->stats->maxareae) return 37;

    for(i=0; i < 6*6; i++) {
      if(i1->stats->S[i] != i2->stats->S[i]) return 38;
      if(i1->stats->Sx[i] != i2->stats->Sx[i]) return 39;
    }

    for(i=0; i < 6; i++) {
      if(i1->stats->S_N[i] != i2->stats->S_N[i]) return 40;
      if(i1->stats->Sx_N[i] != i2->stats->Sx_N[i]) return 41;
    }

    for(i=0; i < 100*255; i++) {
      if(i1->stats->R_N[i] != i2->stats->R_N[i]) return 42;
      if(i1->stats->R[i] != i2->stats->R[i]) return 43;
    }
    
    for(i=0; i < 100; i++) {
      if(i1->stats->RP_N[i] != i2->stats->RP_N[i]) return 42;
      if(i1->stats->RP[i] != i2->stats->RP[i]) return 43;
    }

    for(i=0; i < 255; i++) {
      if(i1->stats->RQ_N[i] != i2->stats->RQ_N[i]) return 42;
      if(i1->stats->RQ[i] != i2->stats->RQ[i]) return 43;
    }


    if (i1->stats->V_N != i2->stats->V_N) return 45;

    for(i=0; i < i1->stats->V_N; i++) {
      if(i1->stats->V[i] != i2->stats->V[i]) return 46;
    }

    if (i1->stats->V_mu != i2->stats->V_mu) return 47;
    if (i1->stats->V_sd != i2->stats->V_sd) return 48;
    if (i1->stats->V_ll != i2->stats->V_ll) return 49;

    if (i1->stats->Vx_N != i2->stats->Vx_N) return 50;

    for(i=0; i < i1->stats->Vx_N; i++) {
      if(i1->stats->Vx[i] != i2->stats->Vx[i]) return 51;
    }
  
    if (i1->stats->Vx_mu != i2->stats->Vx_mu) {
      fprintf(stderr, "vx_mu :%f!=%f:vx_mu\n", i1->stats->Vx_mu, i2->stats->Vx_mu);
      fprintf(stderr, "vx_sd :%f!=%f:vx_sd\n", i1->stats->Vx_sd, i2->stats->Vx_sd);
      fprintf(stderr, "vx_ll :%f!=%f:vx_ll\n", i1->stats->Vx_ll, i2->stats->Vx_ll);
      return 52;
    }
    if (i1->stats->Vx_sd != i2->stats->Vx_sd) return 53;
    if (i1->stats->Vx_ll != i2->stats->Vx_ll) return 54;
  
  }

  return 0;
}

/*-------------------------- bl_matchfileDumpStats ---------------------------
 *    
 * @brief dump the match file stats
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchfileDumpFileStats (void *space, matchfile_t *file)
{

  Uint i, sum;
  Uint *Q_ERR, *P_ERR, *Q_N; 

  Q_N = file->index->Q_N;
  Q_ERR = file->index->Q_ERR;
  P_ERR = file->index->P_ERR;
   
  fprintf(stderr, "Q    " );
  for(i=0; i < QRNGE; i++){
    fprintf(stderr, "%d\t\t\t", MINQUAL+i);
  }

  sum = 0;

  for(i=0; i < QRNGE; i++) { 
    sum += Q_N[i];
  }

  fprintf(stderr, "Q_N  \t");
  for(i=0; i < QRNGE; i++) { 
    fprintf(stderr, "%.1f\t\t\t", (double)Q_N[i]/sum);
  }

  fprintf(stderr, "\nQ_COR\t");
  for(i=0; i < QRNGE; i++) { 
    fprintf(stderr, "%d (%.3f)\t\t\t", Q_N[i]-Q_ERR[i], 
        (double)((double)Q_N[i]-Q_ERR[i])/Q_N[i]);
  }

  fprintf(stderr, "\nQ_ERR\t");
  for(i=0; i < QRNGE; i++) { 
    fprintf(stderr, "%d (%.3f)\t\t\t", Q_ERR[i], (double)Q_ERR[i]/Q_N[i]);
  }

  fprintf(stderr, "\n\nP_ERR\t");
  for(i=0; i < MAXREADLENGTH; i++) { 
    fprintf(stderr, "%.6f\t\t\t", (double)P_ERR[i]/file->index->noofreads);
  }

  fprintf(stderr,"\n");

  return ;
}


/*----------------------- bl_matchfileDumpCrossSection -----------------------
 *    
 * @brief dump cross section
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchfileDumpCrossSection (matchfileCross_t *cs, Uint pos)
{
  
    printf("%c\t", cs->cons);
    printf("%d\t", cs->len);
    printf("%s\t", cs->chars);
    printf("%s\t", cs->quals);
    printf("\n");

	return ;
}

/*------------------------- bl_matchfileAdjustBounds -------------------------
 *    
 * @brief check and adjust frame bounds if violated
 * @author Steve Hoffmann 
 *   
 */
 

Uint
bl_matchfileAdjustBounds(void *space, matchfile_t *file, fasta_t *set,
    Uint setid, Uint matid, Uint start, Uint width, 
    Uint *newstart, Uint *newwidth) {
 
  Uint maxmat=0, maxchr=0, len=0;
 
  *newstart = start;
  *newwidth = width;

  if(set && setid < set->noofseqs) {
    maxchr = bl_fastaGetSequenceLength(set, setid); 
  }

  if(matid < file->index->noofchroms && 
      file->index->noofbins[matid] && file->index->bins[matid]) {
    maxmat = file->index->bins[matid][file->index->noofbins[matid]-1].end;
  }

  len = maxmat;
  if(maxmat) {
    if(start+width > maxmat) {
      if(maxmat > width) {
        *newstart = maxmat-width;     
      } else {
        *newstart = 1;
        *newwidth = maxmat;
      }

    }
  }

  if(maxchr) {
    len = maxchr;
    *newstart = start;
    if(start+width > maxchr) {
      if (maxchr > width ) { 
        *newstart = maxchr - width;
      } else {
        *newstart = 1;
        *newwidth = maxchr;
      }
    }
  }

  return len;
}

/*--------------------------- bl_matchfileGetFrame ---------------------------
 *    
 * @brief get a frame of specified width from indexed matchfile
 * ATTENTION: start is 1-offset
 * @author Steve Hoffmann 
 *   
 */
 
matchfileFrame_t *
bl_matchfileGetFrame(void *space, matchfile_t *file, char *chrname, 
    Uint start, Uint width, fasta_t *set, Uint maxcover, matchfileCross_t *inputcs) {
 
  Uint i=0, k=0, len, adjstart=0, adjwidth=0;
  matchfileFrame_t *frame;
 
  frame = ALLOCMEMORY(space, NULL, matchfileFrame_t, 1);
  frame->chrname = chrname;
  frame->start = start;
  frame->width = width; 
  frame->ref = NULL;
  frame->chrseq = NULL;
  
  k = bl_matchfileGetChromIndexNumber(file->index, chrname);


  if (set) i = bl_fastxFindIDIdx(chrname, set);


  len = bl_matchfileAdjustBounds(space, file, set, i, k, start, width, 
      &adjstart, &adjwidth);

  
  frame->start = adjstart;
  frame->width = adjwidth; 
  frame->chrlen = len;

#ifdef DBGIDX
  DBG("frame [%d,%d] (chrlen:%d)\n", frame->start, 
      frame->start+frame->width, frame->chrlen);
#endif

  if (i < set->noofseqs){
    frame->chrseq = bl_fastaGetSequence(set, i);
    frame->ref = &bl_fastaGetSequence(set, i)[adjstart-1];
    frame->chrname = bl_fastaGetDescription(set, i);
  }

  frame->cs = bl_matchfileRead(space, file, chrname, frame->start, 
      frame->start+adjwidth-1, maxcover, set, 255, inputcs);


  return frame;
}


/*------------------------- bl_matchfileMergeFrames --------------------------
 *    
 * @brief merge two frames, ie. attach frame g to frame f
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchfileMergeFrames (void *space, matchfileFrame_t* f, matchfileFrame_t *g)
{

  assert(!strcmp(f->chrname,g->chrname));
  f->cs = ALLOCMEMORY(NULL, f->cs, matchfileCross_t, f->width+g->width);
  memmove(&f->cs[f->width], g->cs, sizeof(matchfileCross_t)*g->width);
  f->width += g->width;
	
  return ;
}

/*------------------------------ bl_writewiggle ------------------------------
 *    
 * @brief write a wiggle to a file
 * @author Steve Hoffmann 
 *   
 */
 
Uint
bl_writeexpression (void *space, char *filename, Uint filenamelen,
    matchfile_t *file, fasta_t *set, Uint width, Uint maxcover)
{

  char *chromname, *minusfilename, *plusfilename, *basename, *buffer;
  FILE *plus, *minus, *outfile;
  size_t buffersize = 1024, len;
  matchfileCross_t *cs;
  char mheader, pheader;
  Uint i, j, k, s, idx, mstrands, pstrands, basenamelen;


  basename = bl_basename(filename);
  basenamelen = strlen(basename);

  minusfilename = bl_getTempFile(basename, basenamelen);
  plusfilename = bl_getTempFile(basename, basenamelen);

  minus = fopen(minusfilename, "w");
  plus = fopen(plusfilename, "w");


  for(k=0; k < file->index->noofchroms; k++) {
    chromname = file->index->chromnames[k];
    mheader = 0;
    pheader = 0;

    for(j=file->index->matchstart[k]; 
        j < file->index->matchend[k]; j+= width) {

      idx = bl_fastxFindIDIdx(chromname, set);

      if (idx < set->noofseqs){
        chromname = bl_fastaGetDescription(set, idx);
      } else {
        NFO("the chromsome %s was not found in the db. Exit forced!\n", chromname);
        exit(-1);
      }

      NFO("reading interval [%d,%d] on '%s'\n", j, j+width-1, chromname);

      cs = bl_matchfileRead(space, file, chromname, 
          j, j+width-1, maxcover, set, 0, NULL); 

      NFO("dumping interval [%d,%d] to device\n", j, j+width-1);

      for(i=0; i < width; i++) {
        mstrands = 0;
        pstrands = 0;
        for(s=0; s < cs[i].len; s++) {
          if (cs[i].strands[s] == '+')
            pstrands++;
          else
            mstrands++;
        }
        assert(mstrands+pstrands == cs[i].len);
        if(pstrands) {
          if(!pheader) {
            fprintf(plus, "track type=wiggle_0 name=\"%s_(+)\" description=\"%s exp (+)\"\n",
                basename, basename);
            fprintf(plus, "variableStep chrom=%s\n", chromname);
          }
          fprintf(plus, "%d\t%d\n", j+i, pstrands);
          pheader = 1;
        }

        if(mstrands){
          if(!mheader) {
            fprintf(minus,"track type=wiggle_0 name=\"%s_(-)\" description=\"%s exp (-)\"\n",
                basename, basename);
            fprintf(minus, "variableStep chrom=%s\n", chromname);
          }
          fprintf(minus, "%d\t%d\n", j+i, mstrands);
          mheader = 1;
        }
      }

      bl_matchfileDestructCross(space, cs, width);
      FREEMEMORY(space, cs);
    }
  }

  fclose(minus);
  fclose(plus);

  buffer = ALLOCMEMORY(space, NULL, char, buffersize);

  minus = fopen(minusfilename, "rb");
  outfile = fopen(filename, "wb");

  while((len = fread(buffer, 1, buffersize, minus)) > 0) {
    fwrite(buffer, 1, len, outfile);

  }

  fclose(minus);

  plus = fopen(plusfilename, "rb");
  while((len = fread(buffer, 1, buffersize, plus)) > 0) {
    fwrite(buffer, 1, len, outfile);    
  }

  fclose(plus);
  fclose(outfile);
  bl_rm(space, minusfilename);
  bl_rm(space, plusfilename);

  FREEMEMORY(space, buffer);
  FREEMEMORY(space, minusfilename);
  FREEMEMORY(space, plusfilename);

  return 0;
}


/*----------------------------- bl_writeVars2vix -----------------------------
 *    
 * @brief write variants to internal format
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_writeVars2vix ( )
{
	return ;
}

