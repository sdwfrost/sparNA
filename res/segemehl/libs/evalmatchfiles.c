
/*
 *  evalmatchfiles.c
 *  evalutation/statistics
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 06.10.2010 20:26:00 CEST
 *  
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include "sort.h"
#include "alignment.h"
#include "stringutils.h"
#include "basic-types.h"
#include "mathematics.h"
#include "matfile.h"
#include "bitVector.h"
#include "info.h"
#include "vtprogressbar.h"
#include "fileio.h"
#include "matchfilesfields.h"
#include "matchfiles.h"
#include "debug.h"
#include "evalmatchfiles.h"
#include "biofiles.h"
#include "splicesites.h"



/*-------------------------- bl_matchfileSampleCmp ---------------------------
 *    
 * @brief cmp chromosome positions sampled below 
 * @author Steve Hoffmann 
 *   
 */
 
Uint
bl_matchfileSampleCmp (Uint elemA, Uint elemB, void *toSort, void *info)
{
  PairUint *samples;
  samples = (PairUint*) toSort;

  if(samples[elemA].a > samples[elemB].a) 
    return 1;
  if(samples[elemA].a < samples[elemB].a)
    return 2;

  if(samples[elemA].b > samples[elemB].b)
    return 1;
  if(samples[elemA].b < samples[elemB].b)
    return 2;
	
  return 0;
}


/*---------------------- bl_matchfileGetCrossConsError -----------------------
 *    
 * @brief get the error e of a cross section based on the consensus 
 * (not reference)
 * @author Steve Hoffmann 
 *   
 */
 
double
bl_matchfileGetCrossConsError (matchfileFrame_t *frame, Uint pos)
{
  Uint j;
  double e=0;
 
  for(j=0; j < frame->cs[pos].len; j++) {
    if(frame->cs[pos].chars[j] != frame->cs[pos].cons)
      e+=1.0;
  }

  if(frame->cs[pos].len) { 
    e /= (double)frame->cs[pos].len;
  }

  return e;
}

/*---------------------- bl_matchfileGetCrossRefError ------------------------
 *    
 * @brief get the error e of a cross section based on the reference
 * @author Steve Hoffmann 
 *   
 */
 
double
bl_matchfileGetCrossRefError (matchfileFrame_t *frame, Uint pos)
{
  Uint j;
  double e=0;
 
  for(j=0; j < frame->cs[pos].len; j++) {
    if(frame->cs[pos].chars[j] != frame->cs[pos].ref)
      e+=1.0;
  }

  if(frame->cs[pos].len) { 
    e /= (double)frame->cs[pos].len;
  }

  return e;
}

/*------------------------- bl_matchfileGetErrorDensity --------------------------
 *    
 * @brief get sample statistics
 * @author Steve Hoffmann 
 *   
 */

  void
bl_matchfileGetErrorDensity(void *space, matchfileFrame_t *frame, 
    Uint pos, matchfileFrameStats_t *framestats, void *nfo)
{  
  double e, er=.0;
  matchfileSampleStats_t *stats = (matchfileSampleStats_t*) nfo; 
  Uint mincover = stats->mincover;
  Uint maxcover = stats->maxcover;


  if (frame->cs[pos].len < mincover || frame->cs[pos].len > maxcover) { 
    return;
  }

  e = bl_matchfileGetCrossConsError(frame, pos);

  if(e > 0 && stats->e_N < stats->n) {   
    stats->eraw[stats->e_N]=e;
    stats->e[stats->e_N++]=e-er;
  }

  return ;
}


/*------------------------- bl_matchfileGetNTCounts --------------------------
 *    
 * @brief for each cross section in the frame: 
 * Get the count for all nucleotides
 * @author Steve Hoffmann 
 *   
 */

Uint*
bl_matchfileGetNTCounts(matchfileCross_t *cs) {
  Uint j, *cnt;

  cnt = ALLOCMEMORY(space, NULL, Uint, 256);
  memset(cnt, 0, sizeof(Uint)*256);

  for(j=0; j < cs->len; j++) {
    cnt[(int)cs->chars[j]]++;
  }

  return cnt;
}

/*----------------------- bl_matchfileGetPlusCounts -----------------------
 *    
 * @brief for each cross section in the frame: 
 * Get the count for all strands
 * @author Steve Hoffmann 
 *   
 */

Uint*
bl_matchfileGetPlusCounts(matchfileCross_t *cs) {
  Uint j, *cnt;

  cnt = ALLOCMEMORY(space, NULL, Uint, 256);
  memset(cnt, 0, sizeof(Uint)*256);

  for(j=0; j < cs->len; j++) {
    if(cs->strands[j] == '+') cnt[(int)cs->chars[j]]++;
  }

  return cnt;
}


/*------------------------- bl_matchfileGetConsensus -------------------------
 *    
 * @brief calculate the consensus bases in a frame
 * @author Steve Hoffmann 
 *   
 */

void
bl_matchfileGetConsensus(matchfileFrame_t *frame) {
  Uint i, j, *cnt, max;
  double lentropy;
  double rentropy;
  double longentropy;
  char gapalign = 0;
  Uint *tab = NULL;

  tab = ALLOCMEMORY(NULL, NULL, Uint, 256);
  memset(tab, 0, sizeof(Uint)*256);
 
  tab['A'] = 1;
  tab['C'] = 2;
  tab['G'] = 3;
  tab['T'] = 4;
  
  for(i=0; i < frame->width; i++) {
    

  frame->cs[i].entropy = 2.0; 
  frame->cs[i].longentropy = 2.0;

    if (frame->ref) {
      frame->cs[i].ref = frame->ref[i];
      lentropy = 2.0;
      rentropy = 2.0;
      longentropy = 2.0;
 
      if(frame->start+i > 11 && frame->start+i+11 < frame->chrlen) {
        longentropy=
          shannonentropy(NULL, &frame->chrseq[frame->start+i-10], 21, 6, tab);  
      } else { 
        if(frame->start+i > 21) {
        longentropy=
          shannonentropy(NULL, &frame->chrseq[frame->start+i-21], 21, 6, tab); 
        }
        if(frame->start+i+21 < frame->chrlen) { 
        longentropy = 
          shannonentropy(NULL, &frame->chrseq[frame->start+i], 21, 6, tab); 
        }
      }

      if(frame->start+i > 11 && frame->start+i+11 < frame->chrlen) {
        lentropy=
          shannonentropy(NULL, &frame->chrseq[frame->start+i-10], 10, 6, tab); 
        rentropy = 
          shannonentropy(NULL, &frame->chrseq[frame->start+i], 10, 6, tab); 
      } else { 
        if(frame->start+i > 11) {
        lentropy=
          shannonentropy(NULL, &frame->chrseq[frame->start+i-10], 10, 6, tab); 
        }
        if(frame->start+i+11 < frame->chrlen) { 
        rentropy = 
          shannonentropy(NULL, &frame->chrseq[frame->start+i], 10, 6, tab); 
        }
      }
      
      frame->cs[i].entropy = MIN(MIN(lentropy,rentropy),2.0); 
    // frame->cs[i].entropy = MIN(MAX(lentropy,rentropy),2.0); 
      frame->cs[i].longentropy = MIN(longentropy,2.0);
    }

    frame->cs[i].cons = '^';
    cnt = bl_matchfileGetNTCounts(&frame->cs[i]);

    if(frame->cs[i].len){
      max = uarraymax(cnt, 255);
      frame->cs[i].cons = (char)max;
    }

    FREEMEMORY(space, cnt);

    gapalign = 0;
    if(frame->cs[i].noofdels > 1) { 
      for(j=0; j < frame->cs[i].noofdels; j++) {
        if(frame->cs[i].dels[j].len > 1) {
          gapalign = 1;
        }
      }
      if(gapalign) 
      bl_matchfileGapAlign(frame->cs[i].dels, frame->cs[i].noofdels);
    }
  }

  FREEMEMORY(space, tab);
  return;
}


/*-------------------------- bl_matchfileCrossPrint --------------------------
 *    
 * @brief print a cross section to stderr
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchfileCrossPrint(void *space, matchfileFrame_t *frame) {
  Uint i, j, width; //, start;
  matchfileCross_t *cs;

  cs = frame->cs;
  //not used: start = frame->start;
  width = frame->width;
  for(i=0; i < width; i++) {
    fprintf(stderr, "%d: %d\t%d\t%d\t%s\t%s\t", i, cs[i].len, 
        cs[i].starts, cs[i].ends, cs[i].chars, 
        cs[i].quals);
    for(j=0; j < cs[i].len; j++) {
      fprintf(stderr, "%d,", cs[i].row[j]);
    }
    fprintf(stderr, "\t");

    for(j=0; j < cs[i].len; j++) {
      fprintf(stderr, "%d,", cs[i].edist[j]);
    }
    fprintf(stderr, "\n");
  }
}


/*----------------------- bl_matchfileGetScoreSample ------------------------
 *    
 * @brief get a sample of the scores
 * @author Steve Hoffmann 
 *   
 */

  void
bl_matchfileGetScoreSample (void *space, matchfileFrame_t *frame,
    Uint pos, matchfileFrameStats_t *framestats, void *nfo)
{

  matchfileindex_t *idx = (matchfileindex_t *) nfo;
  matchfileSampleStats_t *stats = (matchfileSampleStats_t *) idx->stats;
  Uint mincover = stats->mincover;
  Uint maxcover = stats->maxcover;
   
  if(frame->cs[pos].len < mincover || frame->cs[pos].len > maxcover) {
    return;
  }

  bl_matchfileTest(space, 0, 0, pos, &frame->cs[pos], frame->cs[pos].ref, 
      idx, 0, NULL);

  if(frame->cs[pos].scr_sample > log(0) || frame->cs[pos].scr_cons > log(0)) {
    stats->s[stats->s_N++] = MAX(frame->cs[pos].scr_sample, 
        frame->cs[pos].scr_cons);
   } else {
 //    fprintf(stdout, "%f, %f, %f, rt:%f, rq:%f, mm:%f, minrp:%f, maxrp:%f\n", frame->cs[pos].p_refx, frame->cs[pos].p_ref, log(ecdf(frame->cs[pos].ee, stats->ecdf)), frame->cs[pos].diff_rt, frame->cs[pos].diff_rq, frame->cs[pos].diff_mm, stats->minrp, stats->maxrp);
   }

  return ;
}

/*------------------------- bl_matchfileGetConditional ----------------------
 *    
 * @brief sample values for the scores
 * @author Steve Hoffmann 
 *   
 */

  void
bl_matchfileGetConditionals(void *space, matchfileFrame_t *frame, 
    Uint pos, matchfileFrameStats_t *framestats, void *nfo)
{  
  matchfileSampleStats_t *stats = (matchfileSampleStats_t *) nfo;
  double e,  pxx = stats->pxx; 
  int rpos;
  Uint j;
  
  if(frame->cs[pos].len < stats->mincover 
      || frame->cs[pos].len > stats->maxcover) {
    return;
  }
  
  e = bl_matchfileGetCrossConsError(frame, pos); 

  if(e < pxx) { 

  for(j=0; j < frame->cs[pos].len; j++) { 


    rpos = 
      trunc(((double)(((double)frame->cs[pos].readpos[j]*100.0)
              /((double)frame->cs[pos].readlen[j]))));

    MATRIX2D(stats->R_N, 255, rpos, (int)frame->cs[pos].quals[j])++;

    stats->RP_N[rpos] += 1;


    stats->RQ_N[(int)frame->cs[pos].quals[j]]+=1;

    if(frame->cs[pos].chars[j] != frame->cs[pos].cons) { 
      MATRIX2D(stats->R, 255, rpos, (int)frame->cs[pos].quals[j])++;
      stats->RP[rpos] += 1;
      stats->RQ[(int)frame->cs[pos].quals[j]]++;
    } 
  }
  }

   if(e > 0 && e <pxx) { 
    for(j=0; j < frame->cs[pos].len; j++) { 

    if(frame->cs[pos].edist[j] <= 10)
      stats->RR[(int)frame->cs[pos].edist[j]]++;
    else
      stats->RR[10]++;
    stats->RR_N++;

    if(frame->cs[pos].matchcnt[j] <= 50)
      stats->MM[frame->cs[pos].matchcnt[j]]++;
    else
      stats->MM[50]++;
    stats->MM_N++;

  }
  }

  return ;
}


/*--------------------------- bl_matchfileSmallMap ---------------------------
 *    
 * @brief get a small map to quickly find expressed/covered sites
 * @author Steve Hoffmann 
 *   
 */
 
unsigned char**
bl_matchfileSmallMap (void *space, matchfile_t* file, Uint **mapsize)
{

  FILE *fp = NULL;
  stringset_t *fields = NULL;
  char *buffer = NULL, ch, *curchrom = NULL, *filename;
  unsigned char **map = NULL;
  Uint buffersize = 1024, len = 0, curstart = 0, 
       curend = 0, i, j, bin, lastbin=0, id, lastid=-1, *msz,
       noofseqs =0;

  matchfileindex_t *index;
  unsigned char header = 1;
  unsigned char gzip, fmt;
  struct gzidxfile *gzf = NULL;
  
  filename = file->filename;
  gzip = file->gzip;
  fmt = file->fmt;
  index = file->index;

  noofseqs = index->noofchroms;

  if (gzip) {
    //gzindex = bl_zranGetIndex(filename, &gzlen);
    fp = fopen(filename, "rb");
    gzf = bl_initgzidxfile(fp, file->index->gzindex, 0, LARGECHUNK);
  } else {
    fp = fopen(filename, "r");
  }

  if(fp == NULL) {
    DBGEXIT("Couldn't open file %s. Exit forced!\n", filename);
  }

  buffer = ALLOCMEMORY(space, NULL, char, buffersize);
  map = ALLOCMEMORY(space, NULL, char*, noofseqs);
  memset(map, 0, sizeof(char*)*noofseqs);
  msz = ALLOCMEMORY(space, NULL, Uint, noofseqs);
  memset(msz, 0, sizeof(Uint)*noofseqs);

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
        curchrom = bl_matchfileGetChrom(fields, fmt);
        
        if(curend+1 > 0 && curchrom) {   

          for(id=0; id < index->noofchroms; id++) {
            if(strcmp(curchrom, index->chromnames[id]) == 0) {
              break;
            }

            for(j=0; j < strlen(index->chromnames[id]); j++) {
              if (isspace(index->chromnames[id][j])) break;
            }

            if(strlen(curchrom) <= j && 
                strncmp(curchrom, index->chromnames[id], j) == 0) {
              break;
            }    
          }

          if(id != lastid) {
            lastbin = -1;
            lastid = id;
            NFO("mapping chrom id:%d\n", id);
          }

          if (id >= index->noofchroms) {
            DBGEXIT("reference sequence '%s' not found\n", curchrom);
          }
          
          for(i=curstart; i < curend; i++) {      
            bin = i/255;
            if (bin != lastbin) {
              if(bin >= msz[id]) {
                map[id] = ALLOCMEMORY(space, map[id], unsigned char, bin+1);
                memset(&map[id][msz[id]], 0, sizeof(char)*((bin+1)-msz[id]));
                msz[id] = bin;
              }
            }
            if( ((Uint)map[id][bin])+1 > 0) map[id][bin]++;
            lastbin = bin;
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

  if(gzip) {
    bl_destructgzidxfile(gzf);
    FREEMEMORY(space, gzf);
  }

  FREEMEMORY(space, buffer);
  fclose(fp);

  *mapsize = msz;
  return map;
}

/*------------------ bl_evalmatchfileSampleCrossSections -------------------
 *    
 * @brief sample and execute f on it
 * @author Steve Hoffmann 
 *   
 */
 
Uint
bl_matchfileSampleCrossSections(void *space, 
    matchfile_t *file, fasta_t *set, Uint n, 
    void (*f)(void *, matchfileFrame_t*, Uint, 
      matchfileFrameStats_t *, void *), unsigned char **maps, Uint *mapsize, void *info)
{
  PairUint *samplepos;
  Uint i=0, r=0, j=0, k, *cumchrlen, 
       *order, prev=0, nchr, curchrom=0, curstart=0;
  matchfileFrame_t *frame = NULL;   
  matchfileFrameStats_t *stats = NULL;
  Uint maxcover = 20000;
  Uint setsize = 10000000;
  //Uint setsize = 100000;

  //init random number generator
  srand((unsigned)(time(0))); 
  nchr = file->index->noofchroms;

  samplepos = ALLOCMEMORY(space, NULL, PairUint, n+1);
  memset(samplepos, 0, sizeof(PairUint)*n+1);
  cumchrlen = ALLOCMEMORY(space, NULL, Uint, nchr);
  memset(cumchrlen, 0, sizeof(Uint)*nchr);
/*
  if(mymaps) {
    maps = mymaps;
  } else { 
    MSG("generating small map\n");
    //sum up the length of the references (i.e. chromosomes)
    maps = bl_matchfileSmallMap (space, file, &mapsize);
    MSG("map generated\n");
  }
*/
  cumchrlen[0] = file->index->matchend[0] - file->index->matchstart[0] + 1;
  for(i=1; i < nchr; i++) {
    assert(file->index->matchend[i] >= file->index->matchstart[i]);
    cumchrlen[i] = (file->index->matchend[i] - 
        file->index->matchstart[i]) + cumchrlen[i-1];
/*    NFO("chr %d: length: %d, sum: %d\n", i, (file->index->matchend[i] - 
        file->index->matchstart[i]), cumchrlen[i]);
        */
  }

//  fprintf(stderr, "printint small map\n");
//  for(i=0; i < nchr; i++) {
//    for(j=0; j < mapsize[i]; j++) {
//      fprintf(stdout, "chr %d\t%d\n", i, maps[i][j]);
//    }  
//  }
//  fprintf(stderr, "small map printed\n");


  //randomize n positions across the genome and deterimine their 
  //chromosomes
  i = 0;
  j = 0;
 
  while(i < n && j < setsize) {
     k=0;
     
     samplepos[i].b = 
       (Uint)(((double)cumchrlen[nchr-1]*rand()/RAND_MAX+1))+1;
     
     while(samplepos[i].b > cumchrlen[k] && k+1 < nchr) k++;   
     samplepos[i].a = k;
     
  //     if(k > 0) {
  //     fprintf(stderr, "violation: i=%d, samplepos[i].b=%d, k=%d, 
  //     nchr=%d, cumchrlen[k]=%d\n",
  //          i, samplepos[i].b, k, nchr, cumchrlen[k]);
  //     for(j=0; j < k; k++) {
  //       fprintf(stderr, "cumchrlen[%d]=%d\n", j, cumchrlen[j]);
  //     }
  //     exit(-1);
  //   }

     prev = (k == 0) ? 0 : cumchrlen[k-1];

     if(!maps || (maps[samplepos[i].a] 
        && mapsize[samplepos[i].a] > (samplepos[i].b - prev)/255 
        && maps[samplepos[i].a][(samplepos[i].b - prev)/255] > 200)) {
       i++;
       r++;
     }

     j++;
  }

  NFO("sampling %d positions.\n", i);

  /*
  for(i=0; i < nchr; i++) {
    if(maps[i]) { 
      FREEMEMORY(space, maps[i]);
    }
  }

  FREEMEMORY(space, maps);
  FREEMEMORY(space, mapsize);
 */

  if(j == setsize && r < (int)(0.8*((double)n))) {
    NFO("current sample size %d is below the minimum\n", r);
    /*FREEMEMORY(space, samplepos);
    FREEMEMORY(space, cumchrlen);
    return 0;*/
  }
  
//  for(i=0; i < n; i++) {
//    assert(samplepos[i].a < 1);
//  }

  //sort
  order = quickSort(space, samplepos, n, bl_matchfileSampleCmp, NULL);   

  initProgressBarVT();
 // for(i=0; i < n; i++) {
 //   assert(samplepos[order[i]].a < 1);
 // }
  
  //evaluate
  //to increase speed a frame of size FRAMESIZE is loaded
  for(i=0; i < n; i++) {
 
    progressBarVT("positions sampled.", n, i, 25);
    //is position on a new chromosome or in a new frame?
    if(!frame || samplepos[order[i]].a != curchrom || 
        samplepos[order[i]].b-prev+1 >= frame->start+frame->width) {

      if(frame) { 
        bl_matchfileDestructFrame(space, frame);
        frame = NULL;
        //bl_matchfileDestructFrameStats(space, stats);
      }

      curchrom = samplepos[order[i]].a;
      curstart = samplepos[order[i]].b;
      prev = (samplepos[order[i]].a == 0) ? 0 : cumchrlen[samplepos[order[i]].a-1];

   //   fprintf(stderr, "getting frame for '%s', curstart '%d', prev '%d'\n", 
   //       file->index->chromnames[samplepos[order[i]].a], curstart, prev);
      
      frame = bl_matchfileGetFrame(space, file, 
          file->index->chromnames[samplepos[order[i]].a], 
          curstart-prev+1, 20000, set, maxcover, NULL); 
     
   //   fprintf(stderr, "getting consensus\n" );
      bl_matchfileGetConsensus(frame);
     // stats = bl_matchfileFrameStats (space, frame);
    }
    
  //  fprintf(stderr, "evaluation of %d\n", samplepos[order[i]].b-curstart);
    f(space, frame, samplepos[order[i]].b-curstart, stats, info);   
  }
  fprintf(stderr, "\n");
  NFO("%d positions sampled.\n", n);
    
  if(frame) { 
    bl_matchfileDestructFrame(space,frame);
    frame = NULL;
    //    bl_matchfileDestructFrameStats(space, stats);
  }

  FREEMEMORY(space, order);
  FREEMEMORY(space, samplepos);
  FREEMEMORY(space, cumchrlen);
  
  return 0;
}


/*-------------------------------- phred2prob --------------------------------
 *    
 * @brief phred score to probability
 * @author Steve Hoffmann 
 *   
 */
 
double
phred2prob (char ascii, char offset, char p)
{

  double res, doffset = (double)offset, dascii = (double) ascii;
  res = pow(10,((double)(dascii-doffset)/-10.0));
  
  if(p && res >= 1.0) res = 0.99;

  return res;
}



/*------------------------ bl_matchfileTestNonVariant ------------------------
 *    
 * @brief test variant
 * @author Steve Hoffmann 
 *   
 */

  double
bl_matchfileTestNonVariant ( matchfileCross_t *cs, 
    matchfileSampleStats_t *stats, char ref, 
    matchfileCrossStats_t *css, double minrp, double maxrp, double minrq, 
    double maxrq)
{

  double p = .0, *nv, errors = .0; //, k;
  char *ch, *rq;
  unsigned char *ed;
  Uint len, i, rpos;
  uint32_t *rp, *mc, *rl;
  char usenativequal=1;
 
  len = cs->len;
  ch = cs->chars;
  rp = cs->readpos;
  rq = cs->quals;
  mc = cs->matchcnt;
  ed = cs->edist;
  rl = cs->readlen;
 
  css->mean_rt = .0;
  css->mean_rq = .0;
  css->mean_rr = .0;
  css->mean_mm = .0;

  nv = bl_matchfileGetNTReadPosVar(cs);

  bl_matchfileCrossStatsInit(css, len);

  for(i=0; i < len; i++) { 
    css->var_rt[i] = .0;
    css->var_rq[i] = .0;
    css->var_rr[i] = .0;
    css->var_mm[i] = .0;

    if(ntcode[(int)ch[i]] < 5 ) {

 
      rpos = trunc(((double)(((double)rp[i]*100.0)/((double)rl[i]))));

      if((int)ch[i] == (int)ref) {
    
        //READ POSITION
       //SUGGESTION RED
       css->var_rt[i] =   
        log(1.0-(((double)stats->RP[rpos] + 1.0)/((double) stats->RP_N[rpos] + 1.0))) - log(maxrp);

       //STANDARD 
       // css->var_rt[i] = log((((double)stats->RP[rpos] + 1.0)
       //     /((double) stats->RP_N[rpos] + 1.0))) ;
        
//ORIGINAL
//      css->var_rt[i] = log(1.0-(((double)stats->RP[rpos] + 1.0)
//      /((double) stats->RP_N[rpos] + 1.0))) 
//      - log(maxrp);

        //READ QUALITY
        if(usenativequal) { 
        css->var_rq[i] = log(1.0-phred2prob(rq[i], 64, 1));
        //css->var_rq[i] = log(1.0-pow(10,((double)((double)rq[i]-64.0)/-10.0)));
        } else { 
        css->var_rq[i] = log(1.0-(((double)stats->RQ[(int)rq[i]] + 1.0)
            /((double)stats->RQ_N[(int)rq[i]] + 1.0))) - log(maxrq);
        }

        //READ ERROR
        css->var_rr[i] = log(((double)stats->RR[MIN(ed[i],10)] + 1.0)/
            ((double)stats->RR_N + 1.0));

  
        //MULTIPLE MATCHES
        css->var_mm[i] = .0;
        if( stats->MM_N > stats->MM[MIN(mc[i],10)])
          css->var_mm[i] = MAX(MINMMPROB, 
              log(((double)stats->MM[MIN(mc[i],10)]+1.0)/
              ((double)stats->MM_N + 1.0)));

      } else {

        errors++;

        //READ POSITION
//SUGGESTION RED

       css->var_rt[i] =   
        log((((double)stats->RP[rpos] + 1.0)/((double) stats->RP_N[rpos] + 1.0))) - log(1-minrp);

//ORIGINAL
//        css->var_rt[i] = log(minrp) - log(1.0-(((double)stats->RP[rpos] + 1.0)
//            /((double) stats->RP_N[rpos] + 1.0)));

//STANDARD       
//        css->var_rt[i] =
//          log(1.0-(((double)stats->RP[rpos] + 1.0)
//            /((double) stats->RP_N[rpos] + 1.0)));

        //READ QUALITY
        if(usenativequal) { 
          css->var_rq[i] = log(phred2prob(rq[i], 64, 1));
//          css->var_rq[i] = log(pow(10,((double)((double)rq[i]-64.0)/-10.0)));
        } else { 
          css->var_rq[i] = log(minrq) - 
            log(1.0-(((double)stats->RQ[(int)rq[i]] + 1.0)
            /((double) stats->RQ_N[(int)rq[i]] + 1.0)));
        }

        //READ ERROR
        css->var_rr[i] = log(1.0-((double)stats->RR[MIN(ed[i],10)] + 1.0)/
            ((double)stats->RR_N + 1.0));

 
        //MULTIPLE MATCHES
        css->var_mm[i] = .0;
        if( stats->MM_N > stats->MM[MIN(mc[i],10)]) { 
          css->var_mm[i] = MAX(MINMMPROB, 
              log(1.0-(((double)stats->MM[MIN(mc[i],10)]+1.0)/
                  ((double)stats->MM_N + 1.0))));
        }
      }

      p += css->var_rt[i];
      p += css->var_rq[i];
      p += css->var_rr[i];
      p += css->var_mm[i];
     
      css->mean_rt += css->var_rt[i];
      css->mean_rq += css->var_rq[i];
      css->mean_rr += css->var_rr[i];
      css->mean_mm += css->var_mm[i];

      css->sub[i] = p;
    }   
  }

  css->mean_rt /= (double)len;
  css->mean_rq /= (double)len;
  css->mean_rr /= (double)len;
  css->mean_mm /= (double)len;

  FREEMEMORY(space, nv);
  return p/(double)len;
}


/*------------------------ bl_matchfileTestNonVariant ------------------------
 *    
 * @brief test variant
 * @author Steve Hoffmann 
 *   
 */
 
double
bl_matchfileTestVariant (matchfileCross_t *cs, 
    matchfileSampleStats_t *stats, char ref, matchfileCrossStats_t *css, 
    double minrp, double maxrp, double minrq, double maxrq)
{ 
 

  double px = .0, *nv; //, k;
  char *ch, *rq;
  unsigned char *ed;   
  Uint len, i, curerr, rpos;  
  uint32_t *rp, *mc, *rl;
  char usenativequal=1;


  len = cs->len;
  ch = cs->chars;
  rp = cs->readpos;
  rq = cs->quals;
  mc = cs->matchcnt;
  ed = cs->edist;
  nv = bl_matchfileGetNTReadPosVar(cs);
  rl = cs->readlen;
  
  bl_matchfileCrossStatsInit(css, len);
   
  css->mean_rt = .0;
  css->mean_rq = .0;
  css->mean_rr = .0;
  css->mean_mm = .0;

  /*variant*/
  for(i=0; i < len; i++) { 
   
   css->var_rt[i] = .0;
   css->var_rq[i] = .0;
   css->var_rr[i] = .0;
   css->var_mm[i] = .0;
   
   if(ntcode[(int)ch[i]] < 5 ) {

      rpos = trunc(((double)(((double)rp[i]*100.0)/((double)rl[i]))));
      //ALTERNATIVE RED
      css->var_rt[i] =   
        log(1.0-(((double)stats->RP[rpos] + 1.0)/((double) stats->RP_N[rpos] + 1.0))) - log(maxrp);

      //ORIGINAL
      //css->var_rt[i] = log(1.0-(((double)stats->RP[rpos] + 1.0)
      //      /((double) stats->RP_N[rpos] + 1.0))) - log(maxrp);
 
      //STANDARD
      //css->var_rt[i] = log((((double)stats->RP[rpos] + 1.0)
      //      /((double) stats->RP_N[rpos] + 1.0))) ;

      if(usenativequal) { 
        css->var_rq[i] = log(1.0-phred2prob(rq[i],64,1));
        //css->var_rq[i] = log(1.0-pow(10,((double)((double)rq[i]-64.0)/-10.0)));
      } else { 
        css->var_rq[i] = log(1.0-(((double)stats->RQ[(int)rq[i]] + 1.0)
              /((double)stats->RQ_N[(int)rq[i]] + 1.0))) - log(maxrq);
      }

      curerr = ((int)ch[i] != (int)ref && ed[i]);

      css->var_rr[i] = log(((double)stats->RR[MIN(ed[i]-curerr,10)]+1.0)/
          ((double)stats->RR_N + 1.0));
        

      if( stats->MM_N > stats->MM[MIN(mc[i],10)])
      css->var_mm[i] = MAX(MINMMPROB, log(((double)stats->MM[MIN(mc[i],10)]+1.0)/
          ((double)stats->MM_N + 1.0)));

      px += css->var_rt[i];
      px += css->var_rq[i];
      px += css->var_rr[i];
      px += css->var_mm[i];
     
      css->mean_rt += css->var_rt[i];
      css->mean_rq += css->var_rq[i];
      css->mean_rr += css->var_rr[i];
      css->mean_mm += css->var_mm[i];

      css->sub[i] = px; 
    }   
  }
  
  css->mean_rt /= (double)len;
  css->mean_rq /= (double)len;
  css->mean_rr /= (double)len;
  css->mean_mm /= (double)len;

  FREEMEMORY(space, nv);
  return px/(double)len;
}



/*--------------------------- bl_matchfileTestCons ---------------------------
 *    
 * @brief test the consensus vs. reference
 * @author Steve Hoffmann 
 *   
 */

double
bl_matchfileTestCons (matchfileCross_t *cs, 
    matchfileSampleStats_t *stats, char cons, matchfileCrossStats_t *css, 
    double minrp, double maxrp, double minrq, double maxrq, char usenativequal)
{

  double p = .0, *nv;
  char *ch, *rq;
  unsigned char *ed;
  Uint len, i, rpos;
  uint32_t *rp, *mc, *rl;

  len = cs->len;
  ch = cs->chars;
  rp = cs->readpos;
  rq = cs->quals;
  mc = cs->matchcnt;
  ed = cs->edist;
  rl = cs->readlen;
 
  css->mean_rt = .0;
  css->mean_rq = .0;
  css->mean_rr = .0;
  css->mean_mm = .0;

  nv = bl_matchfileGetNTReadPosVar(cs);
  bl_matchfileCrossStatsInit(css, len);

  for(i=0; i < len; i++) { 

    if(ntcode[(int)ch[i]] < 5 ) {

      rpos = trunc(((double)(((double)rp[i]*100.0)/((double)rl[i]))));

      if((int)ch[i] == (int)cons) {

  
        css->var_rt[i] =   
        log((((double)stats->RP[rpos] + 1.0)/((double) stats->RP_N[rpos] + 1.0))) - log(1-minrp);

/*        css->var_rt[i] = log(minrp) - 
          log(1.0-(((double)stats->RP[rpos] + 1.0)
            /((double) stats->RP_N[rpos] + 1.0)));
*/     
        if(usenativequal) { 
       // css->var_rq[i] = log(1.0-pow(10,phred2prob(rq[i],64,1)));
          css->var_rq[i] = log(phred2prob(rq[i], 64, 1));

        //css->var_rq[i] = log(1.0-pow(10,((double)((double)rq[i]-64.0)/-10.0)));
        } else { 
        css->var_rq[i] = log(minrq) - log(((double)stats->RQ[(int)rq[i]] + 1.0)
            /((double) stats->RQ_N[(int)rq[i]] + 1.0));
        }

        css->var_rq[i]*=QUALFACTOR;

        css->var_rr[i] = log(((double)stats->RR[MIN(ed[i],10)] + 1.0)/
            ((double)stats->RR_N + 1.0));

        css->var_mm[i] = .0;

        if( stats->MM_N > stats->MM[MIN(mc[i],10)])
          css->var_mm[i] = MAX(MINMMPROB, log(((double)stats->MM[MIN(mc[i],10)]+1.0)/
              ((double)stats->MM_N + 1.0)));

      } 

      p += css->var_rt[i];
      p += css->var_rq[i];
      p += css->var_rr[i];
      p += css->var_mm[i];
      css->sub[i] = p;

      css->mean_rt += css->var_rt[i];
      css->mean_rq += css->var_rq[i];
      css->mean_rr += css->var_rr[i];
      css->mean_mm += css->var_mm[i];
    }   
  }
  
  css->mean_rt /= (double)len;
  css->mean_rq /= (double)len;
  css->mean_rr /= (double)len;
  css->mean_mm /= (double)len;
 
  FREEMEMORY(space, nv);

  return p;
}


/*---------------------- bl_matchfileGetStandardization ----------------------
 *    
 * @brief get the standardization
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchfileGetStandardization (void *space, matchfileSampleStats_t *stats)
{
  Uint i;

  stats->MM[10] = (int)((double)stats->MM_N*0.0001)+1;
  stats->RR[10] = (int)((double)stats->RR_N*0.0001)+1;
  
  stats->minrp = 100;
  stats->maxrp = -100;
  stats->maxrq = -100;
  stats->minrq = 10;
  stats->currq = .0;

  for(i=0; i < 100; i++) {
    if(stats->RP[i] > 1) { 
      stats->maxrp = (stats->maxrp < 1.0-(((double)stats->RP[i]+1.0)/(stats->RP_N[i]+1.0))) ? 
        1.0-(((double)stats->RP[i]+1.0)/(stats->RP_N[i]+1.0)) : stats->maxrp;
  
      stats->minrp = (stats->minrp > 1.0-(((double)stats->RP[i]+1.0)/(stats->RP_N[i]+1.0))) ? 
        1.0-(((double)stats->RP[i]+1.0)/(stats->RP_N[i]+1.0)) : stats->minrp;
    }
  }

  for(i=0; i < 255; i++) {
    if(stats->RQ[i] > 10) { 
         
      stats->currq = (((double)stats->RQ[i]+1.0)/(stats->RQ_N[i]+1.0));
      
      stats->maxrq = (stats->maxrq < 1.0-stats->currq && stats->currq < 1.0) ? 
        1.0-stats->currq : stats->maxrq;
      stats->minrq = (stats->minrq > 1.0-stats->currq && stats->currq < 1.0) ? 
        1.0-stats->currq : stats->minrq;

    }
  }

  stats->standardized = 1;

  return ;
}

/*----------------------------- bl_matchfileTest -----------------------------
 *    
 * @brief test
 * @author Steve Hoffmann 
 *   
 */
 
Uint
bl_matchfileTest(void *space, Uint fidx, Uint cidx, Uint pos, 
    matchfileCross_t *cs, char ref, matchfileindex_t *idx, 
    unsigned char show, void *nfo)
{
  matchfileSampleStats_t *stats = idx->stats;
  matchfileCrossStats_t css; 
  Uint *cnt, *strands, conscnt, second, secondcnt, secondplus, secondminus, secondminimum;
    cs->s_cons = 1;
  cs->s_consx = 0;
  cs->s_ref = 1;
  cs->s_refx = 0;
  cs->p_hom = log(0);
  cs->scr_ref = log(0);
  cs->scr_cons = log(0);
  cs->scr_sample = log(0);

  if(!stats || cs->len < stats->mincover 
      || cs->len > stats->maxcover) {
    return 0;
  }
    
  cnt = bl_matchfileGetNTCounts(cs);
  strands = bl_matchfileGetPlusCounts(cs);

  cs->cons = (char) uarraymax(cnt, 255);
  conscnt = cnt[(int) cs->cons];
  second = (char) uarraysecond(cnt, 255, cs->cons);
  secondcnt = cnt[(int) second];
  assert(secondcnt <= conscnt);

  if(secondcnt == 0) { 
    if(cs->cons != ref) {
     secondcnt = conscnt;
     second = cs->cons;
    } else {
     FREEMEMORY(NULL, cnt);
     FREEMEMORY(NULL, strands);
     return 0;
    }
  } 
    
  secondplus = strands[(int) second];
  secondminus = secondcnt -secondplus;
  secondminimum = MIN(secondplus, secondminus);
  secondminimum = secondcnt- secondminimum;

  cs->ee = (double)secondcnt/(double)cs->len;

  FREEMEMORY(space, strands);
  FREEMEMORY(space, cnt);

  
  if(!stats->standardized) bl_matchfileGetStandardization (space, stats);

  cs->p_cons  = bl_matchfileTestNonVariant (cs, stats, cs->cons, &css, stats->minrp, 
      stats->maxrp, stats->minrq, stats->maxrq);

  bl_matchfileCrossStatsDestruct (&css);

  cs->p_consx = bl_matchfileTestVariant (cs, stats, cs->cons, &css, stats->minrp, 
      stats->maxrp, stats->minrq, stats->maxrq);
  
  bl_matchfileCrossStatsDestruct (&css);
   
  
  cs->p_refx = bl_matchfileTestVariant (cs, stats, ref, &css, stats->minrp, 
      stats->maxrp, stats->minrq, stats->maxrq);
  
  cs->diff_mm = 0;
  cs->diff_rt = css.mean_rt;
  cs->diff_rq = css.mean_rq;
  cs->diff_rr = css.mean_rr;
  cs->diff_mm = css.mean_mm;

  bl_matchfileCrossStatsDestruct (&css);
 
  cs->p_ref  = bl_matchfileTestNonVariant (cs, stats, ref, &css, stats->minrp, 
      stats->maxrp, stats->minrq, stats->maxrq);

//  cs->diff_rt -= css.mean_rt;
//  cs->diff_rq -= css.mean_rq;
//  cs->diff_rr -= css.mean_rr;
//  cs->diff_mm -= css.mean_mm;

  bl_matchfileCrossStatsDestruct (&css);
   
  if(cs->cons != ref) 
    cs->ee = (double)conscnt  /(double)cs->len;
  else
    cs->ee = (double)secondcnt/(double)cs->len;

  cs->pee = log(ecdf(cs->ee, stats->ecdf));
  
  cs->secondminimum= secondminimum;
  cs->secondcnt = secondcnt;

  if(stats->strand)
  cs->pbinom = log(pbinom(secondminimum, secondcnt, 0.5, 1));
  else
  cs->pbinom = 0;


  cs->scr_ref = (cs->p_refx - cs->p_ref) + log(ecdf(cs->ee, stats->ecdf)) 
    + cs->pbinom;
  
  cs->ee = (double)secondcnt/(double)cs->len;

  cs->scr_cons = (cs->p_consx - cs->p_cons) + log(ecdf(cs->ee, stats->ecdf)) 
    + cs->pbinom;

  cs->scr_sample = cs->scr_ref;
 
  //cs->scr_sample = (cs->p_refx - cs->p_ref) + log(ecdf(cs->ee, stats->ecdf)) 
  //  + cs->pbinom;



//  if (cs->cons != cs->ref && isinf(p_consx) && isinf(p_refx) 
//      && cs->p_cons-PX+P > cs->p_ref) {
//    cs->p_hom = cs->p_cons;
//  }

  return 0;
}



/*---------------------------- bl_matchfileSNPvcf ----------------------------
 *    
 * @brief get the VCF for a SNP
 * @author Steve Hoffmann 
 *   
 */

void
bl_matchfileSNPvcf (matchfileFrame_t *f, Uint p, char ref, Uint depth, Uint *cnt, Uint phred)
{
  Uint i;
  char id[] = ".", sep;
  char upper[] = {'A','C','G','T','N'};
  char* strings[] = {"A","C","G","T","N"};
  char lower[] = {'a','c','g','t','n'};
  char alt[50];
  char info[1000];
  Uint n = 0;
  Uint altcnt =0;
 
  memset(alt, 0, 50);
  memset(info, 0, 1000);

  n = snprintf(info, 1000, "DP=%d;AC=", depth);

  sep = 0;
  for(i=0; i < 5; i++) {
    if((cnt[(Uint)upper[i]] || cnt[(Uint)lower[i]]) && 
        (upper[i] != ref && lower[i] != ref)) { 
      if (sep == 1) {
        strcat (alt, ",");
        n+=snprintf(info+n, 1000-n, ",");
      }
      strcat(alt, strings[i]);
      altcnt += cnt[(Uint)upper[i]]+cnt[(Uint)lower[i]];
      n+=snprintf(info+n, 1000-n,"%d", cnt[(Uint)upper[i]]+cnt[(Uint)lower[i]]);
      sep=1;
    }
  }
    
  if(altcnt < cnt[(Uint)'-']) {
    n+=snprintf(info+n, 1000-n, ";INDEL");
  }

  printf("%s\t%d\t%s\t",f->chrname, f->start+p, id);
  printf("%c\t%s\t%d\tPASS\t%s\n",ref, alt, phred, info);
  

  return;
}

/*----------------------- bl_matchfileVariationHandler -----------------------
 *    
 * @brief handle variation output
 * @author Steve Hoffmann 
 *   
 */

indelvcf_t*
bl_matchfileVariationHandler (matchfileFrame_t *f, Uint p, double scr_ref, double scr_cons, double cut, 
    char ret, indelvcf_t *indelvcf)
{

 // char upper[] = {'A','C','G','T','N'};     
 // char lower[] = {'a','c','g','t','n'};
  char ref = f->ref[p];
  //Uint prev = f->start+p;
  //Uint nonrefsnps = 0;
  Uint depth = bl_matchfileGetCov(&f->cs[p], 1);
  double phred;
  double maxscr = MAX(scr_cons, scr_ref);
  Uint intphred;




  if(maxscr > cut || ret) { 
    phred = 10*(maxscr/log(10)); 
    intphred = phred;   
    
    Uint *cnt = bl_matchfileGetNTCounts(&f->cs[p]);
    /*
    for(i=0; i < 5; i++) {
      if((cnt[(Uint)upper[i]] || cnt[(Uint)lower[i]]) && 
          (upper[i] != ref && lower[i] != ref)) { 
        nonrefsnps++;
      }
    }

    //very simple check for deletion
    if(cnt[(Uint)'-'] > nonrefsnps) {
      if(indellen == 0) {   
        if(prev > 1) prev -= 2; //zero-offset
        indelvcf = ALLOCMEMORY(NULL, NULL, indelvcf_t, 1);
        indelvcf.ref = ALLOCMEMORY(NULL, indelvcf.ref, char, 2);
        indelvcf.ref[0] = f->refseq[prev];
        indelvcf.ref[1] = 0;
        indelvcf.len = 1;
        indelvcf.alleles = NULL;
        indelvcf.pos = prev+1; //one-offset
      }
      
      indelvcf.ref = ALLOCMEMORY(NULL, indelvcf.ref, char, indelvcf.len+2);
      indelvcf.ref[indelvcf.len] = ref;
      indelvcf.ref[indelvcf.len+1]=0;
      indelvcf.alleles = ALLOCMEMORY(NULL, indelvcf.alleles, char*, indelvcf.len);
      indelvcf.alleles[indelvcf.len-1] = ALLOCMEMORY(NULL, NULL, Uint, 5);
      indelvcf.phreds = ALLOCMEMORY(NULL, indelvcf.phreds, Uint, indelvcf.len);
      indelvcf.phreds[indelvcf.len-1] = intphred;
      for(i=0; i < 5; i++) {
        indelvcf.alleles[i] = cnt[(Uint)upper[i]] + cnt[(Uint)lower[i];
      }
      indelvcf.len += 1;
    
    } else {  
      
      if(indelvcf) {

        printf("%s\t%d\t%s\t",f->chrname, indelvcf.pos, id);
        printf("%c\t%c\t%d\tPASS\t%s\n", indelvcf.ref, indelvcf.ref[0], intphred, info);

        FREEMEMORY(NULL, indelvcf.ref);
        FREEMEMORY(NULL, indelvcf.phred);
        FREEMEMORY(NULL, indelvcf. alleles);
        FREEMEMORY(NULL, indelvcf);
        indelvcf = NULL;
      }
*/
      bl_matchfileSNPvcf (f, p, ref, depth, cnt, intphred);
 //   }

    FREEMEMORY(space, cnt);
  
/*  } else if(indelvcf) {
 *
    FREEMEMORY(NULL, indelvcf.ref);
    FREEMEMORY(NULL, indelvcf.phreds);
    FREEMEMORY(NULL, indelvcf.alleles);
    FREEMEMORY(NULL, indelvcf);
    indelvcf = NULL;*/
  } 

  return indelvcf;
//    return;
}


/*----------------------- bl_matchfileGroupTestsReset ------------------------
 *    
 * @brief reset group tests
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchfileTestGroupsReset (matchfileTestGroups_t *g)
{

  Uint i;
  for(i=0; i < g->noofgroups; i++) { 
    FREEMEMORY(space, g->cnt[i]);
  }
 
  memset(g->cnt, 0, sizeof(Uint*)*g->noofgroups);
  memset(g->p_cons, 0, sizeof(double)*g->noofgroups);
  memset(g->p_ref, 0, sizeof(double)*g->noofgroups);
  memset(g->p_consx, 0, sizeof(double)*g->noofgroups);
  memset(g->p_refx, 0, sizeof(double)*g->noofgroups);
  memset(g->scr_cons, 0, sizeof(double)*g->noofgroups);
  memset(g->scr_ref, 0, sizeof(double)*g->noofgroups);
  memset(g->type, 0, sizeof(char)*g->noofgroups);
	
  return ;
}

/*------------------------ bl_matchfileGroupTestsInit ------------------------
 *    
 * @brief init the group tests
 * @author Steve Hoffmann 
 *   
 */

  void
bl_matchfileTestGroupsInit (void *space, matchfileTestGroups_t *g, Uint n)
{
  g->noofgroups = n;
  g->p_cons = ALLOCMEMORY(space, NULL, double, g->noofgroups);
  g->p_ref = ALLOCMEMORY(space, NULL, double, g->noofgroups);
  g->p_consx = ALLOCMEMORY(space, NULL, double, g->noofgroups);
  g->p_refx = ALLOCMEMORY(space, NULL, double, g->noofgroups);
  g->scr_cons = ALLOCMEMORY(space, NULL, double, g->noofgroups);
  g->scr_ref = ALLOCMEMORY(space, NULL, double, g->noofgroups);
  g->cnt = ALLOCMEMORY(space, NULL, Uint, g->noofgroups);
  g->type = ALLOCMEMORY(sapce, NULL, char, g->noofgroups);
  
  bl_matchfileTestGroupsReset(g);

  return ;
}


/*---------------------- bl_matchfileTestGroupsDestruct ----------------------
 *    
 * @brief destruct the groups
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchfileTestGroupsDestruct (void *space, matchfileTestGroups_t *g)
{
  Uint i;
  for(i=0; i < g->noofgroups; i++) { 
    FREEMEMORY(space, g->cnt[i]);
  }
 
  FREEMEMORY(space, g->scr_cons);
  FREEMEMORY(space, g->scr_ref);
  FREEMEMORY(space, g->cnt);
  FREEMEMORY(space, g->type);
  FREEMEMORY(space, g->p_cons);
  FREEMEMORY(space, g->p_ref);
  FREEMEMORY(sapce, g->p_consx);
  FREEMEMORY(space, g->p_refx);


  return ;
}


/*------------------------ bl_matchfileGroupAddResult ------------------------
 *    
 * @brief add a test result to a group
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchfileTestGroupsAddResult (matchfileTestGroups_t *g, Uint no, 
    matchfileCross_t *cs)
{
  
  Uint *cnt = bl_matchfileGetNTCounts(cs);
  
  g->scr_cons[no] = cs->scr_cons;
  g->scr_ref[no] = cs->scr_ref;
  g->p_cons[no] = cs->p_cons;
  g->p_ref[no] = cs->p_ref;
  g->p_consx[no] = cs->p_consx;
  g->p_refx[no] = cs->p_refx;
  g->cnt[no] = cnt;


  return ;
}



/*----------------------- bl_matchfileTestGroupsPrint ------------------------
 *    
 * @brief print group tests
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchfileTestGroupsPrint (matchfileTestGroups_t *g, Uint no, 
    matchfileFrame_t *f, Uint p)
{
 
  return ;
}

/*---------------------- bl_matchfileEvalCrossSections -----------------------
 *    
 * @brief evaluate all cross sections
 * @author Steve Hoffmann 
 *   
 */


void
bl_matchfileEvalCrossSections (void *space,  matchfile_t **files, 
    int *groups, Uint nooffiles, fasta_t *set, Uint framesize,
    Uint (*f)(void *, Uint fidx, Uint cidx, Uint pos, matchfileCross_t*, char ref, 
    matchfileindex_t *, unsigned char, void *), void *nfo)
{
  Uint i, j, k, ret, n=0, nchr = 0, curchrom = 0, pos=0, maxgroupno=0, chridx, maxchr;
  matchfileFrame_t **frames = NULL;
  matchfileTestGroups_t *g = NULL;
  //unsigned char exclusive = 0;
  Uint maxcover = 20000;
  indelvcf_t *indel;

  for(j=0; j < nooffiles; j++) {
    nchr = MAX(nchr, files[j]->index->noofchroms);
  }

  if(groups) {    
      
    g = ALLOCMEMORY(space, NULL, matchfileTestGroups_t, 1);
    for(i=0; i < nooffiles; i++) {
      maxgroupno = (groups[i]  > maxgroupno) ? groups[i] : maxgroupno;
    }
    maxgroupno+=1;
  }

  frames = ALLOCMEMORY(space, NULL, matchfileFrame_t*, nooffiles);
  memset(frames, 0, sizeof(matchfileFrame_t*)*nooffiles);

  for(k=0; k < nchr; k++) {  

    chridx = bl_fastxFindIDIdx(files[0]->index->chromnames[k], set);
    maxchr = bl_fastaGetSequenceLength(set, chridx);

    for(j=0, n=0; j < nooffiles; j++) {  
      n = MAX(n, files[j]->index->matchend[k]+1);
   }

    NFO("evaluating chr: %d '%s' len:%d\n", chridx, files[0]->index->chromnames[k], maxchr);
    n = MIN(n, maxchr);

    //evaluate
    //to increase speed a frame of size FRAMESIZE is loaded

    for(i=1; i < n; i++) {

      for(j=0; j < nooffiles; j++) {

        if(files[j]->index->matchend[k]>0) { 
          
          if(groups) { 
            bl_matchfileTestGroupsInit(space, g, maxgroupno);
          }

          //is position on a new chromosome or in a new frame?
          if(!frames[j] || k != curchrom || 
              i >= frames[j]->start + frames[j]->width) {

            curchrom = k;

            if(frames[j]) {
              bl_matchfileDestructFrame(space, frames[j]);
              frames[j] = NULL;
            }

            if(files[j]->index->matchend[k] < i ||
                files[j]->index->matchstart[k] > i)
              continue;

            frames[j] = bl_matchfileGetFrame(space, files[j],
                files[j]->index->chromnames[k], i, framesize, set, maxcover, NULL);
            bl_matchfileGetConsensus(frames[j]);

          }

          pos = i - frames[j]->start;
          frames[j]->cs[pos].p_hom = log(0);
          ret = 0;

          if(frames[j]->cs[pos].len) {          
            ret = f(space, j, k, i, &frames[j]->cs[pos], frames[j]->ref[pos], 
                files[j]->index, 0, nfo);   
          }

          double scr_ref  = frames[j]->cs[pos].scr_ref;
          double scr_cons = frames[j]->cs[pos].scr_cons;
          double cut = files[j]->index->stats->cut;

//          if(ret != 0 || 
//              MAX(frames[j]->cs[pos].scr_ref, frames[j]->cs[pos].scr_cons) > 
//              files[j]->index->stats->cut) {

            if(!groups) { 
              //determine type of variation
              indel= bl_matchfileVariationHandler (frames[j], pos, scr_ref, scr_cons, cut, ret, indel);
            }  

//            ccnt++;
//          }

          if(groups && frames[j]->cs[pos].len){ 
            bl_matchfileTestGroupsAddResult (g, groups[j], &frames[j]->cs[pos]);
          }
        }

        //    if(groups) { 

        /*iter the groups*/
        //      for(j=0; j < maxgroupno; j++) {
        //        if(g->s_consx[j] > g->s_cons[j] || 
        //            g->s_refx[j] > g->s_ref[j]) {
        //          exclusive = 1;

        //         for(u=0; u < maxgroupno; u++) {
        //           if(u != j && (g->s_consx[u] > g->s_cons[u] || 
        //                 g->s_refx[u] > g->s_ref[u])) {
        //             exclusive = 0;
        //           }
        //         }

        //         if(exclusive) { 
        //           bl_matchfileTestGroupsPrint (g, j, frames[j], pos);
        //         }
        //       }
        //     }

        //     bl_matchfileTestGroupsDestruct (space, g);
        //   }
      }
    }
  
    for(j=0; j < nooffiles; j++){
        if(frames[j]) bl_matchfileDestructFrame(space,frames[j]);
        frames[j] = NULL;
    }
  }


  if(groups) {
    FREEMEMORY(space, g);
  }

  FREEMEMORY(space, frames);
  return ;
}





/*-------------------------- bl_matchfileConsensus --------------------------
 *    
 * @brief a simple consensus calling based on majority voting,
 *        considering deletions as well (not done by bl_matchfileGetConsensus)
 * @author Christian Otto
 *
 */
Uint
bl_matchfileConsensus ( void *space, Uint fidx, Uint cidx, Uint pos, matchfileCross_t *cs, 
                        char ref, matchfileindex_t *idx, unsigned char show, void *nfo )
{
  Uint i, j, k, max, *cnt = NULL, len = 0;
  char *del, **dels = NULL;

  matchfile_t **files = (matchfile_t **) nfo;

  /* compile deletion string freq table */
  if (cs->noofdels > 0){
    for (i = 0; i < cs->noofdels; i++){
      del = ALLOCMEMORY(space, NULL, char, cs->dels[i].len + 1);
      for (j = 0, k = 0; j < cs->dels[i].len; j++){
        if (cs->dels[i].string[j] != '^')
          del[k++] = cs->dels[i].string[j];
      }
      del[k] = '\0';

      for (j = 0; j < len; j++){
        if (strlen(dels[j]) == strlen(del) &&
            strcmp(dels[j], del) == 0){
          cnt[j]++;
        }            
      }
      if (j == len){
        dels = ALLOCMEMORY(space, dels, char*, len+1);
        cnt = ALLOCMEMORY(space, cnt, Uint, len+1);
        
        dels[len] = ALLOCMEMORY(space, NULL, char, strlen(del)+1);
        memmove(dels[len], del, strlen(del));
        dels[len][strlen(del)] = '\0';        
        cnt[len] = 1;
        len++;
      }
      FREEMEMORY(space, del);
    }

    if (len > 0){
      max = uarraymax(cnt, len);
      
      if (2 * cnt[max] >= cs->len){
        fprintf(stdout, "%s\t%d\t-\t%s\t%d\n", files[fidx]->index->chromnames[cidx], 
                pos, dels[max], cs->len);
      }
    }
  }
  fprintf(stdout, "%s\t%d\t%c\t%c\t%d\n", files[fidx]->index->chromnames[cidx], 
          pos, cs->ref, cs->cons, cs->len);
  
  if (cs->noofdels > 0){
    for (i = 0; i < len; i++){
      FREEMEMORY(space, dels[i]);
    }
    FREEMEMORY(space, dels);
    FREEMEMORY(space, cnt);
  }
  return 0;
}

/*----------------------------- bl_matchfileGetCov -----------------------------
 *    
 * @brief  get coverage in cross section w/ or w/o considering deleted bases
 * @author Christian Otto
 *
 */
Uint
bl_matchfileGetCov(matchfileCross_t *cs, unsigned char allowdel){
  Uint cov, *cnt;
  
  cov = cs->len;
  cnt = bl_matchfileGetNTCounts(cs);
  if (!allowdel){
    cov -= cnt[(Uint)'-'];
  }
  
  free(cnt);  
  return cov;
}

/*--------------------------- bl_matchfileGetContext ---------------------------
 *    
 * @brief  get sequence context of given length w/r/t to given strand
 * @author Christian Otto
 *
 */
char *
bl_matchfileGetContext(fasta_t *fasta, Uint idx, Uint pos, Uint strand, Uint len){
  Uint i, seqlen;
  char *seq, *context;

  context = ALLOCMEMORY(space, NULL, char, len+1);
  memset(context, 'N', len);
  context[len] = '\0';

  if (!strand || idx >= fasta->noofseqs) 
    return context;

  seq = bl_fastaGetSequence(fasta, idx);
  seqlen = bl_fastaGetSequenceLength(fasta, idx);

  if (pos > seqlen)
    return context;

  for (i = 0; i < len; i++){
    if (strand == PLUS_STRAND && pos+i-1 < seqlen){
      context[i] = seq[pos+i-1];
    }
    else if (strand == MINUS_STRAND && pos >= i+1) {
      context[i] = charComplementChar(seq[pos-i-1]);
    }
  }
  return context;
}

