/*
 * manout.c
 * attempt for flexible output of genome mapping w/ SEGEMEHL
 *
 * @author Christian Otto
 * @email christan@bioinf.uni-leipzig.de
 * @date Wed Sep 24 10:56:23 CEST 2008
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <string.h>
#include "basic-types.h"
#include "bitArray.h"
#include "memory.h"
#include "mathematics.h"
#include "sort.h"
#include "info.h"
#include "biofiles.h"
#include "fileio.h"
#include "vtprogressbar.h"
#include "debug.h"
#include "charsequence.h"
#include "manout.h"
#include <assert.h>
#include <pthread.h>
#include "alignment.h"
#include "manoutformats.h"
#include "kdseed.h"
#include "fileBins.h"
#include "segemehl.h"

unsigned char
se_kdMatchListhasMatches(gmatchlist_t *list) {
  return (list->n[0] > 0 || list->n[1] > 0);
}

unsigned char
se_kdMatchListhasMates(gmatchlist_t *list) {
  Uint u,i,k;
  for(u=0; u < 2; u++) {
    for(i=0; i < list->n[u]; i++) {
      for(k=0; k < 4; k++) {
        if(list->matches[u][i].mates[k].al)
          return 1;
      }
    }
  }
  return 0;
}

unsigned char
se_kdMatchHasMates(gmatch_t *match) {
  Uint i;
  for(i=0; i < 4; i++) {
    if (match->mates[i].isset) return 1;
  }
  return 0;
}

void
se_destructMatches(void *space, gread_t *read) {
  FREEMEMORY(space, read->matches[0]);
  FREEMEMORY(space, read->matches[1]);
}

Uint
se_setMatches (void *space, gread_t *read, 
    gmatchlist_t *list, Uint maxedist, segemehl_t *nfo, char rep) {

  Uint i, k, n, noofmatepairs=0, matches=0;
  gmatch_t *m = NULL, *cur;

  read->noofmatepairs = 0;
  read->noofmatches = 0;

  for(k=0; k <= 1; k++) {
    m = NULL;
    n = 0;

//    if(rep) fprintf(stdout, "setMatches: iter list of %d matches\n", list->n[k]);
    for(i=0; i < list->n[k]; i++) {

      cur = &list->matches[k][i];
      if((cur->edist <= maxedist || 
            (se_kdMatchHasMates(cur) && cur->mateminedist+cur->edist <= list->pairminedist)) 
          && !cur->skip) {

//        if(rep) fprintf(stdout, "setMatches: match %d with edist %d < %d maxedist\n", i, cur->edist, maxedist);
        if(se_kdMatchHasMates(cur)) {
          //matemindist of match compared to matemindist of list!
          //if(cur->mateminedist <= cur->mateminedist || !nfo->bestonly) {
//          if(rep) fprintf(stdout, "setMatches: read has Mates\n");
          if(cur->mateminedist+cur->edist <= list->pairminedist || !nfo->bestonly) { 

//            if(rep) fprintf(stdout, "setMatches: added with combinded edist %d+%d=%d <= %d maxedist\n", cur->edist, cur->mateminedist, cur->edist+cur->mateminedist, list->pairminedist);
            noofmatepairs++;
            m = ALLOCMEMORY(space, m, gmatch_t, n+1);
            memmove(&m[n], cur, sizeof(gmatch_t));
            n++;
            matches++;
          } else {

//            if(rep) fprintf(stdout, "setMatches: failed to add because combinded edist %d with edist %d > %d maxedist\n", i, cur->edist+cur->mateminedist, list->pairminedist);
          }
        } else {
          m = ALLOCMEMORY(space, m, gmatch_t, n+1);
          memmove(&m[n], cur, sizeof(gmatch_t));
          n++;
          matches++;
        }
        } else {

//          if(rep) fprintf(stdout, "setMatches: match %d failed because %d > %d or skip:%d\n", i, cur->edist, maxedist, cur->skip);
        }    
      }

      read->noofmatches += n;
      read->n[k] = n;
      read->matches[k] = m;
    }

    read->noofmatepairs = noofmatepairs;
    return matches;
  }


gmatchlist_t*
se_kdMatchListAdd(gmatchlist_t *list, 
    Uint chr_idx, 
    Uint chr_start, 
    Uint chr_end, 
    Uint edist,
    int scr,
    Uint start,
    Uint end,
    double evalue, Alignment *al, Uint u, 
    Uint previdx, Uint prevpos, char prevstrand, 
    Uint nextidx, Uint nextpos, char nextstrand, Uint fragno) {

  Uint n, mat, mis, del, ins, i;
  gmatch_t *match;

  n = list->n[u];
  list->matches[u]=realloc(list->matches[u], sizeof(gmatch_t)*(n+1));
  match = &list->matches[u][n];
  
  initMatch(match);
  countEops(al, &mat, &mis, &ins, &del);
  //changed edist to mis + ins + del
  list->minedist = MIN(list->minedist, mis+ins+del);
  match->scr = scr; 
  match->evalue = evalue;
  //changed edist to mis + ins + del
  match->edist = mis+ins+del;
  match->p = chr_start;
  match->q = chr_end;
  match->i = start; 
  match->j = end; 
  match->mat = mat;
  match->mis = mis;
  match->ins = ins;
  match->del = del;
  match->subject = chr_idx;
  match->al = al;
  match->noofmatematches = 0;
  match->skip = 0;
  match->mateminedist = -1;

  for(i=0; i < 4; i++) {
    match->mates[i].scr = 0;
    match->mates[i].al = NULL;
    match->mates[i].isset = 0;
    match->mates[i].materefdesc = NULL;
    match->mates[i].materefseq = NULL;
    match->mates[i].materefdesclen = 0;
  }

  match->previdx = previdx;
  match->prevpos = prevpos;
  match->nextidx = nextidx;
  match->nextpos = nextpos;
  match->fragno = fragno;
  match->prevflags = 0;
  match->nextflags = 0;
  
  match->prevseqstart = 0;
  match->prevseqrefdesc = NULL;
  match->nextseqstart = 0;
  match->nextseqrefdesc = NULL;
  match->refdesc = NULL;
  match->refdesclen = 0;
  match->refseq = NULL;

  if(prevstrand == '+') {
    match->prevflags |= SPLIT_PREV_PLUS;
  }
  if(nextstrand == '+') {
    match->nextflags |= SPLIT_NEXT_PLUS;
  }

  list->n[u]++;
  return list;
}


gmatchlist_t*
se_kdMatchListSet(void *space,
    gmatchlist_t *list, 
    Uint chr_idx, 
    Uint chr_start, 
    Uint chr_end, 
    Uint edist,
    int scr,
    Uint start,
    Uint end,
    double evalue, Alignment *al, Uint u, Uint n) {

  Uint mat, mis, del, ins, i;
  gmatch_t *match;

  if (n < 0 || n >= list->n[u]){
    return list;
  }

  match = &list->matches[u][n];
  countEops(al, &mat, &mis, &ins, &del);
  list->minedist = MIN(list->minedist, (mis+ins+del));
  match->scr = scr; 
  match->evalue = evalue;
  //changed edist to mis + ins + del
  match->edist = mis+ins+del;
  match->p = chr_start;
  match->q = chr_end;
  match->i = start; 
  match->j = end; 
  match->mat = mat;
  match->mis = mis;
  match->ins = ins;
  match->del = del;
  match->subject = chr_idx;

  if (match->al){
    wrapAlignment(match->al);
    FREEMEMORY(space, match->al);
  }

  match->al = al; 
  match->noofmatematches = 0;
  match->skip = 0;
  match->mateminedist = -1;

  for(i=0; i < 4; i++) {
    if(match->mates[i].al != NULL) {
      wrapAlignment(match->mates[i].al);
    }
    match->mates[i].scr = 0;
    match->mates[i].al = NULL;
    match->mates[i].isset = 0;
  }

  match->fragno = 0;
  match->previdx = -1;
  match->prevpos = -1;
  match->nextidx = -1;
  match->nextpos = -1;

  return list;
}


/*------------------------------- se_kdSetMate -------------------------------
 *    
 * @brief set mate
 * @author Steve Hoffmann 
 *   
 */
 
Uint
se_kdSetMate(void *space, gmatch_t *match, 
    Uint chr_idx, Uint chr_start, Uint chr_end, Uint edist,
    Alignment *al, unsigned char downstream, unsigned char rc) 
{  
  int scr=0;
  unsigned char pos=0;
  gmate_t *mates;
  Uint mat, mis, del, ins;


  assert (downstream < 2 && rc < 2);
  pos = (downstream << 1) | rc;
  countEops(al, &mat, &mis, &ins, &del);

  scr = mat;
  scr -= mis+ins+del;

  mates = match->mates;

  if(mates[pos].isset && mates[pos].al != NULL && mates[pos].scr != 0) {
    if(scr > mates[pos].scr) {
      wrapAlignment(mates[pos].al);
    } else {
      wrapAlignment(al);
      FREEMEMORY(space, al);
      return match->mateminedist;
    }
  }

  mates[pos].isset = 1;
  mates[pos].p = chr_start;
  mates[pos].q = chr_end;
  mates[pos].scr = scr;
  mates[pos].mat = mat;
  mates[pos].mis = mis;
  mates[pos].ins = ins;
  mates[pos].del = del;
  mates[pos].subject = chr_idx;
  //changed edist to mis + ins + del
  mates[pos].edist = mis+ins+del;
  mates[pos].al =  al;

  if((mis+ins+del) < match->mateminedist || !match->noofmatematches) {
//    match->mateminedist = edist;            
//   fprintf(stdout, "mateminedist CHECKPOINT e: %d\n", match->mateminedist);
   match->mateminedist = mis+ins+del;
  }

  return match->mateminedist;
}


/*-------------------------------- initStruct ---------------------------------
 *    
 * @brief inits
 * @author Steve Hoffmann 
 *   
 */

  inline void
initMatch (gmatch_t* m)
{
  m->i = 0; 
  m->p = 0;
  m->q = 0;
  m->scr = 0;
  m->mat = 0;
  m->mis = 0;
  m->del = 0;
  m->ins = 0;
  m->subject = 0;
  m->edist = 0;
  m->al = NULL;
  m->skip = 0;

  memset(&m->mates, 0, sizeof(gmate_t)*4);

  m->fragno = 0;
  m->previdx = -1;
  m->prevpos = -1;
  m->prevflags = 0;
  m->nextidx = -1;
  m->nextpos = -1;
  m->nextflags = 0;

  return ;
}

  void 
initRead(gread_t *r, Uint readid) 
{
  r->id = readid;
  r->n[MINUSSTRAND] = 0;
  r->n[PLUSSTRAND] = 0;
  r->matches[MINUSSTRAND] = NULL;
  r->matches[PLUSSTRAND] = NULL;

  return;
}

  void
initGmap(Gmap *map, MultiCharSeq *mseq, Uint offset) 
{ 
  map->mseq = mseq;
  map->mapoffset = offset;
  map->noofreads = 0;
  map->reads = 0;

  return;
}



inline void
setReads(Gmap *map, gread_t *r, Uint noofreads){
  map->reads = r;
  map->noofreads = noofreads;
}



/*---------------------------- bl_gmatchlistInit -----------------------------
 *    
 * @brief init gmatchlist
 * @author Steve Hoffmann 
 *   
 */

gmatchlist_t*
bl_gmatchlistInit(void *space, int maxedist, int matemaxedist) {

  gmatchlist_t* list;
  list = ALLOCMEMORY(space, NULL, gmatchlist_t, 1);
  list->matches = ALLOCMEMORY(space, NULL, gmatch_t*, 2);
  list->n = ALLOCMEMORY(space, NULL, Uint, 2);
  list->matches[0] = NULL;
  list->matches[1] = NULL;
  list->n[0] = 0;
  list->n[1] = 0;
  list->minedist = maxedist;
  list->mateminedist = matemaxedist;
  list->pairminedist = maxedist+matemaxedist;
 
  return list;

}

/*-------------------------- bl_gmatchlistDestruct ---------------------------
 *    
 * @brief destruct gmatchlist
 * @author Steve Hoffmann 
 *   
 */

void
bl_gmatchlistDestruct(void *space, gmatchlist_t *list){
  Uint i, k, u;

  for(u=0; u < 2; u++) {
    if (list->matches[u]) {
      for(i=0; i < list->n[u]; i++) {
        for(k=0; k < 4; k++) {
          if(list->matches[u][i].mates[k].al) {
            wrapAlignment(list->matches[u][i].mates[k].al);
            FREEMEMORY(space, list->matches[u][i].mates[k].al);
          }
        } 
        if(list->matches[u][i].al) {
          wrapAlignment(list->matches[u][i].al);
          FREEMEMORY(space, list->matches[u][i].al);
          list->matches[u][i].al = NULL;
        }
      }
      FREEMEMORY(space, list->matches[u]);
    }
  }

  FREEMEMORY(space, list->matches);
  FREEMEMORY(space, list->n);
  FREEMEMORY(space, list);
}


/*----------------------------- genericOutput --------------------------------
 *  
 *  @brief write generic output from given string list according to format
 *  @author Christian Otto
 *     
 */

void
genericOutput (FILE *dev, char **list, Uint rep_type, char lf){
  Uint i = 0;
  while(FORMAT[rep_type][i] > -1){
    if (i > 0){
      fprintf(dev, "%s", SEPARATOR);
    }
    fprintf(dev, "%s", list[FORMAT[rep_type][i]]);
    i++;
  }
  fprintf(dev, "%c", lf);
}


/*---------------------------- reportSplicedMatch ----------------------------
 *    
 * @brief reporting spliced matches
 * @author Steve Hoffmann 
 *   
 */

void
reportSplicedMatch(void *space, char* qrydesc, MultiCharSeqAlignment *mcsa, 
    Uint noofaligns, Uint coverage, Uint edist, int score, segemehl_t *nfo) {

  Uint i, j, ulen, vlen, ustart, vstart;
  FILE *dev = NULL;
  bl_fileBin_t *fx = NULL;
  char strands[]= {'+','-'};
  char alignlf= '\n';

  for(j=0; j < noofaligns; j++) {

    ulen = getUalignlen(mcsa[j].al); 
    vlen = getValignlen(mcsa[j].al);
    ustart = mcsa[j].al->uoff;
    vstart = mcsa[j].al->voff;

    if (nfo->splitdev != NULL) {
      dev = nfo->splitdev;
      if (nfo->threadno > 1) {
        pthread_mutex_lock(nfo->mtx3);
      }
    } else {
      fx = bl_fileBinsDomainGetBin(nfo->splitbins, mcsa[j].refdesc, 
          mcsa[j].refstart - mcsa[j].substart + vstart + 1);
      bl_fileBinsLock(fx);
      dev=bl_fileBinsOpen(space, fx, "w");    
    }

    fprintf(dev, "%s\t%d\t%d\t%d\t", qrydesc, edist, coverage, noofaligns);
    fprintf(dev, "%d\t%u\t%u\t%u\t%u\t%u\t%c\t%s\t", 
        j,
        getEdist(mcsa[j].al),
        ustart, ustart+ulen,
        mcsa[j].refstart - mcsa[j].substart + vstart +1, 
        mcsa[j].refstart - mcsa[j].substart + vstart +vlen, 
        strands[mcsa[j].strand], 
        mcsa[j].refdesc);

    for (i = 0; i < noofaligns; i++) {
      ulen = getUalignlen(mcsa[i].al); 
      vlen = getValignlen(mcsa[i].al);
      ustart = mcsa[i].al->uoff;
      vstart = mcsa[i].al->voff;

      if (mcsa[i].strand == 1) {
        ustart = mcsa[i].qrylen-ustart-ulen;
      }

      fprintf(dev, "%u\t%u\t%u\t%u\t%u\t%c\t%s\t", 
          getEdist(mcsa[i].al),
          ustart, ustart+ulen,
          mcsa[i].refstart - mcsa[i].substart + vstart +1, 
          mcsa[i].refstart - mcsa[i].substart + vstart +vlen, 
          strands[mcsa[i].strand], 
          mcsa[i].refdesc);
    }


    if (nfo->align){
      if(nfo->order) {
        alignlf = 7;
        fprintf(dev, "%c", alignlf);
      } else {
        fprintf(dev, "\n");
      }       
      for(i = 0; i < noofaligns; i++) {
        showAlignLF(mcsa[i].al, dev, alignlf);
        fprintf(dev, "%c", alignlf);
      }
    }

    fprintf(dev,"\n");
    fflush(dev);

    if (nfo->splitdev && nfo->threadno > 1) { 
      pthread_mutex_unlock(nfo->mtx3); 
    } 

    if(nfo->splitbins) {
      bl_fileBinsUnlock(fx);
    }  
  }

  return ;
}


/*-------------------------------- getDevice ---------------------------------
 *    
 * @brief get the device
 * @author Steve Hoffmann 
 *   
 */

FILE*
getDevice (void *space, char *chr, Uint pos, bl_fileBin_t **fx, segemehl_t *nfo)
{
  FILE *dev;
  char *bisulfite;

  if(!nfo->bins) {
    dev = nfo->dev;
  } else {
    if (!nfo->bisulfitemerging){
      *fx = bl_fileBinsDomainGetBin(nfo->bins, chr, pos);
      bl_fileBinsLock(*fx);
      dev=bl_fileBinsOpen(space, *fx, "w");    
    } else {
      bisulfite = calloc(MAX_INT_LENGTH+1, 1);
      sprintf(bisulfite, "%u", (nfo->bisulfite + 1) % 2);
      *fx = bl_fileBinsDomainGetBin(nfo->bins, bisulfite, nfo->threadid);
      //DBG("info.bisulfite=%u\tbisulfite=%s\tthreadid=%u\tfilename=%s\tstart=%llu\tend=%llu\n",
      //    nfo->bisulfite, bisulfite, nfo->threadid, (*fx)->fname, (*fx)->id->start, (*fx)->id->end);
      bl_fileBinsLock(*fx);
      dev=bl_fileBinsOpen(space, *fx, "w");
      free(bisulfite);
    }
  } 

  return dev;
}


/*-------------------------------- getSAMTags --------------------------------
 *    
 * @brief get SAM tags
 * @author Steve Hoffmann 
 *   
 */
 

char*
getSAMTags(Gmap *map, Uint ustart, Uint uend, Uint uno, Alignment *al, 
    Uint edist, Uint previdx, Uint prevpos, char prevflags, 
           Uint nextidx, Uint nextpos, char nextflags, Uint matchid,
    Uint noofmatches, Uint noofsplits, char pStatChr, char *meopstr, segemehl_t *nfo, gmatch_t *match) {

  char *tag,
       *md,
       *refdesc=NULL;
  Uint ptr=0,
    len, mis,
       prevseqstart = 0,
       nextseqstart = 0;

  md = mdstring(al, 0);
 

  if(previdx != -1 || nextidx != -1 ) {
    //only one split alignment per read 
    len = snprintf(NULL, 0, "NM:i:%d\tMD:Z:%s\tNH:i:1\tXI:i:%d\tXL:i:%d", edist, md, matchid, noofsplits);
    tag = ALLOCMEMORY(space, NULL, char, len+1);
    snprintf(tag, len+1, "NM:i:%d\tMD:Z:%s\tNH:i:1\tXI:i:%d\tXL:i:%d", edist, md, matchid, noofsplits);
    ptr += len;

  } else { 
 
  len = snprintf(NULL, 0, "NM:i:%d\tMD:Z:%s\tNH:i:%d\tXI:i:%d", edist, md, noofmatches, matchid);
  tag = ALLOCMEMORY(space, NULL, char, len+1);
  snprintf(tag, len+1, "NM:i:%d\tMD:Z:%s\tNH:i:%d\tXI:i:%d", edist, md, noofmatches, matchid);
  ptr += len;
  }

  if (nfo->SAMmeop) {
    len = snprintf(NULL, 0, "\tXE:Z:%s", meopstr);
    tag = ALLOCMEMORY(space, tag, char, ptr+len+1);
    snprintf(&tag[ptr], len+1, "\tXE:Z:%s", meopstr);	  
    ptr += len;
  }

  if(nfo->SAMpairstat) {
    tag = ALLOCMEMORY(space, tag, char, ptr+8);
    snprintf(&tag[ptr], 8, "\tXA:Z:%c", pStatChr);	  
    ptr += 7;
  }

  if (nfo->bisulfiterun == 1){
    len = snprintf(NULL, 0, "\tXB:Z:F%u/CT", nfo->bisulfiteprotocol);
    tag = ALLOCMEMORY(space, tag, char, ptr+len+1);
    snprintf(&tag[ptr], len+1, "\tXB:Z:F%u/CT", nfo->bisulfiteprotocol);
    ptr += len;
  } else if (nfo->bisulfiterun == 2){
    len = snprintf(NULL, 0, "\tXB:Z:F%u/GA", nfo->bisulfiteprotocol);
    tag = ALLOCMEMORY(space, tag, char, ptr+len+1);
    snprintf(&tag[ptr], len+1, "\tXB:Z:F%u/GA", nfo->bisulfiteprotocol);	  
    ptr += len;
  } 

  if (nfo->bisulfite){
    /* 
     * get bisulfite mismatches (hence ones that could 
     * be explained by the bisulfite treatment
     */
    mis = getBisulfiteMismatches(al, nfo->bisulfite);
    len = snprintf(NULL, 0, "\tXD:i:%u", mis);
    tag = ALLOCMEMORY(space, tag, char, ptr+len+1);
    snprintf(&tag[ptr], len+1, "\tXD:i:%u", mis);
    ptr += len;

    /*
     * get wrong-strand bisulfite mismatches (hence ones
     * that could be explained by the bisulfite treatment
     * but ONLY on the other genomic strand (i.e. G/A
     * mismatches in C/T matching run or vice versa)
     * => used to identify wrong strand matches
     */
    mis = getWrongStrandBisulfiteMismatches(al, nfo->bisulfite);
    len = snprintf(NULL, 0, "\tXF:i:%u", mis);
    tag = ALLOCMEMORY(space, tag, char, ptr+len+1);
    snprintf(&tag[ptr], len+1, "\tXF:i:%u", mis);
    ptr += len;
  }

  if(previdx != -1 || nextidx != -1) {
    len = snprintf(NULL, 0, "\tXX:i:%d\tXY:i:%d\tXQ:i:%d", ustart, uend, uno);
    tag = ALLOCMEMORY(space, tag, char, ptr+len+1);
    snprintf(&tag[ptr], len+1, "\tXX:i:%d\tXY:i:%d\tXQ:i:%d", ustart, uend, uno);	  
    ptr += len;
  }

  if(previdx != -1) {
    if(map && map->mseq) { 
    if(previdx > 0)
      prevseqstart = map->mseq->markpos[previdx-1]+1; 
      refdesc = ((CharSequence*)map->mseq->ref[previdx].ref)->description;
    } else { 
      prevseqstart = match->prevseqstart;
      refdesc = match->prevseqrefdesc;
    }
    len = snprintf(NULL, 0, "\tXP:Z:%s\tXU:i:%d\tXS:i:%d", refdesc, prevpos-prevseqstart, prevflags);
    tag = ALLOCMEMORY(space, tag, char, ptr+len+1);
    snprintf(&tag[ptr], len+1,"\tXP:Z:%s\tXU:i:%d\tXS:i:%d", refdesc, prevpos-prevseqstart, prevflags);
    ptr += len;
  }

  if(nextidx != -1) {
    if(map && map->mseq) {
    if(nextidx > 0)
      nextseqstart = map->mseq->markpos[nextidx-1]+1; 
      refdesc = ((CharSequence*) map->mseq->ref[nextidx].ref)->description;
    } else {
      nextseqstart = match->nextseqstart;
      refdesc = match->nextseqrefdesc;
    }
    len = snprintf(NULL, 0, "\tXC:Z:%s\tXV:i:%d\tXT:i:%d", refdesc, nextpos-nextseqstart, nextflags);
    tag = ALLOCMEMORY(space, tag, char, ptr+len+1);
    snprintf(&tag[ptr], len+1,"\tXC:Z:%s\tXV:i:%d\tXT:i:%d", refdesc, nextpos-nextseqstart, nextflags);
    ptr += len;
  }

  tag[ptr] = 0;
  FREEMEMORY(space, md);
  return tag;
}


/*--------------------------- repoertMatchSetFlags ---------------------------
 *    
 * @brief set the flags
 * @author Steve Hoffmann 
 *   
 */ 

void 
reportMatchSetFlags(matchstatus_t pStat, char isMate, Uint noofmatepairs, 
    Uint *flag, Uint *mateflag, char *pStatChr, Uint *noofmatches) {

  
  /*is paired in sequencing*/
  *flag |= 1 ;
  *mateflag |= 1;
  *pStatChr= ' ';

  switch(pStat) {             
    /*query itself is unmapped*/
    case MATE: *flag |= (1 << 3);
               /*second in pair*/
               *flag |= (1 << 7);
               *pStatChr = 'M';
               break;
               /*mate is unmapped*/
    case QUERY: *flag |= (1 << 3);
                /*first in pair*/
                *flag |= (1 << 6);
                *pStatChr = 'Q';
                break;
                /*both mapped but separately*/
    case PAIR_INS: 
                *flag |= (1 << 6);
                *mateflag |= (1 << 7);
                *pStatChr = 'p';
                *flag |= (1 << 1);
                *mateflag |= (1 << 1);
                *noofmatches = noofmatepairs;
                break;              
                /*both mapped in reverse order*/
    case PAIR_REV: 
                *flag |= (1 << 7);
                *mateflag |= (1 << 6); 
                *pStatChr = 'R';
                *flag |= (1 << 1);
                *mateflag |= (1 << 1);
                *noofmatches = noofmatepairs;
                break;
                /*both mapped properly*/        
    case PAIR: 
                *flag |= (1 << 6);
                *mateflag |= (1 << 7);
                *pStatChr = 'P';
                *flag |= (1 << 1);
                *mateflag |= (1 << 1);
                *noofmatches = noofmatepairs;
                break;
                /*query spliced, mate full*/
    case QUERY_SPL_FULL_MATE: 
                if (isMate) {
                  /*second in pair*/
                  *flag |= (1 << 7);
                } else {
                  /*first in pair*/
                  *flag |= (1 << 6);
                  /*split hit*/
                  *flag |= (1 << 8);
                  *noofmatches = 1;
                }
                *pStatChr = 'S';
                break;
                /*query spliced, no mate*/
    case QUERY_SPL_NO_MATE: 
                *flag |= (1 << 3);
                /*first in pair*/
                *flag |= (1 << 6);
                /*split hit*/
                *flag |= (1 << 8);
                *noofmatches = 1;
                *pStatChr = 'T';
                break;
                /*mate spliced, query full*/
    case MATE_SPL_FULL_QUERY: 
                if (isMate) {
                  *flag |= (1 << 7);
                  /*split hit*/
                  *flag |= (1 << 8);
                  *noofmatches = 1;
                } else {
                  *flag |= (1 << 6);
                }
                *pStatChr = 'U';
                break;
                /*mate spliced, no query*/
    case MATE_SPL_NO_QUERY: 
                *flag |= (1 << 3);
                /*second in pair*/
                *flag |= (1 << 7);
                /*split hit*/
                *flag |= (1 << 8);
                *pStatChr = 'V';
                *noofmatches = 1;
                break;
                /*both spliced*/
    case PAIR_SPL: 
                if(isMate) {
                  /*second in pair*/
                  *flag |= (1 << 7);
                } else {
                  /*first in pair*/
                  *flag |= (1 << 6);
                }
                /*split hit*/
                *flag |= (1 << 8);
                *noofmatches = 1;
                *pStatChr = 'X';
                break;
  }
}


/*------------------------------- reportMatch --------------------------------
 *    
 * @brief reports a match to device with different output formats
 * @author Steve Hoffmann 
 *   
 */

inline Uint
reportMatch (void *space, Gmap *map, fasta_t *queries, 
    segemehl_t *nfo, matchstatus_t pStat, 
    unsigned char isMate){

  int i, k, j, l, u, seqlen, qrylen, matedesclen, matelen, matchid,
  mateseqlen, desclen;
  Uint off, clipoff=0, mateclipoff=0, seqstart=0, mateseqstart=0, flag=0, 
    mateflag=0, noofmatepairs=0, noofmatches=0, 
       lclip=0, rclip=0, lclipf=0, rclipf=0, matelclip=0, materclip=0, refdesclen=0, materefdesclen=0, 
       fraglen=0, quallen = 0;

  gread_t *reads, *read;
  gmatch_t *match;
  gmate_t *mates;
  Alignment *al;

  FILE *dev=NULL, *matedev=NULL;
  char *qry,
       *qual,
       *description,
       *refdesc,
       *refseq,
       *materefdesc,
       *matedesc,
       *materefseq,
       *mate,
       *matequal,
       *meopstr,
       *matemeopstr,
       *cigar,
       *matecigar, 
       *tmp,
       *tag,
       strands[2] = {'+','-'},
       strandchr,
       matestrandchr,
       pStatChr = 'Q',
       alignlf = '\n',
       lf = '\n',
       **list;
  unsigned char mateReported=0;
  Uint matereport=0; Uint report=0;
  bl_fileBin_t *fx = NULL, *mfx = NULL;


  if (!nfo->bins && (nfo->threadno > 1))  {
    pthread_mutex_lock(nfo->mtx);
  } 

  if((nfo->order || nfo->bisulfitemerging) && nfo->align) {
    lf = 7;
  }

  off = map->mapoffset;
  reads = map->reads;
  dev = nfo->dev;

  for(i = 0; i < map->noofreads; i++) {
    read = &reads[i];

    noofmatches = read->noofmatches;
    noofmatepairs = read->noofmatepairs;
    matchid = 0;
    //fprintf(stderr, "\nread noofmatches:%d, noofmatepairs:%d\n", noofmatches, noofmatepairs);

    for (k = 0; k <= 1; k++) {

      //fprintf(stdout, "k:%d, read->n[k]=%d\n", k, read->n[k]);
      for(j = 0; j < read->n[k]; j++) {

        mateReported = 0;
        tmp = NULL;
        flag = 0;
        mateflag = 0;
        seqstart = 0;

        match = &read->matches[k][j];

        if(match->subject > 0 && map->mseq) {
          seqstart = map->mseq->markpos[match->subject-1]+1; 
        }
        
        if (!isMate){
          bl_fastaGetClipPos(queries, read->id, &lclip, &rclip);
          // if soft-clipping
          if (!nfo->hardclip) bl_fastaSetClip(queries, read->id, 0, 0);
          
          qry = bl_fastaGetSequence(queries, read->id);
          qual = bl_fastaGetQuality(queries, read->id);
          qrylen = bl_fastaGetSequenceLength(queries, read->id);
        }
        else {
          bl_fastaGetMateClipPos(queries, read->id, &lclip, &rclip);
          // if soft-clipping
          if (!nfo->hardclip) bl_fastaSetMateClip(queries, read->id, 0, 0);
          
          qry = bl_fastaGetMate(queries, read->id);
          qual = bl_fastaGetMateQuality(queries, read->id);
          qrylen = bl_fastaGetMateLength(queries, read->id);         
        }
      
        rclipf = rclip;
        lclipf = lclip;

        //fragment clipping
        if(match->nextpos != -1 || match->prevpos != -1) { 
          rclipf = 0;
          lclipf = 0;
        }


        seqlen = off+match->q - match->p;
        al = match->al;

        list = (char **) calloc(OUTLENGTH+1, sizeof(char *));

        if(!isMate) {
          desclen = bl_fastaGetDescriptionLength(queries, read->id)+1;
          description = bl_fastaGetDescription(queries, read->id);
        } else {
          desclen = bl_fastaGetMateDescriptionLength(queries, read->id)+1;
          description = bl_fastaGetMateDescription(queries, read->id);
        }

        sprintstr(&list[QRY_DESC], description, desclen);
        //fprintf(stdout, "CHECKPOINT 1\n");
        //use alignment instead
        if (seqstart > off+match->p) continue;
        strandchr = strands[k];

        /*default quality string is * */
        if (qual == NULL){
          list[QUAL] = calloc(2, 1);
          list[QUAL][0] = '*';
        } else {
          if(match->previdx != -1 || match->nextidx != -1 ) { // }|| qrylen != match->j - match->i + 1) {
            quallen = match->j - match->i + 1; 
            if(lclip+ match->i + quallen >  qrylen) {  
              DBG("warning: wrong fragment alignment for '%s'! skipping alignment.", description);
              continue;
            }
            
            sprintstr(&list[QUAL], &qual[lclip+match->i], quallen);
          } else { 
            quallen = qrylen; 
            sprintstr(&list[QUAL], qual , quallen);
          }
        } 

        if (strandchr == '-') {
          meopstr = multieopstring(al, rclipf, lclipf, 0);
          cigar = cigarstring(al, rclipf, lclipf, (nfo->hardclip) ? 'H':'S', 0);

          if(match->previdx != -1 || match->nextidx != -1) {
            
            fraglen =match->j - match->i + 1; 
            if(lclip+ match->i + fraglen >  qrylen) { 
              DBG("warning: wrong fragment alignment for '%s'! skipping alignment.", description);
              continue;
            }
            
            tmp = charIUPACcomplement(space, &qry[lclip+match->i], fraglen); 
          } else {
            fraglen = qrylen;
            tmp = charIUPACcomplement(space, qry, fraglen);
          }
                      
          sprintstr(&list[QRY_SEQ], tmp, fraglen);
          free(tmp);

          if (qual != NULL){
            list[QUAL] = strrev(list[QUAL], quallen);
          }

          flag |= (1 << 4);
          mateflag |= (1 << 5);


          clipoff = (nfo->hardclip) ? rclipf : 0;

        } else { 
          
          meopstr = multieopstring(al, lclipf, rclipf, 0);
          cigar = cigarstring(al, lclipf, rclipf, (nfo->hardclip) ? 'H':'S', 0);
          
          if(match->previdx != -1 || match->nextidx != -1 ){  //|| qrylen != match->j - match->i +1) {
            
            fraglen = match->j - match->i + 1;
            if(lclip+ match->i + fraglen >  qrylen) { 
              DBG("warning: wrong fragment alignment for '%s'! skipping alignment.", description);
              continue;
            }

            tmp = ALLOCMEMORY(space, NULL, char, fraglen + 1);
            memmove(tmp, &qry[lclip+match->i], fraglen); 
            tmp[fraglen] = 0;
            sprintstr(&list[QRY_SEQ], tmp, fraglen);
            FREEMEMORY(space, tmp);
          } else {
            sprintstr(&list[QRY_SEQ], qry, qrylen);
          }

          clipoff = (nfo->hardclip) ? lclipf : 0;
        }
        
        /* restore clipping positions for next query match */
        if (!isMate) {
          if (!nfo->hardclip) bl_fastaSetClip(queries, read->id, lclip, rclip);
        }
        else {
          if (!nfo->hardclip) bl_fastaSetMateClip(queries, read->id, lclip, rclip);
        }

        if(map->mseq) { 
          refdesc = ((CharSequence*)
              map->mseq->ref[match->subject].ref)->description;
          refdesclen = ((CharSequence*)
              map->mseq->ref[match->subject].ref)->descrlen;
          refseq = map->mseq->sequences+match->p;
        } else {
          refdesc = match->refdesc;
          refdesclen = match->refdesclen;
          refseq = match->refseq;
        }

        list[MEOP_STR] = meopstr;
        list[SAM_CIGAR] = cigar;
        sprintUint(&list[QRY_LEN], qrylen-1);
        sprintUint(&list[SCR], match->scr);
        sprintflt (&list[EVALUE], match->evalue);
        sprintUint(&list[EDIST], match->edist);
        sprintUint(&list[QRY_S], off+match->i);
        sprintUint(&list[QRY_E], off+match->j);
        sprintUint(&list[SEQ_S], off+match->p-seqstart+clipoff);
        sprintUint(&list[SEQ_E], off+match->q-seqstart+clipoff);
        sprintchar(&list[STRAND], strandchr);
        sprintUint(&list[MAT], match->mat);
        sprintUint(&list[MIS], match->mis);
        sprintUint(&list[INS], match->ins);
        sprintUint(&list[DEL], match->del);      
        sprintstr(&list[REF_SEQ], refseq, seqlen);
        sprintUint(&list[NOOFMATCHES], noofmatches);
        sprintstr (&list[SEQ_DESC], refdesc, refdesclen);     
        sprintstr(&list[SAM_QRY], list[QRY_DESC], desclen);
        //strtok(list[SAM_QRY], "/");
        sprintUint(&list[SAM_MAPQ], 255);
        sprintchar(&list[SAM_MATE_REF], '*');
        sprintUint(&list[MATE_SEQ_S], 0);
        sprintUint(&list[SAM_ISIZE], 0);
        
        //like here ...
        if(bl_fastaHasMate(queries)) { 
          reportMatchSetFlags(pStat, isMate, noofmatepairs, 
              &flag, &mateflag, &pStatChr, &noofmatches); 
        } 

        sprintchar(&list[PAIR_STATUS], pStatChr);

       /* if(noofmatepairs) { 
        //...nofmatches -> noofmatepairs here
        tag = getSAMTags(map, off+match->i, off+match->j, match->fragno, match->al, 
            match->edist, 
            match->previdx, match->prevpos+off, match->prevflags, 
            match->nextidx, match->nextpos+off, match->nextflags,
            noofmatepairs, pStatChr, meopstr, nfo);
        } else {*/
         tag = getSAMTags(map, off+match->i, off+match->j, match->fragno, match->al, 
            match->edist, 
            match->previdx, match->prevpos+off, match->prevflags, 
            match->nextidx, match->nextpos+off, match->nextflags, matchid, 
            noofmatches, read->noofmatches, pStatChr, meopstr, nfo, match);
        //}

        sprintstr(&list[TAG], tag, strlen(tag));
        FREEMEMORY(space, tag); 


        if(bl_fastaHasMate(queries) && 
            (pStat == PAIR  || pStat == PAIR_REV || pStat == PAIR_INS)) {

          if (!isMate){
            bl_fastaGetMateClipPos(queries, read->id, &matelclip, &materclip);

            //if soft-clipping
            if (!nfo->hardclip){
              bl_fastaSetMateClip(queries, read->id, 0, 0);
            }
          }
          else {
            bl_fastaGetClipPos(queries, read->id, &matelclip, &materclip);

            //if soft-clipping
            if (!nfo->hardclip){
              bl_fastaSetClip(queries, read->id, 0, 0);
            }
          }

          mates =  match->mates;
          if (!isMate){
            matelen = bl_fastaGetMateLength(queries, read->id);
          }
          else {
            matelen = bl_fastaGetSequenceLength(queries, read->id);
          }


          if(se_kdMatchHasMates(&read->matches[k][j])) { 

            for (u=0; u < 4 && !mateReported; u++) {
              
              
              if(!mates[u].isset || mates[u].edist >
                  match->mateminedist) continue;

              if(map->mseq) {
                if (mates[u].subject > 0){
                  mateseqstart = map->mseq->markpos[mates[u].subject-1]+1; 
                }
                else {
                  mateseqstart = 0;
                }
                materefseq = map->mseq->sequences + mates[u].p;
              } else { 
                mateseqstart = 0;
                materefseq = mates[u].materefseq;
              }

              if (!isMate){
                matedesclen = 
                  bl_fastaGetMateDescriptionLength(queries, read->id)+1;
                matedesc = bl_fastaGetMateDescription(queries, read->id);
                mate = bl_fastaGetMate(queries, read->id);
                matequal = bl_fastaGetMateQuality(queries, read->id);  
              }
              else {
                matedesclen = 
                  bl_fastaGetDescriptionLength(queries, read->id)+1;
                matedesc = bl_fastaGetDescription(queries, read->id);
                mate = bl_fastaGetSequence(queries, read->id);
                matequal = bl_fastaGetQuality(queries, read->id);                
              }
              mateseqlen = off+mates[u].q - mates[u].p;

              sprintstr(&list[MATE_REF_SEQ], materefseq, mateseqlen);
              //default quality string is *
              if (matequal == NULL){
                list[MATE_QUAL] = calloc(2, 1);
                list[MATE_QUAL][0] = '*';
              } else {
                sprintstr(&list[MATE_QUAL], matequal, matelen);
              }

              /*use alignment instead*/
              if (mateseqstart > off+mates[u].p) continue;              
               
              matestrandchr = strands[k];
              if (u & 1) matestrandchr = strands[(~k)&1];

              if (matestrandchr == '-') {
                matemeopstr = multieopstring(mates[u].al, materclip, matelclip, 0);  
                matecigar = cigarstring(mates[u].al, materclip, matelclip, (nfo->hardclip)?'H':'S', 0);
                tmp =  charIUPACcomplement(space, mate, matelen);
                sprintstr(&list[MATE_QRY_SEQ], tmp, matelen);
                free(tmp);
                if (matequal != NULL){
                  list[MATE_QUAL] = strrev(list[MATE_QUAL], matelen);
                }
                flag  |= (1 << 5);
                mateflag  |= (1 << 4);
                mateclipoff = (nfo->hardclip) ? rclip : 0;
              } else { 
                matemeopstr = multieopstring(mates[u].al, matelclip, materclip, 0);     
                matecigar = cigarstring(mates[u].al, matelclip, materclip, (nfo->hardclip)?'H':'S', 0);
                sprintstr(&list[MATE_QRY_SEQ], mate, matelen);
                mateclipoff = (nfo->hardclip) ? lclip : 0;
              }

              list[MATE_MEOP] = matemeopstr;
              list[SAM_MATE_CIGAR] = matecigar;

              sprintUint(&list[MATE_LEN], matelen);
              sprintUint(&list[MATE_SCR], mates[u].mat - mates[u].edist);
              sprintflt (&list[MATE_EVALUE], 0.0);
              sprintUint(&list[MATE_QRY_S], 1);
              sprintUint(&list[MATE_QRY_E], matelen);

              free(list[MATE_SEQ_S]);

              if(map->mseq) { 
                materefdesc = ((CharSequence*)
                    map->mseq->ref[mates[u].subject].ref)->description;
                materefdesclen = ((CharSequence*)
                    map->mseq->ref[mates[u].subject].ref)->descrlen;
              } else {
                materefdesc = mates[u].materefdesc;
                materefdesclen = mates[u].materefdesclen;
              }

              sprintUint(&list[MATE_SEQ_S], off+mates[u].p - mateseqstart + mateclipoff);
              sprintUint(&list[MATE_SEQ_E], off+mates[u].q - mateseqstart + mateclipoff);
              sprintUint(&list[MATE_MAT], mates[u].mat); 
              sprintUint(&list[MATE_MIS], mates[u].mis);
              sprintUint(&list[MATE_INS], mates[u].ins);
              sprintUint(&list[MATE_DEL], mates[u].del);
              sprintUint(&list[MATE_EDIST], mates[u].edist); 
              sprintchar(&list[MATE_STRAND], matestrandchr);
              sprintstr (&list[MATE_SEQ_DESC], materefdesc, materefdesclen); 
              sprintUint(&list[MATE_NOOFMATCHES], noofmatepairs);
              sprintstr (&list[MATE_DESC], matedesc, matedesclen-1);
              
              sprintstr(&list[SAM_MATE_QRY], list[MATE_DESC], desclen);
              //strtok(list[SAM_MATE_QRY], "/");
              
              sprintUint(&list[SAM_FLAG], flag);
              sprintUint(&list[SAM_MATE_FLAG], mateflag);
              sprintUint(&list[SAM_MATE_MAPQ], 255);

              if(strcmp(list[MATE_SEQ_DESC],list[SEQ_DESC])) {
                free(list[SAM_MATE_REF]);   

                sprintstr(&list[SAM_MATE_REF], list[MATE_SEQ_DESC], 
                    strlen(list[MATE_SEQ_DESC])); 

                sprintstr(&list[SAM_QRY_REF], list[SEQ_DESC], 
                    strlen(list[SEQ_DESC]));
                
                sprintint(&list[SAM_MATE_ISIZE], 0);
              
              } else {
                
                sprintf(list[SAM_MATE_REF],"=");
                sprintchar(&list[SAM_QRY_REF], '=');
                free(list[SAM_ISIZE]);
                
                sprintint(&list[SAM_ISIZE], 
                    (off+mates[u].p - mateseqstart + matelen)-
                    (off+match->p-seqstart));
                
                sprintint(&list[SAM_MATE_ISIZE],
                    (off+match->p-seqstart)-
                    (off+mates[u].p - mateseqstart + matelen));
              }

              //noofmatematches -> noofmatepairs 
              //the number of splits should be 0 because
              //split alignments are always single
              tag = getSAMTags(map, 1, matelen, 1, 
                  mates[u].al, mates[u].edist, 
                  -1, -1, 0, -1, -1, 0, matchid, noofmatepairs, 0, pStatChr, matemeopstr, nfo, NULL);

              sprintstr(&list[MATE_TAG], tag, strlen(tag));
              FREEMEMORY(space, tag);

       
              if(mfx) {
                bl_fileBinsUnlock(mfx);
                mfx = NULL;
              }

              if(nfo->rep_type == 5) {

                dev = getDevice(space, refdesc, off+match->p-seqstart, &fx, nfo);
                genericOutput(dev, list, 5, lf);
                report++;

                if(nfo->align) {
                  if(nfo->order || nfo->bisulfitemerging) {
                    alignlf = 7;
                  }
                  showAlignLF(match->al, dev, alignlf);
                  fprintf(dev, "\n");
                } 

                if(nfo->bins && !nfo->bisulfitemerging) {
                  /*TODO: avoid lock!*/
 
                  mfx = bl_fileBinsDomainGetBin(nfo->bins, materefdesc, 
                      off+mates[u].p - mateseqstart);

                  if(mfx != fx) { 
                    bl_fileBinsUnlock(fx); 
                    fx = NULL;
                    matedev = getDevice(space, materefdesc, 
                        off+mates[u].p - mateseqstart, &mfx, nfo);
                  } else {
                    mfx = NULL;
                    matedev = dev;
                  }
                } else { 
                  matedev = dev; 
                }

                genericOutput(matedev, list, 6, lf);
                matereport++;

                if(nfo->align) {
                  if(nfo->order || nfo->bisulfitemerging) {
                    alignlf = 7;
                  }
                  showAlignLF(mates[u].al, dev, alignlf);
                  fprintf(dev, "\n");
                } 
              }

              if(nfo->rep_type == 12) {

                dev = getDevice(space, refdesc, off+match->p-seqstart, &fx, nfo);
                genericOutput(dev, list, 11, lf);
                report++;

                if(nfo->align) {
                  if(nfo->order || nfo->bisulfitemerging) {
                    alignlf = 7;
                  }
                  showAlignLF(match->al, dev, alignlf);
                  fprintf(dev, "\n");
                } 

                if(nfo->bins && !nfo->bisulfitemerging) {
                  /*TODO: avoid lock!*/
 
                  mfx = bl_fileBinsDomainGetBin(nfo->bins, materefdesc, 
                      off+mates[u].p - mateseqstart);

                  if(mfx != fx) {                   
                    bl_fileBinsUnlock(fx);                    
                    fx = NULL;
                    matedev = getDevice(space, materefdesc, 
                        off+mates[u].p - mateseqstart, &mfx, nfo);
                  } else {
                    mfx = NULL;
                    matedev = dev;
                  }
                } else { 
                  matedev = dev; 
                }

                genericOutput(matedev, list, 13, lf);
                matereport++;

                if(nfo->align) {
                  if(nfo->order || nfo->bisulfitemerging) {
                    alignlf = 7;
                  }
                  showAlignLF(mates[u].al, dev, alignlf);
                  fprintf(dev, "\n");
                } 
              }

              if(nfo->rep_type == 15) { 
                
                dev = getDevice(space, refdesc, off+match->p-seqstart, &fx, nfo);
                genericOutput(dev, list, 15, lf);
                report++;

                if(nfo->align) {
                  if(nfo->order || nfo->bisulfitemerging) {
                    alignlf = 7;
                  }
                  showAlignLF(match->al, dev, alignlf);
                  fprintf(dev, "\n");
                } 
  
                if(nfo->bins && !nfo->bisulfitemerging) {
                  /*TODO: avoid lock!*/
 
                  mfx = bl_fileBinsDomainGetBin(nfo->bins, materefdesc, 
                      off+mates[u].p - mateseqstart);

                  if(mfx != fx) {                    
                    bl_fileBinsUnlock(fx);  
                    fx = NULL;
                    matedev = getDevice(space, materefdesc, 
                        off+mates[u].p - mateseqstart, &mfx, nfo);
                  } else {
                    mfx = NULL;
                    matedev = dev;
                  }
                } else { 
                  matedev = dev; 
                }

                genericOutput(matedev, list, 16, lf); 
                matereport++;

                if(nfo->align) {
                  if(nfo->order || nfo->bisulfitemerging) {
                    alignlf = 7;
                  }
                  showAlignLF(mates[u].al, dev, alignlf);
                  fprintf(dev, "\n");
                }       
              }
              if (match->previdx == -1 && match->nextidx == -1){
                matchid++;
              }

              mateReported = 1;

              for (l=MATE_LEN; l <= MATE_QUAL; l++) {
                if(list[l]) free(list[l]);
                list[l] = NULL;
              }

              for(l=SAM_MATE_FLAG; l <= SAM_MATE_ISIZE; l++) { 
                if(list[l]) free(list[l]);
                list[l] = NULL;
              }
            }
            assert(mateReported);
          } 
          /* restore clipping positions for next query match */
          if (!isMate){
            if (!nfo->hardclip){
              bl_fastaSetMateClip(queries, read->id, matelclip, materclip);
            }
          }
          else {
            if (!nfo->hardclip){
              bl_fastaSetClip(queries, read->id, matelclip, materclip);
            }

          }
        } else {    
          
          
          dev = getDevice(space, refdesc, off+match->p-seqstart, &fx, nfo);
          sprintUint(&list[SAM_FLAG], flag);
          genericOutput(dev, list, nfo->rep_type, lf);  
          report++; 
          if (match->previdx == -1 && match->nextidx == -1){
            matchid++;
          }

          if(match->al) {
            if (nfo->align){
              if (nfo->order || nfo->bisulfitemerging){
                alignlf = 7;
              }
              showAlignLF(match->al, dev, alignlf);
              fprintf(dev, "\n");
            }
          }
        }

        for (l = 0; l < OUTLENGTH; l++){
          if(list[l]) free(list[l]);
        }
        free(list);

        if(nfo->bins) {
          if(fx) { 
            bl_fileBinsUnlock(fx);
            fx = NULL;
          }
          if(mfx) {
            bl_fileBinsUnlock(mfx);
            mfx = NULL;
          }
        }
      }
    }    
  }

  if (!nfo->bins && nfo->threadno > 1){
    pthread_mutex_unlock(nfo->mtx);
  } 

  if(report && matereport) return 3;
  if(!report && matereport) return 2;
  if(report && !matereport) return 1;

  return 0;
}



/*------------------------------- se_SAMHeader -------------------------------
 *    
 * @brief SAM header
 * @author Steve Hoffmann 
 *   
 */

  char*
se_SAMHeader (void *space, char **seq, Uint *seqlen, 
    Uint size, char *cmdline, char sep, char lf, 
    unsigned char sorted)
{

  Uint i,len=1000, curlen=0;
  char *header;

  len += strlen(VERSION);
  if(cmdline)
    len += strlen(cmdline); 


  for(i=0; i < size; i++) {
    len += snprintf(NULL, 0, "@SQ%cSN:%s%cLN:%d%c", sep, seq[i], sep, seqlen[i], lf);
  }

  header = calloc(len, sizeof(char));
  sprintf(header,"@HD%cVN:1.0",sep);
  curlen = strlen(header);

  if(sorted) sprintf(&header[curlen], "%cSO:coordinate", sep);

  curlen = strlen(header);
  sprintf(&header[curlen],"%c",lf);

  for(i=0; i < size; i++) {

    curlen = strlen(header);
    sprintf(&header[curlen],"@SQ%cSN:%s%cLN:%d%c", sep, seq[i], sep, seqlen[i], lf);
  }

  curlen = strlen(header);
  sprintf(&header[curlen],"@PG%cID:segemehl", sep);

  curlen = strlen(header);
  sprintf(&header[curlen],"%cVN:%s", sep, VERSION);

  curlen = strlen(header);
  if(cmdline)
    sprintf(&header[curlen],"%cCL:%s", sep, cmdline);

  curlen = strlen(header);
  sprintf(&header[curlen],"%c",lf);

  return header;
}

/*----------------------------- se_initChromBins -----------------------------
 *    
 * @brief set up bins for chromosomes
 * @author Steve Hoffmann 
 *   
 */

  bl_fileBinDomains_t*
se_createChromDomains (void *space, fasta_t *f, Uint avgbins, Uint maxbins, 
    char *filetemplate, Uint tmplen)
{
  bl_fileBinDomains_t* domains;
  char **desc;
  Uint *size;
  Uint i, no, total=0;

  no = f->noofseqs;
  if(no > maxbins) return NULL;

  desc = ALLOCMEMORY(space, NULL, char*, no);
  size = ALLOCMEMORY(space, NULL, Uint, no);

  for(i=0; i < no; i++) {
    desc[i] = bl_fastaGetDescription(f, i);
    size[i] = bl_fastaGetSequenceLength(f, i);
    total += size[i];
  } 

  domains = bl_fileBinsDomainsInit(space, desc, size, no, total,  
      avgbins, maxbins, filetemplate, tmplen);

  FREEMEMORY(space, desc);
  FREEMEMORY(space, size);

  return domains;
}

/*----------------------------- se_initChromBins -----------------------------
 *    
 * @brief set up bins for chromosomes
 * @author Steve Hoffmann 
 *   
 */

  bl_fileBins_t*
se_createChromBins (void *space, fasta_t *f, int maxbins, char *template, 
    Uint tmplen)
{
  bl_fileBins_t* bins;
  char **desc;
  Uint i, no;

  no = f->noofseqs;
  if(no > maxbins) return NULL;

  bins = ALLOCMEMORY(space, NULL, bl_fileBins_t, 1);
  desc = ALLOCMEMORY(space, NULL, char*, no);
  bl_fileBinsInit(space, bins);

  for(i=0; i < no; i++) {
    desc[i] = bl_fastaGetDescription(f, i);
  }

  bl_fileBinsAdd(space, bins, no, bl_fileBinCClassAssign, desc, NULL, 
      template, tmplen);

  FREEMEMORY(space, desc);
  return bins;
}

/*-------------------------------- se_createBisulifteBins ---------------------
 *
 * @brief set up bin domains for matching runs and threads,
 *        domain names are simply 0...(noofdomains-1) as strings
 * @author Christian Otto
 *
 */

bl_fileBinDomains_t*
se_createBisulfiteBins (void *space, Uint noofdomains,
    Uint threadno, char *filetemplate, Uint tmplen){
  Uint i, j;
  bl_fileBinDomains_t *domains;
  bl_fileBinClass_t *class;

  domains = ALLOCMEMORY(space, NULL, bl_fileBinDomains_t, 1);
  domains->noofdomains = noofdomains;
  domains->exp = 0;
  domains->domain = ALLOCMEMORY(space, NULL, bl_fileBinDomain_t, noofdomains);

  for (i = 0; i < noofdomains; i++){    
    domains->domain[i].domainsize = threadno;    
    domains->domain[i].domainname = ALLOCMEMORY(space, NULL, char, log10(i+1) + 3);
    snprintf(domains->domain[i].domainname, log10(i+1)+2, "%d", i);

    bl_fileBinsInit(space, &domains->domain[i].bins);
    bl_fileBinsAdd(space, &domains->domain[i].bins, threadno, NULL, NULL, NULL,
        filetemplate, tmplen);

    domains->domain[i].bins.noofbins = threadno;

    for (j = 0; j < domains->domain[i].bins.noofbins; j++){
      class = ALLOCMEMORY(space, NULL, bl_fileBinClass_t, 1);
      class->start = j;
      class->end = j;
      class->classname = NULL;
      domains->domain[i].bins.b[j].id = class;
    }
  }
  /*
     DBG("domains: noofdomains=%u\texp=%u\n", domains->noofdomains, domains->exp);
     for (i = 0; i < domains->noofdomains; i++){
     DBG("domain %u: domainname=%s\tdomainsize=%u\tnoofbins=%u\n", i,
     domains->domain[i].domainname, domains->domain[i].domainsize,
     domains->domain[i].bins.noofbins);
     for (j = 0; j < domains->domain[i].bins.noofbins; j++){
     DBG("bin %u: filename=%s\tstart=%llu\tend=%llu\n", j, domains->domain[i].bins.b[j].fname,
     domains->domain[i].bins.b[j].id->start, domains->domain[i].bins.b[j].id->end);
     }
     }*/
  return domains;

}

/*----------------------------- se_defaultHeader ------------------------------
 *    
 * @brief get default header
 * @author Steve Hoffmann 
 *   
 */

  char *
se_defaultHeader (void *space, segemehl_t *info, char sep, char lf)
{
  char *buffer, *timestr; 
  struct tm *timeinfo;
  time_t rawtime;

  time(&rawtime);
  timeinfo = localtime (&rawtime);
  timestr = asctime(timeinfo);
  timestr[strlen(timestr)-1] = 0;

  buffer = ALLOCMEMORY(space, NULL, char, 5000);
  memset(buffer, 0, 5000);
  snprintf(buffer, 5000, "#segemehl %s%c#query: %s (%u seqs)%c#mate: %s%c#subject: %s\
      %c#minsize=%d, diff_seed=%d, jump=%d, acc=%d, maxEvalue: %.5f, hitstrategy: %d\
      %c#splitreads: %s%c%s%c",  
      timestr, lf, info->queryfilename, info->reads->noofseqs, lf, info->matefilename, lf,
      info->idxfilename, lf, info->minsize, info->k_p, info->jump,
      info->accuracy, info->maxevalue, info->bestonly, lf, info->splitfilebasename, lf,
      HEAD[info->rep_type], lf);
  return buffer;
}

/*--------------------------- registerOutputDevice ---------------------------
 *    
 * @brief select the output device for matches and print header
 * @author Steve Hoffmann 
 *   
 */

void
se_registerOutputDevice(void *space, segemehl_t *info) {
  Uint *dmlen, dmno, i;
  char *buffer, **dms;


  if(info->outfile) {
    info->dev=fopen(info->outfile, "w");
    if (info->dev == NULL) {
      fprintf(stderr, "Couldn't open file '%s'. Exit forced.\n", info->outfile);
      exit(EXIT_FAILURE);
    }
  }
  if (info->nohead){
    return;
  }

  if (info->rep_type != 15) {
    if(info->order) {
      buffer = se_defaultHeader(space, info, 8 ,7);
      buffer[strlen(buffer)-1] = 29;
    } else {
      buffer = se_defaultHeader(space, info, '\t', '\n');
    }
  } else {

    dmno = info->fasta->noofseqs;
    dms = ALLOCMEMORY(space, NULL, char**, dmno);
    dmlen = ALLOCMEMORY(space, NULL, Uint*, dmno);

    for(i=0; i < dmno; i++) {
      dms[i] = bl_fastaGetDescription(info->fasta, i); 
      dmlen[i] = bl_fastaGetSequenceLength(info->fasta, i); 
    }

    if(info->order){ 
      buffer = se_SAMHeader(space, dms, dmlen, dmno, info->cmdline, 8, 7, info->order);
      buffer[strlen(buffer)-1] = 29; 
    } else {
      buffer = se_SAMHeader(space, dms, dmlen, dmno, info->cmdline, '\t', '\n', info->order);
    }


    FREEMEMORY(space, dms);
    FREEMEMORY(space, dmlen);
  }

  fprintf(info->dev, "%s", buffer);
  if (info->order){
    fprintf(info->dev, "\n");
  }
  FREEMEMORY(space, buffer);
}




/*------------------------------- se_storeHeader ------------------------------
 *    
 * @brief read and store header (delimitted by '\n') from file
 * @author Steve Hoffmann 
 *   
 */
void
se_storeHeader(void *space, char *filename, char **header, Uint *headerlen){
  FILE *fp;
  int ret;
    
  fp = fopen(filename, "rb+");
  if(!fp) {
    fprintf(stderr,"Couldnt open file '%s'. Exit forced!\n", filename);
    exit(-1);
  }
  
  ret = bl_fgets(space, fp, header);

  if (ret == EOF || ret < 0){
    fprintf(stderr, "Couldn't retrieve header information. End of file reached.\n");
    exit(-1);
  }
  fclose(fp);
  
  *headerlen = (Uint) ret;

  return;
}
