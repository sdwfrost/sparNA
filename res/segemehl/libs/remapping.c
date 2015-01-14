/**
 * remapping.c
 * remapping of unmapped reads
 *
 * @author Christian Otto
 * @email christian@bioinf.uni-leipzig.de
 * @company Bioinformatics, University of Leipzig
 * @date Fri Oct 12 09:55:49 CEST 2012
 */

/*
 * SVN
 * Revision of last commit: $Rev: 407 $
 * Author: $Author: steve $
 * Date: $Date: 2014-02-06 04:55:25 -0500 (Thu, 06 Feb 2014) $
 * Id: $Id: remapping.c 407 2014-02-06 09:55:25Z steve $
 * Url: $URL: http://www2.bioinf.uni-leipzig.de/svn5/segemehl/libs/remapping.c $
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "debug.h"
#include "info.h"
#include "stringutils.h"
#include "basic-types.h"
#include "mathematics.h"
#include "container.h"
#include "list.h"
#include "sort.h"
#include "manout.h"
#include "manopt.h"
#include "multicharseq.h"
#include "charsequence.h"
#include "iupac.h"
#include "alignment.h"
#include "sw.h"
#include "vtprogressbar.h"
#include "segemehl.h"
#include "realign.h"
#include "remapping.h"

void
rm_initStat (remappingstat_t *stat){
  stat->status = NULL;
  stat->aligns = NULL;
  stat->maxdist = NULL;
  stat->extensions = NULL;
  stat->n = 0;
}

void
rm_addStat (void *space, remappingstat_t *stat, 
	    remappingstatus_t status, Uint aligns,
	    Uint extensions, Uint maxdist){
  stat->status = ALLOCMEMORY(space, stat->status, remappingstatus_t, stat->n+1);
  stat->aligns = ALLOCMEMORY(space, stat->aligns, Uint, stat->n+1);
  stat->extensions = ALLOCMEMORY(space, stat->extensions, Uint, stat->n+1);
  stat->maxdist = ALLOCMEMORY(space, stat->maxdist, Uint, stat->n+1);
  stat->status[stat->n] = status;
  stat->aligns[stat->n] = aligns;
  stat->extensions[stat->n] = extensions;
  stat->maxdist[stat->n] = maxdist;
  stat->n++;
}

void
rm_destructStat (void *space, remappingstat_t *stat){
  if (stat->status)
    FREEMEMORY(space, stat->status);
  if (stat->aligns)
    FREEMEMORY(space, stat->aligns);
  if (stat->maxdist)
    FREEMEMORY(space, stat->maxdist);
  if (stat->extensions)
    FREEMEMORY(space, stat->extensions);
  stat->n = 0;
}

void
rm_updateProgressBar(Uint k, remapping_t *nfo) {

  if (!nfo->mute) {
    if (nfo->counter == NULL) {
      progressBarVT("reads remapped.", nfo->reads->noofseqs, k, 25);
    } else {
      (*nfo->counter)++;
    }
  }
  return;
}

void*
remappingworker(void *args) {
  remapping_t *t;
  
  t = (remapping_t*) args;
  bl_remapping(t->space, t->seq, t->reads, t->L, t->R, t);
  return NULL;
}

int cmp_remappingclust_qsort(const void *a, const void* b) {
  remappingclust_t  *aclust = (remappingclust_t*) a;
  remappingclust_t  *bclust = (remappingclust_t*) b;

  if(aclust->cnt > bclust->cnt) return -1;
  if(aclust->cnt <  bclust->cnt) return 1;
  
  return 0;

}

void bl_remappingPrintAlign(FILE *dev, MultiCharSeqAlignment *mcsa, remapping_t *nfo){
  Uint ustart, uend;

  ustart = mcsa->al->uoff;
  if (mcsa->strand == 1) {
    uend = mcsa->al->ulen - ustart - 1;
    ustart = uend - getUalignlen(mcsa->al) + 1;        
  } else {
    uend = ustart + getUalignlen(mcsa->al) - 1;
  }

  fprintf(dev, "chr=%u, strand=%u, refstart=%u, refend=%u, reflen=%u, vstart=%u, vend=%u, vlen=%u, qrylen=%u, ustart=%u, uend=%u, ulen=%u, score=%d, edist=%d\n", mcsa->subidx, mcsa->strand, mcsa->refstart - mcsa->substart, mcsa->refstart - mcsa->substart + mcsa->reflen - 1, mcsa->reflen, mcsa->refstart - mcsa->substart + mcsa->al->voff, mcsa->refstart - mcsa->substart + mcsa->al->voff + getValignlen(mcsa->al) - 1, getValignlen(mcsa->al), mcsa->qrylen, ustart, uend, getUalignlen(mcsa->al), getAlignScore(mcsa->al, nfo->scores, nfo->indel), getEdist(mcsa->al));
  
  showAlign(mcsa->al, dev);
}

void bl_remappingPrintClust(FILE *dev, remappingclust_t *clust, remapping_t *nfo){
  Uint j;

  fprintf(dev, "dist: type=%u, chr=%u, clust=%d", clust->type, clust->chr, clust->clust);
  if (clust->clust != -1){
    fprintf(dev, ", a=%u, b=%u, median=%u", clust->cluster->a - 1, clust->cluster->b - 1, clust->cluster->median - 1);
  }
  fprintf(dev, ", noofmcsa=%u, score=%d\n", clust->noofmcsa, clust->score);

  if (clust->noofmcsa > 0 && clust->mcsa != NULL){
    for (j = 0; j < clust->noofmcsa; j++){
      fprintf(dev, "frag[%u]: ", j);
      bl_remappingPrintAlign(dev, clust->mcsa[j], nfo);
    }
  }
}

remappingclust_t *bl_remappingAlign(void *space, char **seqs, Container *dist, remapping_t *nfo){
  Uint i, j, noofmcsa, ulen, vlen, qrylen = 0, *reflens, *strands;
  int *M, **lmr, **lmv, **lmc, edist, score, purge;
  int bestscore;
  char **refseqs;
  Alignment **aligns;
  remappingclust_t *cur, *best;

  best = NULL;
  bestscore = -1;  

  assert(bl_containerSize(dist) > 0);
  noofmcsa = ((remappingclust_t *) bl_containerGet(dist, 0))->noofmcsa;

  reflens = ALLOCMEMORY(space, NULL, Uint, noofmcsa);
  strands = ALLOCMEMORY(space, NULL, Uint, noofmcsa);
  refseqs = ALLOCMEMORY(space, NULL, char*, noofmcsa);
  aligns = ALLOCMEMORY(space, NULL, Alignment*, noofmcsa);

  for (i = 0; i < bl_containerSize(dist); i++){
    cur = bl_containerGet(dist, i);    
    assert(cur->noofmcsa == noofmcsa);
    assert(noofmcsa < 3);

    for (j = 0; j < cur->noofmcsa; j++){
      if (cur->mcsa[j]->reflen < nfo->minfragmentalignlen)
        continue;
      if (qrylen == 0) 
        qrylen = cur->mcsa[j]->qrylen;
    }


    for (j = 0; j < noofmcsa; j++){
      aligns[j] = cur->mcsa[j]->al;
      refseqs[j] = cur->mcsa[j]->refseq;
      reflens[j] = cur->mcsa[j]->reflen;
      strands[j] = cur->mcsa[j]->strand;
    }    
    
    // no transition penalty
    M = localmultisplicedmatrix(space, seqs[0], seqs[1], qrylen,
                                refseqs, reflens, strands, noofmcsa,
                                nfo->indel, 0, constscr, nfo->scores,
                                &lmv, &lmr, &lmc);
    assert(M != NULL);

    localmultisplicedtraceback(space, M, seqs[0], seqs[1], qrylen, 
      refseqs, reflens, strands, noofmcsa, nfo->indel, 0,
      constscr, nfo->scores, aligns, lmv, lmr, lmc);
    
    purge = 0;
    cur->score = 0;
    for (j = 0; j < noofmcsa; j++){
      ulen = getUalignlen(aligns[j]);
      vlen = getValignlen(aligns[j]);
      score = getAlignScore(aligns[j], nfo->scores, nfo->indel);
      edist = getEdist(aligns[j]);      
      cur->score += score;

      // check if both alignment fragments
      // meet the given requirements
      if (ulen < nfo->minfragmentalignlen ||
          vlen < nfo->minfragmentalignlen ||
          score < nfo->minfragmentalignscore ||
          edist > (ulen - ceil((nfo->accuracy * ulen)/100))) { 
        purge = 1;
        //fprintf(stderr, "purge: ulen=%u, vlen=%u, score=%d, edist=%u\n", ulen,
        //        vlen, score, edist);
      }
    }

    if (!purge){
      /*
      for (j = 0; j < noofmcsa; j++){
        ustart = aligns[j]->uoff;
        vstart = aligns[j]->voff;
        ulen = getUalignlen(aligns[j]);
        vlen = getValignlen(aligns[j]);
        score = getAlignScore(aligns[j], nfo->scores, nfo->indel);
        totalscore += score;
        
        if (strands[j] == 1){
          uend = qrylen - ustart - 1;
          ustart = uend - ulen + 1;
        }
        else {
          uend = ustart + ulen - 1;
        }
        
        }*/

      if (cur->score > bestscore){
        bestscore = cur->score;
        best = cur;
      }
    }

    for (j = 0; j < noofmcsa; j++){
      //wrapAlignment(aligns[j]);
      //FREEMEMORY(space, aligns[j]);
      FREEMEMORY(space, lmv[j]);
      FREEMEMORY(space, lmr[j]);
      FREEMEMORY(space, lmc[j]);
    }
    FREEMEMORY(space, lmv);
    FREEMEMORY(space, lmr);
    FREEMEMORY(space, lmc);
    FREEMEMORY(space, M);
  }
  
  FREEMEMORY(space, reflens);
  FREEMEMORY(space, refseqs);
  FREEMEMORY(space, strands);
  FREEMEMORY(space, aligns);

  return best;
}

void bl_remappingUpdateAlignSeqs(void *space, MultiCharSeq *seq, char **seqs, Uint qrylen, remappingclust_t *cur){
  Uint i, substart, subend;
  MultiCharSeqAlignment *a;

  for (i = 0; i < cur->noofmcsa; i++){
    a = cur->mcsa[i];
    getMultiCharSeqIdxBounds(seq, a->subidx, &substart, &subend);
    a->substart = substart;
    a->subend = subend;
    
    a->refstart += a->substart;
    a->refseq = &seq->sequences[a->refstart];
    a->reflen = (a->subend > (Lint) a->refstart + a->reflen) ? 
      a->reflen : a->subend - a->refstart + 1;
    
    a->query = seqs[a->strand];
    a->qrylen = qrylen;
    a->al = ALLOCMEMORY(space, NULL, Alignment, 1);
    initAlignment(a->al, a->query, a->qrylen, 0, a->refseq, a->reflen, 0);
  }
}

void bl_remappingUpdateAlignOff(void *space, char **seqs, Uint qrylen, remappingclust_t *cur, Uint offset, Uint right){
  Uint i, ulen;
  int edist;
  MultiCharSeqAlignment *mcsa;
  
  for (i = 0; i < cur->noofmcsa; i++){
    mcsa = cur->mcsa[i];

    /* for safety: ulen & edist should not change */
    ulen = getUalignlen(mcsa->al);
    edist = getEdist(mcsa->al);

    //showAlign(mcsa->al, stderr);
    //fprintf(stderr, "%s\n%s\n", seqs[mcsa->strand], mcsa->al->u);
    
    /* 
     * replace query substring by entire one,
     * adjust offsets and lengths
     */
    mcsa->query = seqs[mcsa->strand];
    mcsa->qrylen = qrylen;
    if (mcsa->al->ulen != qrylen){
      mcsa->al->meops = ALLOCMEMORY(space, mcsa->al->meops, Multieop, mcsa->al->vlen + qrylen);
    }
    mcsa->al->ulen = qrylen;
    mcsa->al->u = seqs[mcsa->strand];
    if (right && mcsa->strand == 0) mcsa->al->uoff += offset;
    if (!right && mcsa->strand == 1) mcsa->al->uoff += offset;

    //fprintf(stderr, "offset=%u, uoff=%u, ulen=%u - %u, edist=%d - %d\n",
    //        offset, mcsa->al->uoff, ulen, getUalignlen(mcsa->al), edist,
    //        getEdist(mcsa->al));
    //showAlign(mcsa->al, stderr);
    
    assert(ulen == getUalignlen(mcsa->al) &&
           edist == getEdist(mcsa->al));
  }
}

void bl_remappingDestructClust(void *elem){
  Uint i;
  remappingclust_t *cur = (remappingclust_t *) elem;

  if (cur->mcsa != NULL && cur->noofmcsa > 0){
    for (i = 0; i < cur->noofmcsa; i++){
      if (cur->mcsa[i] != NULL){
        if (cur->mcsa[i]->al != NULL){
          wrapAlignment(cur->mcsa[i]->al);
          free(cur->mcsa[i]->al);
        }
        free(cur->mcsa[i]);
      }
    }
    free(cur->mcsa);
  }
}

void bl_remappingInitAlign(void *space, MultiCharSeqAlignment *mcsa, Uint subidx,
                         Uint start, Uint len, unsigned char strand){
  mcsa->subidx = subidx;
  mcsa->refstart = start;
  mcsa->reflen = len;
  mcsa->strand = strand;
  mcsa->al = NULL;
}

remappingclust_t *bl_remappingUnsplicedAlign(void *space, MultiCharSeq *seq, char **seqs, Uint len,
                                             Uint chr, Uint start, Uint end, unsigned char strand, remapping_t *nfo){
  Container *a;
  remappingclust_t *clust, *ret;

  a = ALLOCMEMORY(space, NULL, Container, 1);
  bl_containerInit(a, 1, sizeof(remappingclust_t));

  clust = ALLOCMEMORY(space, NULL, remappingclust_t, 1);
  clust->chr = 0;
  clust->type = 0;
  clust->clust = -1;
  clust->cluster = NULL;
  clust->mcsa = ALLOCMEMORY(space, NULL, MultiCharSeqAlignment*, 1);
  clust->mcsa[0] = ALLOCMEMORY(space, NULL, MultiCharSeqAlignment, 1);
  clust->noofmcsa = 1;
  bl_remappingInitAlign(space, clust->mcsa[0], chr, start, end-start+1, strand);
  bl_remappingUpdateAlignSeqs(space, seq, seqs, len, clust);
  bl_containerAdd(a, clust);
  
  ret = bl_remappingAlign(space, seqs, a, nfo);
  if (ret == NULL){
    bl_remappingDestructClust(clust);
    FREEMEMORY(space, clust);
  }
  else {
    memmove(clust, ret, sizeof(remappingclust_t));
  }
  
  bl_containerDestruct(a, NULL);
  FREEMEMORY(space, a);

  return(clust);
}

void bl_remappingGetDist(void *space, Container *a, matchsplitsiteclusterlist_t *L, matchsplitsiteclusterlist_t *R,
                         Uint type, Uint chr, Uint clust, unsigned char strand, Uint last, Uint margin, unsigned char right){
  Uint i, j, curstart, curlen, diststart, distlen;
  matchsplitsiteclusterlist_t *list, *other;
  remappingclust_t dist;
    
  if (type == 0){
    list = L;
    other = R;
  }
  else {
    list = R;
    other = L;
  }
  
  for (i = 0; i < list[chr].cluster[clust].noofdistclust; i++){
    dist.chr = list[chr].cluster[clust].distclustchr[i];
    dist.clust = list[chr].cluster[clust].distclust[i];
    if (list[chr].cluster[clust].distclusttype[i] == 1){
      dist.type = 1-type;
      dist.cluster = &other[dist.chr].cluster[dist.clust];
    }
    else {
      dist.type = type;
      dist.cluster = &list[dist.chr].cluster[dist.clust];
    }
    dist.cnt = list[chr].cluster[clust].distclustcnt[i];
    
    /* - 1 offset on median */
    if (type == 0){
      curlen = (last + 1 + 1 >= list[chr].cluster[clust].median) ? 
        last - list[chr].cluster[clust].median + 1 + 1 : 0;
      curstart = list[chr].cluster[clust].median - 1;
    }
    else {
      curlen = (list[chr].cluster[clust].median - 1 + 1 >= last)? 
        list[chr].cluster[clust].median - 1 - last + 1 : 0;
      curstart = last;
    }
    
    if (margin <= curlen || curlen == 0){
      continue;
    }

    /* - 1 offset on median */    
    if (dist.type == 0){
      diststart = dist.cluster->median - 1;
    }
    else {
      diststart = (dist.cluster->median - 1 + curlen + 1 > margin) ?
        dist.cluster->median - 1 - margin + curlen + 1 : 0;
    }
    distlen = margin - curlen;
      
    /*
     * init spliced alignment slides with subidx,
     * start (relative to subidx), length, and
     * strand where the strand of the second slide
     * is determined by the previous strand and the
     * distclusttype: strand stays the same on R2L
     * connections (ie. distclusttype==1) and changes
     * otherwise
     */
    j = (list[chr].cluster[clust].distclusttype[i] == 1) ?
      strand : 1-strand;

    dist.mcsa = ALLOCMEMORY(space, NULL, MultiCharSeqAlignment*, 2);
    dist.mcsa[0] = ALLOCMEMORY(space, NULL, MultiCharSeqAlignment, 1);
    dist.mcsa[1] = ALLOCMEMORY(space, NULL, MultiCharSeqAlignment, 1);
    dist.noofmcsa = 2;
      
    /* order of alignment is dependent on current direction of extension */
    bl_remappingInitAlign(space, dist.mcsa[(right) ? 0 : 1], chr, curstart, curlen, strand);
    bl_remappingInitAlign(space, dist.mcsa[(right) ? 1 : 0], dist.chr, diststart, distlen, j);
    bl_containerAdd(a, &dist);
  }
}

void bl_remappingGetRange(void *space, Container *a, matchsplitsiteclusterlist_t *L, matchsplitsiteclusterlist_t *R,
                          Uint type, Uint chr, Uint start, Uint end, unsigned char strand, Uint maxdist, unsigned char right)
{
  Uint i, u, cur, last;
  matchsplitsiteclusterlist_t *list;

  if (type == 0){
    list = L;
  }
  else {
    list = R;
  }

  u = binarySearch_left(list[chr].cluster, list[chr].noofclusters, &start,
                        cmp_matchsplitsitecluster_bin, NULL);

  for (i = u; i < list[chr].noofclusters; i++){ 
    cur = list[chr].cluster[i].median - 1; /* - 1 offset on median */
    if (cur >= start){
      if (cur > end) break;
      
      last = (type == 0) ? end : start;

      bl_remappingGetDist(space, a, L, R, type, chr, i, strand, last, end-start+1, right);
    }
  }
  qsort(a->contspace, bl_containerSize(a), sizeof(remappingclust_t), cmp_remappingclust_qsort);

  if (maxdist > 0 && bl_containerSize(a) > maxdist){
    for (i = maxdist; i < bl_containerSize(a); i++){
      bl_remappingDestructClust(bl_containerGet(a, i));
    }
    a->nextfree = maxdist;
  }
}

void bl_remappingGetAdjoint(void *space, Container *a, matchsplitsiteclusterlist_t *L, matchsplitsiteclusterlist_t *R,
                            Uint type, Uint chr, Uint clust, Uint margin, unsigned char strand, Uint maxdist, unsigned char right){
  Uint i, j, start, cur, dist, last;
  matchsplitsiteclusterlist_t *list, *other;

  if (type == 0){
    list = L;
    other = R;
  }
  else {
    list = R;
    other = L;
  }
  /* - 1 offset on median */
  start = list[chr].cluster[clust].median - 1;
  for (i = 0; i < list[chr].cluster[clust].noofadjclust; i++){
    j = list[chr].cluster[clust].adjclust[i];
    cur = other[chr].cluster[j].median - 1;
    dist = cur >= start ? cur - start : start - cur;

    if (dist <= margin){      
      last = list[chr].cluster[clust].median - 1;

      bl_remappingGetDist(space, a, L, R, 1-type, chr, j, strand, last, margin, right);
    }
  }
  qsort(a->contspace, bl_containerSize(a), sizeof(remappingclust_t), cmp_remappingclust_qsort);

  if (maxdist > 0 && bl_containerSize(a) > maxdist){
    for (i = maxdist; i < bl_containerSize(a); i++){
      bl_remappingDestructClust(bl_containerGet(a, i));
    }
    a->nextfree = maxdist;
  }
}


unsigned char bl_remappingExtractSeed(void *space, char *desc, Uint desclen, remappingseed_t *seed){
  Uint i;
  int val;
  char *err;
  stringset_t *token, *fields;

  // split description by spaces,
  // the last two entries should contain the
  // seed information
  token = tokensToStringset(space, " ", desc, desclen);
  if (token->noofstrings < 3) { 
    DBG("less than 3 tokens separated by space found (%d)\n", token->noofstrings);
    return 0;
  }

  // split second last entry (query information of seed) by colon
  i = token->noofstrings-2;
  fields = tokensToStringset(space, ":", token->strings[i].str,
                             token->strings[i].len);
  if (fields->noofstrings != 2) {
    DBG("number of fields separated by ':' is wrong (found %d != 2)\n", fields->noofstrings);
    return 0;
  }
  
  val = strtol(fields->strings[0].str, &err, 10);
  if (val < 0 || errno == ERANGE || *err != 0) {
    DBG("ustart out of range: %d\n", val);
    return 0;
  }
  seed->ustart = val;  
  
  val = strtol(fields->strings[1].str, &err, 10);
  if (val < 0 || errno == ERANGE || *err != 0) { 
   DBG("ulen out of range: %d\n", val);
   return 0;
  }
  seed->ulen = val;
  seed->uend = seed->ustart + seed->ulen - 1;
  destructStringset(space, fields);
  
  i++;
  fields = tokensToStringset(space, ":", token->strings[i].str, token->strings[i].len);
  if (fields->noofstrings != 3) {
    DBG("number of fields separated by ':' is wrong (found %d != 3)\n", fields->noofstrings);
    return 0;
  }
  
  val = strtol(fields->strings[0].str, &err, 10);
  if (val < 0 || errno == ERANGE || *err != 0) {
    DBG("chromidx out of range: %d\n", val);
    return 0;
  }
  seed->chromidx = val;

  val = strtol(fields->strings[1].str, &err, 10);
  if (val < 0 || errno == ERANGE || *err != 0){ 
    DBG("vstart out of range: %d\n", val);
    return 0;
  }
  seed->vstart = val;


  val = strtol(fields->strings[2].str, &err, 10);
  if (val < 0 || val > 1 || errno == ERANGE || *err != 0) {
    DBG("strand out of range: %d\n", val);
    return 0;
  }
  seed->strand = (unsigned char) (val & 1);
  destructStringset(space, fields);

  destructStringset(space, token);
  return 1;
}


unsigned char bl_remappingReport(void *space, List *res, fasta_t *reads, Uint k, remapping_t *nfo){
  Uint i, qrylen, ustart, uend, ulen, vstart, vend, vlen, edist, totalcover, totaledist, maxedist;
  int previdx, prevpos, nextidx, nextpos, prevstrand, nextstrand, score, totalscore;
  matchstatus_t pairStatus = QUERY;
  gread_t read;
  Gmap map;
  gmatchlist_t *list=NULL;
  Alignment *alcopy;
  MultiCharSeqAlignment *mcsa, *prev, *next;

  qrylen = bl_fastaGetSequenceLength(reads, k);
  maxedist = qrylen - floor(((double)nfo->accuracy*qrylen)/100.);  
  list = bl_gmatchlistInit(space, maxedist, 0);
  totalcover = 0; totalscore = 0; totaledist = 0;

  for (i = 0; i < bl_listSize(res); i++){
    mcsa = bl_listGetElem(res, i);

    vstart = mcsa->refstart + mcsa->al->voff;
    vlen = getValignlen(mcsa->al);
    vend = vstart + vlen - 1;

    ustart = mcsa->al->uoff;
    ulen = getUalignlen(mcsa->al);
    if (mcsa->strand == 1){
      uend = qrylen - ustart - 1;
      ustart = uend - ulen + 1;
    }
    else {
      uend = ustart + ulen - 1;
    }

    score = getAlignScore(mcsa->al, nfo->scores, nfo->indel);
    edist = getEdist(mcsa->al);

    assert(ulen >= nfo->minfragmentalignlen &&
           vlen >= nfo->minfragmentalignlen &&
           score >= nfo->minfragmentalignscore);

    totalcover += ulen;
    totalscore += score;
    totaledist += edist;

    alcopy = ALLOCMEMORY(space, NULL, Alignment, 1);
    copyAlignment(alcopy, mcsa->al);

    previdx = -1;
    prevpos = -1;
    nextidx = -1;
    nextpos = -1;
    prevstrand = -1;
    nextstrand = -1;

    if (i > 0){
      prev = bl_listGetElem(res, i - 1);
      previdx = prev->subidx;      
      if (prev->strand == 0){
        prevpos = prev->refstart + prev->al->voff + getValignlen(prev->al) - 1;
        prevstrand = '+';        
      }
      else {
        prevpos = prev->refstart + prev->al->voff;
        prevstrand = '-';
      }
    }

    if (i < bl_listSize(res) - 1){
      next = bl_listGetElem(res, i + 1);
      nextidx = next->subidx;
      if (next->strand == 0){
        nextpos = next->refstart + next->al->voff;
        nextstrand = '+';
      }
      else {
        nextpos = next->refstart + next->al->voff + getValignlen(next->al) - 1;
        nextstrand = '-';
      }
    }

    list = se_kdMatchListAdd(list, mcsa->subidx, vstart, vend, edist, score,
                             ustart, uend, .0, alcopy, mcsa->strand,
                             previdx, prevpos, prevstrand, nextidx, nextpos, nextstrand, i);
    
  }

  totalcover *= 100;
  totalcover /= qrylen;

  if (totalscore >= nfo->minsplicedalignscore &&
      totalcover >= nfo->minsplicedaligncover){    

    initRead(&read, k);
    initGmap(&map, nfo->seq, 1);  
    setReads(&map, &read, 1);

    se_setMatches(space, &read, list, maxedist, &nfo->seinfo, 0);
    reportMatch(space, &map, reads, &nfo->seinfo, pairStatus, 0);
    se_destructMatches(space, &read); 
    bl_gmatchlistDestruct(space, list);
    return 1;
  }
  else {    
    bl_gmatchlistDestruct(space, list);
    return 0;
  }  
}

void bl_remappingReportUnmapped (void *space, fasta_t *reads, Uint k, remapping_t *nfo){

  if (!bl_fastaHasMate(reads)) { 
    if (nfo->threadno > 1) pthread_mutex_lock(nfo->seinfo.mtx2);

    if (!bl_fastaHasQuality(reads)){
      fprintf(nfo->seinfo.nomatchdev, ">%s\n%s\n", 
	      bl_fastaGetDescription(reads, k),            
	      bl_fastaGetSequence(reads, k)); 
    } else {
      fprintf(nfo->seinfo.nomatchdev, "@%s\n%s\n+%s\n%s\n",           
	      bl_fastaGetDescription(reads, k), bl_fastaGetSequence(reads, k),
	      bl_fastaGetDescription(reads, k), bl_fastaGetQuality(reads, k));
    }
    fflush(nfo->seinfo.nomatchdev);
    if (nfo->threadno > 1) pthread_mutex_unlock(nfo->seinfo.mtx2);
  }
}

void bl_remapping (void *space, MultiCharSeq *seq, fasta_t *reads,
                   matchsplitsiteclusterlist_t *L, matchsplitsiteclusterlist_t *R, 
                   remapping_t *nfo){

  Uint i, j, k, chr, type, clust, remainder, prevuoff, uoff, ulen,
    vstart, margin, first, right;
  Uint aligns, maxdist, extensions;
  Container *dist;
  List *res;
  remappingclust_t *spliced, *unspliced, *distclust;
  MultiCharSeqAlignment *mcsa;
  unsigned char ret, strand, debug;
  Uint len, desclen;
  Uint maxedist;
  char *seqs[2], *curseqs[2], *desc;
  remappingseed_t seed;
  remappingstatus_t status;
  
  debug = 0;

  for (k = 0; k < reads->noofseqs; k++){
    aligns = 0;
    extensions = 0;
    maxdist = 0;
    
    /* progressbar */
    if(!nfo->mute) rm_updateProgressBar(k, nfo);

    /* clip */

    /* get seed information from descriptions */
    desc = bl_fastaGetDescription(reads, k);
    desclen = bl_fastaGetDescriptionLength(reads, k);

    ret = bl_remappingExtractSeed(space, desc, desclen, &seed);
    if (!ret || seed.chromidx >= seq->numofsequences){
      DBG("Error in parsing seed information: '%s' (ret:%d, seed.chridx:%d[%d,%d] ? %d:numofseqs)\n", desc, ret, seed.chromidx, seed.chromstart, seed.chromend, seq->numofsequences);
      exit(-1);
    }
    getMultiCharSeqIdxBounds(seq, seed.chromidx, &seed.chromstart, &seed.chromend);
    //if (seed.vstart > seed.chromend){
    //  DBG("Error in parsing seed information: position (%u) of chromosome\n", seed.vstart);
    //  exit(-1);
    //}

    if (seed.ulen == 0){
      /* update stats */
      status = NO_SEED;
      if (nfo->stat){
	if (nfo->threadno > 1) pthread_mutex_lock(nfo->seinfo.mtx3);
	rm_addStat(space, nfo->stat, status, aligns, extensions, maxdist);
	if (nfo->threadno > 1) pthread_mutex_unlock(nfo->seinfo.mtx3);
      }
      /* output unmatched */
      if (nfo->seinfo.nomatchdev){
	bl_remappingReportUnmapped(space, reads, k, nfo);
      }
      continue;
    }
    
    seqs[0] = bl_fastaGetSequence(reads, k);
    len = bl_fastaGetSequenceLength(reads, k);
    seqs[1] = charDNAcomplement(space, seqs[0], len);

    /* prefilters (eg entropy, #Ns, coverage of query by seed) */

    maxedist = len-floor(((double)nfo->accuracy*len)/100.);

    first = 1;
    right = 1;

    res = ALLOCMEMORY(space, NULL, List, 1);
    bl_listInit(res, 10, sizeof(MultiCharSeqAlignment));

    /* ============== EXTENSION ============== */
    while (1){
      dist = ALLOCMEMORY(space, NULL, Container, 1);
      bl_containerInit(dist, 100, sizeof(remappingclust_t));

      /* initializations */
      if (first){
        if (debug){
          fprintf(stderr, "seed: seed.ustart=%u, seed.uend=%u, seed.ulen=%u, qrylen=%u, seed.chromidx=%u, seed.vstart=%u, seed.strand=%u\n", seed.ustart, seed.ustart+seed.ulen-1, seed.ulen, len, seed.chromidx, seed.vstart, seed.strand);
        }
        
        chr = seed.chromidx;
        clust = -1;
        prevuoff = 0;
        uoff = 0;
        ulen = len;
        curseqs[0] = seqs[0];
        curseqs[1] = seqs[1];
        unspliced = NULL;
        
        strand = seed.strand;
        
        if (right){
          type = (strand == 0) ? 1 : 0;

          margin = len + maxedist;
          if (strand != 0){
            vstart = (seed.vstart + seed.ulen > margin) ?
               seed.vstart + seed.ulen - margin : 0;
          }
          else {
            vstart = seed.vstart;
          }
        }
        else {
          type = (strand == 0) ? 0 : 1;

          /* update qry & ref info with first elem in res (= seed region) */
          if (bl_listSize(res) > 0){
            mcsa = bl_listGetElem(res, res->first);
            assert(mcsa->strand == strand);

            /*
             *             __________      
             * qry: ------|__________|-----
             *               mcsa         
             *
             *      \----- ulen -----/
             */
            
            if (strand == 0){
              ulen = mcsa->al->uoff + getUalignlen(mcsa->al);
            }
            else {
              ulen = mcsa->al->ulen - mcsa->al->uoff;
            }
            prevuoff = len - ulen;
            curseqs[1] = &seqs[1][prevuoff];
            margin = ulen + maxedist;

            vstart = mcsa->refstart - mcsa->substart + mcsa->al->voff;
            if (strand == 0){
              vstart = (vstart + getValignlen(mcsa->al) > margin) ?
                vstart + getValignlen(mcsa->al) - margin : 0;
            }
            
            /* if entire query is mapped --> no left extension necessary */
            if (ulen - getUalignlen(mcsa->al) == 0){
              break;
            }
            
            /* delete first elem in res (= seed region) */
            mcsa = bl_listUnlink(res, res->first, NULL);  
            wrapMultiCharSeqAlignment(space, mcsa);
            FREEMEMORY(space, mcsa);
          }
          else {
            margin = 2 * (len + maxedist);
            vstart = (seed.vstart > len + maxedist) ? seed.vstart - len - maxedist : 0;
          }
        }

        // possible optimizations:
        // - during right extension: switch to left extension if
        //   remainder of query right of seed is < minfragmentlen
        //   (assumes that there is no split within the seed region)

        
        /*
         * find R or L clusters with median in interval:
         * [seed.vstart, seed.vstart+margin) if seed on plus
         * (seed.vstart-margin, seed.vstart] otherwise
         */
        bl_remappingGetRange(space, dist, L, R, type, chr, vstart, vstart + margin - 1, strand, nfo->maxdist, right);
    
        first = 0;
      }
      else {
        bl_remappingGetAdjoint(space, dist, L, R, type, chr, clust, margin, strand, nfo->maxdist, right);
      }

      if (debug){
        fprintf(stderr, "%s: type=%u, chr=%u, clust=%d, strand=%u, uoff=%u, ulen=%u, vstart=%u, margin=%u, noofdist=%u\n",
                (right) ? "right" : "left", type, chr, clust, strand, uoff, ulen, vstart, margin, bl_containerSize(dist));
      }

      /* reporting stuff */
      if (0){
        for (i = 0; i < bl_containerSize(dist); i++){
          distclust = bl_containerGet(dist, i);
          fprintf(stderr, "dist: %u, type=%u, chr=%u, clust=%d, a=%u, b=%u, median=%u\n", i, 
                  distclust->type, distclust->chr, distclust->clust, distclust->cluster->a - 1,
                distclust->cluster->b - 1, distclust->cluster->median - 1);// - 1 offset due to clusterlist
          for (j = 0; j < distclust->noofmcsa; j++){
            fprintf(stderr, " mcsa[%u]: chr=%u, start=%u, end=%u, len=%u, strand=%u\n", 
                    j, distclust->mcsa[j]->subidx, distclust->mcsa[j]->refstart,
                    distclust->mcsa[j]->refstart + distclust->mcsa[j]->reflen - 1,
                    distclust->mcsa[j]->reflen, distclust->mcsa[j]->strand);
          }
        }
      }

      /* extend by spliced alignment */
      spliced = NULL;      
      if (bl_containerSize(dist) > 0){
        
        /* update ref & query info in alignments */
        for (i = 0; i < bl_containerSize(dist); i++){
          distclust = bl_containerGet(dist, i);
          bl_remappingUpdateAlignSeqs(space, seq, curseqs, ulen, distclust);
        }

        /* calculate spliced alignments => get best one */
        spliced = bl_remappingAlign(space, curseqs, dist, nfo);
        aligns += bl_containerSize(dist);
	if (maxdist < bl_containerSize(dist))
	  maxdist = bl_containerSize(dist);

        /* convert from relative query positions to absolute ones */
        if (spliced != NULL){
          bl_remappingUpdateAlignOff(space, seqs, len, spliced, prevuoff, right);
        }
      }
      
      /* extend by unspliced alignment ==> only necessary once (at the beginning) */
      if (unspliced == NULL){

        /* no valid spliced alignment during initial extension */
        if (spliced == NULL && bl_listSize(res) == 0){
          /*
           * continue with left extension if only unspliced
           * extension is possible since it will be recalculated
           * during left extension anyway
           */
          if (right){
            bl_containerDestruct(dist, bl_remappingDestructClust);
            FREEMEMORY(space, dist);
        
            right = 0;
            first = 1;
            continue;
          }
          /*
           * break if no spliced left extension is possible and no
           * right extension was successful and hence the entire query
           * would only be mapped unspliced which was already done in
           * the initial mapping => do not waste time on this!
           */
          else {
            break;
          }
        }

        unspliced = bl_remappingUnsplicedAlign(space, seq, curseqs, ulen, chr, vstart,
                                             vstart + margin - 1, strand, nfo);
	aligns++;

        if (unspliced != NULL){
          bl_remappingUpdateAlignOff(space, seqs, len, unspliced, prevuoff, right);
        }
      }                            

      /*
       * cases:
       * - spliced set and unspliced not set or is worse than spliced => extend further
       * - unspliced set and spliced not or is worse than unspliced => last extension
       * - unspliced and spliced not set => no extension
       */

      /* if best is set and spliced is better than unspliced */
      if (spliced != NULL && (unspliced == NULL || spliced->score > unspliced->score)){

        if (debug){
          fprintf(stderr, "%s: spliced extension: qryoff=%u, qrylen=%u, vstart=%u, margin=%u, strand=%u\n",
                  (right)? "right" : "left", uoff, ulen, vstart, margin, strand);
          bl_remappingPrintClust(stderr, spliced, nfo);
        }
	extensions++;

        /* 
         * keep information of best split and 
         * update variables
         */
        type = spliced->type;
        chr = spliced->chr;
        clust = spliced->clust;
        if (right) {
          /*
           * in right extension:
           *                                     /-rem-\
           *           __________       ________  
           * qry: ----|__________|-----|________|-------
           *             mcsa[0]         mcsa[1]
           *
           *      \----- off ----/\------- ulen -------/
           */

          strand = spliced->mcsa[1]->strand;

          /* update qry info (with spliced->mcsa[0]) */
          if (spliced->mcsa[0]->strand == 0){
            ulen = spliced->mcsa[0]->al->ulen - spliced->mcsa[0]->al->uoff - getUalignlen(spliced->mcsa[0]->al);
          }
          else {
            ulen = spliced->mcsa[0]->al->uoff;
          }
          uoff = len - ulen;
          curseqs[0] = &seqs[0][uoff];

          margin = ulen + maxedist;
          
          /* calculate remainder of query (with spliced->mcsa[1]) */
          if (spliced->mcsa[1]->strand == 0){
            remainder = spliced->mcsa[1]->al->ulen - spliced->mcsa[1]->al->uoff - getUalignlen(spliced->mcsa[1]->al);
          }
          else {
            remainder = spliced->mcsa[1]->al->uoff;
          }

          /* update ref info (with spliced->mcsa[1]) */
          vstart = spliced->mcsa[1]->refstart - spliced->mcsa[1]->substart + spliced->mcsa[1]->al->voff;
          if (strand != 0){
            vstart = (vstart + getValignlen(spliced->mcsa[1]->al) > margin) ?
              vstart + getValignlen(spliced->mcsa[1]->al) - margin : 0;
          }
        }
        else {

          /*
           * in left extension:
           *      /-rem-\
           *             __________       ________  
           * qry: ------|__________|-----|________|------
           *               mcsa[0]         mcsa[1]
           *
           *      \--------- ulen ------/\---- uoff ----/
           */
          strand = spliced->mcsa[0]->strand;

          /* update qry info (with spliced->mcsa[1]) */
          if (spliced->mcsa[1]->strand == 0){
            ulen = spliced->mcsa[1]->al->uoff;
          }
          else {
            ulen = spliced->mcsa[1]->al->ulen - spliced->mcsa[1]->al->uoff - getUalignlen(spliced->mcsa[1]->al);
          }
          uoff = len - ulen;
          curseqs[1] = &seqs[1][uoff];
          margin = ulen + maxedist;

          /* calculate remainder of query (with spliced->mcsa[0]) */
          if (spliced->mcsa[0]->strand == 0){
            remainder = spliced->mcsa[0]->al->uoff;
          }
          else {
            remainder = spliced->mcsa[0]->al->ulen - spliced->mcsa[0]->al->uoff - getUalignlen(spliced->mcsa[0]->al);
          }

          /* update ref info (with spliced->mcsa[0]) */
          vstart = spliced->mcsa[0]->refstart - spliced->mcsa[0]->substart + spliced->mcsa[0]->al->voff;
          if (strand == 0){
            vstart = (vstart + getValignlen(spliced->mcsa[0]->al) > margin) ?
              vstart + getValignlen(spliced->mcsa[0]->al) - margin : 0;
          }
        }
        

        /* 
         * store first alignment at end in case of right
         * and second alignment at begin in case of left extension
         */
        if (right){
          bl_listInsert(res, res->last, spliced->mcsa[0]);
          FREEMEMORY(space, spliced->mcsa[0]);
          spliced->mcsa[0] = NULL; 
        }
        else {
          bl_listInsert(res, -1, spliced->mcsa[1]);
          FREEMEMORY(space, spliced->mcsa[1]);
          spliced->mcsa[1] = NULL;          
        }

        if (debug){
          fprintf(stderr, "%s: ulen=%u, uoff=%u, remainder=%u\n", (right)? "right" : "left",
                  ulen, uoff, remainder);
        }

        /* store other alignment if remainder is too short for another split */
        if (remainder < nfo->minfragmentalignlen){
          if (right){
            bl_listInsert(res, res->last, spliced->mcsa[1]);
            FREEMEMORY(space, spliced->mcsa[1]);
            spliced->mcsa[1] = NULL;
          }
          else {
            bl_listInsert(res, -1, spliced->mcsa[0]);
            FREEMEMORY(space, spliced->mcsa[0]);
            spliced->mcsa[0] = NULL;
          }   

          if (unspliced != NULL){
            bl_remappingDestructClust(unspliced);
            FREEMEMORY(space, unspliced);
            unspliced = NULL;
          }   

          /* change from right extension to left one */
          if (right){ 
            bl_containerDestruct(dist, bl_remappingDestructClust);
            FREEMEMORY(space, dist);
            
            right = 0;
            first = 1;
            continue;
          }
          
          break;
        }

        /*
         * store second alignment as unspliced for next extension step in case of right
         * and first alignment in case of left extension
         */
        if (unspliced != NULL){
          bl_remappingDestructClust(unspliced);
        }      
        else {
          unspliced = ALLOCMEMORY(space, NULL, remappingclust_t, 1);
          unspliced->chr = 0;
          unspliced->type = 0;
          unspliced->clust = -1;
          unspliced->cluster = NULL;
        }
        unspliced->mcsa = ALLOCMEMORY(space, NULL, MultiCharSeqAlignment*, 1);
        unspliced->noofmcsa = 1;
        if (right){
          unspliced->mcsa[0] = spliced->mcsa[1];
          unspliced->score = getAlignScore(spliced->mcsa[1]->al, nfo->scores, nfo->indel);
          spliced->mcsa[1] = NULL;
        }
        else {
          unspliced->mcsa[0] = spliced->mcsa[0];
          unspliced->score = getAlignScore(spliced->mcsa[0]->al, nfo->scores, nfo->indel);
          spliced->mcsa[0] = NULL;
        }

        /* prepare next extension */
        prevuoff = uoff;
        bl_containerDestruct(dist, bl_remappingDestructClust);
        FREEMEMORY(space, dist);

        continue;
      }   
      else if (unspliced != NULL){
        
        if (debug){
          fprintf(stderr, "%s: unspliced extension: qryoff=%u, qrylen=%u, vstart=%u, margin=%u, strand=%u\n",
                  (right)? "right" : "left", uoff, ulen, vstart, margin, strand);
          bl_remappingPrintClust(stderr, unspliced, nfo);
        }
	extensions++;

        /* store unspliced alignment */
        bl_listInsert(res, (right) ? res->last : -1, unspliced->mcsa[0]);
        FREEMEMORY(space, unspliced->mcsa[0]);
        unspliced->mcsa[0] = NULL;
      }
      else {
        if (debug){
          fprintf(stderr, "%s: no extension: qryoff=%u, qrylen=%u, vstart=%u, margin=%u, strand=%u\n",
                  (right)? "right" : "left", uoff, ulen, vstart, margin, strand);
        }
      }

      if (unspliced != NULL) {
        bl_remappingDestructClust(unspliced);
        FREEMEMORY(space, unspliced);
        unspliced = NULL;
      }

      /* change from right extension to left one */
      if (right){        
        bl_containerDestruct(dist, bl_remappingDestructClust);
        FREEMEMORY(space, dist);

        right = 0;
        first = 1;
        continue;
      }

      break;
    }

    /* report stuff */
    bl_listSweep(res);
    
    if (bl_listSize(res) > 1){
      ret = bl_remappingReport(space, res, reads, k, nfo);
      /* update stats */
      status = ret ? REMAPPED : PURGED;
      if (nfo->stat){
	if (nfo->threadno > 1) pthread_mutex_lock(nfo->seinfo.mtx3);
	rm_addStat(space, nfo->stat, status, aligns, extensions, maxdist);
	if (nfo->threadno > 1) pthread_mutex_unlock(nfo->seinfo.mtx3);
      }
      if (!ret){
	/* output unmatched */
	if (nfo->seinfo.nomatchdev){
	  bl_remappingReportUnmapped(space, reads, k, nfo);
	}	
      }
    }
    else {
      /* update stats */
      status = (aligns == 0) ? NO_DIST : UNSPLICED;
      if (nfo->stat){
	if (nfo->threadno > 1) pthread_mutex_lock(nfo->seinfo.mtx3);
	rm_addStat(space, nfo->stat, status, aligns, extensions, maxdist);
	if (nfo->threadno > 1) pthread_mutex_unlock(nfo->seinfo.mtx3);
      }
      /* output unmatched */
      if (nfo->seinfo.nomatchdev){
	bl_remappingReportUnmapped(space, reads, k, nfo);
      }
    }

    /* clear stuff */    
    bl_containerDestruct(dist, bl_remappingDestructClust);
    FREEMEMORY(space, dist);

    for (i = 0; i < bl_listSize(res); i++){
      mcsa = bl_listGetElem(res, i);
      wrapMultiCharSeqAlignment(space, mcsa);
    }
    bl_listDestruct(res, NULL);
    FREEMEMORY(space, res);
    
    FREEMEMORY(space, seqs[1]);
  }
}

#ifdef REMAPPINGTEST

unsigned char mute = 0;
char *ntcode;
pthread_mutex_t updatemtx;

void*
checkclock(void *args) {
  checkthreadinfo_t *t;

  sleep(2);
  cursorVisible();
  t = (checkthreadinfo_t*) args;
  initProgressBarVT();

  while (pthread_mutex_trylock(&updatemtx) != 0) {
    progressBarVT("reads remapped.", t->noofseqs, (*t->counter), 25);
  }

  cursorVisible();
  fprintf(stderr, "\n");
  return NULL;
}

int main(int argc, char **argv) {
  
  remapping_t nfo, *thnfo;
  
  manopt_optionset optset;
  manopt_arg *unflagged; 
  manopt_arg *queries;
  manopt_arg *dbfilenames;
  manopt_intconstraint threadconstraint;
  manopt_intconstraint accuracyconstraint;

  Uint counter=0, desclen;
  matchsplitsiteclusterlist_t *L = NULL, *R = NULL;
  matchfileindex_t *index;
  matchfile_t **files = NULL; 
  unsigned char gzip = 0, indexreads=1;
  char *version, *desc;
  int i;
  Uint nchr;
  Uint prefixlen=0;
  double difmatch;
  time_t startmatch, endmatch;
  fasta_t  **chopsuey;  
  pthread_t *threads;
  pthread_t clockthread;
  checkthreadinfo_t ch_info;

  threadconstraint.max = 3000;
  threadconstraint.min = 1;
  accuracyconstraint.max = 100;
  accuracyconstraint.min = 0;

  /* set default settings */
  rm_setdefault(&nfo);

  /* init stat struct => comment out to disable */
  //nfo.stat = ALLOCMEMORY(nfo.space, NULL, remappingstat_t, 1);
  //rm_initStat(nfo.stat);

  /* capture command */
  nfo.seinfo.cmdline = malloc(strlen(argv[0])+1);
  strcpy(nfo.seinfo.cmdline, argv[0]);
  
  for(i = 1; i < argc; i++) {    
    nfo.seinfo.cmdline = realloc(nfo.seinfo.cmdline, strlen(nfo.seinfo.cmdline) + 
    strlen(argv[i]) + 2);
    strcat(nfo.seinfo.cmdline," ");
    strcat(nfo.seinfo.cmdline,argv[i]);
  }
  
  initIUPAC(1,1); 
  version = getNiceSVNVersion(VERSION);
  manopt_initoptionset(&optset, argv[0], NULL, 
      "Remapping of unmapped reads\n",
      "SEGEMEHL is free software for non-commercial use \n  (C) 2012 Bioinformatik Leipzig\n",
      version,
      "Please report bugs to christian@bioinf.uni-leipzig.de"); 
  manopt(&optset, LISTOPT, 1, 'd', "database", 
	 "list of path/filename(s) of database sequence(s)", "<file> [<file> ...]", 
	 NULL, NULL);
  manopt(&optset, LISTOPT, 1, 'q', "query", 
	 "path/filename of alignment file", "<file> [<file> ...]", NULL, NULL); 
  manopt(&optset, REQSTRINGOPT, 0, 'o', "outfile", 
      "outputfile", "<string>", NULL, &nfo.outfile);
  manopt(&optset, STRINGOPT, 1, 'r', "remapfilename", 
	 "filename for reads to be remapped", "<file>", NULL, &nfo.remapfile);   
  manopt(&optset, REQSTRINGOPT, 0, 'u', "nomatchfilename", 
      "filename for unmatched reads", "<file>", NULL, &nfo.nomatchfile); 
  manopt(&optset, REQUINTOPT, 0, 't', "threads", 
	 "start <n> threads for remapping", "<n>",
         &threadconstraint, &nfo.threadno);
  manopt(&optset, FLAG, 0, 's', "silent", 
      "shut up!", NULL, NULL, &nfo.mute);
  manopt(&optset, REQUINTOPT, 0, 'A', "accuracy", 
      "min percentage of matches per read in semi-global alignment", "<n>", &accuracyconstraint, &nfo.accuracy);
  manopt(&optset, REQUINTOPT, 0, 'W', "minsplicecover", 
      "min coverage for spliced transcripts", "<n>", &accuracyconstraint, &nfo.minsplicedaligncover);
  manopt(&optset, REQUINTOPT, 0, 'U', "minfragscore", 
      "min score of a spliced fragment", "<n>", NULL, &nfo.minfragmentalignscore);
  manopt(&optset, REQUINTOPT, 0, 'Z', "minfraglen", 
      "min length of a spliced fragment", "<n>", NULL, &nfo.minfragmentalignlen);
  manopt(&optset, REQUINTOPT, 0, 'M', "maxdist", 
      "max number of distant sites to consider, 0 to disable", "<n>", NULL, &nfo.maxdist);

  unflagged = manopt_getopts(&optset, argc, argv);
     
  if(unflagged->noofvalues > 1) { 
    manopt_help(&optset, "unknown argument(s)\n");
  }

  MSG("reading database sequences.\n"); 

  dbfilenames = manopt_getarg(&optset, 'd', "database");
  nfo.fasta = bl_fastxGetSet(nfo.space, dbfilenames->values, 
      dbfilenames->noofvalues, 1, 0, 0, 1);
  
  for(i=0; i < nfo.fasta->noofseqs; i++) {
    desclen = bl_fastaGetDescriptionLength(nfo.fasta, i);
    desc = strclip(nfo.space, bl_fastaGetDescription(nfo.fasta, i), &desclen);
    FREEMEMORY(nfo.space, nfo.fasta->seqs[i]->description);
    nfo.fasta->seqs[i]->description = desc;
    nfo.fasta->seqs[i]->descrlen = desclen;
  }

  NFO("%d database sequences found.\n", nfo.fasta->noofseqs);
  MSG("reading query files.\n");

  queries = manopt_getarg(&optset, 'q', "query");
  if(queries->noofvalues > 30) {
    manopt_help(&optset, "currently no more than 30 query files allowed\n");
  }

  ntcode  = getNTcodekey(nfo.space);
  files = ALLOCMEMORY(nfo.space, NULL, matchfile_t*, queries->noofvalues);
  
  //using index structure only to carry the chr idx
  index = bl_matchfileInitIndex(nfo.space);
  
  nchr = nfo.fasta->noofseqs;
  for(i=0; i < nchr; i++) { 
    bl_matchfileIndexAddChrom(index, bl_fastaGetDescription(nfo.fasta, i));
  }

  for(i=0; i < queries->noofvalues; i++) {

    files[i] = ALLOCMEMORY(nfo.space, NULL, matchfile_t, 1);  
    files[i]->fmt = 0;
    files[i]->index = index;
    files[i]->filename = queries->values[i];

    prefixlen = bl_fileprefixlen(files[i]->filename);

    gzip = 0;
    if(strncmp(&files[i]->filename[prefixlen], ".gz", 3) == 0 || 
        strncmp(&files[i]->filename[prefixlen], ".gzip", 3) == 0) {
      gzip = 1;
    }

    files[i]->gzip = gzip;
  }

  L = ALLOCMEMORY(nfo.space, NULL, matchsplitsiteclusterlist_t,nchr);
  memset(L, 0, sizeof(matchsplitsiteclusterlist_t)*nchr);
  R = ALLOCMEMORY(nfo.space, NULL, matchsplitsiteclusterlist_t, nchr);
  memset(R, 0, sizeof(matchsplitsiteclusterlist_t)*nchr);
 
  for(i=0; i < queries->noofvalues; i++) { 
    bl_matchfileRealignScanFileNew(nfo.space, files[i], NULL, nfo.fasta, 255, &L, &R, &nchr, 1,0); 
  }

  bl_matchLinkAdjoinedCluster(nfo.space, L, R, nchr);
  bl_matchLinkDistCluster (nfo.space, R, L, nchr); 
    
  nfo.seq = concatCharSequences(nfo.space, nfo.fasta->seqs, nfo.fasta->noofseqs, (char)126, (char)127);
  if(indexreads) {
    nfo.reads = bl_fastxGetSet(nfo.space, &nfo.remapfile, 1, 1, 0, 1, nfo.threadno);
  } else {
    nfo.reads = bl_fastxRead(nfo.space, nfo.reads, nfo.remapfile, 1, 0, 0, 0, 0, bl_fastxAdd);
  }
  NFO("%d unmatched read sequences found.\n", nfo.reads->noofseqs);
  
  if (nfo.threadno > nfo.reads->noofseqs) {
    NFO("more threads than unmapped reads. Exit forced\n", NULL);
    exit(EXIT_FAILURE);
  }
  
  nfo.L = L;
  nfo.R = R;
  nfo.seinfo.outfile = nfo.outfile;
  nfo.seinfo.fasta = nfo.fasta;
  nfo.seinfo.rep_type = 15;  
  nfo.seinfo.threadno = nfo.threadno;

  se_registerOutputDevice(nfo.space, &nfo.seinfo);

  if(nfo.nomatchfile != NULL)
    nfo.seinfo.nomatchdev = fopen(nfo.nomatchfile, "w");

  if (nfo.threadno > 1){
    nfo.counter = &counter;
    
    if(indexreads) {
      chopsuey = bl_fastxChopIndex(nfo.space, nfo.reads, nfo.threadno);
    } else {
      chopsuey = bl_fastaChop(nfo.space, nfo.reads, nfo.threadno);
    }
    thnfo = ALLOCMEMORY(nfo.space, NULL, remapping_t, nfo.threadno);
    threads = ALLOCMEMORY(nfo.space, NULL, pthread_t, nfo.threadno);
    ch_info.noofseqs = nfo.reads->noofseqs;
    ch_info.counter = &counter;
    
    if (!nfo.mute) {
      pthread_mutex_init(&updatemtx, NULL);
      pthread_mutex_lock(&updatemtx);
      pthread_create(&clockthread, NULL, checkclock, &ch_info);
    }
    
    for(i=0; i < nfo.threadno; i++) {
      NFO("%d reads in thread %d.\n", chopsuey[i]->noofseqs, i);
    }
    
    time (&startmatch);
    
    for(i=0; i < nfo.threadno; i++) {
      memmove(&thnfo[i], &nfo, sizeof(remapping_t));
      thnfo[i].reads = chopsuey[i];
      thnfo[i].threadid = i;
      pthread_create(&threads[i], NULL, remappingworker, &thnfo[i]);
    }
    
    for(i=0; i < nfo.threadno; i++) {
      pthread_join(threads[i], NULL); 
    } 
      
    if(!nfo.mute) {
      /*notification via mutex - why use signals?*/
      pthread_mutex_unlock(&updatemtx);
      pthread_join(clockthread, NULL);
    }
    
    fflush(nfo.seinfo.dev);
    time (&endmatch);
    difmatch = difftime (endmatch, startmatch);
    NFO("threaded remapping has taken %f seconds.\n", difmatch);
      
    for (i=0; i < nfo.threadno; i++) {
      bl_fastxDestructSequence(nfo.space, chopsuey[i]);
      bl_fastxDestructChunkIndex(nfo.space, chopsuey[i]);
      FREEMEMORY(nfo.space, chopsuey[i]);
    }
      
    FREEMEMORY(nfo.space, chopsuey);
    FREEMEMORY(nfo.space, thnfo);
    FREEMEMORY(nfo.space, threads);
  }
  else {
    initProgressBarVT();
    time (&startmatch); 
    bl_remapping(nfo.space, nfo.seq, nfo.reads, nfo.L, nfo.R, &nfo);
    time (&endmatch);
    difmatch = difftime (endmatch, startmatch);
    NFO("remapping has taken %f seconds.\n", difmatch); 
  }

  if (nfo.stat){
    for (i = 0; i < nfo.stat->n; i++){
      fprintf(stderr, "%d\t%u\t%u\t%u\n", nfo.stat->status[i],
	      nfo.stat->aligns[i], nfo.stat->extensions[i],
	      nfo.stat->maxdist[i]);
    }
    NFO("stats of %d elements\n", nfo.stat->n);
  }

  if (nfo.outfile)
    fclose(nfo.seinfo.dev);

  if(nfo.nomatchfile)
    fclose(nfo.seinfo.nomatchdev);

  bl_fastaDestruct(nfo.space, nfo.reads);
  FREEMEMORY(nfo.space, nfo.reads);
  
  destructMultiCharSeq(nfo.space, nfo.seq);

  if (nfo.stat){
    rm_destructStat(nfo.space, nfo.stat);
    FREEMEMORY(nfo.space, nfo.stat);
  }
  
  FREEMEMORY(nfo.space, nfo.seinfo.mtx);
  FREEMEMORY(nfo.space, nfo.seinfo.mtx2);
  FREEMEMORY(nfo.space, nfo.seinfo.mtx3);
  FREEMEMORY(nfo.space, nfo.seinfo.cmdline);
  
  bl_matchDestructMatchsplitsiteclusterlist(nfo.space, L, nchr); 
  bl_matchDestructMatchsplitsiteclusterlist(nfo.space, R, nchr); 
  FREEMEMORY(nfo.space, L);
  FREEMEMORY(nfo.space, R); 

  bl_fastaDestruct(nfo.space, nfo.fasta);
  FREEMEMORY(nfo.space, nfo.fasta);
  bl_matchfileDestructIndex(nfo.space, index);
  FREEMEMORY(nfo.space, index);

  if(files) {
    for(i=0; i < queries->noofvalues; i++) { 
      FREEMEMORY(nfo.space, files[i]);
    }
    FREEMEMORY(nfo.space, files);
  }

  manopt_destructoptionset(&optset);
  manopt_destructarg(unflagged);
  FREEMEMORY(nfo.space, unflagged);

  FREEMEMORY(nfo.space, ntcode);
  FREEMEMORY(nfo.space, version);

  return EXIT_SUCCESS;
}



#endif
