#ifndef REMAPPING_H
#define REMAPPING_H

/**
 *
 *  remapping.h
 *  remapping of unmapped reads
 * 
 *  @author Christian Otto, christain@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date Fri Oct 12 09:55:49 CEST 2012
 *
 */

/*
 * SVN
 * Revision of last commit: $Rev: 400 $
 * Author: $Author: steve $
 * Date: $Date: 2013-09-11 03:46:50 -0400 (Wed, 11 Sep 2013) $
 * Id: $Id: remapping.h 400 2013-09-11 07:46:50Z steve $
 * Url: $URL: http://www2.bioinf.uni-leipzig.de/svn5/segemehl/libs/remapping.h $
 */

#include <stdlib.h>
#include <stdio.h>
#include "basic-types.h"

typedef struct {
  Uint ustart;
  Uint uend;
  Uint ulen;
  Uint vstart;
  Uint chromstart;
  Uint chromend;
  Uint chromidx;
  unsigned char strand;
} remappingseed_t;

typedef struct {
  matchsplitsitecluster_t *cluster;
  Uint type;
  Uint chr;
  int clust;
  Uint cnt;
  MultiCharSeqAlignment **mcsa;
  Uint noofmcsa;
  int score;
} remappingclust_t;

typedef enum {
  REMAPPED,  /* successfully remapped */
  PURGED,    /* discarded due to error boundaries */
  NO_SEED,   /* no best seed information given */
  NO_DIST,   /* no splice sites in proximity found */
  UNSPLICED  /* no extension possible */
} remappingstatus_t;

typedef struct {
  remappingstatus_t *status;
  Uint *aligns;
  Uint *extensions;
  Uint *maxdist;
  Uint n;
} remappingstat_t;

typedef struct { 
  void *space;
  char *remapfile;
  char *outfile;
  char *nomatchfile;
  matchsplitsiteclusterlist_t *L;
  matchsplitsiteclusterlist_t *R; 
  MultiCharSeq *seq;
  fasta_t *fasta;
  fasta_t *reads;
  Uint threadno;
  Uint threadid;
  unsigned char mute;    
  Uint *counter;
  int scores[2];
  int indel;
  int transition;  
  Uint minfragmentalignlen;
  int  minfragmentalignscore;
  Uint minsplicedaligncover;
  Uint minsplicedalignscore;
  int accuracy;
  Uint maxdist;
  segemehl_t seinfo;
  remappingstat_t *stat;
} remapping_t;


inline static void
rm_setdefault(remapping_t *info) { 
  info->space = NULL;
  info->remapfile = NULL;
  info->outfile = NULL;
  info->nomatchfile = NULL;
  info->L = NULL;
  info->R = NULL;
  info->seq = NULL;
  info->fasta = NULL;
  info->reads = NULL;
  info->threadno = 1;
  info->threadid = 0;
  info->mute = 0;
  info->counter = 0;
  info->scores[0] = 1;
  info->scores[1] = -2;
  info->indel = -2;
  info->transition = -10;
  info->minfragmentalignlen = 5;
  info->minfragmentalignscore = 5;
  info->minsplicedaligncover = 80;
  info->minsplicedalignscore = 5+18;
  info->accuracy = 90;
  info->maxdist = 100;
  se_setdefault(&info->seinfo);
  info->stat = NULL;
}


void bl_remapping (void *space, MultiCharSeq *seq, fasta_t *reads, matchsplitsiteclusterlist_t *L, matchsplitsiteclusterlist_t *R, remapping_t *nfo);
unsigned char bl_remappingReport(void *space, List *res, fasta_t *reads, Uint k, remapping_t *nfo);
void bl_remappingReportUnmapped(void *space, fasta_t *reads, Uint k, remapping_t *nfo);
unsigned char bl_remappingExtractSeed(void *space, char *desc, Uint desclen, remappingseed_t *seed);

void bl_remappingGetRange(void *space, Container *a, matchsplitsiteclusterlist_t *L, matchsplitsiteclusterlist_t *R,
			  Uint type, Uint chr, Uint start, Uint end, unsigned char strand, Uint maxdist, unsigned char right);
void bl_remappingGetAdjoint(void *space, Container *a, matchsplitsiteclusterlist_t *L, matchsplitsiteclusterlist_t *R,
                            Uint type, Uint chr, Uint clust, Uint margin, unsigned char strand, Uint maxdist, unsigned char right);
void bl_remappingGetDist(void *space, Container *a, matchsplitsiteclusterlist_t *L, matchsplitsiteclusterlist_t *R, Uint type, Uint chr, Uint clust, unsigned char strand, Uint last, Uint margin, unsigned char right);

remappingclust_t *bl_remappingUnsplicedAlign(void *space, MultiCharSeq *seq, char **seqs, Uint len, Uint chr, Uint start, Uint end, unsigned char strand, remapping_t *nfo);
remappingclust_t *bl_remappingAlign(void *space, char **seqs, Container *dist, remapping_t *nfo);

void bl_remappingUpdateAlignSeqs(void *space, MultiCharSeq *seq, char **seqs, Uint qrylen, remappingclust_t *cur);
void bl_remappingUpdateAlignOff(void *space, char **seqs, Uint qrylen, remappingclust_t *cur, Uint offset, Uint right);

void bl_remappingPrintAlign(FILE *dev, MultiCharSeqAlignment *mcsa, remapping_t *nfo);
void bl_remappingPrintClust(FILE *dev, remappingclust_t *clust, remapping_t *nfo);
void bl_remappingDestructClust(void *elem);
void bl_remappingInitAlign(void *space, MultiCharSeqAlignment *mcsa, Uint subidx, Uint start, Uint len, unsigned char strand);

#endif /* REMAPPING_H */
