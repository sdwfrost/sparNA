#ifndef REALIGN_H
#define REALIGN_H

/*
 *
 *	realign.h
 *  realignment
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 07/10/2012 10:01:20 PM CEST  
 *
 */

#include "basic-types.h"
#include "matchfiles.h"

typedef struct{
  Uint pos;
  char trans;
} matchsplitsearchkey_t;

typedef struct {
  Uint start;
  Uint end;
  Uint distchr;
  Uint distpos;
  Uint adjoint;
  char trans;
  char *str;
  unsigned char *bookkeeper;
} matchlistelem_t; 

typedef struct{
  Uint chridx;
  Uint pos;
  Uint cnt;
  Uint trans;
} matchsplitsite_t;

typedef struct{
  Uint a;
  Uint b;
  Uint median;
  Uint cnt;
  Uint noofrealigns;
  char *str;
#ifdef CHECKLINKS
  Uint trans;
  Uint lnkcnt;
#endif
  
  Uint *distchr;
  Uint *distpos;
  char *disttrans;
  Uint *distcnt;
  Uint *realigned;
  Uint distsites;
  

  Uint *distclust;
  Uint *distclustchr;
  Uint *distclustcnt;
  Uint *distclustrealigned;
  char *distclusttype;
  Uint noofdistclust;
  double *emat;
  
#if defined PARANOID && defined CHECKLINKS
  char *rpos;
  Uint *linked;
#endif
  
  int16_t* adjoint;
  Uint* adjointcnt;
  Uint noofadjoints;

  Uint *adjclust;
  Uint *adjclustweights;
  Uint noofadjclust;
  
} matchsplitsitecluster_t;

typedef struct{

  matchsplitsitecluster_t* cluster;
  Uint noofclusters;

} matchsplitsiteclusterlist_t;


typedef struct { 
  matchsplitsiteclusterlist_t *L;
  matchsplitsiteclusterlist_t *R; 
  void *space;
  FILE *normdev;
  FILE *transdev;
  FILE *realigndev;
  MultiCharSeq *seq;
  fasta_t *fasta;
  unsigned char mute;
  Uint threadno;
  int scores[2];
  int indel;
  int transition;  
  Uint minfragmentalignlen;
  int  minfragmentalignscore;
  Uint minsplicedaligncover;
  Uint minsplicedalignscore;
  int accuracy;
  char *splitfile;
  char *transfile;
  char *outfile;
  int maxdist;
} realign_t;


inline static void
ra_setdefault(realign_t *info) { 
  info->L = NULL;
  info->R = NULL;
  info->space = NULL;
  info->normdev = stderr;
  info->transdev = stderr;
  info->realigndev = stdout;
  info->seq = NULL;
  info->fasta = NULL;
  info->threadno = 1;
  info->mute = 0;
  info->scores[0] = 1;
  info->scores[1] = -2;
  info->indel = -2;
  info->transition = -10;
  info->minfragmentalignlen = 20;
  info->minfragmentalignscore = 18;
  info->minsplicedaligncover = 60;
  info->minsplicedalignscore = 2*18;
  info->accuracy = 80;
  info->splitfile = NULL;
  info->transfile = NULL;
  info->outfile = NULL;
}

typedef struct {
  int begin;
  int stop;
  int *e2;
  char *seq;
  matchsplitsitecluster_t *T;
  Alignment ***aligns;
  char ***refseq; 
  Uint **reflens; 
  Uint **refstrand; 
  fasta_t *set;
  char *fwd; 
  char *rev;
  Uint ovhglen;
  matchfileRec_t r;
  Uint chromidx;
  Uint seqlen;
  Uint start;
  void *space;/*??*/
  int **lmr;
  int **lmv;
  int **lmc;
  char left;
  Uint locreflen;
  Uint median;
  Uint interval;
  Uint sinterval;
  int oven;
  } realignthread;

typedef struct {
  void *space; 
  FILE *realigndev; 
  List *rlist; 
  matchsplitsitecluster_t *T; 
  char left; 
  fasta_t *set; 
  Uint chromidx; 
  unsigned char fmt; 
  Uint interval;
} splicesitethread;

typedef struct {
  int begin; 
  int stop; 
  matchlistelem_t* arr; 
  unsigned char fmt;  
  matchsplitsitecluster_t *T; 
  char left; 
  fasta_t *set; 
  Uint chromidx;
  Uint noofrealigns/*??*/;
  FILE *realigndev;
  List *rlist;
  Uint interval;
  void *space;
  Uint sinterval;
} readthread;


#define LEFTLIST ((unsigned char) 1)
#define RIGHTLIST ((unsigned char) 1 << 1)
#define LEFTSPLIT ((unsigned char) 1 << 2)
#define RIGHTSPLIT ((unsigned char) 1 << 3)

Uint cmp_matchsplitsitecluster_bin(Uint a, void *data, void *key, void *nfo);
void bl_matchLinkAdjoinedCluster(void *space, matchsplitsiteclusterlist_t *L, matchsplitsiteclusterlist_t *R, Uint nchr);
void bl_matchLinkDistCluster (void *space, matchsplitsiteclusterlist_t *R, matchsplitsiteclusterlist_t *L, Uint n);
Uint bl_matchfileRealignScanFileNew(void *space, matchfile_t *file, FILE *realigndev, fasta_t *set, unsigned char fields, matchsplitsiteclusterlist_t **Lclust, matchsplitsiteclusterlist_t **Rclust, Uint *nchr, int threadno, int maxdist);
void bl_matchDestructMatchsplitsiteclusterlist (void *space, matchsplitsiteclusterlist_t *l, Uint n);
void *threadrealign(int begin, int stop, matchsplitsitecluster_t *T, Uint ovhglen, fasta_t *set, Alignment*** aligns,Uint chromidx,Uint seqlen, char left, matchfileRec_t r,void *space, int *e2, char ***refseqs, Uint **reflens, Uint **refstrand, Uint locreflen, Uint median, Uint interval, Uint start,  Uint end, char *fwd, char *rev, Uint sinterval, int oven);
void *realignthreadstarter(void *args);
void Readrealign(int begin, int stop, matchlistelem_t* arr, unsigned char fmt,  matchsplitsitecluster_t *T, char left, fasta_t *set, Uint chromidx,Uint *noofrealigns/*??*/,FILE *realigndev,  List *rlist, Uint interval,void *space, Uint sinterval);
void *Readthreadstarter(void *args);

#endif

