#ifndef MANOUT_H
#define MANOUT_H

/*
 * manout.h
 * attempt for flexible output of genome mapping w/ SEGEMEHL
 *
 * @author Christian Otto
 * @email christan@bioinf.uni-leipzig.de
 * @date Wed Sep 24 10:56:23 CEST 2008
 *
 */

#include <pthread.h>
#include "basic-types.h"
#include "bitArray.h"
#include "biofiles.h"
#include "charsequence.h"
#include "multicharseq.h"
#include "alignment.h"
#include "kdseed.h"
#include "segemehl.h"

#define MINUSSTRAND 0
#define PLUSSTRAND 1
#define SPLIT_NEXT_PLUS    ((unsigned char) (1 << 5))
#define SPLIT_PREV_PLUS    ((unsigned char) (1 << 6))


typedef enum matchstatus_e {
  QUERY, 
  MATE, 
  PAIR, 
  PAIR_REV, 
  PAIR_INS, 
  QUERY_SPL_NO_MATE, 
  QUERY_SPL_FULL_MATE,
  MATE_SPL_NO_QUERY,
  MATE_SPL_FULL_QUERY,
  PAIR_SPL
} matchstatus_t;

typedef struct gmate_s {

  unsigned char isset;
  Uint p;
  Uint q;
  int scr;
  int mat;
  int mis;
  int ins;
  int del;
  int edist;
  char *materefdesc;
  Uint materefdesclen;
  char *materefseq;

  Alignment *al;
  Uint subject;
 
} gmate_t;

typedef struct gmatch_s{
  Uint subject;
  unsigned char rc;
  Uint i;
  Uint j;
  Uint p;
  Uint q;
  int scr; 
  int mat;
  int mis;
  int ins;
  int del;
  int edist;
  Alignment *al;
  double evalue;

  Uint noofmatematches;
  Uint mateminedist;
  gmate_t mates[4];

  Uint fragno;
  Uint previdx;
  Uint prevpos;
  char prevflags;
  Uint nextidx;
  Uint nextpos;
  char nextflags;

  Uint prevseqstart;
  char *prevseqrefdesc;
  Uint nextseqstart;
  char *nextseqrefdesc;
  char *refdesc;
  Uint refdesclen;
  char *refseq;


  unsigned char skip;
} gmatch_t;


typedef struct gmatchlist_s{

  Uint minedist;
  Uint mateminedist;
  Uint pairminedist;

  Uint *n;
  gmatch_t **matches;

} gmatchlist_t;

typedef struct gread_s{
  Uint id;
  Uint noofmatepairs;
  Uint noofmatches;

  Uint n[2];
  gmatch_t* matches[2];

} gread_t;


typedef struct Gmap_s{

  MultiCharSeq *mseq;
  Uint mapoffset;
  Uint noofreads;
  gread_t *reads;

} Gmap;

void
se_destructMatches(void *space, gread_t *read);

unsigned char
se_hasMatches(gmatchlist_t *list);

unsigned char
se_hasMateMatches(gmatchlist_t *list);

Uint
se_kdSetMate(void *space, gmatch_t *match, 
    Uint chr_idx, Uint chr_start, Uint chr_end, Uint edist,
    Alignment *al, unsigned char downstream, unsigned char rc);

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
    Uint nextidx, Uint nextpos, char nextstrand, Uint fragno);

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
    double evalue, Alignment *al, Uint u, Uint n);
  
void reportSplicedMatch(void *space, char *qrydesc, 
    MultiCharSeqAlignment *mcsa, Uint noofaligns,
    Uint coverage, Uint edist,  int score, segemehl_t *nfo);
Uint se_kdMatchListLength(gmatchlist_t *list, unsigned char strand);
unsigned char se_kdMatchListhasMatches(gmatchlist_t *list);
unsigned char se_kdMatchListhasMates(gmatchlist_t *list);
Uint se_kdMatchListLength(gmatchlist_t *list, unsigned char strand);
Uint se_kdMatchListScore(gmatchlist_t *list);
gmatch_t* se_kdMatchListGet(gmatchlist_t *list, unsigned char strand, 
    Uint elem);
gmate_t* se_kdMatchGetMates(gmatch_t *match);
Uint se_kdMatchGetSubject(gmatch_t *match);
Uint se_kdMatchGetRefStart(gmatch_t *match);
extern void reportMap(FILE*, Gmap *map, Uint level);
extern void initMatch(gmatch_t *);
void initRead(gread_t *, Uint);
void initGmap(Gmap *, MultiCharSeq *, Uint);
extern void setMatches(gread_t*, gmatch_t *, Uint, unsigned char, 
    unsigned char);
extern void setReads(Gmap *, gread_t *, Uint);
extern Uint reportMatch (void *, Gmap *, fasta_t *, segemehl_t *, 
    matchstatus_t pairStatus, unsigned char mate);
Uint se_setMatches(void *space, gread_t *read, gmatchlist_t *list, Uint maxedist, segemehl_t *nfo, char rep);
void matchHeader(FILE* dev, Uint level);
void genericOutput (FILE *dev, char **list, Uint rep_type, char);
void bl_gmatchlistDestruct(void *space, gmatchlist_t *list);
gmatchlist_t* bl_gmatchlistInit(void *space, int maxedist, int matemaxedist);
void se_registerOutputDevice(void *space, segemehl_t *info);
bl_fileBins_t* se_createChromBins (void *space, fasta_t *f, int maxbins, char 
    *template, Uint tmplen);
bl_fileBinDomains_t* se_createChromDomains (void *space, fasta_t *f, 
    Uint minbins, Uint maxbins, char *filetemplate, Uint tmplen);
bl_fileBinDomains_t*
se_createBisulfiteBins (void *space, Uint noofdomains, Uint threadno, char *filetemplate, Uint tmplen);
char* se_SAMHeader (void *space, char **seq, Uint *seqlen, 
    Uint size, char *cmdline, char sep, char lf,
    unsigned char sorted);
char * se_defaultHeader (void *space, segemehl_t *info, char, char);
void se_storeHeader(void *space, char *filename, char **header, Uint *headerlen);

#endif /* MANOUT_H */

