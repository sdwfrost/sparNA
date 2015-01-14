#ifndef KDMATCH_H
#define KDMATCH_H

/*
 *
 *	relaxed.h
 *  Declarations for kdmatch
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 12/26/07 02:41:35 CST  
 *
 *  SVN
 *  Revision of last commit: $Rev: 54 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-09-10 22:13:30 +0200 (Wed, 10 Sep 2008) $
 *
 *  Id: $Id: kdmatch.h 54 2008-09-10 20:13:30Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/trunk/src/kdmatch.h $
 *  
 */

#include "charsequence.h"
#include "kdseed.h"
#include "segemehl.h"
#include "karlin.h"
#include "bitVector.h"

#define SCORE(M, X, I, D)       (M) - (X) - (I) - (D)
#define MATCHES(S, X, I, D)     (S) + (X) + (I) + (D)
#define LEN_Q(S, X, I, D)       (S) + (X) + (I) + (D) + (X) + (I)  
#define LEN_S(S, X, I, D)       (S) + (X) + (I) + (D) + (X) + (D)


typedef struct bestseed_s {
  Uint readstart;
  Uint mat;
  Uint mis;
  Uint ins;
  Uint del;
  Uint len;
  Uint refidx;  
  Uint refpos;
  char refstrand;
  char maxintervalconstraint;
  char maxevalconstraint;

} bestseed_t;

typedef struct split_s {
  Uint subidx;
  char strand;
  Uint start;
  Uint end;
  uint16_t i;
  uint16_t j;
} split_t;

typedef struct spliceevent_s {
  uint8_t noofsites;
  Uint firstreadid;
  char *strand;
  Uint *subidx;
  Uint *start;
  Uint *end;
  uint16_t *i;
  uint16_t *j;
} spliceevent_t;

typedef struct spliceevents_s {
  Uint noofevents;
  spliceevent_t *event;
} spliceevents_t;


typedef struct spliceventmapelem_s{
  unsigned char type; 
  uint8_t site;
  spliceevent_t *ptr;
} spliceeventmapelem_t;

typedef struct spliceeventmap_s{
  Uint size;
  spliceeventmapelem_t *map;
}spliceeventmap_t;

void
se_kdGenomeMatch(void *space, Suffixarray *s, fasta_t *reads, 
    segemehl_t *nfo);

gmatchlist_t*
se_kdMatchStemAlign(void *space, 
    Suffixarray *s,
    MultiCharSeq *seq,
    matchstem_t **stems,
    char **sequences,
    Uint len,
    karlin_t *stats,
    segemehl_t *nfo, 
    Uint *enctab,
    bitvector* D, bestseed_t *best);

  void
bl_kdReportUnmatched (void *space, fasta_t *reads, Uint k, unsigned char matchflag,
    unsigned char matematchflag, bestseed_t* , bestseed_t *, segemehl_t *nfo);

#endif
