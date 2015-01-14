#ifndef KDSEEDS_H
#define KDSEEDS_H

/*
 *
 *	kdseed.h
 *  gettin k-diff seeds
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 04/13/2008 12:05:48 AM CEST  
 *
 *  SVN
 *  Revision of last commit: $Rev: 77 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-11-17 13:16:59 +0100 (Mon, 17 Nov 2008) $
 *
 *  Id: $Id: kdseed.h 77 2008-11-17 12:16:59Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/kdseed.h $
 */

#include "basic-types.h"
#include "vstack.h"
#include "sufarray.h"
#include "karlin.h"

#define KMSTOREBRANCH(B,P)\
                    B[P].mat = m-data.mis;\
                    B[P].q   = m-data.mis;\
                    B[P].mis = data.mis;\
                    B[P].ins = 0;\
                    B[P].del = 0;\
                    B[P].l = data.l;\
                    B[P].r = data.r;\
                    B[P].u = data.u;\
                    B[P].v = data.v


typedef struct branch_s {
    Uint mis;
    Uint mat;
    Uint q;
    Uint ins;
    Uint del;
    Uint l;
    Uint r;
    Uint u;
    Uint v;
} branch_t;

typedef struct matchstem_s {
  
  branch_t *branches;
  Uint noofbranches;
} matchstem_t;



typedef struct {
  Uint sptr;
  Uint qptr;
  Uint l;
  Uint kcnt;
  Uint mat;
  Uint mis;
  Uint ins;
  Uint del;
  PairUint child;
  PairUint parent;
} kdiffm_t;

typedef struct {
  Uint l;
  Uint r;
  Uint u;
  Uint v;
  Uint p;
  Uint mis;
} kmis_t;


extern int kdscore(branch_t *b);

extern  void
dumpkdseeds(Suffixarray *s, matchstem_t *M, Uint m, char strand, Uint T);

double
kd_getBranchEvalue(matchstem_t *stem, Uint i, Uint q, Uint m, Uint n,
    karlin_t *stats);

int
kd_getBranchScore(matchstem_t *M, Uint i, Uint q);

extern int
kd_matchlcp(void *space,
            Suffixarray *s, char *p, unsigned int m, 
            kdiffm_t *data, unsigned int cursufpos, 
            unsigned int lcplen, unsigned int seedlen, VStack *vstack,
            unsigned int maxdiff, int sext, int pmis, int xoff, matchstem_t *b);
extern void
pushkdbranches ( void *space,
		 Suffixarray *s,
		 VStack *vstack, 
		 char *p, 
		 Uint m,
		 Uint kp,
		 kdiffm_t *data, 
		 PairUint pr);

extern void
pushkdlcp ( void *space, 
	    VStack *vstack, 
	    char *p, 
	    Uint m, 
	    Uint kp,
	    kdiffm_t *data);

extern void
kdbest ( void *space,
    Suffixarray *s,
    char *seqs[],
    Uint m,
    Uint sext,
    Uint pmis,
    Uint xoff,
    Uint kp,
    matchstem_t *a[],
    matchstem_t *b0[]);

extern matchstem_t*
kdseeds ( void *space,
    Suffixarray *s,
    char *p,
    Uint m,
    Uint jump,
    Uint sext,
    Uint pmis,
    Uint xoff,
    Uint kp,
    matchstem_t *b0);

extern matchstem_t*
kd_match ( void *space, 
    Suffixarray *s,
    char *p,
    Uint m,
    Uint sext,
    Uint pmis,
    Uint xoff,
    Uint maxdiff,
    Uint iprime,
    Uint jprime,
    Uint kprime,
    Uint lprime,
    Uint sptr,
    Uint qptr );

extern void
matchstemModifyBranch(void *space, 
    matchstem_t *a, Uint k,
    Uint mat, Uint q, 
    Uint mis, Uint ins, Uint del, 
    Uint l, Uint r, Uint u, Uint v);

extern void
matchstemAddBranch(void *space, 
    matchstem_t *a, 
    Uint mat, Uint q, 
    Uint mis, Uint ins, Uint del, 
    Uint l, Uint r, Uint u, Uint v);

extern void 
kd_updatebranches(matchstem_t* b, kdiffm_t *data, unsigned int seedlen);

void matchstemDestruct(void *space, matchstem_t *M);

  branch_t* kmis (void *space, Suffixarray *s, char *P, Uint m, Uint k, Uint *noofmatches);
void kdcompare(TripleSint *a, branch_t *b, Uint m);

void
bl_kdMatchstemDestruct(void *space, matchstem_t* stem, Uint len);

#endif
