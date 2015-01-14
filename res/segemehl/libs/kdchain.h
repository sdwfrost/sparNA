#ifndef KDCHAIN_H
#define KDCHAIN_M

/*
 *
 *	kdchain.h
 *  declarations and marcors for kd chaining
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 03/11/2008 06:40:32 PM CET  
 *
 *  SVN
 *  Revision of last commit: $Rev: 72 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-10-28 18:14:42 +0100 (Tue, 28 Oct 2008) $
 *
 *  Id: $Id: kdchain.h 72 2008-10-28 17:14:42Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/kdchain.h $
 */


#include "manout.h"
#include "container.h"
#include "sufarray.h"
#include "kdseed.h"

#define FSTART_S(a) ((a)->p)
#define FEND_S(a) (int)((a)->p+(a)->mat+(a)->mis+(a)->del-1)
#define RINTERVALSIZE 100
#define FSTART_Q(a) ((a)->i)
#define FEND_Q(a) (int)((a)->i+(a)->mat+(a)->mis+(a)->ins-1)
#define FASSIGN(d, s)        d->i   = s->i;\
                            d->p   = s->p;\
                            d->scr = s->scr;\
                            d->mat = s->mat;\
                            d->mis = s->mis;\
                            d->ins = s->ins;\
                            d->del = s->del

typedef struct {
  Uint size;
  Uint bestscr;
  Uint bestscrpos;
  gmatch_t* chains;
} gchains_t;

typedef struct {
  Uint start;
  Uint end;
  Uint x;
  unsigned char strand;
  Uint branchno;
  Uint substart;
  int score;
  double evalue;
  branch_t *branch;
} branchfragment_t;

typedef struct {
  Uint nooffragments;
  Uint start;
  Uint end;
  int score;
  branchfragment_t **f;
} branchChain_t;

void joinFragments(gmatch_t*, gmatch_t*, gmatch_t*, int err);
gmatch_t* greedyfchain(void *space, gmatch_t* F, Uint n, Uint *scr, 
    Uint *pos, Uint *m);
void wrapChains(void *space, branchChain_t *chains, Uint noofchains);
branchChain_t * branchChain(void *space, Suffixarray *arr, matchstem_t **stems, 
    char **seqs, Uint len, karlin_t *stats, Uint *noofchains, 
    branchfragment_t **, double maxevalue);
void showChains(branchChain_t *chains, Uint noofchains, Suffixarray *arr, 
    FILE *dev, char *, Uint);
int cmp_chainscores(const void *a, const void *b);
extern  void reportfchain(void *space, gchains_t *pi,  gmatch_t *e);
void branch2match(Suffixarray *s, Container *C, branch_t* b, 
    Uint noofbranches);
extern  void reportchaincoords(gmatch_t *dest, gmatch_t *a, gmatch_t *b,
    int ovq, int dmis, int dins, int ddel);
Container* findfchains(void *space, Suffixarray *s, matchstem_t* M,
    Uint m, Uint t, unsigned char strict, int sigmatch, double lambda,
    double H, double K, double maxevalue);
extern int misfrag (gmatch_t *a, gmatch_t *b);
extern int delfrag (gmatch_t *a, gmatch_t *b);
extern int insfrag (gmatch_t *a, gmatch_t *b);
extern void joinfrag(gmatch_t *dest, gmatch_t *a, gmatch_t *b);
int chainscore(branchChain_t *chain, branchfragment_t *f);
Lint chainovl(branchChain_t *chain, branchfragment_t *f);
void chain(branchChain_t *chain, branchfragment_t *f);
branchChain_t* condenseChain (branchChain_t * chains, Uint noofchains, MultiCharSeq *seq, Suffixarray *arr);

#endif
