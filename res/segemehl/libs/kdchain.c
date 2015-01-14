
/*
 *  kdchain.c
 *  implementation of kdchain
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 04/29/2008 07:01:30 PM CEST
 *  
 *  SVN
 *  Revision of last commit: $Rev: 85 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-11-18 15:34:44 +0100 (Tue, 18 Nov 2008) $
 *
 *  Id: $Id: kdchain.c 85 2008-11-18 14:34:44Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/kdchain.c $
 *  
 */

#include "manout.h"
#include "kdchain.h"
#include "mathematics.h"
#include "sufarray.h"
#include "container.h"
#include "kdseed.h"
#include "debug.h"
#include "karlin.h"
#include "bitvectoralg.h"
#include "iupac.h"
#include <assert.h>
#include <limits.h>
#include <unistd.h>
#include <float.h>
#include <math.h>

int
comp_gpos(const void *a, const void *b) {
  gmatch_t *l;
  gmatch_t *r;

  l = (gmatch_t*) a;
  r = (gmatch_t*) b;

  if (l->p < r->p) return -1;
  if (l->p > r->p) return 1;

  return 0;
}

inline void
reportchaincoords(gmatch_t *dest,
                  gmatch_t *a,
                  gmatch_t *b,
                  int ovq,
                  int dmis,
                  int dins,
                  int ddel) {
      
 DBG("a: qcoord=[%d,%d], scoord=[%d,%d] \nscr:%d mat:%d, mis:%d, ins:%d, del:%d\n", 
      FSTART_Q(a), FEND_Q(a), 
      FSTART_S(a), FEND_S(a),
      a->scr, a->mat, a->mis, a->ins, a->del); 
      
  DBG("b: qcoord=[%d,%d], scoord=[%d,%d] \nscr:%d mat:%d, mis:%d, ins:%d, del:%d\n", 
      FSTART_Q(b), FEND_Q(b), 
      FSTART_S(b), FEND_S(b),
      b->scr, b->mat, b->mis, b->ins, b->del); 
  
  DBG("dest: qcoord=[%d,%d], scoord=[%d,%d] \nscr:%d mat:%d, mis:%d, ins:%d, del:%d\n", 
      FSTART_Q(dest), FEND_Q(dest), 
      FSTART_S(dest), FEND_S(dest),
      dest->scr, dest->mat, dest->mis, dest->ins, dest->del); 
      
  DBG ("ovl:%d dmis:%d, dins:%d, ddel:%d\n", ovq, dmis, dins, ddel);

  return;
}

void
initgchains(gchains_t *pi) {
    pi->size = 0;
    pi->bestscr = 0;
    pi->bestscrpos=0;
    pi->chains = NULL;
}

inline int
misfrag (gmatch_t *a, gmatch_t *b) {
    int dq,
        dqa,
        dsa,
        ds,
        ddq,
        dds;
  
  dq = FSTART_Q(b) - FSTART_Q(a);
  ds = FSTART_S(b) - FSTART_S(a);

  dqa = FEND_Q(a) - FSTART_Q(a);
  dsa = FEND_S(a) - FSTART_S(a);
  
  ddq = dq - dqa;
  dds = ds - dsa;

  if (dds == ddq)  return dds;
  if (dds > ddq)   return ddq;
  if (ddq > dds)   return dds;

  return 0;
}


inline int
delfrag (gmatch_t *a, gmatch_t *b) {
    int dq,
        dqa,
        dsa,
        ds,
        ddq,
        dds;
  
  dq = FSTART_Q(b) - FSTART_Q(a);
  ds = FSTART_S(b) - FSTART_S(a);

  dqa = FEND_Q(a) - FSTART_Q(a);
  dsa = FEND_S(a) - FSTART_S(a);
  
  ddq = dq - dqa;
  dds = ds - dsa;

  if (dds > ddq)  return dds - ddq;

  return 0;
}


inline int
insfrag (gmatch_t *a, gmatch_t *b) {
    int dq,
        dqa,
        dsa,
        ds,
        ddq,
        dds;
  
  dq = FSTART_Q(b) - FSTART_Q(a);
  ds = FSTART_S(b) - FSTART_S(a);

  dqa = FEND_Q(a) - FSTART_Q(a);
  dsa = FEND_S(a) - FSTART_S(a);
  
  ddq = dq - dqa;
  dds = ds - dsa;

  if (ddq > dds)  return ddq - dds;

  return 0;
}

inline void
joinfrag(gmatch_t *dest, 
    gmatch_t *a, 
    gmatch_t *b) {

  int dmis=0,
      dins=0,
      ddel=0,
       ovq=0;
  

  if (a->scr > 0 && b->scr > 0) { 
    ovq  = abs(MIN(1, misfrag(a,b))-1);
    dmis = MAX(1, misfrag(a,b))-1;
    dins = insfrag(a,b);
    ddel = delfrag(a,b);
  } else {
    initMatch(dest); 
    return;
  }

  dest->scr = a->scr + b->scr - ovq - dmis - dins -ddel;
  dest->mat = a->mat + b->mat - ovq;
  dest->mis = a->mis + b->mis + dmis;
  dest->ins = a->ins + b->ins + dins;
  dest->del = a->del + b->del + ddel;

  dest->i   = a->i;
  dest->p   = a->p;
  
  if (dest->scr > 0 && FEND_Q(b) != FEND_Q(dest)) { 
    reportchaincoords(dest,a,b,ovq,dmis,dins,ddel);
  }
  
  if (dest->scr > 0 && FEND_S(b) != FEND_S(dest)) { 
    reportchaincoords(dest,a,b,ovq,dmis,dins,ddel);
  }

  dest->j   = FEND_Q(b);
  dest->q   = FEND_S(b);
  dest->subject = a->subject;
  return;
}

inline void
reportfchain(void *space, 
    gchains_t *pi,  
    gmatch_t *e) {

  pi->chains = ALLOCMEMORY(space, pi->chains, gmatch_t, pi->size+1);
  memmove(&pi->chains[pi->size], e, sizeof(gmatch_t));
  if (e->scr > pi->bestscr) {
    pi->bestscr = e->scr;
    pi->bestscrpos = pi->size;
  }
  pi->size++;
  return;
}



/*----------------------------- minDistfragment ------------------------------
 *    
 * @brief minimum distance of fragment hits
 * @author Steve Hoffmann 
 *   
 */
 
Uint
minDistFragmentHits (Suffixarray *arr, branchfragment_t *u, branchfragment_t *v)
{

  Uint i, j, d1idx, d2idx;
  Uint mindist = UINT_MAX;
  Uint d = UINT_MAX;

  for(i=u->branch->l; i <= u->branch->r; i++) {
    for(j=v->branch->l; j <= v->branch->r; j++) {
      d1idx = getMultiCharSeqIndex(arr->seq, &arr->seq->sequences[arr->suftab[i]]);
      d2idx = getMultiCharSeqIndex(arr->seq, &arr->seq->sequences[arr->suftab[j]]);
      if(d1idx == d2idx && (d=llabs((Lint)arr->suftab[i] - arr->suftab[j])) < mindist) {
        mindist = d;
      } 
    }
  }

  return d;
}


/*---------------------------------- chain -----------------------------------
 *    
 * @brief add fragment to a chain
 * @author Steve Hoffmann 
 *   
 */
 
void
chain(branchChain_t *chain, branchfragment_t *f) {
  chain->score = chainscore(chain, f);
  chain->end = f->end;
  chain->nooffragments++;
  chain->f = ALLOCMEMORY(space, chain->f, branchfragment_t*, chain->nooffragments);
  chain->f[chain->nooffragments-1] = f;
}


/*------------------------------- fragmentovl --------------------------------
 *    
 * @brief fragment overlap
 * @author Steve Hoffmann 
 *   
 */
 
int
fragmentovl (branchfragment_t *f1, branchfragment_t *f2)
{
  return (f1->end > f2->start) ? f1->end - f1->start : 0;
}

/*--------------------------------- chainovl ---------------------------------
 *    
 * @brief overlap of a fragment with chain on the query!
 * @author Steve Hoffmann 
 *   
 */
 

Lint
chainovl(branchChain_t *chain, branchfragment_t *f) {
  return ((Lint)chain->end - (Lint)f->start)+1;
}


/*-------------------------------- chainscore --------------------------------
 *    
 * @brief get score of a chain when a fragment is added to it
 * @author Steve Hoffmann 
 *   
 */
 

int
chainscore (branchChain_t *chain, branchfragment_t *f) {
  Lint ovl = 0;
  int score = 0;

  ovl = chainovl(chain, f);
                    //v -- 0 or -ovl
  ovl = (ovl < 0) ? 0 : ovl;

  score = (chain->score + f->score) - ovl;
  return score;
}



/*------------------------------- chainscore2 --------------------------------
 *    
 * @brief get score of chain with two fragments chained
 * @author Steve Hoffmann 
 *   
 */
 
int
chainscore2 (branchChain_t *chain, branchfragment_t *f1, branchfragment_t *f2)
{
  Lint ovl =0;
  Lint ovl2=0;
  int score =0;

  ovl = chainovl(chain, f1);
  ovl = (ovl < 0) ? -ovl : ovl;

  score = (chain->score + f1->score) - ovl;
	
  ovl2 = (Lint)f1->end - (Lint)f2->start + 1;
  ovl2 = (ovl2 < 0) ? -ovl2 : ovl2;

  score += f2->score - ovl2; 
    
  return score;
}

/*----------------------------- cmp_chainscores ------------------------------
 *    
 * @brief compare the scores of a chain (qsort)
 * @author Steve Hoffmann 
 *   
 */
 
int
cmp_chainscores (const void *a, const void *b) {
  branchChain_t *first = (branchChain_t *) a;
  branchChain_t *second = (branchChain_t *) b;

  if(first->score < second->score) return 1;
  if(first->score == second->score) return 0;

  return -1;
}


/*--------------------------- cmp_branchfragments ----------------------------
 *    
 * @brief compare start positions of branch fragmentsi (qsort)
 * @author Steve Hoffmann 
 *   
 */

int
cmp_branchfragments (const void *a, const void *b) {
  branchfragment_t *first = (branchfragment_t*) a;
  branchfragment_t *second = (branchfragment_t*) b;

  if(first->start < second->start) return -1;
  if(first->start == second->start) return 0;

  return 1;
}



/*------------------------------ condenseChain -------------------------------
 *    
 * @brief merge fragments of the chain if they are too close on the reference 
 * practically [u.v][x,y] -> [u,y]; 
 * if(u<x &% u<y) 
 *  if(x<v || x-v < 50) 
 *      merge 
 *  else
 *      dont merge
 * @author Steve Hoffmann 
 *   
 */
 
branchChain_t*
condenseChain (branchChain_t * chains, Uint noofchains, MultiCharSeq *seq, 
    Suffixarray *arr)
{
  Uint i, j, u, v, x, y, k, h, w, len1, len2, strand1, strand2, 
       chr1, chr2, nochain, d1idx=0, d2idx=0, d3idx=0, l, r, ll, rr, p=0, q=0;
  double  mindist = DBL_MAX, d1=0, d2=0, di, dj;
  branchChain_t *newchains = NULL;
  Uint sub_start, sub_end, subidx, **bd, beststart=0;

  l = chains->f[0]->branch->l;
  r = chains->f[0]->branch->r;

  bd = ALLOCMEMORY(space, NULL, Uint*, r-l+1);

  /******
   * first step: minimize distance of fragment hits within chain
   * for all possible start loci in [l,r] select
   * a chain of closest loci
   ******/

  for(i=0, p=l; p <= r; p++, i++) {
    bd[i] = calloc(chains->nooffragments, sizeof(Uint));
    bd[i][0] = p;

    for(di=0, j=1; j < chains->nooffragments; j++) {
      ll = chains->f[j]->branch->l;
      rr = chains->f[j]->branch->r;

      for(dj=0, q=ll; q <= rr; q++) {
        d1idx = getMultiCharSeqIndex(seq, 
            &seq->sequences[arr->suftab[q]]);
        d2idx = getMultiCharSeqIndex(seq, 
            &seq->sequences[arr->suftab[bd[i][j-1]]]);

        if(d1idx != d2idx) {
          d1 = UINT_MAX;
        } else {
          d1 = llabs((LLint) arr->suftab[bd[i][j-1]] - arr->suftab[q]);
        }

        if(bd[i][j]) {
          d3idx = getMultiCharSeqIndex(seq, 
              &seq->sequences[arr->suftab[bd[i][j]]]);
          d2 = llabs((LLint) arr->suftab[bd[i][j-1]] - arr->suftab[bd[i][j]]);
        }

        if(d3idx != d2idx) {
          d2 = UINT_MAX;
        }
        
        if(!bd[i][j] || d1 < d2) {
          bd[i][j] = q;
          dj = d1;
        }       
      }
      di += dj;
    }

    if(di < mindist) {
      beststart = i;
      mindist = di;
    }
  }

  /*decode substarts to real coordinates*/
  for(i=0; i < noofchains; i++) {    
    for(j=0; j < chains[i].nooffragments; j++) {
      chains[i].f[j]->substart = arr->suftab[bd[beststart][j]];
    }  
  }

  /******
   *second step: merge fragments that are close
   ******/

  for(i=0; i < noofchains; i++) {    
    for(j=0; j < chains[i].nooffragments-1; j++) {
      for(k=j+1; k < chains[i].nooffragments; k++) { 
  //   for(k=j+1; k < j+2; k++) { 

        len1 = chains[i].f[j]->end - chains[i].f[j]->start;
        strand1 = chains[i].f[j]->strand;
        h = chains[i].f[j]->end;
        u = arr->suftab[bd[beststart][j]];
        v = u + len1;
        chr1 = getMultiCharSeqIndex(arr->seq, &arr->seq->sequences[u]);

        len2 = chains[i].f[k]->end - chains[i].f[k]->start;
        strand2 = chains[i].f[k]->strand;
        w = chains[i].f[k]->start;
        x = arr->suftab[bd[beststart][k]];
        y = x + len2;
        chr2 = getMultiCharSeqIndex(arr->seq, &arr->seq->sequences[x]);

        // read:           h)--(w
        // reference:   [u,v]  [x,y]
        // if the distance of h and w equals the distance of v and x we assume
        // that the fragments need to be merged
        if(u < x && v < y && chr1 == chr2 && strand1 == strand2 && strand1 == 0) {
          if(x < v || (x-v <= 20 && w-h <= 20) || (w >= h &&  x-v < w-h+20 && x-v+20 > w-h)) { 
            //merge
#ifdef DEBUGTRANSALIGN
            fprintf(stdout, "merging fragment %d with %d. [%d,%d;%d,%d] with [%d,%d;%d,%d]. u=%u : %u=y, u-y=%u \t h=%u : %u=w, w-h=%u\n",
                j, k,
                chains[i].f[j]->start, chains[i].f[j]->end,
                u, v,
                chains[i].f[k]->start, chains[i].f[k]->end,
                x, y,
                u, y,
                u-y,
                h, w,
                w-h);
#endif
            chains[i].f[j]->start = MIN(chains[i].f[j]->start, chains[i].f[k]->start);
            chains[i].f[j]->end = MAX(chains[i].f[j]->end, chains[i].f[k]->end);
            chains[i].f[k]->start =  chains[i].f[j]->start;
            chains[i].f[k]->end = chains[i].f[j]->end;
            bd[beststart][j] = (u < x) ? bd[beststart][j] : bd[beststart][k];
            bd[beststart][k] = bd[beststart][j];
            chains[i].f[j]->substart = arr->suftab[bd[beststart][j]];
            chains[i].f[k]->substart = arr->suftab[bd[beststart][j]];

          } else {
            //dont merge
          }
        }
        //   h)--(w
        //[x,y]  [u,v]
        if(x < u && y < v && chr1 == chr2 && strand1 == strand2 && strand1 == 1) {
          //if(u < y || u-y < ((w >= h) ? w-h+20 : 20)) {
          if(u < y || (u-y <= 20 && w-h <= 20) || (w >= h &&  u-y < w-h+20 && u-y+20 > w-h)) { 
            //merge       
#ifdef DEBUGTRANSALIGN
            fprintf(stdout, "merging fragment %d with %d. [%d,%d;%d,%d] with [%d,%d;%d,%d]. u=%u : %u=y, u-y=%u \t h=%u : %u=w, w-h=%u\n",
                j, k,
                chains[i].f[j]->start, chains[i].f[j]->end,
                u, v,
                chains[i].f[k]->start, chains[i].f[k]->end,
                x, y,
                u, y,
                u-y,
                h, w,
                w-h);
#endif
            chains[i].f[j]->start = MIN(chains[i].f[j]->start, chains[i].f[k]->start);
            chains[i].f[j]->end = MAX(chains[i].f[j]->end, chains[i].f[k]->end);
            chains[i].f[k]->start =  chains[i].f[j]->start;
            chains[i].f[k]->end = chains[i].f[j]->end;
            bd[beststart][j] = (u < x) ? bd[beststart][j] : bd[beststart][k];
            bd[beststart][k] = bd[beststart][j];
            chains[i].f[j]->substart = arr->suftab[bd[beststart][j]];
            chains[i].f[k]->substart = arr->suftab[bd[beststart][j]];

          } else {
            //dont merge
#ifdef DEBUGTRANSALIGN
            fprintf(stdout, "\t not merging fragment %d with %d. [%d,%d;%d,%d] with [%d,%d;%d,%d]. u=%u : %u=y, u-y=%u \t h=%u : %u=w, h-w=%u\n",
                j, k,
                chains[i].f[j]->start, chains[i].f[j]->end,
                u, v,
                chains[i].f[k]->start, chains[i].f[k]->end,
                x, y,
                u, y,
                u-y,
                h, w,
                h-w);
#endif
          }
        }

      }
    }
  }

  newchains = ALLOCMEMORY(space, NULL, branchChain_t, noofchains);

  for(i=0; i < noofchains; i++) {
    newchains[i].nooffragments = 0;
    newchains[i].f = NULL;
    for(j=0; j < chains[i].nooffragments; j++) {
      nochain = 0;
      len1 = chains[i].f[j]->end - chains[i].f[j]->start;
      for (k=j+1; k < chains[i].nooffragments; k++) {
        len2 = chains[i].f[k]->end - chains[i].f[k]->start;
        if(chains[i].f[j]->start  >= chains[i].f[k]->start  && 
           chains[i].f[j]->end    <= chains[i].f[k]->end    && 
           chains[i].f[j]->strand == chains[i].f[k]->strand &&
           chains[i].f[j]->substart >= chains[i].f[k]->substart &&
           chains[i].f[j]->substart+len1 <= chains[i].f[k]->substart+len2) {
          nochain = 1;
        } 
      }

      if(nochain) {
      } else {
        subidx = getMultiCharSeqIndex(seq, &seq->sequences[chains[i].f[j]->substart]);
        getMultiCharSeqIdxBounds(seq, subidx, &sub_start, &sub_end);
#ifdef DEBUGTRANSALIGN
        fprintf(stdout, "adding element [%d,%d] to newchain -> %d (%u)\n", chains[i].f[j]->start, chains[i].f[j]->end, chains[i].f[j]->substart-sub_start, chains[i].f[j]->substart);
#endif
        chain(&newchains[i], chains[i].f[j]); 
      }
    }
  }
 
  for (i=0; i < r-l+1; i++) {
    FREEMEMORY(space, bd[i]);
  }
  FREEMEMORY(space, bd);

  return newchains;
}

/*------------------------------- branchChain --------------------------------
 *    
 * @brief find chain of branches
 * @author Steve Hoffmann 
 *   
 */
 
branchChain_t *
branchChain(void *space, Suffixarray *arr, matchstem_t **stems, char **seqs, 
    Uint len, karlin_t *stats, Uint *noofchains, branchfragment_t **fragments, double maxevalue) {
  
  Uint i, j, u, start, end, substart, 
    k=0, l, r, s, c_prime=0, x, maxocc = 50; //q_prime = 0, c=0, q, v;
  int maxovl = 12, bestscr; //maxgap = 22;
  double minentropy = 1.5, E, H;
//  unsigned char *select =NULL;
  branch_t *branch;
  branchChain_t *chains = NULL;
  branchfragment_t *f = NULL ;// *f_prime = NULL; 
  
  for (i = 0; i < len; i++) {
    for (u = 0; u < 2; u++) {
      x = (u == 0) ? i : len-1-i;

      for (s = 0; s < stems[u][x].noofbranches; s++) {     
        branch = &stems[u][x].branches[s];
//        for (v=branch->l; v <= branch->r; v++) { 
        l = branch->l;
        r = branch->r;

        if (u == 0) {
          start = i;
	      //CHANGED: end position one position too far
          //before: end = i + branch->mat;
          end = i + branch->mat - 1;
          substart = arr->suftab[l]; //l to v
        } else {
          start = i - branch->mat + 1;
          end = i;
          substart = arr->suftab[l] - branch->mat + 1; //l to v
        }

        E = kd_getBranchEvalue(stems[u], x, s, len, arr->numofsuffixes, stats);
        H = minshannonentropy(&seqs[0][start], end-start+1);


#ifdef DEBUGTRANSALIGN       
//      for(q=chains[i].f[j]->branch->l; q <= chains[i].f[j]->branch->r; q++) { 
        Uint sub_idx, sub_start =0, sub_end=0;
        sub_idx = getMultiCharSeqIndex(arr->seq, &arr->seq->sequences[arr->suftab[l]]);
        getMultiCharSeqIdxBounds(arr->seq, sub_idx, &sub_start, &sub_end);

//        fprintf(dev,"%u (chr:%d) -> %u, ",arr->suftab[l], 
//            getMultiCharSeqIndex(arr->seq, &arr->seq->sequences[arr->suftab[l]]),arr->suftab[l]-sub_start);
//      }

        fprintf(stdout, "%d-[%d,%d] -> chr:%d-%d\tx:%d\t", k, start, end, 
            getMultiCharSeqIndex(arr->seq, &arr->seq->sequences[arr->suftab[l]]),arr->suftab[l]-sub_start
            , x);
        fprintf(stdout, "Evalue: %f (max: %f), x:%u, s:%u, len:%u, H: %f (min %f), occ:%d, scr:%d\n", E, maxevalue, x, s, len, H, minentropy, r-l, kdscore(branch));
#endif
        if(E < maxevalue && H > minentropy && l <= r && r - l < maxocc) {
  //        fprintf(stdout, "accepted\n");
          for(j=0; j < k; j++) {
	        //CHANGED:
	        //from (Lint)arr->suftab[branch->l] to (Lint)substart
	        //in order to account for minus strand hits correctly
            //changed (Lint)f[j].start - x to (Lint)f[j].start - (Lint)start
/*
            fprintf(stdout," substart[%d]-substart[cur]:%ld < %ld:start[%d]-start[cur]; substart[%d]:%d, substart[cur]:%d, start[%d]:%d, start[cur]:%d; x:%d\n",
                j, labs((Lint)f[j].substart - (Lint)substart),  
                labs((Lint)f[j].start - x) + maxovl, j,
                j, f[j].substart, substart,
                j, f[j].start, start, x);
//            if( labs((Lint)f[j].substart - (Lint)substart) <
//                labs((Lint)f[j].start - x) + maxovl)
//              break;
*/
            Uint x1 = f[j].start;
            Uint y1 = f[j].end;
            Uint x2 = start;
            Uint y2 = end;

#ifdef DEBUGTRANSALIGN            
             fprintf(stdout, "[%d,%d] vs. [%d,%d]\n", x1, y1, x2, y2);
#endif
            //inclusion
            if(x1 <= x2 && y2 <= y1) { 
#ifdef DEBUGTRANSALIGN
              fprintf(stdout, "inclusion consumption (1)\n");
#endif
              break;
            }
            if(x2 <= x1 && y1 <= y2) { 
#ifdef DEBUGTRANSALIGN
              fprintf(stdout, "inclusion consumption (2)\n");
#endif
              break;
            }
            //overlap
            if(y1 > x2 && y2 > y1 && y1-x2 >= maxovl){ 
#ifdef DEBUGTRANSALIGN
              fprintf(stdout, "overlap consumption (1)\n");
#endif
    //          break;
            }
            if(y2 > x1 && y1 > y2 && y2-x1 >= maxovl) { 
#ifdef DEBUGTRANSALIGN
              fprintf(stdout, "overlap consumption (2)\n");
#endif
  //            break;
            }
          }


          if(j < k) {
            if (kdscore(branch) > f[j].score) { 
  //            fprintf(stdout, "replaced %d\n", k);
              f[j].start = start;  
              f[j].end = end;
              f[j].substart = substart;
              f[j].strand = (unsigned char) u;
              f[j].branchno = s;
              f[j].branch = branch;
              f[j].score = kdscore(branch);
              f[j].x = x;
              f[j].evalue = E;
            }
            continue;
          }

#ifdef DEBUGTRANSALIGN          
          fprintf(stdout, "adding %d [%d,%d]\n", k, start, end);
#endif 
          k++;
          f = ALLOCMEMORY(space, f, branchfragment_t, k);
          f[k-1].start = start ;
          f[k-1].end = end;
          f[k-1].substart = substart;
          f[k-1].strand = (unsigned char) u;
          f[k-1].branchno = s;
          f[k-1].branch = branch;
          f[k-1].score = kdscore(branch);
          f[k-1].x = x;
          f[k-1].evalue = E;


#ifdef FIXINSMALL

          Uint parent_sub_start = 0, parent_sub_end=0;
          Uint parent_sub_idx = getMultiCharSeqIndex(arr->seq, &arr->seq->sequences[arr->suftab[l]]);
          getMultiCharSeqIdxBounds(arr->seq, sub_idx, &parent_sub_start, &parent_sub_end);
          char *fixinseq;
          fixinseq = ALLOCMEMORY(space, NULL, char, 11);
          matchstem_t *b = NULL;
          fixinseq[10] =0;
          if(start > 10) {
            fprintf(stdout, "parent mapped to chr:%d-%d\n", parent_sub_idx, arr->suftab[l]-parent_sub_start);
            for(j=0; j <= start-10; j++) {
              memmove(fixinseq, &seqs[0][j], 10);
              fprintf(stdout, "fixin start %d-%d: %s\n", j, j+10-1, fixinseq);
              b = kd_match(space, arr, fixinseq, 10, 0, 0, 0, 0, 0, arr->numofsuffixes-1, 0, arr->numofsuffixes-1, 0, 0); 
              Uint width = b->branches[0].r - b->branches[0].l;
              //fprintf(stdout, "found interval %d-%d (%d): %d\n", b->branches[0].l, b->branches[0].r, width, b->branches[0].mat);
              for(Uint wi=0; wi < width; wi++) {
                //simple check for distance to save computation of coords
                if(dist_uint(arr->suftab[l], arr->suftab[wi]) < 200000) { 
                Uint sub_idx, sub_start =0, sub_end=0;
                sub_idx = getMultiCharSeqIndex(arr->seq, &arr->seq->sequences[arr->suftab[wi]]);
                getMultiCharSeqIdxBounds(arr->seq, sub_idx, &sub_start, &sub_end);
                if(sub_idx == parent_sub_idx) { 
                fprintf(stdout, "-> chr:%d-%d\n", getMultiCharSeqIndex(arr->seq, &arr->seq->sequences[arr->suftab[wi]]),arr->suftab[wi]-sub_start);
                }
                }
              }
              FREEMEMORY(space, b->branches);
              FREEMEMORY(space, b);
            }
          }
          if((len > 10 && end <= len-10)) {
            for(j=end; j <= len-10; j++) { 
              memmove(fixinseq, &seqs[0][j], 10);
              fprintf(stdout, "fixin end %d-%d: %s\n", j, j+10-1, fixinseq);
            }
          }
          FREEMEMORY(space, fixinseq);
#endif
        } else {
          // fprintf(stdout, "rejected: too short, small evalue\n");
        }
        //      } //new
      }
    }
  }

  qsort(f, k, sizeof(branchfragment_t), cmp_branchfragments);
  chains = ALLOCMEMORY(space, chains, branchChain_t, k);

 // fprintf(stdout, "maxgap: %d, maxovl:%d\n", maxgap, maxovl);
  for(i = 0; i < k; i++) {

    c_prime = i;
    bestscr = 0;

    // search for best precedessor chain
    for (j = 0; j < i; j++){

      // only allow short overlap
      if (chainovl(&chains[j], &f[i]) < maxovl) {

        // update best precessor
        if (bestscr < chainscore(&chains[j], &f[i])){
          bestscr = chainscore(&chains[j], &f[i]);
          c_prime = j;
        }
        if (bestscr == chainscore(&chains[j], &f[i]) && (c_prime == i ||
            minDistFragmentHits(arr, chains[j].f[chains[j].nooffragments-1], &f[i]) <
            minDistFragmentHits(arr, chains[c_prime].f[chains[c_prime].nooffragments-1], &f[i]))) {          
          bestscr = chainscore(&chains[j], &f[i]);
          c_prime = j;
        }
      }
    }

    if (c_prime != i){
      // TODO: add function to do the following
      chains[i].nooffragments = chains[c_prime].nooffragments + 1;
      chains[i].f = ALLOCMEMORY(space, NULL, branchfragment_t*, chains[i].nooffragments);
      memmove(chains[i].f, chains[c_prime].f, sizeof(branchfragment_t*) * chains[c_prime].nooffragments);
      chains[i].f[chains[i].nooffragments-1] = &f[i];
      chains[i].score = bestscr;
      chains[i].end = f[i].end;
    }
    else {
      chains[i].nooffragments = 1;
      chains[i].f = 
        ALLOCMEMORY(space, NULL, branchfragment_t*, chains[i].nooffragments);
      chains[i].f[0] = &f[i];
      chains[i].score = f[i].score;
      chains[i].start = f[i].start;
      chains[i].end = f[i].end;
    }
  }


  (*noofchains) = k;
  (*fragments) = f;

  return chains;
}



/*-------------------------------- showChains --------------------------------
 *    
 * @brief dump the chains
 * @author Steve Hoffmann 
 *   
 */
 
void
showChains(branchChain_t *chains, Uint noofchains, Suffixarray *arr, 
    FILE *dev, char *seq, Uint len) {
  
  Uint i, j, q, subidx, sub_start, sub_end;
  double H;

  for(i=0; i < noofchains; i++) {
    fprintf(dev, "chain %d: %d-%d (%d)\n", i, chains[i].start, 
        chains[i].end, chains[i].score);
    
    for(j=0; j < chains[i].nooffragments; j++) {
           
      fprintf(dev, "fragment %d: %d-%d (%d) (%d:%f) substart:", j, 
          chains[i].f[j]->start, chains[i].f[j]->end, 
          chains[i].f[j]->strand, chains[i].f[j]->score, 
          chains[i].f[j]->evalue);
      
      for(q=chains[i].f[j]->branch->l; q <= chains[i].f[j]->branch->r; q++) {
        subidx = getMultiCharSeqIndex(arr->seq, &arr->seq->sequences[arr->suftab[q]]);
        getMultiCharSeqIdxBounds(arr->seq, subidx, &sub_start, &sub_end);

        fprintf(dev,"%u (chr:%d) -> %u, ",arr->suftab[q], 
            getMultiCharSeqIndex(arr->seq, &arr->seq->sequences[arr->suftab[q]]),arr->suftab[q]-sub_start);
      }
      //CHANGED: from end-start to end-start+1 for length
      //before: H = shannonentropy(NULL, &seq[chains[i].f[j]->start], chains[i].f[j]->end - chains[i].f[j]->start, asize, tab);
      H = minshannonentropy(&seq[chains[i].f[j]->start], chains[i].f[j]->end - chains[i].f[j]->start + 1);

      fprintf(dev, "entropy: %f\n", H);
    }
    fprintf(dev, "\n");
  }
}


/*-------------------------------- wrapChains --------------------------------
 *    
 * @brief remove chains from heap
 * @author Steve Hoffmann 
 *   
 */
 
void
wrapChains(void *space, branchChain_t *chains, Uint noofchains) {
  Uint i;

  for(i=0; i < noofchains; i++) { 
    FREEMEMORY(space, chains[i].f);
  }
  return;
}



gmatch_t*
greedyfchain( void *space,
    gmatch_t *F,
    Uint n,
    Uint *scr,
    Uint *pos,
    Uint *m) {


  int ovl;
  Uint i = 0;
  gmatch_t *h,
           *tmp1,
           *tmp2,
           *f,
           *e;
  gmatch_t *chains;
  gchains_t *pi=NULL;

  f = ALLOCMEMORY(space, NULL, gmatch_t, 1);
  e = ALLOCMEMORY(space, NULL, gmatch_t, 1);
  h = ALLOCMEMORY(space, NULL, gmatch_t, 1);
  tmp1 = ALLOCMEMORY(space, NULL, gmatch_t, 1);
  tmp2 = ALLOCMEMORY(space, NULL, gmatch_t, 1);
  pi = ALLOCMEMORY(space, NULL, gchains_t, 1);

  initgchains(pi);
  initMatch(h);
  memmove(e, &F[0], sizeof(gmatch_t));

  for(i=0; i < n; i++)  {
    memmove(f, &F[i], sizeof(gmatch_t));
    ovl = FEND_S(e) - FSTART_S(f);

    if (FEND_Q(f) > FEND_Q(e) && FSTART_Q(f) < FSTART_Q(e)) {
      memmove(e, f, sizeof(gmatch_t));
    /*case 1: better start*/
    } else if (FEND_S(f) == FEND_S(e) && f->scr >= e->scr) {
      memmove(e, f, sizeof(gmatch_t));
    }
    /*case 2: positive overlap*/
    else if (FEND_S(f) > FEND_S(e) && ovl > 0) {
      joinfrag(tmp1, e, f);
      if(tmp1->scr > h->scr) {
        memmove(h,tmp1,sizeof(gmatch_t));
      } 
    }
    /*case 3: negative/zero overlap*/
    else if (FEND_S(f) > FEND_S(e) && ovl <= 0) {
      joinfrag(tmp1, h,f);
      joinfrag(tmp2, e,f);
      /*case 3a: bridge has better scr*/
      if (h->scr > 0 && FEND_S(f) > FEND_S(h) && tmp1->scr > tmp2->scr) {
        /*case 3aa: bridge of bridges*/
        if (FSTART_S(f) < FEND_S(h)) {
          memmove(e, h, sizeof(gmatch_t));
          memmove(h, tmp1, sizeof(gmatch_t));
        } else {
          memmove(e, tmp1, sizeof(gmatch_t));
          initMatch(h);
        }
      /*case 3b: accept d(e,f)*/
      } else if (tmp2->scr > e->scr) {
        memmove(e, tmp2, sizeof(gmatch_t));
        initMatch(h);
      /*case 3c: chain terminated*/
      } else {
        reportfchain(space, pi, e);
        memmove(e, f, sizeof(gmatch_t));
        initMatch(h);
      }
    } 
  }

  if(h->scr > e->scr) {
    reportfchain(space, pi, h);
  } else {
    reportfchain(space, pi, e);
  }

  FREEMEMORY(space, f);
  FREEMEMORY(space, e);
  FREEMEMORY(space, h);
  FREEMEMORY(space, tmp1);
  FREEMEMORY(space, tmp2);

  (*m) = pi->size;
  (*scr) = pi->bestscr;
  (*pos) = pi->bestscrpos;
  chains = pi->chains;

  FREEMEMORY(space, pi);
  
  return chains;
}

void
branch2match(Suffixarray *s,
    Container *C, 
    branch_t* b, 
    Uint noofbranches) {
  Uint i,
    k,
    idx;
  gmatch_t m;
  
  assert(noofbranches > 0);

  for(i=0; i < noofbranches; i++) {
    for(k=b[i].l; k <= b[i].r; k++) {
      m.p = s->suftab[k];
      m.scr = b[i].mat - b[i].mis - b[i].ins - b[i].del;
      m.mat = b[i].mat;
      m.del = b[i].del;
      m.ins = b[i].ins;
      m.mis = b[i].mis;
      m.i = i;
      m.j = i+b[i].mat+b[i].mis+b[i].ins-1;
      m.q = m.p+b[i].mat+b[i].mis+b[i].del-1;
      m.evalue = 0;
      m.scr = b[i].mat - b[i].mis - b[i].ins - b[i].del;
      idx = getMultiCharSeqIndex(s->seq, s->seq->sequences + m.p);
      m.subject = idx;
      if (b[i].mat > b[i].mis + b[i].ins + b[i].del) {
        bl_containerAdd(C, &m);
      }
    }
  }
  qsort(C->contspace, bl_containerSize(C), sizeof(gmatch_t), comp_gpos);
}

Container*
findfchains(void *space,
    Suffixarray *s,
    matchstem_t* M,
    Uint m,
    Uint t,
    unsigned char strict,
    int sigmatch,
    double lambda,
    double H,
    double K,
    double maxevalue) {

  Uint i,
    j,
    k,
    start,
    range,
    size,
    bestscr,
    bestscrpos,
    chainno,
    idx;

  Container P, 
              *C;
  gmatch_t p,
           *tmp,
           *ptr,
           *last,
           *G; 

  C = (Container *) malloc(sizeof(Container));
  bl_containerInit(C, 1000, sizeof(gmatch_t));
  bl_containerInit(&P, 1000, sizeof(gmatch_t));

  for(i=0; i < m; i++) {
    for(j=0; j < M[i].noofbranches; j++) {
      if(M[i].branches[j].r - M[i].branches[j].l < t) {
        for(k=M[i].branches[j].l; k <= M[i].branches[j].r; k++) {
          p.p = s->suftab[k];
          p.scr = M[i].branches[j].mat - M[i].branches[j].mis - M[i].branches[j].ins - M[i].branches[j].del;
          p.mat = M[i].branches[j].mat;
          p.del = M[i].branches[j].del;
          p.ins = M[i].branches[j].ins;
          p.mis = M[i].branches[j].mis;
          p.i   = i;
          p.j   = p.i+M[i].branches[j].mat+M[i].branches[j].mis+M[i].branches[j].ins-1;
          p.q   = p.p+M[i].branches[j].mat+M[i].branches[j].mis+M[i].branches[j].del-1;
          idx = getMultiCharSeqIndex(s->seq, s->seq->sequences + p.p);
          p.subject = idx;
          if(M[i].branches[j].mat > M[i].branches[j].mis + M[i].branches[j].ins + M[i].branches[j].del) {
            bl_containerAdd(&P, &p);
          }
        }
      }
    }
  }

  if(bl_containerSize(&P)) {
    qsort(P.contspace, bl_containerSize(&P), sizeof(gmatch_t), comp_gpos);

    ptr = (gmatch_t*)bl_containerGet(&P,0);
    last = ptr;
    start = 0;

    for(i=1; i < bl_containerSize(&P); i++) {
      ptr = (gmatch_t*) bl_containerGet(&P,i);

      if ((ptr->p - last->p) > last->scr+RINTERVALSIZE) {
        tmp = bl_containerGet(&P,start);
        range = (last->p - tmp->p)+tmp->scr+1;
        size = i-1-start; 
        DBGL(5, "interval start %d-%d[%d,%d]\n", range, size, start, i-1);

        G = greedyfchain(space, tmp, size, &bestscr, &bestscrpos, &chainno);
        DBGL(5, "sigmatch: %d, bestscr:%d\n", sigmatch, G[bestscrpos].scr); 
        if (G[bestscrpos].scr > 0)
          G[bestscrpos].evalue = evalue(lambda, K, spacemult(m, s->numofsuffixes, H, K), G[bestscrpos].scr);

        if (strict) { 
          if ((G[bestscrpos].scr >= sigmatch 
                && G[bestscrpos].evalue < maxevalue)) { /*G[bestscr] >= 12*/
            bl_containerAdd(C, &G[bestscrpos]);
            DBGL(5, "added match to container %d\n", G[bestscrpos].scr);
          }
        } else {
          if ((G[bestscrpos].scr >= sigmatch 
                || G[bestscrpos].evalue < maxevalue)) {
            bl_containerAdd(C, &G[bestscrpos]);
            DBGL(5, "added match to container %d\n", G[bestscrpos].scr);
          }
        }

        bestscrpos = 0;
        bestscr = 0;
        start = i; 
        free(G);
      }

      DBGL(5, "pos:%d mat:%d mis:%d ins:%d del:%d\n", ptr->p, ptr->scr, 
          ptr->mis, ptr->ins, ptr->del); 
      last = ptr;
    }

    tmp = bl_containerGet(&P,start);
    range = (last->p - tmp->p)+tmp->scr+1;
    size = i-1-start; 
    DBGL(5, "interval start %d-%d[%d,%d]\n", range, size, start, i-1);

    G = greedyfchain(space, tmp, size, &bestscr, &bestscrpos, &chainno);
    DBGL(5, "sigmatch: %d, bestscr:%d\n", sigmatch, G[bestscrpos].scr); 
    if (G[bestscrpos].scr > 0)
      G[bestscrpos].evalue = evalue(lambda, K, spacemult(m, s->numofsuffixes, H, K), G[bestscrpos].scr);

    if (strict) { 
      if (G[bestscrpos].scr > 0 && (G[bestscrpos].scr >= sigmatch 
            && G[bestscrpos].evalue < maxevalue)) { /*G[bestscr] >= 12*/
        bl_containerAdd(C, &G[bestscrpos]);
        DBGL(5, "added match to container %d\n", G[bestscrpos].scr);
      }
    } else {
      if (G[bestscrpos].scr > 0 && (G[bestscrpos].scr >= sigmatch 
            || G[bestscrpos].evalue < maxevalue)) {
        bl_containerAdd(C, &G[bestscrpos]);
        DBGL(5, "added match to container %d\n", G[bestscrpos].scr);
      }
    }


    bestscrpos = 0;
    bestscr = 0;
    start = i; 
    FREEMEMORY(space, G);
  }
  bl_containerDestruct(&P, NULL);
  return C;
}

#ifdef NOCOMPILE
  char *fixinseq;
  fixinseq = ALLOCMEMORY(space, NULL, char, 11);
  fixinseq[10] =0;
  
  for(i=0; i < k; i++) {
    fprintf(stdout, "chain %d: %d-%d\n", i, chains[i].start, chains[i].end);
    if(chains[i].start > 10) {
      for(j=0; j <= chains[i].start-10; j++) {
        memmove(fixinseq, &seqs[0][j], 10);
        fprintf(stdout, "fixin start %d-%d: %s\n", j, j+10-1, fixinseq);
      }
    }
    if((len > 10 && chains[i].end <= len-10)) {
      for(j=chains[i].end; j <= len-10; j++) { 
        memmove(fixinseq, &seqs[0][j], 10);
        fprintf(stdout, "fixin end %d-%d: %s\n", j, j+10-1, fixinseq);
      }
    }
  }
  FREEMEMORY(space, fixinseq);
#endif

