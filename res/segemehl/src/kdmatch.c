
/*
 *  kdmatch.c
 *  routines 4 relaxed alignments
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 11/27/2007 04:08:39 PM CET
 *
 *  SVN
 *  Revision of last commit: $Rev: 103 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-12-10 15:18:18 +0100 (Wed, 10 Dec 2008) $
 *
 *  Id: $Id: kdmatch.c 103 2008-12-10 14:18:18Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/src/kdmatch.c $
 *  
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <unistd.h>
#include <float.h>
#include <math.h>
#include <assert.h>
#include "memory.h"
#include "fileio.h"
#include "stringutils.h"
#include "charsequence.h"
#include "multicharseq.h"
#include "sufarray.h"
#include "mmchar.h"
#include "mathematics.h"
#include "manout.h"
#include "biofiles.h"
#include "vtprogressbar.h"
#include "karlin.h"
#include "sort.h"
#include "basic-types.h"
#include "bitvectoralg.h"
#include "bitVector.h"
#include "kdmatch.h"
#include "bitArray.h"
#include "segemehl.h"
#include "container.h"
#include "kdchain.h"
#include "debug.h"
#include "info.h"
#include "kdseed.h"
#include "alignment.h"
#include "sw.h"
#include "seqclip.h"
#include <pthread.h>
#include "iupac.h"

/*---------------------------- se_kdFindBestMate -----------------------------
 *
 * @brief find the 'best' mate from a list of hits to a given hit
 * @author Christian Otto
 *
 */

gmatchlist_t *
se_kdFindBestMate (void *space, gmatchlist_t *list, gmatch_t *match, Uint maxedist){
  Uint i, u;
  int bestedist = -1;
  unsigned char found = 0;
  Alignment *al;
  PairUint best;
  gmatchlist_t *res;

  res = NULL;
  memset(&best, 0, sizeof(PairUint));

  for (u = 0; u < 2; u++){
    for (i = 0; i < list->n[u]; i++){
      if (list->matches[u][i].edist > maxedist) 
        continue;

      if (list->matches[u][i].subject == match->subject){        
	if (!found || abs((LLint)list->matches[best.a][best.b].p - match->p) > 
	    abs((LLint)list->matches[u][i].p - match->p)){
	  found = 1;
	  best.a = u;
	  best.b = i;
          bestedist = list->matches[u][i].edist;
	}
      }
      else {
        if (!found && (bestedist == -1 || bestedist > list->matches[u][i].edist)){
          best.a = u;
          best.b = i;
          bestedist = list->matches[u][i].edist;
        }
      }
    }
  }
  assert(bestedist != -1);

  al = ALLOCMEMORY(space, NULL, Alignment, 1);
  copyAlignment(al, list->matches[best.a][best.b].al);
  
  res = bl_gmatchlistInit(space, maxedist, 0);
  res = se_kdMatchListAdd(res, list->matches[best.a][best.b].subject,
			  list->matches[best.a][best.b].p,
			  list->matches[best.a][best.b].q,
			  list->matches[best.a][best.b].edist,
			  list->matches[best.a][best.b].scr,
			  list->matches[best.a][best.b].i,
			  list->matches[best.a][best.b].j,
			  list->matches[best.a][best.b].evalue,
			  al, best.a, -1, -1, 0, -1, -1, 0, 0);
    
  return res;
}

/*-------------------------- se_kdFindBestMatePair ---------------------------
 *    
 * @brief find the 'best' pair from two lists of hits for query and mate
 * @author Steve Hoffmann 
 *   
 */

gmatchlist_t*
se_kdFindBestMatePair (void *space, gmatchlist_t *querylist, 
    gmatchlist_t *matelist, Uint maxedist, Uint matemaxedist) {

  Uint u, v, i, j, ucnt=0, vcnt=0, uprime=0, vprime=0;
  unsigned char found = 0, downstream = 0;
  Alignment *al;
  PairUint p,q;
  gmatchlist_t *list;

  list = NULL;
  memset(&p, 0, sizeof(PairUint));
  memset(&q, 0, sizeof(PairUint));

  for(u=0; u < 2; u++) {
    ucnt += querylist->n[u];
    for(i=0; i < querylist->n[u]; i++) {
      uprime = u;
      for(v=0; v < 2; v++) {
        vcnt += matelist->n[v];
        for(j=0; j < matelist->n[v]; j++) {
          vprime = v;
          if(querylist->matches[u][i].subject == 
              matelist->matches[v][j].subject) {
            if(!found || abs((LLint)querylist->matches[p.a][p.b].p - 
                  matelist->matches[q.a][q.b].p) >
                abs((LLint)querylist->matches[u][i].p -
                  matelist->matches[v][j].p)) 
            {
              found = 1;
              p.a = u;
              p.b = i;
              q.a = v;
              q.b = j;
            }
          }
        }
      }
    }
  }


  if(!found && ucnt == 1 && vcnt == 1) { 
    p.a = uprime;
    p.b = 0;
    q.a = vprime;
    q.b = 0;
    found=1;
  }


  if(found) {

    al = ALLOCMEMORY(space, NULL, Alignment, 1);
    copyAlignment(al, querylist->matches[p.a][p.b].al);
   
    list = bl_gmatchlistInit(space, maxedist, matemaxedist);  
    //fprintf(stdout, "adding match with at %d with edist %d\n", querylist->matches[p.a][p.b].p, querylist->matches[p.a][p.b].edist);
    
    list = se_kdMatchListAdd(list,  
        querylist->matches[p.a][p.b].subject,
        querylist->matches[p.a][p.b].p,
        querylist->matches[p.a][p.b].q,
        querylist->matches[p.a][p.b].edist,
        querylist->matches[p.a][p.b].scr,
        querylist->matches[p.a][p.b].i,
        querylist->matches[p.a][p.b].j-1,
        querylist->matches[p.a][p.b].evalue, al, p.a, -1, -1, 0, -1, -1, 0, 0);
        
    if(querylist->matches[p.a][p.b].p >= matelist->matches[q.a][q.b].p) 
      downstream = 1; 
    else 
      downstream = 0; 

    al = ALLOCMEMORY(space, NULL, Alignment, 1);
    copyAlignment(al, matelist->matches[q.a][q.b].al);


   // fprintf(stdout, "adding mate with at %d with edist %d\n",matelist->matches[q.a][q.b].p, matelist->matches[q.a][q.b].edist);
    se_kdSetMate(space, &list->matches[p.a][0], 
//        list->matches[p.a][0].subject, 
        matelist->matches[q.a][q.b].subject,
        matelist->matches[q.a][q.b].p, 
        matelist->matches[q.a][q.b].q, 
        matelist->matches[q.a][q.b].edist, 
        al, downstream, (p.a != q.a));

    if(list->mateminedist > matelist->matches[q.a][q.b].edist) {
      list->mateminedist = matelist->matches[q.a][q.b].edist;
    }
   
  // fprintf(stdout, "list mate min edist %d\n", list->mateminedist);

    if(list->pairminedist > matelist->matches[q.a][q.b].edist + 
        querylist->matches[p.a][p.b].edist) {
      list->pairminedist = matelist->matches[q.a][q.b].edist + 
        querylist->matches[p.a][p.b].edist;
    }
   
  //  fprintf(stdout, "list pair min edist %d\n", list->mateminedist);

    querylist->matches[p.a][p.b].skip = 1;
    matelist->matches[q.a][q.b].skip = 1; 
  }

  return list;
}


/*------------------------------ se_kdAlignMate ------------------------------
 *    
 * @brief find the mate once a sequence was located
 * @author Steve Hoffmann 
 *   
 */


  Uint
se_kdAlignMate(void *space, MultiCharSeq *seq, char **seqs, Uint len, 
    gmatchlist_t *list, Uint maxedist,Uint* enctab, bitvector *D, Uint maxlen) 
{

  PairSint mb;
  Alignment *al;
  bitvector *peq[2];
  char *refseq, *upstreamrefseq;
  Uint u, i, k, p, q;
  Uint idx, refstart, reflen, upstreamreflen, 
       upstreamrefstart, chrstart, chrend;
  Uint edist;

  peq[0] = getpeq(space, seqs[0], len, seq->map, 
      seq->mapsize, enctab);
  peq[1] = getpeq(space, seqs[1], len, seq->map, 
      seq->mapsize, enctab);

  list->mateminedist = maxedist;
  list->pairminedist = list->minedist+maxedist;


  for(u=0; u < 2; u++) {
    for (i=0; i < list->n[u]; i++) {

      idx = list->matches[u][i].subject;
      getMultiCharSeqIdxBounds(seq, idx, &chrstart, &chrend);

      p = list->matches[u][i].p;
      q = list->matches[u][i].q;

      refstart = p;
      reflen = (chrend > (Lint)refstart + maxlen)? maxlen :(chrend-refstart);
      refseq = &seq->sequences[refstart];

      upstreamreflen =((Lint)q-maxlen > chrstart)? maxlen :(Lint)q -chrstart;
      upstreamrefstart = q - upstreamreflen;
      upstreamrefseq = &seq->sequences[upstreamrefstart];

      for(k=0; k < 2; k++) {

        myersbitmatrix(NULL, seqs[k], len, refseq, reflen, 
            seq->map, seq->mapsize, enctab, 
            len-maxedist, peq[k], &mb, D, reflen);
          

        if (mb.a != -1 && mb.b <= maxedist && mb.a < reflen) {  
          al = ALLOCMEMORY(space, NULL, Alignment, 1);

          initAlignment(al, seqs[k], len, 0, refseq, reflen, 0);
          bitvectorbacktrack(al, D, reflen, len, mb.a);

          edist = se_kdSetMate(space, &list->matches[u][i], idx, 
              refstart+al->voff, refstart+mb.a-1, 
              mb.b, al, 1, (u != k));

          if(list->mateminedist > edist) {
            list->mateminedist = edist;
          }

          //pairminedist ... 
          if(list->pairminedist > edist+list->matches[u][i].edist) {
            list->pairminedist = edist+list->matches[u][i].edist;
          } 

          list->matches[u][i].noofmatematches++;
        }

        myersbitmatrix(NULL, seqs[k], len, upstreamrefseq, upstreamreflen, 
            seq->map, seq->mapsize, enctab, 
            len-maxedist, peq[k], &mb, D, upstreamreflen);

        
        if (mb.a != -1 && mb.b <= maxedist && mb.a < upstreamreflen) {  
          al = ALLOCMEMORY(space, NULL, Alignment, 1);

          initAlignment(al, seqs[k], len, 0, upstreamrefseq, 
              upstreamreflen, 0);
          bitvectorbacktrack(al, D, upstreamreflen, len, mb.a);
        
          edist = se_kdSetMate(space, &list->matches[u][i], idx, 
              upstreamrefstart+al->voff, upstreamrefstart+mb.a-1, 
              mb.b, al, 0, (u != k));

          if(list->mateminedist > edist) {
            list->mateminedist = edist;
//            list->noofmatematches = 0;
          }// else if(list->mateminedist == mb.b) {
           // list->noofmatematches++;
          //}
          //pairminedist ... 
          if(list->pairminedist > edist+list->matches[u][i].edist) {
            list->pairminedist = edist+list->matches[u][i].edist;
          }

 
          list->matches[u][i].noofmatematches++;
        }
      }
    }
  }

  for(u=0; u < 2; u++) {
    for(i=0; i < seq->mapsize; i++) {
      FREEMEMORY(space, peq[u][i]);
    }  
    FREEMEMORY(space, peq[u]);
  }

  return 0;
}



/*--------------------------- bl_kdUpdateBestSeed ----------------------------
 *    
 * @brief update the best seed record
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_kdUpdateBestSeed (bestseed_t *best, MultiCharSeq *seq, Suffixarray *s, Uint qrylen,
    Uint readstart, Uint mat, Uint mis, Uint ins, Uint del, Uint l, Uint refstrand)
{
  
  Uint pos, subidx, substart, subend;    
  
  pos = s->suftab[l];
  subidx = getMultiCharSeqIndex(seq, &seq->sequences[pos]);
  getMultiCharSeqIdxBounds(seq, subidx, &substart, &subend);
  pos -= substart;

  if(best->mat < mat) {
    best->maxintervalconstraint = 0;
    best->maxevalconstraint = 0;
    best->readstart = (refstrand == 0) ? readstart : qrylen-readstart-mat-mis-ins;
    assert(qrylen >= readstart + 1);
    best->mat= mat;
    best->mis = mis;
    best->ins = ins;
    best->del = del;
    best->len = mat+mis+ins;
    best->refidx = subidx;
    best->refpos = pos;
    best->refstrand = refstrand;
  }
	
  return ;
}

/*--------------------------- se_kdMatchStemAlign ----------------------------
 *    
 * @brief align the seeds in the matchstem using bv algorithmics
 * @author Steve Hoffmann 
 *   
 */


gmatchlist_t*
se_kdMatchStemAlign(void *space, Suffixarray *s, MultiCharSeq *seq, 
    matchstem_t **stems, char **seqs, Uint len, karlin_t *stats, 
    segemehl_t *nfo, Uint *enctab, bitvector* D, bestseed_t *best) {

  Uint u, k, j, l, r, q, i, pos, mat, mis, ins, del; 
  int maxedist, bestedist, scr, skipmargin=0;
  Uint *check=NULL;
  Uint checklen=0;
  double E;
  bitvector *peq[2];
  PairSint mb;
  MultiCharSeqAlignment mcsa;
  gmatchlist_t *list;

  maxedist = bestedist = len - floor(((double)nfo->accuracy * len)/100.);
  skipmargin = 40*((double)maxedist/100.);
  list = bl_gmatchlistInit(space, maxedist, 0);

  peq[0] = getpeq(NULL, seqs[0], len, seq->map, 
      seq->mapsize, enctab);
  peq[1] = getpeq(NULL, seqs[1], len, seq->map, 
      seq->mapsize, enctab);

  for(u = 0; u < 2; u++) {
    for(i = 0; i < len; i++) {
      for(q = 0; q < stems[u][i].noofbranches; q++) {

        l = stems[u][i].branches[q].l; 
        r = stems[u][i].branches[q].r;
        mat = stems[u][i].branches[q].mat;
        mis = stems[u][i].branches[q].mis;
        ins = stems[u][i].branches[q].ins;
        del = stems[u][i].branches[q].del;

        E = kd_getBranchEvalue(stems[u], i, q, len, s->numofsuffixes, stats);
        
        if(l <= r && E > nfo->maxevalue && best->mat == 0) best->maxevalconstraint=1;
        if(l <= r && (r-l) <= nfo->M && best->mat == 0) best->maxintervalconstraint=1;

        if (l > r || E > nfo->maxevalue || (r-l) > nfo->M) 
          continue;
          
        bl_kdUpdateBestSeed(best, seq, s, len, i, mat, mis, ins, del, l, u);

        for(j = l; j <= r; j++) {
          pos = s->suftab[j];

          if(mat != len || mis+ins+del != 0) {
            initMultiCharSeqAlignment(space, &mcsa, seq, pos, 
                i+maxedist, len+2*(maxedist+1), u, NULL, seqs[u], len);
          } else {
            initMultiCharSeqAlignment(space, &mcsa, seq, pos, 
                0, len, u, NULL, seqs[u], len);
          }

          /*skip or update identical matches*/
          for(k = 0; k < checklen; k++) 
            if (check[k] >= mcsa.refstart-skipmargin && 
                check[k] <= mcsa.refstart+skipmargin)     
              break;	

          if (k < checklen) {
            wrapMultiCharSeqAlignment(space, &mcsa);
            continue;
          }

          check = ALLOCMEMORY(space, check, Uint, checklen+1);
          check[checklen++]= mcsa.refstart;

          if (mat == len && mis+ins+del == 0) {
            scr = kd_getBranchScore(stems[u], i, q);

            for(k=0; k < len; k++) {
              insertEop(mcsa.al, Replacement);
            }

            mb.b = getEdist(mcsa.al);

            if (mb.b <= maxedist && mb.b <= bestedist){	      
              
              list = se_kdMatchListAdd(list, mcsa.subidx, 
                  pos, pos+len-1, mb.b, scr, 0, len-1, E, mcsa.al, 
                  u, -1, -1, 0, -1, -1, 0, 0);

              if(nfo->bestonly) bestedist = list->minedist;
            
            } else {
              wrapMultiCharSeqAlignment(space, &mcsa);
            }
          
          } else {

            myersbitmatrix(NULL, seqs[u], len, mcsa.refseq, mcsa.reflen, 
                seq->map, seq->mapsize, enctab, len-bestedist, peq[u], 
                &mb, D, mcsa.reflen);

            if (mb.a != -1 && mb.b <= maxedist && 
                mb.b <= bestedist && mb.a < mcsa.reflen) {  
              bitvectorbacktrack(mcsa.al, D, mcsa.reflen, len, mb.a);

	      mb.b = getEdist(mcsa.al);

	      if (mb.b <= maxedist && mb.b <= bestedist){
		
		/*skip or update identical matches*/
		for(k = 0; k < list->n[u]; k++) 
		  if (list->matches[u][k].p == mcsa.refstart+mcsa.al->voff)
		    break;	
		
		if (k < list->n[u]) { 
		  if (list->matches[u][k].edist <= mb.b){
		    wrapMultiCharSeqAlignment(space, &mcsa);
		  } else {
		    scr = kd_getBranchScore(stems[u], i, q);
		    
		    list = se_kdMatchListSet(space, list, mcsa.subidx, 
                        mcsa.refstart+mcsa.al->voff, 
                        mcsa.refstart+mb.a-1,
                        mb.b, scr, 0, len-1, E, mcsa.al, u, k);
		  }
		  continue;
		}
		
		scr = kd_getBranchScore(stems[u], i, q); 
		
		list=se_kdMatchListAdd(list, mcsa.subidx, 
                    mcsa.refstart+mcsa.al->voff, 
                    mcsa.refstart+mb.a-1, mb.b, scr, 0, len-1, E, mcsa.al, 
                    u, -1, -1, 0, -1, -1, 0, 0);

		if(nfo->bestonly) bestedist = list->minedist;
            
	      } else { 
		wrapMultiCharSeqAlignment(space, &mcsa);
	      }
	    } else {
	      wrapMultiCharSeqAlignment(space, &mcsa);
	    }
          }
        }
      }
    }

    for(j=0; j < seq->mapsize; j++) {
      FREEMEMORY(space, peq[u][j]);
    }

    FREEMEMORY(space, peq[u]);
    if(check) {
      FREEMEMORY(space, check);
      check = NULL;
      checklen = 0;
    }
  }

  return list;
}


/*----------------------------- se_kdCheckTrans ------------------------------
 *    
 * @brief check trans alignments for better alternatives
 * @author Steve Hoffmann 
 *   
 *
 
MultiCharSeqAlignment*
se_kdCheckTrans(MultiCharSeq *seq, MultiCharSeqAlignment *a,  char **seqs, 
    Uint qrylen, gmatchlist_t *list, int *scores, int indel, 
    Uint *enctab, bitvector *D, unsigned int *noofaligns, segemehl_t *nfo)
{

  bitvector *peq[2];
  PairSint mb;
  Alignment *al;

  char *fseq[2], *refseq, *refseq2;
  unsigned int flen, reflen, reflen2, maxedist, w, idx, pos, refpos, u, foff[2], score;
  PairSint scan[2];
  char nextstrand='+', prevstrand='+', strand;
  Uint i, j, k, prevpos, nextidx, nextpos, vlen = 0,
       ustartj = 0,
       ulen=0, ustart=0, nextustart=0, 
       prevustart=0, uend=0; 
  Uint *subidxs=NULL, *subidxscnt;
  char *strands=NULL, *strandscnt;
  Uint maxstrandcnt = 0, maxsubidxscnt = 0;
  Uint maxsubidx=0;
  char maxstrand=0;

  MultiCharSeqAlignment *b=NULL;
  unsigned int nooffragments = *noofaligns;
  void *space = NULL;

  subidxs = ALLOCMEMORY(space, NULL, Uint, nooffragments);
  strands = ALLOCMEMORY(space, NULL, char, nooffragments);
  subidxscnt = ALLOCMEMORY(space, NULL, Uint, nooffragments);
  strandscnt = ALLOCMEMORY(space, NULL, Uint, nooffragments);

  for(i=0; i < nooffragments; i++) {
    subidxs[i] = a[i].subidx;
    strands[i] = a[i].strand;
  }

  qsort(subidxs, nooffragments, sizeof(Uint), cmp_Uint_qsort); 
  qsort(strands, nooffragments, sizeof(char), cmp_char_qsort); 

  for(i=0; i < nooffragments; i++) {

    if(i==0) { 
      strandscnt[0]=1;
      subidxscnt[0]=1;
      continue;
    }

    if(strands[i-1] == strands[i]) {
      strandscnt[i] = strandscnt[i-1]+1;
    } else {
      strandscnt[i] = 1;
    }
    if(subidxs[i-1] == subidxs[i]) {
      subidxscnt[i] = subidxscnt[i-1]+1;
    } else {
      subidxscnt[i] = 1;
    }

    if(maxstrandcnt < strandscnt[i]) {
      maxstrandcnt = strandscnt[i];
      maxstrand = strands[i]; 
    }

    if(maxsubidxscnt < subidxscnt[i]) {
      maxsubidxscnt = subidxscnt[i];
      maxsubidx = subidxs[i];
    }
  }

  for(k=0, i=0; i < nooffragments; i++) {
    if(a[i].subidx != maxsubidx || a[i].strand != maxstrand) {

      ustart = a[i].al->uoff;
      ulen = getUalignlen(a[i].al);
      idx = a[i].subidx;
      vlen = getValignlen(a[i].al);
      score = getAlignScore(a[i].al, scores, indel);
      
      if(a[i].strand == 1) { 
        strand = '-';
        pos =  a[i].refstart + a[i].al->voff + getValignlen(a[i].al) - 1;
      } else {
        strand = '+';
        pos = a[i].refstart + a[i].al->voff;
      }
      if (a[i].strand == 1) {
        uend = qrylen - ustart - 1;
        ustart = uend - ulen + 1;        
      } else {
        uend = ustart + ulen - 1;
      }

    }
  }

  FREEMEMORY(space, subidxs);
  FREEMEMORY(space, strands);
  FREEMEMORY(space, subidxscnt);
  FREEMEMORY(space, strandscnt);

  return NULL;
}
*/

/*-------------------------------- se_kdFixIn --------------------------------
 *    
 * @brief fix in
 * @author Steve Hoffmann 
 *   
 */

MultiCharSeqAlignment*
se_kdFixIn (MultiCharSeq *seq, MultiCharSeqAlignment *a,  char **seqs, Uint qrylen, 
    gmatchlist_t *list, int *scores, int indel, 
    Uint *enctab, bitvector *D, unsigned int *noofaligns, segemehl_t *nfo)
{

  bitvector *peq[2];
  PairSint mb;
  Alignment *al;

  char *fseq[2], *refseq, *refseq2;
  unsigned int flen, reflen, reflen2, maxedist, w, idx, pos, refpos, u, foff[2], score, chrstart, chrend;
  PairSint scan[2];
  char nextstrand='+', prevstrand='+', strand;
  Uint i, j, k, prevpos, nextidx, nextpos, vlen = 0,
       ustartj = 0,
       ulen=0, ustart=0, nextustart=0, 
       prevustart=0, uend=0; 
  MultiCharSeqAlignment *b=NULL;
  unsigned int nooffragments = *noofaligns;
  void *space = NULL;

 

  for(k=0, i=0; i < nooffragments; i++) {

    b = ALLOCMEMORY(space, b, MultiCharSeqAlignment, k+1);
    memmove(&b[k], &a[i], sizeof(MultiCharSeqAlignment));
    k++;

    ustart = a[i].al->uoff;
    ulen = getUalignlen(a[i].al);
    idx = a[i].subidx;
    vlen = getValignlen(a[i].al);
    score = getAlignScore(a[i].al, scores, indel);


    if(a[i].strand == 1) { 
      strand = '-';
      pos =  a[i].refstart + a[i].al->voff + getValignlen(a[i].al) - 1;
    } else {
      strand = '+';
      pos = a[i].refstart + a[i].al->voff;
    }
    if (a[i].strand == 1) {
      uend = qrylen - ustart - 1;
      ustart = uend - ulen + 1;        
    } else {
      uend = ustart + ulen - 1;
    }
    
    
    if (ulen >= nfo->minfragmentalignlen && 
        vlen >= nfo->minfragmentalignlen &&
        score >= nfo->minfragmentalignscore) { 


    prevpos = -1;
    nextidx = -1;
    nextpos = -1;
    prevstrand = -1;
    nextstrand = -1;
    prevustart =0;
    nextustart = 0;


    for(j=0; j < nooffragments; j++) {

      if(a[j].strand == 1) {
        ustartj = qrylen - a[j].al->uoff - getUalignlen(a[j].al);
      } else {
        ustartj = a[j].al->uoff;
      }

      if (ustartj < ustart &&  (!prevustart || ustartj >= prevustart) &&
          getUalignlen(a[j].al) >= nfo->minfragmentalignlen &&
          getValignlen(a[j].al) >= nfo->minfragmentalignlen &&
          getAlignScore(a[j].al, scores, indel) >= nfo->minfragmentalignscore) { 

        if(a[j].strand == 0) { 
          prevpos = a[j].refstart + a[j].al->voff + getValignlen(a[j].al) - 1;
          prevstrand = '+';
          if(a[j].strand == 1) { 
            prevstrand = '-';
          }
        } else {  
          prevpos = a[j].refstart + a[j].al->voff;
          prevstrand = '-';
          if(a[j].strand == 0) { 
            prevstrand = '+';
          }
        }
        prevustart = ustartj;
      }

      if (ustartj > ustart && (!nextustart || ustartj <= nextustart) &&
          getUalignlen(a[j].al) >= nfo->minfragmentalignlen &&
          getValignlen(a[j].al) >= nfo->minfragmentalignlen &&
          getAlignScore(a[j].al, scores, indel) >= nfo->minfragmentalignscore) { 

        nextidx = a[j].subidx;

        if(a[j].strand == 0) { 
          nextpos = a[j].refstart + a[j].al->voff;
          nextstrand = '+';
        } else { 
          nextpos = a[j].refstart + a[j].al->voff + getValignlen(a[j].al) - 1;
          nextstrand = '-';
        }
        nextustart = ustartj;
      }
    }

    if(nextustart && nextustart > uend+5 && nextidx == idx && strand == nextstrand) {


      flen = nextustart - uend - 1;
      foff[0] = uend+1;
      foff[1] = qrylen-nextustart;
      fseq[0] = &seqs[0][foff[0]];
      fseq[1] = &seqs[1][foff[1]];


      refseq = &seq->sequences[MIN(pos,nextpos)]; 
      reflen = MAX(pos,nextpos) - MIN(pos,nextpos);

      maxedist = flen - floor(((double)nfo->accuracy * flen)/100.);
#ifdef DEBUGFIXIN
      fprintf(stdout, "attempt fix in %u]-[%u into genomic interval %u]-[%u; reflen:%u flen:%u\n", 
          uend, nextustart, pos, nextpos, reflen, flen);
#endif
      peq[0] = getpeq(NULL, fseq[0], flen, seq->map, seq->mapsize, enctab);
      peq[1] = getpeq(NULL, fseq[1], flen, seq->map, seq->mapsize, enctab);

      w = a[i].strand;
      scan[w] = myersbitvector(NULL, fseq[w], flen, refseq, reflen, seq->map, 
          seq->mapsize, enctab, maxedist, peq[w]);

      if(scan[w].a != -1) {
      getMultiCharSeqIdxBounds(seq, idx, &chrstart, &chrend);
#ifdef DEBUGFIXIN
        fprintf(stdout, "found matsch at %u wiff edist %d (strand:%d)\n", 
            scan[w].a, scan[w].b, w);
#endif
        refpos = MIN(pos,nextpos);
        refpos += (scan[w].a > flen) ? scan[w].a-flen : 0;
        refpos += (scan[w].a > flen && scan[w].a-flen > 100) ? -100 : 0;
        refpos = (refpos >= chrstart) ? refpos : chrstart; 
        refseq2 = &seq->sequences[refpos];
        reflen2 = (chrend-(refpos+flen) > 200) ? flen + 200 : chrend-(refpos+flen);

#ifdef DEBUGFIXIN          
        fprintf(stdout, "narrow down to [%u,%u] with length %d\n", 
            refpos, refpos+reflen2, reflen2);
#endif          
        myersbitmatrix(NULL, fseq[w], flen, refseq2, reflen2, 
            seq->map, seq->mapsize, enctab, flen-maxedist, peq[w], &mb, D, reflen2);
#ifdef DEGUGFIXIN
        fprintf(stdout, "aligned matsch at %u wiff edist %d (maxedist:%d, reflen2:%d)\n", 
            mb.a, mb.b, maxedist, reflen2);
#endif
        if (mb.a != -1 && mb.b <= maxedist && mb.a < reflen2) {  
          al = ALLOCMEMORY(space, NULL, Alignment, 1);

          initAlignment(al, seqs[w], qrylen, foff[w], refseq2, reflen2, 0);
          bitvectorbacktrack(al, D, reflen2, flen, mb.a);

          if(getUalignlen(al) > nfo->minfragmentalignlen && 
              getValignlen(al) > nfo->minfragmentalignlen && 
              getAlignScore(al, scores, indel) > nfo->minfragmentalignscore) { 

            b = ALLOCMEMORY(space, b, MultiCharSeqAlignment, k+1);

            initMultiCharSeqAlignment(space, &b[k], seq, refpos, 
                0, reflen2, w, NULL, seqs[w], qrylen);

            wrapAlignment(b[k].al);
            FREEMEMORY(space, b[k].al);
            b[k].al = al;
            k++;
#ifdef DEGUBFIXIN    
            showAlign(al, stdout);
#endif
          } else { 
            wrapAlignment(al);
            FREEMEMORY(space, al);
          }
        } 
      }

      for(w=0; w < 2; w++) {
        for(u=0; u < seq->mapsize; u++) {
          FREEMEMORY(space, peq[w][u]);
        }  
        FREEMEMORY(space, peq[w]);
      }
    }

    if(prevpos && prevstrand == nextstrand && prevstrand != strand) {
      //fprintf(stdout, "double check frag\n");
    }
  }
  }

  FREEMEMORY(space, a);

  *noofaligns = k;
  return b;

}


/*------------------------- se_kdAlignEvalSplitAlign -------------------------
 *    
 * @brief post processing of the multisplitalignment
 * @author Steve Hoffmann 
 *   
 */
 
char
se_kdAlignEvalSplitAlign (MultiCharSeq *seq, MultiCharSeqAlignment *a,  char **seqs,  
    Uint qrylen,  gmatchlist_t *list, Uint *totalcover, int *totalscore, 
    unsigned char *trans, int *scores, int indel, Uint *enctab, bitvector *D, 
    unsigned int noofaligns, segemehl_t *nfo)
{

  Alignment *alcopy; 

  unsigned char laststrand=0, purge = 0;
  char nextstrand='+', prevstrand='+';
  Uint i, j, k, previdx, prevpos, nextidx, nextpos, 
       lastsubidx =0 , totaledist = 0, ustartj = 0,
       ulen=0, vlen=0, vstart=0, ustart=0, nextustart=0, 
       prevustart=0, uend=0, edist=0, fragno = 0; 
  int score; 
#ifdef DEBUGTRANSALIGN
  Uint sub_start, sub_end;
#endif
  for(k=0, i=0; i < noofaligns; i++) {
 
    //remove trailing indels
    //if(i == cur->nooffragments && clean5prime) {
    //  clean5prime(a[i].al);
    //}

    ustart = a[i].al->uoff;
    vstart = a[i].al->voff;
    ulen = getUalignlen(a[i].al);
    vlen = getValignlen(a[i].al);
    score = getAlignScore(a[i].al, scores, indel);
    edist = getEdist(a[i].al);


#ifdef DEBUGTRANSALIGN
        getMultiCharSeqIdxBounds(seq, a[i].subidx, &sub_start, &sub_end);

        fprintf(stdout, "frag:%d, [%u,%u], off:%d, [%u,%u], voff:%d, strand:%d, score:%d, edist:%d, ulen:%d, vlen:%d\n", 
        i, a[i].qrystart, a[i].qrystart+a[i].qrylen-1, ustart, 
        a[i].refstart-sub_start, a[i].refstart-sub_start+a[i].reflen-1, 
        vstart, a[i].strand,  score, edist, ulen, vlen);
        
        showAlign(a[i].al, stdout);
#endif

   
    
    if(edist > (ulen - floor(((double)nfo->accuracy * ulen)/100.))) {  
#ifdef DEBUGTRANSALIGN
      fprintf(stdout, "purging!\n");
#endif
      purge = 1;
    }
    
    if (ulen >= nfo->minfragmentalignlen && 
        vlen >= nfo->minfragmentalignlen &&
        score >= nfo->minfragmentalignscore) { 

      *totalcover += ulen;
      *totalscore += score;
      totaledist += edist;
      k++;
      
      alcopy = ALLOCMEMORY(space, NULL, Alignment, 1);
      copyAlignment(alcopy, a[i].al);     
      if (a[i].strand == 1) {
        uend = qrylen - ustart - 1;
        ustart = uend - ulen + 1;        
      } else {
        uend = ustart + ulen - 1;
      }

      previdx = -1;
      prevpos = -1;
      nextidx = -1;
      nextpos = -1;
      prevstrand = -1;
      nextstrand = -1;
      prevustart =0;
      nextustart = 0;


      for(j=0; j < noofaligns; j++) {

        if(a[j].strand == 1) {
          ustartj = qrylen - a[j].al->uoff - getUalignlen(a[j].al);
        } else {
          ustartj = a[j].al->uoff;
        }

        if (ustartj < ustart &&  (!prevustart || ustartj >= prevustart) &&
            getUalignlen(a[j].al) >= nfo->minfragmentalignlen &&
            getValignlen(a[j].al) >= nfo->minfragmentalignlen &&
            getAlignScore(a[j].al, scores, indel) >= nfo->minfragmentalignscore) { 
          
          previdx = a[j].subidx;

          if(a[j].strand == 0) { 
            prevpos = a[j].refstart + a[j].al->voff + getValignlen(a[j].al) - 1;
            prevstrand = '+';
            if(a[j].strand == 1) { 
              prevstrand = '-';
            }
          } else {  
            prevpos = a[j].refstart + a[j].al->voff;
            prevstrand = '-';
            if(a[j].strand == 0) { 
              prevstrand = '+';
            }
          }
          prevustart = ustartj;
        }
   
        if (ustartj > ustart && (!nextustart || ustartj <= nextustart) &&
            getUalignlen(a[j].al) >= nfo->minfragmentalignlen &&
            getValignlen(a[j].al) >= nfo->minfragmentalignlen &&
            getAlignScore(a[j].al, scores, indel) >= nfo->minfragmentalignscore) { 
          
          nextidx = a[j].subidx;
          
          if(a[j].strand == 0) { 
            nextpos = a[j].refstart + a[j].al->voff;
            nextstrand = '+';
          } else { 
            nextpos = a[j].refstart + a[j].al->voff + getValignlen(a[j].al) - 1;
            nextstrand = '-';
          }
          nextustart = ustartj;
        }
      }

      list = se_kdMatchListAdd(list, a[i].subidx, 
          a[i].refstart + vstart, 
          a[i].refstart + vstart + vlen - 1, 
          edist, score, ustart, //ustart + ulen - 1, 
          uend, .0, alcopy, a[i].strand, 
          previdx, prevpos, prevstrand, 
          nextidx, nextpos, nextstrand, fragno);

      fragno++;

      if (k > 1 && (laststrand != a[i].strand || lastsubidx != a[i].subidx)) {
        *trans = 1;
      }

      laststrand = a[i].strand;
      lastsubidx = a[i].subidx;
    } else {
      //fprintf(stdout, "minalignlen and score failure\n");
      //purge =1;
    } 
  }

  if(fragno < 2) purge = 1;
  
  return purge;
}

/*--------------------------- se_kdAlignSplitChain ---------------------------
 *    
 * @brief align a chain of fragments using local multi spliced alignment
 * @author Steve Hoffmann 
 *   
 */

gmatchlist_t*
se_kdAlignSplitChain (void *space, branchChain_t *chains, Uint noofchains,
    Suffixarray *arr, MultiCharSeq *seq, char *querydesc, matchstem_t **stems,
    char **seqs, Uint qrylen, int *scores, int indel, int transition, 
    spliceevents_t *events, Uint *enctab, bitvector *D, segemehl_t *nfo) {

  Uint i, j, start, floff = 0, flen =0, 
    maxedist,  
    *strands, *starts, *ends, *tstarts, *tends, *lengths, *reflens, totalcover = 0,
    sub_start, sub_end;

  unsigned int margin=50, maxmargin=100, noofaligns; //50;
  unsigned char trans=0, purge=0;
  char **refseqs;
  char ***K = NULL;
  PairUint *bestscr;
  int ***M, **lmr, **lmv, **lmc, totalscore=0;
  branchChain_t *cur;
  MultiCharSeqAlignment *a;
  Alignment **aligns; 
  gmatchlist_t *list=NULL;
  PairUint *diag;
#ifdef DEBUGKBAND
  char ***B;
#endif



  maxedist = qrylen - floor(((double)nfo->accuracy * qrylen)/100.);
  list = bl_gmatchlistInit(space, maxedist, 0);

  //DBG("kdalignsplitchain with %d chains", noofchains);

  if(noofchains == 0) return list;

  qsort(chains, noofchains, sizeof(branchChain_t), cmp_chainscores); 
  cur = &chains[0];

#ifdef DEBUGTRANSALIGN  
  showChains(chains, noofchains, arr, stdout, seqs[1], qrylen);
  fprintf(stdout, "nooffrags: %d; scr1:%d scr:2:%d\n", cur->nooffragments, 
      chains[0].score, chains[1].score);
#endif

 

  if(cur->nooffragments <= 1) return list; 

  a = ALLOCMEMORY(space, NULL, MultiCharSeqAlignment, cur->nooffragments);
  reflens = ALLOCMEMORY(space, NULL, Uint, cur->nooffragments);
  strands = ALLOCMEMORY(space, NULL, Uint, cur->nooffragments);
  starts = ALLOCMEMORY(space, NULL, Uint, cur->nooffragments);
  ends = ALLOCMEMORY(space, NULL, Uint, cur->nooffragments);
  tstarts = ALLOCMEMORY(space, NULL, Uint, cur->nooffragments);
  tends = ALLOCMEMORY(space, NULL, Uint, cur->nooffragments);
  lengths = ALLOCMEMORY(space, NULL, Uint, cur->nooffragments);
  refseqs = ALLOCMEMORY(space, NULL, char*, cur->nooffragments);
  aligns =  ALLOCMEMORY(space, NULL, Alignment*, cur->nooffragments);
  diag = ALLOCMEMORY(space, NULL, PairUint, cur->nooffragments);



 /*attention cur is a new chain and includes the beststarts already. 
   * no need for bd[beststart] beyond this call*/
  cur = condenseChain(cur, 1, seq, arr);
  
  for(i=0; i < cur->nooffragments; i++) {

    Uint uloff;
    Uint uroff;
 
    if(i > 0) {
       uloff = (cur->f[i-1]->end < cur->f[i]->start) ? cur->f[i]->start - cur->f[i-1]->end : 0; 
    } else {
       uloff = cur->f[i]->start;
    }

    if(i < cur->nooffragments-1) {
       uroff = (cur->f[i]->end < cur->f[i+1]->start) ? cur->f[i+1]->start - cur->f[i]->end : 0; 
    } else {
       uroff = qrylen-cur->f[i]->end;
    }

    if(cur->f[i]->strand) {
       floff = uroff + MIN(maxedist + margin, maxmargin);
       flen = floff + (cur->f[i]->end - cur->f[i]->start) + uloff + MIN(maxmargin, maxedist + margin);
    } else { 
       floff = maxedist + uloff + margin;
       flen = floff + (cur->f[i]->end - cur->f[i]->start) + uroff + MIN(maxmargin, maxedist + margin);
    }

#ifdef DEBUGTRANSALIGN
    fprintf(stdout, "strand:%d floff:%d\tflen:%d\t[%d,%d]->%d\n", 
        cur->f[i]->strand, floff, flen, cur->f[i]->start, 
        cur->f[i]->end, (i < cur->nooffragments-1) ? cur->f[i+1]->start : qrylen);
#endif

    start = cur->f[i]->substart;
    initMultiCharSeqAlignmentOpt(space, &a[i], seq, start,
        querydesc, seqs[cur->f[i]->strand], cur->f[i]->start, cur->f[i]->end,
        qrylen, floff, flen, uloff, uroff, MIN(maxmargin, maxedist+margin), cur->f[i]->strand);

    aligns[i]  = a[i].al;
    refseqs[i] = a[i].refseq;
    reflens[i] = a[i].reflen; 
    strands[i] = cur->f[i]->strand;
    lengths[i] = a[i].qrylen;
    tstarts[i] = a[i].qrystart;
    tends[i] = tstarts[i]+lengths[i]-1;
    
    if(strands[i]==0) { 
      starts[i] = a[i].qrystart;
      ends[i] = starts[i]+lengths[i]-1;
    } else {
      starts[i] = qrylen - (a[i].qrystart + lengths[i]);
      assert(qrylen >= a[i].qrystart+lengths[i]);
      ends[i] = starts[i]+lengths[i]-1;
      assert(ends[i] <=  qrylen);
    }

    if(strands[i]==0) { 
    diag[i].b = a[i].floff;
    } else {
       if(a[i].reflen >= MIN(maxmargin, maxedist+margin) + uroff + (cur->f[i]->end - cur->f[i]->start) - 1){ 
         diag[i].b = a[i].reflen - MIN(maxmargin, maxedist+margin) - uroff - (cur->f[i]->end - cur->f[i]->start) + 1;
       } else {
         diag[i].b = 0;
       }
    } 
    
    diag[i].a = cur->f[i]->start - a[i].qrystart;


// -DDEBUGMULTISPLICEOPT -DDEBUGTRANSALIGN  
#ifdef DEBUGTRANSALIGN 
    fprintf (stdout, "query sequence of fragment %d [%d->%d]\n", i, starts[i], ends[i]);
   
    Uint h=0;
    for(h=0; h < lengths[i]; h++) {
      if(h && (h%60) == 0) fprintf(stdout, "\n");
      fprintf(stdout, "%c",  a[i].query[starts[i]+h]);
    }
    fprintf(stdout,"\n");
   
    fprintf (stdout, "reference sequence of fragment %d\n", i);
    for(h=0; h < reflens[i]; h++) {
      if(h && (h%60) == 0) fprintf(stdout, "\n");
      fprintf(stdout, "%c",  refseqs[i][h]);
    }
    fprintf(stdout,"\n");

    fprintf(stdout, "%s\n qrylen:%d, fragment:%d, start:%d, strand:%d, curstart:%d, curend:%d, maxedist:%d mapping to [%d,%d]\n", 
            querydesc, qrylen, i, start, strands[i],  starts[i], ends[i], maxedist, a[i].refstart, a[i].refstart+a[i].reflen -1);
    fprintf(stdout, "\n");
#endif
  }

  M = localmultisplicedmatrixopt(space, seqs[0], seqs[1], qrylen, lengths,
      refseqs, reflens, strands, starts, ends, tstarts, tends, cur->nooffragments, indel, transition,
      constscr, scores, &lmv, &lmr, &lmc, &bestscr, &K, diag);


  if(M == NULL) {
    fprintf(stderr, "empty matrix returned for seqs: '%s'/'%s' (%d)\n", 
        seqs[0], seqs[1], qrylen);

    for(i=0; i < cur->nooffragments; i++) {

      getMultiCharSeqIdxBounds(seq, a[i].subidx, &sub_start, &sub_end);
      fprintf(stderr, "fragment %d: %d in %d[%d,%d] '", 
          i, 1 /*arr->suftab[bd[beststart][i]]*/, a[i].subidx, sub_start, sub_end);
      for(j=0; j< qrylen; j++) fprintf(stderr, "%c", refseqs[i][j]);
      fprintf(stderr, "'(%d) strand:%d\n", reflens[i], strands[i]);
    }
    return list;
  }
  
#ifdef DEBUGKBAND 
  B = 
#endif 
      localmultisplicedtracebackopt(space, M, seqs[0], seqs[1], qrylen, lengths, 
      refseqs, reflens, strands, starts, ends, tstarts, tends, 
      cur->nooffragments, indel, transition, constscr, scores, 
      aligns, lmv, lmr, lmc, bestscr);

  for(i=0; i < cur->nooffragments; i++) {
    FREEMEMORY(space, lmv[i]);
    FREEMEMORY(space, lmr[i]);
    FREEMEMORY(space, lmc[i]);


#ifdef DEBUGKBAND
    Uint uloff;
    Uint uroff;

    if(i > 0) {
       uloff = (cur->f[i-1]->end < cur->f[i]->start) ? cur->f[i]->start - cur->f[i-1]->end : 0; 
    } else {
       uloff = cur->f[i]->start;
    }

    if(i < cur->nooffragments-1) {
       uroff = (cur->f[i]->end < cur->f[i+1]->start) ? cur->f[i+1]->start - cur->f[i]->end : 0; 
    } else {
       uroff = qrylen-cur->f[i]->end;
    }
   
    fprintf(stderr, "matrix %d of %d\n", i, cur->nooffragments);
    fprintf (stderr, "query sequence of fragment %d (%d,%d)[%d->%d], starts:%u ends:%u fragstart:%d fragend:%d, uloff:%d, uroff:%d, floff:%d, reflen:%d\n", i, cur->f[i]->start, cur->f[i]->end, starts[i], ends[i], diag[i].a, diag[i].b, cur->f[i]->start, cur->f[i]->end, uloff, uroff, a[i].floff, a[i].reflen);

    dumprecursionmatrix2D(stderr, M[i], B[i], K[i], lengths[i], reflens[i], &diag[i]);
#endif

    for(j=0; j < lengths[i]+1; j++) {
#ifdef DEBUGKBAND
      FREEMEMORY(space, B[i][j]);
      FREEMEMORY(space, K[i][j]);
#endif
      FREEMEMORY(space, M[i][j]);
    }
#ifdef DEBUGKBAND
    FREEMEMORY(space, B[i]);
    FREEMEMORY(space, K[i]);
#endif
    FREEMEMORY(space, M[i]);
  }
#ifdef DEBUGKBAND
  FREEMEMORY(space, B);
  FREEMEMORY(space, K);
#endif

  FREEMEMORY(space, M);
  FREEMEMORY(space, bestscr);
  FREEMEMORY(space, diag);

  noofaligns = cur->nooffragments;

  a = se_kdFixIn(seq, a, seqs, qrylen, list, scores, indel, enctab, 
      D, &noofaligns, nfo);

  purge = se_kdAlignEvalSplitAlign(seq, a, seqs, qrylen, list, &totalcover, 
      &totalscore, &trans, scores, indel, enctab, D, noofaligns, nfo);

  totalcover *= 100;
  totalcover /= qrylen;

#ifdef DEBUGTRANSALIGN
   fprintf(stdout, "qrylen:%d, totalcover %d, totalscore %d, noofchains %d, mincover:%d, mintotalscore:%d \n", 
      qrylen, totalcover, totalscore, noofchains, nfo->minsplicedaligncover, nfo->minsplicedalignscore);
#endif
  
 
  if(totalscore >= nfo->minsplicedalignscore && 
     totalcover >= nfo->minsplicedaligncover /*&& k > 1*/ 
     && !purge) {
    /*restrictive policy for reporting trans splicing events*/
//    if(!trans || noofchains == 1) {
//      reportSplicedMatch(space, querydesc, b, k, 
//        totalcover, totaledist, totalscore, nfo);
//    }
//    store splice sites internally for online remapping
//    bl_storeSpliceEvent (space, seq, list, events, 0, 100, seqs, querydesc);

  } else { 
    
//    fprintf(stdout, "destructing list\n");
    bl_gmatchlistDestruct(space, list);
    list = bl_gmatchlistInit(space, maxedist, 0);
  }

  for(i=0; i < noofaligns; i++) { 
    wrapMultiCharSeqAlignment(space, &a[i]);
  }

  wrapChains(space, cur, 1);
  FREEMEMORY(space, cur);


  FREEMEMORY(space, reflens);
  FREEMEMORY(space, refseqs);
  FREEMEMORY(space, strands);
  FREEMEMORY(space, starts);
  FREEMEMORY(space, ends);
  FREEMEMORY(space, tstarts);
  FREEMEMORY(space, tends);
  FREEMEMORY(space, lengths);

 
  FREEMEMORY(space, lmv);
  FREEMEMORY(space, lmr);
  FREEMEMORY(space, lmc);
  FREEMEMORY(space, aligns);
  FREEMEMORY(space, a);

  return list;
}


/*------------------------------ se_kdSplitRead ------------------------------
 *    
 * @brief find the splits of a chimeric reads from matchstem data
 * @author Steve Hoffmann 
 *   
 */

gmatchlist_t*
se_kdSplitRead(void *space, Suffixarray *arr, MultiCharSeq *seq, 
    char *querydesc, matchstem_t **stems, char **seqs, Uint len, 
    karlin_t *stats, spliceevents_t *events, Uint *enctab, bitvector *D, segemehl_t *nfo) 
{
  int indel = -2;
  int transition = -10;
  int scores[]={1, -2};
  Uint noofchains;
  gmatchlist_t* list;
  branchfragment_t* fragments;
  branchChain_t *chains;


  chains = branchChain(space, arr, stems, seqs, len, stats, 
      &noofchains, &fragments, nfo->maxsplitevalue);


  list = se_kdAlignSplitChain (space, chains, noofchains,
      arr, seq, querydesc, stems, seqs, len, scores, indel, transition, events, enctab, D, nfo);

  wrapChains(space, chains, noofchains);
  FREEMEMORY(space, fragments);
  FREEMEMORY(space, chains);

  return list;
}


/*--------------------------------- se_clip ----------------------------------
 *    
 * @brief clipping sequences
 * @author Steve Hoffmann 
 *   
 */
 
void
se_clip (void *space, fasta_t *reads, Uint elem, segemehl_t *nfo)
{

  if(nfo->hardclip3Prime || nfo->hardclip5Prime) {
    bl_fastaHardClip(space, reads, elem, nfo->hardclip5Prime, 
        nfo->hardclip3Prime);
    if(bl_fastaHasMate(reads)) {
      bl_fastaMateHardClip(space, reads, elem, nfo->hardclip5Prime, 
          nfo->hardclip3Prime);
    } 
  } 

  if(nfo->softclip3Prime || nfo->softclip5Prime) {

    bl_fastaSoftClip(space, reads, elem, 
        nfo->softclip5Prime, nfo->softclip5PrimeLen, nfo->minclipscr5,
        nfo->softclip3Prime, nfo->softclip3PrimeLen, nfo->clipacc, nfo->polyAlen);
    if(bl_fastaHasMate(reads)) {
      bl_fastaMateSoftClip(space, reads, elem, 
          nfo->softclip5Prime, nfo->softclip5PrimeLen, nfo->minclipscr5,
          nfo->softclip3Prime, nfo->softclip3PrimeLen, nfo->clipacc, nfo->polyAlen);
    }
  }

  return ;
}

/*----------------------------- se_kdGenomeMatch -----------------------------
 *    
 * @brief map reads to the genome
 * @author Steve Hoffmann 
 *   
 */

void
se_kdGenomeMatch(void *space, Suffixarray *s, fasta_t *reads, 
    segemehl_t *nfo) {

  unsigned char matchflag, matematchflag;
  matchstatus_t pairStatus = QUERY;
  char *seqs[2], *mateseqs[2], rep=0;
  Uint k, i, u, *enctab, dim, wordno, len, matelen=0, jump, maxedist, 
       matemaxedist=0; //, setmatches;
  karlin_t stats;
  bitvector *D, *Mv;
  Gmap map;
  gread_t read;
  matchstem_t *stems[2] = {NULL, NULL}, *matestems[2] = {NULL, NULL},
              *b0[2], *mateb0[2];
  gmatchlist_t *list=NULL, *matelist=NULL, *templist, 
               *bestpairlist=NULL, *slist=NULL, *slist2=NULL;
  spliceevents_t *events;
  bestseed_t best, bestmate;
  PairSint frag;

  events = ALLOCMEMORY(space, NULL, spliceevents_t, 1);
  events->noofevents = 0;
  events->event = NULL;

  enctab = encodetab(nfo->seq->map, nfo->seq->mapsize);
  dim = reads->maxlen+1000;

  if(bl_fastaHasMate(reads)) {
    dim += nfo->maxinsertsize;
  }

  dim += 2*((reads->maxlen-floor(((double)nfo->accuracy*reads->maxlen)/100.))+4);
  wordno = (reads->maxlen/BITVECTOR_WORDSIZE)+1; 

  D = ALLOCMEMORY(space, NULL, bitvector, 2*(dim+1));
  Mv = &D[dim+1];

  for(i=0; i <= dim; i++) {
    D[i] = initbitvector(space, wordno*BITVECTOR_WORDSIZE);
    Mv[i] = initbitvector(space, wordno*BITVECTOR_WORDSIZE);
  }  

  karlinunitcostpp(space, &stats.lambda, &stats.H, &stats.K);

  for (k=0; k < reads->noofseqs; k++) {
    pairStatus = QUERY;
    matchflag = 0;
    matematchflag = 0;

    best.mat = 0;
    best.mis = 0;
    best.ins =0;
    best.del = 0;
    best.len = 0;
    best.readstart = 0;
    best.refidx = 0;
    best.refpos = 0; 
    best.refstrand = 0;
    best.maxevalconstraint = 0;
    best.maxintervalconstraint = 0;

    bestmate.mat = 0;
    bestmate.mis = 0;
    bestmate.ins = 0;
    bestmate.del = 0;
    bestmate.len = 0;
    bestmate.readstart = 0;
    bestmate.refidx = 0;
    bestmate.refpos = 0; 
    bestmate.refstrand = 0;
    bestmate.maxevalconstraint = 0;
    bestmate.maxintervalconstraint = 0;



    if(!nfo->mute) se_updateProgressBar(k, nfo);
    se_clip(space, reads, k, nfo);

  
    seqs[0] = bl_fastaGetSequence(reads, k);
    len = bl_fastaGetSequenceLength(reads, k);

#ifdef HASHING
    if (bl_fastaGetQuantity(reads, k) == 1){
      DBG("%u: %s\t%u\n",  k, bl_fastaGetSequence(reads, k), bl_fastaGetQuantity(reads, k));
    }

    //    fprintf(nfo->dev,"@%s\n%s\n+\n%s\n", 
    //        bl_fastaGetDescription(reads,k), seqs[0], bl_fastaGetQuality(reads,k));
    //    continue;
    continue;
#endif
    //    pthread_mutex_lock(nfo->mtx2);
    //    fprintf(nfo->dev, "%s\n", bl_fastaGetDescription(reads,k));
    //    fprintf(nfo->dev, "%s\n", bl_fastaGetMateDescription(reads,k));
    //    pthread_mutex_unlock(nfo->mtx2);


    if(len >= nfo->minsize) {  
      seqs[1] = charIUPACcomplement(space, seqs[0], len);

      /* convert for seed search */
      if (nfo->bisulfite){    
        seqs[0] = ALLOCMEMORY(space, NULL, char, len+1);
        memmove(seqs[0], bl_fastaGetSequence(reads, k), len+1);
        bl_convertBisulfite(seqs[0], len, nfo->bisulfite, 1);
        bl_convertBisulfite(seqs[1], len, nfo->bisulfite, 1);
      }

      initGmap(&map, nfo->seq, 1);
      initRead(&read, k);

      if (nfo->jump == 0) {
        jump = floor(len/75) * 2;
        jump = (jump > 0) ? jump : 1;
        jump = MIN(jump, 15); //limit jump to 15
      } else {
        jump = nfo->jump;
      }

      stems[0] = NULL; stems[1] = NULL;
      b0[0] = NULL; b0[1] = NULL;

      /* restrict search to one strand */
      for (u = 0; u < 2; u++){
        /* nfo->strand == 1 : search only on plus strand
         * => init stems[1] as empty
         * nfo->strand == 2 : search only on minus strand
         * => init stems[0] as empty
         * Note: empty but initialized stems are ignored
         * in function kdbest
         */
        if (nfo->strand == 2 - u){	 
          stems[u] = ALLOCMEMORY(space, NULL, matchstem_t, len);
          for (i = 0; i < len; i++){
            stems[u][i].branches = NULL;
            stems[u][i].noofbranches = 0;
          }
        }
      }

      /*
       * try to find full match, only possible
       * if there are no more than k_p
       * unmatchable characters are in read
       */
      if (nfo->bestonly && countNonMatchingChars(seqs[0], len) <= nfo->k_p){
        kdbest(space, s, seqs, len, nfo->s_ext, nfo->p_mis,
            nfo->Xoff, nfo->k_p, stems, b0);
      }

      if (stems[0] == NULL){
        stems[0]=kdseeds(space, s, seqs[0], len, jump, nfo->s_ext, nfo->p_mis,
            nfo->Xoff, nfo->k_p, b0[0]);
      }
      if (stems[1] == NULL){
        stems[1]=kdseeds(space, s, seqs[1], len, jump, nfo->s_ext, nfo->p_mis,
            nfo->Xoff, nfo->k_p, b0[1]);
      }

      /* convert for alignment */
      if (nfo->bisulfite){
        FREEMEMORY(space, seqs[1]);
        memmove(seqs[0], bl_fastaGetSequence(reads, k), len+1);
        seqs[1] = charIUPACcomplement(space, seqs[0], len);
        bl_convertBisulfite(seqs[0], len, nfo->bisulfite, 0); 
        bl_convertBisulfite(seqs[1], len, nfo->bisulfite, 0); 
      }

      list = se_kdMatchStemAlign(space, s, nfo->seq, stems, seqs,
          len, &stats, nfo, enctab, D, &best);

      if (bl_fastaHasMate(reads)) {
        bestpairlist = NULL;

        mateseqs[0] = bl_fastaGetMate(reads, k);
        matelen = bl_fastaGetMateLength(reads, k); 

        mateseqs[1] = charIUPACcomplement(space, mateseqs[0], matelen); 

        /* convert for direct mate alignment */
        if (nfo->bisulfite){
          mateseqs[0] = ALLOCMEMORY(space, NULL, char, matelen+1);
          memmove(mateseqs[0], bl_fastaGetMate(reads, k), matelen+1);
          bl_convertBisulfite(mateseqs[0], matelen, nfo->bisulfite, 0);
          bl_convertBisulfite(mateseqs[1], matelen, nfo->bisulfite, 0);
        }   

        if(se_kdMatchListhasMatches(list)) {
          matemaxedist = matelen - floor(((double)nfo->accuracy * matelen)/100.);
          se_kdAlignMate(space, nfo->seq, mateseqs, matelen, 
              list, matemaxedist, enctab, D, nfo->maxinsertsize);
        }

        if (se_kdMatchListhasMatches(list) &&
            se_kdMatchListhasMates(list)) {
          /*pair is fully matched*/
          pairStatus = PAIR;
        } else { 
          /*try to find mate first*/
          if (nfo->jump == 0) {
            jump = floor(matelen/75) * 2;
            jump = (jump > 0) ? jump : 1;
          } else {
            jump = nfo->jump;
          }

          /* convert for mate seed search */
          if (nfo->bisulfite){
            FREEMEMORY(space, mateseqs[1]);
            memmove(mateseqs[0], bl_fastaGetMate(reads, k), matelen+1);
            mateseqs[1] = charIUPACcomplement(space, mateseqs[0], matelen);
            bl_convertBisulfite(mateseqs[0], matelen, nfo->bisulfite, 1);
            bl_convertBisulfite(mateseqs[1], matelen, nfo->bisulfite, 1);
          } 

          matestems[0] = NULL; matestems[1] = NULL;
          mateb0[0] = NULL; mateb0[1] = NULL;

          /* restrict search to one strand */
          for (u = 0; u < 2; u++){
            /* nfo->strand == 1 : search only on plus strand
             * => search for mate only on minus strand
             * => init stems[0] as empty
             * nfo->strand == 2 : search only on minus strand
             * => search for mate only on plus strand
             * => init stems[1] as empty
             * Note: empty but initialized stems are ignored
             * in function kdbest
             */
            if (nfo->strand == u + 1){
              matestems[u] = ALLOCMEMORY(space, NULL, matchstem_t, matelen);
              for (i = 0; i < matelen; i++){
                matestems[u][i].branches = NULL;
                matestems[u][i].noofbranches = 0;
              }
            }
          }

          /* 
           * try to find full match, only possible
           * if there are no more than k_p
           * unmatchable characters are in read
           */
          if (nfo->bestonly && countNonMatchingChars(mateseqs[0], matelen) <= nfo->k_p){
            kdbest(space, s, mateseqs, matelen, nfo->s_ext, nfo->p_mis,
                nfo->Xoff, nfo->k_p, matestems, mateb0);
          }

          if (matestems[0] == NULL){
            matestems[0]=kdseeds(space, s, mateseqs[0], matelen, 
                jump, nfo->s_ext, nfo->p_mis,
                nfo->Xoff, nfo->k_p, mateb0[0]);
          }
          if (matestems[1] == NULL){
            matestems[1]=kdseeds(space, s, mateseqs[1], matelen, 
                jump, nfo->s_ext, nfo->p_mis,
                nfo->Xoff, nfo->k_p, mateb0[1]);
          }

          /* convert for mate alignment */
          if (nfo->bisulfite){
            FREEMEMORY(space, mateseqs[1]);
            memmove(mateseqs[0], bl_fastaGetMate(reads, k), matelen+1);
            mateseqs[1] = charIUPACcomplement(space, mateseqs[0], matelen);
            bl_convertBisulfite(mateseqs[0], matelen, nfo->bisulfite, 0);
            bl_convertBisulfite(mateseqs[1], matelen, nfo->bisulfite, 0);
          }   

          matelist = se_kdMatchStemAlign(space, s, nfo->seq, matestems, 
              mateseqs, matelen, &stats, nfo, enctab, D, &bestmate);

          maxedist = len - floor(((double)nfo->accuracy * len)/100.);

          se_kdAlignMate(space, nfo->seq, seqs, len, matelist, 
              maxedist, enctab, D, nfo->maxinsertsize);

          if (se_kdMatchListhasMatches(matelist) && 
              !se_kdMatchListhasMates(matelist) &&
              !se_kdMatchListhasMatches(list)) {
            /*query remains unmatched*/
            pairStatus = MATE;
          }

          if(!se_kdMatchListhasMatches(matelist) &&
              se_kdMatchListhasMatches(list)) {
            /*mate remains unmatched*/
            pairStatus = QUERY;
          }

          if(se_kdMatchListhasMatches(list) &&
              se_kdMatchListhasMatches(matelist) && 
              !se_kdMatchListhasMates(matelist)) {
            /*pair not aligned properly but we have hits (long indel!)*/
            maxedist = len - floor(((double)nfo->accuracy * len)/100.);
            matemaxedist = matelen - floor(((double)nfo->accuracy * matelen)/100.);
            bestpairlist = se_kdFindBestMatePair(space, list, matelist, maxedist, matemaxedist);
            pairStatus = PAIR_INS;
          } 

          if (se_kdMatchListhasMatches(matelist) && 
              se_kdMatchListhasMates(matelist)) {
            /*pair is fully matched in reverse order*/
            templist = list;
            list = matelist;
            matelist = templist;
            pairStatus = PAIR_REV;
          }
        }
      }

      if (nfo->bestonly) {
        maxedist = list->minedist;
        if(matelist) matemaxedist = matelist->minedist;
      } else {
        maxedist = len - floor(((double)nfo->accuracy * len)/100.);
        if(matelist) matemaxedist = matelen - floor(((double)nfo->accuracy * matelen)/100.);
      }

//      if(rep) fprintf(nfo->dev, "pair status %d\n", pairStatus);

      matchflag = 0;
      matematchflag = 0;
      setReads(&map, &read, 1);

      /*report: single ends, fully matched pairs*/
      if(!bl_fastaHasMate(reads) || pairStatus == PAIR_REV || 
          pairStatus == PAIR) {

//        if(rep) fprintf(nfo->dev, "PAIR edist:%d mateedist:%d pairedist:%d (Pair:%d, Rev:%d)\n", maxedist, matemaxedist, list->pairminedist, pairStatus == PAIR, pairStatus == PAIR_REV);
        //if (list->n[0] || list->n[1]) matchflag = 1;
        se_setMatches(space, &read, list, maxedist, nfo, rep);
        matchflag = reportMatch(space, &map, reads, nfo, pairStatus, pairStatus == PAIR_REV);
        se_destructMatches(space, &read); 
      }

      /*report: spliced single ends */
      if(nfo->split && !bl_fastaHasMate(reads) && 
          !se_kdMatchListhasMatches(list)) {

//        if(rep) fprintf(nfo->dev, "SINGLE\n");
        slist = se_kdSplitRead(space, s, nfo->seq, 
            bl_fastaGetDescription(reads, k), 
            stems, seqs, len, &stats, events, enctab, D, nfo);

        se_setMatches(space, &read, slist, maxedist, nfo, rep);
        matchflag = reportMatch(space, &map, reads, nfo, pairStatus, 0);
        se_destructMatches(space, &read); 
        bl_gmatchlistDestruct(space, slist);
      }

      /*report: bestpair from two separately calculated match lists*/
      if(pairStatus == PAIR_INS && bestpairlist) {

//        if(rep) fprintf(nfo->dev, "PAIRINS\n");
        se_setMatches(space, &read, bestpairlist, maxedist, nfo, rep);
        //matchflag = 1;
        //matematchflag = 1;
        matchflag = reportMatch(space, &map, reads, nfo, pairStatus, 0);
        se_destructMatches(space, &read); 
        bl_gmatchlistDestruct(space, bestpairlist);
      }

      /*report: spliced unmatched mate pairs*/
      if(bl_fastaHasMate(reads) && 
          (pairStatus == MATE || pairStatus == QUERY)) {

//        if(rep) fprintf(nfo->dev, "MATE OR QUERY\n");
        if (nfo->split) {

          slist = NULL;
          slist2 = NULL;

          if(!se_kdMatchListhasMatches(list)) {
            slist = se_kdSplitRead(space, s,  nfo->seq,
                bl_fastaGetDescription(reads, k), 
                stems, seqs, len, &stats, events, enctab, D, nfo);

//            if(rep) fprintf(nfo->dev, "ATTEMPT SPLICING FOR QUERY %d\n", (slist->n[0] || slist->n[1]) );
          } 



          if(!se_kdMatchListhasMatches(matelist)) {
            slist2 = se_kdSplitRead(space, s, nfo->seq,
                bl_fastaGetMateDescription(reads, k),
                matestems, mateseqs, matelen, &stats, events, enctab, D, nfo);

//            if (rep) fprintf(nfo->dev, "ATTEMPT SPLICING FOR MATE: %d\n", (slist2->n[0] || slist2->n[1]) );
          } 



          if(slist && se_kdMatchListhasMatches(slist) 
              && (!slist2 || !se_kdMatchListhasMatches(slist2))) {
            /*spliced query full mate*/ 
//            if(rep)     fprintf(nfo->dev, "SPLICED QUERY FULL MATE\n");

            if (se_kdMatchListhasMatches(matelist)) {
              /*report the full mate match*/
//              if(rep)       fprintf(nfo->dev, "MATELIST HAS MATCHES\n");
              pairStatus = QUERY_SPL_FULL_MATE;

	      /* select best mate to spliced query */
	      if (matelist->n[0] + matelist->n[1] > 1){

                /* get last fragment */
                frag.a = frag.b = -1;
                for (u = 0; u < 2; u++){
                  for (i = 0; i < slist->n[u]; i++){                    
                    if (slist->matches[u][i].fragno == 
                        slist->n[0]+slist->n[1]-1){
                      frag.a = u;
                      frag.b = i;
                    }
                  }
                }
                assert(frag.a != -1 && frag.b != -1);
		bestpairlist = se_kdFindBestMate(space, matelist, &slist->matches[frag.a][frag.b], matemaxedist);
		bl_gmatchlistDestruct(space, matelist);
		matelist = bestpairlist;
	      }
              se_setMatches(space, &read, matelist, matemaxedist, nfo, rep);
              matematchflag = reportMatch(space, &map, reads, nfo, pairStatus, 1);
              se_destructMatches(space, &read); 
              /*spliced query no mate*/
            } else {
              pairStatus = QUERY_SPL_NO_MATE;
            }

            se_setMatches(space, &read, slist, maxedist, nfo, rep);
            matchflag = reportMatch(space, &map, reads, nfo, pairStatus, 0);  

            if(matematchflag && matchflag) matchflag = 3;
            else if(matchflag && !matematchflag) matchflag = 1;
            else if(!matchflag && matematchflag) matchflag = 2;

            se_destructMatches(space, &read); 
          }

          if(slist2 && se_kdMatchListhasMatches(slist2) && 
              (!slist || !se_kdMatchListhasMatches(slist))) { 
            /*spliced mate full query*/ 
//            if(rep)    fprintf(nfo->dev, "SPLICED MATE FULL QUERY\n");

            if (se_kdMatchListhasMatches(list)) {
//              if(rep)       fprintf(nfo->dev, "QUERY LIST HAS MATCHES\n");
              pairStatus = MATE_SPL_FULL_QUERY;
              
	      /* select best query to spliced mate */
	      if (list->n[0] + list->n[1] > 1){

                /* get first fragment */
                frag.a = frag.b = -1;
                for (u = 0; u < 2; u++){
                  for (i = 0; i < slist2->n[u]; i++){
                    if (slist2->matches[u][i].fragno == 0){
                      frag.a = u;
                      frag.b = i;
                    }
                  }
                }
                assert(frag.a != -1 && frag.b != -1);
		bestpairlist = se_kdFindBestMate(space, list, &slist2->matches[frag.a][frag.b], maxedist);
		bl_gmatchlistDestruct(space, list);
		list = bestpairlist;
	      }
              se_setMatches(space, &read, list, maxedist, nfo, rep);
              matchflag = reportMatch(space, &map, reads, nfo, pairStatus, 0);
              se_destructMatches(space, &read); 
              /*spliced query no mate*/
            } else {
              pairStatus = MATE_SPL_NO_QUERY;
            }

            se_setMatches(space, &read, slist2, matemaxedist, nfo, rep);
            matematchflag = reportMatch(space, &map, reads, nfo, pairStatus, 1);

            if(matematchflag && matchflag) matchflag = 3;
            else if(matchflag && !matematchflag) matchflag = 1;
            else if(!matchflag && matematchflag) matchflag = 2;

            se_destructMatches(space, &read); 
          }

          if(slist && se_kdMatchListhasMatches(slist) && slist2 && 
              se_kdMatchListhasMatches(slist2)) {
            /*both spliced*/
//            if(rep)     fprintf(nfo->dev, "BOTH SPLICED\n");
            pairStatus = PAIR_SPL;

            se_setMatches(space, &read, slist, maxedist, nfo, rep);
            matchflag = reportMatch(space, &map, reads, nfo, pairStatus, 0);
            se_destructMatches(space, &read); 

            se_setMatches(space, &read, slist2, matemaxedist, nfo, rep);
            matematchflag = reportMatch(space, &map, reads, nfo, pairStatus, 1);
            se_destructMatches(space, &read); 

            if(matematchflag && matchflag) matchflag = 3;
            else if(matchflag && !matematchflag) matchflag = 1;
            else if(!matchflag && matematchflag) matchflag = 2;

          }

          if((!slist || !se_kdMatchListhasMatches(slist)) && se_kdMatchListhasMatches(matelist)) {
            pairStatus = MATE;

//            if(rep )     fprintf(nfo->dev, "ONLY MATE matemaxedist %d\n", matemaxedist);
            se_setMatches(space, &read, matelist, matemaxedist, nfo, rep);
            matematchflag = reportMatch(space, &map, reads, nfo, pairStatus, 1);
            se_destructMatches(space, &read);
            if(matematchflag) matchflag = 2; else matchflag = 0;
          }

          if((!slist2 || !se_kdMatchListhasMatches(slist2)) && se_kdMatchListhasMatches(list)) {
            pairStatus = QUERY;

//            if(rep )     fprintf(nfo->dev, "ONLY QUERY\n");
            se_setMatches(space, &read, list, maxedist, nfo, rep);
            matchflag = reportMatch(space, &map, reads, nfo, pairStatus, 0);
            se_destructMatches(space, &read); 
            if(matchflag) matchflag = 1; else matchflag = 0;
          }


          if(slist)  bl_gmatchlistDestruct(space, slist);
          if(slist2) bl_gmatchlistDestruct(space, slist2);

        } else {

//          if(rep )     fprintf(nfo->dev, "DISJOINT\n");
          matchflag = 0;
          matematchflag = 0;

          if(list && se_kdMatchListhasMatches(list)) {
            // matchflag = 1;
            se_setMatches(space, &read, list, maxedist, nfo, rep);
            matchflag = reportMatch(space, &map, reads, nfo, pairStatus, 0);
            se_destructMatches(space, &read); 
          }

          if(matelist && se_kdMatchListhasMatches(matelist)) {
            // matematchflag = 1;
            se_setMatches(space, &read, matelist, matemaxedist, nfo, rep);
            matematchflag = reportMatch(space, &map, reads, nfo, pairStatus, 1);
            se_destructMatches(space, &read); 

          }

          if(matematchflag && matchflag) matchflag = 3;
          else if(matchflag && !matematchflag) matchflag = 1;
          else if(!matchflag && matematchflag) matchflag = 2;


        }
      }

      bl_kdMatchstemDestruct(space, stems[0], len);
      bl_kdMatchstemDestruct(space, stems[1], len);
      if(matestems[0]) {
        bl_kdMatchstemDestruct(space, matestems[0], matelen);
        bl_kdMatchstemDestruct(space, matestems[1], matelen); 
        matestems[0] = NULL;
        matestems[1] = NULL;
      }

      bl_gmatchlistDestruct(space, list);
      if (nfo->bisulfite){
        FREEMEMORY(space, seqs[0]);
      }
      FREEMEMORY(space, seqs[1]);  

      if (bl_fastaHasMate(reads)) { 
        if (nfo->bisulfite){
          FREEMEMORY(space, mateseqs[0]);
        }
        FREEMEMORY(space, mateseqs[1]);
      }
      if (matelist) {
        bl_gmatchlistDestruct(space, matelist);
        matelist = NULL;
      }
    }

    bl_kdReportUnmatched(space, reads, k, matchflag, matematchflag, &best, &bestmate, nfo);
  }

  wrapBitmatrix(space, D, 2*(dim+1));
  FREEMEMORY(space, D);
  FREEMEMORY(space, enctab);
  FREEMEMORY(space, events);
  return;
}



/*--------------------------- bl_kdReportUnmatched ---------------------------
 *    
 * @brief dump the unmatched sequences to a device
 * @author Steve Hoffmann 
 *   
 */

void
bl_kdReportUnmatched (void *space, fasta_t *reads, Uint k, 
    unsigned char matchflag, unsigned char matematchflag, 
    bestseed_t *best, bestseed_t *bestmate, 
    segemehl_t *nfo)
{
   
  if(nfo->nomatchdev) { 

    if (matchflag == 0 && !bl_fastaHasMate(reads)) { 
      if (nfo->threadno > 1) pthread_mutex_lock(nfo->mtx2);

      if (!bl_fastaHasQuality(reads)){
        fprintf(nfo->nomatchdev, ">%s ef:%d;if:%d %d:%d %d:%d:%d\n%s\n", 
            bl_fastaGetDescription(reads, k), 
            best->maxevalconstraint, best->maxintervalconstraint,
            best->readstart, best->len,
            best->refidx, best->refpos, best->refstrand,
            
            bl_fastaGetSequence(reads, k)); 
      } else {	   
        fprintf(nfo->nomatchdev, "@%s ef:%d;if:%d %d:%d %d:%d:%d\n%s\n+%s\n%s\n",		
            bl_fastaGetDescription(reads, k),  
            best->maxevalconstraint, best->maxintervalconstraint,

            best->readstart, best->len,
            best->refidx, best->refpos, best->refstrand, 
            
            bl_fastaGetSequence(reads, k),
            bl_fastaGetDescription(reads, k), bl_fastaGetQuality(reads, k));
      }

      fflush(nfo->nomatchdev);
      if (nfo->threadno > 1) pthread_mutex_unlock(nfo->mtx2);
    }


    if ((matchflag < 3) && bl_fastaHasMate(reads)) {
      if (nfo->threadno > 1) pthread_mutex_lock(nfo->mtx2);


      if(matchflag == 0 || matchflag ==2) {
        if (!bl_fastaHasQuality(reads)){
          fprintf(nfo->nomatchdev, ">%s ef:%d;if:%d %d:%d %d:%d:%d\n%s\n", 
              bl_fastaGetDescription(reads, k), 
              best->maxevalconstraint, best->maxintervalconstraint,

              best->readstart, best->len,
              best->refidx, best->refpos, best->refstrand, 

              bl_fastaGetSequence(reads, k)); 
        } else {	   
          fprintf(nfo->nomatchdev, "@%s ef:%d;if:%d %d:%d %d:%d:%d\n%s\n+%s\n%s\n",		
              bl_fastaGetDescription(reads, k), 
              best->maxevalconstraint, best->maxintervalconstraint,
              
              best->readstart, best->len,
              best->refidx, best->refpos, best->refstrand,

              bl_fastaGetSequence(reads, k),
              bl_fastaGetDescription(reads, k), bl_fastaGetQuality(reads, k));
        }
      }

      if(matchflag == 0 || matchflag == 1) {
        if (!bl_fastaHasQuality(reads)){
          fprintf(nfo->nomatchdev, ">%s ef:%d;if:%d %d:%d %d:%d:%d\n%s\n", 
              bl_fastaGetMateDescription(reads, k),               
              bestmate->maxevalconstraint, bestmate->maxintervalconstraint,
              
              bestmate->readstart, bestmate->len,
              bestmate->refidx, bestmate->refpos, bestmate->refstrand, 
              
              bl_fastaGetMate(reads, k)); 
        } else {	   
          fprintf(nfo->nomatchdev, "@%s ef:%d;if:%d %d:%d %d:%d:%d\n%s\n+%s\n%s\n",		
              bl_fastaGetMateDescription(reads, k),               
              bestmate->maxevalconstraint, bestmate->maxintervalconstraint,

              bestmate->readstart, bestmate->len,
              bestmate->refidx, bestmate->refpos, bestmate->refstrand, 
              
              bl_fastaGetMate(reads, k),
              bl_fastaGetMateDescription(reads, k), bl_fastaGetMateQuality(reads, k));
        }
      }

      fflush(nfo->nomatchdev);
      if (nfo->threadno > 1) pthread_mutex_unlock(nfo->mtx2);
    }

  }

  return ;
}

