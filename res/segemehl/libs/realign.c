/*
 *  realign.c
 *  realigning
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 26.06.2012 12:55:18 CEST
 *   
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <pthread.h>
#include "alignment.h"
#include "debug.h"
#include "stringutils.h"
#include "basic-types.h"
#include "mathematics.h"
#include "manout.h"
#include "sort.h"
#include "matfile.h"
#include "bitVector.h"
#include "info.h"
#include "zran.h"
#include "sw.h"
#include "matchfiles.h"
#include "evalmatchfiles.h"
#include "manout.h"
#include "matchfilesfields.h"
#include "matepairs.h"
#include "matchfiles.h"
#include "manopt.h"
#include "iupac.h"
#include "info.h"
#include "fqueue.h"
#include "list.h"
#include "realign.h"
#include "multicharseq.h"
#include "charsequence.h"
#include "container.h"
#include "sw.h"
#include "segemehl.h"

pthread_mutex_t inout;
pthread_mutex_t mutkuh; /*?*/

/*------------------------------ destructBuffer ------------------------------
 *    
 * @brief to clear the realignmentlists' string buffersi
 * @author Steve Hoffmann 
 *   
 */
 
void
destructBuffer (void *elem)
{
 
  matchlistelem_t* myelem = (matchlistelem_t*) elem;
  FREEMEMORY(NULL, myelem->str);
	
  return ;
}

void*
bl_insertArray(void *arr, size_t nmemb, size_t size, void *elem, Uint pos) {
  void *dst, *src;

  arr = realloc(arr, size*(nmemb+1));
  if(pos>=nmemb) {
    dst = ((char*)arr)+(size*nmemb);
    memmove(dst, elem, size);
  } else {  
    src = ((char*)arr)+(size*pos);
    dst = ((char*)arr)+(size*(pos+1));
    memmove(dst, src, size*(nmemb-pos));
    memmove(src, elem, size);
  } 

  return arr;
}

void*
bl_removeArray(void *arr, size_t nmemb, size_t size, void *elem, Uint pos) {
  void *dst, *src;

  if(pos>=nmemb-1) {
    arr = realloc(arr, size*(nmemb-1));
  } else {  
    dst = ((char*)arr)+(size*pos);
    src = ((char*)arr)+(size*(pos+1));
    memmove(dst, src, size*(nmemb-pos));
    arr = realloc(arr, size*(nmemb-1));
  } 

  return arr;
}
/*------------------------------ cmp functions -------------------------------
 *    
 * @brief compare functions
 * @author Steve Hoffmann 
 *   
 */
 
int cmp_matchsplitsitecluster_qsort(const void *a, const void *b) {
  matchsplitsitecluster_t  *l = (matchsplitsitecluster_t*) a;
  matchsplitsitecluster_t  *r = (matchsplitsitecluster_t*) b;

  if(l->a < r->a) return -1;
  if(l->a > r->a) return 1;
  if(l->b < r->b) return -1;
  if(l->b > r->b) return 1;

  return 0;

}


Uint cmp_matchsplitsites_quickSort(Uint i, Uint j, void *arr, void* nfo) {
  matchsplitsite_t  *sites = (matchsplitsite_t*) arr;

  if(sites[i].cnt <  sites[j].cnt) return 1;
  if(sites[i].cnt > sites[j].cnt) return 2;
  
  return 0;

}

 Uint cmp_distsites_quickSort(Uint i, Uint j, void *arr, void* nfo) {
  matchsplitsearchkey_t *key = (matchsplitsearchkey_t*) nfo;
  matchlistelem_t  *sites = (matchlistelem_t*) arr;
  Uint pos1, pos2;

  if(sites[i].distchr >  sites[j].distchr) return 1;
  if(sites[i].distchr < sites[j].distchr) return 2;
  
  if(sites[i].distpos > sites[j].distpos) return 1; 
  if(sites[i].distpos <  sites[j].distpos) return 2;

  /*key->trans = left*/
  if(key->trans) {
    pos1 = sites[i].start - key->pos;
    pos2 = sites[j].start - key->pos;
  } else {
    pos1 = sites[i].end - key->pos;
    pos2 = sites[j].end - key->pos;
  }

  if(pos1 > pos2) return 1;
  if(pos1 < pos2) return 2;

  if(sites[i].trans > sites[j].trans) return 2;
  if(sites[i].trans < sites[j].trans) return 1;

  return 0;
}



Uint cmp_matchsplitsitecluster_bin(Uint a, void *data, void *key, void *nfo) {
	matchsplitsitecluster_t *d = (matchsplitsitecluster_t*) data;
	Uint *k = (Uint*) key;

	if (d[a].a > *k) return 1;
	if (d[a].b < *k) return 2;

	return 0;
 }


/*---------------------- bl_matchJoinSplitSiteClusters -----------------------
 *    
 * @brief 
 * @author Steve Hoffmann 
 *   
 */
 
matchsplitsitecluster_t *
bl_matchJoinSplitSiteClusters (void *space, matchsplitsitecluster_t *dst, 
    matchsplitsitecluster_t *src, Uint dstelems, Uint srcelems)
{
  dst = ALLOCMEMORY(space, dst, matchsplitsitecluster_t, dstelems+srcelems);
  memmove(&dst[dstelems], src, sizeof(matchsplitsitecluster_t)*(srcelems));

  return dst;
}

/*------------------------ bl_matchfileDestructMatch -------------------------
 *    
 * @brief destroy a match
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchfileDestructMatchListElem (void *elem)
{
  matchlistelem_t *match = (matchlistelem_t*) elem;
  FREEMEMORY(NULL, match->bookkeeper);
  FREEMEMORY(NULL, match->str);

  return ;
}
/*------------------------- bl_matchfileEnlistMatch --------------------------
 *    
 * @brief enlist match alignment to unsorted list
 * @author Steve Hoffmann 
 *   
 */
 
List *
bl_matchfileEnlistMatch (void *space, List *l, Uint start, Uint end, 
    Uint distchr, Uint distpos, Uint adjointpos, char trans, char *str, 
			 unsigned char* bookkeeper)
{
 
  matchlistelem_t elem;
  elem.start = start;
  elem.end = end;
  elem.str = str;
  elem.distchr = distchr;
  elem.distpos = distpos;
  elem.trans = trans;
  elem.adjoint = adjointpos;
  elem.bookkeeper = bookkeeper;
 
  bl_listInsert(l, l->last, &elem);
  return l;
}

/*------------------------- bl_matchfileUnlistMatch -------------------------
 *    
 * @brief unlist a match from unsorted list if it is mapped to [a,b]
 * @author Steve Hoffmann 
 *   
 */

  matchlistelem_t*
bl_matchfileUnlistMatches (void *space, List *l, Uint a, Uint b, Uint *k, 
    char left)
{
  matchlistelem_t *elem, *arr=NULL;
  Uint n=0;
  Uint pos=0;
  int cur, next;

  cur = l->first;
  while(cur != -1) { 
    next = l->nodes[cur].next;
    elem = (matchlistelem_t*) bl_listGetElem (l, cur);  
    pos = (left) ? elem->start : elem->end;
    if(elem && pos >= a && pos <=b) {
      elem = (matchlistelem_t*) bl_listUnlink(l, cur, NULL);
      arr = ALLOCMEMORY(space, arr, matchlistelem_t, n+1);
      memmove(&arr[n], elem, sizeof(matchlistelem_t));
      FREEMEMORY(space, elem);
      n++;
    }
    
    cur = next;
  }

  *k = n;
  return arr;
}


/*------------------------ bl_matchfileClusterSplits -------------------------
 *    
 * @brief cluster the split sites w/ greedy lr extension within r-l+1=interval
 * @author Steve Hoffmann 
 *   
 */

  matchsplitsitecluster_t*
bl_matchfileClusterSplits (void *space, matchsplitsite_t *sites, Uint n, 
    Uint interval, Uint *noofclusters)
{
  Uint *idx;
  char *mrk;
  matchsplitsitecluster_t *C=NULL;
  Uint k=0, i, j, unmarked = n;
  Uint med, mid, sum=0, *cum=NULL;
  
  mrk = ALLOCMEMORY(space, NULL, char, n);
  memset(mrk, 0, n);
  idx = quickSort(space, sites, n, cmp_matchsplitsites_quickSort, NULL);

  while(unmarked > 0) { 
    C = ALLOCMEMORY(space, C, matchsplitsitecluster_t, k+1);
    C[k].a =0;
    C[k].b =0;
    C[k].cnt = 0;

#ifdef CHECKLINKS
    C[k].trans = 0;
    C[k].lnkcnt = 0;
#endif

    sum =0;
    for(i=0; i < n; i++) {
      if(!mrk[idx[i]]) { 
        //left expand
        if(C[k].a == 0 ||
            (C[k].a > sites[idx[i]].pos && 
             sites[idx[i]].pos+interval >= C[k].b)) {
          C[k].a = sites[idx[i]].pos;
          mrk[idx[i]]=1;
        }
        //right expand
        if(C[k].b == 0 ||
            (C[k].b < sites[idx[i]].pos && 
             sites[idx[i]].pos <= C[k].a+interval)) {
          C[k].b = sites[idx[i]].pos;
          mrk[idx[i]]=1;
        }
        //inclusion
        if(C[k].a <= sites[idx[i]].pos && C[k].b >= sites[idx[i]].pos) {

          for(j=0; j < sites[idx[i]].cnt; j++) { 
          cum = bl_insertArray(cum, sum, sizeof(Uint), 
              &sites[idx[i]].pos, sum);
          sum++;
          }

          C[k].cnt += sites[idx[i]].cnt;
#ifdef CHECKLINKS
          C[k].trans += sites[idx[i]].trans;
#endif
          mrk[idx[i]]=1;
        }

        if(mrk[idx[i]]) {
          unmarked--;
        }
      }
    }
  
    med = C[k].a;
    if(sum) {
      qsort(cum, sum, sizeof(Uint), cmp_Uint_qsort);
      mid = sum/2 + sum%2;
      med = cum[mid-1];
      FREEMEMORY(space, cum);
    }

    assert(med >= C[k].a && C[k].b >= med);
    C[k].median = med;
     /*  C[k].median*/
    k++;
  }

  FREEMEMORY(space, idx);
  FREEMEMORY(space, mrk);
  *noofclusters = k;

  return C;
}


/*----------------------------- bl_matchDistPos ------------------------------
 *    
 * @brief get the condensed dist pos list
 * @author Steve Hoffmann 
 *   
 */

  void
bl_matchsplitdistpos (void *space, matchsplitsitecluster_t *C, 
    matchlistelem_t *arr, Uint n, char left)  
{

  Uint i, j;

#if defined PARANOID && defined CHECKLINKS
  Uint pos;
#endif

  matchlistelem_t *elem;
  Uint *idx;
  Uint noofdists=0;
  int16_t *adjoint = NULL;
  Uint inelem =0;
  int16_t inelem1 =0;
  Uint *adjointcnt = NULL;
  int16_t relpos  =0;
  Uint noofadjoints = 0;

  matchsplitsearchkey_t key;

  key.trans = left;
  key.pos = C->a;

  /*generate sort idx to for unsorted list arr*/
  idx = quickSort(space, arr, n, cmp_distsites_quickSort, &key);

  C->distchr = ALLOCMEMORY(space, NULL, Uint, noofdists+1);
  C->distpos = ALLOCMEMORY(space, NULL, Uint, noofdists+1);
  C->disttrans = ALLOCMEMORY(space, NULL, char, noofdists+1);
  C->distcnt = ALLOCMEMORY(space, NULL, Uint, noofdists+1);
  
#if defined PARANOID && defined CHECKLINKS
  C->rpos = ALLOCMEMORY(space, NULL, char, noofdists+1);
  pos = (left) ? arr[idx[0]].start - C->a : arr[idx[0]].end - C->a;
  C->rpos[0] = pos;
#endif

  C->distchr[0] = arr[idx[0]].distchr;
  C->distpos[0] = arr[idx[0]].distpos;
  C->disttrans[0] = arr[idx[0]].trans;
  C->distcnt[0] = 1;
  noofdists++;

  /*iter arr using the sort index idx*/
  for(i=1; i < n; i++) {
    elem = &arr[idx[i]];

#if defined PARANOID && defined CHECKLINKS
    pos = (left) ? arr[idx[i]].start - C->a : arr[idx[i]].end - C->a;
#endif   
    if( elem->distchr != C->distchr[noofdists-1] || 
        elem->distpos != C->distpos[noofdists-1] ||
#if defined PARANOID && defined CHECKLINKS
                  pos != C->rpos[noofdists-1]    || 
#endif
          elem->trans != C->disttrans[noofdists-1]
      ) {

      C->distchr = ALLOCMEMORY(space, C->distchr, Uint, noofdists+1);
      C->distpos = ALLOCMEMORY(space, C->distpos, Uint, noofdists+1);
      C->disttrans = ALLOCMEMORY(space, C->disttrans, char, noofdists+1);
      C->distcnt = ALLOCMEMORY(space, C->distcnt, Uint, noofdists+1);
#if defined PARANOID && defined CHECKLINKS
      C->rpos = ALLOCMEMORY(space, C->rpos, char, noofdists+1);
      C->rpos[noofdists] = pos;
#endif
      C->distchr[noofdists] = elem->distchr;
      C->distpos[noofdists] = elem->distpos;
      C->disttrans[noofdists] = elem->trans;
      C->distcnt[noofdists] = 1;
      noofdists++;
    } else { 
      C->distcnt[noofdists-1] += 1;
    }
    
    /*keep a sorted list for adjoints and their counts*/
    if(elem->adjoint != -1) {
      relpos = elem->adjoint - C->median;
      for(j=0; j < noofadjoints; j++) {
        if(adjoint[j] == relpos) {
            adjointcnt[j]++;
            break;
        }
        if(adjoint[j] > relpos) {
          break;
        }
      }

      if(noofadjoints==j || adjoint[j] > relpos) {
        inelem1 = relpos;
        adjoint = bl_insertArray(adjoint, noofadjoints, 
            sizeof(int16_t), &inelem1, j);
        inelem = 1;
        adjointcnt = bl_insertArray(adjointcnt, noofadjoints, sizeof(Uint), 
            &inelem, j);
        noofadjoints++;
      }
    }
  }


  C->distsites = noofdists;
  C->adjoint = adjoint;
  C->adjointcnt = adjointcnt;
  C->noofadjoints = noofadjoints;
  C->noofadjclust = 0;
  C->adjclust = NULL;
  C->adjclustweights = NULL;
  C->distclust = NULL;
  C->distclustchr = NULL;
  C->distclustcnt = NULL;
  C->distclusttype = NULL;
  C->realigned = NULL;
  C->emat = NULL;
  C->distclustrealigned = NULL;
  C->noofdistclust = 0;

#if defined PARANOID && defined CHECKLINKS
  C->linked = calloc(noofdists, sizeof(Uint));
#endif

  FREEMEMORY(space, idx);

  return ;
}

/*------------------------ bl_matchfileScanRightList ------------------------
 *    
 * @brief scan the list to accumulate splits
 * @author Steve Hoffmann 
 *   
 */

  matchsplitsitecluster_t *
bl_matchfileScanList (void *space, List *l, matchsplitsite_t *sites, 
    Uint n, Uint interval, Uint curstart, char left, Uint *noofclusters)
{

  Uint i, nstr, k;
  matchlistelem_t* arr=NULL;
  matchsplitsitecluster_t *C;

  /*first determine the cluster intervals and the medians in sites array*/
  C = bl_matchfileClusterSplits (space, sites, n, interval, &k);

  for(i=0; i < k; i++) {
    /*then unlist the matches in unsorted list*/
    arr = bl_matchfileUnlistMatches (space, l, C[i].a, C[i].b, &nstr, left);
    /*add the elems for unsorted list to cluster*/
    bl_matchsplitdistpos (space, &C[i], arr, nstr, left) ;

    if(nstr >0) FREEMEMORY(space, arr);
  }

#ifdef CHECKLINKS
  Uint j;
  for(i=0; i < k; i++) {
    nstr = 0;
    for(j=0; j < C[i].distsites; j++) {
      if(C[i].disttrans[j])
        nstr += C[i].distcnt[j];
    }
    assert(nstr == C[i].trans);
  }
#endif
  

  qsort(C, k, sizeof(matchsplitsitecluster_t), cmp_matchsplitsitecluster_qsort);

  *noofclusters = k;
  return C;
}

/*------------------------ bl_matchfileUpdateSplitSites------------------------
 *    
 * @brief update and sort the lsites and rsites arrays. 
 * @author Steve Hoffmann 
 *   
 */

  matchsplitsite_t*
bl_matchfileUpdateSplitSites (void *space, Uint e,  
    matchsplitsite_t *sites, Uint *noofsites, Uint interval, char trans)
{
  Uint i, n;
  matchsplitsite_t elem;
  n = *noofsites;

  for(i=0; i < n; i++) { 
    if(sites[i].pos == e) {
      sites[i].cnt++;
      sites[i].trans+=trans;
    }
    if(sites[i].pos >= e) {
      break;
    }
  }
    
  if(i==n || sites[i].pos > e) {
    elem.pos = e;
    elem.cnt = 1;
    elem.trans = trans;
    sites = bl_insertArray(sites, n, sizeof(matchsplitsite_t), &elem, i);
    n++;
  }

  *noofsites = n;
  return sites;
}


/*------------------------- bl_matchfileGetOverhang --------------------------
 *    
 * @brief get the overhanging sequence, boundaries and score
 * @author Steve Hoffmann 
 *   
 */

int
bl_matchfileSplitAlignment (char *read, char *ref, char *aln, 
    Uint allen, Uint splitpos, int *scores, 
    Alignment **leftaln, Alignment **rightaln, char left)
{ 
  int xscr =0, yscr=0;
  Uint  p=0, q=0, k=0, xqlen=0, xrlen=0, yqlen=0, yrlen =0, yendS=0;
  Alignment *xal, *yal;
  Uint   x0; 
   
  xal = ALLOCMEMORY(space, NULL, Alignment, 1);
  initAlignment(xal, read, strlen(read), 0, ref, strlen(read), 0);
  yal = ALLOCMEMORY(space, NULL, Alignment, 1);
  initAlignment(yal, read, strlen(read), 0, ref, strlen(read), 0);

  //lookahead
  while(aln[k] == 'S' || aln[k] =='H') {
    if(aln[k] != 'H') q++; 
    k++;
  }

  x0 = q;

  for(; k < allen; k++) {
    switch(aln[k]) {
      case 'M':
      case 'X': 
      case '=':
        if(p <= splitpos) {
          xscr += scores[!matchIUPAC(ref[p],read[q])];/*[(ref[p]!=read[q])];*/
           insertEop(xal, Replacement);
          xrlen++; xqlen++;
            } else {
            yscr += scores[!matchIUPAC(ref[p],read[q])];/*[(ref[p]!=read[q])];*/
          insertEop(yal, Replacement);
          yrlen++; yqlen++;
        }
        p++; q++;
        break;
      case 'D':
          if(p <= splitpos) {
          xscr += scores[1];
           insertEop(xal, Deletion);
          xrlen++;
         } else {
           yscr += scores[1];
          insertEop(yal, Deletion);
          yrlen++; 
        }
        p++;
        break;
      case 'I':
        //left
        if(p <= splitpos+left) {
          xscr += scores[1];
          insertEop(xal, Insertion);
          xqlen++;
        } else {
          yscr += scores[1];
          insertEop(yal, Insertion);
          yqlen++; 
        }
        q++;
        break;
    case 'S':
      if(p <= splitpos) {
      }
      else {
	yendS++;
      }
      break;
    case 'N':
    case 'P':
          p++;
      default:
        break;
    }
  }

  xal->u = read;
  xal->uoff = x0;
  xal->ulen = xqlen;
  
  xal->v = ref;
  xal->voff= 0;
  xal->vlen= xrlen;
  
  yal->v = read;
  yal->uoff = strlen(read)-(yqlen+yendS);
  yal->ulen = yqlen;

  yal->v = ref;
  yal->voff = splitpos+1;
  yal->vlen= yrlen;
 
  *leftaln  = xal;
  *rightaln = yal;
 
  assert(getAlignScore(xal, scores, -1) == xscr);
  assert(getAlignScore(yal, scores, -1) == yscr);

  return q;
}


/*------------------------ bl_matchfileRealignWriteSAM -----------------------
 *    
 * @brief write the realignment in SAM
 * @author Steve Hoffmann 
 *   
 */

  char*
bl_matchfileRealignWriteSAM (matchfileRec_t *orig, char *qual,  
    char *rname, Uint pos, char strand, Alignment* al, Uint ustart, Uint uend,  
    Uint uno, Uint off, char *rnext, Uint pnext, char fnext)
{


  char *tag = NULL, *md, *cigar, flg=0, *newqual, *newseq;
  Uint len=0,  ptr=0, flag, allen =0,   edist=0,
       noofsplits=0, resplits=0;
  int addx=0;


  flag = orig->flag;
  flag |= (1 << 8);
  edist = getEdist(al);
  allen  = getUalignlen(al);

  if(strand == '-') { 
    flag |= (1 << 4);
  } else {
    flag &= ~(1 << 4);
  }
  
  newqual = ALLOCMEMORY(space, NULL, char, strlen(qual)+1); 

  memmove(newqual, qual, strlen(qual));
  if(orig->strand != strand) {
    newqual = strrev(newqual, strlen(qual));
  }

  newseq = ALLOCMEMORY(space, NULL, char, allen+1);  
  memmove(newseq, &al->u[al->uoff], allen);
  memmove(newqual, &newqual[al->uoff], allen);
  newseq[allen] = 0;
  newqual[allen] = 0;


  noofsplits = orig->noofsplits;
  resplits = 2;

  cigar = cigarstring(al, 0, 0, 'S', 0);
  
  md = mdstring(al,0);

  len = snprintf(NULL, 0, "%s\t%d\t%s\t%d\t255\t%s\t*\t0\t0\t%s\t%s\t", 
      orig->curname, flag, rname, pos, cigar, newseq, newqual);
  tag = ALLOCMEMORY(space, tag, char, ptr+len+1);
  snprintf(&tag[ptr], len+1, "%s\t%d\t%s\t%d\t255\t%s\t*\t0\t0\t%s\t%s\t", 
      orig->curname, flag, rname, pos, cigar, newseq, newqual);
  ptr += len; 

  if(orig->noofsplits) { 
    len = snprintf(NULL, 0, 
        "NM:i:%d\tMD:Z:%s\tNH:i:%d\tXI:i:%d\tXL:i:%d\tXR:i:%d\tXO:i:%d\t", 
		   edist, md, orig->curcnt, orig->identity,noofsplits, resplits, orig->xno);

    tag = ALLOCMEMORY(space, tag, char, ptr+len+1);
    snprintf(&tag[ptr], len+1,
        "NM:i:%d\tMD:Z:%s\tNH:i:%d\tXI:i:%d\tXL:i:%d\tXR:i:%d\tXO:i:%d\t", 
        edist, md, orig->curcnt, orig->identity, noofsplits, resplits, orig->xno);
    ptr += len;

  } else {
    len = snprintf(NULL, 0, 
        "NM:i:%d\tMD:Z:%s\tNH:i:%d\tXI:i:%d\tXR:i:%d\t", 
        edist, md, orig->curcnt, orig->identity, resplits);

    tag = ALLOCMEMORY(space, tag, char, ptr+len+1);    
    snprintf(&tag[ptr], len+1,
        "NM:i:%d\tMD:Z:%s\tNH:i:%d\tXI:i:%d\tXR:i:%d\t", 
        edist, md, orig->curcnt, orig->identity, resplits);
    ptr += len;

  }
  if(orig->xstart>0){
    addx=orig->xstart-1;
  }
  len = snprintf(NULL, 0, "XX:i:%d\tXY:i:%d\tXQ:i:%d\t", ustart+addx, uend+addx, uno);
  tag = ALLOCMEMORY(space, tag, char, ptr+len+1);
  snprintf(&tag[ptr], len+1, "XX:i:%d\tXY:i:%d\tXQ:i:%d\t", ustart+addx, uend+addx, uno);	  
  ptr += len;

  if(uno == 0) { /*??*/ 
    if(orig->donorchr) { 
      len = snprintf(NULL, 0, "XP:Z:%s\tXU:i:%d\tXS:i:%d\t",orig->donorchr, 
		     orig->donorpos, orig->donorflg); /**/
      tag = ALLOCMEMORY(space, tag, char, ptr+len+1);
      snprintf(&tag[ptr], len+1,"XP:Z:%s\tXU:i:%d\tXS:i:%d\t",orig->donorchr, 
          orig->donorpos, orig->donorflg); /**/
      ptr += len;
    }

    flg = (fnext == '-') ? 0 : SPLIT_NEXT_PLUS;

    len = snprintf(NULL, 0, "XC:Z:%s\tXV:i:%d\tXT:i:%d", rnext, pnext, flg);      /**/
    tag = ALLOCMEMORY(space, tag, char, ptr+len+1);
    snprintf(&tag[ptr], len+1,"XC:Z:%s\tXV:i:%d\tXT:i:%d",rnext, pnext, flg);  /**/
    ptr += len;

  } else {
    flg = (fnext == '-') ? 0 : SPLIT_PREV_PLUS;
    len = snprintf(NULL, 0, "XP:Z:%s\tXU:i:%d\tXS:i:%d",rnext, pnext, flg);  /**/
    tag = ALLOCMEMORY(space, tag, char, ptr+len+1);
    snprintf(&tag[ptr], len+1,"XP:Z:%s\tXU:i:%d\tXS:i:%d",rnext, pnext, flg);  /**/
    ptr += len;

    if(orig->acceptorchr) { 
      len = snprintf(NULL, 0, "\tXC:Z:%s\tXV:i:%d\tXT:i:%d",orig->acceptorchr, 
          orig->acceptorpos, orig->acceptorflg); 
      tag = ALLOCMEMORY(space, tag, char, ptr+len+1);
      snprintf(&tag[ptr], len+1,"\tXC:Z:%s\tXV:i:%d\tXT:i:%d",orig->acceptorchr, 
          orig->acceptorpos, orig->acceptorflg);
      ptr += len;
    }
  }


  FREEMEMORY(space, newseq);
  FREEMEMORY(space, newqual);
  FREEMEMORY(space, cigar);
  FREEMEMORY(space, md);


  return tag;
}

/*------------------------- bl_matchfileRealign -------------------------
 *    
 * @brief realign right and left sides of alignments
 *
 *              median
 *   -------------|-----> |
 *   <------------|------ |
 *                      median + interval
 *
 *
 *  realign left sides of alignments
 *                  median
 *            | <-----|------------
 *  median+interval
 *
 *
 *
 * @author Steve Hoffmann 
 *   
 */

  Uint
bl_matchfileRealign (void *space, FILE *realigndev, List *rlist, matchsplitsitecluster_t *T, 
		     Uint n, char left, fasta_t *set, Uint chromidx, unsigned char fmt, Uint interval, int threadno, Uint sinterval,int MAXSPLICE)
{

  Uint  i, q,  nstr; /* *siteidx,*/
  matchlistelem_t* arr=NULL;
  Uint noofrealigns=0;
  pthread_t *threads;
  readthread *ReadT;

  /*this would be the point if we boss it*/
  for(i=0; i < n; i++) {
    if(T[i].noofdistclust<MAXSPLICE) {
    if(left) {
      arr = bl_matchfileUnlistMatches(space, rlist, T[i].median-interval, T[i].median, 
          &nstr, 1);  
    } else { 
      arr = bl_matchfileUnlistMatches(space, rlist, T[i].median, T[i].median+interval, 
          &nstr, 0);
    }
    //    fprintf(stderr,"%d,%d\n",T[i].distsites,nstr); 
    if(threadno<2||nstr<2) {
      Readrealign(0, nstr-1, arr, fmt, &(T[i]),left,set,chromidx, &noofrealigns,realigndev, rlist,interval,space,sinterval);
    }
    else {
      ReadT=ALLOCMEMORY(space, NULL, realignthread , threadno);/*N_threads??*/
      threads = ALLOCMEMORY(space, NULL, pthread_t, threadno);  
      ReadT[0].T=&(T[i]);;
      ReadT[0].arr=arr;
      ReadT[0].fmt=fmt;
      ReadT[0].left=left;
      ReadT[0].set=set;
      ReadT[0].chromidx=chromidx;
      ReadT[0].noofrealigns=0;
      ReadT[0].realigndev=realigndev;
      ReadT[0].rlist=rlist;
      ReadT[0].interval=interval;
      ReadT[0].space=space;
      ReadT[0].sinterval=sinterval;
      //pthread_create(&threads[0]/*???*/, NULL, Readthreadstarter, &RT[0]);
      if(nstr<threadno) {
	ReadT[0].begin=0;
	ReadT[0].stop=0;	
	pthread_create(&threads[0]/*???*/, NULL, Readthreadstarter, &ReadT[0]);
     
	for(q=1; q < nstr; q++) {
	  memmove(&ReadT[q], &ReadT[0], sizeof(readthread));
	  ReadT[q].begin=q;
	  ReadT[q].stop=q; 
	  pthread_create(&threads[q]/*???*/, NULL, Readthreadstarter, &ReadT[q]);
	}
	for(q=0; q < nstr; q++) {
	  pthread_join(threads[q], NULL); 
	}
	for(q=0; q < nstr; q++) {
	  noofrealigns+=ReadT[q].noofrealigns;
	}
	
      }
      else {
	int chunk;
	int overhead;
	int where=0;
	chunk=nstr/threadno;
	overhead=nstr%threadno;
	chunk--;
	if (where<overhead) {
	  q=1+chunk;
	}
	else {
	  q=chunk;
	}
	ReadT[where].begin=0;
	ReadT[where].stop=q;
	assert(q<nstr);
	pthread_create(&threads[where]/*???*/, NULL, Readthreadstarter, &ReadT[where]);
	for(where=1; where<threadno; where++) {
	  q++;
	  memmove(&ReadT[where], &ReadT[0], sizeof(readthread));
	  ReadT[where].begin=q;
	  if (where<overhead) {
	    q++;
	  }
	  q+=chunk;
	  ReadT[where].stop=q;
	  assert(q<nstr);
	  pthread_create(&threads[where]/*???*/, NULL, Readthreadstarter, &ReadT[where]);
	}  
	for(q=0; q < threadno; q++) {
	  pthread_join(threads[q], NULL); 
	}
	for(q=0; q < threadno; q++) {
	  noofrealigns+=ReadT[q].noofrealigns;
	}
      }
      FREEMEMORY(space, threads);
      FREEMEMORY(space, ReadT);
    }

    FREEMEMORY(space, arr);
  }
}
  return noofrealigns;
}


 
/*----------------------- bl_matchfileDumpUnRealigned ------------------------
 *    
 * @brief dump the un-realigned reads to device and clean list
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchfileDumpUnRealigned (void *space, FILE* realigndev, List *llist, List *rlist, 
    Uint curstart)
{
  Uint cur, next;
  matchlistelem_t* elem=NULL;
  cur = llist->first;

  while(cur != -1) { 
    next = llist->nodes[cur].next;
    elem = (matchlistelem_t*) bl_listGetElem (llist, cur);  
 
    if(elem->end < curstart) { 
   
      elem = (matchlistelem_t*) bl_listUnlink(llist, cur, NULL);
      if(!(elem->bookkeeper[0] & LEFTSPLIT) && 
          !(elem->bookkeeper[0] & RIGHTSPLIT)) {
        elem->bookkeeper[0] |= LEFTSPLIT;
        fprintf(realigndev, "%s\n", elem->str);
      }
      
      if(!(elem->bookkeeper[0] & RIGHTLIST)) 
        FREEMEMORY(space, elem->bookkeeper);

      FREEMEMORY(space, elem->str);
      FREEMEMORY(space, elem);
    }
    cur = next;
  }

  cur = rlist->first;
  while(cur != -1) { 
    next = rlist->nodes[cur].next;
    elem = (matchlistelem_t*) bl_listGetElem (rlist, cur);  
 
    if(elem->end < curstart) { 
   
      elem = (matchlistelem_t*) bl_listUnlink(rlist, cur, NULL);
      if(!(elem->bookkeeper[0] & LEFTSPLIT) && 
          !(elem->bookkeeper[0] & RIGHTSPLIT)) {
        elem->bookkeeper[0] |= RIGHTSPLIT;
        fprintf(realigndev, "%s\n", elem->str);
      }
      
      FREEMEMORY(space, elem->bookkeeper);
      FREEMEMORY(space, elem->str);
      FREEMEMORY(space, elem);
    }
    cur = next;
  }

   return ;
}

/*----------------------- bl_matchfileRealignScanFile ------------------------
 *    
 * @brief read all matches from start to end on chromname
 * @author Steve Hoffmann 
 *   
 */

Uint
bl_matchfileRealignScanFileNew(void *space, matchfile_t *file, FILE *realigndev,
    fasta_t *set, unsigned char fields, matchsplitsiteclusterlist_t **Lclust,   
			       matchsplitsiteclusterlist_t **Rclust, Uint *nchr, int threadno, int maxdist) {

  stringset_t *token;

  Uint buffersize=1024, //startbin, //endbin, 
       len=0, k=-1, 
       curstart=0, curend=0, 
    acceptorpos=0, donorpos=0, /*xstart, xend, xno,*/  
       acceptorchridx, donorchridx, 
       acceptorflg = 0, donorflg = 0, curchromidx=0, adjoint=0
         //,counter=0, pnext = 0, curcnt = 0
         ; 
  char *buffer = NULL, *buffer2 =  NULL, *buffer3 = NULL, 
       ch,
       //*curseq=NULL, *curqual=NULL, 
       //*curaln, 
       *filename, strand,
       *acceptorchr = NULL, *donorchr = NULL, 
       //*rnext, 
       *curchrom;
  unsigned char header = 1;
    //readlen=0, u, 
    /*int  edist=0;*/

  matchfileindex_t *index;
  unsigned char gzip, fmt;/*, curedist=0;*/
  FILE *fp = NULL;
  struct gzidxfile *gzf = NULL;
  matchsplitsitecluster_t *T=NULL;
  matchsplitsiteclusterlist_t *L=NULL, *R=NULL;
  int gzlen;
  Uint temp = 0;
  List rightl, leftl, rightrealign, leftrealign;  
  Uint interval=10;
  matchsplitsite_t *rsites=NULL;
  matchsplitsite_t *lsites=NULL;
  Uint nooflsites=0;
  Uint noofrsites=0;
  Uint trans = 0;
  Uint noofleftrealigns = 0, noofrightrealigns = 0;
  unsigned char *bookkeeper=NULL;

  gzip = file->gzip;
  fmt = file->fmt; 
  filename = file->filename;
  index = file->index;
  buffer = ALLOCMEMORY(space, NULL, char, buffersize);

  if (gzip) {  
    index->gzindex = bl_zranGetIndex(filename, &gzlen);
    fp = fopen(filename, "rb");
    gzf = bl_initgzidxfile(fp, index->gzindex, 0, CHUNK);
  } else {
    fp = fopen(filename, "r");
  }

  if(fp == NULL) {
    DBGEXIT("Couldn't open file %s. Exit forced!\n", filename);
  }

  //initialize buffer list
  bl_listInit(&rightl, 100, sizeof(matchlistelem_t));
  bl_listInit(&leftl, 100, sizeof(matchlistelem_t));
  bl_listInit(&rightrealign, 100, sizeof(matchlistelem_t));
  bl_listInit(&leftrealign, 100, sizeof(matchlistelem_t));

  L = *Lclust;
  R = *Rclust;

  while((ch = (gzip) ? bl_getgzidxc(gzf) : getc(fp)) != EOF) {    

    if(len == buffersize-1) {
      buffersize = 2*buffersize+1;
      buffer = ALLOCMEMORY(space, buffer, char, buffersize);
    }

    if(ch == '\n' && len > 0) {

      buffer = ALLOCMEMORY(space, buffer, char, len+1); 
      buffer[len] = '\0';

#ifdef DBGIDX
      DBG("buffer: %s\n", buffer);
#endif

      if(header) header = bl_matchfileIsHeader(buffer, len, fmt);

      if(!header) { 
        token = tokensToStringset(space, "\t", buffer, len);

        curchrom = bl_matchfileGetChrom(token, fmt);
        curstart = bl_matchfileGetStartPos(token, fmt);
        curend   = bl_matchfileGetEndPos(token, fmt);
        curchromidx = bl_matchfileGetChromIndexNumber(index, curchrom);
        if(k == -1) {
          k = curchromidx;
         // fprintf(stderr, "chrom:%s %d\n", curchrom, curchromidx);
        }

        if(k != curchromidx) {
          //finalize stuff

          T = bl_matchfileScanList (space, &leftl, 
              lsites, nooflsites, interval, -1, 1, &temp);
          
          if(realigndev) { 
            noofleftrealigns += bl_matchfileRealign(space, realigndev, &leftrealign, T, temp, 1, set, 
						    k, fmt, 25,threadno, interval, maxdist);
            bl_listSweep(&leftrealign);
          }

          L[k].cluster = 
            bl_matchJoinSplitSiteClusters(space, L[k].cluster, T, L[k].noofclusters, temp);
          L[k].noofclusters += temp;

          FREEMEMORY(space, T);
          FREEMEMORY(space, lsites);
          nooflsites = 0;

          T = bl_matchfileScanList (space, &rightl, 
              rsites, noofrsites, interval, -1, 0, &temp);

          if(realigndev) { 
            noofrightrealigns += bl_matchfileRealign(space, realigndev, &rightrealign, T, temp, 0, set, k, fmt, 25,threadno,interval,maxdist);
            bl_listSweep(&rightrealign);
            bl_matchfileDumpUnRealigned(space, realigndev, &leftrealign, &rightrealign, -1);
          }

          R[k].cluster = 
            bl_matchJoinSplitSiteClusters(space, R[k].cluster, T, R[k].noofclusters, temp);
          R[k].noofclusters += temp;

          FREEMEMORY(space, T);
          FREEMEMORY(space, rsites);
          noofrsites = 0;

          k = curchromidx;
   
          bl_listDestruct(&rightrealign, NULL);
          bl_listDestruct(&leftrealign, NULL);
          bl_listInit(&rightrealign, 100, sizeof(matchlistelem_t));
          bl_listInit(&leftrealign, 100, sizeof(matchlistelem_t));

         // fprintf(stderr, "chrom:%s %d\n", curchrom, curchromidx);
        }


        //fprintf(stderr, "reading %s: %d - %d\n", curchrom, curstart, curend);
        /*last condition to avoid inclusion of 0-alignments in BAM files*/
        if (curstart != curend && curend+1 > 0) {

	  /*  edist    = bl_matchfileGetEdist(token, fmt);*/
          strand   = bl_matchfileGetStrand(token, fmt);

          if(fields & MFREAD_SPLITS) { 
            acceptorpos = bl_matchfileGetNextPos(token, fmt);
            acceptorchr = bl_matchfileGetNextChr(token, fmt);
            acceptorflg = bl_matchfileGetNextFlag(token, fmt);
            donorpos = bl_matchfileGetPrevPos(token, fmt);
            donorchr = bl_matchfileGetPrevChr(token, fmt);
            donorflg = bl_matchfileGetPrevFlag(token, fmt);
            /*xstart = bl_matchfileGetSplitStart(token, fmt);
	      xend = bl_matchfileGetSplitEnd(token, fmt);
	      xno = bl_matchfileGetSplitNumber(token, fmt);*/
          }

	  /*          if(edist > 255) {
            curedist = 255; 
	    } else curedist = edist;*/

          donorchridx = -1;
          acceptorchridx = -1;


          /*
           * rsites (lsites) holds split positions that occur at 
           * rightmost (leftmost) 
           * sites of split alignments. Both arrays are sorted and the 
           * largest position is last. 
           * The actual alignment information is stored in two unsorted lists. 
           * They are scanned if the current alignment
           * start is at least interval nt larger then the largest element 
           * in lsites or rsites. 
           * In this case there is no further split alignment (sorted file)
           * that belongs to a cluster of splits in the list. The scan
           * returns an array of split site clusters (interval, median ...) 
           * which is appended to the L or R split site cluster array
           * 
           */

          /*
           * cleaning up
           *
           */

          if(lsites && lsites[nooflsites-1].pos+interval < curstart) {

            T = bl_matchfileScanList (space, &leftl, 
                lsites, nooflsites, interval, curstart, 1, &temp);

            if(realigndev) { 
              noofleftrealigns += bl_matchfileRealign(space, realigndev, &leftrealign, T, temp, 1, set, 
						      curchromidx, fmt, 25,threadno, interval, maxdist);
              bl_listSweep(&leftrealign);
            }

            L[k].cluster = 
              bl_matchJoinSplitSiteClusters(space, L[k].cluster, T, L[k].noofclusters, temp);
            L[k].noofclusters += temp;

            bl_listSweep(&leftl);
            FREEMEMORY(space, T);
            FREEMEMORY(space, lsites);
            nooflsites = 0;
          }

          if(rsites && rsites[noofrsites-1].pos+interval < curstart) {

            T = bl_matchfileScanList (space, &rightl, 
                rsites, noofrsites, interval, curstart, 0, &temp);

            if(realigndev) { 
              noofrightrealigns += bl_matchfileRealign(space, realigndev, &rightrealign, T, temp, 0, set, 
						       curchromidx, fmt, 25,threadno, interval, maxdist);
              bl_listSweep(&rightrealign);
              bl_matchfileDumpUnRealigned(space, realigndev, &leftrealign, &rightrealign, curstart);
            }

            R[k].cluster = 
              bl_matchJoinSplitSiteClusters(space, R[k].cluster, T, R[k].noofclusters, temp);      
            R[k].noofclusters += temp;

            bl_listSweep(&rightl);
            FREEMEMORY(space, T);
            FREEMEMORY(space, rsites);
            noofrsites = 0;
          }

          /*
           * building up lists
           *
           */

          if(acceptorchr) {
            acceptorchridx = bl_matchfileGetChromIndexNumber(index, acceptorchr);
          }

          if(realigndev) { 
            if(!acceptorchr || !donorchr) {
              bookkeeper = ALLOCMEMORY(space, NULL, unsigned char, 1);
              bookkeeper[0] = 0;

              buffer2 = ALLOCMEMORY(space, NULL, char, len+1);
              memmove(buffer2, buffer, len);
              buffer2[len] = 0;
            }

            if(!acceptorchr && !donorchr) {
              buffer3 = ALLOCMEMORY(space, NULL, char, len+1);
              memmove(buffer3, buffer, len);
              buffer3[len] = 0;
            }
          }

          if(donorchr) {

            trans = 0;
            donorchridx = bl_matchfileGetChromIndexNumber(index, donorchr);

            if(strand == '-') {
              trans = (donorflg & SPLIT_PREV_PLUS) ? 1 : 0; 
              adjoint = (acceptorchridx != -1) ? curstart : -1;
	      /*Do not take into account very short introns (size<10)*/
	      if(trans || donorchridx != curchromidx || (donorpos-curend)*(donorpos-curend)>100) { 
		bl_matchfileEnlistMatch (space, &rightl, curstart, curend, 
					 donorchridx, donorpos, adjoint, trans, NULL, NULL);
		
		rsites = bl_matchfileUpdateSplitSites (space, curend, 
						       rsites, &noofrsites, interval, trans);
	      }
	      /*    else {
		fprintf(stderr,"NOT included in splicelist: %s\n",buffer);
		}*/
              //store read for realignment
              if(!acceptorchr && realigndev) { 
                bookkeeper[0] |= LEFTLIST;
                bl_matchfileEnlistMatch (space, &leftrealign, curstart, curend, 
                    donorchridx, donorpos, adjoint, trans, buffer2, bookkeeper);
              }

            } else {
              trans = (!(donorflg & SPLIT_PREV_PLUS)) ? 1 : 0;

              adjoint = (acceptorchridx != -1) ? curend : -1;
	      if(trans || donorchridx != curchromidx || (donorpos-curstart)*(donorpos-curstart)>100) {
              bl_matchfileEnlistMatch (space, &leftl, curstart, curend, 
                  donorchridx, donorpos, adjoint, trans, NULL, NULL);

              lsites = bl_matchfileUpdateSplitSites(space, curstart,
                  lsites, &nooflsites, interval, trans);
	      }
	      /*  else {
		fprintf(stderr,"NOT included in splicelist: %s\n",buffer);
		}*/
              //store read for realignment
              if(!acceptorchr && realigndev) {
                bookkeeper[0] |= RIGHTLIST;
                bl_matchfileEnlistMatch (space, &rightrealign, curstart, curend, 
                    donorchridx, donorpos, adjoint, trans, buffer2, bookkeeper);
              }
            }
          }
 
          if(acceptorchr) {

            if(strand == '-') {
              trans = (acceptorflg & SPLIT_NEXT_PLUS) ? 1 : 0; 
              adjoint = (donorchridx != -1) ? curend : -1;
	      if(trans || acceptorchridx != curchromidx || (acceptorpos-curstart)*(acceptorpos-curstart)>100) {
		bl_matchfileEnlistMatch (space, &leftl, curstart, curend, 
					 acceptorchridx, acceptorpos, adjoint, trans, NULL, NULL);
		
		lsites = bl_matchfileUpdateSplitSites(space, curstart,
						      lsites, &nooflsites, interval, trans);
	      }
	      /*  else {
		fprintf(stderr,"NOT included in splicelist: %s\n",buffer);
		}*/
              //store read for realignment
              if(!donorchr && realigndev) { 
                bookkeeper[0] |= RIGHTLIST;
                bl_matchfileEnlistMatch (space, &rightrealign, curstart, curend, 
                    acceptorchridx, acceptorpos, adjoint, trans, buffer2, bookkeeper);
              }

            } else {
              trans = (!(acceptorflg & SPLIT_NEXT_PLUS)) ? 1 : 0;
              adjoint = (donorchridx != -1) ? curstart : -1;
	      if(trans || acceptorchridx != curchromidx || (acceptorpos-curend)*(acceptorpos-curend)>100) { 
              bl_matchfileEnlistMatch (space, &rightl, curstart, curend, 
                  acceptorchridx, acceptorpos, adjoint, trans, NULL, NULL);

              rsites = bl_matchfileUpdateSplitSites(space, curend, 
                  rsites, &noofrsites, interval, trans);
	      }
	      /* else {
		fprintf(stderr,"NOT included in splicelist: %s\n",buffer);
		}*/
              //store read for realignment
              if(!donorchr && realigndev) { 
                bookkeeper[0] |= LEFTLIST; 
                bl_matchfileEnlistMatch (space, &leftrealign, curstart, curend, 
                    acceptorchridx, acceptorpos, adjoint, trans, buffer2, bookkeeper);
              }
            }

          }

          //store read for realignment
          if(!acceptorchr && !donorchr && realigndev) {
            bookkeeper[0] |= LEFTLIST;
            bookkeeper[0] |= RIGHTLIST;

            bl_matchfileEnlistMatch (space, &leftrealign, curstart,
				     curend, acceptorchridx, acceptorpos, adjoint, trans, buffer2, bookkeeper);

            bl_matchfileEnlistMatch (space, &rightrealign, curstart,
				     curend, acceptorchridx, acceptorpos, adjoint, trans, buffer3, bookkeeper);
          }
	  if(acceptorchr && donorchr && realigndev) {
	    // printf( "Was here!\n");
	     if (realigndev != NULL) {fprintf(realigndev,"%s\n",buffer);} /**/
	}
        }  

        destructStringset(space, token);
      }
      else {
	if (realigndev != NULL) {fprintf(realigndev,"%s\n",buffer);}
      }

      buffer = ALLOCMEMORY(space, buffer, char, buffersize);
      len = 0;

    } else {
      if(ch != '\n') buffer[len++] = ch;
    }
  }

  T = bl_matchfileScanList (space, &leftl, lsites, nooflsites, interval, -1, 1, &temp);
  
  if(realigndev) { 
    noofleftrealigns += bl_matchfileRealign(space, realigndev, &leftrealign, T, temp, 1, set, 
					    curchromidx, fmt, 25,threadno, interval, maxdist);
    bl_listSweep(&leftrealign);
  }

  L[k].cluster = bl_matchJoinSplitSiteClusters(space, L[k].cluster, T, L[k].noofclusters, temp);
  L[k].noofclusters += temp;

  FREEMEMORY(space, T);
  FREEMEMORY(space, lsites);

  T = bl_matchfileScanList (space, &rightl, rsites, noofrsites, interval, -1, 0, &temp);

  if(realigndev) { 
    noofrightrealigns += bl_matchfileRealign(space, realigndev, &rightrealign, T, temp, 0, set, 
					     curchromidx, fmt, 25,threadno, interval, maxdist);

    bl_listSweep(&rightrealign);
    bl_matchfileDumpUnRealigned(space, realigndev, &leftrealign, &rightrealign, -1);
  }

  R[k].cluster = 
    bl_matchJoinSplitSiteClusters(space, R[k].cluster, T, R[k].noofclusters, temp);
  R[k].noofclusters += temp;

  FREEMEMORY(space, T);
  FREEMEMORY(space, rsites);


  bl_listDestruct(&rightl, NULL);
  bl_listDestruct(&leftl, NULL); 
  bl_listDestruct(&rightrealign, NULL);
  bl_listDestruct(&leftrealign, NULL);

  *Lclust = L;
  *Rclust = R;

  FREEMEMORY(space, buffer);
  
  fclose(fp);
  if(gzip) bl_destructgzidxfile(gzf);
  FREEMEMORY(space, gzf);

  return noofrightrealigns+noofleftrealigns;
}

/*-------------------------- bl_matchAddDistClusterLink --------------------------
 *    
 * @brief add a dist cluster link to a cluster
 * @author Steve Hoffmann 
 *   
 */

  void
bl_matchAddDistClusterLink(matchsplitsitecluster_t *e1, Uint chr, Uint link, 
    Uint cnt, Uint realigned, char type)
{

  Uint l;
e1->cnt +=realigned;
  for(l=0; l < e1->noofdistclust; l++) {
    if(e1->distclustchr[l]==chr && 
        e1->distclust[l] == link && e1->distclusttype[l] == type) { 
      e1->distclustcnt[l] += cnt;
      e1->distclustrealigned[l] += realigned;
      /*NEW*/
      e1->distclustcnt[l]+= realigned;
      
      /*end(NEW)*/
      break;
    }
  }

  if(l==e1->noofdistclust) {   
    e1->distclust = bl_insertArray(e1->distclust, e1->noofdistclust, 
        sizeof(Uint), &link, e1->noofdistclust);

    e1->distclustchr = bl_insertArray(e1->distclustchr, e1->noofdistclust, 
        sizeof(Uint), &chr, e1->noofdistclust);

    e1->distclustcnt = bl_insertArray(e1->distclustcnt, e1->noofdistclust, 
        sizeof(Uint), &cnt, e1->noofdistclust);

    e1->distclustrealigned = bl_insertArray(e1->distclustrealigned, e1->noofdistclust, 
        sizeof(Uint), &realigned, e1->noofdistclust);

    e1->distclusttype = bl_insertArray(e1->distclusttype, e1->noofdistclust, 
        sizeof(char), &type, e1->noofdistclust);

    e1->noofdistclust++;
    e1->distclustcnt[l]+=realigned;
    
  }
  return ;
}


/*---------------------- bl_matchAddAdjoinedClusterLink ----------------------
 *    
 * @brief add a link to adjoined clusters
 * @author Steve Hoffmann 
 *   
 */

  void
bl_matchAddAdjoinedClusterLink (matchsplitsitecluster_t *e1, Uint link, Uint cnt)
{
  Uint l;

  for(l=0; l < e1->noofadjclust; l++) {
    if(link == e1->adjclust[l]) {
      e1->adjclustweights[l] += cnt;
      break;
    }
  }
  if(l==e1->noofadjclust) {

    e1->adjclust = bl_insertArray(e1->adjclust, e1->noofadjclust, 
        sizeof(Uint), &link, e1->noofadjclust);

    e1->adjclustweights = bl_insertArray(e1->adjclustweights, e1->noofadjclust, 
        sizeof(Uint), &cnt, e1->noofadjclust);

    e1->noofadjclust++;
  }

  return ;
}

/*----------------------- bl_matchLinkAdjoinedCluster -----------------------
 *    
 * @brief link all clusters that are adjoined by two splice sites
 * @author Steve Hoffmann 
 *   
 */

  void
bl_matchLinkAdjoinedCluster(void *space, matchsplitsiteclusterlist_t *L,
    matchsplitsiteclusterlist_t *R, Uint nchr)
{
  Uint i, j, k, pos, u;
  matchsplitsitecluster_t *e1, *e2;


  for(i=0; i < nchr; i++) {
    for(j=0; j < R[i].noofclusters; j++) {
      e1 = &R[i].cluster[j];
      for(k=0; k < e1->noofadjoints; k++) {

        if(e1->adjoint[k] < 0) {  

          pos = e1->median + e1->adjoint[k];
          u = binarySearch_left(L[i].cluster, L[i].noofclusters, &pos, 
              cmp_matchsplitsitecluster_bin, NULL);

          if(u < L[i].noofclusters) { 
            e2 = &L[i].cluster[u];
            if(e2->a <= pos && pos <= e2->b) {

              bl_matchAddAdjoinedClusterLink (e2, j, e1->adjointcnt[k]);
              bl_matchAddAdjoinedClusterLink (e1, u, e1->adjointcnt[k]);

            } 
          }
        }
      }
    }
  }

  return ;
}


/*--------------------------- bl_matchLinkDistCluster ----------------------------
 *    
 * @brief get loc links from cluster lists
 * @author Steve Hoffmann 
 *   
 */

  void
bl_matchLinkDistCluster (void *space, matchsplitsiteclusterlist_t *R,  
    matchsplitsiteclusterlist_t *L, Uint n)
{
  Uint i, j, k, u, q, trans, cnt, pos, chr;
  matchsplitsitecluster_t *e1, *e2;
  Uint realigned;


#if defined CHECKLINKS
  Uint temp;
#endif
#if defined PARANOID && defined CHECKLINKS
  Uint dpos, v;
#endif

  //all chroms
  for(i=0; i < n; i++) {
    //all clusters
    for(j=0; j < R[i].noofclusters; j++) { 
      //all dist sites R2L and R2R
      for(k=0; k < R[i].cluster[j].distsites; k++){ 
        e1 = &R[i].cluster[j];
        pos = R[i].cluster[j].distpos[k];
        chr = R[i].cluster[j].distchr[k];
        cnt = R[i].cluster[j].distcnt[k];
        trans = R[i].cluster[j].disttrans[k];

        if(R[i].cluster[j].realigned)
          realigned = R[i].cluster[j].realigned[k];
        else
          realigned = 0;

        if(!trans) { 

          u = binarySearch_left(L[chr].cluster, L[chr].noofclusters, &pos, 
              cmp_matchsplitsitecluster_bin, NULL);

          if(u < L[chr].noofclusters) { 
            e2 = &L[chr].cluster[u];
            if(pos >= e2->a && e2->b >= pos) { 

              if(e2->realigned) {
                for(q=0; q < e2->distsites; q++) {
                  if(e2->distchr[q] == chr && e2->distpos[q]>=e1->a && e1->b>=e2->distpos[q])
                    break;
                }
                if(q< e2->distsites) {
                  realigned += e2->realigned[q];
                }
              } 

              bl_matchAddDistClusterLink(e1, chr, u, cnt, realigned, 1);
	      bl_matchAddDistClusterLink(e2, i, j, cnt, realigned, 1);
	      

#ifdef CHECKLINKS        
              R[i].cluster[j].lnkcnt += cnt;
              e2->lnkcnt += cnt;
#endif
#if defined PARANOID && defined CHECKLINKS
              if(R[i].cluster[j].linked[k] < cnt) { 
                dpos =  R[i].cluster[j].a + R[i].cluster[j].rpos[k];
                for(v=0; v < e2->distsites; v++) {

                  if(e2->distchr[v] == i && dpos == e2->distpos[v] && 
                      e2->rpos[v] +e2->a == pos &&
                      e2->disttrans[v] == trans) {
                    assert(cnt == e2->distcnt[v]);
                    e2->linked[v] += cnt;
                    R[i].cluster[j].linked[k] += e2->distcnt[v];

                  }
                }
              }
#endif
            } else {
              fprintf(stderr, "cluster not found (range check) [%d,%d] looking for %d chr %d\n",e1->a, e1->b, pos,chr);
            }          
          } else {
            fprintf(stderr,"cluster not found (search): [%d,%d] looking for %d, chr %d\n",
		    e1->a, e1->b, pos,chr);
            for(q=0; q < L[chr].noofclusters; q++) {
              if(pos >= L[chr].cluster[q].a && L[chr].cluster[q].b >= pos) {
                fprintf(stderr, "found in linear scan\n");
                break;
              }
            } 
            if(q == L[chr].noofclusters) {
              fprintf(stderr, "not found in linear scan: pos:%d:%d\n", chr, pos);
            }
          } 

        } else { 

          u = binarySearch_left(R[chr].cluster, R[chr].noofclusters, &pos, 
              cmp_matchsplitsitecluster_bin, NULL);

          if(u < R[chr].noofclusters) { 
            e2 = &R[chr].cluster[u];
            if(pos >= e2->a && e2->b >= pos) { 

              //TOO MUCH IS LINKED HERE! TRANS1 -> TRANS2 AND TRANS2->TRANS1 
              //- BELoW (PaRANOID) THIS IS AVOIDED BY RIGOROURS POS CHECKING
              //HENCE LINK ONLY IF ELEM >   

              if(chr > i || (chr == i && e2 >= &R[i].cluster[j])) { 
                if(e2->realigned) {
                  for(q=0; q < e2->distsites; q++) {
                    if(e2->distchr[q] == chr && e2->distpos[q]>=e1->a && e1->b>=e2->distpos[q])
                      break;
                  }
                  if(q< e2->distsites) {
		    if(e1 != e2) { /*do not increase twice if splicing back onto itself*/
		      realigned += e2->realigned[q];
		    }
                  }
                } 

                bl_matchAddDistClusterLink(e1, chr, u, cnt, realigned, 2);
		if(e1 != e2) { /*do not increase twice if splicing back onto itself*/
		  bl_matchAddDistClusterLink(e2, i, j, cnt, realigned, 2);
		}
              }
#ifdef CHECKLINKS
              if(chr > i || (chr == i && e2 >= &R[i].cluster[j])) { 
                R[i].cluster[j].lnkcnt += cnt;

                if(chr > i || j != u) {  
                  e2->lnkcnt += cnt;
                }
              }
#endif
#if defined PARANOID && defined CHECKLINKS
              if(R[i].cluster[j].linked[k] < cnt) { 

                dpos =  R[i].cluster[j].a + R[i].cluster[j].rpos[k];

                for(v=0; v < e2->distsites; v++) {
                  if(e2->distchr[v] == i &&  dpos == e2->distpos[v]  && 
                      e2->rpos[v]+e2->a == pos &&
                      e2->disttrans[v] == trans && 
                      e2->linked[v] < e2->distcnt[v]) {

                    assert(cnt == e2->distcnt[v]);
                    if(R[i].cluster[j].linked != e2->linked || k != v)
                      e2->linked[v] += cnt;
                    R[i].cluster[j].linked[k] += e2->distcnt[v];
                  }
                }
              }
#endif
            } else {
              fprintf(stderr, "trans R2R cluster not found (range check)[%d,%d] looking for %d chr %d\n",e1->a, e1->b, pos,chr);
            } 
          } else {
            fprintf(stderr, "trans R2R cluster not found (search)[%d,%d] looking for %d chr %d\n",e1->a, e1->b, pos,chr);
          }
        }
      }
    }

    //L2L
    for(j=0; j < L[i].noofclusters; j++) {

      for(k=0; k < L[i].cluster[j].distsites; k++){ 

        e1 = &L[i].cluster[j];

        pos = L[i].cluster[j].distpos[k];
        chr = L[i].cluster[j].distchr[k];
        cnt = L[i].cluster[j].distcnt[k];
        trans = L[i].cluster[j].disttrans[k];
        if(L[i].cluster[j].realigned)
          realigned = L[i].cluster[j].realigned[k];
        else
          realigned = 0;

        if(trans) { 

          u = binarySearch_left(L[chr].cluster, L[chr].noofclusters, &pos, 
              cmp_matchsplitsitecluster_bin, NULL);

          if(u < L[chr].noofclusters) { 
            e2 = &L[chr].cluster[u];
            if(pos >= e2->a && e2->b >= pos) { 

              if(chr > i || (chr == i && e2 >= &L[i].cluster[j])) { 

                if(e2->realigned) {
                  for(q=0; q < e2->distsites; q++) {
                    if(e2->distchr[q] == chr &&  e2->distpos[q]>=e1->a && e1->b>=e2->distpos[q])
                      break;
                  }
                  if(q< e2->distsites) {
                    realigned += e2->realigned[q];
                  }
                } 

                bl_matchAddDistClusterLink(e1, chr, u, cnt, realigned, 3);
                bl_matchAddDistClusterLink(e2, i, j, cnt, realigned, 3);
              }
#ifdef CHECKLINKS
              if(chr > i || (chr == i && e2 >= &L[i].cluster[j])) { 
                L[i].cluster[j].lnkcnt += cnt;
                if(chr > i || j != u)   
                  e2->lnkcnt += cnt;
              }
#endif
#if defined PARANOID && defined CHECKLINKS
              if(trans && L[i].cluster[j].linked[k] < cnt)  { 
                dpos =  L[i].cluster[j].a + L[i].cluster[j].rpos[k];
                for(v=0; v < e2->distsites; v++) {
                  if(e2->distchr[v] == i &&  dpos == e2->distpos[v] 
                      && pos == e2->rpos[v]+e2->a  
                      && e2->disttrans[v] == trans) {
                    assert(cnt == e2->distcnt[v]);
                    if(L[i].cluster[j].linked != e2->linked || k != v)
                      e2->linked[v] += cnt;
                    L[i].cluster[j].linked[k] += e2->distcnt[v];

                  }         
                }
              }
#endif
            }
          }
        }
      }
    }
  }

#ifdef CHECKLINKS
  for(i=0; i < n; i++) {
    //all clusters
    for(j=0; j < R[i].noofclusters; j++) { 
      //all dist sites R2L and R2R
      temp = 0;

      if(R[i].cluster[j].cnt != R[i].cluster[j].lnkcnt) { 
        fprintf(stderr, "R linkcount wrong: %d <> %d (trans:%d)\n", 
            R[i].cluster[j].cnt, R[i].cluster[j].lnkcnt, R[i].cluster[j].trans);
      }

      for(k=0; k < R[i].cluster[j].distsites; k++){
        temp += R[i].cluster[j].distcnt[k];

#if defined PARANOID && defined CHECKLINKS        
        if(R[i].cluster[j].distcnt[k] != R[i].cluster[j].linked[k]) {
          fprintf(stderr, 
              "R wrong trans: (%d, %d, %d): %d:[%d]-%d <> %d (trans %d)\n", 
              i, j, k,R[i].cluster[j].distchr[k], R[i].cluster[j].distpos[k], 
              R[i].cluster[j].distcnt[k], R[i].cluster[j].linked[k], 
              R[i].cluster[j].disttrans[k]);
        }
#endif
      }
      assert(temp == R[i].cluster[j].cnt);
    }
  }

  for(i=0; i < n; i++) {
    //all clusters
    for(j=0; j < L[i].noofclusters; j++) { 
      //all dist sites L2L 

      if(L[i].cluster[j].cnt != L[i].cluster[j].lnkcnt) { 
        fprintf(stderr, "L linkcount wrong: %d <> %d (trans:%d)\n", 
            L[i].cluster[j].cnt, L[i].cluster[j].lnkcnt, L[i].cluster[j].trans);
      }

      temp = 0;

      for(k=0; k < L[i].cluster[j].distsites; k++){
        temp += L[i].cluster[j].distcnt[k];
#if defined PARANOID && defined CHECKLINKS
        if(L[i].cluster[j].distcnt[k] != L[i].cluster[j].linked[k]) {
          fprintf(stderr, 
              "L wrong trans: (%d, %d, %d): %d:[%d]-%d <> %d (trans %d)\n", 
              i, j, k, L[i].cluster[j].distchr[k], L[i].cluster[j].distpos[k], 
              L[i].cluster[j].distcnt[k], L[i].cluster[j].linked[k], 
              L[i].cluster[j].disttrans[k]);
        }
#endif
      }
      assert(temp == L[i].cluster[j].cnt);
    }
  }

#endif

  return ;
}


/*------------------ bl_matchCompareLinkedClusterSequences -------------------
 *    
 * @brief this function compares the sequences of linked clusters. 
 * The more different the sequences are, the higher the credibility of 
 * the links
 * 
 * @author Steve Hoffmann 
 *   
 */

void
bl_matchCompareLinkedClusterSequences (void *space, fasta_t *set,
				       matchsplitsiteclusterlist_t *R,  
				       matchsplitsiteclusterlist_t *L, Uint n, int maxdist)
{
  Uint i,j,k,q;
  Uint lmrgn, rmrgn, lmrgn2, rmrgn2, chr, pos, distclust, distclust2,
    distchr, distpos, distchr2, 
    distpos2, reflen, reflen2, reflen3, len, len2, noofdistclust;
  Uint edist, minlen;
  char *seq, *seq2, *seq3, disttrans, disttrans2, *rm=NULL, *rm2=NULL;
  int *M, scores[] = {1,-1};
  double accuracy, *emat;
  Uint range=25;
  Alignment al;


  for(i=0; i < n; i++) {
    for(j=0; j < L[i].noofclusters; j++) { 

      //fprintf(stderr, "chrom:%d cluster:%d\n", i, j);
      pos = L[i].cluster[j].median;
      chr = i;
      noofdistclust =  L[i].cluster[j].noofdistclust;
      if (!maxdist || noofdistclust<maxdist) { /**/
	emat = 
	  ALLOCMEMORY(space, NULL, double, ((noofdistclust+1)*(noofdistclust)));
	memset(emat, 0, sizeof(double)*((noofdistclust+1)*(noofdistclust)));
    	//compare splits with left site
	reflen = bl_fastaGetSequenceLength(set, chr);
	lmrgn = (pos > range) ? range : pos;
	rmrgn = (pos+range < reflen) ? range : reflen - pos;
	seq = &bl_fastaGetSequence(set, chr)[pos-lmrgn];

	  for(k=0; k < L[i].cluster[j].noofdistclust; k++) {
	    distchr = L[i].cluster[j].distclustchr[k];
	    distclust = L[i].cluster[j].distclust[k];

	    if(L[i].cluster[j].distclusttype[k] == 1) {   
	      distpos = R[distchr].cluster[distclust].median;
	      disttrans = 0;
	    } else {
	      distpos = L[distchr].cluster[distclust].median;
	      disttrans = 1;
	    }

	    reflen2 = bl_fastaGetSequenceLength(set, distchr);

	    lmrgn = (distpos > lmrgn) ? lmrgn : distpos;
	    rmrgn = (distpos+rmrgn < reflen2) ? rmrgn : reflen2 - distpos;
	    len = rmrgn + lmrgn;

	    seq2 =&bl_fastaGetSequence(set, distchr)[distpos-lmrgn];
	    if(disttrans) {
	      rm = charIUPACcomplement(space, seq2, len);
	      seq2 = rm;
	    }

	    M = swmatrix(space, seq, len, seq2, len, -1, constscr, scores); 
	    initAlignment(&al, seq, len, 0, seq2, len, 0);
	    swtraceback(space, M, seq, len, seq2, len, -1, constscr, scores, &al);
	    MATRIX2D(emat, noofdistclust, 0, k) = .0;

	    minlen = MIN(getUalignlen(&al),getValignlen(&al));
	    edist = getEdist(&al);
	    accuracy = 1.0 - (double)edist/(double)minlen;

	    if(minlen >= 20) {
	      MATRIX2D(emat, noofdistclust, 0, k) = accuracy;
	    }
	    /*
	      if(minlen >= 20 && accuracy >= 0.75) { 
	      fprintf(stderr, "comparing %d w/ %d (%d)\n", pos, distpos, disttrans);
	      fprintf(stderr, "accuracy: %f, minlen:%d\n", accuracy, minlen);
	      showAlign(&al, stderr);
	      }
	    */
	    wrapAlignment(&al);
	    FREEMEMORY(space, M);

	    //compare the splits with each other
	    for(q=k+1; q < L[i].cluster[j].noofdistclust; q++) {

	      distchr2 = L[i].cluster[j].distclustchr[q]; 
	      distclust2 = L[i].cluster[j].distclust[q];

	      if(L[i].cluster[j].distclusttype[q] == 1) {   
		distpos2 = R[distchr2].cluster[distclust2].median;
		disttrans2 = 0;
	      } else {
		distpos2 = L[distchr2].cluster[distclust2].median;
		disttrans2 = 1;
	      }

	      reflen3 = bl_fastaGetSequenceLength(set, distchr2);

	      lmrgn2 = (distpos > range) ? range : distpos;
	      lmrgn2 = (distpos2 > lmrgn2) ? lmrgn2 : distpos2;         
	      rmrgn2 = (distpos+range < reflen2) ? range : reflen2 - distpos;
	      rmrgn2 = (distpos2+rmrgn2 < reflen3) ? rmrgn2 : reflen3 - distpos2;

	      seq3 = &bl_fastaGetSequence(set, distchr2)[distpos2-lmrgn2];
	      len2 = rmrgn2 + lmrgn2;

	      if(disttrans2) {
		rm2 = charIUPACcomplement(space, seq3, len2);
		seq3 = rm2;
	      }

	      M = swmatrix(space, seq2, len, seq3, len2, -1, constscr, scores); 
	      initAlignment(&al, seq2, len, 0, seq3, len2, 0);
	      swtraceback(space, M, seq2, len, seq3, len2, -1, constscr, scores, &al);

	      MATRIX2D(emat, noofdistclust, k+1, q) = .0;
	      MATRIX2D(emat, noofdistclust, q+1, k) = .0;

	      minlen = MIN(getUalignlen(&al),getValignlen(&al));
	      edist = getEdist(&al);
	      accuracy = 1.0 - (double)edist/(double)minlen;

	      if(minlen >= 20) {
		MATRIX2D(emat, noofdistclust, k+1, q) = accuracy;
		MATRIX2D(emat, noofdistclust, q+1, k) = accuracy;
	      }
	      /*
		if(minlen >= 20 && accuracy >= 0.75) { 
		fprintf(stderr, "comparing %d (%d) w/ %d (%d)\n", distpos, disttrans, distpos2, disttrans2);
		fprintf(stderr, "accuracy: %f, minlen:%d\n", accuracy, minlen);
		showAlign(&al, stderr);
		}
	      */
	      wrapAlignment(&al);
	      FREEMEMORY(space, M);
	      if(disttrans2) {
		FREEMEMORY(space, rm2);
	      }
	    }

	    if(disttrans) {
	      FREEMEMORY(space, rm);
	    }
	  }
   
	L[i].cluster[j].emat = emat;
      }
      else {
	L[i].cluster[j].emat = NULL /*geht das? will it free?*/;
      }
    }
  }

  return ;
}

/*-------------------- bl_matchShowMatchsplitsiteclusterlist -----------------
 *    
 * @brief shows the lists
 * @author Steve Hoffmann 
 *   
 */

  void
bl_matchShowMatchsplitsiteclusterlist (void *space, 
    matchsplitsiteclusterlist_t *l, matchsplitsiteclusterlist_t *r, 
    Uint n, char type)
{
  Uint i, k, j, u;

  for(i=0; i < n; i++) {
    k = l[i].noofclusters;
    fprintf(stderr, "chromosome %d  with %d clusters\n",i, k);
    
    for(j=0; j < k; j++) { 
      fprintf(stderr,"cluster %d\t[%d,%d]", j, 
          l[i].cluster[j].a, l[i].cluster[j].b);
#ifdef CHECKLINKS
      fprintf(stderr,"-%d",l[i].cluster[j].cnt);
#endif
      fprintf(stderr, "\n");

      fprintf(stderr, "\tdistant loci:\n");
      for(u=0; u < l[i].cluster[j].distsites; u++) {
        fprintf(stderr, "\t\tto: %d:%d  cnt:%d (trans:%d)",
            l[i].cluster[j].distchr[u], l[i].cluster[j].distpos[u], 
            l[i].cluster[j].distcnt[u], l[i].cluster[j].disttrans[u]);
        if(l[i].cluster[j].realigned) {
            fprintf(stderr, "realigned: %d\n",l[i].cluster[j].realigned[u]);
        } else {
          fprintf(stderr, "\n");
        }
      }
      
      fprintf(stderr, "\tdistant cluster:\n");
      for(u=0; u < l[i].cluster[j].noofdistclust; u++) {
        fprintf(stderr, "\t\tto: %d:%d  cnt:%d (type:%d) realigned:%d\n",
            l[i].cluster[j].distclustchr[u], l[i].cluster[j].distclust[u],
            l[i].cluster[j].distclustcnt[u], l[i].cluster[j].distclusttype[u], 
            l[i].cluster[j].distclustrealigned[u]);
      }
      
     fprintf(stderr, "\tadjoint loci:\n");
      for(u=0; u < l[i].cluster[j].noofadjoints; u++) {
        fprintf(stderr, "\t\tto: %d  cnt:%d\n",
            l[i].cluster[j].adjoint[u], 
            l[i].cluster[j].adjointcnt[u]);
      }

 
      fprintf(stderr, "\tadjoint cluster:\n");
      for(u=0; u < l[i].cluster[j].noofadjclust; u++) {
        fprintf(stderr, "\t\tto: %d cnt:%d [%d,%d]\n",
            l[i].cluster[j].adjclust[u], l[i].cluster[j].adjclustweights[u], 
            r[i].cluster[l[i].cluster[j].adjclust[u]].a, 
            r[i].cluster[l[i].cluster[j].adjclust[u]].b);
      }
      fprintf(stderr, "\n");
    }
  }
  return ;
}



/*----------------- bl_matchGetMatchsplitsiteclusterlistBED ------------------
 *    
 * @brief get a bed from matchsplitsitecluster
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchGetMatchsplitsiteclusterlistBED (void *space, fasta_t *set, FILE *normdev, FILE *transdev,
					 matchsplitsiteclusterlist_t *l, matchsplitsiteclusterlist_t *r, 
					 Uint n, int maxdist)
{

  Uint i, k, j, u, q, v, distclust, distchr, noofdistclust, noofdistclust2;
  char flag;
  double *emat, max, max2;
  Lint distance;

  //iter chromosomes
  for(i=0; i < n; i++) {
    k = l[i].noofclusters;
    //iter clusters
    for(j=0; j < k; j++) { 
      //iter dist clusters

      noofdistclust = l[i].cluster[j].noofdistclust;
      // no fancy stuff if cluster is promiscuitive
      if (!maxdist || noofdistclust<maxdist) { 
	for(u=0; u < noofdistclust; u++) {

	  distchr = l[i].cluster[j].distclustchr[u];
	  distclust = l[i].cluster[j].distclust[u];
	  emat = l[i].cluster[j].emat;

	  max = MATRIX2D(emat, noofdistclust, 0, u);
	
	  for(q=0; q < noofdistclust; q++) {
	    max = MAX(max, MATRIX2D(emat, noofdistclust, q+1, u));
	  }

	  flag = 'P';
	  if(max >= 0.75) flag = 'F';
        
	  distance = (l[i].cluster[j].distclusttype[u] == 1) ?
	    labs((Lint)l[i].cluster[j].median - r[distchr].cluster[distclust].median) :          
	    labs((Lint)l[i].cluster[j].median - l[distchr].cluster[distclust].median) ;


	  if( l[i].cluster[j].distclusttype[u] == 1 &&
	      distchr == i &&
	      distance <= 200000) { 

	    emat = r[distchr].cluster[distclust].emat;

	    noofdistclust2 = r[distchr].cluster[distclust].noofdistclust;
	    if (!maxdist || noofdistclust2<maxdist) { //second cluster is promiscuitive
	      for(v=0; v < noofdistclust2; v++) {
		if(r[distchr].cluster[distclust].distclustchr[v] == i &&
		   r[distchr].cluster[distclust].distclust[v] == j)
		  break;
	      }
	      assert(v < noofdistclust2);
	      max2 = MATRIX2D(emat, noofdistclust2, 0, v);
	      for(q=0; q < noofdistclust2;  q++) {
		max2 = MAX(max2, MATRIX2D(emat, noofdistclust2, q+1, v));
	      }
	      flag ='P';
	      if(max >= 0.75 || max2 >= 0.75) flag = 'F'; 
	    }  
	    else {
	      flag = 'M';//second cluster is promiscuitive
	    }
	   
	    if(l[i].cluster[j].median<r[distchr].cluster[distclust].median) {
	      /*here: include information about circulars*/
	      fprintf(normdev, "%s\t%d\t%d\tsplits:%d:%d:%d:C:%c\t%d\t+\n",
		      bl_fastaGetDescription(set, i), l[i].cluster[j].median,
		      r[distchr].cluster[distclust].median, 
		      l[i].cluster[j].distclustcnt[u], l[i].cluster[j].cnt,
		      r[distchr].cluster[distclust].cnt, flag, l[i].cluster[j].distclustrealigned[u]);
	    }
	    else {
	      fprintf(normdev, "%s\t%d\t%d\tsplits:%d:%d:%d:N:%c\t%d\t+\n",
		      bl_fastaGetDescription(set, i),  r[distchr].cluster[distclust].median,
		      l[i].cluster[j].median,
		      l[i].cluster[j].distclustcnt[u], r[distchr].cluster[distclust].cnt, 
		      l[i].cluster[j].cnt, flag, l[i].cluster[j].distclustrealigned[u]);
	    }
	  
	  } else {

	    if(distance > 200000 || i != l[i].cluster[j].distclustchr[u]) { 
	      if(l[i].cluster[j].distclusttype[u] > 1){ 

		if(j <= distclust) { 
		  emat = l[distchr].cluster[distclust].emat;
		  noofdistclust2 = l[distchr].cluster[distclust].noofdistclust;
		  if (!maxdist || noofdistclust2<maxdist) {  //second cluster is promiscuitive
		    for(v=0; v < noofdistclust2; v++) {
		      if(l[distchr].cluster[distclust].distclustchr[v] == i &&
			 l[distchr].cluster[distclust].distclust[v] == j)
			break;
		    }
		    assert(v < noofdistclust2);
		    max2 = MATRIX2D(emat, noofdistclust2, 0, v);
		    for(q=0; q < noofdistclust2;  q++) {
		      max2 = MAX(max2, MATRIX2D(emat, noofdistclust2, q+1, v));
		    }

		    flag ='P';
		    if(max >= 0.75 || max2 >= 0.75) flag = 'F';
		  } 
		  else {
		    flag = 'M'; //second cluster is promiscuitive
		  }
		 
		  //L2L
		  fprintf(transdev, "%s\t%d\t%d\tdiststrandsplice:%s:%d:%d:%d:%d:L:%c\t%d\t+\n",
			  bl_fastaGetDescription(set, i), l[i].cluster[j].median,
			  l[i].cluster[j].median, 
			  bl_fastaGetDescription(set, distchr),
			  l[distchr].cluster[distclust].median,
			  l[i].cluster[j].distclustcnt[u], 
			  l[i].cluster[j].cnt,
			  l[distchr].cluster[distclust].cnt, flag, 
			  l[i].cluster[j].distclustrealigned[u]);

		  fprintf(transdev, "%s\t%d\t%d\tdiststrandsplice:%s:%d:%d:%d:%d:L:%c\t%d\t+\n",
			  bl_fastaGetDescription(set, l[i].cluster[j].distclustchr[u]),
			  l[distchr].cluster[distclust].median,
			  l[distchr].cluster[distclust].median,
                  
			  bl_fastaGetDescription(set, i), l[i].cluster[j].median,
			  l[i].cluster[j].distclustcnt[u],
			  l[distchr].cluster[distclust].cnt,
			  l[i].cluster[j].cnt, flag,
			  l[i].cluster[j].distclustrealigned[u]);
              	 
		}
	      } else {

 
		emat = r[distchr].cluster[distclust].emat;

		noofdistclust2 = r[distchr].cluster[distclust].noofdistclust;
		if (!maxdist || noofdistclust2<maxdist) {//second cluster is promiscuitive
		  for(v=0; v < noofdistclust2; v++) {
		    if(r[distchr].cluster[distclust].distclustchr[v] == i &&
		       r[distchr].cluster[distclust].distclust[v] == j)
		      break;
		  }
		  assert(v < noofdistclust2);
		  max2 = MATRIX2D(emat, noofdistclust2, 0, v);
		  for(q=0; q < noofdistclust2;  q++) {
		    max2 = MAX(max2, MATRIX2D(emat, noofdistclust2, q+1, v));
		  }

		  flag ='P';
		  if(max >= 0.75 || max2 >= 0.75) flag = 'F';
		} 
		else {
		  flag = 'M'; //second cluster is promiscuitive
		}
		/*  here if i want to put out promiscuitives*/
		//L2R
		fprintf(transdev, "%s\t%d\t%d\tdistsplice:%s:%d:%d:%d:%d:L:%c\t%d\t+\n",
			bl_fastaGetDescription(set, i), l[i].cluster[j].median,
			l[i].cluster[j].median, 
			bl_fastaGetDescription(set, l[i].cluster[j].distclustchr[u]),
			r[distchr].cluster[distclust].median,
                  
			l[i].cluster[j].distclustcnt[u], 
			l[i].cluster[j].cnt,
			r[distchr].cluster[distclust].cnt, 
			flag,
			l[i].cluster[j].distclustrealigned[u]);

		fprintf(transdev, "%s\t%d\t%d\tdistsplice:%s:%d:%d:%d:%d:R:%c\t%d\t+\n",
			bl_fastaGetDescription(set, l[i].cluster[j].distclustchr[u]),
			r[distchr].cluster[distclust].median,
			r[distchr].cluster[distclust].median,
			bl_fastaGetDescription(set, i), l[i].cluster[j].median,
                  
			l[i].cluster[j].distclustcnt[u],
			r[distchr].cluster[distclust].cnt,
			l[i].cluster[j].cnt, 
			flag,
			l[i].cluster[j].distclustrealigned[u]);
	
	      }
	    } else {
	      if(j <= distclust) { 

		emat = l[distchr].cluster[distclust].emat;

		noofdistclust2 = l[distchr].cluster[distclust].noofdistclust;
		if (!maxdist || noofdistclust2<maxdist) {//second cluster is promiscuitive
		  for(v=0; v < noofdistclust2; v++) {
		    if(l[distchr].cluster[distclust].distclustchr[v] == i &&
		       l[distchr].cluster[distclust].distclust[v] == j)
		      break;
		  }
		  assert(v < noofdistclust2);
		  max2 = MATRIX2D(emat, noofdistclust2, 0, v);
		  for(q=0; q < noofdistclust2;  q++) {
		    max2 = MAX(max2, MATRIX2D(emat, noofdistclust2, q+1, v));
		  }

		  flag ='P';
		  if(max >= 0.75 || max2 >= 0.75) flag = 'F';
		} 
		else {
		  flag = 'M'; //second cluster is promiscuitive
		}
		/*	   here if i want to put out promiscuitives*/
	
		//L2L
		fprintf(transdev, "%s\t%d\t%d\tstrandsplice:%d:%d:%d:L:%c\t%d\t+\n",
			bl_fastaGetDescription(set, i), l[i].cluster[j].median,
			l[distchr].cluster[distclust].median,

			l[i].cluster[j].distclustcnt[u], 
			l[i].cluster[j].cnt,
			l[distchr].cluster[distclust].cnt, 
			flag,
			l[i].cluster[j].distclustrealigned[u]);
	      }
	    }
	  }
	}
      }
      else {/*put out all multiple splice sites also*/
	for(u=0; u < noofdistclust; u++) {

	  distchr = l[i].cluster[j].distclustchr[u];
	  distclust = l[i].cluster[j].distclust[u];
	  //	  emat = l[i].cluster[j].emat;

	  //max = MATRIX2D(emat, noofdistclust, 0, u);
	
	  //for(q=0; q < noofdistclust; q++) {
          //max = MAX(max, MATRIX2D(emat, noofdistclust, q+1, u));
	  //}

	  flag = 'M';
	  //if(max >= 0.75) flag = 'F';
        
	  distance = (l[i].cluster[j].distclusttype[u] == 1) ?
	    labs((Lint)l[i].cluster[j].median - r[distchr].cluster[distclust].median) :          
	    labs((Lint)l[i].cluster[j].median - l[distchr].cluster[distclust].median) ;


	  if( l[i].cluster[j].distclusttype[u] == 1 &&
	      distchr == i &&
	      distance <= 200000) { 

	    // emat = r[distchr].cluster[distclust].emat;

	    /* noofdistclust2 = r[distchr].cluster[distclust].noofdistclust;
	       if (noofdistclust2<maxdist) {
	       for(v=0; v < noofdistclust2; v++) {
	       if(r[distchr].cluster[distclust].distclustchr[v] == i &&
	       r[distchr].cluster[distclust].distclust[v] == j)
	       break;
	       }
	       assert(v < noofdistclust2);
	       max2 = MATRIX2D(emat, noofdistclust2, 0, v);
	       for(q=0; q < noofdistclust2;  q++) {
	       max2 = MAX(max2, MATRIX2D(emat, noofdistclust2, q+1, v));
	       }

	       flag ='P';
	       if(max >= 0.75 || max2 >= 0.75) flag = 'F'; 
	       }  
	       else {
	       flag = 'M';
	       }*/
	 
	    if(l[i].cluster[j].median<r[distchr].cluster[distclust].median) {
	      fprintf(normdev, "%s\t%d\t%d\tsplits:%d:%d:%d:L:%c\t%d\t+\n",
		      bl_fastaGetDescription(set, i), l[i].cluster[j].median,
		      r[distchr].cluster[distclust].median, 
		      l[i].cluster[j].distclustcnt[u], l[i].cluster[j].cnt,
		      r[distchr].cluster[distclust].cnt, flag, l[i].cluster[j].distclustrealigned[u]);
	    }
	    else {
	      fprintf(normdev, "%s\t%d\t%d\tsplits:%d:%d:%d:R:%c\t%d\t+\n",
		      bl_fastaGetDescription(set, i),  r[distchr].cluster[distclust].median,
		      l[i].cluster[j].median,
		      l[i].cluster[j].distclustcnt[u], r[distchr].cluster[distclust].cnt, 
		      l[i].cluster[j].cnt, flag, l[i].cluster[j].distclustrealigned[u]);
	    }
	  
	  } else {

	    if(distance > 200000 || i != l[i].cluster[j].distclustchr[u]) { 
	      if(l[i].cluster[j].distclusttype[u] > 1){ 

		if(j <= distclust) { 
		  /*
		    emat = l[distchr].cluster[distclust].emat;
		    noofdistclust2 = l[distchr].cluster[distclust].noofdistclust;
		    if (noofdistclust2<maxdist) {
		    for(v=0; v < noofdistclust2; v++) {
		    if(l[distchr].cluster[distclust].distclustchr[v] == i &&
		    l[distchr].cluster[distclust].distclust[v] == j)
		    break;
		    }
		    assert(v < noofdistclust2);
		    max2 = MATRIX2D(emat, noofdistclust2, 0, v);
		    for(q=0; q < noofdistclust2;  q++) {
		    max2 = MAX(max2, MATRIX2D(emat, noofdistclust2, q+1, v));
		    }

		    flag ='P';
		    if(max >= 0.75 || max2 >= 0.75) flag = 'F';
		    } 
		    else {
		    flag = 'M';
		    }*/
		  /*   here if i want to put out promiscuitives*/
		  //L2L
		  fprintf(transdev, "%s\t%d\t%d\tdiststrandsplice:%s:%d:%d:%d:%d:L:%c\t%d\t+\n",
			  bl_fastaGetDescription(set, i), l[i].cluster[j].median,
			  l[i].cluster[j].median, 
			  bl_fastaGetDescription(set, distchr),
			  l[distchr].cluster[distclust].median,
			  l[i].cluster[j].distclustcnt[u], 
			  l[i].cluster[j].cnt,
			  l[distchr].cluster[distclust].cnt, flag, 
			  l[i].cluster[j].distclustrealigned[u]);

		  fprintf(transdev, "%s\t%d\t%d\tdiststrandsplice:%s:%d:%d:%d:%d:L:%c\t%d\t+\n",
			  bl_fastaGetDescription(set, l[i].cluster[j].distclustchr[u]),
			  l[distchr].cluster[distclust].median,
			  l[distchr].cluster[distclust].median,
                  
			  bl_fastaGetDescription(set, i), l[i].cluster[j].median,
			  l[i].cluster[j].distclustcnt[u],
			  l[distchr].cluster[distclust].cnt,
			  l[i].cluster[j].cnt, flag,
			  l[i].cluster[j].distclustrealigned[u]);
		}
	      } else {

 
		/*  emat = r[distchr].cluster[distclust].emat;

		    noofdistclust2 = r[distchr].cluster[distclust].noofdistclust;
		    if (noofdistclust2<maxdist) {
		    for(v=0; v < noofdistclust2; v++) {
		    if(r[distchr].cluster[distclust].distclustchr[v] == i &&
		    r[distchr].cluster[distclust].distclust[v] == j)
                    break;
		    }
		    assert(v < noofdistclust2);
		    max2 = MATRIX2D(emat, noofdistclust2, 0, v);
		    for(q=0; q < noofdistclust2;  q++) {
		    max2 = MAX(max2, MATRIX2D(emat, noofdistclust2, q+1, v));
		    }

		    flag ='P';
		    if(max >= 0.75 || max2 >= 0.75) flag = 'F';
		    } 
		    else {
		    flag = 'M';
		    }*/
		/*  here if i want to put out promiscuitives*/
		//L2R
		fprintf(transdev, "%s\t%d\t%d\tdistsplice:%s:%d:%d:%d:%d:L:%c\t%d\t+\n",
			bl_fastaGetDescription(set, i), l[i].cluster[j].median,
			l[i].cluster[j].median, 
			bl_fastaGetDescription(set, l[i].cluster[j].distclustchr[u]),
			r[distchr].cluster[distclust].median,
                  
			l[i].cluster[j].distclustcnt[u], 
			l[i].cluster[j].cnt,
			r[distchr].cluster[distclust].cnt, 
			flag,
			l[i].cluster[j].distclustrealigned[u]);

		fprintf(transdev, "%s\t%d\t%d\tdistsplice:%s:%d:%d:%d:%d:R:%c\t%d\t+\n",
			bl_fastaGetDescription(set, l[i].cluster[j].distclustchr[u]),
			r[distchr].cluster[distclust].median,
			r[distchr].cluster[distclust].median,
			bl_fastaGetDescription(set, i), l[i].cluster[j].median,
                  
			l[i].cluster[j].distclustcnt[u],
			r[distchr].cluster[distclust].cnt,
			l[i].cluster[j].cnt, 
			flag,
			l[i].cluster[j].distclustrealigned[u]);
		
	      }
	    } else {
	      if(j <= distclust) { 
		/*
		  emat = l[distchr].cluster[distclust].emat;

		  noofdistclust2 = l[distchr].cluster[distclust].noofdistclust;
		  if (noofdistclust2<maxdist) {
		  for(v=0; v < noofdistclust2; v++) {
		  if(l[distchr].cluster[distclust].distclustchr[v] == i &&
		  l[distchr].cluster[distclust].distclust[v] == j)
		  break;
		  }
		  assert(v < noofdistclust2);
		  max2 = MATRIX2D(emat, noofdistclust2, 0, v);
		  for(q=0; q < noofdistclust2;  q++) {
		  max2 = MAX(max2, MATRIX2D(emat, noofdistclust2, q+1, v));
		  }

		  flag ='P';
		  if(max >= 0.75 || max2 >= 0.75) flag = 'F';
		  } 
		  else {
		  flag = 'M';
		  }*/
		/*	   here if i want to put out promiscuitives*/
	
		//L2L
		fprintf(transdev, "%s\t%d\t%d\tstrandsplice:%d:%d:%d:L:%c\t%d\t+\n",
			bl_fastaGetDescription(set, i), l[i].cluster[j].median,
			l[distchr].cluster[distclust].median,

			l[i].cluster[j].distclustcnt[u], 
			l[i].cluster[j].cnt,
			l[distchr].cluster[distclust].cnt, 
			flag,
			l[i].cluster[j].distclustrealigned[u]);
	    
	      }
	    }
	  }
	}
      }
    }
  }
  


  for(i=0; i < n; i++) {
    k = r[i].noofclusters;
    //iter clusters
    for(j=0; j < k; j++) { 
      //iter dist clusters
      noofdistclust = r[i].cluster[j].noofdistclust;
      if (!maxdist || noofdistclust<maxdist) {
	for(u=0; u < noofdistclust; u++) {
	  if (r[i].cluster[j].distclusttype[u] == 2) { 

	    distchr = r[i].cluster[j].distclustchr[u];
	    distclust = r[i].cluster[j].distclust[u];
	    emat = r[i].cluster[j].emat;

	    max = MATRIX2D(emat, noofdistclust, 0, u);
	    for(q=0; q < noofdistclust; q++) {
	      max = MAX(max, MATRIX2D(emat, noofdistclust, q+1, u));
	    }
      
	    emat = r[distchr].cluster[distclust].emat;

	    noofdistclust2 = r[distchr].cluster[distclust].noofdistclust;
	    if (!maxdist || noofdistclust2<maxdist) {
	      for(v=0; v < noofdistclust2; v++) {
		if(r[distchr].cluster[distclust].distclustchr[v] == i &&
		   r[distchr].cluster[distclust].distclust[v] == j)
		  break;
	      }
	      assert(v < noofdistclust2);
          
	      max2 = MATRIX2D(emat, noofdistclust2, 0, v);
	      for(q=0; q < noofdistclust2;  q++) {
		max2 = MAX(max2, MATRIX2D(emat, noofdistclust2, q+1, v));
	      }

	      flag ='P';
	      if(max >= 0.75 || max2 >= 0.75) flag = 'F';
	    } 
	    else {
	      flag = 'M';
	    }
	    /*here if i want to put out promiscuitives*/
	
	      distance = 
		labs((Lint)r[i].cluster[j].median - r[distchr].cluster[distclust].median);


	      if(distance <= 200000 && i == distchr) {
    
		if(j <= distclust) { 
		  //R2R
		  fprintf(transdev, "%s\t%d\t%d\tstrandsplice:%d:%d:%d:R:%c\t%d\t+\n",
			  bl_fastaGetDescription(set, i), r[i].cluster[j].median,
			  r[i].cluster[distclust].median,
			  r[i].cluster[j].distclustcnt[u], r[i].cluster[j].cnt,
			  r[i].cluster[distclust].cnt, 
			  flag, 
			  r[i].cluster[j].distclustrealigned[u]);
		}
	      } else {
           
		if(j <= distclust) { 
		  fprintf(transdev, "%s\t%d\t%d\tdiststrandsplice:%s:%d:%d:%d:%d:R:%c\t%d\t+\n",
			  bl_fastaGetDescription(set, i), r[i].cluster[j].median,
			  r[i].cluster[j].median, 
                  
			  bl_fastaGetDescription(set, distchr),
			  r[distchr].cluster[distclust].median,
			  r[i].cluster[j].distclustcnt[u], 
			  r[i].cluster[j].cnt,
			  r[distchr].cluster[distclust].cnt, 
			  flag,
			  r[i].cluster[j].distclustrealigned[u]);

		  fprintf(transdev, "%s\t%d\t%d\tdiststrandsplice:%s:%d:%d:%d:%d:R:%c\t%d\t+\n",
			  bl_fastaGetDescription(set, r[i].cluster[j].distclustchr[u]),
			  r[distchr].cluster[distclust].median,
			  r[distchr].cluster[distclust].median,
                  
			  bl_fastaGetDescription(set, i), 
			  r[i].cluster[j].median,
			  r[i].cluster[j].distclustcnt[u],
			  r[distchr].cluster[distclust].cnt, 
			  r[i].cluster[j].cnt, 
			  flag,
			  r[i].cluster[j].distclustrealigned[u]);
		}
	      }
	    }
      
	  }
	}
      
      else {
	flag = 'M';
	for(u=0; u < noofdistclust; u++) {
	  if (r[i].cluster[j].distclusttype[u] == 2) { 

	    distchr = r[i].cluster[j].distclustchr[u];
	    distclust = r[i].cluster[j].distclust[u];
	    distance = 
	      labs((Lint)r[i].cluster[j].median - r[distchr].cluster[distclust].median);


	    if(distance <= 200000 && i == distchr) {
    
		if(j <= distclust) { 
		  //R2R
		  fprintf(transdev, "%s\t%d\t%d\tstrandsplice:%d:%d:%d:R:%c\t%d\t+\n",
			  bl_fastaGetDescription(set, i), r[i].cluster[j].median,
			  r[i].cluster[distclust].median,
			  r[i].cluster[j].distclustcnt[u], r[i].cluster[j].cnt,
			  r[i].cluster[distclust].cnt, 
			  flag, 
			  r[i].cluster[j].distclustrealigned[u]);
		}
	      } else {
           
	      if(j <= distclust) { 
		  fprintf(transdev, "%s\t%d\t%d\tdiststrandsplice:%s:%d:%d:%d:%d:R:%c\t%d\t+\n",
			  bl_fastaGetDescription(set, i), r[i].cluster[j].median,
			  r[i].cluster[j].median, 
                  
			  bl_fastaGetDescription(set, distchr),
			  r[distchr].cluster[distclust].median,
			  r[i].cluster[j].distclustcnt[u], 
			  r[i].cluster[j].cnt,
			  r[distchr].cluster[distclust].cnt, 
			  flag,
			  r[i].cluster[j].distclustrealigned[u]);

		  fprintf(transdev, "%s\t%d\t%d\tdiststrandsplice:%s:%d:%d:%d:%d:R:%c\t%d\t+\n",
			  bl_fastaGetDescription(set, r[i].cluster[j].distclustchr[u]),
			  r[distchr].cluster[distclust].median,
			  r[distchr].cluster[distclust].median,
                  
			  bl_fastaGetDescription(set, i), 
			  r[i].cluster[j].median,
			  r[i].cluster[j].distclustcnt[u],
			  r[distchr].cluster[distclust].cnt, 
			  r[i].cluster[j].cnt, 
			  flag,
			  r[i].cluster[j].distclustrealigned[u]);
		}
	      }
	    }
	}
      }
    }
  }

  return ;
}

/*---------------- bl_matchDestructMatchsplitsiteclusterlist -----------------
 *    
 * @brief destruct the lists
 * @author Steve Hoffmann 
 *   
 */

  void
bl_matchDestructMatchsplitsiteclusterlist (void *space, 
    matchsplitsiteclusterlist_t *l, Uint n)
{
  Uint i,j;
  for(i=0; i < n;i++) {
    for(j=0; j < l[i].noofclusters; j++) { 
      FREEMEMORY(space, l[i].cluster[j].distpos);  
      FREEMEMORY(space, l[i].cluster[j].distchr);
      FREEMEMORY(space, l[i].cluster[j].disttrans);
#if defined PARANOID && defined CHECKLINKS
      FREEMEMORY(space, l[i].cluster[j].rpos);
      FREEMEMORY(space, l[i].cluster[j].linked);
#endif
      FREEMEMORY(space, l[i].cluster[j].distcnt);
      if(l[i].cluster[j].adjoint) { 
        FREEMEMORY(space, l[i].cluster[j].adjoint);
        FREEMEMORY(space, l[i].cluster[j].adjointcnt);
      }
      if(l[i].cluster[j].realigned) { 
        FREEMEMORY(space, l[i].cluster[j].realigned);
      }

      if(l[i].cluster[j].emat) {
        FREEMEMORY(space, l[i].cluster[j].emat);
      }

      if(l[i].cluster[j].adjclust) {
         FREEMEMORY(space, l[i].cluster[j].adjclust);
         FREEMEMORY(space, l[i].cluster[j].adjclustweights);
      }

      if(l[i].cluster[j].distclust) {
        FREEMEMORY(space, l[i].cluster[j].distclust);
        FREEMEMORY(space, l[i].cluster[j].distclustchr);
        FREEMEMORY(space, l[i].cluster[j].distclustcnt);
        FREEMEMORY(space, l[i].cluster[j].distclusttype);
        FREEMEMORY(space, l[i].cluster[j].distclustrealigned);
      }
    }

    l[i].noofclusters = 0;
    FREEMEMORY(space, l[i].cluster);
    l[i].cluster = NULL;
  }
  return ;
}

#ifdef REALIGNTEST

unsigned char mute = 0;
char *ntcode;

int main(int argc, char **argv) {
  
  realign_t nfo;
  
  manopt_optionset optset;
  manopt_arg *unflagged; 
  manopt_arg *queries;
  manopt_arg *dbfilenames;
  manopt_intconstraint threadconstraint;

  Uint realigned=0, desclen;
  char norealign=0;
  matchsplitsiteclusterlist_t *L = NULL, *R = NULL;
  matchfileindex_t *index;
  matchfile_t **files = NULL; 
  unsigned char gzip = 0;
  char version[]="0.1", expand=0, *desc;
  int i;
  char verbose = 0;
  Uint nchr;
  Uint prefixlen=0;
  
  threadconstraint.max = 3000;
  threadconstraint.min = 1;
  ra_setdefault(&nfo);
  nfo.splitfile="splicesites.bed";
  nfo.transfile="transrealigned.bed";
  nfo.maxdist=100;
  initIUPAC(1,1); 
  manopt_initoptionset(&optset, argv[0], NULL, 
      "Heuristic mapping of short sequences\n",
      "SEGEMEHL is free software for non-commercial use \n  (C) 2008 Bioinformatik Leipzig\n",
      version,
      "Please report bugs to steve@bioinf.uni-leipzig.de"); 
  manopt(&optset, LISTOPT, 1, 'd', "database", 
	 "list of path/filename(s) of database sequence(s)", "<file> [<file> ...]", 
	 NULL, NULL);
  manopt(&optset, LISTOPT, 1, 'q', "query", 
	 "path/filename of alignment file", "<file> [<file> ...]", NULL, NULL); 
  manopt(&optset, FLAG, 0, 'E', "expand", 
	 "expand", NULL, NULL, &expand);
  manopt(&optset, FLAG, 0, 'v', "verbose", 
	 "verbose", NULL, NULL, &verbose);
  manopt(&optset, FLAG, 0, 'n', "norealign", 
	 "do not realign", NULL, NULL, &norealign);
  manopt(&optset, REQUINTOPT, 0, 't', "threads", 
	 "start <n> threads for realigning", "<n>",
         &threadconstraint, &nfo.threadno);
  manopt(&optset, REQSTRINGOPT, 0, 'U', "splitfile", 
      "path/filename of the split bedfile", "<file>", &nfo.splitfile,  &nfo.splitfile);
  manopt(&optset, REQSTRINGOPT, 0, 'T', "transfile", "path/filename of bed files containing trans-split", "<file>",&nfo.transfile , &nfo.transfile);
  manopt(&optset, REQSTRINGOPT, 0, 'o', "outfile", "path/filename of output sam file", "<file>",NULL , &nfo.outfile);
  manopt(&optset, REQUINTOPT, 0, 'M', "maxdist", 
      "max number of distant sites to consider, 0 to disable", "<n>", NULL, &nfo.maxdist);

  //open file for bed output 
  
  unflagged = manopt_getopts(&optset, argc, argv);
 
  if(unflagged->noofvalues > 1) { 
    manopt_help(&optset, "unknown argument(s)\n");
  }
  nfo.transdev=fopen(nfo.transfile,"w");
  nfo.normdev=fopen(nfo.splitfile,"w");
  if(nfo.outfile != NULL) {
    nfo.realigndev = fopen(nfo.outfile,"w");
  }
  
  if(norealign) { 
    nfo.realigndev = NULL;
  }
 
    pthread_mutex_init(&inout, NULL);
    pthread_mutex_init(&mutkuh, NULL);
  
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
    realigned += 
      bl_matchfileRealignScanFileNew(nfo.space, files[i], nfo.realigndev, nfo.fasta, 255, &L, &R, &nchr,nfo.threadno, nfo.maxdist); 
  }
  /*  printf("realigned:%d\n",realigned);*/
  bl_matchLinkAdjoinedCluster(nfo.space, L, R, nchr);
  bl_matchLinkDistCluster (nfo.space, R, L, nchr); 
  bl_matchCompareLinkedClusterSequences (nfo.space, nfo.fasta, R, L, nchr, nfo.maxdist);
  bl_matchCompareLinkedClusterSequences (nfo.space, nfo.fasta, L, R, nchr, nfo.maxdist);
  bl_matchGetMatchsplitsiteclusterlistBED (nfo.space, nfo.fasta, nfo.normdev, nfo.transdev, L, R, nchr, nfo.maxdist);

  if(verbose) { 
    fprintf(stderr, "LEFT -----------------\n");
    bl_matchShowMatchsplitsiteclusterlist(nfo.space, L, R, nchr, 0);

    fprintf(stderr, "RIGHT -----------------\n");
    bl_matchShowMatchsplitsiteclusterlist(nfo.space, R, L, nchr, 1);
        fprintf(stderr, "realigned:%d\n",realigned);
  }

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
}

#endif

/* T[i].distsites;*/
void *threadrealign(int begin, int stop, matchsplitsitecluster_t *T, Uint ovhglen, fasta_t *set, Alignment*** aligns,Uint chromidx,Uint seqlen, char left, matchfileRec_t r,void *space, int *e2, char ***refseqs, Uint **reflens, Uint **refstrand, Uint locreflen, Uint median, Uint interval, Uint start,  Uint end, char *fwd, char *rev, Uint sinterval, int oven) { 
  int k, p, q, scores[] = {1,-1};; 
  Uint distpos, distchr,  distreflen,  l;  /* **reflens, **refstrand,*/
  int *M, **lmr, **lmv, **lmc, tempe, check;
  char disttrans;/*  ***refseqs,*/
  int overhang;
  int startleft;
  int check2; /*reinsert deletion penalties?*/
  overhang=2+ovhglen+ovhglen*.4;
  /* if(overhang>interval){
    overhang=interval;
    }*/
    startleft=(sinterval<overhang)? sinterval:overhang;
  /*  startleft=overhang;*/
    for(k=begin; k <= stop; k++) {
    e2[k]=0;
    distpos = T->distpos[k];/*->?*/
    distchr = T->distchr[k];
    disttrans = T->disttrans[k];
    distreflen = bl_fastaGetSequenceLength(set, distchr);    
    /*  mrgn = (distpos > interval) ? interval : distpos;  */
    refseqs[k] = ALLOCMEMORY(NULL, NULL, char*, 2);
    reflens[k] = ALLOCMEMORY(NULL, NULL, Uint, 2);
    refstrand[k] = ALLOCMEMORY(NULL, NULL, Uint, 2);
    
    aligns[k] = NULL;

    if(distchr != chromidx && (ovhglen < 10 || oven < ovhglen-(0.35*ovhglen))) continue; /*??always true?*/
    if(T->median+seqlen+1 >= locreflen) continue;
    if(distpos+overhang/*mrgn*/+1 >= distreflen) continue;
        
    if(left) { 
      if ((distchr == chromidx)&&(distpos>T->median) &&(T->median+seqlen/*??*/+1>distpos+ovhglen)) continue;
      if(r.strand=='+') {
	p =0; q =1;
      } else {
	p =1; q =0;
      }
      if ((distpos<startleft)) {
	startleft=distpos;
      }
      if (overhang>distpos){
	overhang=distpos;
      }
      refseqs[k][p] = &bl_fastaGetSequence(set, distchr)[distpos-overhang]; 
      reflens[k][p] = startleft+overhang;
      refseqs[k][q] = &bl_fastaGetSequence(set, chromidx)[T->median-1];
      reflens[k][q] =end-T->median+1 ;/*seqlen; */
      refstrand[k][q] = (r.strand == '+') ? 0 : 1;
      if(disttrans) 
	refstrand[k][p] = (r.strand == '+') ? 1 : 0;
      else
	refstrand[k][p] = (r.strand == '+') ? 0 : 1;  
    } else { 
      if ((distchr == chromidx)&&(distpos<T->median)&&(T->median-seqlen/*??*/-1<distpos-ovhglen)) continue; /*strand??*/
      if(r.strand=='+') {
	p =0; q =1;
      } else {
	p =1; q =0;
      }
      if (startleft>distpos) {
	startleft=distpos;
      }
      refseqs[k][p] = &bl_fastaGetSequence(set, chromidx)[start-1];
      reflens[k][p] = median+5;
      refseqs[k][q] = &bl_fastaGetSequence(set, distchr)[distpos-/*overhang*/startleft]; 
      reflens[k][q] = startleft+overhang;
      refstrand[k][p] = (r.strand == '+') ? 0 : 1;
      if(disttrans) 
	refstrand[k][q] = (r.strand == '+') ? 1 : 0;
      else
	refstrand[k][q] = (r.strand == '+') ? 0 : 1;
     }

    aligns[k] = ALLOCMEMORY(space, NULL, Alignment*, 2);  
    aligns[k][0] = ALLOCMEMORY(space, NULL, Alignment, 1);
    aligns[k][1] = ALLOCMEMORY(space, NULL, Alignment, 1);
    
    initAlignment(aligns[k][p], (refstrand[k ][p]==0)? fwd : rev, seqlen, 
		  0, refseqs[k ][p], reflens[k ][p], 0);
    initAlignment(aligns[k][q], (refstrand[k ][q]==0)? fwd : rev, seqlen, 
		  0, refseqs[k ][q], reflens[k ][q], 0);
    
    M = localmultisplicedmatrix(space, fwd, rev, seqlen,
				refseqs[k ], reflens[k ], refstrand[k ], 2, -1, 0, 
				constscr, scores, &lmv, &lmr, &lmc);
    
    localmultisplicedtraceback(space, M, fwd, rev, seqlen, 
			       refseqs[k ], reflens[k ], refstrand[k ], 2, -1, 0,
			       constscr, scores, aligns[k], lmv, lmr, lmc);
    e2[k]=getAlignScore(aligns[k][0], scores, -1);
    tempe= getAlignScore(aligns[k][1], scores, -1);
    /*how many bases are not aligned, add deletion scores*/
    check2=getUalignlen(aligns[k][q])+getUalignlen(aligns[k][p])-strlen(fwd);
    /*check whether at least 50% of the overhang are aligned*/
    if (left){
      check=getAlignScore(aligns[k][p], scores, -1);
      check*=2;
      check-=(seqlen-getUalignlen(aligns[k][q]));
      if(getUalignlen(aligns[k][p])<4) {
	check=-1;
      }
      //    if((getValignlen(aligns[k][q])+aligns[k][q]->voff)!=end) { /*v or u??*/
      //	check=-1;
      //}
    }
    else {
      check=getAlignScore(aligns[k][q], scores, -1)*2;
      check-=(seqlen-getUalignlen(aligns[k][p]));
      if(getUalignlen(aligns[k][q])<4) {
	check=-1;
      }
     // if (aligns[k][p]->voff != start) { /*v or u??*/
      //	check=-1;
      // }
    }
    /*if one of the parts is below 4??, set both zero*/
    if((check<0)||(e2[k] < 4) || (tempe < 4)) {
      e2[k]=0;/*in threading, e2=0?:*/ 
    }
    else {
      e2[k] += tempe;
      e2[k] += check2;
    }
    
    for(l=0; l < 2; l++) { 
      FREEMEMORY(space, lmv[l]);
      FREEMEMORY(space, lmr[l]);
      FREEMEMORY(space, lmc[l]);
    }
    FREEMEMORY(space, lmv);
    FREEMEMORY(space, lmr);
    FREEMEMORY(space, lmc);
    FREEMEMORY(space, M);
  }
  return NULL;
}


void Readrealign(int begin, int stop, matchlistelem_t* arr, unsigned char fmt,  matchsplitsitecluster_t *T, char left, fasta_t *set, Uint chromidx,Uint *noofrealigns/*??*/,FILE *realigndev,  List *rlist, Uint interval,void *space, Uint sinterval) {
  int j,p,k,l, rest, overhang;
  int scores[] = {1,-1};
  unsigned char *bookkeeper;
  Uint start;
  stringset_t *token;
  matchfileRec_t r;
  Alignment *laln, *raln, ***aligns;
  Uint bestdistalign=0, bestdistpos=0, bestlocalign =0, bestlocpos=0, bestlocstrand=0, 
    bestaligns=0, bestdistspliceoff=0, bestlocspliceoff=0, bestlocustart=0, bestlocuend=0, 
    bestdistustart=0, bestdistuend=0, bestdistchr=0, seqlen,locreflen=0; /*, distreflen=0,  q, v,*/
  Uint  distpos, distchr, ovhglen, median, end;
  int  bestscore, origscore, best;
  char *seq, *qual, *out, *rnext, *rm, *fwd, *rev, *str, *ref, ***refseqs,
       disttrans, bestdiststrand = '+', realigned = 0;
  int noofrealignsT=0, oven;
  Uint **reflens, **refstrand, mrgn;
  for(j=begin; j <= stop; j++) { 
      str = arr[j].str;
      start = arr[j].start;
      end = arr[j].end;
      bookkeeper = arr[j].bookkeeper; 

      if(!(bookkeeper[0] & LEFTSPLIT) && !(bookkeeper[0] & RIGHTSPLIT)) { 
	realigned = 0;
	token = tokensToStringset(space, "\t", str, strlen(str));
        bl_matchfileGetMatchFileRec(&r, 255, token, fmt); 
	//get overhanging sequence
        ref = &bl_fastaGetSequence(set, chromidx)[start-1]; 
	//both with 1-offset median points to position with 0-offset
        median = T->median - start; 
	//get boundaries and score of overhanging sequence
        rest = bl_matchfileSplitAlignment (r.curseq, ref, r.curaln, strlen(r.curaln), 
            median, scores, &laln, &raln, left);

        bestscore = getAlignScore(laln, scores, -1) + getAlignScore(raln, scores, -1);
        origscore = bestscore;

        if(left) {
          ovhglen = getUalignlen(laln) ;
          oven = getAlignScore(laln, scores, -1);
	  if(arr[j].end<T->median) {
	    ovhglen=0;
	  }
        } else { 
          ovhglen = getUalignlen(raln);    
          oven = getAlignScore(raln, scores, -1);      
	  if(arr[j].start>T->median) {
	    ovhglen=0;
	  }
        }
	
        seqlen = rest - laln->uoff;
        qual = ALLOCMEMORY(space, NULL, char, seqlen+1);
        seq = ALLOCMEMORY(space, NULL, char, seqlen+1);
        memmove(seq, &r.curseq[laln->uoff], seqlen);
        memmove(qual, &r.curqual[laln->uoff], seqlen);
        seq[seqlen] = 0;
        qual[seqlen] = 0;

        if(oven < ovhglen-(0.25*ovhglen) && ovhglen >= 4) { 
	  int *e2/*??*/;
	  int q;
          rm = charIUPACcomplement(space, seq, seqlen);
	  
/* 
          fprintf(stderr, "splitting: %s\n", str);
          fprintf(stderr, "ovhglen:%d, e:%d\n", ovhglen, e);
          fprintf(stderr, "distsites:%d\n", T[i].distsites);
          fprintf(stderr, "left:\n");
          showAlign(laln,stderr);
          fprintf(stderr, "right:\n");
          showAlign(raln,stderr);
          fprintf(stderr, "rco: %s\n", rm);
*/
          //align overhanging sequence to distant splits
          aligns = ALLOCMEMORY(space, NULL, Alignment**, T->distsites);
          refseqs = ALLOCMEMORY(space, NULL, char*, T->distsites);
          reflens = ALLOCMEMORY(space, NULL, Uint*, T->distsites);
          refstrand = ALLOCMEMORY(space, NULL, Uint*, T->distsites);
          locreflen = bl_fastaGetSequenceLength(set, chromidx);
	  e2= ALLOCMEMORY(space, NULL, int, T->distsites+1);
          //ensure original read direction
          if(r.strand == '-') {
            fwd = rm;
            rev = seq;
          } else {
            fwd = seq;
            rev = rm;
          }
	  assert(strlen(fwd) == seqlen);
	  
	  best=0; /*da??*/
	  /*e2???*/
	  
	  threadrealign(0,T->distsites-1 , T, ovhglen, set, aligns,chromidx,seqlen, left,r,space, e2, refseqs, reflens, refstrand, locreflen, median, interval, start, end, fwd, rev, sinterval,oven);
	  if(left) { 
	    if(r.strand=='+') {
	      p =0; q =1;
	    } else {
	      p =1; q =0;
	    }
	  } else { 
              if(r.strand=='+') {
                p =0; q =1;
              } else {
                p =1; q =0;
              }
	  }
	  
	  for (k=0; k< T->distsites; k++) {
	    if (e2[k]< bestscore) {
	      continue;
	    }
	    if (e2[k]> bestscore) {
	      if((left && aligns[k][q]->voff != 0) || 
		 (!left && aligns[k][p]->voff+getValignlen(aligns[k][p]) != median+1) ||
		 (left && ((T->median+aligns[k][q]->voff+getValignlen(aligns[k][q])) != end+1)) || /*?*/
		 (!left && aligns[k][p]->voff != 0) /*+1?*/ ) {
		//  fprintf(stderr, "not at median\n");
	      }
	      else {
		bestscore=e2[k];
		best=k;
	      }
	      continue;
	    }
	    if (T->distcnt[best]<T->distcnt[k]/*anzahl reads groesser*/) {
	      if((left && aligns[k][q]->voff != 0) || 
		 (!left && aligns[k][p]->voff+getValignlen(aligns[k][p]) != median+1) ||
		 (left && T->median+aligns[k][q]->voff+getValignlen(aligns[k][q]) != end+1) || /*???*/
		 (!left && aligns[k][p]->voff != 0)) {
		//  fprintf(stderr, "not at median\n");
	      }
	      else {
		bestscore=e2[k];
		best=k;
	      }
	    }
	  }
	 
	
	  if(bestscore > origscore) {
	    /*  if(e2[best] >= bestscore) { */

	    if((left && aligns[best][q]->voff != 0) || 
	       (!left && aligns[best][p]->voff+getValignlen(aligns[best][p]) != median+1)||
	       (left && T->median+aligns[best][q]->voff+getValignlen(aligns[best][q]) != end+1) || /*???*/
	       (!left && aligns[best][p]->voff != 0)) {
	   	  
	      //  fprintf(stderr, "not at median\n");
	    } 
	    else{
	      //  if(disttrans) fprintf(stderr, "distant split %d -> %d\n",refstrand[k][0], refstrand[k][1]);
	      //  else fprintf(stderr, "regular split %d -> %d\n", refstrand[k][0], refstrand[k][1]);
	      //  fprintf(stderr, "at median with score %d\n", e2);
	      bestaligns = best;
	      bestscore = e2[best];            
	      distpos=T->distpos[best]; 
	      distchr=T->distchr[best];
	      disttrans=T->disttrans[best];
	      overhang=2+ovhglen+ovhglen*.4;
	      mrgn=(distpos > overhang) ? overhang : distpos; 
	      if(left) { 
		bestdistalign = p; 
		bestlocspliceoff = 0;
		/*if(sinterval<mrgn){
		  mrgn=sinterval;
		  } this seems to be wrongly applied..*/
	      } else {  
		bestdistalign = q;
		if(sinterval<mrgn){
		  mrgn=sinterval;
		}
	      }
	     
	      bestdistpos = distpos-mrgn + aligns[bestaligns][bestdistalign]->voff + 1;
	      bestdistchr = distchr ;
	    
	      if (!disttrans) {
		bestdiststrand = r.strand;
	      } else {
		bestdiststrand = (r.strand == '+') ? '-' : '+';
	      }
	    
	      bestdistustart = (refstrand[bestaligns][bestdistalign] == 0) ? 
		aligns[bestaligns][bestdistalign]->uoff+1 : 
		seqlen-aligns[best][bestdistalign]->uoff - 
		getUalignlen(aligns[bestaligns][bestdistalign])+1;
	    
	      bestdistuend = bestdistustart+getUalignlen(aligns[bestaligns][bestdistalign])-1;
	    
	      if (left) { 
		bestlocalign = q; 
		bestlocpos = T->median;
	      
		bestdistspliceoff = (disttrans) ? 0 : 
		  getValignlen(aligns[bestaligns][bestdistalign])-1;

	      } else {  
		bestlocalign = p; 
		bestlocpos = start;
	      
		bestdistspliceoff = (!disttrans) ? 0 : 
		  getValignlen(aligns[bestaligns][bestdistalign])-1;

		bestlocspliceoff = getValignlen(aligns[bestaligns][bestlocalign]) - 1;
	      }

	      bestlocstrand = r.strand;
	    
	      bestlocustart = (refstrand[bestaligns][bestlocalign] == 0) ? 
		aligns[bestaligns][bestlocalign]->uoff+1 : 
		seqlen-aligns[bestaligns][bestlocalign]->uoff - 
		getUalignlen(aligns[bestaligns][bestlocalign])+1;

	      bestlocuend = bestlocustart+getUalignlen(aligns[bestaligns][bestlocalign])-1;
	      /*	  }*/
	      /*          fprintf(stderr, "aligns[0]:\n");
			  showAlign(aligns[k][0],stderr);
			  fprintf(stderr, "aligns[1]:\n");
			  showAlign(aligns[k][1],stderr);
	      */
	      /* }*/
	      /*
		FREEMEMORY(space, siteidx);
	      */
	      
	      /*  }*/ 
	      /*    if(bestscore > origscore) { */
	      bookkeeper[0] |= (left) ? LEFTSPLIT : RIGHTSPLIT;  
	      realigned = 1;
	      pthread_mutex_lock(&mutkuh);
	      
	      if(!T->realigned) {
		T->realigned = ALLOCMEMORY(space, NULL, Uint, T->distsites);
		memset(T->realigned, 0, sizeof(Uint)*T->distsites);
	      }
	      T->realigned[bestaligns]++;
	      /*    T->real++;*/
	      /*New for counting once only*/
	      /*	      T->cnt++;
			      T->distcnt[bestaligns]++;*/
	  
	      /*end new for counting once only*/
	      
	      pthread_mutex_unlock(&mutkuh);
	    
	      rnext = bl_fastaGetDescription(set, bestdistchr);
	      pthread_mutex_lock(&inout);
	    
	      out = bl_matchfileRealignWriteSAM (&r, qual, r.curchrom, bestlocpos, bestlocstrand, 
						 aligns[bestaligns][bestlocalign], bestlocustart, bestlocuend, 
						 (bestlocustart<bestdistustart?/*new*/0:1), 0, rnext, bestdistpos+bestdistspliceoff, bestdiststrand);
	      //	    pthread_mutex_lock(&inout);
	      fprintf(realigndev, "%s\n", out);
	      //	    pthread_mutex_unlock(&inout);
            
	      //          for(l=0; l < getValignlen(aligns[bestaligns][bestlocalign]); l++) {
	      //    fprintf(stderr, "%c",((char*)&bl_fastaGetSequence(set, chromidx)[bestlocpos-1])[l]);
	      //          }
	      //   fprintf(stderr, "\n");

	      FREEMEMORY(space, out);

	      out = bl_matchfileRealignWriteSAM (&r, qual, rnext, bestdistpos, bestdiststrand, 
						 aligns[bestaligns][bestdistalign], bestdistustart, bestdistuend, 
						 (bestlocustart<bestdistustart?/*new*/1:0)/*new*/, 0, r.curchrom, bestlocpos+bestlocspliceoff, bestlocstrand);

	      fprintf(realigndev, "%s\n", out);
	      pthread_mutex_unlock(&inout);

	      // for(l=0; l < getValignlen(aligns[bestaligns][bestdistalign]); l++) {
	      //   fprintf(stderr, "%c", ((char*)&bl_fastaGetSequence(set, bestdistchr)[bestdistpos-1])[l]);
	      // }
	      // fprintf(stderr, "\n");

	      FREEMEMORY(space, out);
	    }
	  }
	  FREEMEMORY(space, e2);
          for(k=0; k < T->distsites; k++) {
            if(aligns[k]) { 
              for(l=0; l < 2; l++) {  
                wrapAlignment(aligns[k][l]);
                FREEMEMORY(space, aligns[k][l]);
              }
              FREEMEMORY(space, aligns[k]);
            }
            FREEMEMORY(space, refseqs[k]);
            FREEMEMORY(space, reflens[k]);
            FREEMEMORY(space, refstrand[k]);
          }

          //enlist all unmarked sequences to be used in next iter
          FREEMEMORY(space, aligns);
          FREEMEMORY(space, refseqs);
          FREEMEMORY(space, reflens);
          FREEMEMORY(space, refstrand);
          FREEMEMORY(space, rm);

        }

        wrapAlignment(laln);
        wrapAlignment(raln);

        FREEMEMORY(space, laln);
        FREEMEMORY(space, raln);
        FREEMEMORY(space, r.curaln); 
        FREEMEMORY(space, r.diff);
        FREEMEMORY(space, seq);
        FREEMEMORY(space, qual);
        destructStringset(space, token);

        if(realigned) {
          noofrealignsT++;
        }
      } 
      
      pthread_mutex_lock(&mutkuh);

      bl_matchfileEnlistMatch (space, rlist, arr[j].start, arr[j].end, 
			       arr[j].distchr, arr[j].distpos, arr[j].adjoint, 
			       arr[j].trans, arr[j].str, arr[j].bookkeeper);
      pthread_mutex_unlock(&mutkuh);

     
  } 
  *noofrealigns=noofrealignsT;
  return;
}
      
void *Readthreadstarter(void *args) {
  readthread *t;
  t=(readthread*)args;
  Readrealign(t->begin, t->stop, t->arr, t->fmt, t->T,t->left,t->set,t->chromidx,&t->noofrealigns,t->realigndev,t->rlist, t->interval,t->space, t->sinterval);
  return NULL;
}


