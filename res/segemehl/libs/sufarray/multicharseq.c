
/*
 *  multiseq.c
 *  some functions to handle multiseqs (type char)
 *
 *  @author Steve Hoffmann
 *  @email shoffmann@zbh.uni-hamburg.de
 *  @date 12/15/06 11:42:53 CET
 *  
 *  SVN
 *  Revision of last commit: $Rev: 66 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-10-02 13:38:05 +0200 (Thu, 02 Oct 2008) $
 *
 *  Id: $Id: multicharseq.c 66 2008-10-02 11:38:05Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/sufarray/multicharseq.c $
 */

 #include <stdio.h>
 #include <stdlib.h>
 #include <string.h>
 #include "basic-types.h"
 #include "memory.h"
 #include "debug.h"
 #include "info.h"
 #include "charsequence.h"
 #include "vtprogressbar.h"
 #include "multicharseq.h"
 #include "alignment.h"
 #include "mathematics.h"
 #include "sort.h"
 
/*---------------------------- concatCharSequences ----------------------------
 *    
 * concatenates CharSequences using a given Uint delimiter
 * and stores them in a MultiCharSeq container.
 * 
 */
 
MultiCharSeq *
concatCharSequences (void *space, CharSequence **s, Uint len, 
					char delim, char sentinel)
{
    char *buf=NULL;
    char *map = NULL;
    Uint i, j, k=0, 
		 totallength=0, 
		 *markpos;
	MultiCharSeq *mseq;

	mseq = ALLOCMEMORY(space, NULL, MultiCharSeq, 1);
	markpos = ALLOCMEMORY(space, NULL, Uint, len);
	mseq->ref = ALLOCMEMORY(space, NULL, SeqReference, len);
    map = ALLOCMEMORY(space, NULL, char, 257);
    memset(map, 0, 256);
    mseq->delim = delim;

	for(i=0; i < len; i++) {

        mseq->ref[i].ref = s[i];
	  
        totallength += (s[i]->length+1);
		buf = ALLOCMEMORY(space, buf, char, totallength+1);
		if (buf==NULL) {
          NFO("allocation of %d bytes failed: exiting\n", totallength);
          exit(-1);
        }

		for(j=0; j < s[i]->length; j++) {
			buf[k] = s[i]->sequence[j];
                        if ((Uint)buf[k] == 0){
                          NFO("invalid character (NUL) in database sequences. Exit forced\n", NULL);
                          exit(-1);
                        }
            map[(Uint)buf[k]]=buf[k];
            k++;
		}
		/*separate sequences or finalize*/
		if (i == (len-1)) {
		  buf[k] = sentinel;
          map[(Uint)buf[k]]=buf[k];
		  markpos[i] = k;
		  k++;
          buf[k]='\0';

		} else {
		  buf[k] = delim;
		  map[(Uint)buf[k]]=buf[k];
          markpos[i] = k;
		  k++;
		}
        
        /*FREEMEMORY(space, s[i]->sequence);*/
	}
	mseq->totallength = totallength;
	mseq->numofsequences = len;
	mseq->sequences = buf;
	mseq->markpos = markpos;

    for(i=0; i < 256; i++) {
        if(map[i]==0) {
           j=i+1;
           while(j<256 && map[j]==0) j++;
           if (j < 256) {
             map[i]=map[j];
             map[j]=0;
           } else {
             break;
           }
        }
    }

    map = ALLOCMEMORY(space, map, char, i+1);
    mseq->map = map;
    mseq->mapsize = i;


	return mseq;
}



/*----------------------------- destructMultiSeq -----------------------------
 *    
 * destructs a MultiSeq structure
 * 
 */

void
destructMultiCharSeq (void *space, MultiCharSeq *mseq)
{
    
	FREEMEMORY(space, mseq->sequences);
	if (mseq->markpos != NULL)
      FREEMEMORY(space, mseq->markpos);
	if (mseq->map != NULL)
      FREEMEMORY(space, mseq->map);
    if (mseq->ref != NULL)
      FREEMEMORY(space, mseq->ref);
    FREEMEMORY(space, mseq);
	return ;
}


/*------------------------------- cmp_markpos --------------------------------
 *    
 * compare function for getMultiSeqIndex
 * 
 */
 
Uint
cmp_markpos (Uint a, void *data, void *key, void *info)
{
    Uint *d = (Uint*) data;
	Uint *k = (Uint*) key;
	
	if (d[a] > *k) {
		if (a > 0) {
			if (d[a-1] < *k) {
				return 0;
			} else {
				return 1;
			}
		} else {
			return 0;
		}
	}
	
    if (d[a] < *k) return 2;
	return 0;
}

/*-------------------------- getMultiSeqIndex --------------------------
 *    
 * returns index of a sequence in multiseq addressed by a pointer
 * 
 */
 
Uint
getMultiCharSeqIndex (MultiCharSeq *mseq, char *ptr)
{	
	Uint pos, i;
	
	if (mseq->numofsequences == 1){
	  return 0;
	}
	pos = (ptr - mseq->sequences); 
    if (mseq->numofsequences < MSEQ_BSEARCH_THRESHOLD) {
      i=binarySearch(mseq->markpos, mseq->numofsequences, &pos, 
          cmp_markpos, NULL);
    } else {
      for (i=0; i < mseq->numofsequences; i++) {
        if (mseq->markpos[i] > pos) break;
      }
    }

	return i;
}

/*------------------------- getMultiCharSeqIdxBounds -------------------------
 *    
 * @brief return start and end of idx 
 * @author Steve Hoffmann 
 *   
 */
 
void
getMultiCharSeqIdxBounds(MultiCharSeq *mseq, Uint idx, Uint *start, Uint *end)
{

  *start = (idx > 0) ? mseq->markpos[idx-1]+1 : 0;
  *end = mseq->markpos[idx]; 
  return ;
}


/*---------------------------- nextMultiSeqDelim -----------------------------
 *    
 * returns positions of next delimiter in multiseq (ie. end of current seq)
 * 
 */

Uint
nextMultiSeqDelim (MultiCharSeq *mseq, char *ptr)
{	
  return mseq->markpos[getMultiCharSeqIndex(mseq,ptr)];
}

 

/*---------------------------- getMultiSeqRelPos -----------------------------
 *    
 * returns the relative position of a pointer to multiseq
 * with respect to the addressed sequence.
 * 
 */
 
Uint
getMultiCharSeqRelPos (MultiCharSeq *mseq, char *ptr)
{
  Uint idx;
  CharSequence *seq;
  idx = getMultiCharSeqIndex(mseq, ptr);
  seq = getCharSequence(mseq, idx);
  return (ptr - seq->sequence);
}


/*------------------------------- dumpMultiSeq -------------------------------
 *    
 * dumps a multiseq to the screen
 * 
 */

void
dumpMultiCharSeq (MultiCharSeq *mseq)
{
  	Uint i;

	for(i=0; i < mseq->totallength; i++) {
		printf("%c-", mseq->sequences[i]);	
	}

	printf("\n");
	return ;
}


/*----------------------------- getCharSequence ------------------------------
 *    
 * @brief return the CharSequence at index idx
 * @author Steve Hoffmann 
 *   
 */
 
CharSequence*
getCharSequence (MultiCharSeq *mseq, Uint idx)
{
	return (CharSequence*) mseq->ref[idx].ref;
}


/*------------------------ initMultiCharSeqAlignment -------------------------
 *    
 * @brief initalize an alignment in the multichar seq
 *
 * @author Steve Hoffmann 
 *   
 */
 
int
initMultiCharSeqAlignment(
    void *space, MultiCharSeqAlignment* a, MultiCharSeq *seq, Uint pos, 
    Uint loff, Uint len, unsigned char strand, 
    char *qrydesc, char *query, Uint qrylen) 
{
  Uint sub_start, 
       sub_end;
  //Uint     i;

  a->subidx = getMultiCharSeqIndex(seq, &seq->sequences[pos]);
  getMultiCharSeqIdxBounds(seq, a->subidx, &sub_start, &sub_end);
  a->substart = sub_start;
  a->subend = sub_end;
   
  a->refstart = MAX(sub_start, (Lint)pos-loff);
  
  if(a->refstart > sub_end) {
    fprintf(stderr, "refstart > substart: skiping MultiCharSeqAlignment\n");
    return 0;
  }

  a->reflen = (sub_end > (Lint)a->refstart + len)?len:(sub_end - a->refstart)+1;
  a->refseq = &seq->sequences[a->refstart];
  a->refdesc = ((CharSequence*)seq->ref[a->subidx].ref)->description;
  a->qrydesc = qrydesc;
  a->query = query;
  a->strand = strand;
  a->qrylen = qrylen;
/*
  for(i=0; i < a->reflen; i++) {
    fprintf(stdout, "%c", a->refseq[i]);
  }
  fprintf(stdout, "\n");
*/  
  a->al = ALLOCMEMORY(space, NULL, Alignment, 1);
  initAlignment(a->al, query, qrylen, 0, a->refseq, a->reflen, 0);
	
  return 1;
}

void
wrapMultiCharSeqAlignment(void *space, MultiCharSeqAlignment *a) {
  wrapAlignment(a->al);
  FREEMEMORY(space, a->al);
}
/*----------------------- initMultiCharSeqAlignmentOpt -----------------------
 *    
 * @brief init a mcsa with query bounds
 * @author Steve Hoffmann 
 *   
 */
 
int
initMultiCharSeqAlignmentOpt(
    void *space, MultiCharSeqAlignment* a, MultiCharSeq *seq, Uint pos, 
    char *qrydesc, char *query, Uint start, Uint end, 
    Uint qrylen, Uint floff, Uint flen, Uint uloff, Uint uroff, Uint maxoff, unsigned char strand) 
{

   Uint sub_start, 
       sub_end,
       rlen, qstart, qend, qlen;
   
  //get bounds and length of reference sequence
  a->subidx = getMultiCharSeqIndex(seq, &seq->sequences[pos]);
  getMultiCharSeqIdxBounds(seq, a->subidx, &sub_start, &sub_end);
  a->substart = sub_start;
  a->subend = sub_end;
  a->refstart = MAX(sub_start, (Lint)pos-floff); //maxoff
  a->floff = pos - a->refstart;

  //this should not happen  
  if(a->refstart > sub_end) {
    fprintf(stderr, "refstart > substart: skiping MultiCharSeqAlignment\n");
    return 0;
  }

  rlen = flen;
  a->refseq = &seq->sequences[a->refstart];
  a->reflen = (sub_end > (Lint)a->refstart + rlen) ? rlen : (sub_end - a->refstart)+1;

 
  //get bounds and length of query sequence
  qstart = (start > maxoff+uloff) ? start-(maxoff+uloff) : 0;
  qend = (end + uroff + maxoff < qrylen) ? end + uroff + maxoff : qrylen;
  qlen = qend - qstart;
  a->query = query;
  a->qrystart = qstart;
  a->qrylen = qlen;
  a->strand = strand;

  //descriptions
  a->refdesc = ((CharSequence*)seq->ref[a->subidx].ref)->description;
  a->qrydesc = qrydesc;
  //init alignment and return
  a->al = ALLOCMEMORY(space, NULL, Alignment, 1);
  initAlignment(a->al, query, qlen, 0, a->refseq, a->reflen, 0);
	
  return 1;
}


