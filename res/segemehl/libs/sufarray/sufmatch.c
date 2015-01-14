
/*
 *  sufmatch.c
 *  functions for matching in suffixarrays
 *
 *  @author Steve Hoffmann
 *  @email shoffmann@zbh.uni-hamburg.de
 *  @date 12/19/06 15:13:18 CET
 *  
 *  SVN
 *  Revision of last commit: $Rev: 19 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-05-14 15:43:29 +0200 (Wed, 14 May 2008) $
 *
 *  Id: $Id: sufmatch.c 19 2008-05-14 13:43:29Z steve $
 *  Url: $URL: file:///homes/bierdepot/steve/svn/segemehl/trunk/libs/sufarray/sufmatch.c $
 */


 #include <stdlib.h>
 #include <stdio.h>
 #include <string.h>
 #include <math.h>
 #include "basic-types.h"
 #include "memory.h"
 #include "mathematics.h"
 #include "sufarray.h"
 #include "sufmatch.h"
 #include "mm.h"
 #include "intsequence.h"
 #include "list.h"
 #include "depictseqs.h"
 #include "gnuplot_i.h"
 #include "dpalign.h"
 #include "cantor.h"
 #include "wurstimbiss.h"
 #include "falphabet.h"

/*------------------------------- suffixscore --------------------------------
 *    
 * scoring function for sw alignments of structure sequences 
 * 
 */
 
int
suffixscore (symtype a, symtype b, void* info)
{
   double* scores;

   scores = (double*) info;
   if (a == b) return 3;
   
	return -2;
}


/*------------------------------- sufSubstring -------------------------------
 *    
 * retrieve the longest substring matches in a suffixarray
 * 
 */
 
PairSint*
sufSubstring (void *space, Suffixarray *arr, Uint *pattern, 
			  Uint len, Uint sublen)
{
  	Uint i;
	PairSint *res, d;

	if (len <= sublen) 
	{
	  return NULL;
	}
	
	res = ALLOCMEMORY(space, NULL, PairSint, len-sublen);
	for(i=0; i < len-sublen; i++) 
	{
		d=mmsearch(arr, &pattern[i], sublen, 0, 0, arr->numofsuffixes-1);
		res[i].a=d.a;
	  	res[i].b=d.b;	  
	}
	
	return res;
}


/*------------------------------ reportSufmatch ------------------------------
 *    
 * returns a beautified match string
 * 
 */
 
void
reportSufmatch (Suffixarray *a, PairSint *matches, 
				Uint len, Uint threshold,
				IntSequence **s)
{
    Uint i, j, idx;
	char *info;

	for (i=0; i < len; i++) {
	   if (matches[i].b >= ((matches[i].a)+threshold)) {
		   /*valid matches*/
		  for (j=matches[i].a; j <= matches[i].b; j++) {	
             idx = getMultiSeqIndex(a->seq, &a->seq->sequences[a->suftab[j]]);
			 info = s[idx]->description;
			 printf("[%d]:\t %s\n", j, info);
		  }
	   }
    }
	
    return ;
}


/*---------------------------------- cmp_* -----------------------------------
 *    
 * compare functions for clibs's qsort and bsearch
 * 1. cmp_suffixno : used in rankSufmatch to sort sequence numbers
 * 2. cmp_ranks    : used in rankSufmatch to sort sequence ranks
 * 3. cmp_ranks_ptr: used in rankSufmatchList to sort sequence ranks
 * 
 */

int
cmp_blast(const void *a, const void *b) {
    Matchtype *first = (Matchtype*)a;
    Matchtype *second =(Matchtype*)b;	
    double frac_first, frac_second;

	
	frac_first = (double) first->blast ;
	frac_second = (double) second->blast ;

	if(frac_first > frac_second) return 1;
	if(frac_first < frac_second) return -1;
	
	return 0;
}


int
cmp_swscore(const void *a, const void *b) {
    Matchtype *first = (Matchtype*)a;
    Matchtype *second =(Matchtype*)b;	
    double frac_first, frac_second;

	if (first->swscore == 0 && second->swscore == 0) {
		if(first->count > second->count) return 1;
		if(first->count < second->count) return -1;
	}
	
	frac_first = (double) first->swscore ;
	frac_second = (double) second->swscore ;

	if(frac_first > frac_second) return 1;
	if(frac_first < frac_second) return -1;
	
	return 0;
}

int
cmp_score(const void *a, const void *b) {
    Matchtype *first = (Matchtype*)a;
    Matchtype *second =(Matchtype*)b;	
    double frac_first, frac_second;

	if (first->score == 0 && second->score == 0) {
		if(first->count > second->count) return 1;
		if(first->count < second->count) return -1;
	}
	
	frac_first = (double) first->score ;
	frac_second = (double) second->score;

	if(frac_first > frac_second) return 1;
	if(frac_first < frac_second) return -1;
	
	return 0;
}


int 
cmp_suffixno(const void *a, const void* b) {
    Matchtype *first = (Matchtype*)a;
    Matchtype *second =(Matchtype*)b;	
 
	
	if(first->id > second->id) return 1;
	if(first->id < second->id) return -1;

	return 0;
}


int 
cmp_ranks(const void *a, const void* b) {
    Matchtype *first = (Matchtype*)a;
    Matchtype *second =(Matchtype*)b;	
  
	if(first->count > second->count) return 1;
	if(first->count < second->count) return -1;

	return 0;
}

int 
cmp_rank_ptr(const void *a, const void* b) {
    Matchtype **first =  (Matchtype**)a;
    Matchtype **second = (Matchtype**)b;	
  
	if(first[0]->count > second[0]->count) return 1;
	if(first[0]->count < second[0]->count) return -1;

	return 0;
}

/*------------------------------ freeMatchtype -------------------------------
 *    
 * a function to delete Matchtypes in a list
 * 
 */
 
void
freeMatchtype (void *space, void *data)
{
    Matchtype *d = (Matchtype*)data;

	FREEMEMORY(space, d->pos);
	FREEMEMORY(space, data);
}


/*---------------------------------- subscr ----------------------------------
 *    
 * a function that assigns scores from a substitution matrix for matches 
 * and mismatches given in info matrix
 * 
 */
 
int
subscr (symtype a, symtype b, void *info )
{
  	double* M;	
	FAlphabet *alphabet;
	imbissinfo *p;
	int val;
	
	
	p = (imbissinfo*) info;
	alphabet = p->alphabet;
	M = (double*) p->sub;

/*	ad = alphabet->mapdomain[(Uint)a];
	bd = alphabet->mapdomain[(Uint)b];
	printf("decoding %d and %d", ad, bd);
	tupela = decodeCantor(NULL, ad, 2);
	tupelb = decodeCantor(NULL, bd, 2);	
	printf("... ended ... to %d and %d \n", VECTOR(tupela,0), VECTOR(tupelb,0));
	val =  MATRIX2D(M, 309, VECTOR(tupela, 0), VECTOR(tupelb,0));
*/
	val = MATRIX2D(M, 308, a ,b);
	if (val < -10) val = -10;
	
	return val;
}



/*-------------------------------- getEntropy --------------------------------
 *    
 * calculates the entropy of a sequence, given probabilities.
 * 
 */
 
double
getEntropy(void *space, Uint* sequence, Uint l, double* prob) {
	int i;
	double sum=0;
	
	for(i=0; i < l; i++) {
		sum += prob[sequence[i]]*log2(prob[sequence[i]]);
	}

	return sum;
}


			/*local alignment*/
			/*printf("max sw score: %f\n", occ[i-1].swscore);
			
			
	  		swres = swmatrix(space, queryseq->sequence, queryseq->length, 
				s[occ[i-1].id]->sequence, s[occ[i-1].id]->length, 
				-5, constscr, swscores);	    
		
			
			printf("max sw score: %d\n", swres[arraymax(swres, 
	  			(queryseq->length+1)*(s[occ[i-1].id]->length+1))]);
			
			align = swgaplesstraceback(space, swres,  
					queryseq->sequence, queryseq->length, 
					s[occ[k].id]->sequence, s[occ[k].id]->length,
					//suffixscore,((imbissinfo*)info)->score 
					-5, 
					constscr, swscores,
					&alignsize);
			 
			if (depictsw) {
			   	alignstr = printAlignment(space, align, alignsize,
				  		queryseq, s[occ[i-1].id], 80);
				printf("%s\n", alignstr);
				FREEMEMORY(space, alignstr);
				FREEMEMORY(space, align); 
			}*/


Matchtype*
selectBlastScoreSWconst(void *space, Matchtype *m, Uint k, 
			IntSequence *a, IntSequence **s, void *info) {
 
  	Uint l, i;
	int *swres;
	imbissinfo *imbiss;

	imbiss = (imbissinfo*) info;
		
	qsort(m, k, sizeof(Matchtype), cmp_blast);
	
	l=0;
	for (i=k; i > 0 && l < 1000; i--) {
		if (m[i-1].count >= imbiss->minseeds) {
					
					swres = swgapless(space, 
							a->sequence, a->length, 
							s[m[i-1].id]->sequence, s[m[i-1].id]->length, 
							constscr, imbiss->swscores
							/*subscr, info*/
						);	    
				
					m[i-1].swscore= swres[arraymax(swres, 
	  						(a->length+1)*(s[m[i-1].id]->length+1))];
	
			 		FREEMEMORY(space, swres);
		} else {
			m[i-1].swscore = 0;
		}
		l++;
    }

	qsort(m, k, sizeof(Matchtype), cmp_swscore);
	return m;
}


Matchtype*
selectScoreSWconst(void *space, Matchtype *m, Uint k, 
			IntSequence *a, IntSequence **s, void *info) {
 
  	Uint l, i;
	int *swres;
	imbissinfo *imbiss;

	imbiss = (imbissinfo*) info;
		
	qsort(m, k, sizeof(Matchtype), cmp_score);
	
	l=0;
	for (i=k; i > 0 && l < 1000; i--) {
		if (m[i-1].count >= imbiss->minseeds) {
					
					swres = swgapless(space, 
							a->sequence, a->length, 
							s[m[i-1].id]->sequence, s[m[i-1].id]->length, 
							constscr, imbiss->swscores
							/*subscr, info*/
						);	    
				
					m[i-1].swscore= swres[arraymax(swres, 
	  						(a->length+1)*(s[m[i-1].id]->length+1))];
	
			 		FREEMEMORY(space, swres);
		} else {
			m[i-1].swscore = 0;
		}
		l++;
    }

	qsort(m, k, sizeof(Matchtype), cmp_swscore);
	return m;
}


Matchtype*
selectSW (void *space, Matchtype *m, Uint k, 
		IntSequence *a, IntSequence **s, void *info) {
  qsort(m, k, sizeof(Matchtype), cmp_swscore);	
  return m;
}

Matchtype*
selectBlastScore (void *space, Matchtype *m, Uint k, 
		IntSequence *a, IntSequence **s, void* info) {
  
  qsort(m, k, sizeof(Matchtype), cmp_blast);
  return m;
}

Matchtype*
selectScore (void *space, Matchtype *m, Uint k,
		IntSequence *a, IntSequence **s, void* info) {

  qsort(m, k, sizeof(Matchtype), cmp_score);
  return m;
}


double
scorefilter (void *space, Matchtype *m, IntSequence *a, IntSequence *b,
					Uint *ptr, Uint len, Uint pos, void *info) {
  Uint l;
  double temp = 0;
  double sum = 0; 
  imbissinfo *imbiss;
  
  imbiss=(imbissinfo*) info;
  
  m->count++;
  m->pos = ALLOCMEMORY(space, m->pos, Uint, m->count);
  m->org = ALLOCMEMORY(space, m->org, Uint, m->count);
  m->pos[(m->count)-1]=pos;
  m->org[(m->count)-1]=pos;
					
  for (l=0; l < len; l++){
	temp = ((imbissinfo*)info)->score[(Uint)*ptr];
	sum += temp;
	m->score += temp;	
	ptr++;
  }
  
  m->blast = m->blast > sum ? m->blast : sum;

  imbiss->consensus[pos] += (Uint) 1
					/*((double)imbiss->lambda*sum)*/;
  
  return sum > 0 ? sum : 0;
}



double
swconstfilter(void *space, Matchtype *m, IntSequence *a, IntSequence *b,
					Uint *ptr, Uint len, Uint pos, void *info) {
 	
    imbissinfo *imbiss; 
	int *swres;
	double t;

	imbiss = (imbissinfo*) info;
	t=scorefilter(space, m, a, b, ptr, len, pos, info);
	
	if (m->count == imbiss->minseeds) {
	  	swres = swgapless(space, a->sequence, a->length, 
			  	b->sequence, b->length, 
			  	constscr, imbiss->swscores);	    
		m->swscore= swres[arraymax(swres, (a->length+1)*(b->length+1))];
		
		FREEMEMORY(space, swres);
	}
	
	return t;	
}



/*------------------------------ initMatchtype -------------------------------
 *    
 * init a Matchtype struct
 * 
 */
 
void
initMatchtype (Matchtype *m, Uint id)
{
   
	m->id = id;
	m->count = 0;
	m->pos = NULL;
	m->org = NULL;
	m->m=0;
	m->blast=0;
	m->score = 0;
	m->swscore = 0;
	
	return ;
}

/*-------------------------------  rankSufmatch  ------------------------------
 *    
 * ranks matches 
 * given in an array of PairSint of length len. Sorting is done
 * by several calls to clib's qsort. For each item of the sorted
 * array a handler-function is invoked.
 * 
 */

void
rankSufmatch ( 	void *space, 
				Suffixarray *a, 
				PairSint *matches, 
				Uint len, 	
				Uint maxmatches, 
				Uint S, 
				IntSequence **s,
				Uint noofseqs,
				double (*fltr)(		void *, 
				  					Matchtype *, 
									IntSequence *, 
									IntSequence *, 
	  				        		Uint *, 
									Uint, 
									Uint, 
									void *),
				Matchtype* (*sel)(	void *, 
				  					Matchtype *, 
									Uint, 
									IntSequence *,
									IntSequence **,
									void *),
				int (*handler)(	void *, 
				  					Matchtype *, 
									IntSequence **, 
									Uint, 
									Uint, 
									void *), 	
				IntSequence *queryseq, 
				void *info, 
				double *scores, 
				unsigned char depictsw
)

{
  
	Matchtype key, *cur, *occ=NULL;
	int i=0, j=0, k=0, l=0, r=0, retval=0;	
	int *hashTable;
	double t;
	Uint *ptr;	

	hashTable = ALLOCMEMORY(space, NULL, int, (noofseqs+1));
	memset(hashTable, -1, sizeof(int)*(noofseqs+1));
	  	
	for(i=0; i < len; i++) {
		if(matches[i].b >= ((matches[i].a))) {
		
		  	for(j=matches[i].a; j <= matches[i].b; j++) {
		        key.id = getMultiSeqIndex(a->seq, &a->seq->sequences[a->suftab[j]]);	
				r = hashTable[key.id];
					
				if (r == -1) {
				    occ = ALLOCMEMORY(space, occ, Matchtype, k+1);
					hashTable[key.id]=(&occ[k])-(occ);
					initMatchtype(&occ[k], key.id);
					cur = &occ[k];
					k++; 
				} else {
					cur = ((Matchtype*)(occ+r));
				}
				
				/*score the matches if no < maxmatches*/	
				if ((matches[i].b-matches[i].a) < maxmatches) {
					ptr = &a->seq->sequences[a->suftab[j]]);
					t=fltr(space, cur, queryseq, s[cur->id], ptr, S, i, info);
					
					if (t == -1) break;
				}	
			}
		}
	}
	
	occ = sel(space, occ, k, queryseq, s, info);

	l=0;
	for (i=k; i > 0; i--) {
			retval = handler(space, &occ[i-1], s, len, l, info); 
			if (retval) l++;
			if (retval == -1) break;
	}
   
	FREEMEMORY(space, hashTable);

	for(i=0; i < k; i++) {
		FREEMEMORY(space, occ[i].pos);
		FREEMEMORY(space, occ[i].org);
	}

	FREEMEMORY(space, occ);    
}

