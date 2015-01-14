
#ifndef SUF_MATCH
#define SUF_MATCH

/*
 *
 *	sufmatch.h
 *  declarations for matching functions in suffixarrays
 * 
 *  @author Steve Hoffmann, shoffmann@zbh.uni-hamburg.de
 *  @company Center for Bioinformatics, Hamburg 
 *  @date 12/19/06 15:16:31 CET  
 *
 *  SVN
 *  Revision of last commit: $Rev: 19 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-05-14 15:43:29 +0200 (Wed, 14 May 2008) $
 *
 *  Id: $Id: sufmatch.h 19 2008-05-14 13:43:29Z steve $
 *  Url: $URL: file:///homes/bierdepot/steve/svn/segemehl/trunk/libs/sufarray/sufmatch.h $
 */

 #include "basic-types.h"
 #include "intsequence.h"
 #include "sufmatch.h"
 #include "sufarray.h" 

typedef struct {
	Uint id;
	Uint count;
	Uint *pos;
	Uint *org;
	Uint m;
	float score;
	float swscore;
	double blast;
} Matchtype;

PairSint* sufSubstring (void *, Suffixarray *, Uint *, Uint, Uint);
void reportSufmatch    (Suffixarray *, PairSint *, Uint, Uint, 
						IntSequence **);


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
);

Matchtype*
selectScoreSWconst(void *space, Matchtype *m, Uint k, 
	IntSequence *a, IntSequence **s, void *info);
 
Matchtype*
selectSW (void *space, Matchtype *m, Uint k, IntSequence *a, 
		IntSequence **s, void* info);

Matchtype*
selectScore (void *space, Matchtype *m, Uint k, IntSequence *a, 
		IntSequence **s, void* info);

double
swconstfilter(void *space, Matchtype *m, IntSequence *a, IntSequence *b,
					Uint *ptr, Uint len, Uint pos, void *info); 

double
scorefilter(void *space, Matchtype *m, IntSequence *a, IntSequence *b,
					Uint *ptr, Uint len, Uint pos, void *info);


Matchtype*
selectBlastScore (void *space, Matchtype *m, Uint k, 
		IntSequence *a, IntSequence **s, void* info);


Matchtype*
selectBlastScoreSWconst(void *space, Matchtype *m, Uint k, 
			IntSequence *a, IntSequence **s, void *info); 

#endif

