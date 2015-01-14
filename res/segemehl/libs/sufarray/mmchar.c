
/*
 *  mmchar.c
 *  implementations for searches manber myers style
 *  on enhanced suffix arrays (type char)
 *
 *  @author Steve Hoffmann
 *  @email shoffmann@zbh.uni-hamburg.de
 *  @date 12/22/06 19:13:27 CET
 *  
 *  SVN
 *  Revision of last commit: $Rev: 19 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-05-14 15:43:29 +0200 (Wed, 14 May 2008) $
 *
 *  Id: $Id: mmchar.c 19 2008-05-14 13:43:29Z steve $
 *  Url: $URL: file:///homes/bierdepot/steve/svn/segemehl/trunk/libs/sufarray/mmchar.c $
 */

 #include <stdlib.h>
 #include <stdio.h>
 #include "basic-types.h"
 #include "memory.h"
 #include "mathematics.h"
 #include "sufarray.h"
 #include "mmchar.h"

/*---------------------------------- mmleft ----------------------------------
 *    
 * part of the manber-myers pattern search algorithm
 * 
 */
 
int
mmleft(Suffixarray *arr, char *pattern, Uint len, int h, int l, int r) 
{
  	PairSint lbound, rbound, ibound;
	int mid;

	lbound = mmcompare (arr, pattern, len, l, h);
	if(lbound.b <= 0) {
		return l;
	}
	rbound = mmcompare (arr, pattern, len, r, h);
	if(rbound.b > 0) {
		return r+1;
	}	
	while (r > l+1) {
		mid = (l+r)/2;
		ibound = mmcompare (arr, pattern, len, mid, MIN(lbound.a, rbound.a));
		if(ibound.b <= 0) {
			rbound.a = ibound.a;
			r = mid;
		} else {
			lbound.a = ibound.a;
			l = mid;
		}
	}
	return r;
}

/*--------------------------------- mmright ----------------------------------
 *    
 * part of the manber-myers pattern search algorithm
 * 
 */
 
int
mmright (Suffixarray *arr, char *pattern, Uint len, int h, int l, int r)
{
  	PairSint lbound, rbound, ibound;
	int mid;
	
	lbound = mmcompare(arr, pattern, len, l, h);
	if (lbound.b < 0) {
		return -1; 
	}
	rbound = mmcompare(arr, pattern, len, r, h);
	if (rbound.b >= 0) {
		return r;
	}
	
	while (r > l+1) {
		mid = (l+r)/2;
		ibound = mmcompare(arr, pattern, len, mid, MIN(lbound.a, rbound.a));
		if (ibound.b >= 0 ) {
			lbound.a = ibound.a;
			l = mid;
		} else {
			rbound.a = ibound.a;
			r = mid;
		}
	}	
	return l;
}

/*-------------------------------- mmcompare ---------------------------------
 *    
 * manber 'n' myers compare for suffixarrays
 * 
 */
 
PairSint
mmcompare (Suffixarray *arr, char *pattern, Uint len, int idx, int start)
{
    char *sufptr;
	int t, margin;
	PairSint res;
	
	sufptr = &arr->seq->sequences[arr->suftab[idx]];
	t = start;
						/*length = totallength - relativeposition*/
	margin = MIN((len),(arr->seq->totallength-(sufptr-arr->seq->sequences)));

	while (t < margin) {
		if(pattern[t] < sufptr[t]) {
			res.a=t;
			res.b=-1;
			
			return res;
		} else {
			if(pattern[t] > sufptr[t]) {
				res.a=t;
				res.b=1;
				
				return res;
			} else {
				t++;
			}
		}
	}
	
	if (t == len) {
		  res.a=t;
		  res.b=0;
		  return res;
	} 
	
	res.a=t;
	res.b=-1;
	
	return res;
}


/*--------------------------------- mmsearch ---------------------------------
 *    
 * search manber'n'myers style
 * 
 */
 
PairSint
mmsearch (Suffixarray *arr, char *pattern, Uint len, Uint h, Uint l, Uint r)
{
    int p, q;
	PairSint ps;

	p = mmleft(arr, pattern, len, h, l, r);
	q = mmright(arr, pattern, len, h, l, r);

	if (p <= q) {
	    ps.a = p;
		ps.b = q;
		
		return ps;
	}
	
	ps.a = 1;
	ps.b = 0;
	
	return ps;
}


