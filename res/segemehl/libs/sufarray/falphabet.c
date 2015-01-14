
/*
 * falphabet.c
 * implmentations for a flexible alphabet
 *
 *  SVN
 *  Revision of last commit: $Rev: 39 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-09-02 12:22:02 +0200 (Tue, 02 Sep 2008) $
 *
 *  Id: $Id: falphabet.c 39 2008-09-02 10:22:02Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/trunk/libs/sufarray/falphabet.c $
 */

 #include <stdlib.h>
 #include "basic-types.h"
 #include "falphabet.h"
 #include "sort.h"
 #include "memory.h"
 #include "debug.h"

 Uint lookupChar(FAlphabet* alphabet, Uint mapped) {
	Uint ch;
	
	ch = binarySearch(alphabet->mapdomain, alphabet->mapsize, &mapped, cmp_int_bin, NULL);
	return ch;
 	
 }

 void destructAlphabet (void *space, FAlphabet *alphabet) {
	
   	FREEMEMORY(space, alphabet->characters);
	FREEMEMORY(space, alphabet->mapdomain);
 	FREEMEMORY(space, alphabet);
 }

 Uint lookupMapping(FAlphabet* alphabet, Uint ch) {
        Uint i;
        i = binarySearch(alphabet->characters, alphabet->domainsize, &ch, cmp_int_bin, NULL);
        return i;     
 }

 Uint cmp_map(FAlphabet* f, Uint a, Uint b){
 	Uint amap = lookupMapping(f, a);
	Uint bmap = lookupMapping(f, b);
 	if (amap < 0 || amap >= f->mapsize || bmap < 0 || amap >= f->mapsize) {
		DBG("Found char could not be mapped. Exit forced.\n",NULL);
		exit(-1);
	}
	else
		return f->mapdomain[amap] == f->mapdomain[bmap];
 }
