 #ifndef FALPHABET_H
 #define FALPHABET_H

/*
 * alphabet.h
 * declarations for a flexible alphabet
 *
 *  SVN
 *  Revision of last commit: $Rev: 40 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-09-03 10:42:48 +0200 (Wed, 03 Sep 2008) $
 *
 *  Id: $Id: falphabet.h 40 2008-09-03 08:42:48Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/trunk/libs/sufarray/falphabet.h $
 */

 #include "basic-types.h"

 typedef struct {

	Uint *characters,
		 *mapdomain;
	
	Uint domainsize,
		 mapsize,

		 mappedwildcards,
		 undefsymbol,
		 *symbolmap;

 } FAlphabet;


 /*from mapdomain to character*/
 Uint lookupChar(FAlphabet *, Uint);	
 void destructAlphabet(void *space, FAlphabet *);
 Uint cmp_map(FAlphabet*, Uint, Uint);	
 Uint lookupMapping(FAlphabet*, Uint);

#endif
