
/*
 *
 *	bitvectoralg.h
 *  
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 05/26/2008 12:26:03 PM CEST  
 *
 */
#include "basic-types.h"
#include "alignment.h"
#include "bitVector.h"

PairSint myersbitvector( void *space,
    char *query, 
    Uint qlen, 
    char *subject, 
    Uint slen, 
    char *alphabet, 
    Uint asize,
    Uint *enctab,
    Uint k,
    bitvector *peq);

bitvector*
myersbitmatrix( void *space,
    char *query, 
    Uint qlen, 
    char *subject, 
    Uint slen, 
    char *alphabet, 
    Uint asize,
    Uint *enctab,
    Uint k,
    bitvector *peq,
    PairSint *res,
    bitvector *D,
    Uint dim);

bitvector*
myersblockbitmatrix(
    void *space,
    char *query, 
    Uint qlen, 
    char *subject, 
    Uint slen, 
    char *alphabet, 
    Uint asize,
    Uint *enctab,
    Uint k,
    bitvector *peq,
    PairSint *res,
    bitvector *D,
    Uint dim);

bitvector*
getpeq(void *space,
    char *query, 
    Uint qlen,
    char *alphabet,
    Uint asize,
    Uint *enctab);

Uint*
encodetab(char *alphabet, Uint asize) ;


char*
getstringalphabet (void *space, char *string, Uint len, Uint *asize);

Alignment*
bitvectorbacktrack(Alignment *al, bitvector *D, Uint dim, Uint k, Uint l);

