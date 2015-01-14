#ifndef SEQCLIP_H
#define SEQCLIP_H

/*
 *
 *	seqclip.h
 *  
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 24.04.2010 22:32:24 CEST  
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include "basic-types.h"


Uint
bl_seqclipPolyA(void *space, char *sequence, Uint len, char* clp, Uint clen);

Uint
bl_seqclipSoft3Prime(void *space, char *sequence, Uint len, 
    char *toClip, Uint clipsize, Uint minclipacc, Uint pAlen);

Uint
bl_seqclipSoft5Prime(void *space, char *s, Uint len, 
    char *C, Uint clen, Uint minclipscr);

Uint
bl_seqclipHard3Prime(Uint len, Uint clipsize);

char*
bl_seqclipHard5Prime(char *s, Uint len, Uint clen);

char*
bl_seqclipFind3Prime (void *space, fasta_t *set, Uint samplesize, Uint fs, int ws);


#endif
