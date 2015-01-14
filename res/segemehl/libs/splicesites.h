#ifndef SPLICESITES_H
#define SPLICESITES_H

/*
 *
 *	splicesites.h
 *  
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 06/15/2011 11:49:01 PM CEST  
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "basic-types.h"
#include "matchfiles.h"
#include "biofiles.h"

void bl_matchfileInitSplitMap (void *space, splitmap_t *splitmap, annotationtrack_t *bed, matchfile_t **files, fasta_t*);
Uint bl_matchfileSplit (void *space, Uint fidx, Uint cidx, Uint pos, matchfileCross_t *cs, char ref, matchfileindex_t *idx, unsigned char show, void *nfo);
splicemap_t* bl_matchfileSpliceMap (void *space, splitmap_t *map, Uint interval, Uint minsplitno);
void printsplits (void *space, char *chr, Uint pos, matchfileCross_t *cs, splitmap_t *map);
void printsplice (void *space, splicemap_t *sm, FILE *out);
void printsplicebed (void *space, splicemap_t *sm, Uint minsplitno, char *title, FILE *out, FILE *transout);
void bl_matchfileDestructSpliceMap (void *space, splicemap_t *sm);
void bl_matchfileDestructSplitMap (void *space, splitmap_t *splitmap);
distsplitsites_t* bl_matchfileGetDistantSplitSites (void *space, matchfileCross_t *cs, Uint pos,
    Uint cidx, char type, Uint *noofsites, Uint *checkptr);
void bl_matchfileSpliceAnnotation (void *space, splicemap_t* sm, annotationtrack_t *track);

#endif
