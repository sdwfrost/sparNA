#ifndef MATEPAIRS_H
#define MATEPAIRS_H

/*
 *
 *	matepairs.h
 *  
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 11/09/2011 05:07:25 PM CET  
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "alignment.h"
#include "debug.h"
#include "stringutils.h"
#include "basic-types.h"
#include "mathematics.h"
#include "sort.h"
#include "matfile.h"
#include "bitVector.h"
#include "info.h"
#include "zran.h"
#include "nw.h"
#include "matchfiles.h"
#include "evalmatchfiles.h"
#include "manout.h"
#include "matchfilesfields.h"
#include "matchfiles.h"
#define MATEPAIRDIST 500
#define MATEBINSIZE 200
#define MATEPAIRSEARCHBINMARGIN 1

Uint bl_matchfileSearchMateLink (void *space, matebinmap_t *map, 
    Uint frompos, Uint fromref, Uint topos, Uint toref);
matelink_t * bl_matchfileMateScan (void *space, matchfileCross_t *cs);
void bl_matchfileAddMate (void *space, matchfileCross_t *cs, Uint refidx, Uint refpos);
void bl_matchfileInitMateMap (void *space, matemap_t *map);
void bl_matchfileDestructMateMap (void *space, matemap_t *map);
void bl_matchfileAddDistMateToMap (void * space, matemap_t *map, matchfileCross_t *cs,
    Uint cidx, Uint pos, char ref);
void bl_matchfileInitMateBinMap (void *space, matebinmap_t *map);
void bl_matchfileDestructMateBinMap (void *space, matebinmap_t *map);
void bl_matchfileAddToMateBinMap (void * space, matebinmap_t *map, 
    matchfileCross_t *cs, Uint cidx, Uint pos, char ref);



#endif
