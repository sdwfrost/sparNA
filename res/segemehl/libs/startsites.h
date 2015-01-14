
/*
 *
 *	startsites.h
 *  get the start sites
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 07/02/2011 01:00:18 AM CEST  
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include "sort.h"
#include "alignment.h"
#include "stringutils.h"
#include "basic-types.h"
#include "mathematics.h"
#include "matfile.h"
#include "bitVector.h"
#include "info.h"
#include "vtprogressbar.h"
#include "fileio.h"
#include "matchfilesfields.h"
#include "matchfiles.h"
#include "debug.h"
#include "evalmatchfiles.h"
#include "list.h"
#include "biofiles.h"

Uint
bl_matchfileStartSites(void *space, Uint fidx, Uint cidx, Uint pos, 
    matchfileCross_t *cs, char ref, matchfileSampleStats_t *stats, 
    unsigned char show, void *nfo);
Uint
bl_coverage(void *space, Uint fidx, Uint cidx, Uint pos, 
    matchfileCross_t *cs, char ref, matchfileindex_t *idx, 
    unsigned char show, void *nfo);

