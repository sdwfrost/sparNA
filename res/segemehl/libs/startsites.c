
/*
 *  startsites.c
 *  finding start sites
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 07/01/2011 09:53:24 PM CEST
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
#include "startsites.h"



/*----------------------------- bl_getStartSites -----------------------------
 *    
 * @brief getting the startsites
 * @author Steve Hoffmann 
 *   
 */
 

Uint
bl_matchfileStartSites(void *space, Uint fidx, Uint cidx, Uint pos, 
    matchfileCross_t *cs, char ref, matchfileSampleStats_t *stats, 
    unsigned char show, void *nfo) {

  Uint *cntr = (Uint *) nfo;

  if(cs->len > 5) {
    if(cs->starts < 255) cntr[(int)((double)(cs->starts*100.0)/cs->len)]++;
  }

  return 0;
}

Uint
bl_coverage(void *space, Uint fidx, Uint cidx, Uint pos, 
    matchfileCross_t *cs, char ref, matchfileindex_t *idx, 
    unsigned char show, void *nfo) {

  fprintf(stdout, "%d\n", cs->len);

  return 0;
}

