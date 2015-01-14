
/*
 *  matepairs.c
 *  
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 11/09/2011 05:16:34 PM CET
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
#include "matepairs.h"

/*------------------------- bl_matchfileInitMateMap --------------------------
 *    
 * @brief init the mate map
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchfileInitMateMap (void *space, matemap_t *map)
{
  map->noofmates = 0;
  map->refpos = NULL;
  map->refidx = NULL;
  map->materefpos = NULL;
  map->materefidx = NULL;
  map->matecount = NULL;

  return ;
}

/*----------------------- bl_matchfileDestructMateMap ------------------------
 *    
 * @brief destruct the mate map
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchfileDestructMateMap (void *space, matemap_t *map)
{

  if(map->refpos) FREEMEMORY(space, map->refpos);
  if(map->refidx) FREEMEMORY(space, map->refidx);
  if(map->materefpos) FREEMEMORY(space, map->materefpos);
  if(map->materefidx) FREEMEMORY(space, map->materefidx);
  if(map->matecount) FREEMEMORY(space, map->matecount);

  return ;
}

/*--------------------------- bl_matchfileMateScan ---------------------------
 *    
 * @brief get the mate locations for a map interval
 * @author Steve Hoffmann 
 *   
 */
 
matelink_t *
bl_matchfileMateScan (void *space, matchfileCross_t *cs)
{
  matelink_t *mates = NULL;



  return mates;
}

/*--------------------------- bl_matchfileAddMate ----------------------------
 *    
 * @brief add a mate to the cross section
 * @author Steve Hoffmann 
 *   
 */

void
bl_matchfileAddMate (void *space, matchfileCross_t *cs, Uint refidx, 
    Uint refpos)
{

  Uint i,k;

  for(i=0; i < cs->noofmatelinks; i++) {
    if(cs->matelinks[i].refpos == refpos && 
        cs->matelinks[i].refidx == refidx) {
      cs->matelinks[i].noofmates++;
      break;
    }
  }

  if(i == cs->noofmatelinks){ 
    k = cs->noofmatelinks;
    cs->matelinks = ALLOCMEMORY(space, cs->matelinks, matelink_t, k+1);
    cs->matelinks[k].refidx = refidx;
    cs->matelinks[k].refpos = refpos;
    cs->matelinks[k].noofmates = 1;
    cs->noofmatelinks++;
  }

  return ;
}


/*------------------------- bl_matchfileAddMateToMap -------------------------
 *    
 * @brief adding a mate to mate map
 * @author Steve Hoffmann 
 *   
 */

  void
bl_matchfileAddDistMateToMap (void * space, matemap_t *map, 
    matchfileCross_t *cs, Uint cidx, Uint pos, char ref)
{

  Uint k,i;
  k = map->noofmates;

  for(i=0; i < cs->noofmatelinks; i++) {  
    if(cs->matelinks[i].refidx != cidx ||       
        llabs(((Lint)cs->matelinks[i].refpos-(Lint)pos)) > MATEPAIRDIST) { 
   
      map->refpos = ALLOCMEMORY(space, map->refpos, Uint, k+1);
      map->refidx = ALLOCMEMORY(space, map->refidx, Uint, k+1);
      map->materefpos = ALLOCMEMORY(space, map->materefpos, Uint, k+1);
      map->materefidx = ALLOCMEMORY(space, map->materefidx, Uint, k+1);
      map->matecount = ALLOCMEMORY(space, map->matecount, Uint, k+1);

      map->refpos[k] = pos;
      map->refidx[k] = cidx;
      map->materefpos[k] = cs->matelinks[i].refpos;
      map->materefidx[k] = cs->matelinks[i].refidx;
      map->matecount[k] = cs->matelinks[i].noofmates;
      k++;
    }
  }

  map->noofmates = k;
  return ;
}


/*------------------------- bl_matchfileInitMateBin --------------------------
 *    
 * @brief init mate bin
 * @author Steve Hoffmann 
 *   
 */
 

void
bl_matchfileInitMateBin (void *space, matebin_t *bin, Uint pos, Uint ref)
{

  bin->binpos = pos;
  bin->binref = ref;
  bin->noofmates = 0;
  bin->matebinpos = NULL;
  bin->matebinref = NULL;
  bin->matebincnt = NULL;

  return ;
}

  
/*------------------------ bl_matchfileInitMateBinMap ------------------------
 *    
 * @brief init mate bin map
 * @author Steve Hoffmann 
 *   
 */
 

void
bl_matchfileInitMateBinMap (void *space, matebinmap_t *map)
{
  map->noofbins = 0;
  map->bins = NULL;

  return ;
}


/*------------------------- bl_matchfileScanMateBins -------------------------
 *    
 * @brief scan to find the bin at a given position
 * @author Steve Hoffmann 
 *   
 */
 
Uint
bl_matchfileScanBins (void *space, matebinmap_t *map, Uint pos, Uint ref)
{
  Uint binpos, startbinpos, endbinpos, i;
  char found = 0;

  binpos = startbinpos = pos - (pos % MATEBINSIZE);

  if (binpos > MATEBINSIZE) {
    binpos -= MATEBINSIZE;
  }
  endbinpos = binpos + MATEBINSIZE;

  
  for(i=0; i < map->noofbins; i++) {
    if(map->bins[i].binpos >= startbinpos && 
       map->bins[i].binpos <= endbinpos &&
       map->bins[i].binref == ref) {
      found = 1;        
      break;
    }
    
    if((map->bins[i].binpos > endbinpos &&  
        map->bins[i].binref == ref) ||
        map->bins[i].binref > ref) {
      break;
    }
  }


  if(found) return i;

  return -1;
}


/*-------------------------- bl_matchfileSearchMate --------------------------
 *    
 * @brief find a mate link from pos to dpos
 * @author Steve Hoffmann 
 *   
 */
 
Uint
bl_matchfileSearchMateLink (void *space, matebinmap_t *map, 
    Uint frompos, Uint fromref, Uint topos, Uint toref)
{
  Uint i, j, k, no=0, tobinpos;

  k = bl_matchfileScanBins(space, map, frompos, fromref);
  
  if (k == -1) {  
    return 0;
  } else if (k > MATEPAIRSEARCHBINMARGIN) { 
    k -= MATEPAIRSEARCHBINMARGIN; 
  } else if (k > 0) { 
    k = 0;
  } 

  tobinpos = topos - (topos % MATEBINSIZE);

  for(i=k; i < map->noofbins && i < k+MATEPAIRSEARCHBINMARGIN+1; i++) {
    for(j=0; j < map->bins[i].noofmates; j++) {
      if(map->bins[i].matebinpos[j] == tobinpos &&
          map->bins[i].matebinref[j] == toref) {
        no += map->bins[i].matebincnt[j];
      }
    }
  }

  return no;
}


/*------------------------- bl_matchfileAddToMateBin -------------------------
 *    
 * @brief adding a mate to mate map
 * @author Steve Hoffmann 
 *   
 */

  void
bl_matchfileAddToMateBinMap (void *space, matebinmap_t *map, 
    matchfileCross_t *cs, Uint cidx, Uint pos, char ref)
{

  Uint i, j, k, binpos, curbinpos = 0, curbinref = 0, 
       matebinpos = 0, matebinref = 0;
  
  k = map->noofbins;

  binpos = pos - (pos % MATEBINSIZE);
  curbinpos = (k > 0) ? map->bins[k-1].binpos : 0;
  curbinref = (k > 0) ? map->bins[k-1].binref : 0;

  if(k == 0 || curbinpos != binpos || cidx != curbinref) { 
    map->bins = ALLOCMEMORY(space, map->bins, matebin_t, k+1);
    bl_matchfileInitMateBin(space, &map->bins[k], binpos, cidx);
    map->noofbins++;
    k++;
  }

  for(i=0; i < cs->noofmatelinks; i++) {  

    matebinpos = 
      cs->matelinks[i].refpos - (cs->matelinks[i].refpos % MATEBINSIZE);
    matebinref = cs->matelinks[i].refidx;

    for(j=0; j < map->bins[k-1].noofmates; j++) {
       if(map->bins[k-1].matebinpos[j] == matebinpos &&
          map->bins[k-1].matebinref[j] == matebinref) { 
        break;
      }
    }

    if(j ==  map->bins[k-1].noofmates) {
      map->bins[k-1].matebinpos = 
        ALLOCMEMORY(space, map->bins[k-1].matebinpos, Uint, j+1);
      map->bins[k-1].matebinref = 
        ALLOCMEMORY(space, map->bins[k-1].matebinref, Uint, j+1);
      map->bins[k-1].matebincnt =
        ALLOCMEMORY(space, map->bins[k-1].matebincnt, Uint, j+1);
      
      map->bins[k-1].matebinpos[j] = matebinpos;
      map->bins[k-1].matebinref[j] = matebinref;
      map->bins[k-1].matebincnt[j] = 0;
      map->bins[k-1].noofmates++;
    } 
      
    map->bins[k-1].matebincnt[j]++;
  }

  return ;
}


/*---------------------- bl_matchfileDestructMateBinMap ----------------------
 *    
 * @brief destruct mate bin map
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchfileDestructMateBinMap (void *space, matebinmap_t *map)
{
  Uint i;
  for(i=0; i < map->noofbins; i++) {
    if(map->bins[i].matebinpos) {
      FREEMEMORY(space, map->bins[i].matebinpos);
      FREEMEMORY(space, map->bins[i].matebinref);
      FREEMEMORY(space, map->bins[i].matebincnt);
    }
  }
	
  if(map->bins) FREEMEMORY(space, map->bins);
  return ;
}

