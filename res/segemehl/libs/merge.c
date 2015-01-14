/*
 * merge.c
 * functions to merge matches
 *
 *  SVN
 *  Revision of last commit: $Rev: 348 $
 *  Author: $Author: steve $
 *  Date: $Date: 2012-08-24 08:46:52 -0400 (Fri, 24 Aug 2012) $
 *
 *  Id: $Id: merge.c 348 2012-08-24 12:46:52Z steve $
 *  Url: $URL: http://www2.bioinf.uni-leipzig.de/svn5/segemehl/libs/merge.c $
 *  
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "debug.h"
#include "info.h"
#include "basic-types.h"
#include "stringutils.h"
#include "mathematics.h"
#include "biofiles.h"
#include "fileBins.h"
#include "matchfilesfields.h"
#include "merge.h"



/*------------------------- bl_mergefilesInit ---------------------------------
 *    
 * @brief init container for multiple merge files
 * @author Christian Otto
 *   
 */
void bl_mergefilesInit(void *space, bl_mergefiles_t *files, Uint nooffiles){
    files->f = ALLOCMEMORY(space, NULL, bl_mergefile_t, nooffiles);
    files->nooffiles = nooffiles;
}


/*----------------------- bl_mergefilesDestruct -------------------------------
 *    
 * @brief destruct container for multiple merge files
 * @author Christian Otto
 *   
 */
void bl_mergefilesDestruct(void *space, bl_mergefiles_t *files){
  if (files->f != NULL){
    FREEMEMORY(space, files->f);
    files->f = NULL;
  }
  files->nooffiles = 0;
}

/*--------------------------- bl_mergefileInit --------------------------------
 *    
 * @brief init merge file container
 * @author Christian Otto
 *   
 */
void bl_mergefileInit(void *space, bl_mergefile_t *file, FILE *fp){
  file->fp = fp;
  file->eof = 0;
  file->complete = 0;
  file->entry = ALLOCMEMORY(space, NULL, bl_mergefilematch_t, 1);
  bl_mergefilematchInit(space, file->entry);
}


/*------------------------ bl_mergefileDestruct -------------------------------
 *    
 * @brief destruct merge file container
 * @author Christian Otto
 *   
 */
void bl_mergefileDestruct(void *space, bl_mergefile_t *file){
  file->fp = NULL;
  file->eof = 0;
  file->complete = 0;
  bl_mergefilematchDestruct(space, file->entry);
  FREEMEMORY(space, file->entry);
  file->entry = NULL;
}


/*------------------------ bl_mergefilematchInit ------------------------------
 *    
 * @brief init merge file entry
 * @author Christian Otto
 *   
 */
void bl_mergefilematchInit(void *space, bl_mergefilematch_t *entry){
  entry->match = NULL;
  entry->matchlen = 0;
  entry->matematch = NULL;
  entry->matematchlen = 0;  
  entry->key = NULL;
  entry->keylen = 0;
  entry->flag = 0;
  entry->mateflag = 0;
  entry->edist = 0;
  entry->mateedist = 0;
  entry->noofmatches = 0;
  entry->noofmatematches = 0;
  entry->rname = NULL;
  entry->rstart = 0;
}


/*---------------------- bl_mergefilematchDestruct ----------------------------
 *    
 * @brief destruct merge file entry
 * @author Christian Otto
 *   
 */
void bl_mergefilematchDestruct(void *space, bl_mergefilematch_t *entry){
  if (entry->match != NULL){
    FREEMEMORY(space, entry->match);
    entry->match = NULL;
  }
  entry->matchlen = 0;
  if (entry->matematch != NULL){
    FREEMEMORY(space, entry->matematch);
    entry->matematch = NULL;
  }
  entry->matematchlen = 0;
  if (entry->key != NULL){
    FREEMEMORY(space, entry->key);
    entry->key = NULL;
  }
  entry->keylen = 0;
  entry->flag = 0;
  entry->mateflag = 0;
  entry->edist = 0;
  entry->mateedist = 0;
  entry->noofmatches = 0;
  entry->noofmatematches = 0;
  if (entry->rname != NULL){
    FREEMEMORY(space, entry->rname);
    entry->rname = NULL;
  }
  entry->rstart = 0;
}

/*------------------------- bl_mergefileCompare -------------------------------
 *    
 * @brief compare two merge file entries regarding SAM flag in case of
 *        paired-end reads (i.e. best pairing state), returns -1 if
 *        first is better, 1 if second is better, 0 otherwise @author
 *        Christian Otto
 *   
 */
int bl_mergefilematchCompareFlag(bl_mergefilematch_t *i, bl_mergefilematch_t *j){
  int tmpi, tmpj;

  /* flag compare */
  if (i->flag & 1) {
    /* if qry/mate unmapped */
    tmpi = (i->flag >> 3) & 1;
    tmpj = (j->flag >> 3) & 1;
    if (tmpi != tmpj){
      return tmpi - tmpj;
    }
    /* if proper pair */
    tmpi = ((i->flag >> 1) & 1);
    tmpj = ((j->flag >> 1) & 1);
    if (tmpi != tmpj){
      return -1 * (tmpi - tmpj);
    }
  }
  return 0;
}

/*---------------------- bl_mergefileCompareEdist -----------------------------
 *    
 * @brief compare two merge file entries regarding edit distance (or
 *        pair edit distance in case of paired-end reads), returns -1
 *        if first is better, 1 if second is better, 0 otherwise
 *        @author Christian Otto
 *   
 */
int bl_mergefilematchCompareEdist(bl_mergefilematch_t *i, bl_mergefilematch_t *j){
  int tmpi, tmpj;

  /* edist compare */
  tmpi = i->edist + i->mateedist;
  tmpj = j->edist + j->mateedist;
  if (tmpi != tmpj){
    if (tmpi < tmpj){
      return -1;
    }
    else {
      return 1;
    }
  }
  else {
    return 0;
  }
}


/*------------------------- bl_mergefileFastaIDCompare -------------------------------
 *    
 * @brief compare two fasta descriptions if they contain the same fasta ID,
 *        in case of paired-end data, it disregards /1 or /2 at the end or
 *        any  differences after the first white space
 *        returns 1 if both descriptions contain the same ID,  0 otherwise
 * @author Christian Otto
 *   
 */
unsigned char
bl_mergefileFastaIDCompare(char *desc1, Uint desc1len, char *desc2, Uint desc2len) {

  char *id, *id2, *tok1, *tok2;
  unsigned char res;

  id = ALLOCMEMORY(space, NULL, char, desc1len+2); 
  id2 = ALLOCMEMORY(space, NULL, char, desc2len+2); 

  strcpy(id, desc1);
  strcpy(id2, desc2);

  tok1 = strtok(id, "/");
  tok2 = strtok(id2, "/");
  res = (strcmp(tok1, tok2)==0);

  if(!res) { 
    FREEMEMORY(space, id);
    FREEMEMORY(space, id2);

    id = ALLOCMEMORY(space, NULL, char, desc1len+2); 
    id2 = ALLOCMEMORY(space, NULL, char, desc2len+2); 

    strcpy(id, desc1);
    strcpy(id2, desc2);

    tok1 = strtok(id, " ");
    tok2 = strtok(id2, " ");
    res = (strcmp(tok1, tok2)==0);
  }

  FREEMEMORY(space, id);
  FREEMEMORY(space, id2);
  return res;
}


/*------------------------- bl_mergeParseLine ---------------------------------
 *    
 * @brief parses a SAM-formatted line (single or paired-end) and 
 *        inserts the data in the given container
 *        NOTE: split reads not supported up to now
 * @author Christian Otto
 *   
 */
unsigned char bl_mergeParseLine(void *space, bl_mergefilematch_t *entry, char *line, Uint len){
  unsigned char complete = 0;
  char *tmp;
  Uint tmplen, flag;
  stringset_t *fields = NULL;

  fields = tokensToStringset(space, "\t\007", line, len);
  flag = bl_matchfileGetFlag(fields, SAM);

  if ((flag >> 8) & 1){
    DBG("Split reads not supported yet. Exit forced.\n", NULL);
    exit(-1);
  }
  /* ignore unmapped fragments */
  if ((flag >> 2) & 1){    
    destructStringset(space, fields);
    return complete;
  }

  /* 
   * reading match
   * (simply first fragment in
   * output, not necessarily
   * the read match)
   */
  if (entry->match == NULL){
    entry->match = line;
    entry->matchlen = len;

    tmp = bl_matchfileGetQname(fields, SAM);
    if (tmp == NULL){
      DBG("Error in parsing line. Exit forced.\n", NULL);
      exit(-1);
    }
    entry->keylen = strlen(tmp);
    entry->key = ALLOCMEMORY(space, NULL, char, entry->keylen + 1);
    memmove(entry->key, tmp, entry->keylen);
    entry->key[entry->keylen] = '\0';
    
    entry->flag = flag;
    entry->edist = bl_matchfileGetEdist(fields, SAM);
    entry->noofmatches = bl_matchfileGetMatchCnt(fields, SAM);
    
    tmp = bl_matchfileGetChrom(fields, SAM);
    tmplen = strlen(tmp);
    entry->rname = ALLOCMEMORY(space, NULL, char, tmplen + 1);
    memmove(entry->rname, tmp, tmplen);
    entry->rname[tmplen] = '\0';
    entry->rstart = bl_matchfileGetStartPos(fields, SAM);

    /* abort if non-valid flags if paired (either first or second in pair) */
    if ((entry->flag & 1) && !(((entry->flag >> 6) & 1) ^ 
			       ((entry->flag >> 7) & 1))){
      DBG("Incorrect flag information in paired-end match. Exit forced.\n", NULL);
      exit(-1);
    }

    /* match complete if unpaired or unmapped other fragment */
    if (!(entry->flag & 1) ||
	((entry->flag >> 3) & 1)){
      complete = 1;
    }
  }
  /* 
   * reading mate match 
   * (simply second fragment in
   * output, not necessarily
   * the mate match)
   */
  else {
    entry->matematch = line;
    entry->matematchlen = len;


    /* abort if mateflag already set */
    if (entry->mateflag & 1){
      DBG("Error in reading paired-end matches. Exit forced.\n", NULL);
      exit(-1);
    }
    entry->mateflag = flag;

    /* abort if not equal query name */
    tmp = bl_matchfileGetQname(fields, SAM);
    if (tmp == NULL){
      DBG("Error in parsing line. Exit forced.\n", NULL);
      exit(-1);
    }
    if (! bl_mergefileFastaIDCompare(entry->key, entry->keylen, tmp, strlen(tmp))){
      DBG("Error in reading paired-end matches. Exit forced.\n", NULL);
      exit(-1);
    }

    /* 
     * abort with non-valid flags
     * (mate unpaired, both/none first/second in pair)
     */
    if (!(entry->mateflag & 1) ||
	!(((entry->flag >> 6) & 1) ^ ((entry->mateflag >> 6) & 1)) ||
	!(((entry->flag >> 7) & 1) ^ ((entry->mateflag >> 7) & 1))) {
      DBG("Incorrect flag information in paired-end matches. Exit forced.\n", NULL);
      exit(-1);
    }
    entry->mateedist = bl_matchfileGetEdist(fields, SAM);
    entry->noofmatematches = bl_matchfileGetMatchCnt(fields, SAM);

    complete = 1;
  }
  destructStringset(space, fields);
  return complete;
}


/*-------------------------- bl_mergeReadNext ---------------------------------
 *    
 * @brief read next match (and mate match) entry in mergefile
 * @author Christian Otto
 *   
 */
void bl_mergeReadNext(void *space, bl_mergefile_t *file){
  Uint buffersize = 1024, len = 0;
  char *buffer = NULL, ch;

  if (!file->complete && !file->eof){
    buffer = ALLOCMEMORY(space, NULL, char, buffersize);
    len = 0;
    while((ch = getc(file->fp)) != EOF) {
      /* extend buffer */
      if(len == buffersize-1) {
	buffersize = 2*buffersize+1;
	buffer = ALLOCMEMORY(space, buffer, char, buffersize);
      }
      /* process buffer */
      if(ch == '\n' && len > 0) {
	buffer[len++] = ch;
	buffer = ALLOCMEMORY(space, buffer, char, len+1);  
	buffer[len] = '\0';

	file->complete = bl_mergeParseLine(space, file->entry, buffer, len);
	  
	if (file->complete){
	  break;
	}
	else {
	  buffer = ALLOCMEMORY(space, NULL, char, buffersize);
	  len = 0;
	  continue;
	}
      } else {
	if(ch != '\n') buffer[len++] = ch;
      }
    }
    /* set end of file */
    if (!file->eof && ch == EOF){	
      file->eof = 1;

      if (len > 0 && !file->complete){
	DBG("%u:%s\n", len, buffer);
	DBG("Incomplete read matching entry at end of file. Exit forced.\n", NULL);
	exit(-1);
      }
    }
    if (!file->complete){
      FREEMEMORY(space, buffer);
    }
  }
}


/*------------------------- bl_mergeUpdateTag ---------------------------------
 *    
 * @brief replaces noofmatches in SAM tag in given input line,
 *        input line and length is given as references
 * @author Christian Otto
 *   
 */
void bl_mergeUpdateTag(void *space, char **line, Uint *len, Uint noofmatches){
  Uint i, totallen;
  char *pch, *res;

  /* search for NH:i: tag in read match */
  pch = strstr(*line, "NH:i:");
  if (pch == NULL){
    DBG("Error in updating NH TAG. Exit forced.\n", NULL);
    exit(-1);
  }
  
  /* find end of tag */
  for (i = 1; i < strlen(pch); i++){
    if (ISWHITESPACE(pch[i]) || pch[i] == '\007'){
      break;
    }
  }
  pch[0] = '\0';
  totallen = snprintf(NULL, 0, "%sNH:i:%u%s", *line, noofmatches, pch+i);
  res = ALLOCMEMORY(space, NULL, char, totallen + 1);
  sprintf(res, "%sNH:i:%u%s", *line, noofmatches, pch+i);
  res[totallen] = '\0';
  FREEMEMORY(space, *line);
  *len = totallen;
  *line = res;
}


/*------------------------- bl_mergeBisulfiteBins -----------------------------
 *    
 * @brief  merging of bisulfite bins according to given read order between
 *         different bisulfite matching runs for each bin separately
 *         NOTE: this may be threaded if necessary later
 * @author Christian Otto
 *   
 */
void
se_mergeBisulfiteBins (void *space, bl_fileBinDomains_t *bsdomains, fasta_t **reads,
		       FILE *dev, bl_fileBinDomains_t *chrdomains, unsigned char remove,
                       Uint bestonly){
  Uint i, j, k, l, curlen, noofbins,
    noofbest, allocbest = 1000;    
  int cmp;
  char *curkey;
  FILE *fp;
  bl_mergefiles_t files;
  bl_mergefilematch_t **best;

  assert(bsdomains->noofdomains > 0);
  noofbins = bsdomains->domain[0].bins.noofbins;
  for (i = 1; i < bsdomains->noofdomains; i++){
    if (bsdomains->domain[i].bins.noofbins != noofbins){
      DBG("Inconsistent noofbins in domains. Exit forced.\n", NULL);
      exit(-1);
    }
  }

  best = ALLOCMEMORY(space, NULL, bl_mergefilematch_t **, allocbest);

  for (i = 0; i < noofbins; i++){
    NFO("Merging bisulfite bin %d.\n", i);
    /* init and open files */
    bl_mergefilesInit(space, &files, bsdomains->noofdomains);
    for (k = 0; k < files.nooffiles; k++){
      bl_mergefileInit(space, &files.f[k], bl_fileBinsOpen(space, &bsdomains->domain[k].bins.b[i], "r"));
    }

    /* perform merging */
    for (j = 0; j < reads[i]->noofseqs; j++){
      noofbest = 0;

      /* get next key */
      curkey = bl_fastaGetDescription(reads[i], j);
      curlen = bl_fastaGetDescriptionLength(reads[i], j);
      //DBG("queryfile:id=%d\tkey=%s\n", j, curkey);
      /* 
       * find match(es) with current key,
       * best pairing state, minimal edist (qry edist + mate edist)
       */
      for (k = 0; k < files.nooffiles; k++){
	while (1){
	  /* read next entry */
	  if (!files.f[k].complete){
	    bl_mergeReadNext(space, &files.f[k]);
	  }
          //DBG("files.f[%d]: curkey=%s\nmatch=%s\n", k, files.f[k].entry->key, files.f[k].entry->match);
	  /* 
           * end of file reached or next match with different key
           * Note: allow for partial matches (match one to the end)
           * due to clipping of /1 or /2 in paired-end data
           */
	  if (files.f[k].eof || ! bl_mergefileFastaIDCompare(curkey, curlen, files.f[k].entry->key,files.f[k].entry->keylen)){
	    break;
	  }
	  
	  /* compare current with previous best match (pairing state & edist) */
	  if (noofbest > 0){
	    cmp = bl_mergefilematchCompareFlag(best[0], files.f[k].entry);

            /* compare edit distance only in case of best-only */
            if (cmp == 0 && bestonly){
              cmp = bl_mergefilematchCompareEdist(best[0], files.f[k].entry);
            }
	  }
	  else {
	    cmp = 0;
	  }
          //DBG("cmp=%d\n", cmp);
	  
	  if (cmp >= 0){
	    if (cmp > 0){
	      /* new best match found => destruct previous ones */
	      for (l = 0; l < noofbest; l++){
		bl_mergefilematchDestruct(space, best[l]);
		FREEMEMORY(space, best[l]);
	      }
	      noofbest = 0;
	    }
	    /* extend best buffer */
	    if (noofbest == allocbest - 1){
	      allocbest *= 2;
	      best = ALLOCMEMORY(space, best, bl_mergefilematch_t **, allocbest);
	    }
	    /* append current match to best */
	    best[noofbest++] = files.f[k].entry;
	    
	    files.f[k].entry = ALLOCMEMORY(space, NULL, bl_mergefilematch_t, 1);
	    bl_mergefilematchInit(space, files.f[k].entry);
	  }
	  /* better match already found => clear data */
	  else {
	    bl_mergefilematchDestruct(space, files.f[k].entry);
	  }
	  files.f[k].complete = 0;
	}
      }

      for (k = 0; k < noofbest; k++){
	/* updating the 'number of matches' tag in SAM format (if necessary) */
	if (best[k]->noofmatches != noofbest){
	  bl_mergeUpdateTag(space, &best[k]->match, &best[k]->matchlen, noofbest);
	}
        if (best[k]->matematch != NULL && best[k]->noofmatematches != noofbest){
	  bl_mergeUpdateTag(space, &best[k]->matematch, &best[k]->matematchlen, noofbest);          
        }
	
	/* select output device */
	if (chrdomains == NULL) {
	  fp = dev;
	}
	else {
	  fp = bl_fileBinsOpen(space, bl_fileBinsDomainGetBin(chrdomains, best[k]->rname,
							      best[k]->rstart), "w");
	}
	
	/* output found match */
	fprintf(fp, "%s", best[k]->match);
	//fprintf(stderr, "best: %s", best[k]->match);

	if (best[k]->matematch != NULL){
          if (chrdomains != NULL){
            fp = bl_fileBinsOpen(space, bl_fileBinsDomainGetBin(chrdomains, best[k]->matername,
                                                                best[k]->materstart), "w");
          }
	  fprintf(fp, "%s", best[k]->matematch);
	  //fprintf(stderr, "%s", best[k]->matematch);    
	}

	/* clear data */
	bl_mergefilematchDestruct(space, best[k]);
	FREEMEMORY(space, best[k]);
      }
    }

    for (k = 0; k < files.nooffiles; k++){
      /* check whether match file is entirely processed */
      bl_mergeReadNext(space, &files.f[k]);
      if (!files.f[k].eof || files.f[k].complete){
	DBG("Files not yet entirely processed. Exit forced.\n", NULL);
        //DBG("files.f[%d]: key=%s\nmatch=%s\n", k, files.f[k].entry->match, files.f[k].entry->key);
	exit(-1);
      }
      /* destruct */
      bl_mergefileDestruct(space, &files.f[k]);
      /* close file */
      bl_fileBinsClose(&bsdomains->domain[k].bins.b[i]);
      if (remove){
	bl_rm(space, bsdomains->domain[k].bins.b[i].fname);
	bsdomains->domain[k].bins.b[i].unlinked = 1;
      }
    }
    bl_mergefilesDestruct(space, &files);
  }
  FREEMEMORY(space, best);  
}
