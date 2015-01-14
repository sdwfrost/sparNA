/*
 *  fileBins.c
 *  segemehl
 *
 *  Created by Steve Hoffmann on 09.02.10.
 *  Copyright 2010 University Leipzig. 
 *  All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
//#include <malloc/malloc.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <pthread.h>

#include "memory.h"
#include "fileio.h"
#include "basic-types.h"
#include "radixsort.h"
#include "fileBins.h"
#include "info.h"
#include "debug.h"


#define _FILE_OFFSET_BITS 64


/*---------------------------- bl_fileBinsGetInfo ----------------------------
 *    
 * @brief get info on file bins
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_fileBinsGetInfo (bl_fileBins_t *fb)
{   
  Uint i;

  if(!fb) {
      DBG("fileBins not initialized:\n", NULL);
      return;
    }

    NFO("total number of filebins: %d\n", fb->noofbins);
    for(i=0; i < fb->noofbins; i++) {
      NFO("Bin[%d] %s (classname: %s, range:%lld-%lld)\n", i, fb->b[i].fname, fb->b[i].id->classname, fb->b[i].id->start, fb->b[i].id->end);
    }
	return ;
}


/*------------------------- bl_fileBinsDomainGetInfo -------------------------
 *    
 * @brief get file bins domain info
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_fileBinsDomainGetInfo (bl_fileBinDomains_t *dms) {
  Uint i;

  NFO("total number of domains: %d\n", dms->noofdomains);
  for(i=0; i < dms->noofdomains; i++) {
    NFO("Domain[%d] %s, domainsize: %d\n",
        i, dms->domain[i].domainname, dms->domain[i].domainsize);
      bl_fileBinsGetInfo(&dms->domain[i].bins);
  }

  return;
}



/*--------------------------- bl_fileBinsWriteLine ---------------------------
 *    
 * @brief write a line
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_fileBinsWriteLn (void *space, bl_fileBin_t *fx, char *line)
{
  fprintf(fx->fp,"%s\n", line);
  fx->lines++;
  return ;
}


/*---------------------------- bl_fileBinsIsOpen -----------------------------
 *    
 * @brief returns 1 if file is open, 0 if not
 * @author Steve Hoffmann 
 *   
 */
 
unsigned char
bl_fileBinsIsOpen (bl_fileBin_t *fx)
{
  	return (fx->fp != NULL);
}

/*----------------------------- bl_fileBinsClose ------------------------------
 *    
 * @brief close a bin
 * @author Steve Hoffmann 
 *   
 */
 

int
bl_fileBinsClose(bl_fileBin_t *fx){
  int ret;

  assert(fx->fp);
  ret = fclose(fx->fp);
  fx->fp = NULL; 

  return ret;
}



/*------------------------------ fileBinsUnlock ------------------------------
 *    
 * @brief unlock file bin
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_fileBinsUnlock (bl_fileBin_t *bin)
{

  int ret;
  assert(bin);
  ret = pthread_mutex_trylock(bin->mtx);
  assert(ret == EBUSY);
  pthread_mutex_unlock(bin->mtx);
	
  return ;
}


/*------------------------------ bl_fileBinsLock -----------------------------
 *    
 * @brief file bin set lock
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_fileBinsLock (bl_fileBin_t *bin)
{
  assert(bin);
  pthread_mutex_lock(bin->mtx);
  return ;
}

/*----------------------------- bl_fileBinsOpen ------------------------------
 *    
 * @brief open a bin
 * @author Steve Hoffmann 
 *   
 */
 

FILE*
bl_fileBinsOpen(void *space, bl_fileBin_t* bin, const char *mode){

 /*empty files?*/ 
  if (!bin->fp) bin->fp = fopen(bin->fname, mode);

  if(!bin->fp) {
    DBG("filebins couldnt open file %s in mode '%s'. Exit forced.\n", bin->fname, mode);
   printf( "Error opening file: %s\n", strerror( errno ) ); 
    exit(-1);
  }

  return bin->fp;
}

/*----------------------------- bl_fileBinsFind -----------------------------
 *    
 * @brief return bin 
 * @author Steve Hoffmann 
 *   
 */
 

bl_fileBin_t *
bl_fileBinsFind (void *space, bl_fileBins_t* bins, 
    int (*selector)(void *id, void *nfo), void *nfo) 
{
	
  Uint i;

  for(i=0; i < bins->noofbins; i++) {
    if(selector(bins->b[i].id, nfo)) {
      return &bins->b[i];
    }
  }
  return NULL;
}


/*----------------------------- bl_fileBinsInit ------------------------------
 *    
 * @brief initialize file bins
 * @author Steve Hoffmann 
 *   
 */

void
bl_fileBinsInit(void *space, bl_fileBins_t *bins) {

  bins->b = NULL;
  bins->noofbins = 0;

  return;
}



/*----------------- bl_fileBinsDomainsGetNames(space, desc) ------------------
 *    
 * @brief pass list of k domain names and their length 
 * @author Steve Hoffmann 
 *   
 */
 
Uint 
bl_fileBinsDomainsGetList(void *space, bl_fileBinDomains_t *domains, 
    char **domainnames[], Uint **domainsizes) 
{
  Uint i=0, k=0;
  char **list=NULL;
  Uint *ll = NULL;

  for(i=0; i < domains->noofdomains; i++) {
    if(list == NULL || 
        list[k-1] != domains->domain[i].domainname) {
      list = ALLOCMEMORY(space, list, char*, k+1);
      ll   = ALLOCMEMORY(space, ll, Uint, k+1);
      list[k] = domains->domain[i].domainname;
      ll[k] = domains->domain[i].domainsize;
      k++;
    }
  }
	
  *domainnames = list;
  *domainsizes = ll;
  return k;
}

/*--------------------------- bl_fileBinsCloseAll ----------------------------
 *    
 * @brief closes all file bins
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_fileBinsCloseAll (bl_fileBins_t *bins)
{
  Uint i;
  assert(bins); 
  for(i=0; i < bins->noofbins; i++) {
    if(bl_fileBinsIsOpen(&bins->b[i])) bl_fileBinsClose(&bins->b[i]);
  }
  return ;
}

/*------------------------------ bl_fileBinsAdd ------------------------------
 *    
 * @brief add file bins to file bin container
 * @author Steve Hoffmann 
 *   
 */
 

void
bl_fileBinsAdd (void *space, bl_fileBins_t *bins, Uint add, 
    bl_fileBinClass_t* (*assigner)(void *, int, void *), void *nfo, char **filenames,
    char *template, Uint tmplen) {
  Uint i;
  char *fname;

  bins->b = ALLOCMEMORY(bins->b, NULL, bl_fileBin_t, bins->noofbins+add);
  bins->noofbins += add;

  for(i=0; i < bins->noofbins; i++) {
    if(filenames == NULL) {
      fname = bl_getTempFile(template, tmplen);
    } else {
      fname = filenames[i];
    }
    
    bins->b[i].unlinked = 0;
    bins->b[i].fname = fname;
    
    if (assigner) 
      bins->b[i].id = assigner(space, i, nfo);
    else 
      bins->b[i].id = NULL;
    
    bins->b[i].lines=0;
    bins->b[i].sorted = 0;
    bins->b[i].fp = NULL;
    bins->b[i].mtx = NULL;
    bins->b[i].mtx = ALLOCMEMORY(space, NULL, pthread_mutex_t, 1);
    pthread_mutex_init(bins->b[i].mtx, NULL);
  }

  return;
}


/*--------------------------- bl_fileBinsDestruct ----------------------------
 *    
 * @brief destruct file bins
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_fileBinsDestruct (void *space, bl_fileBins_t *bins)
{
  Uint i;

  for(i=0; i < bins->noofbins; i++) {
    FREEMEMORY(space, bins->b[i].fname);
    FREEMEMORY(space, bins->b[i].id);
    FREEMEMORY(space, bins->b[i].mtx);
  }

  FREEMEMORY(space, bins->b);
  bins->noofbins = 0;
  bins->b = NULL;
  return;
}


/*--------------------------- bl_fileBinDomainsDestruct ----------------------------
 *    
 * @brief destruct file bin domains
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_fileBinDomainsDestruct (void *space, bl_fileBinDomains_t *dms)
{
  Uint i;

  for(i=0; i < dms->noofdomains; i++) {
    bl_fileBinsDestruct(space, &dms->domain[i].bins);
    FREEMEMORY(space, dms->domain[i].domainname);
  }

  FREEMEMORY(space, dms->domain);

  dms->noofdomains = 0;
  dms->domain = NULL;
  return;
}

/*----------------------.-- bl_fileBinDomainsCloseAll ------------------------
 *    
 * @brief closes all file bins
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_fileBinDomainsCloseAll (bl_fileBinDomains_t *dms)
{
  Uint i;
  assert(dms->domain);

  for(i=0; i < dms->noofdomains; i++) {
    bl_fileBinsCloseAll(&dms->domain[i].bins);
  }
  return ;
}


/*-------------------------- bl_fileBinDomainsInit ---------------------------
 *    
 * @brief initalize file bin domains 
 *        find next highest power of two; make binsize a power of two
 * @author Steve Hoffmann 
 *   
 */

bl_fileBinDomains_t* 
bl_fileBinsDomainsInit(void *space, char **domainnames, Uint *domainsizes, 
    Uint noofdomains, Uint totalsize, Uint avgbins, Uint maxbins, 
    char *filenametemplate, Uint tmplen){

  Uint i, j, noofbins, e=0, binsize, maxbinperdomain, 
  maxdomainsize=0, est=0;
  bl_fileBinDomains_t* dms;
  bl_fileBinClass_t *ptr;


  if(noofdomains > maxbins || maxbins == 0 || avgbins > maxbins) {
    DBG("bl_fileBinDomainsInit: maxbins=%u < %u=noofdomains\n", maxbins, noofdomains);
    return NULL;
  }

  for (i=0; i < noofdomains; i++) maxdomainsize = domainsizes[i];

  dms = ALLOCMEMORY(space, NULL, bl_fileBinDomains_t, 1);
  dms->noofdomains = noofdomains;
  binsize = ceil(totalsize/avgbins);

  while (((binsize-1) >> ++e) >= 1);
  
  if (e >= 31) {
    DBG("bl_fileBinDomainsInit: binsize 2^%u is out of range.\n", e);
    return NULL;
  }
  
  binsize = (1 << e);

  for(i=0; i < noofdomains; i++) {
    est += ceil((double)domainsizes[i]/binsize);
  }

  if(est >= maxbins) {
    maxbinperdomain = floor(maxbins/noofdomains);
    binsize = maxdomainsize/maxbinperdomain + 1; 
    while (((binsize-1) >> ++e) >= 1);
    binsize = (1 << e);
     if (e >= 31) {
        DBG("bl_fileBinDomainsInit: binsize 2^%u is out of range.\n", e);
        return NULL;
    }
  }
  
  dms->exp = e;
  dms->domain = ALLOCMEMORY(space, NULL, bl_fileBinDomain_t, noofdomains);

  for (i=0; i < noofdomains; i++) {
    dms->domain[i].domainsize = domainsizes[i];
    dms->domain[i].domainname = ALLOCMEMORY(space, NULL, char, strlen(domainnames[i])+1);
    memmove(dms->domain[i].domainname, domainnames[i], strlen(domainnames[i]));
    dms->domain[i].domainname[strlen(domainnames[i])] = '\0';

    noofbins = (domainsizes[i] >> e) + ((domainsizes[i] & (binsize-1)) > 0); 

    bl_fileBinsInit(space, &dms->domain[i].bins);
    bl_fileBinsAdd (space, &dms->domain[i].bins, noofbins, NULL, NULL, NULL,
        filenametemplate, tmplen);

    dms->domain[i].bins.noofbins = noofbins;

    for (j=0; j < noofbins; j++) {
      ptr = ALLOCMEMORY(space, NULL, bl_fileBinClass_t, 1);      
      ptr->start = j*binsize;
      ptr->end = (j+1)*binsize -1;
      ptr->classname = NULL;
      dms->domain[i].bins.b[j].id = ptr;
    }
  }

  return dms;
}


/*------------------------- bl_fileBinDomainsFindBin -------------------------
 *    
 * @brief find domain and return bin
 * @author Steve Hoffmann 
 *   
 */
 
bl_fileBin_t *
bl_fileBinsDomainGetBin (bl_fileBinDomains_t *dms, char *domainname, Uint pos)
{
  Uint i;

  for(i=0; i < dms->noofdomains; i++) {
    //fprintf(stderr,"compare %s == %s\n", dms->domain[i].domainname, domainname);
    if(strcmp(dms->domain[i].domainname, domainname) == 0) break;
  }

  if (i == dms->noofdomains) 
    return NULL;

  return &dms->domain[i].bins.b[(pos >> dms->exp)]; 
}



/*--------------------------- bl_fileBinsUnixSort ----------------------------
 *    
 * @brief start unix sort tool
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_fileBinsUnixSort (void *space, bl_fileBins_t *fb, const char *fldstr, const char delim)
{
    Uint i;
    //int ret;
    char *filename;

    for(i=0; i < fb->noofbins; i++) {
      NFO("sorting file '%s'.\n", fb->b[i].fname);
      filename = fb->b[i].fname;
      //not used: ret = bl_UnixSort(space, filename, fldstr);
      bl_UnixSort(space, filename, fldstr, delim);
    }

    return ;
}




/*-------------------------- fileBinsDomainUnixSort --------------------------
 *    
 * @brief sort all bins of all domains
 * @author Steve Hoffmann 
 *   
 */
 
void 
bl_fileBinDomainsUnixSort (void *space, bl_fileBinDomains_t *dms, const char *fldstr, const char delim) 
{

    Uint i, j;
    //int ret;
    char *filename;

    for(i=0; i < dms->noofdomains; i++) {
      NFO("sorting domain %d.\n", i);
      for(j=0; j < dms->domain[i].bins.noofbins; j++) {
        filename = dms->domain[i].bins.b[j].fname;
        //not used: ret = bl_UnixSort(space, filename, fldstr);
        bl_UnixSort(space, filename, fldstr, delim);
      }
    }

	return ;
}

void
bl_fileBinDomainsSortMerge(void *space, bl_fileBinDomains_t *dms,
    char *bname, Uint bnamelen,
    char *suf, Uint suflen,
    const char *fldstr, const char delim,
    unsigned char remove) 
{

  char *cname, *ccname, *newname;
  char **filenames;
  Uint i,j,cnamelen;



  for(i=0; i < dms->noofdomains; i++) {

    filenames = ALLOCMEMORY(space, NULL, char*, dms->domain[i].bins.noofbins);
    cname = dms->domain[i].domainname; 
    cnamelen = strlen(cname);
    ccname = bl_replacenonalphanum(cname, cnamelen);
    newname = ALLOCMEMORY(space, NULL, char, cnamelen + bnamelen + suflen + 4);
    sprintf(newname, "%s_%s.%s", bname, ccname, suf);
    

    for(j=0; j < dms->domain[i].bins.noofbins; j++) {
      filenames[j] = dms->domain[i].bins.b[j].fname;
    }

    bl_UnixSortMerge(space, filenames, dms->domain[i].bins.noofbins, fldstr, delim, newname);

    FREEMEMORY(space, ccname);
    FREEMEMORY(space, newname);
  }

}


/*------------------------- bl_fileBinsDomainsMerge --------------------------
 *    
 * @brief merge all bins of all domains
 * @author Steve Hoffmann 
 *   
 */

  void
bl_fileBinDomainsMerge (void *space, bl_fileBinDomains_t *dms, 
    char *bname, Uint bnamelen, 
    char *suf, Uint suflen, char **header,
    unsigned char remove)
{

  FILE *outfile=NULL, *fp; 
  size_t buffersize = 1024, len; 
  off_t fs;
  char *buffer, *cname, *ccname, *newname;
  Uint i,j,cnamelen; 

  buffer = ALLOCMEMORY(space, NULL, char, buffersize);

  for(i=0; i < dms->noofdomains; i++) {

    cname = dms->domain[i].domainname; 
    cnamelen = strlen(cname);
    ccname = bl_replacenonalphanum(cname, cnamelen);
    newname = ALLOCMEMORY(space, NULL, char, cnamelen + bnamelen + suflen + 4);
    sprintf(newname, "%s_%s.%s", bname, ccname, suf);
    

    if(header && header[i]) {
      outfile = fopen(newname,"w");
      fprintf(outfile, "%s", header[i]);
      fclose(outfile);
    }
    
    outfile = fopen(newname, "ab");
    if (!outfile) {
      DBG("Opening of file %s failed. Exit forced.\n", newname);
      exit(EXIT_FAILURE);
    }

    for(j=0; j < dms->domain[i].bins.noofbins; j++) {

      fp = fopen(dms->domain[i].bins.b[j].fname, "rb");

      if (fp == NULL){
        DBG("Opening of file %s failed. Exit forced.\n", 
            dms->domain[i].bins.b[j].fname);
        exit(EXIT_FAILURE);
      }

      fseek (fp , 0 , SEEK_END);
      fs = ftello(fp);
      rewind (fp);

      while((len = fread(buffer, 1, buffersize, fp)) > 0) {
        fwrite(buffer, 1, len, outfile);
        fs -= len;
      }

      if(fs > 0) {
        DBG("Could not read %s entirely (fs:%zu)\n", 
            dms->domain[i].bins.b[j].fname, fs);
      }

      fclose(fp);
      if (remove) {
        bl_rm(space,  dms->domain[i].bins.b[j].fname);
        dms->domain[i].bins.b[j].unlinked = 1;
      }
    }

    FREEMEMORY(space, newname);
    FREEMEMORY(space, ccname);
  }

  FREEMEMORY(space, buffer);
  fclose(outfile);
  return ;
}





/* helpers 
 *
 *
 *
 *
 *
 */


/*--------------------------- bl_fileBinsSortLine ----------------------------
 *    
 * @brief sort lines in a bin (in memory or directly in the file)
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_fileBinsSortLine(void *space, bl_fileBins_t* bins, 
    unsigned char fileSort, char *filename, unsigned char remove,
    LLint (*keyaccess)(char *, void*), void* nfo) {

  Uint i, k, j=0, len;
  bl_fileSort_t *data;
  LLint res;
  char *line, *tmpname=NULL;
  off_t fs;
  int fseekres;
  size_t fwriteres, freadres;
  FILE *fp;
  FILE *outfile;
  struct stat st;
  char *template = "filebinsort";

  if (filename) {
    outfile = fopen(filename, "w"); 
  } else {
    tmpname = bl_getTempFile(template, 11);
    outfile = fopen(tmpname, "w");
  }

  if (outfile == NULL){
    fprintf(stderr, "Opening temp file failed. Exit forced.\n");
    exit(EXIT_FAILURE);
  }

  for(i=0; i < bins->noofbins; i++) {

    data = malloc(sizeof(bl_fileSort_t)*bins->b[i].lines); 
    if(!data) {
      bins->b[i].sorted = 0;
      fprintf(stderr,"warning: not enough memory for fileBins. Try unix sort.");
      continue;
    } else {
      fileSort=1;
    }

    if (stat(bins->b[i].fname, &st) == 0) {
      line = malloc(st.st_size);
      if(!line) { 
        bins->b[i].sorted = 0;
        fprintf(stderr,"warning: not enough memory for fileBins. Try sort.");
        continue;
      } 
      free(line);
    } else {
      continue;
    }

    fp = fopen(bins->b[i].fname, "r");
    if (fp == NULL){
      fprintf(stderr, "Opening file %s failed. Exit forced.\n", 
          bins->b[i].fname);
      exit(EXIT_FAILURE);
    }

    fs = ftello(fp); 
    if (fs == -1) {
      fprintf(stderr,"File access error for %s. Exit forced.\n", 
          bins->b[i].fname);
      exit(EXIT_FAILURE);
    }
    j = 0;

    while((len = bl_fgets(space, fp, &line)) != EOF) {
      res = keyaccess(line, nfo);

      data[j].key = res;
      data[j].ptr = fs;
      data[j].len = len;

      if (fileSort) {
        FREEMEMORY(space, line);
        data[j].line = NULL;
        fs = ftello(fp);

        if (fs == -1) {
          fprintf(stderr,"File access error for %s. Exit forced.\n", 
              bins->b[i].fname);
          exit(EXIT_FAILURE);
        }
      } else {
        data[j].line = line;
      }
      j += 1;
    }

    if(fileSort) FREEMEMORY(space, line); 
    bl_radixSortKeyFirst(space, data, sizeof(bl_fileSort_t), j, 16);

    for(k=0; k < j; k++) {
      if(fileSort) {
        line = ALLOCMEMORY(space, NULL, sizeof(char), data[k].len+1);
        fseekres = fseeko(fp, data[k].ptr, 0);
        if (fseekres == -1) {
          fprintf(stderr,"File access error for %s. Exit forced.\n", 
              bins->b[i].fname);
          exit(EXIT_FAILURE);
        }
        freadres = fread(line, sizeof(char), data[k].len+1, fp);
        if (freadres != data[k].len+1) {
          fprintf(stderr,"File access error for %s. Exit forced.\n", 
              bins->b[i].fname);
          exit(EXIT_FAILURE);
        }
        fwriteres = fwrite(line, sizeof(char), data[k].len+1, outfile);
        if (fwriteres != data[k].len+1) {
          fprintf(stderr,"File access error for %s. Exit forced.\n", tmpname);
          exit(EXIT_FAILURE);
        }
        FREEMEMORY(space, line);
      } else {
        fprintf(outfile, "%s\n", data[k].line);
        FREEMEMORY(space, data[k].line);
      }
    }

    fclose(fp);

    if (filename == NULL) {
      fclose(outfile);
      unlink(bins->b[i].fname);
      rename(tmpname, bins->b[i].fname);
      outfile = fopen(tmpname, "w");
      if (outfile == NULL){
        fprintf(stderr, "Opening temp file failed. Exit forced.\n");
        exit(EXIT_FAILURE);
      }
    } else if(remove) {
      unlink(bins->b[i].fname);
      bins->b[i].unlinked = 1;
    }

    FREEMEMORY(space, data); 
  }

  fclose(outfile);
  return;	
}


/*----------------------------- bl_fileBinsMerge -----------------------------
 *    
 * @brief merge file bins
 * @author Steve Hoffmann 
 *   
 */

void
bl_fileBinsMerge(void *space, char *filename, bl_fileBins_t* bins, 
    unsigned char delete) {

  FILE *outfile, *fp; 
  char *line;
  Uint i, len;

  outfile = fopen(filename, "w");
  for(i=0; i < bins->noofbins; i++) {    
    fp = fopen(bins->b[i].fname, "r");
    if (fp == NULL){
      fprintf(stderr, "Opening of file %s failed. Exit forced.\n", 
          bins->b[i].fname);
      exit(EXIT_FAILURE);
    }
    fprintf(stderr, "start file\n");

    while((len = bl_fgets(space, fp, &line)) != EOF) {
      fprintf(outfile, "%s\n", line);
    }
    //unlink file
  }
}


/*------------------------- bl_fileBinsCClassSelect --------------------------
 *    
 * @brief select classname
 * @author Steve Hoffmann 
 *   
 */
 
int
bl_fileBinsCClassSelect (void *id, void *nfo)
{
  bl_fileBinClass_t *elem;
  char *toSelect;

  elem = (bl_fileBinClass_t*) id;
  toSelect = (char*) nfo;
  
  if (strcmp(elem->classname, toSelect) == 0 || elem->classname[0]=='*') {
    return 1;
  }

  return 0;
}


/*------------------------- bl_fileBinsCClassRename --------------------------
 *    
 * @brief rename to classname
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_fileBinsCClassRename (void *space, bl_fileBins_t *fb, 
    char *bname, Uint bnamelen, char *suf, Uint suflen)
{
    char *newname, *cname, *ccname;
    Uint i, cnamelen;
    int ret;

    for(i=0; i < fb->noofbins; i++) {
      cname = fb->b[i].id->classname; 
      cnamelen = strlen(cname);
      ccname = bl_replacenonalphanum(cname, cnamelen);
      newname = ALLOCMEMORY(space, NULL, char, cnamelen + bnamelen + suflen + 4);
      sprintf(newname, "%s_%s.%s", bname, ccname, suf);

      ret = rename(fb->b[i].fname, newname);
      assert(ret != -1);

      FREEMEMORY(space, ccname);
      FREEMEMORY(space, newname);
    }

	return ;
}



/*-------------------------- bl_fileBinCClassAssign --------------------------
 *    
 * @brief assign classname
 * @author Steve Hoffmann 
 *   
 */
 
bl_fileBinClass_t*
bl_fileBinCClassAssign (void *space, int id, void *nfo)
{
  bl_fileBinClass_t *ptr;
  char **classnames;

  ptr = ALLOCMEMORY(space, NULL, bl_fileBinClass_t, 1);
  classnames = (char**) nfo;
  ptr->classname = classnames[id];

  return ptr;
}



