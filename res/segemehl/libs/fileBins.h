#ifndef FILEBINS_H
#define FILEBINS_H

/*
 *  fileBins.h
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
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <assert.h>
#include <sys/types.h>
#include "basic-types.h"
#include <pthread.h>


#ifndef HAVE_FSEEKO
    int fseeko(FILE *stream, off_t offset, int whence);
#endif

#ifndef HAVE_FTELLO
    off_t ftello(FILE *stream);
#endif

typedef struct bl_fileBinClass_s{
  char *classname;
  int classno;
  LLint start;
  LLint end;
} bl_fileBinClass_t;

typedef struct fileBin_s{
  FILE *fp;
  bl_fileBinClass_t *id;
  char *fname;
  unsigned char unlinked;
  off_t maxsize;
  pthread_mutex_t* mtx;
  unsigned char sorted;
  unsigned long long int lines;
} bl_fileBin_t;

typedef struct fileBins_s {
  Uint noofbins;
  bl_fileBin_t *b;
} bl_fileBins_t;

typedef struct fileBinDomain_s {
  char *domainname;
  Uint domainsize;
  bl_fileBins_t bins;
} bl_fileBinDomain_t;

typedef struct fileBinDomains_t {
  Uint noofdomains;
  Uint exp; //to the base of two
  bl_fileBinDomain_t *domain;
} bl_fileBinDomains_t;

typedef struct fileSort_s {
  LLint key;
  Uint len;
  char *line;
  off_t ptr;
} bl_fileSort_t;


void
bl_fileBinsUnlock (bl_fileBin_t *bin);

void
bl_fileBinsLock (bl_fileBin_t *bin);

void
bl_fileBinsCloseAll (bl_fileBins_t *bins);

void 
bl_fileBinsAdd(void *space, bl_fileBins_t* bins, Uint add,
    bl_fileBinClass_t* (*assigner)(void *, int, void *), void *nfo, char** names,
    char *template, Uint tmplen);

void
bl_fileBinsSortLine(void *space, bl_fileBins_t* bins, 
    unsigned char fileSort, char *filename, unsigned char ulink,
    LLint (*key)(char *, void*), void* nfo);

int
bl_fileBinsCClassSelect (void *id, void *nfo);

bl_fileBinClass_t*
bl_fileBinCClassAssign (void *space, int id, void *nfo);

int
bl_fileBinsClose(bl_fileBin_t *fx);

void
bl_fileBinsDestruct (void *space, bl_fileBins_t *bins);

void
bl_fileBinsInit(void *space, bl_fileBins_t *bins);

bl_fileBin_t *
bl_fileBinsFind (void *space, bl_fileBins_t* bins, 
    int (*selector)(void *id, void *nfo), void *nfo);

FILE*
bl_fileBinsOpen(void *space, bl_fileBin_t *bin, const char *mode);

void
bl_fileBinsWriteLn(void *space, bl_fileBin_t *fx, char *line);

unsigned char
bl_fileBinsIsOpen (bl_fileBin_t *fx);

char *
bl_fileBinsGetTemp(char *tmp, Uint tmplen);

void
bl_fileBinsGetInfo(bl_fileBins_t *);

void
bl_fileBinsCClassRename (void *space, bl_fileBins_t *fb, 
    char *bname, Uint bnamelen, char *suf, Uint suflen);

void
bl_fileBinsUnixSort (void *space, bl_fileBins_t *fb, const char *fieldstring, const char delim);

bl_fileBinDomains_t* 
bl_fileBinsDomainsInit(void *space, char **domainnames, Uint *domainsizes, 
    Uint noofdomains, Uint total, Uint avgbins, Uint maxbins, char *filenametemplate, Uint tmplen);

bl_fileBin_t *
bl_fileBinsDomainGetBin (bl_fileBinDomains_t *dms, char *domainname, Uint pos);

void
bl_fileBinDomainsDestruct (void *space, bl_fileBinDomains_t *dms);


void
bl_fileBinsDomainGetInfo (bl_fileBinDomains_t *dms);

void
bl_fileBinDomainsCloseAll (bl_fileBinDomains_t *dms);

void 
bl_fileBinDomainsUnixSort (void *space, bl_fileBinDomains_t *dms, const char *fldstr, const char delim); 


void
bl_fileBinDomainsMerge (void *space, bl_fileBinDomains_t *dms, 
    char *bname, Uint bnamelen, 
    char *suf, Uint suflen, char **header,
    unsigned char remove);

void
bl_fileBinDomainsSortMerge(void *space, bl_fileBinDomains_t *dms,
    char *bname, Uint bnamelen,
    char *suf, Uint suflen,
    const char *fldstr, const char delim,
    unsigned char remove);

Uint 
bl_fileBinsDomainsGetList(void *space, bl_fileBinDomains_t *domains, 
    char **domainnames[], Uint **domainsizes);

void
bl_fileBinDomainsDestruct (void *space, bl_fileBinDomains_t *dms);


#endif
