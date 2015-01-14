#ifndef SEGEMEHL_H
#define SEGEMEHL_H

/*
 *
 *	segemehl.h
 *  declarations for threaded segemehl
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 01/02/2008 10:12:46 AM CET  
 *
 *  SVN
 *  Revision of last commit: $Rev: 101 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-12-08 02:29:27 +0100 (Mon, 08 Dec 2008) $
 *
 *  Id: $Id: segemehl.h 101 2008-12-08 01:29:27Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/src/segemehl.h $
 *  
 */

#include <stdlib.h>
#include <stdio.h>
#include "biofiles.h"
#include "fileio.h"
#include "stringutils.h"
#include "charsequence.h"
#include "multicharseq.h"
#include "sufarray.h"
#include "fileBins.h"
#include <pthread.h>
#include "version.h"
#include "vtprogressbar.h"

#define MAXFILEBINS 50

typedef struct segemehl_s {
    void *space;
    char *outfile;
    char *queryfilename;
    char *matefilename;
    char *idxfilename;
    char *idx2filename;
    char *dbfilename;
    char *nomatchname;
    char *matelessfilename;
    char *softclip3Prime;
    char *softclip5Prime;
    char *filebinbasename;
    char *splitfilebasename;
    FILE *dev;
    FILE *nomatchdev;
    FILE *splitdev;
    bl_fileBinDomains_t* bins;
    bl_fileBinDomains_t *splitbins;
    int (*slct)(void *, void *);
    Suffixarray *sarray;
    MultiCharSeq *seq;
    fasta_t *reads;
    fasta_t *fasta;
    char * cmdline;
    Uint hardclip;
    Uint hardclip3Prime;
    Uint hardclip5Prime;
    Uint softclip3PrimeLen;
    Uint softclip5PrimeLen;
    Uint totallength;
    Uint minsize;
    Uint *mapmatches;
    Uint *counter;
    Uint M;
    Uint jump;
    Uint s_ext;
    Uint p_mis;
    Uint Xoff;
    Uint k;
    Uint k_p;
    Uint rep_type;
    Uint check;
    Uint kMis;
    Uint threadno;
    Uint threadid;
    Uint bestonly;
    Uint maxinsertsize;
    Uint polyAlen;
    Uint minclipscr3;
    Uint minclipscr5;
    Uint bisulfiterun;
    Uint bisulfiteprotocol;
    Uint bisulfitemerging;
    Uint bisulfite;
    Uint strand;
    Uint minfragmentalignlen;
    int  minfragmentalignscore;
    Uint minsplicedaligncover;
    int minsplicedalignscore;
    unsigned char split;
    unsigned char index;
    unsigned char bining;
    unsigned char matchingstat;
    unsigned char align;
    unsigned char mute;
    unsigned char oldversion;
    unsigned char showMateless;
    unsigned char polyA;
    unsigned char order;
    unsigned char sam;
    unsigned char autoclip;
    unsigned char SAMmeop;
    unsigned char SAMpairstat;
    unsigned char nohead;
    double maxevalue;
    double maxsplitevalue;
    int accuracy;
    int clipacc;
    Uint bedist;
    Uint fusion;
    double chainscorescale;
    pthread_mutex_t *mtx;
    pthread_mutex_t *mtx2;
    pthread_mutex_t *mtx3;
} segemehl_t;

typedef struct checkthread_s {
    Uint noofseqs;
    Uint *counter;
} checkthreadinfo_t;

inline static void
se_setdefault(segemehl_t *info) {
  info->bins = NULL;
  info->split = 0;
  info->splitbins = NULL;
  info->bining = 0;
  info->slct = bl_fileBinsCClassSelect;
  info->dev = stdout;
  info->idxfilename = NULL;
  info->idx2filename = NULL;
  info->dbfilename = NULL;
  info->queryfilename = NULL;
  info->matefilename=NULL;
  info->splitfilebasename = NULL;
  info->sarray = NULL;
  info->seq = NULL;
  info->fasta = NULL;
  info->polyAlen = 0;
  info->reads = NULL;
  info->outfile=NULL;
  info->totallength=0;
  info->counter=0;
  info->minsize = 12;
  info->k_p = 1;
  info->index = 1;
  info->threadno = 1;
  info->threadid = 0;
  info->accuracy = 90;
  info->clipacc = 70;
  info->maxsplitevalue = 50;
  info->M = 100;
  info->jump = 0;
  info->s_ext = 2;
  info->p_mis = 4;
  info->Xoff = 8;
  info->rep_type = 12;
  info->kMis = 0;
  info->mute = 0;
  info->maxinsertsize = 5000;
  info->matchingstat = 0;
  info->bestonly = 1;
  info->check = 0;
  info->maxevalue = 5;
  info->space = NULL;
  info->nomatchname = NULL;
  info->nomatchdev = NULL;
  info->splitdev = NULL;
  info->showMateless = 1;
  info->bedist=0;
  info->autoclip=0;
  info->align=0;
  info->order=0;
  info->sam = 1;
  info->bisulfiterun = 0;
  info->bisulfiteprotocol = 0;
  info->bisulfitemerging = 0;
  info->bisulfite = 0;
  info->strand = 0;
  info->oldversion=0;
  info->fusion=0;
  info->hardclip = 0;
  info->hardclip3Prime = 0;
  info->hardclip5Prime = 0;
  info->softclip3Prime = NULL;
  info->softclip5Prime = NULL;
  info->filebinbasename = NULL;
  info->polyA = 0;
  info->softclip3PrimeLen = 0;
  info->softclip5PrimeLen = 0;
  info->minclipscr3 = 5;
  info->minclipscr5 = 5;
  info->SAMmeop = 0;
  info->SAMpairstat = 1;
  info->chainscorescale = 1.0;
  info->minfragmentalignlen = 20;
  info->minfragmentalignscore = 18;
  info->minsplicedaligncover = 80;
  info->cmdline = NULL;
  info->minsplicedalignscore = 2*18;
  info->mtx = ALLOCMEMORY(space, NULL, pthread_mutex_t, 1);
  info->mtx2 = ALLOCMEMORY(space, NULL, pthread_mutex_t, 1);
  info->mtx3 = ALLOCMEMORY(space, NULL, pthread_mutex_t, 1);
  pthread_mutex_init(info->mtx, NULL);
  pthread_mutex_init(info->mtx2, NULL);
  pthread_mutex_init(info->mtx3, NULL);
}

extern const char *SORT[];
extern const char *SORTBIN[];
extern const char SORTDELIM;
void se_updateProgressBar(Uint k, segemehl_t *nfo);

#endif
