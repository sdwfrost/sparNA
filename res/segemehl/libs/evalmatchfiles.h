#ifndef EVALMATCHFILE_H
#define EVALMATCHFILE_H

/*
 *
 *	evalmatchfiles.h
 *  evalutate matchfiles
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 10/14/2010 12:07:57 AM CEST  
 *
 */

/*
 * SVN
 * Revision of last commit: $Rev: 411 $
 * Author: $Author: steve $
 * Date: $Date: 2014-06-18 08:38:48 -0400 (Wed, 18 Jun 2014) $
 * Id: $Id: evalmatchfiles.h 411 2014-06-18 12:38:48Z steve $
 * Url: $URL: http://www2.bioinf.uni-leipzig.de/svn5/segemehl/libs/evalmatchfiles.h $
 */



#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "basic-types.h"
#include "matchfiles.h"
#include "biofiles.h"


#define MAX_STARTSITES 20000
#define EPSILON_NV 0.0000001
#define SDFRACTION 3
#define MINSUBPROB -4.0
#define MINRVPROB -3.0
#define MINMMPROB -3.0
#define MINEEPROB -7.0
#define STRANDPENALTY -3.0
#define QUALFACTOR 1.0
#define PLUS_STRAND 1
#define MINUS_STRAND (1 << 1)
#define BOTH_STRANDS (PLUS_STRAND | MINUS_STRAND)
#define EEWEIGHT .0
#define SPWEIGHT 0.75
#define RQWEIGHT 1.0
#define RRWEIGHT -1.0
#define RTWEIGHT 1500.0
#define MMWEIGHT 1.5 


typedef struct {
  Uint noofgroups;
  char *type;
  double *scr_cons;
  double *scr_ref;
  double *p_consx;
  double *p_cons;
  double *p_refx;
  double *p_ref;
  Uint **cnt;
} matchfileTestGroups_t;


typedef struct {
  char *ref;
  char *phreds;
  char len;
  char **alleles;
  Uint pos;
 
} indelvcf_t;

extern char * ntcode;
extern char * ntdecode;

  void
bl_matchfileEvalCrossSections (void *space,  matchfile_t **file, int *gropus, Uint nooffiles, fasta_t *set, Uint maxframesize, 
    Uint (*f)(void *, Uint fidx, Uint cidx, Uint pos, matchfileCross_t*, char, matchfileindex_t *, unsigned char, void *), void *nfo);

matchfileFrameStats_t *bl_matchfileFrameStats (void *space, matchfileFrame_t *frame);
void bl_matchfileDestructFrame(void *space, matchfileFrame_t *frame);
void bl_matchfileGetConsensus(matchfileFrame_t *frame);
Uint* bl_matchfileGetNTCounts(matchfileCross_t *cs);
double* bl_matchfileGetNTError(matchfileFrame_t *frame, Uint pos);
double* bl_matchfileGetNTRedundancy(matchfileFrame_t *frame, Uint pos);
double* bl_matchfileGetNTEdist(matchfileFrame_t *frame, Uint pos);
double* bl_matchfileGetNTReadPos(matchfileFrame_t *frame, Uint pos);
double* bl_matchfileGetNTReadPosVar(matchfileCross_t *);
matchfileFrameStats_t* bl_matchfileFrameStats(void *space, matchfileFrame_t *frame);
void bl_matchfileDestructFrameStats(void *space, matchfileFrameStats_t *stats);
void bl_matchfileRSSGNUPLOT(void *space, matchfileFrame_t *frame, matchfileFrameStats_t *stats);
void bl_matchfileCOVGNUPLOT(void *space, matchfileFrame_t *frame);
extern FILE *popen( const char *command, const char *modes);
extern int pclose(FILE *stream);
Uint bl_matchfileSampleCmp (Uint elemA, Uint elemB, void *toSort, void *info);
void bl_matchfileGetErrorDensity(void *space, matchfileFrame_t *frame, Uint pos, matchfileFrameStats_t *, void *nfo);
Uint bl_matchfileTest(void *space, Uint fidx, Uint cidx, Uint pos, matchfileCross_t *cs, char ref, matchfileindex_t *, unsigned char, void *nfo);
Uint bl_matchfileSampleCrossSections(void *space, matchfile_t *file, fasta_t *set, Uint n, 
    void (*f)(void *, matchfileFrame_t*, Uint, matchfileFrameStats_t *, void *), unsigned char **maps, Uint *mapsize, void *info);
void bl_matchfileGetConditionals (void *space, matchfileFrame_t *frame,
    Uint pos, matchfileFrameStats_t *stats, void *nfo);
matchfileSampleStats_t*
bl_matchfileInitSampleStats (void *space, Uint maxsample, Uint maxcover, Uint mincover, double minfrac, char entropyfilter, 
    Uint areasize, double maxareae);
void bl_matchfileDestructSampleStats (void *space, matchfileSampleStats_t *stats);
void bl_matchfileDumpSampleStats (matchfileSampleStats_t *stats);
void bl_matchfileFitGEV (void *space, matchfileSampleStats_t *stats);
unsigned char** bl_matchfileSmallMap (void *space, matchfile_t* file, Uint **mapsize);
Uint bl_matchfileSimpleGEV (void *space, Uint fidx, Uint cidx, Uint pos, matchfileCross_t *cs, 
    char ref, matchfileindex_t *idx, unsigned char show, void *nfo);
Uint bl_matchfileSimpleGATK ( void *space, Uint fidx, Uint cidx, Uint pos, matchfileCross_t *cs, 
    char ref, matchfileindex_t *idx, unsigned char show, void *nfo);
Uint bl_matchfileConsensus ( void *space, Uint fidx, Uint cidx, Uint pos, matchfileCross_t *cs, 
    char ref, matchfileindex_t *idx, unsigned char show, void *nfo);
Uint bl_matchfileGetCov(matchfileCross_t *cs, unsigned char allowdel);
char *bl_matchfileGetContext(fasta_t *fasta, Uint idx, Uint pos, Uint strand, Uint len);
void bl_matchfileGetScoreSample (void *space, matchfileFrame_t *frame, Uint pos, matchfileFrameStats_t *framestats, void *nfo);

void
bl_matchfileCrossStatsDestruct (matchfileCrossStats_t *css);
void
bl_matchfileCrossStatsInit (matchfileCrossStats_t *css, Uint len);

#endif
