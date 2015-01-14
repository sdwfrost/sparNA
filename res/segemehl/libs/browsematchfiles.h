#ifndef BROWSE_MAT_FILE_H
#define BROWSE_MAT_FILE_H

/*
 *
 *	browsematchfiles.h
 *  a small browser
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 06.10.2010 01:48:30 CEST  
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "manout.h"
#include "matchfiles.h"
#include "evalmatchfiles.h"
#include <ncurses.h>

#define MAXPADLINES 100


#define WHITEONBLACK 1
#define BLACKONWHITE 2
#define REDONBLACK 3
#define BLUEONBLACK 4
#define WHITEONBLUE 5
#define BLUEONGREEN 6
#define BLACKONBLUE 7
#define BLUEONYELLOW 9
#define REDONYELLOW 10
#define BLUEONMAGENTA 11
#define WHITEONGREEN 12
#define REDONGREEN 13
#define WHITEONYELLOW 14
#define REDONWHITE 15
#define BLUEONWHITE 16
#define GREENONWHITE 17
#define MAGENTAONWHITE 18
#define CYANONWHITE 19
#define GREENONYELLOW 20
#define REDONBLUE 21


typedef struct {
  matchfile_t *file;
  fasta_t *set;
  matchfileFrame_t *curframe;
  matchfileFrameStats_t *curframestats;
  matchfileSampleStats_t *stats;
  matchfileindex_t *idx;
 
  annotationtrack_t *annotation; 
  Uint annotationoffset;

  WINDOW* annotationpad;
  Uint annotationscroll;
  WINDOW* pad;
  Uint scroll;
  Uint *map;
  Uint *imap;
  Uint offset;
  Uint width;

} matchfileView_t;

typedef struct {

  Uint noofviews;
  Uint activeview;
  matchfileView_t **views;
  WINDOW *activebox;
  int *pminrow;
  int *pmincol;
  int *sminrow;
  int *smincol;
  int *smaxrow;
  int *smaxcol;

} matchfilePanel_t;

void
bl_matchfileViewer(void *space, matchfile_t **files, Uint nooffiles,
    fasta_t *set,  annotationtrack_t *annotation, Uint start, Uint width ); 

void
bl_matchfileDestructView(void *space, matchfileView_t *view);

#endif
