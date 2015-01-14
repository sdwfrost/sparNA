#ifndef _PLOTMATCHFILES_H
#define _PLOTMATCHFILES_H

/*
 *
 *	plotmatchfiles.h
 *  gnuplot routines
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 10/22/2010 05:23:11 PM CEST  
 *
 */

void bl_matchfileQERRGNUPLOT(void *space, matchfileindex_t *index);
void bl_matchfilePERRGNUPLOT(void *space, matchfileindex_t *index);
void bl_matchfileCOVGNUPLOT(void *space, matchfileFrame_t *frame);
void bl_matchfileRSSGNUPLOT(void *space, matchfileFrame_t *frame, matchfileFrameStats_t *stats);
void bl_matchfileSUBGNUPLOT(void *space, matchfileindex_t *index);

#endif
