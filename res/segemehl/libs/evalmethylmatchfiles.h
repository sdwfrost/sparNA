#ifndef EVALMETHYLMATCHFILES_H
#define EVALMETHYLMATCHFILES_H

/**
 * evalmethylmatchfiles.c
 * evaluation and statistics of matchfiles from methylC-seq
 *
 * @author Christian Otto & Helene Kretzmer
 * @email christian@bioinf.uni-leipzig.de
 * @company Bioinformatics, University of Leipzig
 * @date Thu May  2 10:07:27 EDT 2013
 *  
 */

/*
 * SVN
 * Revision of last commit: $Rev: 408 $
 * Author: $Author: steve $
 * Date: $Date: 2014-06-12 07:10:00 -0400 (Thu, 12 Jun 2014) $
 * Id: $Id: evalmethylmatchfiles.h 408 2014-06-12 11:10:00Z steve $
 * Url: $URL: http://www2.bioinf.uni-leipzig.de/svn5/segemehl/libs/evalmethylmatchfiles.h $
 */

matchfileCross_t* bl_matchfileGetBSCross(matchfileCross_t *cs);
Uint bl_matchfileGetBSStrand(char ref);
char bl_matchfileGetBSBase(Uint strand, unsigned char conv);
Uint bl_matchfileGetBSCount(matchfileCross_t *bscs, Uint strand, unsigned char conv);
double bl_matchfileGetBSRateSimple(matchfileCross_t *bscs, Uint strand, unsigned char allowdel);
Uint bl_matchfileCallMethylSimple ( void *space, Uint fidx, Uint cidx, Uint pos, matchfileCross_t *cs, 
    char ref, matchfileindex_t *idx, unsigned char show, void *nfo);
Uint bl_matchfileCalcMethylBias ( void *space, Uint fidx, Uint cidx, Uint pos, matchfileCross_t *cs, 
    char ref, matchfileindex_t *stats, unsigned char show, void *nfo);
void bl_matchfileSampleStatsBS(void *space, matchfileFrame_t *frame, Uint pos, matchfileFrameStats_t *, void *nfo);
int bl_matchfileSampleCrossSectionsBS(void *space, matchfile_t *file, fasta_t *set, Uint n, 
    void (*f)(void *, matchfileFrame_t*, Uint, matchfileFrameStats_t *, void *), void *info);
int *bl_matchfileGetReadQualBS (void *space, matchfileFrame_t *frame, Uint pos);
int *bl_matchfileGetReadPosBS (void *space, matchfileFrame_t *frame, Uint pos);

#endif
