#ifndef MATCHFILESFIELDS_H
#define MATCHFILESFIELDS_H

/*
 *
 *	matchfilesfields.h
 *  
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 01.03.2011 17:40:56 CET  
 *
 */
#include "matchfiles.h"

Uint bl_matchfileGetPNext (stringset_t *fields, unsigned char fmt);
char* bl_matchfileGetRNext ( stringset_t *fields, unsigned char fmt);
char* bl_matchfileGetQname(stringset_t *fields, unsigned char fmt);
Uint bl_matchfileGetFlag(stringset_t *fields, unsigned char fmt);
unsigned char bl_matchfileIsHeader(char *buffer, Uint len, unsigned char fmt);
Uint bl_matchfileGetStartPos(stringset_t *fields, unsigned char fmt);
Uint bl_matchfileGetEndPos(stringset_t *fields, unsigned char fmt);
char* bl_matchfileGetRead(stringset_t *fields, unsigned char fmt);
char* bl_matchfileGetQual(stringset_t *fields, unsigned char fmt);
char bl_matchfileGetStrand(stringset_t *fields, unsigned fmt);
char* bl_matchfileGetChrom(stringset_t *fields, unsigned fmt);
Uint bl_matchfileGetMatchCnt(stringset_t *fields, unsigned char fmt);
Uint bl_matchfileGetEdist(stringset_t *fields, unsigned char fmt);
Uint bl_matchfileGetBisulfite(stringset_t *fields, unsigned char fmt);
char* bl_matchfileGetAln(stringset_t *fields, unsigned char fmt);
char* bl_matchfileGetDiffString(stringset_t *fields, unsigned char fmt);
Uint bl_matchfileGetPrevPos(stringset_t *fields, unsigned char fmt);
Uint bl_matchfileGetNextPos(stringset_t *fields, unsigned char fmt);
char* bl_matchfileGetNextChr(stringset_t *fields, unsigned char fmt);
char* bl_matchfileGetPrevChr(stringset_t *fields, unsigned char fmt);
Uint bl_matchfileGetSplitStart(stringset_t *fields, unsigned char fmt);
Uint bl_matchfileGetSplitEnd(stringset_t *fields, unsigned char fmt);
Uint bl_matchfileGetSplitNumber(stringset_t *fields, unsigned char fmt);
Uint bl_matchfileGetPrevFlag(stringset_t *fields, unsigned char fmt);
Uint bl_matchfileGetNextFlag(stringset_t *fields, unsigned char fmt);
Uint bl_matchfileGetMappingID(stringset_t *fields, unsigned char fmt);
matchfileRec_t * bl_matchfileGetMatchFileRec(matchfileRec_t *rec, Uint fields, stringset_t *token, Uint fmt);

#endif
