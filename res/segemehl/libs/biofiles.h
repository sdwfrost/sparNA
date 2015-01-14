#ifndef _BIOFILES_
#define _BIOFILES_

/*
 *
 *	biofiles.h
 *  declarations
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 07/10/2007 02:32:29 PM CEST  
 *
 *  SVN
 *  Revision of last commit: $Rev: 19 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-05-14 15:43:29 +0200 (Wed, 14 May 2008) $
 *
 *  Id: $Id: biofiles.h 19 2008-05-14 13:43:29Z steve $
 *  Url: $URL: file:///homes/bierdepot/steve/svn/segemehl/trunk/libs/biofiles.h $
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "stringutils.h"
#include "basic-types.h"
#include "charsequence.h"
#include "zran.h"

#define ID 1
#define IDEND 2
#define SEQEND 3
#define QUALIDEND 4
#define END 5

#define GFFITEM 0
#define BEDITEM 1
#define SNPITEM 2

typedef struct seqaccesspoint_s {
  off_t offset;
  Uint fileid;
  Uint noofseqs;
  Uint cumnoofseqs;
} seqaccesspoint_t;


typedef struct fastxseqindex_s {
  seqaccesspoint_t *ap;
  Uint size;
  Uint allocated;
} fastxseqindex_t;


typedef struct fastxfileindex_s {
  seqaccesspoint_t *ap;
  Uint size;
  Uint allocated;
} fastxfileindex_t;


typedef struct fasta_s {

  CharSequence** seqs;
  CharSequence** quals;

  Uint *matestart;
  Uint noofseqs;
  Uint active_noofseqs;
  Uint active_noofmates;
  Uint minlen;
  Uint maxlen;
  char minqual;
  char maxqual;
  Uint curchunk;
  Uint offset;

  unsigned char lower;
  unsigned char upper;
  unsigned char gzip;
  unsigned char hasMates;
  unsigned char hasIndex;
  unsigned char chunkIsActive;
 
  Uint nooffiles;
  Uint *filetotal; 
  char **filenames;
  char **matefilenames;  

  fastxseqindex_t *chunkindex;
  fastxseqindex_t *matechunkindex;

  fastxfileindex_t **findex;
  fastxfileindex_t **matefindex;
  
  struct access **gzindex;
  struct access **mategzindex;

} fasta_t;




typedef struct {

  unsigned char type;
  /*
   *
   * chromname is gff seqname
   *
   */
  char *chromname; 
  Uint chromnamelen;

  /*
   * BED und GFF: 1-offset
   * for personalSNP: start with 0-offset
   * end base is not part of the feature ie. if 
   * end = 100 last feature base is 99. see below.
   *
   */

  Uint start;
  Uint end;

  /* 
   * GFF: name is the feature key
   * BED: name is the name of the acutal feature
   * for personalSNP track name is alleles A,C,T,G separated by '/'. 
   * Leading '-' is indel: insertion if start-end=0 
   *
   * */
  
  char *name;     
  Uint namelen;
  double score;
  unsigned char strand;

  /*GFF fields*/
  unsigned char frame; 
  char *source; 
  Uint sourcelen;
  Uint noofattributes;
  char **attributes;
  Uint *attributelen;

  /*BED fields*/
  Uint thickStart;
  Uint thickEnd;
  Uint *itemRgb;
  Uint blockCount;
  Uint* blockSizes;
  Uint* blockStarts;
  Uint noofovl;
  Uint firstovl;
  Uint level;
  /*extension*/
  char **blockRefseqs;
  char *blockStrands;

  /*SNPitem*/
  Uint alleleCount;//number of alleles in name
  Uint *alleleFreq;//from comma separated list of number of observed alleles - if unkowns 0
  Uint *alleleScores;//from a comma separated list - if unkown 0

} annotationitem_t;



typedef struct {

  char *trackname;
  Uint tracknamelen;
  char *description;
  Uint descriptionlen;
  Uint noofitems;
  annotationitem_t *items;

} annotationtrack_t;


typedef struct {
  Uint start;
  Uint end;
} cds_t;

typedef struct {
  Uint start;
  Uint end;
  char *refchr;
  char strand;
  Uint noofcds;
  cds_t *cds;
} exon_t;

typedef struct {
  char *id;
  char direction;
  Uint noofexons;
  exon_t *exons;
  Uint startcodon;
  Uint stopcodon;
} gene_t;


typedef struct {
  Uint noofgenes;
  gene_t *genes;
} geneset_t;


fasta_t**
bl_fastxChopIndex(void *space, fasta_t *f, Uint pieces);
void bl_fastaGetClipPos (fasta_t *f, Uint elem, Uint *p5, Uint *p3);
void bl_fastaGetMateClipPos (fasta_t *f, Uint elem, Uint *p5, Uint *p3); 
Uint bl_fastaSoftClip (void *space, fasta_t *f, Uint elem, 
    char *p5, Uint p5len, Uint p5scr, char *p3, Uint p3len, Uint p3scr, Uint pAlen);
Uint bl_fastaMateSoftClip (void *space, fasta_t *f, Uint elem, 
    char *p5, Uint p5len, Uint p5scr, char *p3, Uint p3len, Uint p3scr, Uint pAlen);
Uint bl_fastaHardClip (void *space, fasta_t *f, Uint elem, 
    Uint p5, Uint p3);
Uint bl_fastaMateHardClip (void *space, fasta_t *f, Uint elem, 
    Uint p5, Uint p3);
void bl_fastaDestruct(void *space, fasta_t* f);
void bl_fastxDestructSequence(void *space, fasta_t* f);
int bl_fastxGetChunk (fasta_t *fasta, Uint k);
fasta_t** bl_fastaChop(void *space, fasta_t* f, Uint pieces);
Uint bl_fastaGetMateDescriptionLength(fasta_t *f, Uint elem);
char* bl_fastaGetMateDescription(fasta_t *f, Uint elem);
unsigned char bl_fastaHasQuality(fasta_t *f);
fasta_t* bl_fastaInit(void *);
void bl_fastaAdd(void *space, fasta_t*, char *desc, Uint, 
    char* sequence, Uint, Uint);
fasta_t* bl_fastaRead(void *space, fasta_t*, char* filename, 
    unsigned char upper, unsigned char lower, unsigned int n, 
    void (*handler) (void *, fasta_t*, char*, Uint, char*, Uint, Uint));
fasta_t* bl_fastxGetSet(void *space, char **filenames, unsigned int nooffiles,
    unsigned char upper, unsigned char lower, unsigned char index, Uint pieces);
fasta_t* bl_fastaGetSet(void *space, char **filenames, unsigned int nooffiles,
    unsigned char upper, unsigned char lower);
Uint bl_fastaGetMateStart(fasta_t *f, Uint elem);
void bl_fastxDestructChunkIndex (void *space, fasta_t *);
fasta_t* bl_fastxGetMateSet(void *space, fasta_t* set, char** filenames, 
    unsigned int nooffiles, unsigned char upper, unsigned lower, 
    unsigned char index, Uint pieces);
void bl_fastxDestructSequence(void *space, fasta_t* f);
fasta_t* bl_fastxgzRead (void *space, fasta_t *fasta, char *filename, struct access *idx,
    unsigned char upper, unsigned char lower, off_t offset, Uint startseq, Uint lastseq,
    void (*handler)(void *, fasta_t*, char *, Uint, char *, char *, Uint, Uint));
fasta_t* bl_fastxRead(void *space, fasta_t* fasta, char* filename, 
    unsigned char upper, unsigned char lower, off_t offset, Uint startseq, Uint lastseq, 
    void (*handler) (void *, fasta_t*, char *, Uint, char *, char *, Uint, Uint));
void bl_fastxAdd(void *space, fasta_t *f, char *desc, Uint descrlen,
    char *sequence, char *quality, Uint seqlen, Uint sNo);
void bl_fastxAddMate(void *space, fasta_t *f, char *desc, Uint desclen,
    char *sequence, char *quality, Uint seqlen, Uint sNo);
void bl_fastxDump( void *space, fasta_t *fasta, char *desc, Uint desclen, 
    char *sequence, char *quality, Uint quallen, Uint sNo);
Uint bl_fastxGetChunkElem (void *space, fasta_t *f, Uint k);
fasta_t* bl_fastxIndex(void *space, fasta_t *f, char **filenames, Uint nooffiles, 
    unsigned char isMate, unsigned char gzip, Uint pieces);
Uint bl_fastaGetDescriptionLength(fasta_t *f, Uint elem);
char* bl_fastaGetDescription(fasta_t *f, Uint elem);
Uint bl_fastaGetSequenceLength(fasta_t *f, Uint elem);
char* bl_fastaGetSequence(fasta_t *f, Uint elem);
char* bl_fastaGetQuality(fasta_t* f, Uint elem);
#ifdef HASHING
Uint bl_fastaGetQuantity(fasta_t* f, Uint elem);
void bl_fastaSetQuantity(fasta_t* f, Uint elem, Uint quantity);
#endif
unsigned char bl_fastaHasMate(fasta_t *f);
Uint bl_fastaGetMateLength(fasta_t *f, Uint elem);
char* bl_fastaGetMate(fasta_t *f, Uint elem);
char* bl_fastaGetMateQuality(fasta_t *f, Uint elem);
fasta_t* bl_fastxgzIndex(void *space, char *gzfilename);
fastxseqindex_t* bl_fastxChunkIndex (void *space, char **filenames, struct access **gzindex, 
    fastxfileindex_t **findex, Uint *n, Uint nooffiles, Uint total, Uint k);
void bl_fastaSetMateClip (fasta_t *f, Uint elem, Uint p5, Uint p3);
void bl_fastaSetClip (fasta_t *f, Uint elem, Uint p5, Uint p3);
int bl_rm(void *space, char *filename);
Uint bl_fastxFindIDIdx (char *id, fasta_t *set);
annotationtrack_t* bl_BEDread (void *space, char *filename);
void bl_BEDwrite (annotationtrack_t *track, FILE *dev);
void bl_annotationtrackDestruct (void *space, annotationtrack_t *track);
annotationtrack_t* bl_GFFread (void *space, char *filename);
Uint bl_annotationitem_cmp_track (Uint item, void *track, void *elem, void *nfo);
int bl_fastxIDcmp (char *a, char *b);
void bl_GFFAddAttribute (void *space, annotationitem_t *item, char *attr, Uint len);
void bl_GFFwrite(char *filename, annotationtrack_t *set);
Uint bl_annotationtrackGetStats (void *space, annotationtrack_t *track);

#endif
