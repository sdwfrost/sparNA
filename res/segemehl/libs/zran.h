#ifndef ZRAN_H
#define ZRAN_H
/*
 *
 *	zran.h
 *  
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 07/06/2010 03:24:28 PM CEST  
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include "zlib.h"

#define SPAN 1048576L       /* desired distance between access points */
#define WINSIZE 32768U      /* sliding window size */
#define CHUNK 16384         /* file input buffer size */
#define LARGECHUNK 1638400000
#define MEDIUMCHUNK 1638400


/* access point entry */
struct point {
    off_t out;          /* corresponding offset in uncompressed data */
    off_t in;           /* offset in input file of first full byte */
    int bits;           /* number of bits (1-7) from byte at in - 1, or 0 */
    unsigned char window[WINSIZE];  /* preceding 32K of uncompressed data */
};

/* access point list */
struct access {
    int have;           /* number of list entries filled in */
    int size;           /* number of list entries allocated */
    struct point *list; /* allocated list */
};

struct gzidxfile {
  FILE *fp;
  struct access *index;
  off_t curap;
  int mychunk;
  unsigned char *buf;
  unsigned char *pos;
  int len;
};

/* Deallocate an index built by build_index() */
void free_index(struct access *index);

int extract(FILE *in, struct access *index, off_t offset,
                  unsigned char *buf, int len);

struct access* bl_zranGetIndex(char *filename, int *len);

void bl_destructgzidxfile(struct gzidxfile *file);

off_t bl_ftellgzidx(struct gzidxfile *f);

int build_index(FILE *in, off_t span, struct access **built);

struct gzidxfile* bl_initgzidxfile(FILE *fp, struct access *index, off_t offset, int len);


int bl_getgzidxc (struct gzidxfile *f);


#endif
