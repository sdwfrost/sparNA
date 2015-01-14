#ifndef MATFILE_H
#define MATFILE_H

/*
 *
 *	matfile.h
 *  
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 08/25/2010 03:47:54 PM CEST  
 *
 */

#include "biofiles.h"
#include "matchfiles.h"

typedef struct {
  FILE *dev;
  fasta_t *fasta;
  matchfile_t **files;
} matfile_t;


#endif
