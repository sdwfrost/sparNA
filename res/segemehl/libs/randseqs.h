#ifndef _RANDSEQS_H_
#define _RANDSEQS_H_

/*
 *
 *	randseqs.h
 *  
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 07/15/2010 10:42:14 AM CEST  
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "stringutils.h"
#include "basic-types.h"
#include "charsequence.h"
#include "randseqs.h"

  geneset_t *
bl_getGeneModelFromBEDtrack (void *space, annotationtrack_t *track);

void
bl_copyGene (void *space, gene_t *to, gene_t *from);

annotationtrack_t *
bl_getTrackFromGeneModel (void *space, geneset_t *set);

  geneset_t*
bl_simulateTransSplicing (void *space, geneset_t *set, char type, Uint n);

  void
bl_printSplicingEdges (void *space, FILE *dev, geneset_t *set);

  char*
bl_getGeneSequence(void *space, fasta_t *reference, gene_t *gene);

  void
bl_simulateGeneSequencing (void *space, FILE *dev, fasta_t *reference, gene_t *gene, Uint readlen, 
    Uint cov, char *alphabet, Uint alphabetsize, Uint minqual, Uint maxqual,
    double acc, double Pmis, double Pins);

  void
bl_simulateGeneSetSequencing (void *space, FILE *dev, fasta_t *reference, geneset_t *genes, 
    Uint readlen, Uint cov, char *alphabet, Uint alphabetsize, Uint minqual, Uint maxqual,
    double acc, double Pmis, double Pins);

Uint
bl_fastxScramble (char *buffer, char *quality, 
    char *template, Uint len, 
    double acc, double Pmis, double Pins, 
    Uint uoff, Uint voff, 
    char *alphabet, Uint alphabetsize, 
    Uint minqual, Uint maxqual,
    char *editstring, Uint *editstringlen, Uint *readerrcnt);

 void
bl_fastxSimulateSpliceSites (void *space, char *sequence,
    Uint seqlen, Uint n, Uint maxchildren, 
    char *alphabet, Uint alphabetsize, double acc, 
    double Pmis, double Pins, double Pdel,
    Uint minqual, Uint maxqual,
    double Pcis, Uint mincisdist, Uint maxcisdist, 
    double Pstrandswitch, Uint readlen);

void 
bl_fastxPrintRandomReads(FILE *ref, char *sequence, Uint reflen, Uint n, 
    Uint minlen, Uint maxlen, char *alphabet, Uint alphabetsize,
    double acc, double Pmis, double Pins, double Pdel, unsigned char fastq, 
    Uint minqual, Uint maxqual, 
    char *five, Uint fivelen, char *three, Uint threelen, Uint polyAlen);


void
bl_fastxPrintRandomSplitReads (FILE *dev, char *sequence, Uint seqlen, Uint n, 
    Uint minspltlen, Uint maxspltlen, 
    char *alphabet, Uint alphabetsize,
    double acc,
    double Pmis, double Pins, double Pdel,
    unsigned char fastq,
    Uint minqual, Uint maxqual, 
    char *five, Uint fivelen,
    char *three, Uint threelen, Uint polyAlen);

  
void
bl_fastxPrintRandomMatePairs(FILE *dev, FILE *matedev, 
    char *sequence, Uint seqlen, Uint n, 
    Uint minlen, Uint maxlen, Uint mindist, Uint maxdist,  
    char *alphabet, Uint alphabetsize,
    double acc,
    double Pmis, double Pins, double Pdel,
    unsigned char fastq,
    Uint minqual, Uint maxqual,
    char *five, Uint fivelen,
    char *three, Uint threelen, Uint polyAlen);

void 
bl_fastxPrintRandomBisulfiteReads(FILE *ref, char *sequence, char *sequencerates, Uint reflen, 
    Uint n, Uint minlen, Uint maxlen, char *alphabet, Uint alphabetsize,
    double acc, double Pmis, double Pins, double Pdel, unsigned char fastq, 
    Uint minqual, Uint maxqual, 
    char *five, Uint fivelen, char *three, Uint threelen, Uint polyAlen, char *prefix);

Uint
bl_fastxBisulfiteScramble (char *buffer, char *quality, 
                           char *template, char *rates, Uint len, 
                           double acc, double Pmis, double Pins, 
                           Uint uoff, Uint voff, 
                           char *alphabet, Uint alphabetsize, 
                           Uint minqual, Uint maxqual,
                           char *editstring, Uint *editstringlen, Uint *readerrcnt);

#endif
