

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

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include "basic-types.h"
#include "info.h"
#include "debug.h"
#include "sort.h"
#include "vtprogressbar.h"
#include "mathematics.h"
#include "biofiles.h"
#include "container.h"
#include "matchfilesfields.h"
#include "matchfiles.h"
#include "evalmatchfiles.h"
#include "evalmethylmatchfiles.h"
#include "matfile.h"


/*-------------------------- bl_matchfileGetBSCross ----------------------------
 *    
 * @brief for each cross section in the frame: 
 * Get the count for all nucleotides
 * @author Christian Otto
 *   
 */
matchfileCross_t *
bl_matchfileGetBSCross(matchfileCross_t *cs) {
  Uint i, j, *space=NULL, *cnt, max;
  matchfileCross_t *bscs;

  bscs = ALLOCMEMORY(space, NULL, matchfileCross_t, 1);
  memset(bscs, 0, sizeof(matchfileCross_t));

  if (cs->ref != 'C' && cs->ref != 'G')
    return bscs;

  bscs->chars = ALLOCMEMORY(space, NULL, char, cs->len);
  bscs->quals = ALLOCMEMORY(space, NULL, char, cs->len);
  bscs->strands = ALLOCMEMORY(space, NULL, char, cs->len);
  bscs->feat = ALLOCMEMORY(space, NULL, char, cs->len);
  bscs->readpos = ALLOCMEMORY(space, NULL, uint32_t, cs->len);
  bscs->readlen = ALLOCMEMORY(space, NULL, uint32_t, cs->len);
  bscs->row = ALLOCMEMORY(space, NULL, uint32_t, cs->len);
  bscs->matchcnt = ALLOCMEMORY(space, NULL, uint32_t, cs->len);
  bscs->edist = ALLOCMEMORY(space, NULL, unsigned char, cs->len);
  bscs->bisulfite = ALLOCMEMORY(space, NULL, uint32_t, cs->len);

  /*
   * currently: no copy of dels (would be possible) & 
   *  splits (is useless) &
   *  matelinks (would need further modifications
   *  on input reading method)
   */
  for (i = 0, j = 0; j < cs->len; j++){
    if ((cs->ref == 'C' && cs->bisulfite[j] == 1) ||
        (cs->ref == 'G' && cs->bisulfite[j] == 2)){
      
      bscs->chars[i] = cs->chars[j];
      bscs->quals[i] = cs->quals[j];
      bscs->strands[i] = cs->strands[j];
      bscs->feat[i] = cs->feat[j];
      /* update: starts, ends */
      if (cs->feat[j] == '*'){
        bscs->starts++;
      }
      else if (cs->feat[j] == '$'){
        bscs->ends++;
      }
      bscs->readpos[i] = cs->readpos[j];
      bscs->readlen[i] = cs->readlen[j];
      bscs->row[i] = cs->row[j];
      /* update: maxrow */
      if(bscs->maxrow < cs->row[j]) {
        bscs->maxrow = cs->row[j];
      }
      bscs->edist[i] = cs->edist[j];
      bscs->bisulfite[i] = cs->bisulfite[j];
      i++;
    }
  }
  /* update: ref, len */
  bscs->ref = cs->ref;
  bscs->len = i;
  if (!bscs->len){    
    bl_matchfileDestructCross(space, bscs, 1);
    memset(bscs, 0, sizeof(matchfileCross_t));

    return bscs;
  }
  /* get consensus */
  bscs->cons = '^';
  cnt = bl_matchfileGetNTCounts(bscs);
  
  max = uarraymax(cnt, 256);
  bscs->cons = (char)max;

  /* realloc space */
  bscs->chars = ALLOCMEMORY(space, bscs->chars, char, bscs->len);
  bscs->quals = ALLOCMEMORY(space, bscs->quals, char, bscs->len);
  bscs->strands = ALLOCMEMORY(space, bscs->strands, char, bscs->len);
  bscs->feat = ALLOCMEMORY(space, bscs->feat, char, bscs->len);
  bscs->readpos = ALLOCMEMORY(space, bscs->readpos, uint32_t, bscs->len);
  bscs->readlen = ALLOCMEMORY(space, bscs->readlen, uint32_t, bscs->len);
  bscs->row = ALLOCMEMORY(space, bscs->row, uint32_t, bscs->len);
  bscs->matchcnt = ALLOCMEMORY(space, bscs->matchcnt, uint32_t, bscs->len);
  bscs->edist = ALLOCMEMORY(space, bscs->edist, unsigned char, bscs->len);
  bscs->bisulfite = ALLOCMEMORY(space, bscs->bisulfite, uint32_t, bscs->len);

  FREEMEMORY(space, cnt);
  
  return bscs;
}

/*--------------------------- bl_matchfileGetBSStrand --------------------------
 *    
 * @brief  get bisulfite-relevant strand with given reference base
 * @author Christian Otto
 *
 */
Uint
bl_matchfileGetBSStrand(char ref){
  Uint strand = 0;

  if (ref == 'C'){
    strand = PLUS_STRAND;
  }
  else if (ref == 'G'){
    strand = MINUS_STRAND;
  }
  else {
    strand = 0;
  }
  return strand;
}


/*--------------------------- bl_matchfileGetBSBase ----------------------------
 *    
 * @brief  get converted or unconverted bisulfite-related bases
 *         w/r/t given strand
 * @author Christian Otto
 *
 */
char
bl_matchfileGetBSBase(Uint strand, unsigned char conv){
  if (!strand)
    return 'N';

  if (!conv){
    return (strand == PLUS_STRAND) ? 'C' : 'G';    
  }
  else {
    return (strand == PLUS_STRAND) ? 'T' : 'A';
  }
}


/*--------------------------- bl_matchfileGetBSCount ---------------------------
 *    
 * @brief  get count of converted or unconverted bisulfite-related bases
 *         w/r/t given strand
 * @author Christian Otto
 *
 */
Uint
bl_matchfileGetBSCount(matchfileCross_t *bscs, Uint strand, unsigned char conv){
  Uint *cnt, ret;
  char ch;

  if (!strand) return 0;

  cnt = bl_matchfileGetNTCounts(bscs);
  ch = bl_matchfileGetBSBase(strand, conv);
  ret = cnt[(Uint) ch];

  free(cnt);
  return ret;
}

/*------------------------- bl_matchfileGetBSRateSimple ------------------------
 *    
 * @brief  get simple methylation rate estimate
 * @author Christian Otto
 *
 */
double
bl_matchfileGetBSRateSimple(matchfileCross_t *bscs, Uint strand, unsigned char allowdel){
  Uint i, j, max, meth, methch, unmeth, unmethch, *cnt;
  double rate = -1;

  if (!strand) return rate;

  cnt = bl_matchfileGetNTCounts(bscs);
  if (!allowdel) cnt[(int)'-'] = 0;

  /* get unconverted/methylated and converted/unmethylated base */
  methch = (Uint) bl_matchfileGetBSBase(strand, 0);
  meth = cnt[methch];
  unmethch = (Uint) bl_matchfileGetBSBase(strand, 1);
  unmeth = cnt[unmethch];

  /* calculate methylation rate */
  if (meth + unmeth > 0){
    i = uarraymax(cnt, 255);
    max = cnt[i];
    cnt[i] = 0;
    j = uarraymax(cnt, 255);
      
    /* Two valid cases:
     * case 1: 
     * unique maximal occurring character in cross
     * which is either the methylated or unmethylated
     * character
     * case 2:
     * ambiguous maximal occurring character in cross
     * which must be the methylation and unmethylated
     * character
     */       
    if ((max > cnt[j] && (i == methch || i == unmethch)) ||
        (max == cnt[j] && ((i == methch && j == unmethch)||
                           (i == unmethch && j == methch)))){
      rate = (double) meth/((double) meth + (double) unmeth);
    }
  }
  free(cnt);
  return rate;
}

/*------------------------ bl_matchfileCallMethylSimple ------------------------
 *    
 * @brief  Simple methylation caller
 * @author Christian Otto
 *
 */
Uint
bl_matchfileCallMethylSimple (void *space, Uint fidx, Uint cidx, Uint pos, matchfileCross_t *cs, 
                              char ref, matchfileindex_t *idx, unsigned char show, void *nfo )
{
  Uint i, strand, cov, methcov, meth, unmeth, other;
  double rate;
  char strandch, missing, *chrom, *context;
  matchfileCross_t *bscs;
  matfile_t *matfile = (matfile_t *) nfo;

  /* basic settings */
  chrom = matfile->files[fidx]->index->chromnames[cidx];
  missing = strandch = '.';

  /* get bisulfite conversion strand */
  strand = bl_matchfileGetBSStrand(ref);
  if (strand == PLUS_STRAND){
    strandch = '+';
  }
  else if (strand == MINUS_STRAND){
    strandch = '-';
  }

  /* do not process non-methylation sites */
  if (!strand){
    return 0;
  }

  /* get bisulfite cross */
  bscs = bl_matchfileGetBSCross(cs);

  /* get sequence context */
  i = bl_fastxFindIDIdx(chrom, matfile->fasta);
  context = bl_matchfileGetContext(matfile->fasta, i, pos, strand, 2);
  
  /* get coverage (DP) and methylation coverage (MDP) */
  cov = bl_matchfileGetCov(cs, 0);
  methcov = bl_matchfileGetCov(bscs, 0);

  /* get base counts in bisulfite cross */
  meth = unmeth = other = 0;
  rate = -1;

  if (methcov > 0){
    assert(strand);

    /* count unconverted/methylated and converted/unmethylated bases */
    meth = bl_matchfileGetBSCount(bscs, strand, 0);
    unmeth = bl_matchfileGetBSCount(bscs, strand, 1);
    other = methcov - meth - unmeth;

    /* calculate methylation rate */
    rate = bl_matchfileGetBSRateSimple(bscs, strand, 0);
  }

  /* report output */
  if (rate != -1){
    /* output first fields */
    fprintf(matfile->dev, "%s\t%d\t%c\t%c\t%c\t%c\t%c", chrom, pos, missing, ref, missing, missing, missing);
    fprintf(matfile->dev, "\t");
    /* info field */
    fprintf(matfile->dev, "CS=%c;CC=%s;NS=1;MMR=%.2f;DMR=.", strandch, context, rate);
    fprintf(matfile->dev, "\t");
    /* format field */
    fprintf(matfile->dev, "DP:MDP:MDP3:MRDP:CM:CU:MR");
    fprintf(matfile->dev, "\t");
    /* data field */
    fprintf(matfile->dev, "%d:%d:%d,%d,%d:%d:%d:%d:%.2f", cov, methcov, meth, unmeth, other,
            meth + unmeth, meth, unmeth, rate);
    fprintf(matfile->dev, "\n");
  }

  /* destruct everything */
  bl_matchfileDestructCross(space, bscs, 1);
  FREEMEMORY(space, bscs);
  FREEMEMORY(space, context);
  
  return 0;
}

Uint
bl_matchfileCalcMethylBias (void *space, Uint fidx, Uint cidx, Uint pos, matchfileCross_t *cs, 
                            char ref, matchfileindex_t *idx, unsigned char show, void *nfo )
{
  //matchfileSampleStats_t *stats = idx->stats;
  Uint i, strand, methcov, methbase, unmethbase, base;
  char *chrom, *context;
  matchfileCross_t *bscs;
  matfile_t *matfile = (matfile_t *) nfo;

  /* get bisulfite conversion strand */
  strand = bl_matchfileGetBSStrand(ref);

  /* do not process non-methylation sites */
  if (!strand){
    return 0;
  }

  /* get bisulfite cross */
  bscs = bl_matchfileGetBSCross(cs);

  /* get methylation coverage (MDP) */
  methcov = bl_matchfileGetCov(bscs, 0);

  if (methcov > 0){
    
    /* get sequence context */
    chrom = matfile->files[fidx]->index->chromnames[cidx];
    i = bl_fastxFindIDIdx(chrom, matfile->fasta);
    context = bl_matchfileGetContext(matfile->fasta, i, pos, strand, 2);
  
    /* TODO: parameter for sequence context filtering (given via nfo) */
    /* TODO: coverage filter */
    /* do not process nonCpG sites */
    if (strcmp(context, "CG") == 0){

      methbase = bl_matchfileGetBSBase(strand, 0);
      unmethbase = bl_matchfileGetBSBase(strand, 1);

      for (i = 0; i < bscs->len; i++){
  
        base = bscs->chars[i];

        /* exclude deleted chars (TODO: same as in BSmooth?) */
        if (base != '-'){

          /* report output */
          fprintf(matfile->dev, "%u\t%u\t%u\t%u\n", //pos, ref, 
                  bscs->readpos[i],
                  (base == methbase) ? 1 : 0,
                  (base == unmethbase) ? 1 : 0,
                  (base != methbase && base != unmethbase) ? 1 : 0);                  
        }
      }     
    }
    FREEMEMORY(space, context);
  }
    
  /* destruct everything */  
  bl_matchfileDestructCross(space, bscs, 1);
  FREEMEMORY(space, bscs);

  return 0;
}

/*------------------ bl_evalmatchfileSampleCrossSectionsBS -------------------
 *    
 * @brief sample and execute f on it
 * @author Steve Hoffmann 
 *   
 */
 
int
bl_matchfileSampleCrossSectionsBS(void *space, 
				matchfile_t *file, fasta_t *set, Uint n, 
				void (*f)(void *, matchfileFrame_t*, Uint, 
					  matchfileFrameStats_t *, void *), void *info)
{
  PairUint *samplepos;
  Uint i=0, r=0, j=0, k, *cumchrlen, 
    *order, prev=0, nchr, curchrom=0, curstart=0,
    *mapsize = NULL;
  matchfileFrame_t *frame = NULL;   
  matchfileFrameStats_t *stats = NULL;
  Uint maxcover = 20000;
  Uint setsize = 10000000;
  unsigned char **maps;

  char *sequence;
  char *pch;
  char next;
  Container *positions;
  Lint pos;


  //init random number generator
  srand((unsigned)(time(0))); 
  nchr = file->index->noofchroms;

  samplepos = ALLOCMEMORY(space, NULL, PairUint, n+1);
  memset(samplepos, 0, sizeof(PairUint)*n+1);
  cumchrlen = ALLOCMEMORY(space, NULL, Uint, nchr);
  memset(cumchrlen, 0, sizeof(Uint)*nchr);

  MSG("generating small map\n");
  //sum up the length of the references (i.e. chromosomes)
  maps = bl_matchfileSmallMap (space, file, &mapsize);
  MSG("map generated\n");

  cumchrlen[0] = file->index->matchend[0] - file->index->matchstart[0] + 1;
  fprintf(stderr, "%u: cumchrlen %u\n", 0, cumchrlen[0]);
  for(i=1; i < nchr; i++) {
    assert(file->index->matchend[i] >= file->index->matchstart[i]);
    cumchrlen[i] = (file->index->matchend[i] - 
		    file->index->matchstart[i]) + cumchrlen[i-1];
    fprintf(stderr, "%u: cumchrlen %u\n", i, cumchrlen[i]);
  }

  //  HELENE
  positions = ALLOCMEMORY(space, NULL, Container, 1);
  bl_containerInit(positions, 10000, sizeof(Lint));

  for (i=0; i<nchr; i++){
    sequence = bl_fastaGetSequence(set, i);
    printf("%s\t%d\n", bl_fastaGetDescription(set, i), bl_fastaGetSequenceLength(set, i));
    
    //position C on plus
    pch = strchr(sequence, 'C');
    while (pch!=NULL){

      next = sequence[pch-sequence+1];
      if (next == 'G'){
	if (i == 0){
	  pos = pch-sequence;
	  bl_containerAdd(positions, &pos);
	}
	else{
	  pos = pch-sequence+cumchrlen[i-1];
	  bl_containerAdd(positions, &pos);
	}
      }      
      pch=strchr(pch+1,'C');
    }
    //    printf("Start-End: %u\t%u\n", file->index->matchstart[i],  file->index->matchend[i]);
  }

  //  printf ("Anzahl CG: %u\n", bl_containerSize(positions));

  //randomize n positions across the genome and deterimine their 
  //chromosomes
  i = 0;
  j = 0;
 
  while(i < n && j < setsize) {
    k=0;
    
    pos = RANDINT(bl_containerSize(positions)-1);
    samplepos[i].b = * ((Uint *) bl_containerGet(positions, pos));
    //    printf("%d\n", samplepos[i].b);
    
    while(samplepos[i].b > cumchrlen[k] && k+1 < nchr) k++;   
    samplepos[i].a = k;
    //    printf("%d\n\n", samplepos[i].a);

    prev = (k == 0) ? 0 : cumchrlen[k-1];

    if(maps[samplepos[i].a] 
       && mapsize[samplepos[i].a] > (samplepos[i].b - prev)/255 
       && maps[samplepos[i].a][(samplepos[i].b - prev)/255] > 200) {
      i++;
      r++;
    }

    j++;
  }

  NFO("\n selected %d positions for sampling %d %d\n", i, j, n);
 
  for(i=0; i < nchr; i++) {
    if(maps[i]) { 
      FREEMEMORY(space, maps[i]);
    }
  }

  bl_containerDestruct(positions, NULL);
  FREEMEMORY(space, positions);
  FREEMEMORY(space, maps);
  FREEMEMORY(space, mapsize);
 
  if(j == setsize && r < (int)(0.8*((double)n))) {
    DBG("current sample size %d is below the minimum\n", r);
    FREEMEMORY(space, samplepos);
    FREEMEMORY(space, cumchrlen);
    return 0;
  }

  //sort
  order = quickSort(space, samplepos, n, bl_matchfileSampleCmp, NULL);   

  initProgressBarVT();
  
  //evaluate
  //to increase speed a frame of size FRAMESIZE is loaded
  for(i=0; i < n; i++) {

  
    progressBarVT("positions sampled.", n, i, 25);
    //is position on a new chromosome or in a new frame?
    if(!frame || samplepos[order[i]].a != curchrom || 
       samplepos[order[i]].b-prev+1 >= frame->start+frame->width) {

      if(frame) { 
        bl_matchfileDestructFrame(space, frame);
        frame = NULL;
        //bl_matchfileDestructFrameStats(space, stats);
      }

      curchrom = samplepos[order[i]].a;
      curstart = samplepos[order[i]].b;
      prev = (samplepos[order[i]].a == 0) ? 0 : cumchrlen[samplepos[order[i]].a-1];

      fprintf(stderr, "getting frame for '%s', curstart '%d', prev '%d'\n", 
      	      file->index->chromnames[samplepos[order[i]].a], curstart, prev);
      
      frame = bl_matchfileGetFrame(space, file, 
				   file->index->chromnames[samplepos[order[i]].a], 
				   curstart-prev+1, 20000, set, maxcover, NULL); 
     
      fprintf(stderr, "getting consensus\n" );
      bl_matchfileGetConsensus(frame);
      // stats = bl_matchfileFrameStats (space, frame);
    }
    
    fprintf(stderr, "evaluation of %d\n", samplepos[order[i]].b-curstart);
    f(space, frame, samplepos[order[i]].b-curstart, stats, info);   
  }

  NFO("\n %d positions sampled.\n", n);
    
  if(frame) { 
    bl_matchfileDestructFrame(space,frame);
    frame = NULL;
    //    bl_matchfileDestructFrameStats(space, stats);
  }

  FREEMEMORY(space, order);
  FREEMEMORY(space, samplepos);
  FREEMEMORY(space, cumchrlen);
  return 1;
}


/*------------------------- bl_matchfileSampleStatsBS --------------------------
 *    
 * @brief get bisulfite sample statistics
 * @author Helene 
 *   
 */

void
bl_matchfileSampleStatsBS(void *space, matchfileFrame_t *frame, 
                          Uint pos, matchfileFrameStats_t *framestats, void *nfo)
{  
  double rate;//, er=.0, b;
  //Uint i;
  Uint methcov, strand;
  matchfileCross_t *bscs;
//  matchfileSampleStats_t *stats = (matchfileSampleStats_t*) nfo; 
 
  /* get methylation coverage */
  bscs = bl_matchfileGetBSCross(&frame->cs[pos]);
  methcov = bl_matchfileGetCov(bscs, 0);

  printf("BS_cov\t%d\n", methcov);

  if (methcov < 15) 
    return;

  /* calculate methylation rate */
  strand = bl_matchfileGetBSStrand(bscs->ref);
  rate = bl_matchfileGetBSRateSimple(bscs, strand, 0);
  printf("BS_rate\t%lf\n", rate);

  /* TODO: rewrite functions
  int *positions = bl_matchfileGetReadPosBS(space, frame, pos);
  printf("BS_pos\t");
  for (i = 0; i < methcov; i++){
    printf("%d\t", positions[i]);
  }
  printf("\n");
  FREEMEMORY(space, positions);

  int *qualities = bl_matchfileGetReadQualBS(space, frame, pos);
  printf("BS_quals\t");
  for (i = 0; i < methcov; i++){
    printf("%d\t", qualities[i]);
  }
  printf("\n");
  FREEMEMORY(space, qualities);
  */

/*    b = bl_matchfileGetStrandBias(frame, pos);

  if(e > 0 && stats->e_N < stats->n) {   
    stats->entropy[stats->e_N] = frame->cs[pos].longentropy;
    stats->eraw[stats->e_N]=e;
    stats->b[stats->e_N]=b;
    stats->e[stats->e_N++]=e-er;
    }
*/
  return;
}

/*---------------------- bl_matchfileGetReadPosBS -----------------------
 *    
 * @brief get the positions of a read that cover a cross section
 * (not reference)
 * @author Helene
 *   
 */
/*
double*
bl_matchfileGetReadPosBS (void *space, matchfileFrame_t *frame, Uint pos)
{
  Uint j;
  Uint i=0;
  double *r;
  matchfileCross_t *bscs, *save;
  //TODO: rewrite function with use of following function
  bscs = bl_matchfileGetBSCross(&frame->cs[pos]);
  
  
  positions = ALLOCMEMORY(space, NULL, int, frame->cs[pos].len);

  if (frame->cs[pos].len){
    for(j=0; j < frame->cs[pos].len; j++) {
      positions[i] = -1;
      if (frame->cs[pos].ref == 'C'){
	if (frame->cs[pos].bisulfite[j] == 1){
	  if (frame->cs[pos].chars[j] == 'C' || frame->cs[pos].chars[j] == 'T'){
	    positions[i] = frame->cs[pos].readpos[j];
	    i++;
	  }
	}
      }
      else if (frame->cs[pos].ref == 'G'){
	if (frame->cs[pos].bisulfite[j] == 2){
	  if (frame->cs[pos].chars[j] == 'G' || frame->cs[pos].chars[j] == 'A'){
	    positions[i] = frame->cs[pos].readpos[j];
	    i++;
	  }
	}
      }
    }
  }

  return positions;
}
*/

/*---------------------- bl_matchfileGetReadPosBS -----------------------
 *    
 * @brief get the qualities of the reads that cover a cross section
 * (not reference)
 * @author Helene
 *   
 */
/*
int *
bl_matchfileGetReadQualBS (void *space, matchfileFrame_t *frame, Uint pos)
{
  Uint j;
  Uint i=0;
  int *quality;
  quality = ALLOCMEMORY(space, NULL, int, frame->cs[pos].len);

  //TODO: rewrite function with use of following function
  bscs = bl_matchfileGetBSCross(&frame->cs[pos]);

  if (frame->cs[pos].len){
    for(j=0; j < frame->cs[pos].len; j++) {
      quality[i] = -1;
      if (frame->cs[pos].ref == 'C'){
	if (frame->cs[pos].bisulfite[j] == 1){
	  if (frame->cs[pos].chars[j] == 'C' || frame->cs[pos].chars[j] == 'T'){
	    quality[i] = frame->cs[pos].quals[j];
	    i++;
	  }
	}
      }
      else if (frame->cs[pos].ref == 'G'){
	if (frame->cs[pos].bisulfite[j] == 2){
	  if (frame->cs[pos].chars[j] == 'G' || frame->cs[pos].chars[j] == 'A'){
	    quality[i] = frame->cs[pos].quals[j];
	    i++;
	  }
	}
      }
    }
  }

  return quality;
}
*/
