
/*
 *  matfile.c
 *  match files
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 08/25/2010 03:42:37 PM CEST
 *  
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include "alignment.h"
#include "stringutils.h"
#include "basic-types.h"
#include "mathematics.h"
#include "bitVector.h"
#include "matchfiles.h"
#include "browsematchfiles.h"
#include "matfile.h"
#include "iupac.h"
#include "info.h"
#include "manopt.h"
#include "evalmatchfiles.h"
#include "evalmethylmatchfiles.h"
#include "splicesites.h"
#include "startsites.h"
#include "snvsplines.h"

unsigned char mute = 0;
char *ntcode;



/*---------------------------------- stats -----------------------------------
 *    
 * @brief get matchfile stats
 * @author Steve Hoffmann 
 *   
 */

  void
getstats (void *space, matchfile_t *file, fasta_t *fasta, Uint mincover, 
    Uint maxcover, double minfrac, char entropyfilter, char strandbias, char usenativequal, Uint samplecond, Uint samplescr)
{
  matchfileSampleStats_t *stats;
  unsigned char **maps=NULL;
  Uint *mapsize=NULL, i;
  double mx, sx, kx;
  Uint nchr = file->index->noofchroms;

  stats = bl_matchfileInitSampleStats(space, MAX(samplescr,samplecond)+20000, mincover, maxcover, minfrac, entropyfilter, 10, .15);
  stats->strand = strandbias;
  file->index->stats = stats;
  
  
  MSG("generating small map\n");
  maps = bl_matchfileSmallMap (space, file, &mapsize);
  
  
  MSG("evaluating error distribution.\n");
  bl_matchfileSampleCrossSections(space, file, fasta, 100000, 
      bl_matchfileGetErrorDensity, maps, mapsize, stats);
  
  qsort(stats->e, stats->e_N, sizeof(double), cmp_dbl_qsort);
  gevLmoment(stats->e, stats->e_N, &mx, &sx, &kx);
  gevmle(NULL, stats->e, stats->e_N, &mx, &sx, &kx, 100000, stats->e[0], stats->e[stats->e_N-1]);
  stats->pxx = mx;
 
  NFO("setting maximum error for sampling to:%f.\n", stats->pxx);

  
  MSG("sampling parameter.\n");
  bl_matchfileSampleCrossSections(space, file, fasta, samplecond, 
      bl_matchfileGetConditionals, maps, mapsize, stats);


  MSG("sampling scores.\n");
  stats->ecdf = ecdf_init(stats->eraw, stats->e_N); 
  
  bl_matchfileSampleCrossSections(space, file, fasta, samplescr, 
      bl_matchfileGetScoreSample, maps, mapsize, file->index);

//  getcutoff(stats, NULL, NULL, NULL, NULL, NULL);

 
  for(i=0; i < nchr; i++) {
    if(maps[i]) { 
      FREEMEMORY(space, maps[i]);
    }
  }

  FREEMEMORY(space, maps);
  FREEMEMORY(space, mapsize);

  return ;
}

/*----------------------------------- view -----------------------------------
 *    
 * @brief view the matchfiles
 * @author Steve Hoffmann 
 *   
 */

void
view (void *space, matchfile_t **files, Uint nooffiles, fasta_t *fasta, 
    annotationtrack_t *bed)
{
 
  MSG("starting viewer.\n");
  bl_matchfileViewer(space, files, nooffiles, fasta, bed, 100, 1000);
  return ;
}

/*----------------------------------- eval -----------------------------------
 *    
 * @brief evaluation of matchfiles
 * @author Steve Hoffmann 
 *   
 */
 
void
eval (void *space, matchfile_t **files, int *groups, Uint nooffiles, fasta_t *fasta, Uint maxframesize, double cut)
{

  MSG("evaluating x-sections.\n");

  if(cut != -1) {
    MSG("override cutoff\n");
    files[0]->index->stats->cut = cut;
  } else { 
    MSG("calculating cutoff\n");
    getcutoff(  files[0]->index->stats , NULL, NULL, NULL, NULL, NULL);
  }

  bl_matchfileEvalCrossSections(space, files, groups, nooffiles, fasta, maxframesize, 
      bl_matchfileTest, NULL);
 
  return ;
}
/*----------------------------------- simpleeval -----------------------------------
 *    
 * @brief simple evaluation of matchfiles
 * @author Steve Hoffmann 
 *   
 */
 
void
simpleeval (void *space, matchfile_t **files, int *groups, Uint nooffiles, fasta_t *fasta, Uint maxframesize, char usenativequal)
{

  
  MSG("simple evaluating x-sections.\n");
  
  bl_matchfileEvalCrossSections(space, files, groups, nooffiles, fasta, maxframesize, 
      bl_matchfileSimpleGEV, &usenativequal);
 
  return ;
}

/*----------------------------------- gatkeval -----------------------------------
 *    
 * @brief simple evaluation of matchfiles
 * @author Steve Hoffmann 
 *   
 */
 
void
gatkeval (void *space, matchfile_t **files, int *groups, Uint nooffiles, fasta_t *fasta, Uint maxframesize, char usenativequal)
{

  
  MSG("gatk evaluating x-sections.\n");
  
  bl_matchfileEvalCrossSections(space, files, groups, nooffiles, fasta, maxframesize, 
       bl_matchfileSimpleGATK , &usenativequal);
 
  return ;
}

/*---------------------------------- consensus ---------------------------------
 *    
 * @brief simple evaluation of matchfiles
 * @author Christian Otto
 *   
 */

void
evalconsensus (void *space, matchfile_t **files, int *groups, Uint nooffiles, fasta_t *fasta, Uint maxframesize)
{

  
  MSG("consensus evaluating x-sections.\n");
  
  bl_matchfileEvalCrossSections(space, files, groups, nooffiles, fasta, maxframesize, 
       bl_matchfileConsensus , files);
 
  return ;
}

/*------------------------------ callMethylSimple ------------------------------
 *    
 * @brief simple evaluation of matchfiles
 * @author Christian Otto
 *   
 */

void
callMethylSimple (void *space, matchfile_t **files, int *groups, Uint nooffiles, fasta_t *fasta, Uint maxframesize)
{
  assert(nooffiles == 1);

  Uint i, len;
  char *name, *basename;
  matfile_t matfile;
  matfile.dev = stdout;
  matfile.fasta = fasta;
  matfile.files = files;  
  
  /* prepare sample name from file name */
  basename = bl_basename(files[0]->filename);
  len = strlen(basename);
  name = ALLOCMEMORY(space, NULL, char, len + 1);
  memmove(name, basename, len + 1);

  if (files[0]->gzip){
    if (strncmp(&name[len-3], ".gz", 3) == 0){
      name[len-3] = '\0';
    }
    else {
      if (strncmp(&name[len-5], ".gzip", 5) == 0) {
        name[len-5] = '\0';
      }
    }  
    assert(strlen(name) != len);
    len = strlen(name);
  }
  len = bl_fileprefixlen(name);
  name[len] = '\0';
  
  MSG("evaluating x-sections.\n");
  
  /* output VCF header */
  fprintf(matfile.dev, "##fileformat=VCFv4.1\n");
  
  for (i = 0; i < fasta->nooffiles; i++){
    fprintf(matfile.dev, "##reference=file://%s\n",
            fasta->filenames[i]);
  }

  for (i = 0; i < fasta->noofseqs; i++){
    fprintf(matfile.dev, "##contig=<ID=%s,length=%d>\n", 
            bl_fastaGetDescription(fasta, i), 
            bl_fastaGetSequenceLength(fasta,i));
  }
  //fileDate? program call?
  fprintf(matfile.dev, "##INFO=<ID=CS,Number=1,Type=Character,Description=\"Bisulfite conversion strand, i.e., strand of cytosine relative to reference\">\n");
  fprintf(matfile.dev, "##INFO=<ID=CC,Number=1,Type=String,Description=\"Sequence context of two bases on the reference relative to bisulfite conversion strand, missing value if context is beyond reference boundaries\">\n");
  fprintf(matfile.dev, "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">\n");
  fprintf(matfile.dev, "##INFO=<ID=MMR,Number=1,Type=Float,Description=\"Mean methylation rate among samples with data\">\n");
  fprintf(matfile.dev, "##INFO=<ID=DMR,Number=1,Type=Float,Description=\"Difference in methylation rate (Sample1 - Sample2).\">\n");
  fprintf(matfile.dev, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth, not filtered by mapping quality, base quality, or bisulfite conversion information\">\n");
  fprintf(matfile.dev, "##FORMAT=<ID=MDP,Number=1,Type=Integer,Description=\"Number of reads after filtering by bisulfite conversion strand\">\n");
  fprintf(matfile.dev, "##FORMAT=<ID=MDP3,Number=3,Type=Integer,Description=\"Number of reads with 1) unconverted cytosines, 2) converted cytosines or thymines, and 3) other bases after filtering by bisulfite conversion strand\">\n");
  fprintf(matfile.dev, "##FORMAT=<ID=MRDP,Number=1,Type=Integer,Description=\"Number of reads with unconverted (methylated) and converted (unmethylated) cytosines after filtering by bisulfite conversion strand\">\n");
  fprintf(matfile.dev, "##FORMAT=<ID=CM,Number=1,Type=Integer,Description=\"Number of reads with unconverted (methylated) cytosines after filtering by bisulfite conversion strand\">\n");
  fprintf(matfile.dev, "##FORMAT=<ID=CU,Number=1,Type=Integer,Description=\"Number of reads with converted (unmethylated) cytosines after filtering by bisulfite conversion strand\">\n");
  fprintf(matfile.dev, "##FORMAT=<ID=MR,Number=1,Type=Float,Description=\"Estimated methylation rate\">\n");
  fprintf(matfile.dev, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n", name);

  bl_matchfileEvalCrossSections(space, files, groups, nooffiles, fasta, maxframesize, 
       bl_matchfileCallMethylSimple, &matfile);

  FREEMEMORY(space, name); 
  return ;
}

void
calcMethylBias (void *space, matchfile_t **files, int *groups, Uint nooffiles, fasta_t *fasta, Uint maxframesize)
{
  matfile_t matfile;
  matfile.dev = stdout;
  matfile.fasta = fasta;
  matfile.files = files;

  bl_matchfileEvalCrossSections(space, files, groups, nooffiles, fasta, maxframesize, 
       bl_matchfileCalcMethylBias, &matfile);
 
  return ;
}


/*-------------------------------- evalstarts --------------------------------
 *    
 * @brief get start sites
 * @author Steve Hoffmann 
 *   
 */
 
void
evalstarts (void *space, matchfile_t **files, int *groups, Uint nooffiles, 
    fasta_t *fasta, Uint maxframesize, annotationtrack_t *bed )
{
 
  Uint *cntr;
//  Uint i;
  cntr = ALLOCMEMORY(space, NULL, Uint, 255);
  memset(cntr, 0, sizeof(Uint)*255);
 
  MSG("start sites.\n");
/*
  bl_matchfileEvalCrossSections(space, files, groups, nooffiles, fasta, 
      bl_matchfileStartSites, cntr);

  for(i=0; i < 255; i++) {
    printf("%d\t%d\n", i, cntr[i]);
  }
*/
  
  MSG("coverage.\n");
  bl_matchfileEvalCrossSections(space, files, groups, nooffiles, fasta, maxframesize, 
      bl_coverage, cntr);

  return ;
}

/*-------------------------------- writewiggle -------------------------------
 *    
 * @brief get start sites
 * @author Steve Hoffmann 
 *   
 */
 
void
writeexpression (void *space, matchfile_t **files, int *groups, Uint nooffiles, 
    fasta_t *fasta, annotationtrack_t *bed)
{
 
  Uint i;
  char *filename;
  MSG("write expression file.\n");

  for(i=0; i < nooffiles; i++) { 
    filename = ALLOCMEMORY(space, NULL, char, strlen(files[i]->filename)+5);
    memmove(filename, files[i]->filename, strlen(files[i]->filename));
    sprintf(&filename[strlen(files[i]->filename)], ".wig");

    bl_writeexpression(space, filename, strlen(filename), files[i], fasta, 
        10000000, 1000000);
  }

  return ;
}



/*---------------------------------- splice ----------------------------------
 *    
 * @brief dump the splice sites
 * @author Steve Hoffmann 
 *   
 */
 
void
evalsplice (void *space, matchfile_t **files, int *groups, Uint nooffiles, fasta_t *fasta, Uint maxframesize, 
    annotationtrack_t *bed, char *filename, char *transfilename, Uint minsplitno)
{
  Uint i, cidx, pos, len;
  char *chr, *chr2, *ptr, *base;
  splitmap_t map;
  splicemap_t *sm;
  FILE *fp1, *fp2;

  bl_matchfileInitSplitMap(space, &map, bed, files, fasta);
  MSG("eval splice sites.\n");
  bl_matchfileEvalCrossSections(space, files, groups, nooffiles, fasta, maxframesize, 
      bl_matchfileSplit, &map);
  MSG("condensing sites.\n");
  sm = bl_matchfileSpliceMap (space, &map, 10, minsplitno);
  MSG("writing splice sites to stdout.\n");



  fp1 = fopen(filename, "w");
  if (fp1 == NULL) {
    fprintf(stderr, "Couldnt open %s for reading. Exit forced.\n", filename);
    exit(-1);
  } 


  fp2 = fopen(transfilename, "w");
  if (fp2 == NULL) {
    fprintf(stderr, "Couldnt open %s for reading. Exit forced.\n", transfilename);
    exit(-1);
  } 


  base = bl_basename(filename);

  //printsplice(space, sm, stdout);
  printsplicebed(space, sm, minsplitno, base, fp1, fp2);
 
  MSG("writing splice sites to gff.\n");
  if(bed) {
    bl_annotationtrackGetStats (space, bed);
    bl_matchfileSpliceAnnotation(space, sm, bed);
    bl_GFFwrite("splice.gff", bed);
  }

  for(i=0; i < map.noofsplits; i++) {
    cidx = map.cidx[i];
    chr = bl_fastaGetDescription(fasta, cidx);
    pos = map.pos[i];
 
    len = strlen(chr);
    chr2 = ALLOCMEMORY(space, NULL, char, len+1);
    memmove(chr2, chr, len);
    chr2[len] = 0;
    ptr = strtok(chr2, " ");
    
    // dummy if to avoid not used warnings
    if (0){
      printsplits(space, ptr, pos, &map.cs[i], &map);
    }
    FREEMEMORY(space, chr2);
  }

  bl_matchfileDestructSplitMap(space, &map);
  bl_matchfileDestructSpliceMap (space, sm);  
  FREEMEMORY(space, sm);

  fclose(fp1);
  fclose(fp2);

  return ;
}

/*----------------------------------- main -----------------------------------
 *    
 * @brief the main
 * @author Steve Hoffmann 
 *   
 */
 

int
main(int argc, char **argv) {
  
  void *space = NULL;
  
  manopt_optionset optset;
  manopt_arg *unflagged; 
  manopt_arg *dbfilenames;
  manopt_arg *queries;
  manopt_arg *indices;
  manopt_arg *grouplist;
  manopt_arg *bedfilenames;
  matchfile_t **files = NULL;
  annotationtrack_t *track = NULL;
  fasta_t *fasta = NULL;
//  matchfileindex_t *idx2 = NULL;
  Uint prefixlen=0; 
  Uint splicebasenamelen;
  Uint scoreplotbasenamelen;
  Uint minsplitno = 4;
  Uint mincover = 6;
  Uint maxqual = 64;
  Uint maxcover = 200;
  char entropyfilter = 0;
  double minfrac = 0.1;
  unsigned char gzip = 0, browse=0, call=0, saveindex=0, stats=0, dumpstats=0, starts=0, wiggle=0, strandbias=0, consensus=0, methylcall=0, methylbias=0;
  char version[]="0.1";
  char *splicebasename = NULL;
  char *scoreplotbasename = NULL;
  char usenativequal = 0, simplegev=0, gatk=0;
  int *groups = NULL;
  int i;
  Uint maxframesize = 100000;

  FILE *beddev = NULL;
  char *filename;
  char *transfilename;
  double minscore = -1;
  Uint samplescr = 200000;
  Uint samplecond = 60000;

  initIUPAC(1,1); 
  manopt_initoptionset(&optset, argv[0], NULL, 
      "Heuristic mapping of short sequences\n",
      "SEGEMEHL is free software for non-commercial use \n  (C) 2008 Bioinformatik Leipzig\n",
      version,
      "Please report bugs to steve@bioinf.uni-leipzig.de"); 
  manopt(&optset, LISTOPT, 1, 'd', "database", 
      "list of path/filename(s) of database sequence(s)", "<file> [<file> ...]", 
      NULL, NULL);
  manopt(&optset, LISTOPT, 1, 'q', "query", 
      "path/filename of alignment file", "<file> [<file> ...]", NULL, NULL); 
  manopt(&optset, LISTOPT, 0, 'i', "index", 
      "path/filename of db index", "[<file> ... ]", NULL, NULL);
  manopt(&optset, LISTOPT, 0, 'a', "annotation", 
      "path/filename of bed annotation", "[<bedfile> ... ]", NULL, NULL);
  manopt(&optset, LISTOPT, 0, 'x', "generate", 
      "generate db index and store to disk", "<file>", NULL, NULL); 
  manopt(&optset, LISTOPT, 0, 'g', "group", 
      "group number for all nput files (1-offset and at most #files groups)", "[<int> ...]", NULL, NULL); 
  manopt(&optset, REQSTRINGOPT, 0, 's', "splice", 
      "dump splice sites to <basename>", NULL, NULL, &splicebasename);
  manopt(&optset, FLAG, 0, 'S', "starts", 
      "dump start sites", NULL, NULL, &starts);
  manopt(&optset, FLAG, 0, 'b', "browse", 
      "start browser", NULL, NULL, &browse);
  manopt(&optset, FLAG, 0, 'c', "call",
      "variant caller", NULL, NULL, &call);
  manopt(&optset, FLAG, 0, 'w', "expression",
      "generate a expression graph file from the match files", NULL, NULL, &wiggle);
  manopt(&optset, FLAG, 0, 'Z', "dumpstats",
      "dump data stats", NULL, NULL, &dumpstats);
  manopt(&optset, REQSTRINGOPT, 0, 'R', "Rplots",
      "write files for R score plots <basename>", NULL, NULL, &scoreplotbasename);
  manopt(&optset, FLAG, 0, 'H', "stats", 
      "stats calculation", NULL, NULL, &stats);
  manopt(&optset, REQUINTOPT, 0, 'm', "minsplitno", 
      "minimum number of splits required to call a splice site", "<n>",
      NULL, &minsplitno);
  manopt(&optset, REQUINTOPT, 0, 'N', "mincover", 
      "minimum coverage to call sites", "<n>",
      NULL, &mincover);
  manopt(&optset, REQUINTOPT, 0, 'M', "maxcover", 
      "minimum coverage to call sites", "<n>",
      NULL, &maxcover);
  manopt(&optset, REQUINTOPT, 0, 'K', "maxframesize", 
      "maximum size of memory resident frame", "<n>",
      NULL, &maxframesize);
  manopt(&optset, REQDBLOPT, 0, 'F', "minfrac", 
      "minimum fraction of ALT alleles to call sites", "<f>",
      NULL, &minfrac);
  manopt(&optset, REQDBLOPT, 0, 'u', "minscore", 
      "override minimum score", "<f>", NULL, &minscore);
  manopt(&optset, FLAG, 0, 'E', "entropyfilter",
      "turn on entropy filter", NULL, NULL, &entropyfilter);
  manopt(&optset, FLAG, 0, 'G', "strand",
      "use pbinom strand", NULL, NULL, &strandbias);
  manopt(&optset, FLAG, 0, 'Q', "nativequal",
      "use native quals", NULL, NULL, &usenativequal);
  manopt(&optset, FLAG, 0, 'V', "simplegev",
      "use simple gev", NULL, NULL, &simplegev);
  manopt(&optset, FLAG, 0, 'A', "gatk", 
      "use GATK model", NULL, NULL, &gatk);
  manopt(&optset, FLAG, 0, 'C', "consensus", 
      "use simple consensus calling", NULL, NULL, &consensus);
  manopt(&optset, FLAG, 0, 'B', "methylcall", 
      "use simple methylation calling", NULL, NULL, &methylcall);
  manopt(&optset, FLAG, 0, 'Y', "methylbias", 
      "calculate methylation bias", NULL, NULL, &methylbias);
  manopt(&optset, REQUINTOPT, 0, 'O', "qoffset", 
      "quality offset", "<n>", NULL, &maxqual);
 manopt(&optset, REQUINTOPT, 0, 'U', "samplecond", 
      "sample size for conditionals", "<n>", NULL, &samplecond);
 manopt(&optset, REQUINTOPT, 0, 'X', "samplescr", 
      "sample size for scores", "<n>", NULL, &samplescr);

 

  unflagged = manopt_getopts(&optset, argc, argv);
  saveindex = manopt_isset(&optset, 'x', NULL);
  
  if(!(!manopt_isset(&optset, 'i', NULL) ^ !manopt_isset(&optset, 'x', NULL))) {
    manopt_help(&optset, "please give index filename using -i XOR -x option\n");
  } else if(unflagged->noofvalues > 1) { 
    manopt_help(&optset, "unknown argument(s)\n");
  }

  MSG("reading database sequences.\n"); 
  NFO("minsplitno set to %d\n", minsplitno);

  dbfilenames = manopt_getarg(&optset, 'd', "database");
  fasta = bl_fastxGetSet(space, dbfilenames->values, 
      dbfilenames->noofvalues, 1, 0, 0, 1);

  NFO("%d database sequences found.\n", fasta->noofseqs);
  MSG("reading query files.\n");

  queries = manopt_getarg(&optset, 'q', "query");
  if(queries->noofvalues > 30) {
    manopt_help(&optset, "currently no more than 30 query files allowed\n");
  }

  grouplist = manopt_getarg(&optset, 'g', "group");
  if(grouplist) {
    if(grouplist->noofvalues != queries->noofvalues) 
      manopt_help(&optset, "please provide a group name for each input file");

    groups = ALLOCMEMORY(space, NULL, Uint, grouplist->noofvalues);
    
    for(i=0; i < grouplist->noofvalues; i++) {
      groups[i] = atoi(grouplist->values[i]);
      if(groups[i] == 0)
        manopt_help(&optset, "please provide group numbers (int) > 0");
      if(groups[i] > queries->noofvalues)
        manopt_help(&optset, "please provide groupnumbers <= number of input files");
      NFO("found group number %d\n", groups[i]);
    }
  }

  if (methylcall || methylbias){
    if (queries->noofvalues > 1){
      manopt_help(&optset, "multiple query files are not supported for methylation analyses\n");
    }
  }

  if(saveindex) {
    indices = manopt_getarg(&optset, 'x', "generate");
  } else {
    indices = manopt_getarg(&optset, 'i', "index");
  }

  if(indices->noofvalues != queries->noofvalues) {
    manopt_help(&optset, "please provide an index file name for each query file\n");
  }

  ntcode  = getNTcodekey(space);
  files = ALLOCMEMORY(space, NULL, matchfile_t*, queries->noofvalues);

  for(i=0; i < queries->noofvalues; i++) {

    files[i] = ALLOCMEMORY(space, NULL, matchfile_t, 1);  
    files[i]->fmt = 0;
    files[i]->index = NULL;
    files[i]->filename = queries->values[i];

    prefixlen = bl_fileprefixlen(files[i]->filename);

    if(strncmp(&files[i]->filename[prefixlen], ".gz", 3) == 0 || 
        strncmp(&files[i]->filename[prefixlen], ".gzip", 5) == 0) {
      gzip = 1;
    }

    files[i]->gzip = gzip;

    if(saveindex) {

      bl_matchfileIndex(space, files[i], fasta);
      bl_matchfileWriteIndex(files[i]->index, indices->values[i]);  
   /*   if(stats) { 
        getstats(space, files[i], fasta, mincover, maxcover, minfrac, entropyfilter, strandbias, usenativequal); 
      }

      bl_matchfileWriteIndex(files[i]->index, indices->values[i]);  
      idx2 = bl_matchfileReadIndex(space, indices->values[i]);

      fprintf(stderr, "compare index (%p:%p):%d\n", 
          (void*)files[i]->index, (void *)idx2, 
          bl_compareIndices(files[i]->index, idx2));
      bl_matchfileDestructIndex(space, idx2); 
      FREEMEMORY(space, idx2);
   */
    } else if(indices->values[i]) {
    
      MSG("reading index file\n");
      files[i]->index = bl_matchfileReadIndex(space, indices->values[i]);
    }
   

    if(stats) { 
       
      /*if(files[i]->index->stats) {
         bl_matchfileDestructSampleStats(NULL, files[i]->index->stats); 
       }
      */ 
       getstats(space, files[i], fasta, mincover, maxcover, minfrac, entropyfilter, strandbias, usenativequal, samplecond, samplescr); 
       bl_matchfileWriteIndex(files[i]->index, indices->values[i]);   
      
       /*
       idx2 = bl_matchfileReadIndex(space, indices->values[i]);
       fprintf(stderr, "compare index (%p:%p):%d\n", 
           (void*)files[i]->index, (void *)idx2, bl_compareIndices(files[i]->index, idx2));
       bl_matchfileDestructIndex(space, idx2); 
       FREEMEMORY(space, idx2);
       */
    }

   
    /*
     *  typically stats should be present if stats
     *  failed however we are not allowed to set
     *  this one
     */
    
    if(files[i]->index->stats) { 
        if(files[i]->index->stats->mincover != mincover || 
            files[i]->index->stats->maxcover != maxcover || 
            files[i]->index->stats->minfrac != minfrac) {
          MSG("WARNING: resetting minfrac, mincover or maxcover!");
          NFO("mincover: %d, maxcover:%d, minfrac:%f\n",files[i]->index->stats->mincover, files[i]->index->stats->maxcover, files[i]->index->stats->minfrac);
        }
        files[i]->index->stats->mincover = mincover;
        files[i]->index->stats->maxcover = maxcover; 
        files[i]->index->stats->minfrac = minfrac;
        files[i]->index->stats->entropyfilter = entropyfilter;
        files[i]->index->stats->usenativequal = usenativequal;
        files[i]->index->stats->usegev = 0;
        files[i]->index->stats->strand = strandbias;
        if(stats && files[i]->index->stats->ecdf) {
            ecdf_destruct(files[i]->index->stats->ecdf);
            FREEMEMORY(space, files[i]->index->stats->ecdf);
        }
        files[i]->index->stats->ecdf = 
          ecdf_init(files[i]->index->stats->e, files[i]->index->stats->e_N);

    }

    //fprintf(stderr, "maxcover: %d\n",files[i]->index->stats->maxcover);

    /*
     * if stats are not present in any file
     * neither gev, dumpstats, nor call
     * can be called
     */
    if (!files[i]->index->stats &&
        (dumpstats || scoreplotbasename || call || starts)){
      manopt_help(&optset, "please generate an index with stats (option -H) if the options -Z, -c, -S are used\n");
    }
  }


  if(manopt_isset(&optset, 'a', "annotation")) { 
    bedfilenames = manopt_getarg(&optset, 'a', "annotation");
    for(i=0; i < bedfilenames->noofvalues; i++) { 
  
      prefixlen = bl_fileprefixlen(bedfilenames->values[i]);
      if( strncmp(&bedfilenames->values[i][prefixlen], ".bed", 3) == 0 ||
        strncmp(&bedfilenames->values[i][prefixlen], ".BED", 3) == 0) {
        track = bl_BEDread(space, bedfilenames->values[i]);

        beddev = fopen("sorted.bed", "w");
        if(beddev  == NULL) {
          fprintf(stderr, "could not open file %s. Exit forced.", "sorted.bed");
          exit(-1);
        }

        bl_BEDwrite(track, beddev);
      
      } else if( strncmp(&bedfilenames->values[i][prefixlen], ".gff", 3) == 0 ||
        strncmp(&bedfilenames->values[i][prefixlen], ".GFF", 3) == 0) {
        track = bl_GFFread(space, bedfilenames->values[i]);
        bl_BEDwrite(track, beddev);
      
      } else {
        manopt_help(&optset, "please provide files with .GFF or .BED extension\n");
      }
    }
  }

  /*
   * TODO: stats for all files, what if stats failed!
   */
  
  if(gatk) {fprintf(stderr, "using GATK model\n"); gatkeval(space, files, groups, queries->noofvalues, fasta, maxframesize, usenativequal); }
  if(consensus){ fprintf(stderr, "calling consensus\n"); evalconsensus(space, files, groups, queries->noofvalues, fasta, maxframesize); }
  //if(GEV)  { fprintf(stderr, "fitting GEV\n"); bl_matchfileFitGEV (NULL, files[0]->index->stats); }
  if(dumpstats) bl_matchfileDumpSampleStats(files[0]->index->stats);
  if(call) { 

    eval(space, files, groups, queries->noofvalues, fasta, maxframesize, minscore);

  }
  if(simplegev) simpleeval(space, files, groups, queries->noofvalues, fasta, maxframesize, usenativequal);
  if(methylcall) callMethylSimple(space, files, groups, queries->noofvalues, fasta, maxframesize);
  if(methylbias) calcMethylBias(space, files, groups, queries->noofvalues, fasta, maxframesize);

  if(manopt_isset(&optset, 'R', "scoreplot")) {
   

    if(!scoreplotbasename) { 
      scoreplotbasename = bl_basename(files[0]->filename);
    }
    NFO("plotting sample scores to basename: %s", scoreplotbasename);     
    
    scoreplotbasenamelen  = strlen(scoreplotbasename);

     char *histofilename = ALLOCMEMORY(space, NULL, char, scoreplotbasenamelen+7+4+1);
     sprintf(histofilename, "%s.histo.dat", scoreplotbasename);

     char *scorefilename = ALLOCMEMORY(space, NULL, char, scoreplotbasenamelen+7+4+1);
     sprintf(scorefilename, "%s.score.dat", scoreplotbasename);
 
     char *cutfilename = ALLOCMEMORY(space, NULL, char, scoreplotbasenamelen+7+4+1);
     sprintf(cutfilename, "%s.cut.dat", scoreplotbasename);
     
     char *splinefilename = ALLOCMEMORY(space, NULL, char, scoreplotbasenamelen+7+4+1);
     sprintf(splinefilename, "%s.spline.dat", scoreplotbasename);
     
     char *estimatefilename = ALLOCMEMORY(space, NULL, char, scoreplotbasenamelen+7+4+1);
     sprintf(estimatefilename, "%s.estim.dat", scoreplotbasename);

     getcutoff(files[0]->index->stats, histofilename, scorefilename, cutfilename, splinefilename, estimatefilename);

     FREEMEMORY(space, histofilename);
     FREEMEMORY(space, scorefilename);
     FREEMEMORY(space, cutfilename);
     FREEMEMORY(space, splinefilename);
     FREEMEMORY(space, estimatefilename);

  }


  if(manopt_isset(&optset, 's', "splice") && 
      manopt_isset(&optset, 'q', "query")) {

      if(!splicebasename) { 
        splicebasename = bl_basename(files[0]->filename);
      }
  
      splicebasenamelen  = strlen(splicebasename);

      filename = ALLOCMEMORY(space, NULL, char, splicebasenamelen+7+4+1);
      sprintf(filename, "%s.splice.bed", splicebasename);

      transfilename = ALLOCMEMORY(space, NULL, char, splicebasenamelen+7+4+1);
      sprintf(transfilename, "%s.trans.bed", splicebasename);

      // not used: splice = 1;
      evalsplice(space, files, groups, queries->noofvalues, fasta, maxframesize, track, filename, transfilename, minsplitno);

      FREEMEMORY(space, filename);
      FREEMEMORY(space, transfilename);
  }

  
  if(browse) view(space, files, queries->noofvalues, fasta, track);
  if(starts) evalstarts(space, files, groups, queries->noofvalues, fasta, maxframesize, track);
  if(wiggle) writeexpression(space, files, groups, queries->noofvalues, fasta, track);

  bl_fastaDestruct(space, fasta);
  FREEMEMORY(space, fasta);
  
  if(files) {
    for(i=0; i < queries->noofvalues; i++) { 
      if (files[i]->index) {
        bl_matchfileDestructIndex(space, files[i]->index);
        FREEMEMORY(space, files[i]->index);
      }
      FREEMEMORY(space, files[i]);
    }
    FREEMEMORY(space, files);
  }

  if(groups) {
    FREEMEMORY(space, groups);
  }


  if(track) { 
    bl_annotationtrackDestruct(space, track);
    FREEMEMORY(space, track);
  }


  MSG("Goodbye.\n");

  manopt_destructoptionset(&optset);
  manopt_destructarg(unflagged);
  FREEMEMORY(space, unflagged);

  FREEMEMORY(space, ntcode);

  return 0;
}

