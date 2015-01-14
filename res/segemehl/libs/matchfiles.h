#ifndef MATCHFILES_H
#define MATCHFILES_H

/*
 *
 *	readmatchfiles.h
 *  
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 06.10.2010 01:10:47 CEST  
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <sys/types.h>
#include "zran.h"
#include "biofiles.h"
#include "basic-types.h"
#include "mathematics.h"

#define SAM 0
#define SEG 1
#define MAXREADLENGTH 1500
#define MINQUAL 33
#define MAXQUAL 104
#define QRNGE (MAXQUAL-MINQUAL+1)

#define GZ_IDX_STORED ((unsigned char) (1 << 0))
#define STATS_TAB_STORED ((unsigned char) (1 << 1))
#define IDXMD5_STORED ((unsigned char ) (1 << 2))
#define STATS_GEV_STORED ((unsigned char ) (1 << 3))
#define MOTIF_STORED ((unsigned char ) (1 << 4))
#define STRAND_STORED ((unsigned char ) (1 << 5))
#define PXX_STORED ((unsigned char ) (1 << 6))

#define MFREAD_FEAT  ((unsigned char) (1 << 0))
#define MFREAD_QUAL ((unsigned char) (1 << 1))
#define MFREAD_DELS ((unsigned char) (1 << 2))
#define MFREAD_SPLITS ((unsigned char) (1 << 3))
#define MFREAD_RLEN ((unsigned char) (1 << 4))
#define MFREAD_MCNT ((unsigned char) (1 << 5))
#define MFREAD_ROW ((unsigned char) (1 << 6))
#define MFREAD_RPOS ((unsigned char) (1 << 7))


typedef struct {
  Uint len;
  char *string;
  char *quals;
  Uint row;
  Uint edist;
  Uint readpos;
  Uint curstart;
  Uint matchcnt;
  Uint bisulfite;
} matchfileDeletion_t;


typedef struct {
  Uint spliceid;
  char edgetype;                                        
  char strand;  /* strand of the exonedge                                   */ 
  char trans;
  Uint xno;
  Uint xstart;  /* start of the split in the read                           */
  Uint xend;    /* end of the split in the read                             */
  Uint edgechridx;
  Uint edge;
  Uint adjoint;
} matchfileSplit_t;

typedef struct {
  Uint refidx;
  Uint refpos;
  Uint noofmates;
} matelink_t;

typedef struct {
  Uint len;
  char ref; 
  char cons;
  char *chars;
  char *quals;
  char *strands;
  char *feat;
  uint32_t *readpos;
  uint32_t *readlen;
  uint32_t *row;
  uint32_t *matchcnt; 
  uint32_t starts;
  uint32_t ends;
  Uint noofmatelinks;
  matelink_t *matelinks;
  Uint maxrow; 
  uint32_t noofdels; 
  matchfileDeletion_t *dels;
  uint32_t noofsplits;
  matchfileSplit_t *splits;
  unsigned char *edist;
  uint32_t *bisulfite;
  double s_ref;
  double s_refx;
  double p_ref;
  double p_refx;
  double s_cons;
  double s_consx;
  double p_cons;
  double p_consx;
  double p_hom;
  double entropy;
  double longentropy;
  double pentropy;
  double ee;
  double scr_ref;
  double scr_cons;
  double scr_sample;
  double diff_rt;
  double diff_rq;
  double diff_rr;
  double diff_mm;
  double pee;
  double pbinom;
  Uint secondminimum;
  Uint secondcnt;
} matchfileCross_t;

typedef struct {
  Uint start;
  Uint width;
  Uint maxheight;
  Uint chrlen;
  char *chrname;
  char *chrseq;
  char *ref;
  matchfileCross_t *cs;
} matchfileFrame_t;


typedef struct { 
  Uint n;
  double px;
  double pxx;
  Uint maxcover;
  Uint mincover;
  double minfrac;
  char entropyfilter;
  /*errordensity*/
  Uint e_N;
  double *eraw;
  double *entropy;
  double *b;
  double b_mu;
  double b_sd;
  double b_ll;
  double *e;
  double *e_mu;
  double *e_sd;
  double e_ll;
  double *gev_mu;
  double *gev_si;
  double *gev_xi;
  double *gev_ll;
  Uint P;
  Uint X;
  Uint N;
  /*readerror at non variant positions: 0,1,2,3,4,5,>5*/
  Uint RR_N;
  Uint *RR;
  /*multiple mapping at non variant positions: 0,1,2,3,4,5,>5*/
  Uint MM_N;
  Uint *MM;
  /*norm area*/
  Uint areasize;
  double maxareae;
  /*substitution*/
  Uint *S_N;
  double *S;
  Uint *Sx_N; 
  double *Sx;
  /*noise*/
  Uint *R_N;
  Uint *R;
  Uint *RP_N;
  Uint *RP;
  Uint *RQ_N;
  Uint *RQ;
  Uint *MO_N;
  Uint *MO;
  /*readvariance*/
  Uint V_N;
  double *V;
  double V_mu;
  double V_sd;
  double V_ll;
  Uint Vx_N;
  double *Vx;
  double Vx_mu;
  double Vx_sd;
  double Vx_ll; 
  double entropydensitystep;
  Uint entropydensitylen;
  double *entropydensity;
  char usenativequal;
  char usegev;
  /*standardization*/
  char standardized;
  double maxrp;
  double minrp;
  double maxrq;
  double minrq;
  double currq;
  ecdf_t *ecdf;
  double *s;
  Uint s_N;
  double cut;
  char strand;
  unsigned char minqual;
  unsigned char maxqual;
} matchfileSampleStats_t;

typedef struct {
  Uint start;
  Uint end;
  Uint matches;
  off_t offset;
  off_t endoff;
} matchfileBin_t;

typedef struct {
  struct access *gzindex;
  Uint md5;
  Uint noofchroms;
  Uint noofreads;
  Uint exp; //power of two
  char **chromnames;
  Uint *matchstart;
  Uint *matchend;
  Uint *noofbins;
  matchfileBin_t **bins;
  Uint maxreadlen;
  double *submatrix;
  double mean_coverage;
  double mean_qual;
  Uint *P_ERR;
  Uint *Q_ERR;
  Uint *Q_N;  
  matchfileSampleStats_t *stats;
  unsigned char minqual;
  unsigned char maxqual;
} matchfileindex_t;

typedef struct {
  unsigned char fmt;
  unsigned char gzip;
  char *filename;
  matchfileindex_t *index;
} matchfile_t;

typedef struct {
  double *mean_err;
  double *mean_sde;
  double *mean_pos;
  double *mean_mul;
  double *mean_dis;
  Uint rss;
  Uint *char_frq;
  Uint *dist_rss;
  Uint dist_rss_N;
  Uint dist_rss_ymax;
  Uint dist_rss_ymax_1;
  Uint dist_rss_xmax;
  double mean_cov;
  double prime5;
  double prime3;
  Uint *ntcnt;
} matchfileFrameStats_t;

typedef struct {
  Uint len;
  double var_ee;
  double* var_s;
  double* var_rt;
  double* var_rq;
  double* var_rr;
  double* var_rv;
  double* var_mm;
  double* sub;
  double pentropy;
  double strandpenalty;
  double mean_rt;
  double mean_rq;
  double mean_rr;
  double mean_mm;

} matchfileCrossStats_t;

typedef struct{
  Uint noofsplits;
  Uint distpos;
  Uint distcidx;
  matchfileCross_t *cs;
  Uint pos;
  Uint cidx;
  Uint splicesite;
  Uint acceptor;
  Uint donor;
  Uint transsplits;
  char seen;
} distsplitsites_t;

typedef struct{ 
  Uint binpos;
  Uint binref;
  Uint noofmates;
  Uint *matebinpos;
  Uint *matebinref;
  Uint *matebincnt;
} matebin_t;

typedef struct{
  Uint noofbins;
  matebin_t *bins;
} matebinmap_t;

typedef struct {
  Uint noofmates;
  Uint *refpos;
  Uint *refidx;
  Uint *materefpos;
  Uint *materefidx;
  Uint *matecount;
} matemap_t;

typedef struct {
  annotationtrack_t *bed;
  matchfile_t **files;
  fasta_t *set;
  Uint noofsplits;
  Uint *pos;
  Uint *cidx;
  matchfileCross_t *cs;
  matebinmap_t matemap;
} splitmap_t;

typedef struct{  
  char *chromname;
  Uint cidx;            /*the idx of the chromsome*/
  Uint start;           /*the start of the interval*/
  Uint end;             /*the end of the interval*/        
  Uint median;          /*median split site = splice site*/

  /* *
   *
   * information on splitsites and their splits
   * 
   * */

  Uint noofsplitsites;  /*no of splitsites*/
  Uint *splitsites;     /*splitsites contributing to this splice site*/
  Uint *noofsplits;     /*no of splits for each splitsite*/

                        /*the respective cross sections*/
  matchfileCross_t **cs;

  Uint totalsplits;     /*total no of splits*/
  Uint dtypes;          /*no of donor splits*/
  Uint atypes;          /*no of accep splits*/
  Uint mstrands;        /*no of minus splits*/
  Uint pstrands;        /*no of plus splits*/
  Uint transsplits;     /*no of trans splits*/

  char type;            /*consensus type (A or D)*/
  char strand;          /*consensus strand (- or +)*/
 
  Uint noofrightsites;
                        /*pointers into this very list ( )*/
  Uint *rightsiteidx;
  Uint *rightedgeweight;
  uint16_t *rightedgeacceptor;
  uint16_t *rightedgedonor;
  uint16_t *righttranssplits;
  uint16_t *leftmatesupport;
  uint16_t *rightmatesupport;

  Uint noofleftsites;
                        /*pointers into this very list ( )*/
  Uint *leftsiteidx;
  Uint *leftedgeweight;
  uint16_t *leftedgeacceptor;
  uint16_t *leftedgedonor;
  uint16_t *lefttranssplits;

} splicesite_t;


typedef struct {
  
  Uint noofsplicesites;
  splicesite_t *map;
 
  //stats
  Uint interval;
  Uint *histogram;
  Uint **charhist;
  Uint **charhistA;
  Uint *chrcnt;

} splicemap_t;

typedef struct{
  char *curname;
  char *curchrom;
  Uint flag;
  char *flgs;
  Uint curstart;
  Uint curend;
  Uint curchromidx;
  char *curseq;
  char *curqual;
  char *diff;
  Uint curcnt;
  char *curaln;
  Uint curalnlen;
  int edist;
  Uint bisulfite;
  char strand;
  char *rnext;
  Uint pnext;
  Uint acceptorpos;
  char *acceptorchr;
  Uint acceptorflg;
  Uint donorpos;
  char *donorchr;
  Uint donorflg;
  Uint xstart;
  Uint xend;
  Uint xno;
  Uint noofsplits;
  Uint identity;
} matchfileRec_t;


Uint
bl_matchfileGetChromIndexNumberDB (fasta_t *set, matchfileindex_t *index, 
    char *chr);
matchfileCross_t* bl_matchfileRead(void *space, matchfile_t *file, 
    char *chromname, Uint start, Uint end, Uint maxcover, fasta_t *set, unsigned char fields, matchfileCross_t*);
void bl_matchfileIndex(void *space, matchfile_t *file, fasta_t *set);
void bl_matchfileDestructCross(void *space, matchfileCross_t *cs, Uint len);
Uint bl_matchfileGetChromIndexNumber(matchfileindex_t *index, char *chromname);
void bl_matchfileDestructIndex(void *space, matchfileindex_t *index);
matchfileCross_t* bl_matchfileCopyCross(void *space, matchfileCross_t *xs, matchfileCross_t *cs);
void bl_matchfileGapAlign(matchfileDeletion_t *dels, Uint noofdels);
void bl_matchfileDumpFileStats (void *space, matchfile_t *file);
void bl_matchfileDumpCrossSection (matchfileCross_t *cs, Uint len);
matchfileindex_t * bl_matchfileReadIndex (void *space, char *filename);
void bl_matchfileWriteIndex(matchfileindex_t *idx, char *filename); 
Uint bl_compareIndices(matchfileindex_t *i1, matchfileindex_t *i2);
matchfileindex_t* bl_matchfileInitIndex(void *space); 
Uint bl_matchfileAdjustBounds(void *space, matchfile_t *file, fasta_t *set,
    Uint setid, Uint matid, Uint start, Uint width, Uint *newstart, Uint *newwidth);
matchfileFrame_t * bl_matchfileGetFrame(void *space, matchfile_t *file, 
    char *chrname, Uint start, Uint width, fasta_t *set, Uint maxcover, matchfileCross_t*);
Uint bl_writeexpression (void *space, char *filename, Uint filenamelen,
    matchfile_t *file, fasta_t *set, Uint width, Uint maxcover);
Uint bl_matchfileIndexAddChrom(matchfileindex_t *index, char *chromname);
void bl_matchfileMergeFrames (void *space, matchfileFrame_t* f, matchfileFrame_t *g);


#endif
