
/*
 *  evalmatchfileshelper.c
 *  helper functions 
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 12/24/2013 02:48:45 AM CET
 *  
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include "sort.h"
#include "alignment.h"
#include "stringutils.h"
#include "basic-types.h"
#include "mathematics.h"
#include "matfile.h"
#include "bitVector.h"
#include "info.h"
#include "vtprogressbar.h"
#include "fileio.h"
#include "matchfilesfields.h"
#include "matchfiles.h"
#include "debug.h"
#include "evalmatchfiles.h"
#include "biofiles.h"
#include "splicesites.h"


/*-------------------------- bl_matchfileSimpleGEV ---------------------------
 *    
 * @brief simple GEV calling
 * @author Steve Hoffmann 
 *   
 */

Uint
bl_matchfileSimpleGEV (void *space, Uint fidx, Uint cidx, Uint pos, matchfileCross_t *cs, 
    char ref, matchfileindex_t *idx, unsigned char show, void *nfo )
{

  double p = .0, px =.0, e, pxP =.0, pP = .0;
  char *ch, *rq, *usequal;
  matchfileSampleStats_t *stats = idx->stats;
  Uint len, i, errors=0;

  usequal = (char*) nfo;

  len = cs->len;
  if(len < 2) return 0;

  ch = cs->chars;
  rq = cs->quals;
  
  
  

  for(i=0; i < len; i++) { 

    if((int)ch[i] != (int)ref) {
      errors++;
      //error probability P
      //if P is high -> bad for px, good for p        
      //if P is low  -> bad for p, good for px

      if(*usequal) { 
        pP  += log(pow(10.0,((double)((double)rq[i]-64.0)/-10.0)))+log(0.9999);
        pxP += log(1.0-(pow(10.0,((double)((double)rq[i]-64.0)/-10.0))))+log(0.0001);
      }
    } 
  }   

  if(errors) { 
    e = (double)errors/(double)len;
    p  = log(0.999)+log(1-gevcdf(e,stats->gev_mu[0],stats->gev_si[0],stats->gev_xi[0]))+pP;
    px = log(0.001)+log(gevcdf(e,stats->gev_mu[0],stats->gev_si[0],stats->gev_xi[0]))+pxP;

    if(px > p) fprintf(stdout, "chr21\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%s\n", pos, len, e, p, px, pP, pxP, cs->chars);
    fflush(stdout);
  }

  return 0;
}



/*-------------------------- bl_matchfileSimpleGATK --------------------------
 *    
 * @brief a simple GATK implementation w/o fragment stuff
 * @author Steve Hoffmann 
 * @info modeled after https://github.com/broadgsa/gatk/blob/master/public/java/src/org/broadinstitute/sting/gatk/walkers/genotyper/DiploidSNPGenotypeLikelihoods.java
 *
 *
 * Suppose we have bases b1, b2, ..., bN with qualities scores q1, q2, ..., qN.  This object
 * calculates:
 *
 * P(G | D) = P(G) * P(D | G)
 *
 * where
 *
 * P(D | G) = sum_i log10 P(bi | G)
 *
 * and
 *
 * P(bi | G) = 1 - P(error | q1) if bi is in G
 *           = P(error | q1) / 3 if bi is not in G
 *
 * for homozygous genotypes and for heterozygous genotypes:
 *
 * P(bi | G) = 1 - P(error | q1) / 2 + P(error | q1) / 6 if bi is in G
 *           = P(error | q1) / 3 if bi is not in G
 *
 * for each of the 10 unique diploid genotypes AA, AC, AG, .., TT
 *
 * Everything is stored as arrays indexed by DiploidGenotype.ordinal() values in log10 space.
 */

Uint
bl_matchfileSimpleGATK ( void *space, Uint fidx, Uint cidx, Uint pos, matchfileCross_t *cs, 
    char ref, matchfileindex_t *idx, unsigned char show, void *nfo )
{

  char *ch, *rq, n1, n2;
  double *gtslik = NULL, bl;
  double HETLIK = log10(10e-3);
  char *gts[] = {"AA", "AC", "AG", "AT", "CC", "CG", "CT", "GG", "GT", "TT"};
  char *obs;
  Uint len, i, j, errors=0, bj=0;

  len = cs->len;
  if(len < 2) return 0;

  ch = cs->chars;
  rq = cs->quals;

  gtslik = ALLOCMEMORY(space, NULL, double, 10);
  memset(gtslik, 0, sizeof(double)*10);
  gtslik[1] = HETLIK;
  gtslik[2] = HETLIK;
  gtslik[3] = HETLIK;
  gtslik[5] = HETLIK;
  gtslik[6] = HETLIK;
  gtslik[8] = HETLIK;

  obs = ALLOCMEMORY(space, NULL, char, 256);
  memset(obs, 0, sizeof(char)*256);
  

  //AA, AC  AG, AT, CC, CG, CT, GG, GT, TT
  // 0,  1,  2,  3,  4,  5,  6,  7,  8,  9

  for(j=0; j < 10; j++) {
    n1 = gts[j][0];
    n2 = gts[j][1];

    for(i=0; i < len; i++) { 
      obs[(int)ch[i]] = 1;
      if(n1 == n2) { 
        if((int)ch[i] != (int)n1) {
          errors++;
          //error probability P
          //if P is high -> bad for px, good for p        
          //if P is low  -> bad for p, good for px
          gtslik[j] += log10(
                   pow(10.0,((double)((double)rq[i]-64.0)/-10.0)))-log10(3.0);
        } else {
          gtslik[j] += log10(
              1.0-(pow(10.0,((double)((double)rq[i]-64.0)/-10.0))));
        }

      } else {
        if((int)ch[i] != (int)n1 && (int)ch[i] != (int)n2) {
          errors++;
          //error probability P
          //if P is high -> bad for px, good for p        
          //if P is low  -> bad for p, good for px
          gtslik[j] += log10(pow(10.0,((double)((double)rq[i]-64.0)/-10.0)))-log10(3.0);
        } else {
          gtslik[j] += log10(
              1.0-(pow(10.0,((double)((double)rq[i]-64.0)/-10.0))/2) + 
                  (pow(10.0,((double)((double)rq[i]-64.0)/-10.0))/6));
        }
      }
    }   
  }

  switch(ref) {
    case 'A':
    case 'a':
      bj = 0;
      break;
    case 'C':
    case 'c':
      bj = 4;
      break;
    case 'G':
    case 'g':
      bj = 7;
      break;
    case 'T':
    case 't':
      bj = 9;
      break;
    default:
      bj = 0;
  }

  bl = gtslik[bj];

//  fprintf(stdout, "--- \n");
  for(j=0; j < 10; j++) {
    if(pos == 9413199)   fprintf(stdout, "GT:%s\t%f\n", gts[j], gtslik[j]);
    if(gtslik[j] > bl && gtslik[j] != 0.0 && obs[(int)gts[j][0]] && obs[(int)gts[j][1]]) {
      bl = gtslik[j];
      bj = j;
    }
  }

  n1 = gts[bj][0];
  n2 = gts[bj][1];

  if(ref != (int) 'N' && (n1 != (int)ref || n2 != (int)ref)) {
    fprintf(stdout, "CS[%d]:%s\n", len, ch);
    fprintf(stdout, "%d\t%d\t%c\t%c\t%c\n", cidx, pos, ref, n1, n2);
    if(pos == 9413199) exit(-1);
   // exit(-1);
  }

  FREEMEMORY(space, gtslik);
  return 0;
}

/*----------------- bl_matchfileGetCrossConsErrorQualScaled -----------------
 *    
 * @brief get the error e of a cross section based on the consensus 
 * (not reference)
 * @author Steve Hoffmann 
 *   
 */
 
double
bl_matchfileGetCrossConsErrorQualScaled (matchfileFrame_t *frame, Uint pos)
{
  Uint j;
  double e=.0, mat = .0, mis = .0;
 
  for(j=0; j < frame->cs[pos].len; j++) {
    if(frame->cs[pos].chars[j] != frame->cs[pos].cons) { 
       mis += pow(10,((double)((double)frame->cs[pos].quals[j]-64.0)/-10.0)); 
    } else {
       mat += 1.0-pow(10,((double)((double)frame->cs[pos].quals[j]-64.0)/-10.0));        
    }
  }
   
  if(frame->cs[pos].len && mat) { 
    e = mis/mat;
  }

  return e;
}



/*------------------------ bl_matchfileGetStrandBias -------------------------
 *    
 * @brief get the strand bias of a cross section
 * @author Steve Hoffmann 
 *   
 */
 
double
bl_matchfileGetStrandBias (matchfileFrame_t *frame, Uint pos)
{ 
  
  Uint j;
  double e=0;
 
  for(j=0; j < frame->cs[pos].len; j++) {
    if(frame->cs[pos].strands[j] == '+')
      e++;
  }

  if(frame->cs[pos].len) { 
    e /= frame->cs[pos].len;
  }

  return e;
}
/*----------------------- bl_matchfileGetNTRedundancy ------------------------
 *    
 * @brief get for each NT average redundancy of reads (multiple read hit count) 
 * in the cross section
 * @author Steve Hoffmann 
 *   
 */

double*
bl_matchfileGetNTRedundancy(matchfileFrame_t *frame, Uint pos) {
  Uint j, *c;
  double *r;

  r = ALLOCMEMORY(space, NULL, double, 256);
  c = ALLOCMEMORY(space, NULL, Uint, 256);
  memset(r, 0, sizeof(double)*256);
  memset(c, 0, sizeof(Uint)*256);

  for(j=0; j < frame->cs[pos].len; j++) {
    r[(int)frame->cs[pos].chars[j]] 
      += (int)frame->cs[pos].matchcnt[j];
    c[(int)frame->cs[pos].chars[j]]++;
  }

  for(j=0; j < 256; j++) {
    if(c[j]) r[j] /= c[j];
  }

  FREEMEMORY(space, c);
  return r;
}



/*-------------------------- bl_matchfileGetNTError --------------------------
 *    
 * @brief return average quality score, ie error probability,
 * in a cross section for each nucleotide present in the cross section
 * @author Steve Hoffmann 
 *   
 */

double*
bl_matchfileGetNTError(matchfileFrame_t *frame, Uint pos) {
  Uint j;
  Uint *c;
  double *p;


  p = ALLOCMEMORY(space, NULL, double, 256);
  c = ALLOCMEMORY(space, NULL, Uint, 256);
  memset(c, 0, sizeof(Uint)*256);

  for(j=0; j < 256; j++) {
    p[j] = log10(0);
  }

  for(j=0; j < frame->cs[pos].len; j++) {
    p[(int)frame->cs[pos].chars[j]] = 
      log10add(p[(int)frame->cs[pos].chars[j]],
          (((double)frame->cs[pos].quals[j])/-10.)); 
    c[(int)frame->cs[pos].chars[j]]++;
  }

  for(j=0; j < 256; j++) {
    if(c[j]) p[j] -= log10(c[j]);
  }

  FREEMEMORY(space, c);
  return p;
}

/*--------------------------- bl_matchfileGetNTQual ----------------------------
 *    
 * @brief get for each NT average qual in a cross section
 * @author Steve Hoffmann 
 *   
 */
 

double*
bl_matchfileGetNTQual(matchfileFrame_t *frame, Uint pos) {
  Uint j, *c;
  double *r;

  r = ALLOCMEMORY(space, NULL, double, 256);
  c = ALLOCMEMORY(space, NULL, Uint, 256);
  memset(r, 0, sizeof(Uint)*256);
  memset(c, 0, sizeof(Uint)*256);

  for(j=0; j < frame->cs[pos].len; j++) {
    r[(int)frame->cs[pos].chars[j]] 
      += frame->cs[pos].quals[j];
    c[(int)frame->cs[pos].chars[j]]++;
  }

  for(j=0; j < 256; j++) {
    if(c[j]) r[j] /= c[j];
  }

  FREEMEMORY(space, c);
  return r;
}


/*------------------------- bl_matchfileGetNTReadPos -------------------------
 *    
 * @brief get for each NT average read position in a cross section
 * @author Steve Hoffmann 
 *   
 */
 

double*
bl_matchfileGetNTReadPos(matchfileFrame_t *frame, Uint pos) {
  Uint j, *c;
  double *r;

  r = ALLOCMEMORY(space, NULL, double, 256);
  c = ALLOCMEMORY(space, NULL, Uint, 256);
  memset(r, 0, sizeof(Uint)*256);
  memset(c, 0, sizeof(Uint)*256);

  for(j=0; j < frame->cs[pos].len; j++) {
    r[(int)frame->cs[pos].chars[j]] 
      += frame->cs[pos].readpos[j];
    c[(int)frame->cs[pos].chars[j]]++;
  }

  for(j=0; j < 256; j++) {
    if(c[j]) r[j] /= c[j];
  }

  FREEMEMORY(space, c);
  return r;
}


/*----------------------- bl_matchfileGetNTReadPosVar ------------------------
 *    
 * @brief get for each NT read position variance in a cross section pos
 * @author Steve Hoffmann 
 *   
 */
 

double*
bl_matchfileGetNTReadPosVar(matchfileCross_t *cs) {
  double *v, **r;
  int *c, j;

  r = ALLOCMEMORY(space, NULL, double*, 256);
  c = ALLOCMEMORY(space, NULL, Uint, 256);
  v = ALLOCMEMORY(space, NULL, double, 256);
  memset(r, 0, sizeof(double*)*256);
  memset(c, 0, sizeof(Uint)*256);
  memset(v, 0, sizeof(double)*256);

  for(j=0; j < cs->len; j++) {

    r[(int)cs->chars[j]] =
      ALLOCMEMORY(space, 
          r[(int)cs->chars[j]], double, 
          c[(int)cs->chars[j]]+1);
 
    r[(int)cs->chars[j]][c[(int)cs->chars[j]]]
      = (double)cs->readpos[j] / (double)cs->readlen[j];  

    c[(int)cs->chars[j]]++;
  }

  for(j=0; j < 256; j++) {
    if(c[j]) { 
      v[j] = var(r[j],c[j]);
    }
    FREEMEMORY(space, r[j]);
  }

  FREEMEMORY(space, r);
  FREEMEMORY(space, c);
  return v;
}




/*-------------------------- bl_matchfileGetNTEdist --------------------------
 *    
 * @brief for each cross section in the frame:
 * Get the average edist of reads in the cross section
 * @author Steve Hoffmann 
 *   
 */
 
double*
bl_matchfileGetNTEdist(matchfileFrame_t *frame, Uint pos) {
  Uint j, *c;
  double *e;

  e = ALLOCMEMORY(space, NULL, double, 256);
  c = ALLOCMEMORY(space, NULL, Uint, 256);
  memset(e, 0, sizeof(Uint)*256);
  memset(c, 0, sizeof(Uint)*256);


  for(j=0; j < frame->cs[pos].len; j++) {
    e[(int)frame->cs[pos].chars[j]] 
      += (int)frame->cs[pos].edist[j];
    c[(int)frame->cs[pos].chars[j]]++;
  }

  for(j=0; j < 256; j++) {
    if(c[j]) e[j] /= c[j];
  }

  FREEMEMORY(space, c);
  return e;
}



/*---------------------- bl_matchfileGetRegionMotif ----------------------
 *    
 * @brief get the error e of a of a region of cross section based 
 * (not reference)
 * @author Steve Hoffmann 
 *   
 */

Uint
bl_matchfileGetRegionMotif (matchfileFrame_t *frame, 
    Uint pos, Uint left, Uint right)
{
  Uint i, l, r, motif=0;
  

  assert(frame->width > right);

  if(left >= pos) {
    return -1;
  } else { 
    if(pos + right >= frame->width-1) {
      return -1;
    } else { 
      l = pos - left;
      r = pos + right;
    }
  }


  for(i=r; i >= l; i--) { 
    switch(frame->cs[i-1].cons) { 
    case 'A':
      // motif += 0* pow(4,(r-i))
      break;
    case 'C':
      motif += 1*pow(4,(r-i));
      break;
    case 'G':
      motif += 2*pow(4,(r-i));
      break;
    case 'T':
      motif += 3*pow(4,(r-i));
      break;
    default:
      return -1;
      break;
    }
  }
  
  return motif;
}

/*---------------------- bl_matchfileGetRegionConsError ----------------------
 *    
 * @brief get the error e of a of a region of cross section based 
 * (not reference)
 * @author Steve Hoffmann 
 *   
 */

  double
bl_matchfileGetRegionConsError (matchfileFrame_t *frame, 
    Uint pos, Uint range)
{
  Uint i, j, l, r;
  double e=0, ef=0;

  assert(frame->width > range);

  if(range/2 >= pos) {
    l = 0;
    r = range;
  } else { 
    if(range/2 + pos >= frame->width-1) {
      l = frame->width-1-range;
      r = frame->width-1;
    } else { 
      l = pos - range/2;
      r = pos + range/2-1;
    }
  }


  for(i=l; i < r; i++) { 

    for(j=0, e=0; j < frame->cs[i].len; j++) {
      if(frame->cs[i].chars[j] != frame->cs[i].cons)
        e++;
    }

    if(frame->cs[i].len) { 
      e /= (double)frame->cs[i].len;
      ef += e;
    }
  }
  
  return ef/(double) range;
}


/*---------------------------- bl_matchfileFitGEV ----------------------------
 *    
 * @brief fit a generalized extreme value distribution to raw and adjust errors
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchfileFitGEV (void *space, matchfileSampleStats_t *stats) 
{

  double m=.0, s=.0, k=.0;
  double mx=.0, sx=.0, kx=.0;

  //NFO("evaluating sample stats with %d values (raw error).\n", stats->e_N); 
  qsort(stats->eraw, stats->e_N, sizeof(double), cmp_dbl_qsort);
  
  gevLmoment(stats->eraw, stats->e_N, &m, &s, &k);
  NFO("lmoment m:%f, s:%f, xi:%f \n", m, s, -k);
  NFO("log-likelihood:%f\n", gevll( stats->eraw, stats->e_N, m, s, k)); 
  gevmle(NULL, stats->eraw, stats->e_N, &m, &s, &k, 10000, stats->eraw[0], 
      stats->eraw[stats->e_N-1]);
  NFO("gev m:%f, s:%f, xi:%f \n", m, s, -k);
  NFO("log-likelihood:%f\n", gevll( stats->eraw, stats->e_N, m, s, k));
  NFO("cdf example for mu%f =%f\n", mx, gevcdf(m, m, s, -k ));
  NFO("cdf example for 0.1 %f\n", gevcdf(0.1, m, s, -k ));
  NFO("cdf example for 0.2 %f\n", gevcdf(0.2, m, s, -k ));
  NFO("cdf example for 0.3 %f\n", gevcdf(0.3, m, s, -k ));
  NFO("cdf example for 0.4 %f\n", gevcdf(0.4, m, s, -k ));
  NFO("cdf example for 0.5 %f\n", gevcdf(0.5, m, s, -k ));
  NFO("cdf example for 0.6 %f\n", gevcdf(0.6, m, s, -k ));
 
  stats->gev_mu[0] = m;
  stats->gev_si[0] = s;
  stats->gev_xi[0] = -k;
  stats->gev_ll[0] = gevll(stats->eraw, stats->e_N, mx, sx, kx);

  NFO("evaluating sample stats with %d values (adj error).\n", stats->e_N); 
  qsort(stats->e, stats->e_N, sizeof(double), cmp_dbl_qsort);
  gevLmoment(stats->e, stats->e_N, &mx, &sx, &kx);

  NFO("lmoment m:%f, s:%f, xi:%f \n", mx, sx, -kx);
  NFO("log-likelihood:%f\n", gevll(stats->e, stats->e_N, mx, sx, kx));
  
  gevmle(NULL, stats->e, stats->e_N, &mx, &sx, &kx, 10000, stats->e[0], 
      stats->e[stats->e_N-1]);
  NFO("gev m:%f, s:%f, xi:%f \n", mx, sx, -kx);
  NFO("log-likelihood:%f\n", gevll( stats->e, stats->e_N, mx, sx, kx));
  NFO("variance:%f (ll:%f)\n",gevvar(mx,sx,-kx), gevcdf(gevvar(mx,sx,-kx)+mx, 
        mx, sx, -kx));
  NFO("cdf example for mu%f =%f\n", mx, gevcdf(mx, mx, sx, -kx ));
  NFO("cdf example for 0.1 %f\n", gevcdf(0.1, mx, sx, -kx ));
  NFO("cdf example for 0.2 %f\n", gevcdf(0.2, mx, sx, -kx ));
  NFO("cdf example for 0.3 %f\n", gevcdf(0.3, mx, sx, -kx ));
  NFO("cdf example for 0.4 %f\n", gevcdf(0.4, mx, sx, -kx ));
  NFO("cdf example for 0.5 %f\n", gevcdf(0.5, mx, sx, -kx ));
  NFO("cdf example for 0.6 %f\n", gevcdf(0.6, mx, sx, -kx ));
  
  stats->gev_mu[1] = mx;
  stats->gev_si[1] = sx;
  stats->gev_xi[1] = -kx;
  stats->gev_ll[1] = gevll(stats->e, stats->e_N, mx, sx, kx);


  return ;
}


/*----------------------- bl_matchfileInitSampleStats ------------------------
 *    
 * @brief initialize sample stats
 * @author Steve Hoffmann 
 *   
 */
 
matchfileSampleStats_t*
bl_matchfileInitSampleStats (void *space, Uint maxsample, Uint mincover, 
    Uint maxcover, double minfrac, char entropyfilter, Uint areasize, double maxareae)
{

  matchfileSampleStats_t *stats;

  stats = ALLOCMEMORY(space, NULL, matchfileSampleStats_t, 1);
  stats->standardized = 0;
  stats->n = maxsample;
  stats->px = .0;
  stats->pxx = .0;
  stats->b_ll = .0;
  stats->e = ALLOCMEMORY(space, NULL, double, maxsample);
  stats->b = ALLOCMEMORY(space, NULL, double, maxsample);
  stats->s= ALLOCMEMORY(space, NULL, double, maxsample);
  stats->eraw = ALLOCMEMORY(space, NULL, double, maxsample);
  stats->entropy = ALLOCMEMORY(space, NULL, double, maxsample);

  stats->entropydensity = NULL;
  stats->entropydensitystep = 0;
  stats->entropydensitylen = 0;
  stats->mincover = mincover;
  stats->maxcover = maxcover;
  stats->areasize = areasize;
  stats->maxareae = maxareae;
  stats->minfrac = minfrac;
  stats->entropyfilter = entropyfilter;
  stats->e_N = 0;
  stats->s_N= 0;
  stats->px = .0;
  stats->V_N = 0;
  stats->V_mu = .1;
  stats->V_sd = .1;
  stats->V_ll = 0;
  stats->Vx_N = 0;
  stats->Vx_mu = .1;
  stats->Vx_sd = .1;
  stats->Vx_ll = 0;
  stats->P = 0;
  stats->X = 0;
  stats->N = 0;
    
  stats->RR_N = 0;
  stats->RR = ALLOCMEMORY(space, NULL, Uint, 11);
  memset(stats->RR, 0, sizeof(Uint)*11);

  stats->MM_N = 0;
  stats->MM = ALLOCMEMORY(space, NULL, Uint, 51);
  memset(stats->MM, 0, sizeof(Uint)*51);
  
  stats->e_mu = ALLOCMEMORY(space, NULL, double, 2);
  stats->e_sd = ALLOCMEMORY(space, NULL, double, 2);
  stats->gev_mu = ALLOCMEMORY(space, NULL, double, 2);
  stats->gev_si = ALLOCMEMORY(space, NULL, double, 2); 
  stats->gev_xi = ALLOCMEMORY(space, NULL, double, 2);
  stats->gev_ll = ALLOCMEMORY(space, NULL, double, 2);

  stats->e_mu[0] = 0.1;
  stats->e_mu[1] = 0.6;
  stats->e_sd[0] = 0.1;
  stats->e_sd[1] = 0.6;
  stats->e_ll =.0;

  stats->gev_mu[0] = 0.044763;
  stats->gev_mu[1] = 0.020171;
  stats->gev_si[0] = 0.022864;
  stats->gev_si[1] = 0.031077;
  stats->gev_xi[0] = 0.212219;
  stats->gev_xi[1] = -0.041355; 
  stats->gev_ll[0] = 6291.208397;
  stats->gev_ll[1] = 5908.074411;


  stats->S = ALLOCMEMORY(space, NULL, double, 6*6);
  stats->S_N = ALLOCMEMORY(space, NULL, Uint, 6);
  stats->Sx = ALLOCMEMORY(space, NULL, double, 6*6);
  stats->Sx_N = ALLOCMEMORY(space, NULL, double, 6);
  stats->R = ALLOCMEMORY(space, NULL, Uint, 100*255);
  stats->R_N = ALLOCMEMORY(space, NULL, Uint, 100*255);
  stats->V = ALLOCMEMORY(space, NULL, double, maxsample);
  stats->Vx = ALLOCMEMORY(space, NULL, double, maxsample); 
  stats->RP = ALLOCMEMORY(space, NULL, Uint, 100);
  stats->RP_N = ALLOCMEMORY(space, NULL, Uint, 100);
  stats->RQ = ALLOCMEMORY(space, NULL, Uint, 255);
  stats->RQ_N = ALLOCMEMORY(space, NULL, Uint, 255);
  stats->MO = ALLOCMEMORY(space, NULL, Uint, 1024);
  stats->MO_N = ALLOCMEMORY(space, NULL, Uint, 1024);
  memset(stats->S, 0, sizeof(double)*(6*6));
  memset(stats->S_N, 0, sizeof(Uint)*6);
  memset(stats->Sx, 0, sizeof(double)*(6*6));
  memset(stats->Sx_N, 0, sizeof(Uint)*6);
  memset(stats->R, 0, sizeof(Uint)*(100*255));
  memset(stats->RP, 0, sizeof(Uint)*(100));
  memset(stats->RQ, 0, sizeof(Uint)*(255));
  memset(stats->R_N, 0, sizeof(Uint)*(100*255));
  memset(stats->RP_N, 0, sizeof(Uint)*(100));
  memset(stats->RQ_N, 0, sizeof(Uint)*(255));
  memset(stats->MO, 0, sizeof(Uint)*(1024));
  memset(stats->MO_N, 0, sizeof(Uint)*(1024));

  memset(stats->V, 0, sizeof(double)*maxsample);
  memset(stats->Vx, 0, sizeof(double)*maxsample);
  memset(stats->s, 0, sizeof(double)*maxsample);
  return stats;
}


/*--------------------- bl_matchfileDestructSampleStats ----------------------
 *    
 * @brief destruct sample stats
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchfileDestructSampleStats (void *space, matchfileSampleStats_t *stats)
{
  
  FREEMEMORY(space, stats->entropydensity);
  FREEMEMORY(space, stats->e);

  //if(stats->b) FREEMEMORY(space, stats->b);
  //if(stats->s) FREEMEMORY(space, stats->s);
  
  FREEMEMORY(space, stats->eraw);
  FREEMEMORY(space, stats->entropy);
  FREEMEMORY(space, stats->e_mu);
  FREEMEMORY(space, stats->e_sd);
  FREEMEMORY(space, stats->S);
  FREEMEMORY(space, stats->S_N);
  FREEMEMORY(space, stats->Sx);
  FREEMEMORY(space, stats->Sx_N);
  FREEMEMORY(space, stats->R);
  FREEMEMORY(space, stats->R_N);
  FREEMEMORY(space, stats->V);
  FREEMEMORY(space, stats->Vx);
  FREEMEMORY(space, stats->RP);
  FREEMEMORY(space, stats->RQ);
  FREEMEMORY(space, stats->RP_N);
  FREEMEMORY(space, stats->RQ_N);
  FREEMEMORY(space, stats->RR);
  FREEMEMORY(space, stats->MM);
  FREEMEMORY(space, stats->MO);
  FREEMEMORY(space, stats->MO_N);
  
  if(stats->ecdf) {  
    ecdf_destruct(stats->ecdf);
    FREEMEMORY(space, stats->ecdf);
  }

  if(stats->gev_mu) { 
    FREEMEMORY(space, stats->gev_mu);
    FREEMEMORY(space, stats->gev_si);
    FREEMEMORY(space, stats->gev_xi);
    FREEMEMORY(space, stats->gev_ll);
  }

  return ;
}

/*-------------------------- bl_matchfileFrameStats --------------------------
 *    
 * @brief get descriptive statistics for a frame
 * @author Steve Hoffmann 
 *   
 */
 

matchfileFrameStats_t *
bl_matchfileFrameStats (void *space, matchfileFrame_t *frame) {
  Uint j, pos=0, ch, qu, rp, mc, rss, *c, *f, *n, 
       *D, D_ymax=0, D_ymax_1=0, D_xmax=0, noofstarts=0;
  matchfileFrameStats_t* stats;
  int **s;
  double *e, *r, *v, *y, *d, x=0;

  stats = ALLOCMEMORY(space, NULL, matchfileFrameStats_t, 1);
  v = ALLOCMEMORY(space, NULL, double, 256);
  d = ALLOCMEMORY(space, NULL, double, 256);
  e = ALLOCMEMORY(space, NULL, double, 256);
  c = ALLOCMEMORY(space, NULL, double, 256);
  r = ALLOCMEMORY(space, NULL, double, 256);
  y = ALLOCMEMORY(space, NULL, double, 256);
  f = ALLOCMEMORY(space, NULL, Uint, 256);
  D = ALLOCMEMORY(space, NULL, Uint, MAX_STARTSITES);

  memset(v, 0, sizeof(double)*256);
  memset(c, 0, sizeof(double)*256);
  memset(r, 0, sizeof(double)*256);
  memset(y, 0, sizeof(double)*256);
  memset(d, 0, sizeof(double)*256);
  memset(f, 0, sizeof(Uint)  *256);
  memset(D, 0, sizeof(Uint)*MAX_STARTSITES);

  for(j=0; j < 256; j++) e[j] = log10(0);

  for(pos=0; pos < frame->width; pos++) {

    x += frame->cs[pos].len;
    s = ALLOCMEMORY(space, NULL, int*, 256);
    n = ALLOCMEMORY(space, NULL, Uint, 256);
    memset(s, 0, sizeof(int*)*256);
    memset(n, 0, sizeof(Uint)*256);
    f[(int)frame->cs[pos].ref]++;
    rss = frame->cs[pos].starts;
    noofstarts += rss;

    if(rss >= MAX_STARTSITES) {
      rss = MAX_STARTSITES-1;
      D_xmax = MAX_STARTSITES;
    } else {
      D_xmax = (rss > D_xmax) ? rss: D_xmax;
    }
    
    D[rss]++;
    if(D[rss] > D_ymax) {
      D_ymax = D[rss];
    }

    if(rss > 0 && D[rss] >D_ymax_1) {
      D_ymax_1 = D[rss];
    }

    for(j=0; j < frame->cs[pos].len; j++) {
      ch = frame->cs[pos].chars[j];
      qu = frame->cs[pos].quals[j];
      rp = frame->cs[pos].readpos[j];
      mc = frame->cs[pos].matchcnt[j];


      if(ch != frame->cs[pos].ref) {
        d[ch]++;
      }

      s[ch] = ALLOCMEMORY(space, s[ch], Uint, n[ch]+1);
      s[ch][n[ch]] = rp;
      e[ch] = log10add(e[ch],(qu/-10.));
      r[ch] += rp;
      y[ch] += mc;
      c[ch]++;
      n[ch]++;
    }

    for(j=0; j < 256; j++) {
      if(n[j]) v[j] += sqrt(var_int(s[j],n[j]));
      FREEMEMORY(space, s[j]);
    }
    
    FREEMEMORY(space, n);
    FREEMEMORY(space, s);
  }

  for(j=0; j < 256; j++) {
    if(c[j]) e[j] -= log10(c[j]);
    if(c[j]) r[j] /= c[j];
    if(c[j]) y[j] /= c[j];
    if(c[j]) v[j] /= f[j];
    if(c[j]) d[j] /= c[j];
  }


  stats->ntcnt = c;
  stats->mean_err = e;
  stats->mean_sde = v;
  stats->mean_pos = r;
  stats->mean_mul = y;
  stats->mean_dis = d;
  stats->char_frq = f;
  stats->mean_cov = x/frame->width;
  stats->dist_rss_xmax = D_xmax;
  stats->dist_rss_ymax = D_ymax;
  stats->dist_rss_ymax_1 = D_ymax_1;
  stats->dist_rss = D;
  stats->rss = noofstarts;

  return stats;
}



/*------------------------ bl_matchfileDestructFrame -------------------------
 *    
 * @brief wrap the frame
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchfileDestructFrame(void *space, matchfileFrame_t *frame)
{
  
  bl_matchfileDestructCross(space, frame->cs, frame->width);
  FREEMEMORY(space, frame->cs);
  FREEMEMORY(space, frame);

  return ;
}

/*---------------------- bl_matchfileDestructFrameStats ----------------------
 *    
 * @brief destruct the frame statistics structure
 * @author Steve Hoffmann 
 *   
 */
 

void
bl_matchfileDestructFrameStats(void *space, matchfileFrameStats_t *stats) {
  FREEMEMORY(space, stats->mean_err);
  FREEMEMORY(space, stats->mean_sde);
  FREEMEMORY(space, stats->mean_pos);
  FREEMEMORY(space, stats->mean_mul);
  FREEMEMORY(space, stats->mean_dis);
  FREEMEMORY(space, stats->dist_rss);
  FREEMEMORY(space, stats->char_frq);
  FREEMEMORY(space, stats->ntcnt);
}

/*------------------------ bl_matchfileDumpSampleStats -------------------------
 *    
 * @brief dump the stats
 * @author Steve Hoffmann 
 *   
 */

  void
bl_matchfileDumpSampleStats (matchfileSampleStats_t *stats)
{

  double mx =.0 , sx = .0, kx = .0;
  //double ws[]={0.95, 0.05};
  Uint i;
  //Uint j;
  //
/*
  stats->V_mu = .1;
  stats->V_sd = .1;
  for(i=0; i < stats->entropydensitylen; i++) {
    fprintf(stderr, "%.6f\t%.6f\n", (double)i*0.05, stats->entropydensity[i]);
  }
stats->V_ll = 0;
  stats->Vx_mu = .1;
  stats->Vx_sd = .1;
  stats->Vx_ll = 0;

  ws[0] = 1.0;
  fprintf(stderr, "fitting %d\n", stats->V_N);

  stats->V_ll=gmm(NULL, stats->V, stats->V_N, 1, 
      &stats->V_mu, &stats->V_sd, ws, 1, 100000); 

  stats->Vx_ll=gmm(NULL, stats->Vx, stats->Vx_N, 1, 
      &stats->Vx_mu, &stats->Vx_sd, ws, 1, 100000);

*/
  for(i=0; i < stats->e_N; i++) {
    fprintf(stderr, "%.6f\n", stats->e[i]);
  }
  
  fprintf(stderr, "-------------\n");

  for(i=0; i < stats->s_N; i++) {
    fprintf(stderr, "%.6f\n", stats->s[i]);
  }

 // return;


  fprintf(stderr, "px: %f\n",stats->px);
  fprintf(stderr, "pxx: %f\n",stats->pxx);

  fprintf(stderr, "edensity:\n");

  fprintf(stderr, "mu={%f }, sd=%f, ll=%f\n", 
      stats->e_mu[0], stats->e_sd[0], stats->e_ll);
 
  fprintf(stderr, "mu={%f }, sd=%f, ll=%f\n", 
      stats->e_mu[1], stats->e_sd[1], stats->e_ll);


  /*
     fprintf(stderr, "noise matrix\n");

     for(i=30; i < 100; i++) {
     for(j=0; j < 100; j++) {
     if(MATRIX2D(stats->R_N, 255, j, i))
     fprintf(stderr, "%d %d %f\n", i, j, 
     (double)MATRIX2D(stats->R, 255, j, i)/MATRIX2D(stats->R_N, 255, j, i));
     }
     }
     */

  fprintf(stderr, "P=%d, X=%d, N=%d\n", stats->P, stats->X, stats->N);
  fprintf(stderr, "P(X)=%f log:%f\n",(double)stats->X/stats->N, log((double)stats->X/stats->N));
  fprintf(stderr, "P(N)=%f log:%f\n",(double)stats->P/stats->N, log((double)stats->P/stats->N));

  if(stats->MO) { 
    fprintf(stderr, "motif");
    for(i=0; i < 1024; i++) {
      if(stats->MO_N[i]) {
        fprintf(stderr, "%d %d %d %f\n", i, stats->MO[i], stats->MO_N[i], log((double)stats->MO[i]/(double)stats->MO_N[i]));
      }
    }
  }

  fprintf(stderr, "readpos\n");
  for(i=0; i < 100; i++) {
    if(stats->RP_N[i])
      fprintf(stderr, "%d %d %d %f\n", i, stats->RP[i], stats->RP_N[i], log((double)stats->RP[i]/(double)stats->RP_N[i]));
  }

  fprintf(stderr, "readqual\n");
  for(i=0; i < 255; i++) {
    if(stats->RQ_N[i])
      fprintf(stderr, "%d %d %d %f\n", i, stats->RQ[i], stats->RQ_N[i], (double)stats->RQ[i] / stats->RQ_N[i]);
  }

  fprintf(stderr, "readerror\n");
  for(i=0; i < 11; i++) {
    fprintf(stderr, "%d %d %d %f\n", i, stats->RR[i], stats->RR_N, (double)stats->RR[i] / stats->RR_N);
  }

  fprintf(stderr, "multiple matches\n");
  for(i=0; i < 50; i++) {
    fprintf(stderr, "%d %d %d %f\n", i, stats->MM[i], stats->MM_N, (double)stats->MM[i] / stats->MM_N);
  }

  fprintf(stderr, "readstartvar\n");
  for(i=0; i < stats->V_N; i++) {
    fprintf(stderr, "%f\n", stats->V[i]);
  }
 
  fprintf(stderr, "readstartvar gaussian model\n");

  fprintf(stderr, "mu=%f, sd=%f, ll=%f\n", 
      stats->V_mu, stats->V_sd, stats->V_ll);
  
  fprintf(stderr, "readstartvar X");

  for(i=0; i < stats->Vx_N; i++) {
    fprintf(stderr, "%f\n", stats->Vx[i]);
  }
 
  fprintf(stderr, "readstartvar X gaussian model\n");

  fprintf(stderr, "mu=%f, sd=%f, ll=%f\n", 
      stats->Vx_mu, stats->Vx_sd, stats->Vx_ll);

  if(stats->b) { 
    fprintf(stderr, "strand bias\n");
    for(i=0; i < stats->e_N; i++) {
      fprintf(stderr, "%.6f\n", stats->b[i]);
    }
  }

  fprintf(stderr, "raw error\n");
 
  for(i=0; i < stats->e_N; i++) {
    fprintf(stderr, "%.6f\n", stats->eraw[i]);
  }
  
  fprintf(stderr, "adjust error\n");
 
  for(i=0; i < stats->e_N; i++) {
    fprintf(stderr, "%.6f\n", stats->e[i]);
  }

  fprintf(stderr, "entropy\n");

  for(i=0; i < stats->e_N; i++) {
    fprintf(stderr, "%.6f\n", stats->entropy[i]);
  }

  fprintf(stderr, "entropydensity: %d\n", stats->entropydensitylen);

  for(i=0; i < stats->entropydensitylen; i++) {
    fprintf(stderr, "%.6f\t%.6f\n", (double)i*0.05, stats->entropydensity[i]);
  }
 
  qsort(stats->eraw, stats->e_N, sizeof(double), cmp_dbl_qsort);
  gevLmoment(stats->eraw, stats->e_N, &mx, &sx, &kx);

  NFO("lmoment m:%f, s:%f, xi:%f \n", mx, sx, -kx);
  NFO("log-likelihood:%f\n", gevll(stats->eraw, stats->e_N, mx, sx, kx));
 
  gevmle(NULL, stats->eraw, stats->e_N, &mx, &sx, &kx, 10000, stats->eraw[0], stats->eraw[stats->e_N-1]);
  NFO("gev m:%f, s:%f, xi:%f \n", mx, sx, -kx);
  NFO("log-likelihood:%f\n", gevll(stats->eraw, stats->e_N, mx, sx, kx));

  stats->px = mx;
  while(gevcdf(stats->px, mx, sx, -kx) < 0.99) 
    stats->px += 0.00001;

  fprintf(stderr, "eraw px: %f\n", stats->px);

  qsort(stats->e, stats->e_N, sizeof(double), cmp_dbl_qsort);
  gevLmoment(stats->e, stats->e_N, &mx, &sx, &kx);

  NFO("lmoment m:%f, s:%f, xi:%f \n", mx, sx, -kx);
  NFO("log-likelihood:%f\n", gevll(stats->e, stats->e_N, mx, sx, kx));
  gevmle(NULL, stats->e, stats->e_N, &mx, &sx, &kx, 10000, stats->e[0], stats->e[stats->e_N-1]);
  NFO("gev m:%f, s:%f, xi:%f \n", mx, sx, -kx);
  NFO("log-likelihood:%f\n", gevll(stats->e, stats->e_N, mx, sx, kx));

  stats->px = mx;
  while(gevcdf(stats->px, mx, sx, -kx) < 0.99) 
    stats->px += 0.00001;

  fprintf(stderr, "99: %f\n", stats->px);
   
  stats->px = mx;
  while(gevcdf(stats->px, mx, sx, -kx) < 0.95) 
    stats->px += 0.00001;

  fprintf(stderr, "95: %f\n", stats->px);
  
  stats->px = mx;
  while(gevcdf(stats->px, mx, sx, -kx) < 0.90) 
    stats->px += 0.00001;

  fprintf(stderr, "90: %f\n", stats->px);


  return ;
}
/*------------------------ bl_matchfileTabulateFull -------------------------
 *    
 * @brief full tabulation of all factors for the calculation of snvs
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchfileTabulateFull (FILE *dev, Uint pos, matchfileCrossStats_t* cvcss, matchfileCrossStats_t* cnvcss,  matchfileCrossStats_t* rvcss, matchfileCrossStats_t* rnvcss,   
    matchfileCross_t *cs, matchfileSampleStats_t *stats, char ref, char cons, double p_cons, double p_ref, double p_consx, double p_refx)
{
  double *nv;
  char *ch, *rq, *st;
  unsigned char *ed;
  uint32_t *rp, *mc, i, rpos, *rl, errors=0;
  double e;
  double *rvrt, *rvrq, *rvrr, *rvrv, *rvmm, rvsp;
  double *rnrt, *rnrq, *rnrr, *rnrv, *rnmm, rnsp;
  double *cvrt, *cvrq, *cvrr, *cvrv, *cvmm, cvsp;
  double *cnrt, *cnrq, *cnrr, *cnrv, *cnmm, cnsp;
  double rvee, cvee, rnee, cnee;
 
  ch = cs->chars;
  st = cs->strands;
  rq = cs->quals;
  ed = cs->edist;
  rp = cs->readpos;
  nv = bl_matchfileGetNTReadPosVar(cs);
  mc = cs->matchcnt;
  rl = cs->readlen;

  rvrt = rvcss->var_rt;
  double rvrt_sum = .0;
  rvrq = rvcss->var_rq;
  double rvrq_sum = .0;
  rvrr = rvcss->var_rr;
  double rvrr_sum = .0;
  rvrv = rvcss->var_rv;
  double rvrv_sum = .0;
  rvmm = rvcss->var_mm;
  double rvmm_sum = .0;
  rvee = rvcss->var_ee;
  rvsp = rvcss->strandpenalty;

  cvrt = cvcss->var_rt;
  cvrq = cvcss->var_rq;
  cvrr = cvcss->var_rr;
  cvrv = cvcss->var_rv;
  cvmm = cvcss->var_mm;
  cvee = cvcss->var_ee;
  cvsp = cvcss->strandpenalty; 
   
  rnrt = rnvcss->var_rt;
  double rnrt_sum = .0;
  rnrq = rnvcss->var_rq;
  double rnrq_sum = .0;
  rnrr = rnvcss->var_rr;
  double rnrr_sum = .0;
  rnrv = rnvcss->var_rv;
  double rnrv_sum = .0;
  rnmm = rnvcss->var_mm;
  double rnmm_sum = .0;
  rnee = rnvcss->var_ee;
  rnsp = .0;

  cnrt = cnvcss->var_rt;
  cnrq = cnvcss->var_rq;
  cnrr = cnvcss->var_rr;
  cnrv = cnvcss->var_rv;
  cnmm = cnvcss->var_mm;
  cnee = cnvcss->var_ee;
  cnsp = .0;



  fprintf(dev, "%d\t%c\t%c\t%f\t%f\t%f\t%f\t", pos, ref, cons, p_cons, 
      p_ref, p_consx, p_refx);

  for(i=0; i < cs->len; i++) {  
    fprintf(dev, "%c",  ch[i]);
      if((int)ch[i] != (int)ref && ch[i] != 'N') {
        errors++;
      }
  }
  fprintf(dev, "\t");

  for(i=0; i < cs->len; i++) {  
    fprintf(dev, "%c",  st[i]);
  }
  fprintf(dev, "\t");


  e = (double)errors/(double)cs->len;

  fprintf(dev, "%d\t%f\t", errors, (double)errors/(double)cs->len);
  
  for(i=0; i < cs->len; i++) fprintf(dev, "%c",  rq[i]);
  fprintf(dev, "\t");
  
  for(i=0; i < cs->len-1; i++) fprintf(dev, "%d,",  ed[i]);
  fprintf(dev, "%d\t", ed[cs->len-1]);
  
  for(i=0; i < cs->len-1; i++) fprintf(dev, "%d,",  rp[i]);
  fprintf(dev, "%d\t", rp[cs->len-1]);

  for(i=0; i < cs->len-1; i++) fprintf(dev, "%d,",  mc[i]);
  fprintf(dev, "%d\t", mc[cs->len-1]);

  for(i=0; i < cs->len-1; i++) { 
    rpos = trunc(((double)(((double)rp[i]*100.0)/((double)rl[i]))));
    fprintf(dev, "%d,", rpos);
  }
  rpos = trunc(((double)(((double)rp[cs->len-1]*100.0)/((double)rl[cs->len-1]))));
  fprintf(dev, "%d\t",  rpos);

  for(i=0; i < cs->len-1; i++) fprintf(dev, "%f,", nv[(int)ch[i]]);
  fprintf(dev, "%f\t", nv[(int)ch[i]]);


  for(i=0; i < cs->len-1; i++) { 
    rvrt_sum += rvrt[i];
    fprintf(dev, "%f,",  rvrt[i]);
  }
  rvrt_sum += rvrt[i];
  fprintf(dev, "%f\t", rvrt[cs->len-1]);
 
  for(i=0; i < cs->len-1; i++) { 
    rvrq_sum += rvrq[i];
    fprintf(dev, "%f,",  rvrq[i]);
  }
  rvrq_sum += rvrq[i];
  fprintf(dev, "%f\t", rvrq[cs->len-1]);
 
  for(i=0; i < cs->len-1; i++) { 
    rvrr_sum += rvrr[i];
    fprintf(dev, "%f,",  rvrr[i]);
  }
  rvrr_sum += rvrr[i];
  fprintf(dev, "%f\t", rvrr[cs->len-1]);
 
  for(i=0; i < cs->len-1; i++) { 
    rvrv_sum += rvrv[i]; 
    fprintf(dev, "%f,",  rvrv[i]);
  }
  rvrv_sum += rvrv[i];
  fprintf(dev, "%f\t", rvrv[cs->len-1]);
 
  for(i=0; i < cs->len-1; i++) { 
    rvmm_sum += rvmm[i];
    fprintf(dev, "%f,",  rvmm[i]);
  }
  rvmm_sum += rvmm[i];
  fprintf(dev, "%f\t", rvmm[cs->len-1]);

  fprintf(dev, "%f\t", rvsp);

  fprintf(dev, "%f\t", rvee);

  for(i=0; i < cs->len-1; i++) { 
    rnrt_sum += rnrt[i];
    fprintf(dev, "%f,",  rnrt[i]);
  }
  rnrt_sum += rnrt[i];
  fprintf(dev, "%f\t", rnrt[cs->len-1]);
 
  for(i=0; i < cs->len-1; i++) { 
    rnrq_sum += rnrq[i];
    fprintf(dev, "%f,",  rnrq[i]);
  }
  rnrq_sum += rnrq[i];
  fprintf(dev, "%f\t", rnrq[cs->len-1]);
 
  for(i=0; i < cs->len-1; i++) { 
    rnrr_sum += rnrr[i];
    fprintf(dev, "%f,",  rnrr[i]);
  }
  rnrr_sum += rnrr[i];
  fprintf(dev, "%f\t", rnrr[cs->len-1]);
 
  for(i=0; i < cs->len-1; i++) { 
    rnrv_sum += rnrv[i];
    fprintf(dev, "%f,",  rnrv[i]);
  }
  rnrv_sum += rnrv[i];
  fprintf(dev, "%f\t", rnrv[cs->len-1]);
 
  for(i=0; i < cs->len-1; i++) { 
    rnmm_sum += rnmm[i];
    fprintf(dev, "%f,",  rnmm[i]);
  }
  rnmm_sum += rnmm[i];
  fprintf(dev, "%f\t", rnmm[cs->len-1]);
 
  fprintf(dev, "%f\t", rnsp);
  
  fprintf(dev, "%f\t", rnee);
  
  for(i=0; i < cs->len-1; i++) fprintf(dev, "%f,",  cvrt[i]);
  fprintf(dev, "%f\t", cvrt[cs->len-1]);
 
  for(i=0; i < cs->len-1; i++) fprintf(dev, "%f,",  cvrq[i]);
  fprintf(dev, "%f\t", cvrq[cs->len-1]);
 
  for(i=0; i < cs->len-1; i++) fprintf(dev, "%f,",  cvrr[i]);
  fprintf(dev, "%f\t", cvrr[cs->len-1]);
 
  for(i=0; i < cs->len-1; i++) fprintf(dev, "%f,",  cvrv[i]);
  fprintf(dev, "%f\t", cvrv[cs->len-1]);
 
  for(i=0; i < cs->len-1; i++) fprintf(dev, "%f,",  cvmm[i]);
  fprintf(dev, "%f\t", cvmm[cs->len-1]);

  fprintf(dev, "%f\t", cvsp);

  fprintf(dev, "%f\t", cvee);

  for(i=0; i < cs->len-1; i++) fprintf(dev, "%f,",  cnrt[i]);
  fprintf(dev, "%f\t", cnrt[cs->len-1]);
 
  for(i=0; i < cs->len-1; i++) fprintf(dev, "%f,",  cnrq[i]);
  fprintf(dev, "%f\t", cnrq[cs->len-1]);

  for(i=0; i < cs->len-1; i++) fprintf(dev, "%f,",  cnrr[i]);
  fprintf(dev, "%f\t", cnrr[cs->len-1]);
 
  for(i=0; i < cs->len-1; i++) fprintf(dev, "%f,",  cnrv[i]);
  fprintf(dev, "%f\t", cnrv[cs->len-1]);
 
  for(i=0; i < cs->len-1; i++) fprintf(dev, "%f,",  cnmm[i]);
  fprintf(dev, "%f\t", cnmm[cs->len-1]);

  fprintf(dev, "%f\t", cnsp);

  fprintf(dev, "%f\t", cnee);

  for(i=0; i < cs->len-1; i++) fprintf(dev, "%f,", (double)stats->RP[rpos]);
  fprintf(dev, "%f\t", (double)stats->RP[rpos]);
   
  for(i=0; i < cs->len-1; i++) fprintf(dev, "%f,", (double)stats->RP_N[rpos]);
  fprintf(dev, "%f\t", (double)stats->RP_N[rpos]);

  for(i=0; i < cs->len-1; i++) fprintf(dev, "%f,", (double)stats->RR[MIN(ed[i],10)]);
  fprintf(dev, "%f\t", (double)stats->RR[MIN(ed[i],10)]);

  fprintf(dev, "%f\t", (double)stats->RR_N);
 
  for(i=0; i < cs->len-1; i++) fprintf(dev, "%f,", 
      univarnormcdf(nv[(int)ch[i]], 
        stats->V_mu-(stats->V_sd/SDFRACTION), stats->V_sd));
  fprintf(dev, "%f\t", univarnormcdf(nv[(int)ch[i]], 
        stats->V_mu-(stats->V_sd/SDFRACTION), stats->V_sd));

  for(i=0; i < cs->len-1; i++) fprintf(dev, "%f,", (double)stats->MM[MIN(mc[i],10)]);
  fprintf(dev, "%f\t", (double)stats->MM[MIN(mc[i],10)]);

  fprintf(dev, "%f\t", (double)stats->MM_N);

  fprintf(dev, "%f\n", gevcdf(e, stats->gev_mu[0], stats->gev_si[0], stats->gev_xi[0]));


  return ;
}


/*------------------------- bl_matchfilePrintVariant -------------------------
 *    
 * @brief print variant
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchfileCrossStatsPrint (FILE *fp, matchfileCrossStats_t *css, matchfileCross_t *cs, 
    matchfileSampleStats_t *stats, char ref)
{

  double *nv;
  char *ch, *rq;
  unsigned char *ed;
  uint32_t *rp, *mc, i, rpos, *rl;

 
  rp = cs->readpos;
  rq = cs->quals;
  ed = cs->edist;
  ch = cs->chars;
  nv = bl_matchfileGetNTReadPosVar(cs);
  mc = cs->matchcnt;
  rl = cs->readlen;

  fprintf(fp, "----stats----\n");
  
  for(i=0; i < cs->len; i++) { 

    rpos = trunc(((double)(((double)rp[i]*100.0)/((double)rl[i]))));
    
    fprintf(fp, "------------------------------");
    fprintf(fp, "nucleotide %i\n", i);
    fprintf(fp, "P(%c -> %c) = %f\n", ref, ch[i], css->var_s[i]);
    fprintf(fp, "RP(%d)=%f\n", (int)rpos, css->var_rt[i]);
    fprintf(fp, "RQ(%d)=%f (p:%f, %f)\n", (int)rq[i], css->var_rq[i], pow(10,((double)((double)rq[i]-64.0)/-10.0)), log(pow(10,((double)((double)rq[i]-64.0)/-10.0))) );
//    fprintf(fp, "(1-RQ + 1-RP)/2=%f\n", logadd(css->var_rt[i], css->var_rq[i]) - log(2));
    fprintf(fp, "RR(%d)=%f\n", ed[i], css->var_rr[i]);
    fprintf(fp, "RV(%f)=%f (mu: %f sd:%f)\n", nv[(int)ch[i]], css->var_rv[i],
        stats->V_mu-(stats->V_sd/SDFRACTION), stats->V_sd);
    fprintf(fp, "MM(%d)=%f (%d/%d )\n", mc[i], css->var_mm[i], 
        stats->MM[MIN(mc[i],10)], stats->MM_N);
    fprintf(fp, "subtotal: %f\n", css->sub[i]);
  }
 
  fprintf(fp, "------------------------------");
  fprintf(fp, "pentropy %f\n", css->pentropy);
  fprintf(fp, "strandpenalty %f\n", css->strandpenalty);

  return ;
}

/*----------------------------- bl_matchfileTestCheck -----------------------------
 *    
 * @brief test
 * @author Steve Hoffmann 
 *   
 *
 
Uint
bl_matchfileTestCheck(void *space, Uint fidx, Uint cidx, Uint pos, matchfileCross_t *cs, 
    char ref, matchfileSampleStats_t *stats, unsigned char show, void *nfo)
{
  matchfileCrossStats_t cnvcss, cvcss, rnvcss, rvcss; 
  FILE *fp = NULL;
  double p_cons = .0, p_consx = .0, p_ref = .0, p_refx = .0, PX, P; 
  Uint *cnt, i;
  
  cs->s_cons = 1;
  cs->s_consx = 0;
  cs->s_ref = 1;
  cs->s_refx = 0;
  cs->p_hom = log(0);
  
  //override
  //stats->maxcover = 10000;
  
  fp = stdout;

  if(!stats || cs->len < stats->mincover 
      || cs->len > stats->maxcover) { 
    return 0;
  }

  if(cs->len && !cs->cons){
    cnt = bl_matchfileGetNTCounts(cs);
    cs->cons = (char) uarraymax(cnt, 255);
    FREEMEMORY(space, cnt);
  }
  
  for(i=0; i < cs->len; i++) {
    if(cs->chars[i] != cs->cons || 
        cs->chars[i] != cs->ref) break;
  }

  if(i == cs->len) {
    return 0;
  }
  
  if(!stats->standardized) bl_matchfileGetStandardization (space, stats);

  PX = log((double)stats->X/stats->N);
  P = log((double)stats->P/stats->N);
  //PX = MAX(PX, -3.5);
  //P =  MAX(P, -0.03);

  PX = log(0.01);
  P = log(0.99);
 
  p_cons  = bl_matchfileTestNonVariant (cs, stats, cs->cons, &cnvcss, stats->minrp, stats->maxrp, stats->minrq, stats->maxrq);
  p_consx = bl_matchfileTestVariant (cs, stats, cs->cons, &cvcss, stats->minrp, stats->maxrp, stats->minrq, stats->maxrq);
  p_ref  = bl_matchfileTestNonVariant (cs, stats, ref, &rnvcss, stats->minrp, stats->maxrp, stats->minrq, stats->maxrq);
  p_refx = bl_matchfileTestVariant (cs, stats, ref, &rvcss, stats->minrp, stats->maxrp, stats->minrq, stats->maxrq);


  cs->p_cons = P+p_cons;
  cs->p_consx = PX+p_consx;

  cs->p_ref = P+p_ref;
  cs->p_refx = PX+p_refx;


  bl_matchfileTabulateFull (fp, pos, &cvcss, &cnvcss, &rvcss, &rnvcss, cs, stats, 
      ref, cs->cons, cs->p_cons, cs->p_ref, cs->p_consx, cs->p_refx);

  bl_matchfileCrossStatsDestruct (&cnvcss);
  bl_matchfileCrossStatsDestruct (&cvcss); 
  bl_matchfileCrossStatsDestruct (&rnvcss);
  bl_matchfileCrossStatsDestruct (&rvcss); 
 
  
  return 0;

}

*/



/*------------------------ bl_matchfileCrossStatsInit ------------------------
 *    
 * @brief initalize cross stats
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchfileCrossStatsInit (matchfileCrossStats_t *css, Uint len)
{

  css->len = len;
  css->var_s = ALLOCMEMORY(space, NULL, double, len);
  css->var_rt = ALLOCMEMORY(space, NULL, double, len);
  css->var_rq = ALLOCMEMORY(space, NULL, double, len);
  css->var_rr = ALLOCMEMORY(space, NULL, double, len);
  css->var_rv = ALLOCMEMORY(space, NULL, double, len);
  css->var_mm = ALLOCMEMORY(space, NULL, double, len);
  css->sub = ALLOCMEMORY(space, NULL, double, len);
  css->var_ee = .0;
  css->pentropy = .0;
  css->strandpenalty = .0;

  return ;
}

/*---------------------- bl_matchfileDestructCrossStats ----------------------
 *    
 * @brief destruct cross stats
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchfileCrossStatsDestruct (matchfileCrossStats_t *css)
{

  if (css->var_s) FREEMEMORY(space, css->var_s);
  if (css->var_rt) FREEMEMORY(space, css->var_rt);
  if (css->var_rq) FREEMEMORY(space, css->var_rq);
  if (css->var_rr) FREEMEMORY(space, css->var_rr);
  if (css->var_rv) FREEMEMORY(space, css->var_rv);
  if (css->var_mm) FREEMEMORY(space, css->var_mm);
  if (css->sub) FREEMEMORY(space, css->sub);

  return ;
}

/*----------------------- bl_matchfilePrintTestResult ------------------------
 *    
 * @brief dump the test result
 * @author Steve Hoffmann 
 *   
 */

  void
bl_matchfileTestPrint (matchfileFrame_t *f, Uint p)
{

  char type;
  
  if(!isinf(f->cs[p].p_hom)) {
    type = 'H';
  } else if(f->cs[p].p_consx > f->cs[p].p_cons && 
      f->cs[p].p_refx > f->cs[p].p_ref) { 
    type = 'B';
  } else if (f->cs[p].p_consx > f->cs[p].p_cons) { 
    type = 'C';
  } else{ 
    type = 'R';
  }

  printf("%s\t%d\t%c\t%c\t%c\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%d\t%f\t%f\t%f\n", 
      f->chrname, f->start+p, type, 
      f->ref[p], f->cs[p].cons, 
      f->cs[p].chars, f->cs[p].diff_rt, f->cs[p].diff_rq, 
      f->cs[p].diff_rr, f->cs[p].diff_mm, 
      f->cs[p].ee,  f->cs[p].pee, f->cs[p].secondminimum, f->cs[p].secondcnt, f->cs[p].pbinom, f->cs[p].scr_cons, f->cs[p].scr_ref);

  return ;
}


