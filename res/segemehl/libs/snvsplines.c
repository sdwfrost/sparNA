
/*
 *  snvsplines.c
 *  
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 06/18/14 14:37:12 CEST
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
#include "snvsplines.h"
#include "splines.h"

/*-------------------------------- getcutoff ---------------------------------
 *    
 * @brief calculate the cutoff for the calculation
 * @author Steve Hoffmann 
 *   
 */

void
getcutoff (matchfileSampleStats_t *stats, char *histofilename, 
    char *scorefilename, char *cutfilename, char* splinefilename, 
    char* estimatefilename)
{

  Uint nbins = 30;
  double start = 0.0;
  double end = 3.0;
  double p[] = {0.25, 0.50, 0.75};
  double step = 0.01;
  double *X = stats->s;
  double *Xhi, *bins, *y, binsize, *q, maxval, minval;
  double xi, dyi, dyim1, ddyi, ddyim1, yerr;
  double yi;
  double infval = 0;
  FILE *histofp=NULL, *scorefp = NULL, *cutfp = NULL, *splinefp=NULL, *estimatefp=NULL;
  Uint n = stats->s_N, m, i, *cnt;
  splinefit_t fit;

  gsl_bspline_deriv_workspace *dbw = gsl_bspline_deriv_alloc(4);
  gsl_matrix *dB = gsl_matrix_alloc(7, 5); 
  gsl_matrix *cov = gsl_matrix_alloc(7, 7);
  gsl_vector *ddv = gsl_vector_alloc(7);
  gsl_vector *dv = gsl_vector_alloc(7);

  if(histofilename) histofp = fopen(histofilename, "w");
  if(scorefilename) scorefp = fopen(scorefilename, "w");
  if(cutfilename) cutfp = fopen(cutfilename, "w");
  if(splinefilename) splinefp = fopen(splinefilename, "w");
  if(estimatefilename) estimatefp = fopen(estimatefilename, "w");
 

  ecdf_t * secdf = ecdf_init(X, n);
  fprintf(stderr, "n=%d, e->n=%d", n , secdf->n);
  fprintf(stderr, "%f = %f\n", 0.1, ecdf(0.1, secdf));
  fprintf(stderr, "%f = %f\n", 0.25, ecdf(0.25, secdf));
  fprintf(stderr, "%f = %f\n", 0.33, ecdf(0.33, secdf));
  fprintf(stderr, "%f = %f\n", 0.50, ecdf(0.50, secdf));
  fprintf(stderr, "%f = %f\n", 0.66, ecdf(0.66, secdf));
  fprintf(stderr, "%f = %f\n", 0.75, ecdf(0.75, secdf));
  fprintf(stderr, "%f = %f\n", 0.80, ecdf(0.80, secdf));
  fprintf(stderr, "%f = %f\n", 0.90, ecdf(0.90, secdf));
  fprintf(stderr, "%f = %f\n", 0.95, ecdf(0.95, secdf));
  fprintf(stderr, "%f = %f\n", 0.98, ecdf(0.98, secdf));
  fprintf(stderr, "%f = %f\n", 0.99, ecdf(0.99, secdf));
  fprintf(stderr, "%f = %f\n", 0.999, ecdf(0.999, secdf));
  fprintf(stderr, "%f = %f\n", 1.0, ecdf(1.0, secdf));
  fprintf(stderr, "%f = %f\n", 1.5, ecdf(1.5, secdf));
  fprintf(stderr, "%f = %f\n", 1.7, ecdf(1.7, secdf));
  fprintf(stderr, "%f = %f\n", 1.8, ecdf(1.8, secdf));
  fprintf(stderr, "%f = %f\n", 1.9, ecdf(1.9, secdf));
  fprintf(stderr, "%f = %f\n", 2.0, ecdf(2.0, secdf));
  fprintf(stderr, "%f = %f\n", 2.1, ecdf(2.1, secdf));
  fprintf(stderr, "%f = %f\n", 2.2, ecdf(2.2, secdf));
  fprintf(stderr, "%f = %f\n", 2.3, ecdf(2.3, secdf));
  fprintf(stderr, "%f = %f\n", 2.4, ecdf(2.4, secdf));
  fprintf(stderr, "%f = %f\n", 2.5, ecdf(2.5, secdf));
  
  double minrp=100, maxrp=-100; 
  for(i=0; i < 100; i++) {
    if(stats->RP[i] > 1) { 
      maxrp = (maxrp < 1.0-(((double)stats->RP[i]+1.0)/(stats->RP_N[i]+1.0))) ? 
        1.0-(((double)stats->RP[i]+1.0)/(stats->RP_N[i]+1.0)) : maxrp;
  
      minrp = (minrp > 1.0-(((double)stats->RP[i]+1.0)/(stats->RP_N[i]+1.0))) ? 
        1.0-(((double)stats->RP[i]+1.0)/(stats->RP_N[i]+1.0)) : minrp;
    }
  }

  fprintf(stderr, "minrp: %f, maxrp:%f\n", minrp, maxrp);

  for(i=0; i < 100; i++) {
    fprintf(stderr, "%d\t%d\t%f\t%f\t%f\t%f\t%f\n", stats->RP_N[i], stats->RP[i], (((double)stats->RP[i] + 1.0)
            /((double) stats->RP_N[i] + 1.0)), log(((double)stats->RP[i] + 1.0)
            /((double) stats->RP_N[i] + 1.0)), log(1-(((double)stats->RP[i] + 1.0)
            /((double) stats->RP_N[i] + 1.0))), 
          log(minrp) - log(1.0-(((double)stats->RP[i] + 1.0)/((double) stats->RP_N[i] + 1.0))),         
          log(1.0-(((double)stats->RP[i] + 1.0)/((double) stats->RP_N[i] + 1.0))) - log(maxrp));
  }
  
  fprintf(stderr, "\n");


  double t = 0.5;

//  while(ecdf(t, secdf) > 0.96) t+=0.01;
  while(ecdf(t, secdf) <  0.95) t+=0.01;
  
  fprintf(stderr, "%f = %f\n", t, ecdf(t, secdf));
  

  qsort(X, n, sizeof(double), cmp_dbl_qsort);
  end = MAX(X[n-1], X[0]);
  
  Xhi = ALLOCMEMORY(NULL, NULL, double, n);
  memset(Xhi, 0, sizeof(double)*n);

  for(i=0, m=0; i < n; i++) {
    if(X[i] > start && X[i] <= end) {
      Xhi[m] = X[i];
      m++;
    }
    if(scorefp) fprintf(scorefp, "%f\n", X[i]);
  }

  Uint binmin = 0;
  fprintf(stderr, "m %d - %d n\n", m, n);
  cnt = bin(Xhi, m, &bins, &nbins);
  y = ALLOCMEMORY(NULL, NULL, double, nbins);
  memset(y, 0, sizeof(double)*nbins);
  binsize = (bins[nbins-1]-bins[0])/((double)nbins-1.0);
  fprintf(stderr, "%f - %f / %d = %f\n",bins[nbins-1], bins[0], nbins-1, binsize);
  for(i=0; i < nbins; i++) {
    bins[i] += binsize/2.0;
    y[i] = (double)cnt[i];
    if(histofp) fprintf(histofp, "%f\t%f\n", bins[i], y[i]);
    if(!binmin) {
      if(i > 0 && i < nbins-1)
      NFO("y[%d]=%f (of nbins:%d), %f < %f = %d, %f < %f = %d\n", i, y[i], nbins, y[i-1], y[i], y[i-1] > y[i], y[i], y[i+1], y[i] < y[i+1]);
      if(i > 0 && i < nbins-1 && cnt[i-1] > cnt[i] && cnt[i] < cnt[i+1]) {
        binmin = i-1;
      }
    }
  }

  q = quantiles(bins, nbins, p, 3);

  //GSL voodoo
  histsplineglm(bins, y, nbins, q, 3, &fit);
 
  if(estimatefp){ 
      
    for(i = 0; i < fit.c->size; i++) {
       fprintf(estimatefp, "%d\t%f\n", i, gsl_vector_get(fit.c, i));
    }    
  }

  if(splinefp) {     
    for (xi = bins[0]; xi <= bins[nbins-1]; xi += 0.1)
    {

    gsl_bspline_deriv_eval(xi, 0, dB, fit.bw, dbw); 
    gsl_bspline_eval(xi, fit.B, fit.bw);  
    gsl_vector_set(fit.B, 0, 1.0); //set the intercept

      gsl_multifit_linear_est(fit.B, fit.c, cov, &yi, &yerr);
      fprintf(splinefp, "%f\t%f\n", xi, yi);
      fprintf(stderr, "%f\t%f\n", xi, exp(yi));
    }
  }

  //go get the rightmost maximum in spline
  maxval = bins[nbins-1];

  for (xi = bins[nbins-1]; xi > bins[0]+step; xi -= step) {

    gsl_bspline_deriv_eval(xi, 2, dB, fit.bw, dbw);  
    gsl_matrix_get_col (dv, dB, 1);
    gsl_vector_set(dv,  0, 0.0);
    gsl_multifit_linear_est(dv, fit.c, fit.cov, &dyi, &yerr);

    gsl_bspline_deriv_eval(xi-step, 2, dB, fit.bw, dbw);  
    gsl_matrix_get_col (dv, dB, 1);
    gsl_vector_set(dv,  0, 0.0);
    gsl_multifit_linear_est(dv, fit.c, fit.cov, &dyim1, &yerr);

 
    if(dyi <  0 && dyim1 > 0 && maxval == bins[nbins-1]) { 
      maxval = xi; 
    }
  }

  minval = 0;
  infval = 0;
  //go get the minimum before the rightmost maximum
  for (xi = bins[0]+0.01; xi < maxval; xi += 0.01)
  {

      gsl_bspline_deriv_eval(xi, 2, dB, fit.bw, dbw);  
      gsl_matrix_get_col(dv, dB, 1);
      gsl_matrix_get_col(ddv, dB, 2);
      gsl_vector_set(dv, 0, 0.0);
      gsl_vector_set(ddv, 0, 0.0);

      gsl_multifit_linear_est(dv, fit.c, fit.cov, &dyi, &yerr);
      gsl_multifit_linear_est(ddv, fit.c, fit.cov, &ddyi, &yerr);

      gsl_bspline_deriv_eval(xi-step, 2, dB, fit.bw, dbw);  
      gsl_matrix_get_col (dv, dB, 1);
      gsl_matrix_get_col(ddv, dB, 2);
      gsl_vector_set(dv, 0, 0.0);
      gsl_vector_set(ddv, 0, 0.0);

      gsl_multifit_linear_est(dv, fit.c, fit.cov, &dyim1, &yerr);
      gsl_multifit_linear_est(ddv, fit.c, fit.cov, &ddyim1, &yerr);
      
      if (dyi > 0 && dyim1 < 0) {
        minval = xi;
      }
   }


  if(minval > bins[0]+0.11) {
    //maximum and minimum found - now find minimum and the inflection point
    NFO("selected minimum-inflection method (%f)\n", maxval );
    for (xi = minval; xi < maxval; xi += step) {

      gsl_bspline_deriv_eval(xi, 2, dB, fit.bw, dbw);  
      gsl_matrix_get_col(dv, dB, 1);
      gsl_matrix_get_col(ddv, dB, 2);
      gsl_vector_set(dv, 0, 0.0);
      gsl_vector_set(ddv, 0, 0.0);

      gsl_multifit_linear_est(dv, fit.c, fit.cov, &dyi, &yerr);
      gsl_multifit_linear_est(ddv, fit.c, fit.cov, &ddyi, &yerr);

      gsl_bspline_deriv_eval(xi-step, 2, dB, fit.bw, dbw);  
      gsl_matrix_get_col (dv, dB, 1);
      gsl_matrix_get_col(ddv, dB, 2);
      gsl_vector_set(dv, 0, 0.0);
      gsl_vector_set(ddv, 0, 0.0);

      gsl_multifit_linear_est(dv, fit.c, fit.cov, &dyim1, &yerr);
      gsl_multifit_linear_est(ddv, fit.c, fit.cov, &ddyim1, &yerr);

 /*     if (dyi < 0 && dyim1 > 0) {
        //minimum
        minval = xi;
      }
 */    
      if (ddyi < 0 && ddyim1 > 0) {
        //negative inflection
        infval = xi;
        break;
      }
    }
  } else {

/*    minval = 0;
    infval = 0;
    NFO("inflection-inflection method selected (%f) - PLEASE CHECK SCORE DISTRIBUTION!\n", maxval);
    //maximum not found - use the first positive and the last negative inflection point
    for (xi = bins[0]+step; xi < bins[nbins-1]; xi += step) {

      gsl_bspline_deriv_eval(xi, 2, dB, fit.bw, dbw);  
      gsl_matrix_get_col(dv, dB, 1);
      gsl_matrix_get_col(ddv, dB, 2);
      gsl_vector_set(dv, 0, 0.0);
      gsl_vector_set(ddv, 0, 0.0);

      gsl_multifit_linear_est(dv, fit.c, fit.cov, &dyi, &yerr);
      gsl_multifit_linear_est(ddv, fit.c, fit.cov, &ddyi, &yerr);

      gsl_bspline_deriv_eval(xi-step, 2, dB, fit.bw, dbw);  
      gsl_matrix_get_col (dv, dB, 1);
      gsl_matrix_get_col(ddv, dB, 2);
      gsl_vector_set(dv, 0, 0.0);
      gsl_vector_set(ddv, 0, 0.0);

      gsl_multifit_linear_est(dv, fit.c, fit.cov, &dyim1, &yerr);
      gsl_multifit_linear_est(ddv, fit.c, fit.cov, &ddyim1, &yerr);
     
      if (ddyi > 0 && ddyim1 < 0 && minval == 0) {
        //first positive inflection
        minval = xi;
      }
     
      if (ddyi < 0 && ddyim1 > 0 && minval != 0) {
        //second negative inflection
        infval = xi;
        break;
      }
    }*/
  
    stats->cut = t;
  }
    
 
  NFO("minimum bin: %d\n", binmin);
  if(minval < bins[0]+0.11 || infval < bins[0]+0.11) {  
    if(binmin) {
      stats->cut = bins[binmin];
      NFO("setting cutoff: %f [%f,%f]\n", stats->cut, minval, infval);
    } else { 
      stats->cut = t;
      NFO("ATTENTION! NO CUTOFF FOUND! PLEASE CHECK SCORE DISTRIBUTION - SETTING TO %f\n", stats->cut);
    }
  } else { 
  //  stats->cut = (minval + infval)/2.0;
   
    stats->cut = (minval + minval)/2.0;
    stats->cut = MIN(10.0, stats->cut);
    NFO("setting cutoff: %f [%f,%f]\n", stats->cut, minval, infval);
  }

  if(cutfp) {
    fprintf(cutfp, "%f\t", minval);
    fprintf(cutfp, "%f\t", infval);
    fprintf(cutfp, "%f\n", stats->cut);
  }

  if(histofilename) fclose(histofp);
  if(scorefilename) fclose(scorefp);
  if(cutfilename) fclose(cutfp);
  if(splinefilename) fclose(splinefp);
  if(estimatefilename) fclose(estimatefp);

  destructSplineFit(&fit);  

  gsl_bspline_deriv_free(dbw);
  gsl_matrix_free(dB);
  gsl_matrix_free(cov);
  gsl_vector_free(dv);
  gsl_vector_free(ddv);

  FREEMEMORY(NULL, bins);
  FREEMEMORY(NULL, q);
  FREEMEMORY(NULL, cnt);
  FREEMEMORY(NULL, y);
  FREEMEMORY(NULL, Xhi);
  return ;
}
