
/*
 *  splines.c
 *  
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 06/18/14 14:06:03 CEST
 *  
 */
#include "mathematics.h"
#include "sort.h"
#include <float.h>
#include <string.h>
#include <limits.h>
#include "708.h"
#include <math.h>
#include <complex.h>
#include "splines.h"
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_bspline.h>


/*---------------------------- destructSplineFit -----------------------------
 *    
 * @brief destruct the splinefit 
 * @author Steve Hoffmann 
 *   
 */
 
void
destructSplineFit (splinefit_t* fit)
{
   
  gsl_vector_free(fit->B);
  gsl_vector_free(fit->c);
  gsl_matrix_free(fit->cov);
  gsl_bspline_free(fit->bw);

   return ;
}

/*----------------------------- binsplinepoisson -----------------------------
 *    
 * @brief fits a spline basis to histogram using glm poisson
 * @author Steve Hoffmann 
 *   
 */
 
gsl_vector*
histsplineglm(double *myx, double *myy, int myn, double* q, Uint myqn, 
    splinefit_t *fit) {
  size_t n = myn; 
  const size_t K = 4;
  const size_t ncoeffs = 7;
  const size_t nbreak = ncoeffs-2;
  size_t i, j;
  gsl_bspline_workspace *bw;
  gsl_vector *B;
  gsl_rng *r;
  gsl_vector *c, *w;
  gsl_vector *x, *y;
  gsl_vector *breaks;
  gsl_matrix *X, *cov;
  gsl_multifit_linear_workspace *mw;

  gsl_rng_env_setup();
  r = gsl_rng_alloc(gsl_rng_default);

  /* alloc a cubic bspline workspace (k = 4) */ 
  bw = gsl_bspline_alloc(K, nbreak);
  B = gsl_vector_alloc(ncoeffs);

  x = gsl_vector_alloc(n);
  y = gsl_vector_alloc(n);
  X = gsl_matrix_alloc(n, ncoeffs);
  c = gsl_vector_alloc(ncoeffs);
  breaks = gsl_vector_alloc(nbreak);
  w = gsl_vector_alloc(n);
  cov = gsl_matrix_alloc(ncoeffs, ncoeffs);
  mw = gsl_multifit_linear_alloc(n, ncoeffs);

  
  for(i=0; i < myn; i++) {
    gsl_vector_set(x, i, myx[i]);
    gsl_vector_set(y, i, myy[i]);
    gsl_vector_set(w, i, 1.0);
  }

  gsl_vector_set(breaks, 0, myx[0]); //myx[0] 
  gsl_vector_set(breaks, 1, q[0]);
  gsl_vector_set(breaks, 2, q[1]);
  gsl_vector_set(breaks, 3, q[2]);
  gsl_vector_set(breaks, 4, myx[myn-1]); //myx[myn-1]

  gsl_bspline_knots(breaks, bw);

  /* construct the fit matrix X */
  for (i = 0; i < n; ++i){
      double xi = gsl_vector_get(x, i);

      /* compute B_j(xi) for all j */
      gsl_bspline_eval(xi, B, bw);
  

      /* fill in row i of X */
      for (j = 0; j < ncoeffs; ++j)
        {
          double Bj = gsl_vector_get(B, j);
          gsl_matrix_set(X, i, j, Bj);
        }
  }
 
  /* do the fit */ 
  irls(X, w, y, c, cov, mw, n, ncoeffs); 

  fit->B = B;
  fit->c = c;
  fit->cov = cov;
  fit->bw = bw;

  gsl_rng_free(r);
  gsl_vector_free(x);
  gsl_vector_free(breaks);
  gsl_vector_free(y);
  gsl_matrix_free(X);
  gsl_vector_free(w);
  gsl_multifit_linear_free(mw);

  return NULL;

}


/*----------------------------------- irls -----------------------------------
 *    
 * @brief iterative least square method to fit poisson family function
 * @author Steve Hoffmann
 *
 * TODO: implement SVD and matrix methods to remove gsl lib dependency
 *   
 */
 
gsl_vector *
irls(gsl_matrix *X,  gsl_vector *w, gsl_vector *y, gsl_vector *c, gsl_matrix *cov, 
    gsl_multifit_linear_workspace *mw,  Uint n, Uint ncoeffs) {

  Uint i;
  size_t rank;
  double oldchisq =-1, chisq, yi, mui, zi;
  gsl_vector *z = gsl_vector_calloc(n);
  gsl_vector *mu = gsl_vector_calloc(n);

  //init the mu
  for(i=0; i < n; i++) {
    yi = gsl_vector_get(y,i);
    if(yi == 0) { 
      gsl_vector_set(mu, i, 0.01);
    } else { 
      gsl_vector_set(mu, i, yi);
    }
    //remove intercept
    gsl_matrix_set(X, i, 0, 1.0); 
  }

  //IRLS iteration
  while(1){
    
    for(i=0; i < n; i++) {
      //calculate Z - log transformed response value
      yi = gsl_vector_get(y,i);
      mui = gsl_vector_get(mu, i);
      zi = log(mui) + (1.0/mui)*(yi-mui); 
      gsl_vector_set(z, i, zi);
      //calculate W
      gsl_vector_set(w, i, mui);
    }

    gsl_multifit_wlinear_svd(X, w, z, GSL_DBL_EPSILON, &rank, c, cov, &chisq, mw);
    
    if(fabs(chisq - oldchisq) < 1e-12){ 
      break;
    }

    oldchisq = chisq;
  
    //compute new mu
    gsl_matrix *B = gsl_matrix_calloc(ncoeffs,1);
    for(i=0; i < ncoeffs; i++) {
      gsl_matrix_set(B, i, 0, gsl_vector_get(c, i));
    }

    gsl_matrix *C = gsl_matrix_calloc(n, 1);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, X, B, 0, C); //C = XB

    for(i=0; i < n; i++) {
      mui = exp(gsl_matrix_get(C, i, 0));
      gsl_vector_set(mu, i, mui);
    }

    gsl_matrix_free(B);
    gsl_matrix_free(C);
  }


  gsl_vector_free(mu);
  gsl_vector_free(z);
  return c;
}

