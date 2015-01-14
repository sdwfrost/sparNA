
/*
 *
 *	splines.h
 *  
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 06/18/14 14:11:41 CEST  
 *
 */

#include <gsl/gsl_statistics.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_bspline.h>

typedef struct
{
    size_t k;
    gsl_matrix *A; 
    gsl_matrix *dB; 
} gsl_bspline_deriv_workspace;

int
gsl_bspline_deriv_eval(const double x,
                       const size_t nderiv,
                       gsl_matrix *dB,
                       gsl_bspline_workspace *w,
                       gsl_bspline_deriv_workspace *dw);

gsl_bspline_deriv_workspace * gsl_bspline_deriv_alloc(const size_t k);
void gsl_bspline_deriv_free(gsl_bspline_deriv_workspace *w);


typedef struct{
  gsl_matrix *cov;
  gsl_vector *B;
  gsl_vector *c;
  gsl_bspline_workspace *bw;
} splinefit_t;

gsl_vector * irls(gsl_matrix *X,  gsl_vector *w, gsl_vector *y, gsl_vector *c, gsl_matrix *cov, 
    gsl_multifit_linear_workspace *mw,  Uint n, Uint ncoeffs);
double *quantiles(double *x, Uint n, double* q, Uint k);
gsl_vector *histsplineglm(double *myx, double *myy, int myn, double* q, Uint myqn, 
    splinefit_t *fit) ;
void destructSplineFit (splinefit_t* fit);
