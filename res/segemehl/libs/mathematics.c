
/*
 * mathematics.c
 * implemtation of various mathematical functions
 *
 * @author Steve Hoffmann
 * @date Wed 22 Nov 2006
 *
 *  SVN
 *  Revision of last commit: $Rev: 54 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-09-10 22:13:30 +0200 (Wed, 10 Sep 2008) $
 *
 *  Id: $Id: mathematics.c 54 2008-09-10 20:13:30Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/trunk/libs/mathematics.c $
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

int* intrev(int* n, Uint len){
  int end = len-1;
  int start = 0;

  while (start<end) {
    n[start] ^= n[end];
    n[end] ^= n[start];
    n[start] ^= n[end];
    start++;
    end--;
  }
  return n;
}

 void *initArray(void *space, int size, size_t datatype) {
	void *ptr=NULL;

	/*dirty trick: sizeof(char) == 1*/
	ptr = ALLOCMEMORY(space, ptr, char, size*datatype);
	return ptr;
 }


void appendvector(void *space, vector_t *v, vectorelem elem) { 

  	 v->elements = (vectorelem*) ALLOCMEMORY(space, v->elements, vectorelem, (v->length+1));
	 v->elements[v->length]=elem;
	 v->length++;
}




/*--------------------------------- mindist ----------------------------------
 *    
 * @brief expects a sorted vector to find the minimum distance between 
 * vec[i] and vec[j]
 * @author Steve Hoffmann 
 *   
 */

Uint
minvecdist(void *space, vector_t *vec, Uint i, Uint j) {
  Uint k, 
       size_i,
       size_j;
  int range,
      dist = INT_MAX,
      l;
  vectorelem *e_i,
             *e_j;

  size_j = LENGTHVEC(&vec[j]);
  size_i = LENGTHVEC(&vec[i]);

  if (size_i == 0 || size_j == 0) 
    return 0;

  e_j = &vec[j].elements[0];
  for(k=0; k < size_j ; k++, e_j++) {
    e_i = &vec[i].elements[0];        
    for(l=0; l < size_i; l++, e_i++) {
      range = abs((int)*e_j - (int)*e_i);
      if (range < dist) {
        dist = range;
      }
    }
  }

  return dist;
}


Uint
minvecdist2(void *space, vector_t *vec1, vector_t *vec2, Uint *which) {
  Uint k, 
       size_i,
       size_j;
  int range,
      dist = INT_MAX,
      l;
  vectorelem *e_i,
             *e_j;

  size_j = LENGTHVEC(vec2);
  size_i = LENGTHVEC(vec1);

  if (size_i == 0 || size_j == 0) 
    return 0;

  e_j = &vec2->elements[0];
  for(k=0; k < size_j ; k++, e_j++) {
    e_i = &vec1->elements[0];        
    for(l=0; l < size_i; l++, e_i++) {
      range = abs((int)*e_j - (int)*e_i);
      if (range < dist) {
        dist = range;
        *which = l;
      }
    }
  }

  return dist;
}



void dumpMatrix_int(int *M, int m, int n) {
	int i,j;

	for (i=0; i < m; i++) {
		for (j=0; j < n; j++){
			printf("%d ", MATRIX2D(M,n,i,j));
		}
			printf("\n");
	}
 }

Uint uarraymax(Uint *arr, Uint l) {
	Uint i;
	Uint max =0;
	
  	for(i=0; i < l; i++) {
		if (arr[i]>arr[max]) max=i;
	}
	
	return max;
}

Uint uarraysecond(Uint *arr, Uint l, Uint max) {
	Uint i;
	Uint second =0;
	
  	for(i=0; i < l; i++) {
		if (arr[i]>arr[second] && i!=max) 
          second=i;
	}
	
	return second;
}

int arraymax(int *arr, int l) {
	int i;
	int max =0;
	
  	for(i=0; i < l; i++) {
		if (arr[i]>arr[max]) max=i;
	}
	
	return max;
}

void dumpMatrix_Uint(Uint *M, Uint m, Uint n) {
	Uint i,j;

	for (i=0; i < m; i++) {
		for (j=0; j < n; j++){
			printf("%d ", MATRIX2D(M,n,i,j));
		}
			printf("\n");
	}
 }


void dumpMatrix_dbl(double *M, Uint m, Uint n) {
	Uint i,j;

	for (i=0; i < m; i++) {
		for (j=0; j < n; j++){
			printf("%f ", MATRIX2D(M,n,i,j));
		}
			printf("\n");
	}
 }


 void dumpMatrix3D_int(int *M, int m, int n, int l) {
	int i,j,k;

	for (i=0; i < m; i++) {
		for (j=0; j < n; j++){
			for (k=0; k < l; k++) {
				printf("%d ", MATRIX3D(M,n,l,i,j,k));
			}
		printf(";");
		}
	printf("\n");
	}
 }

void dumpVector(vector_t *v) {

	int i;
	for (i=0; i < v->length; i++) {
		printf("%d ", v->elements[i]);
	}

	printf("\n");
}


void destructVector(void *space, vector_t *v) {
	
    if (v!=NULL) {
    	if (v->elements) FREEMEMORY(space, v->elements);
		FREEMEMORY(space, v);
	}
}

void reverseVector(Uint a, Uint b, vector_t *v) {
	Uint i;
	
	for (i=0; i < (b-a); i++) {
		SWAPVEC(a+i,b-i,v);
	}
}

int nextPermutation(vector_t *v) {
	Uint i,j; 
	vectorelem *e=v->elements;

	for (i=(v->length)-1; i > 0; i--)
		if(e[i-1]<=e[i]) break;
	
	if (i==0) return 0;

	for (j=i+1; j < (Uint) v->length; j++ )
		if(e[i-1]>=e[j]) break;
	
	SWAPVEC(i-1, j-1, v);
	REVERSEVEC(i, (v->length)-1, v);

	return 1;
}

/*----------------------------------- norm -----------------------------------
 *    
 * @brief normalize
 * @author Steve Hoffmann 
 *   
 */
 
void
normalize (double *a, Uint n)
{
  Uint i;
  double sum = 0;

  for(i=0; i < n; i++) {
    sum += a[i];
  }
	  
  for(i=0; i < n; i++) {
    a[i] /= sum;
  }
 
  return ;
}



/*----------------------------------- gcd ------------------------------------
 *    
 * calculate the greatest common divisor of two integer values
 * 
 */
 
int
gcd (int a, int b)
{
    int val;

	b = abs(b);
	
	if (b > a)
	  val=a, a=b, b=val;

	while (b != 0) {
		val = a%b;
		a = b;
		b = val;
	}
	
	return a;
}


/*---------------------------------- power -----------------------------------
 *    
 * the power may be with you! 
 * 
 */
 
double
power (double x, int n)
{
  	double y = 1.;

	if(n==0)
	  return 1;
	if(x==0) {
		if(n < 0) {
			return MAX_DOUBLE;
		}
		return 0;
	}

	if (n < 0) {
		x = 1./x;
		n = -n;
	}

	while(n > 0) {
		if (n & 1) {
			y *= x;
		}
		n /= 2;
		x *= x;
	}

	return y;
}




/*----------------------------------- fak ------------------------------------
 *    
 * @brief get the factorial (works only for n <=10!!)
 * @author Steve Hoffmann 
 *   
 */
 

Uint fak(Uint n) {
  Uint i,x=n;
  
  for(i=x-1; i > 0; i--) {
  	x *= i;
  }

  return x;
}


/*--------------------------------- uniroot ----------------------------------
 *    
 * getting the zero-root of a given function
 * 
 * according to G. Forsythe, M. Malcom et al.
 * Computer methods for mathematical computations, 1980
 * 
 */
 
double
uniroot (double start, double end, double (*f)(double, void*), 
    double tolx, void* info)
{	
  	double a, b, c;
	double fa, fb, fc;
	double prev;
	double currenttol;
	double p, q, new_step;
	double cb, t1, t2;

	a = start; b= end; fa = (*f)(a,info); fb=(*f)(b,info);
	c = a; fc = fa;
	
	if ((fa > (double) 0 && fb > (double) 0) 
        || (fa < (double)0 && fb < (double)0)) {
		printf("mooep!\n");	
	  /*return 0;*/
	} 
	
	while(1) {

	  	prev = b-a;
		
		if (fabs(fc) < fabs(fb)) {
			a=b; b=c; c=a;
			fa=fb; fb=fc; fc=fa;
		}
		currenttol = 2 * FLT_EPSILON * fabs(b) + tolx/2;
		new_step = (c-b)/2;
		if (fabs(new_step) <= currenttol || fb == (double)0) {
			return b;
		}
		
		if ( fabs(prev) >= currenttol && fabs(fa) > fabs(fb) ) {
			cb = c-b;
			if(a==c) {
				t1 = fb/fa;
				p = cb*t1;
				q = 1.0 - t1;
			} else {
				q = fa/fc;
				t1 = fb/fc;
				t2 = fb/fa;
				p = t2 * ( cb * q * (q-t1) - (b-a)*(t1-1.0) );
				q = (q-1.0) * (t1-1.0) * (t2-1.0);	
			}
			if ( p > (double)0) {
				q = -q;
			} else {
				p = -p;
			}

			if(p < (0.75 * cb * q - fabs(currenttol*q)/2) 
				&& p < fabs(prev * q/2) ) {
				new_step = p/q;
			}
		}

		if (fabs(new_step) < currenttol ) {
			if(new_step > (double)0) {
				new_step = currenttol;
			} else {
				new_step = -currenttol;
			}
		}
		
		a=b; fa=fb;
		b+= new_step;
		fb = (*f)(b,info);
		if( (fb>0 && fc>0) || (fb < 0 && fc < 0) ) {
			c=a; fc=fa;
		}	
	}
	
	return 0;
}


/*---------------------------------- coldel ----------------------------------
 *    
 * @brief delete column for matrix
 * @author Steve Hoffmann 
 *   
 */
 
double*
coldel (void *space, double *a, Uint m, Uint n, Uint d) {
	
	double *t;
	Uint	i,
			j=-1,
			k=0,
			l=0;

  t = (double*) INITMATRIX2D(space, m, (n-1), sizeof(double));

  for(i=0; i < m*n; i++) {
	if(i % n == 0) { 
	  j++; k=0; l=0;
	} 	
	if(k++ != d) {
	  MATRIX2D(t, n-1, j, l++) = a[i];
	}
  }
	
  FREEMEMORY(space, a);
  return t;
}


/*---------------------------------- rowdel ----------------------------------
 *    
 * @brief delete row from matrix
 * @author Steve Hoffmann 
 *   
 */
 
double*
rowdel (void *space, double *a, Uint m, Uint n, Uint d) {
	
	double *t;
	Uint	i,
			j=-1,
			k=0,
			l=-1;

  t = (double*) INITMATRIX2D(space, (n-1), m, sizeof(double));

  for(i=0; i < m*n; i++) {
	if(i % n == 0) { 
	  j++; k=0;
	  l = (j != d) ? l+1 : l;
	} 	
	if(j != d) {
	  MATRIX2D(t, n, l, k++) = a[i];
	}
  }

  FREEMEMORY(space, a);
  return t;
}


/*---------------------------------- xprod -----------------------------------
 *    
 * @brief calculate the cross product of two vectors
 * @author Steve Hoffmann 
 *   
 */
 
double*
xprod (void *space, double* x, Uint m, double *y, Uint n) {
	double *p;
	Uint 	i,	
			j;

	p = (double*) INITMATRIX2D(space, m, n, sizeof(double));

	for (i=0; i < m; i++) {
		for(j=0; j < n; j++) {
			MATRIX2D(p, n, i, j) = x[i]*y[i];
		}
	}
	return p;
}


/*-------------------------------- transpose ---------------------------------
 *    
 * @brief transpose a matrix $a$ of dimensions $m x n$
 * @author Steve Hoffmann 
 *   
 */

double*
transpose (void* space, double *a, Uint m, Uint n) {
  double *t;
  Uint	i,
        j=-1,
        k=0;

  t = (double*) INITMATRIX2D(space, n, m, sizeof(double));

  for(i=0; i < m*n; i++) {
    if(i % n == 0) { j++; k=0;} 	
    MATRIX2D(t, m, k, j) = a[i];
    k++;
  }

  FREEMEMORY(space, a);
  return t;
}



/*--------------------------------- simpson ----------------------------------
 *    
 * @brief implementation of the simpson algorithm to determine the integral
 * of $f(x)$ in the interval [$a$,$b$]. sdiv denotes number of subdivisions
 * @author Steve Hoffmann 
 *   
 */

  double
simpson( double a, double b, int sdiv, 
    double (*f) (double, void*),
    void* info) 
{

  double 	k,
            sum1=0,
            sum2=0;
  int 	i;

  k = ((double) b-a)/((double)2*sdiv);

  for (i=1; i < sdiv; i++) {
    sum1+=f(a + k*2*i, info);
    sum2+=f(a + k*(2*i-1), info);
  }

  sum2+=f(a + k*(2*i-1), info);

  return ((double)(k/3)*(f(a, info)+f(b, info)+2*sum1+4*sum2));
}


/*-------------------------------- simpson1D ---------------------------------
 *    
 * @brief helper function for simpson2D
 * @author Steve Hoffmann 
 *   
 */

  double
simpson1D(double x, int sdiv, 
    double (*f) (double, double, void*),
    double (*c) (double, void*),
    double (*d) (double, void*),
    void* info) 
{

  double 	k,
            sum1=0,
            sum2=0,
            ca,
            da;
  int 	    i;

  ca = c(x, info);
  da = d(x, info);

  k = ((double) da-ca)/((double)2*sdiv);

  for (i=1; i < sdiv; i++) {
    sum1+=f(x, ca + k*2*i, info);
    sum2+=f(x, ca + k*(2*i-1), info);
  }

  sum2+=f(x, ca + k*(2*i-1), info);


  return ((double)(k/3)*(f(x, ca, info)+f(x, da, info)+2*sum1+4*sum2));
}


/*-------------------------------- simpson2D ---------------------------------
 *    
 * @brief calculates the 2-dim integral of function $f$ given the interval
 * [$a$,$b$] in the first and [$c(x)$,$d(x)$] in the second dimension
 * sdiv, sdiv2 denote the subdivisions in the first and second dimension
 * @author Steve Hoffmann 
 *   
 */

double 
simpson2D(double a, double b, int sdiv, int sdiv2, 
    double (*f) (double, double, void*), 
    double (*c) (double, void*),
    double (*d) (double, void*),
    void *info) {

  double 	h,	
  sum1=0,
  sum2=0;
  int		i;

  h = ((double)b-a)/((double)2*sdiv);

  for (i=1; i < sdiv; i++) {
    sum1 += simpson1D((a+h*2*(i)), sdiv2, f, c, d, info);
    sum2 += simpson1D((a+h*(2*i-1)), sdiv2, f, c, d, info);
  }

  sum2 += simpson1D((a+h*(2*i-1)), sdiv2, f, c, d, info);

  return ((double)(h/3) * (simpson1D(a, sdiv2, f, c, d, info) + 
        simpson1D(b, sdiv2, f, c, d, info) + 2*sum1 + 4*sum2 ));
}



/*--------------------------------- myMinor ----------------------------------
 *    
 * @brief helper function for the laplacian algorithm used in det()
 * @author Steve Hoffmann 
 *   
 */
double* 
myMinor(void *space, double* M, Uint m, Uint n, Uint i, Uint j) {

  double *t;

  t = (double*) ALLOCMEMORY(space, NULL, double, m*n);
  memmove(t, M, sizeof(double)*(m*n));

  t = rowdel(NULL, t, m, n, i);
  t = coldel(NULL, t, m-1, n, j);  

  return t;
}


/*----------------------------------- det ------------------------------------
 *    
 * @brief calculates the determinant of a square matrix m of size n x n using 
 * Laplacian algorithm (recursive implementation)
 * @author Steve Hoffmann 
 *   
 */

double
det(void *space, double *M, int n) {
  double sum=0,
         *t=NULL;
  int		j;

  if (n==1) {	
    return MATRIX2D(M, n, 0, 0);
  }

  for(j=0; j < n; j++) {
    t = myMinor(space, M, n, n, 0, j);
    sum += pow(-1.0, (j+2))*MATRIX2D(M,n,0,j)*det(space, t, n-1);
    FREEMEMORY(space, t);
  }

  return sum;
}

/*---------------------------------- invert ----------------------------------
 *    
 * @brief invert a matrix
 * @author Steve Hoffmann 
 *   
 */
 
double*
invert3D (void *space, double *M)
{
  
  double a, b, c, d, e, f, g, h, k;
  double detM;

  if((detM=det(space, M, 3)) !=0) {
    
    a = MATRIX2D(M, 3, 0, 0);
    b = MATRIX2D(M, 3, 0, 1);
    c = MATRIX2D(M, 3, 0, 2);
    d = MATRIX2D(M, 3, 1, 0);
    e = MATRIX2D(M, 3, 1, 1);
    f = MATRIX2D(M, 3, 1, 2);
    g = MATRIX2D(M, 3, 2, 0);
    h = MATRIX2D(M, 3, 2, 1);
    k = MATRIX2D(M, 3, 2, 2);


    MATRIX2D(M, 3, 0, 0) = (e*k-f*h)/detM;
    MATRIX2D(M, 3, 0, 1) = (f*g-d*k)/detM;
    MATRIX2D(M, 3, 0, 2) = (d*h-e*g)/detM;
    
    MATRIX2D(M, 3, 1, 0) = (c*h-b*k)/detM;
    MATRIX2D(M, 3, 1, 1) = (a*k-c*g)/detM;
    MATRIX2D(M, 3, 1, 2) = (g*b-a*h)/detM;
    
    MATRIX2D(M, 3, 2, 0) = (b*f-c*e)/detM;
    MATRIX2D(M, 3, 2, 1) = (c*d-a*f)/detM;
    MATRIX2D(M, 3, 2, 2) = (a*e-b*d)/detM;

    M = transpose(space, M, 3, 3);
    return M;
  }

   return NULL;
}



/*----------------------------------- add ------------------------------------
 *    
 * @brief componentwise addition of a to a vector of length m
 * @author Steve Hoffmann 
 *   
 */

double*
add(double *x, Uint m, double a) {
  Uint i;

  for(i=0; i < m; i++) {
    x[i] += a;
  }
  return x;
}


/*----------------------------------- mean -----------------------------------
 *    
 * @brief calculate the arithmetic mean for a vector of length m
 * @author Steve Hoffmann 
 *   
 */

double
mean (double *x, Uint m) {
  Uint i;
  double sum=0;

  for (i=0; i < m; i++) {
    sum += x[i];
  }

  return sum /= m; 
}

/*----------------------------------- mean -----------------------------------
 *    
 * @brief calculate the arithmetic mean for a vector of length m
 * @author Steve Hoffmann 
 *   
 */

double
mean_int (int *x, Uint m) {
  Uint i;
  double sum=0;

  for (i=0; i < m; i++) {
    sum += x[i];
  }

  return sum /= m; 
}



/*---------------------------------- scalar ----------------------------------
 *    
 * @brief calculate the scalar product of two vectors of length m
 * @author Steve Hoffmann 
 *   
 */

double
scalar (double* x, double *y, Uint m) {
  double  p=0;
  Uint 	i;

  for (i=0; i < m; i++) {
    p += x[i]*y[i];
  }
  return p;
}


/*----------------------------------- cov ------------------------------------
 *    
 * @brief get the covariance matrix (2x2) for two vectors of length m
 * @author Steve Hoffmann 
 *   
 */

double*
cov (void *space, double *x, double *y, Uint m) {
  double *c,
         xm,
         ym;

  c = (double*) INITMATRIX2D(space, 2, 2, sizeof(double));
  xm = mean(x, m);
  ym = mean(y, m);

  /*center*/
  add(x, m, (-1)*xm);
  add(y, m, (-1)*ym);

  MATRIX2D(c, 2, 0, 0) = (double) scalar(x,x,m)/(m-1);
  MATRIX2D(c, 2, 0, 1) = MATRIX2D(c, 2, 1, 0) = (double) scalar(x,y,m)/(m-1);
  MATRIX2D(c, 2, 1, 1) = (double) scalar(y,y,m)/(m-1);

  return c;
}




/*----------------------------------- var ------------------------------------
 *    
 * @brief get the sample variance
 * @author Steve Hoffmann 
 *   
 */
 
double
samplevar (double *x, double *p, double n)
{   
    int i;
    double m, r, sum=0;

    m=mean(x, n);
    for (i=0; i < n; i++) {
      r = x[i]-m;
      sum += (r*r)*p[i];
    }

	return sum/n;
}


/*----------------------------------- var ------------------------------------
 *    
 * @brief get the variance
 * @author Steve Hoffmann 
 *   
 */
 
double
var_int (int *x, Uint n)
{   
    int i;
    double m, r, sum=0;

    m=mean_int(x, n);
    for (i=0; i < n; i++) {
      r = x[i]-m;
      sum += (r*r);
    }

	return sum/n;
}

/*--------------------------------- poisson ----------------------------------
 *    
 * @brief the <.><< distribution
 * @author Steve Hoffmann 
 *   
 */
 
double
poisson(double lambda, double x) {
  assert(x >= 0);
  return (pow(lambda,x)/tgamma(x+1))*exp(-lambda);
}


/*-------------------------------- logpoisson --------------------------------
 *    
 * @brief poisson in log space
 * @author Steve Hoffmann 
 *   
 */
 
double
logpoisson (double loglambda, double logx)
{
  //assert(logx >= 0);
  return (exp(logx)*loglambda) - log(tgamma(exp(logx)+1)) - exp(loglambda);   
	
}

/*----------------------------------- var ------------------------------------
 *    
 * @brief get the variance
 * @author Steve Hoffmann 
 *   
 */
 
double
var (double *x, Uint n)
{   
    int i;
    double m, r, sum=0;

    m=mean(x, n);
    for (i=0; i < n; i++) {
      r = x[i]-m;
      sum += (r*r);
    }

	return sum/n;
}


/*---------------------------------- stddev ----------------------------------
 *    
 * @brief get the standard deviation
 * @author Steve Hoffmann 
 *   
 */
 
double
stddev (double *x, double n)
{
	
    return sqrt(var(x, n));
}


/*----------------------------------- rho ------------------------------------
 *    
 * @brief calculate correlation $\rho$ for two vectors of length m
 * @author Steve Hoffmann 
 *   
 */

double
rho (void *space, double *x, double *y, Uint m) {
  double *cv;

  cv = cov(space, x, y, m); 
  return (MATRIX2D(cv, 2, 0, 1)/sqrt(MATRIX2D(cv, 2, 0, 0)*MATRIX2D(cv, 2, 1, 1)));
}

/*-------------------------------- univarnorm --------------------------------
 *    
 * @brief pdf gaussian
 * @author Steve Hoffmann 
 *   
 */
 
double
univarnorm (double x, double mu, double sd)
{
    double d = (x-mu);
    double sdsq = sd * sd;

	return exp(-0.5*(d*d/sdsq))/sqrt(2*M_PI*sdsq);
}


/*------------------------------ univarnormcdf -------------------------------
 *    
 * @brief cdf gaussian
 * @author Steve Hoffmann 
 *   
 */
 
double
univarnormcdf (double x, double mu, double sd)
{
	return 0.5 * (1+erf((x-mu)/(sd*M_SQRT2)));
}


/*------------------------------ randunivarnorm ------------------------------
 *    
 * @brief algorithm adapted from Dr. Everett (Skip) Carter, Jr.
 * @author Steve Hoffmann 
 *   
 */
 
double
randunivarnorm (double mu, double sd)
{
  double x1, x2, w, y2;// y1;
  
  do{
    x1 = 2.0 * (((double)rand())/((double)RAND_MAX)) - 1.0;
    x2 = 2.0 * (((double)rand())/((double)RAND_MAX)) - 1.0;
    w = x1 * x1 + x2 * x2; 
  } while (w >= 1.0 || w == .0);

  w = sqrt((-2.0*log(w))/w);
  // not used: y1 = x1 *w;
  y2 = x2 *w;

  return y2*sd+mu;
}

/*-------------------------------- bivarcond ---------------------------------
 *    
 * @brief conditional bivar. norm. distrib. f(y|x) given location parameter
 * $mu1$, $mu2$ and covariance matrix $cv$ of size (2x2)
 * @author Steve Hoffmann 
 *   
 */

double
bivarcond(double x, double y, double mu1, double mu2, double *cv) {
  double rho,
         s1,
         s1sq,
         s2,
         s2sq,
         m,
         e;

  s1sq = MATRIX2D(cv, 2, 0, 0);
  s2sq = MATRIX2D(cv, 2, 1, 1);
  s1 = sqrt(s1sq);
  s2 = sqrt(s2sq);
  rho = MATRIX2D(cv, 2, 0, 1)/sqrt(s1sq*s2sq);

  m  = 1/sqrt((2*M_PI*s2sq*(1-(rho*rho))));
  e = (y-mu2-rho*(s2/s1)*(x-mu1));
  e *= e;
  e /= s2sq*(1-(rho*rho));

  return(m*exp(-0.5*e));
}


/*-------------------------------- bivarnorm ---------------------------------
 *    
 * @brief bivariate normal distribution f(x,y) given location parameter
 * $mu1$, $mu2$ and covariance matrix $cv$ of size (2x2)
 * @author Steve Hoffmann 
 *   
 */


double
bivarnorm(double x, double y, double mu1, double mu2, double* cv) {
  double rho,
         s1,
         s1sq,
         s2,
         s2sq,
         m,
         e1,
         e2;

  s1sq = MATRIX2D(cv, 2, 0, 0);
  s2sq = MATRIX2D(cv, 2, 1, 1);
  s1 = sqrt(s1sq);
  s2 = sqrt(s2sq);
  rho = MATRIX2D(cv, 2, 0, 1)/sqrt(s1sq*s2sq);

  m = 1/(2*M_PI*s1*s2*sqrt(1-(rho*rho)));

  e1 = (-1)/(2*(1-(rho*rho)));
  e2 = ((x-mu1)*(x-mu1))/s1sq  
    - (2*rho*(x-mu1)*(y-mu2))/(s1*s2)
        + ((y-mu2)*(y-mu2))/s2sq; 

  return m*exp(e1*e2);
}

/*------------------------------ multivarnorm -------------------------------
 *    
 * @brief n-dimensional gaussian probability density function wo correlation term
 * i.e. orthogonal variation!
 * @author Steve Hoffmann 
 *   
 */
 
double
multivarnorm (double *pt, double *mu, double *sd, Uint n)
{
  Uint i;
  double det = 1; 
  double exponent = 0;
  double perturb = 0;

  for(i=0; i < n; i++) {
    exponent += pow((pt[i]-mu[i]),2) * (1/pow(MAX(sd[i],perturb),2));
    det *= pow(sd[i],2);
  }

  //fprintf(stderr, "exponent: %f, det: %f, exp(expo):%f, norm:%f\n", exponent, det,
  //    exp(exponent*-0.5), sqrt(pow((2*M_PI),n)*det));

  return exp(exponent*-0.5)/sqrt(pow((2*M_PI),n)*det);	
}


/*this is ncbi intellectual property*/

double BLAST_Expm1(double x)
{
  double	absx = ABS(x);

  if (absx > .33)
	return exp(x) - 1.;

  if (absx < 1.e-16)
	return x;

  return x * (1. + x *
	  (1./2. + x * 
	   (1./6. + x *
		(1./24. + x * 
		 (1./120. + x *
		  (1./720. + x * 
		   (1./5040. + x *
			(1./40320. + x * 
			 (1./362880. + x *
			  (1./3628800. + x * 
			   (1./39916800. + x *
				(1./479001600. + 
				 x/6227020800.))))))))))));
}


/*---------------------------------- log10 -----------------------------------
 *    
 * @brief logarithm to the base 10
 * @author Steve Hoffmann 
 *   
 */
 
double
log10(double x) {
  return (log(x)/log(10));
}

/*----------------------------------- log2 -----------------------------------
 *    
 * @brief logarithm to the base 2
 * @author Steve Hoffmann 
 *   
 */
 

double
log2(double x) {
  return (log(x)/log(2));
}


/*--------------------------------- log10add ---------------------------------
 *    
 * @brief addition in log10 space
 * @author Steve Hoffmann 
 *   
 */
 
double 
log10add(double a, double b) {
  double max, min;

  if(a > b) {
    if (b == log10(0)) return a;
    else {
      max=a;
      min=b;
    }
  } else {
    if (a == log10(0)) return b;
    else {
      max=b;
      min=a;
    }
  }
  return max + log10(1+pow(10.,min-max));  
}

/*---------------------------------- logadd ----------------------------------
 *    
 * @brief addition in logarithmus naturalis space
 * @author Steve Hoffmann 
 *   
 */
 

double 
logadd(double a, double b) {
  double max, min;

  if(a > b) {
    if (b == log(0)) return a;
    else {
      max=a;
      min=b;
    }
  } else {
    if (a == log(0)) return b;
    else {
      max=b;
      min=a;
    }
  }
  return max + log(1+exp(min-max));  
}

/*-------------------------------- seqentropy --------------------------------
 *    
 * @brief zero order sequence entropy
 * @author Steve Hoffmann 
 *   
 */
 

double
shannonentropy(void *space, char *seq, Uint len, Uint asize, Uint *encodetab) {
  //CHANGED: count should start with 0
  //--> possibly just use len instead of norm 
  //Uint i, norm=1;
  Uint i, norm=0;
  double *p, H=0;


  p = ALLOCMEMORY(space, NULL, double, asize);
  memset(p, 0, sizeof(double)*asize);

  for(i=0; i < len; i++) {
    p[encodetab[(Uint)seq[i]]]++;
    norm++;
  }

  for(i=0; i < asize; i++) {
    p[i]/=norm;
  }

  for(i=0; i < asize; i++) {
    if(p[i] > 0)
    H += p[i] * log2(p[i]); 
  }

  FREEMEMORY(space, p);
  return -1*H;
}

/*-------------------------------- smoothavg ---------------------------------
 *    
 * @brief smoothing average
 * @author Steve Hoffmann 
 *   
 */
 
double*
smoothavg (double *y, Uint n, Uint r)
{
  double *ys;
  Uint i, j;

  ys = ALLOCMEMORY(space, NULL, double, n);

  for(i=r; i < n; i++) {
    for(j=0; j < 2*r; j++) {
      ys[r+j] += y[i-r+j] / (2*r+1);
    }
  }

  return ys;
}

/*----------------------------------- gmm ------------------------------------
 *    
 * @brief gaussian mixture model for m n-dimensional data points and g gaussians
 * @author Steve Hoffmann 
 *   
 */
 
double
gmm(void *space, double *pt, Uint m, Uint n, 
    double *mu, double *sd, double *w, Uint g, Uint maxiter) {

  double *ms, *mu_, *sd_, no, dd, ll, oll, epsilon=0.0000000000000001;
  Uint i, j, k, l=0;

  ms = ALLOCMEMORY(space, NULL, double, m*g);
  mu_ = ALLOCMEMORY(space, NULL, double, n);
  sd_ = ALLOCMEMORY(space, NULL, double, n);
  memset(sd_, 0, sizeof(double)*n);
  memset(mu_, 0, sizeof(double)*n);
  memset(ms, 0, sizeof(double)*m*g);

  ll = 1;

  do {

    oll = ll;

    /*expectation step*/
    for(ll=0, i=0; i < m; i++) {
      for(j=0; j < g; j++) {
        MATRIX2D(ms, g, i, j) = w[j] * 
          multivarnorm(&MATRIX2D(pt, n, i, 0), &MATRIX2D(mu, n, j, 0),
              &MATRIX2D(sd, n, j, 0), n);
        ll += log(MATRIX2D(ms, g, i, j));
      }
      normalize(&MATRIX2D(ms, g, i, 0), g);
    }

    /*maximization step - weighted normalized mean*/
    for(j=0; j < g; j++) {
      memset(mu_, 0, sizeof(double)*n);

      for(i=0, no=0; i < m; i++) {
        no += MATRIX2D(ms, g, i, j);
        for(k=0; k < n; k++) {
          mu_[k] += MATRIX2D(pt, n, i, k) * MATRIX2D(ms, g, i, j);
        }
      }

      /*update mean*/   
      for(k=0; k < n; k++) { 
        mu_[k] /= no;
        MATRIX2D(mu, n, j, k) = mu_[k];
      }

      for(i=0; i < m; i++) {
        for(k=0; k < n; k++) {
          dd = MATRIX2D(pt, n, i, k) - mu_[k];
          sd_[k] += (dd * dd) * MATRIX2D(ms, g, i, j);
        }
      }

      /*update sd*/
      for(k=0; k < n; k++) { 
        sd_[k] /= no;
        MATRIX2D(sd, n, j, k) = sqrt(sd_[k]);
      }

      /*update weight*/
      w[j] = no/((double)m);

    }
   
  } while(l++ < maxiter && fabs((double)ll-oll) > -1.0*ll*epsilon);

   
  FREEMEMORY(space, mu_);
  FREEMEMORY(space, sd_);
  FREEMEMORY(space, ms);

  return ll;
}



/*--------------------------------- bl_RSSi ----------------------------------
 *    
 * @brief calculation of a RSS for a vector of length n from u to v for intercept
 * @author Steve Hoffmann 
 *   
 */
 
double*
bl_RSS (void *space, double *x, Uint n, Uint u, Uint v)
{
  Uint i;
  double *cum, *y;
  
  assert(v>u);
  
  cum = ALLOCMEMORY(space, NULL, double, v-u+1);
  y = ALLOCMEMORY(space, NULL, double, v-u+1);

  cum[0] = x[u];
  y[0] = x[u] * sqrt(2);

  for(i=u+1; i < v; i++) {
    cum[i-u] = x[i] + cum[i-u-1]; 
    y[i-u] = x[i] - cum[i-u]/((double)i-u+1);
    y[i-u] *= sqrt(1.0+(1.0/((double)i-u)));
  }
      
  cum[0] = 0;

  for(i=1; i < v-u+1; i++) {
    cum[i] = cum[i-1] + y[i]*y[i];
  }

  //diagonal
  cum[0] = 0;

  FREEMEMORY(space, y);
  return cum;
}


/*------------------------------- bl_RSSmatrix -------------------------------
 *    
 * @brief compute the triangular RSS matrix 
 * @author Steve Hoffmann 
 *   
 */
 
breakpoints_t*
bl_RSSmatrix (void *space, double *x, Uint n, Uint h, Uint noofbreaks)
{
  Uint i, j, minarg, *POS, m, resc;  //m breaks
  breakpoints_t *breaks;
  double *M, *RSS, *temp, minval;
  double *y;

  /*calculation of rss matrix*/
  M = ALLOCMEMORY(space, NULL, double, (n-h+1)*n);
  memset(M, 0, sizeof(double)*((n-h+1)*n));

  for(i=0; i < n-h+1; i++) {
    y = bl_RSS(space, x, n, i, n);
    memmove(&MATRIX2D(M, n, i, i), y, sizeof(double)*(n-i));
    FREEMEMORY(space, y);
  }
  
  /*DP in the rss matrix*/
  RSS = ALLOCMEMORY(space, NULL, double, (noofbreaks+1)*(n+1));
  POS = ALLOCMEMORY(space, NULL, Uint, (noofbreaks+1)*(n+1));
  memset(RSS, 0, (noofbreaks+1)*(n+1)*sizeof(double));
  memset(POS, 0, (noofbreaks+1)*(n+1)*sizeof(Uint));

  for(i=0; i < n; i++) {
    MATRIX2D(RSS, noofbreaks+1, i, 1) = MATRIX2D(M, n, 0, i);
    MATRIX2D(POS, noofbreaks+1, i, 1) = i;
  }

  for(m=2; m <= noofbreaks; m++) { 
    for(i=m*h-1; i < n-h; i++) {
      temp = ALLOCMEMORY(space, NULL, double, (i-h)-((m-1)*h-1)+1);
      /*Recursion*/
      for(j=(m-1)*h-1; j < i-h+1; j++) {
        temp[j-(((m-1)*h)-1)] = MATRIX2D(RSS, noofbreaks+1, j, m-1) + MATRIX2D(M, n, j+1, i);
      }
      /*find the minimum*/
      minarg = 0;
      minval = temp[0];
      for(j=1; j < (i-h)-((m-1)*h-1)+1; j++) {
        if(temp[j] < minval) {
          minarg = j;
          minval = temp[j];
        }
      }
      /*register*/
      MATRIX2D(RSS, noofbreaks+1, i, m) = temp[minarg];
      MATRIX2D(POS, noofbreaks+1, i, m) = (m-1)*h-1+minarg;

      FREEMEMORY(space, temp);
    }
  }

  
  breaks = ALLOCMEMORY(space, NULL, breakpoints_t, noofbreaks+1);
  breaks[0].noofbreaks = 0;
  breaks[0].RSS = MATRIX2D(RSS, noofbreaks+1, n-1, 1);
  breaks[0].LL = -0.5 * n * (log(breaks[0].RSS) + 1 - log(n) + log(2*M_PI));
  breaks[0].BIC = -2*breaks[0].LL + 2*log(n);
  breaks[0].breaks = NULL;


  for(m=noofbreaks; m >= 1; m--) { 
    /*calling breaks and backtrace*/
    temp = ALLOCMEMORY(space, NULL, double, n);

    for(i=h-1; i < n-h; i++) {
      if (MATRIX2D(RSS, noofbreaks+1, i, m) == 0.0 || 
          MATRIX2D(M, n, i+1, n-1) == 0.0) {
        temp[i]= -1.0;
      } else {  
        temp[i] = MATRIX2D(RSS, noofbreaks+1, i, m) + MATRIX2D(M, n, i+1, n-1);
      }
    }

    minarg = 0;
    minval = temp[h-1];
    for(i=h-1; i < n-h; i++) {
      if((temp[i] > 0.0 && temp[i] < minval) || minval == -1.0) {
        minarg = i;
        minval = temp[i];
      }
    }

    breaks[m].noofbreaks = m; 
    breaks[m].breaks = ALLOCMEMORY(space, NULL, Uint, m);
    breaks[m].RSS = temp[minarg];
    breaks[m].LL =  -0.5 * n * (log(temp[minarg]) + 1 - log(n) + log(2*M_PI));
    breaks[m].BIC =  -2*breaks[m].LL + (2*(m+1))*log(n);
    breaks[m].breaks[0] = minarg;


    for(i=m, j=1; i >= 2; i--, j++) {
       breaks[m].breaks[j] = MATRIX2D(POS, noofbreaks+1, breaks[m].breaks[j-1], i);
    }
 
    //reverse
    for(i=0; i < m/2; i++) { 
      resc = breaks[m].breaks[m-i-1];
      breaks[m].breaks[m-i-1] = breaks[m].breaks[i];
      breaks[m].breaks[i] = resc;
    }
    FREEMEMORY(space, temp);
  }
   
  FREEMEMORY(space, RSS);
  FREEMEMORY(space, M);
  
  return breaks;
}


/*---------------------------------- gevcdf ----------------------------------
 *    
 * @brief cumulative extreme value distribution
 * @author Steve Hoffmann 
 *   
 */
 
double
gevcdf (double x, double mu, double si, double xi)
{
  double t, r=.0;

    if(xi != 0.0) { 
      t = 1.0 + xi*((x-mu)/si);
      t = pow(t,(-1.0/xi));
    } else {
      t = exp(-(x-mu)/si);
    }

    r = exp(-t);
    return r;
}

/*---------------------------------- llgev -----------------------------------
 *    
 * @brief log likelihood of a gev Prescott Walden Formulation
 *        l(mu, sigma, kappa) = 
 *              -\Sum log(sigma) - (1-k) \Sum x_i - \Sum exp(-x_i)
 *        x_i = -(1/k) log[1 - k((y_i - mu)/sigma)]
 *        NOTE: k = -\xi !!!
 * @author Steve Hoffmann 
 *   
 */
 
double
gevll (double *data, Uint n, double mu, double sigma, double kappa)
{
  Uint i;
  double sum1 = .0, sum2 =.0, sum3=.0, x, z;

  sum1 =-1.0*log(sigma) * (double)n;

  for(i=0; i < n; i++) {
    z = 1.0-(kappa*((data[i]-mu)/sigma));
    if(z >= .0)
        x = -1.0/kappa * log(z);
    else
        x = 1/1e-5;
    sum2 += x;
    sum3 += exp(-1.0*x);
  }

  sum1 = sum1 - ((1.0-kappa)*sum2) - sum3;
  return sum1;
}


/*-------------------------------- gev1stdev ---------------------------------
 *    
 * @brief first derivative of gev location mu
 * @author Steve Hoffmann 
 *   
 */
 
double*
gevinvhessian(void *space, double *y, Uint n, double m, double s, double k, 
        double* dm, double* ds, double* dk, Uint *excl)

{
  //the likelihood is given by:
  //l(m,s,k) = \sum \log(s) - (1-k) \Sum x_i - \Sum exp (-x_i)
  //where x_i is given by
  //x_i = (-1/k) \log(1-k\frac{y-m}{s})
  //we denote q = -(1-k) \Sum x_i
  //we denote z = 1-k\frac{y-m}{s}
  //we denote r = exp(-1*x_i)

  Uint i, exclude=0;
  double x=0, ex=0, kmys, kmys2, logs, k2, k3, s2, my, my2, logkmys;
  double cdqm =.0, cdqs = .0, cdqk =.0;
  double cdqmm =.0, cdqss =.0, cdqkk =.0;
  double cdqms =.0, cdqmk =.0, cdqsk =.0;
  double cdrm =.0, cdrs = .0, cdrk =.0;
  double cdrmm =.0, cdrss =.0, cdrkk =.0;
  double cdrms =.0, cdrmk =.0, cdrsk =.0;
  double dxm =.0, dxs =.0, dxk =.0; 
  double dxmm, dxss, dxkk, dxms, dxmk, dxsk;
  double *H;
  double fulldm = .0;

  logs = log(s);
  k2 = k*k;
  k3 = k*k*k;
  s2 = s*s;

  for(i = 0; i < n; i++) {

    my = m - y[i];
    my2 = my * my; 
    kmys = k*my+s;   //kmys = z * s;
    kmys2 = kmys*kmys;  
    logkmys = log(kmys);  

    if(kmys > (double).0) { 

      x = (-1.0/k)*log(1.0-(k*(1.0/s)*((y[i]-m))));
      ex = exp(-1.0*x);

      dxm = -1.0/(kmys);
      dxs = my/(s*kmys);
      dxk  = -(1.0/k2) * ((k*my)/kmys - logkmys + logs);

      dxmm = k/kmys2;
      dxss = -my*(k*my+2.0*s)/(s2*kmys2);
      dxkk = (1.0/k3) * (((k2*my2)/(kmys2)) + 2.0*((k*my)/kmys -logkmys + logs));

      //dxms = 1.0/((k*(-m+y[i])-s)*(k*(-m+y[i])-s));
      dxms = 1.0/(kmys2);
      dxmk = my/kmys2;
      dxsk = -1.0*my2/(s*kmys2);

      //sums of derivatives for q direct or product rule
      cdqm += dxm;
      cdqs += dxs;
      cdqk += -1.0*x + (1.0-k)*dxk; //product rule
      cdqmm += dxmm;
      cdqss += dxss;
      cdqkk += (-1.0* dxk) - dxk + (1.0-k)*dxkk; //2x product rule
      cdqms += dxms;
      cdqmk += (m+s-y[i])/kmys2; //direct derivative
      cdqsk += -(my*(m+s-y[i]))/(s*kmys2); //direct derivative

      //sums of derivatives for r using chain, product rule
      cdrm += ex * -1.0*dxm;
      cdrs += ex * -1.0*dxs;
      cdrk += ex * -1.0*dxk;
      cdrmm += ex * ((dxm*dxm) + -1.0*dxmm);         
      cdrss += ex * ((dxs*dxs) + -1.0*dxss); 
      cdrkk += ex * ((dxk*dxk) + -1.0*dxkk);
      cdrms += ex * ((dxm*dxs) + -1.0*dxms);
      cdrmk += ex * ((dxm*dxk) + -1.0*dxmk);
      cdrsk += ex * ((dxs*dxk) + -1.0*dxsk);       

    fulldm += -(exp((1/k)*log((k*m-k*y[i]+s)/s))+k-1) * (1/kmys);

    }  else {
        exclude++;
    }
    
  }

  cdqm *= -(k-1.0); cdqs *= -(k-1.0);
  cdqmm *= -(k-1.0); cdqss *= -(k-1.0); cdqms *= -(k-1.0);

  *dm = (- cdqm - cdrm);
  *ds = (-1.0*n*(1/s) - cdqs - cdrs);
  *dk = (-cdqk - cdrk);

  H = ALLOCMEMORY(space, NULL, double, 9);

  H[0] = (-cdqmm - cdrmm) * -1.0;
  H[4] = (n*(1.0/s2) - cdqss - cdrss) * -1.0;
  H[8] = (-cdqkk - cdrkk) * -1.0;
  H[1] = H[3] = (-cdqms - cdrms) * -1.0;
  H[2] = H[6] = (-cdqmk - cdrmk) * -1.0;
  H[5] = H[7] = (-cdqsk - cdrsk) * -1.0;

//  fprintf(stderr, "excluded: %d\n", exclude);
//  fprintf(stderr, "du:%f, da:%f, dg:%f, fulldu:%f\n", *dm, *ds, *dk, fulldm);
//  fprintf(stderr, "duu:%f, dua:%f, daa:%f, dug:%f, dag:%f, dgg:%f\n", H[0], H[1], H[4], H[2], H[5], H[8]);

  H = invert3D(space, H);

  if(H) {
//    fprintf(stderr, "du:%e, da:%e, dg:%e\n", *dm, *ds, *dk);
//    fprintf(stderr, "duu:%e, dua:%e, daa:%e, dug:%e, dag:%e, dgg:%e\n", H[0], H[3], H[4], H[6], H[7], H[8]);
  }

  *excl = exclude;
  return H;
}

/*--------------------------------- gevmle -----------------------------------
 *    
 * @brief fit a generalized extreme value distribution
 * @author Steve Hoffmann 
 *   
 */
 
double
gevmle(void *space, double *y, Uint n,  
    double *m, double *s, double *k, Uint maxiter, double ymin, double ymax) {

  Uint l=0, maxrcnt = 20, rcnt=0, exclude;
  double *H, oll=-1.0*DBL_MAX, ll;
  double vm=.0, vs=.0, vk=.0;
  double dm, ds, dk;
  double lm, ls, lk;
  double sm=0.5, ss=0.25, sk=0.02, acc=1e-8, rho=0.25, smax;
  double z;
  
  lm = *m;
  ls = *s;
  lk = *k;

  //adjust values to avoid failure of invHessian (ie. kmys <= 0)
  if(fabs(lk) < 1e-4) lk = 1e-4;
  if(ls < 0) ls = 1;
  
  if(lk <= 0) {
    if(ymin < lm) {
      z = ls/(ymin-lm);
      if(lk <= z) {
        lk = z + 1e-4;
        if(lk >= 0) lk = 0.5 * z;
      }
    }
  } else {
    if(ymax > lm) {
      z = ls/(ymax-lm);
      if(lk >= z) {
        lk = z - 1e-4;
        if(lk <= 0) lk = 0.5 *z;
      }
    }
  }


  do {
    
    H = gevinvhessian(space, y, n, lm, ls, lk, &dm, &ds, &dk, &exclude);
    ll = gevll(y, n, lm, ls, lk);
    if(exclude) ll *= 0.9;

   // fprintf(stderr, "-------------\n iter:%d m:%f, s:%f, k:%f -> ll:%.10f oll:%.10f\n", 
  //      l, lm, ls, lk, ll, oll);

    if(H && (H[0] >= 0 && H[4] >= 0 && H[8] >= 0) && ll > oll) {
      oll = ll;
      
   //   fprintf(stderr, "likelihood has increased\n");

      vm = H[0]*dm + H[3]*ds + H[6]*dk;
      vs = H[4]*ds + H[3]*dm + H[7]*dk;
      vk = H[8]*dk + H[6]*dm + H[7]*ds;
   
      FREEMEMORY(space, H);

      smax = (fabs(vm)/(sm*ls) > fabs(ds)/(ss*ls)) ? fabs(vm)/(sm*ls)  : fabs(ds)/(ss*ls);
      smax = (smax > fabs(vk)/sk) ? smax : fabs(vk)/sk;

      if(smax < 1) {  
        smax = 1/smax;
        vm *= smax;
        vs *= smax;
        vk *= smax;
      }

    } else {

      if(ll <= oll) {
   //     fprintf(stderr, "last move was wrong");
        lm -= vm;
        ls -= vs;
        lk -= vk;
      } else { 
   //     fprintf(stderr, "steepest ascent\n");
        double tmp1 = 1e37;
        double tmp2 = 1e37;
        double tmp3 = 1e37;

        if(dm != 0.0) tmp1 = sm/(lm*fabs(dm));
        if(ds != 0.0) tmp2 = ss/(ls*fabs(ds));
        if(dk != 0.0) tmp3 = sk/(fabs(dk));
        smax = MAX3(tmp1, tmp2, tmp3);

        vm = dm * lm * lm * smax;
        vs = ds * ls * ls * smax;
        vk = dk * smax;
      }
    }

     rcnt = 0;
     do {  
   //   fprintf(stderr, "rescale (vm:%f, vs:%f, vk:%f)\n", vm, vs, vk);
      vm *= rho;
      vs *= rho;
      vk *= rho;   
    }  while ((1-(lk+vk)*(ymin-(lm+vm))/(ls+vs) < 0 || 
               (1-(lk+vk)*(ymax-(lm+vm))/(ls+vs) < 0)) && rcnt++ <maxrcnt);
    
    lm += vm;
    ls += vs;
    lk += vk;
   
 //   fprintf(stderr, "vm:%.8f, vs:%.8f, vk:%.8f\n", vm, vs, vk);
 //   fprintf(stderr, "lm:%.8f, ls:%.8f, lk:%.8f\n", lm, ls, lk);

  } while(l++ < maxiter && fabs(vm)>acc*ls && fabs(vs)>acc*ls && fabs(vk)>acc);

  *m = lm;
  *s = ls;
  *k = lk;

  return .0;
}



/*---------------------------------- gevvar ----------------------------------
 *    
 * @brief get the variance for gev
 * @author Steve Hoffmann 
 *   
 */
 
double
gevvar (double mu, double s, double k)
{
  if(mu != .0 && k < 1) {
    return (s*s)*(exp(lgamma(1.0-2.0*k)) - exp(lgamma(1.0-k))*exp(lgamma(1.0-k)))/(k*k);
  }

  if(mu == .0) {
    return s*s*((M_PI*M_PI)/k);
  }

  return -1*log(0);
}


/*---------------------------------- stirling ----------------------------------
 *    
 * @brief stirlings's approximation for the gamma function (numerical recipies)
 * @author Steve Hoffmann 
 *   
 */
 
double
gammaln (double x)
{
    Uint i;
    double y,tmp;
    double ser =  1.000000000190015;
    double coef[6] = { 76.18009172947146, -86.50532032941677, 
                       24.01409824083091, -1.231739572450155,  
                       0.1208650973866179e-2, -0.5395239384953e-5}; 
    y = x;
    
    for(i=0; i < 6; i++) {
      ser += coef[i]/++y;
    }

    tmp = (x+5.5)-(x+0.5)*log(x+5.5);
    return -tmp+log(2.5066282746310005*ser/x);
}



/*-------------------------------- gevmoment ---------------------------------
 *    
 * @brief L moment estimators (Hosking 1985)for a generalized extreme value 
 *        distribution using probability weighted moments PWM 
 *        (Wallis 1980, 1985; cf. Martins Stedinger 2000); 
 *        data is ordered x1 < x2
 * @author Steve Hoffmann 
 *   
 */
 
void
gevLmoment (double *data, Uint n, double *m, double *s, double *k)
{
  
    Uint i;
    double beta[3], lambda[4], c, kappa, alpha, xi, q, g;

    beta[0] = .0;
    beta[1] = .0;
    beta[2] = .0;

    for(i = 1 ; i <= n; i++) {
      beta[0] += data[i-1];
      q = (((double)(i-1))/(n-1));
      beta[1] += q * data[i-1];
      beta[2] += (q*(((double)(i-2))/(n-2))) * data[i-1];
    }
    
    beta[0] /= n; 
    beta[1] /= n; 
    beta[2] /= n;
    
    lambda[1] = beta[0];
    lambda[2] = 2.0*beta[1] - beta[0];
    lambda[3] = 6.0*beta[2] - 6.0*beta[1] + beta[0];

    //calculation of theta3 = lambda3 / lambda2
    c = 2.0/(3.0+(lambda[3]/lambda[2])) - log(2)/log(3);
    kappa = 7.8590*c + 2.9554*(c*c);
    g = exp(gammaln(1.0+kappa));
    alpha = (lambda[2]*kappa)/((1-pow(2.0,-1.0*kappa)) * g); 
    xi = lambda[1] - alpha*(1.0-exp(gammaln(1.0+kappa)))/kappa;

    *m = xi;
    *s = alpha;
    *k = kappa;

	return ;
}


/*------------------------------- dumpmatrix2D -------------------------------
 *    
 * @brief dump a 2-Dmatrix
 * @author Steve Hoffmann 
 *   
 */
 
void
dumprecursionmatrix2D (FILE *dev, int **M, char**B, char **K, Uint m, Uint n, 
    PairUint *mark)
{
  Uint i,j;


  fprintf(dev, " \t");
  for(j=0; j < n-1; j++) {
    fprintf(dev, "  %d    \t", j);
  }
  fprintf(dev, "\n");

  for(i=0; i < m; i++) { 
    fprintf(dev, "%d\t", i);
    for(j=0; j < n-1; j++) {
      if(mark->a == i && mark->b ==j)
        fprintf(dev, "^");
      if(B[i][j] && K[i][j])
        fprintf(dev, "-*%u*-\t", M[i][j]);
      else if(B[i][j])
        fprintf(dev, " *%u* \t", M[i][j]);
      else if(K[i][j])
        fprintf(dev, "- %u -\t", M[i][j]);
      else
        fprintf(dev, "  %u  \t", M[i][j]);
    }
      if(B[i][j] && K[i][j])
        fprintf(dev, "-*%u*-\n", M[i][j]);
      else if(B[i][j])
        fprintf(dev, " *%u* \n", M[i][j]);
      else if(K[i][j])
        fprintf(dev, "- %u -\n", M[i][j]);
      else
        fprintf(dev, "  %u  \n", M[i][j]);
  }
	
  return ;
}


/*-------------------------------- upperbound --------------------------------
 *    
 * @brief get the first element in a sorted list greater x
 * @author Steve Hoffmann 
 *   
 */
 
Uint
upperbound (double *l, Uint lo, Uint hi, double x)
{
  Uint cur=lo, begin = lo, size=hi-lo, step =0; 

  while(size > 0) {
    step=size/2; 
    cur = begin+step;
    if(x >= l[cur]) {
      begin = ++cur;
      size -= (step+1);
    } else {
      size = step;
    }
  }
	
  return begin;
}


/*-------------------------------- lowerbound --------------------------------
 *    
 * @brief get the first element in a sorted list that equals x
 * @author Steve Hoffmann 
 *   
 */
 
Uint
lowerbound (double *l, Uint lo, Uint hi, double x)
{
  Uint cur=lo, begin = lo, size=hi-lo, step =0; //rename count to size

  while(size > 0) {
    step=size/2; 
    cur = begin+step;
    if(x > l[cur]) {
      begin = ++cur;
      size -= (step+1);
    } else {
      size = step;
    }
  }
	
  return begin;
}


/*----------------------------------- ecdf -----------------------------------
 *    
 * @brief calculate the probability for a given ecdf
 * @author Steve Hoffmann 
 *   
 */
 
double
ecdf (double x, ecdf_t *e)
{
  Uint pos = upperbound(e->l, 0, e->n-1, x);
  return ((double)pos)/((double)e->n);
}


/*-------------------------------- ecdf_init ---------------------------------
 *    
 * @brief initalize the ecdf
 * @author Steve Hoffmann 
 *   
 */
 
ecdf_t*
ecdf_init (double *x, Uint n)
{
  ecdf_t* e;
  double *l;

  l = ALLOCMEMORY(NULL, NULL, double, n);
  memmove(l, x, sizeof(double)*n);
  e = ALLOCMEMORY(NULL, NULL, ecdf_t, 1);
  
  qsort(l, n, sizeof(double), cmp_dbl_qsort);

  e->l = l;
  e->n = n;

  return e;
}


/*------------------------------ ecdf_destruct -------------------------------
 *    
 * @brief destruct ecdf
 * @author Steve Hoffmann 
 *   
 */
 
void
ecdf_destruct (ecdf_t *e)
{
  FREEMEMORY(NULL, e->l);
  return ;
}

/*--------------------------------- incbeta ---------------------------------
 *    
 * @brief 
 * @author Steve Hoffmann 
 *   
 */
 
double
incbeta (double x, double a, double b, char lower)
{
 
  if(x <= 0) return .0;
  if(x >= 1) return 1.;
  double x1 = 0.5 - x + 0.5, w, wc;
  long int ierr;

  bratio_(&a, &b, &x, &x1, &w, &wc, &ierr);
  return lower ? w : wc;
}

/*---------------------------------- pbinom ----------------------------------
 *    
 * @brief binominal distribution
 * @author Steve Hoffmann 
 *   
 */
 
double
pbinom (double x, double n, double p, char lower)
{
  if(n<0 || p<0 || p>1) return 0;
  if (x < 0) return 0.0;
  x = floor(x + 1e-7);
  if (n <= x) return 1.0;
  return incbeta(p, x + 1, n - x, !lower); 
}


/*-------------------------------- quantiles ---------------------------------
 *    
 * @brief get k quantiles from x
 * @author Steve Hoffmann 
 *   
 */
 
double*
quantiles (double *x, Uint n, double *p, Uint k)
{
  Uint i;
  double *q;

  //sort the array
  qsort(x, n, sizeof(double), cmp_dbl_qsort);

  q = ALLOCMEMORY(NULL, NULL, double, k);
    
  for(i=0; i < k; i++) {
    double h = ((double)n-1)*p[i] + 1.0;
    Uint hlo = floor(h); 
    q[i] = x[hlo-1] + (h - ((double)hlo))*(x[hlo] - x[hlo-1]);
  }

  return q;
}


/*------------------------- choleskyTriDiagArrowFact -------------------------
 *    
 * @brief cholesky tri-diagnal arrow factorization for splines
 * @author Steve Hoffmann 
 *  
 * DEBUG! 
 * http://www.math.ethz.ch/education/bachelor/lectures/hs2012/math/nummath_cse/sol5.pdf
 */

void
choleskyTriDiagArrowFact(double *dia, double *off, double *bot, Uint n) {

  Uint i;
  double sum =0;

  assert(n > 3);
  dia[0] = sqrt(dia[0]);  //L
  off[0] = off[0]/dia[0]; //B = M
  bot[0] = bot[0]/dia[0]; //E

  for(i=1; i < n-3; i++){

    dia[i] = dia[i] - off[i-1]*off[i-1];
    assert(dia[i] >= 0); // pos def matrix
    dia[i] = sqrt(dia[i]);
    off[i] = off[i]/dia[i]; //off = M

    bot[i] = bot[i] - bot[i-1]*off[i-1] / dia[i];
    sum += bot[i-1]*bot[i-1]; //bot = E
  }

  dia[n-3] = dia[n-3] - off[n-4]*off[n-4];
  assert(dia[n-3] >= 0);
  dia[n-3] = sqrt(dia[n-3]);
  off[n-3] = off[n-3] - bot[n-4]*off[n-4]/dia[n-3];
  bot[n-3] = bot[n-3] - bot[n-4]*off[n-4]/dia[n-3];

  dia[n-2] = dia[n-2] - off[n-3]*off[n-3] - sum;
  assert(dia[i] >= 0);
  dia[n-2] = sqrt(dia[n-2]);
  
  return;
}

 
/*----------------------------- splines_periodic -----------------------------
 *    
 * @brief http://svn.r-project.org/R/trunk/src/library/stats/src/splines.c
 * @author Steve Hoffmann 
 *   
 */
 

double*
splines_periodic(double *x, double *y, Uint n) {
  int i;
  double sum = .0;
  double *dia, *off, *bot, *c;

  dia = ALLOCMEMORY(NULL, NULL, double, n);
  off = ALLOCMEMORY(NULL, NULL, double, n);
  bot = ALLOCMEMORY(NULL, NULL, double, n);
  c = ALLOCMEMORY(NULL, NULL, double, n);

  off[0]  = x[1] - x[0];
  off[n-2]= x[n-1] - x[n-2];
  dia[0] = 2.0 * (off[0] + off[n-2]);
  bot[0] = x[n-1] - x[n-2];
  c[0] = (y[1] - y[0])/off[0] - (y[n-1] - y[n-2])/off[n-2];

  for(i = 1; i < n-1; i++) {
    off[i] = x[i+1] - x[i];
    dia[i] = 2.0 * (off[i] + off[i-1]);
    bot[i] = 0;
    c[i] = (y[i+1] - y[i])/off[i] - (y[i] - y[i-1])/off[i-1];
  }

  choleskyTriDiagArrowFact(dia, off, bot, n);

  //do forward substititution
  c[0] = c[0]/dia[0];
  for(i=1; i < n-2; i++) {
    c[i] = (c[i] - off[i-1] * c[i-1])/dia[i];
    //for periodic boundary condition
    sum += bot[i-1] * c[i-1];
  }
  c[n-2] = (c[n-2] - off[n-3] * c[n-3] - sum) / dia[n-2];

  //now backward substitution
  c[n-2] = c[n-2]/dia[n-2];
  //for periodic boundary condition
  c[n-3] = (c[n-3] - off[n-3] * c[n-2])/dia[n-3];

  for(i=n-4; i >= 0; i--) { 
    c[i] = (c[i] - off[i]*c[i] - bot[i])/dia[i];
  }

  c[n-1] = c[0];

  for(i=0; i < n-1; i++) {
    sum = x[i+1] - x[i];
    dia[i] = (y[i+1]-y[i])/sum - sum*(c[i+1]+2.0*c[i]);
    off[i] = (c[i+1]-c[i])/sum;
    c[i] = 3.0*c[i];
  }
  
  dia[n-1] = dia[0];
  off[n-1] = off[0];
  c[n-1] = c[0];

  return c;
}


/*--------- function pretty was adapted from the R project pretty.c ----------
 *    
 * @brief 
 * @author Steve Hoffmann 
 *   
 */
 

double *
prettyarray(double min, double max, Uint n, Uint *r)
{

  double *arr = NULL;
  double eps = 1e-7;
  double cell, U, base, unit, ns, nu;
  double shrink = 0.75;
  double h = 1.5; 
  double h5 = .5 + 1.5*h;
  double d = max-min;
  double delta;
  int min_n =  n / 3;
  char i_small = 0;
  int k, newn;

  if(max == 0 && min == 0) {
    cell = 1;
    i_small = 1;
  } else {
    cell = MAX(fabs(min),fabs(max));
    U = (1 + (h5 >= 1.5*h+.5)) ? 1/(1+h) : 1.5/(1+h5);
    i_small = d < cell * U * MAX(1,n) * DBL_EPSILON *3; 
  }

//  fprintf(stderr, "small:%d, cell:%f, U:%f, d:%f\n", i_small, cell, U, d); 

  if(i_small) {
    if(cell > 10) cell = 9 + cell/10;
    cell *= shrink;
    if(min_n > 1) cell /= (double)min_n;
  } else {
    cell = d;
    if(n > 1) cell /= (double)n;
  }

//  fprintf(stderr, "cell: %f\n", cell);

  if(cell < 20*DBL_MIN) {
     cell = 20*DBL_MIN;
  } else if(cell * 10 > DBL_MAX) { 
     cell = .1*DBL_MAX;
  }

//  fprintf(stderr, "cell: %f\n", cell);
  
  base = pow(10., floor(log10(cell)));
  unit = base;
    if((U = 2*base)-cell <  h*(cell-unit)) { 
        unit = U;
        if((U = 5*base)-cell < h5*(cell-unit)) { 
            unit = U;
            if((U =10*base)-cell <  h*(cell-unit)) 
                unit = U; 
        }
    }
   
//  fprintf(stderr, "unit: %f\n", unit);

  ns = floor(min/unit+eps);
  nu = ceil (max/unit-eps);
 
// fprintf(stderr, "initial ns-nu [%f,%f] - %f %f %f\n", ns, nu, ns*unit, min, eps*unit);

  while(ns*unit > min + eps*unit) ns--;
  while(nu*unit < max - eps*unit) nu++;

  k = (int)(0.5 + nu - ns);
//  fprintf(stderr, "itered ns-nu [%f,%f] - %d, %d\n", ns, nu, k, min_n); 

  if(k < min_n) {
    k = min_n - k;
    if(ns >= 0.) {
        nu += k/2;
        ns -= k/2 + k%2;/* ==> nu-ns = old(nu-ns) + min_n -k = min_n */
    } else {
        ns -= k/2;
        nu += k/2 + k%2;
    }
    newn = min_n;
  } else {
    newn = k;
  }

  arr = ALLOCMEMORY(NULL, NULL, double, newn+2);
  delta = (nu-ns)*unit/((double)newn);
  //fprintf(stderr, "from %f to %f, delta:%f, newn:%d\n", nu, ns, delta, newn);

  for(int i=0; i < newn; i++) {
    arr[i] = ns*unit + (((double)i)*delta);
  }
  
  *r = newn;
  return arr;
}


Uint *
bin(double *x, Uint n, double **breaks, Uint *nbins) {

  double *bins, min, max;
  Uint i, newn, curbin = 0, *cnt;

  assert(n > 0);

  qsort(x, n, sizeof(double), cmp_dbl_qsort);
  min = x[0];
  max = x[n-1]; 
 
  if(*nbins == 0) *nbins = ceil(log2(n) + 1); 

  fprintf(stderr, "\n[%f,%f]\n", min, max);
  bins = prettyarray(min, max, *nbins, &newn);

  for(i=0; i < newn; i++) fprintf(stderr, "%d %f\n", i, bins[i]);
  cnt = ALLOCMEMORY(space, NULL, Uint, newn);
  memset(cnt, 0, sizeof(Uint)*newn);
  cnt[0] = 0;

  for(i=0; i < n; i++) {
    while(curbin+1 < newn && x[i] >= bins[curbin+1]) {
        cnt[curbin+1] = 0;
        curbin++;
    }
    if(curbin > 10) fprintf(stderr, "%f -> bin[%f]\n", x[i], bins[curbin]);
    cnt[curbin]++;
  }

  *nbins = newn;
  *breaks = bins;
  
  return cnt;
}

/*-------------------------------- dist_uint ---------------------------------
 *    
 * @brief absolute distances
 * @author Steve Hoffmann 
 *   
 */
 
Uint
dist_uint (Uint a, Uint b)
{
  if(a > b) {
    return a-b;

  }
return b-a;
}
