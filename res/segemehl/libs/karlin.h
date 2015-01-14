#ifndef KARLIN_H
#define KARLIN_H
/*
 *
 *	karlin.h
 *  declaration for karlin.c
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 05/06/2008 08:55:40 PM CEST  
 *
 *  SVN
 *  Revision of last commit: $Rev: 19 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-05-14 15:43:29 +0200 (Wed, 14 May 2008) $
 *
 *  Id: $Id: karlin.h 19 2008-05-14 13:43:29Z steve $
 *  Url: $URL: file:///homes/bierdepot/steve/svn/segemehl/trunk/libs/karlin.h $
 */


typedef struct karlin_s {
  double lambda;
  double H;
  double K;
} karlin_t;


int karlinunitcostpp(void *space, double *lambda, double*H, double *K);
double significance (double lambda,double K,double multiplier, int score);
double evalue (double lambda,double K,double multiplier, int score);
double explength(Uint m, Uint n, double H, double K);
double effSubjectLength(Uint m, Uint n, double H, double K);
double effQueryLength(Uint m, Uint n, double H, double K);
double spacemult(Uint m, Uint n, double H, double K);
#endif
