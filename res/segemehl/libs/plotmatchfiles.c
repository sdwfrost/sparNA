
/*
 *  plotmatchfiles.c
 *  
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 10/22/2010 05:17:21 PM CEST
 *  
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <sys/types.h>
#include <unistd.h>
#include "alignment.h"
#include "stringutils.h"
#include "basic-types.h"
#include "mathematics.h"
#include "matfile.h"
#include "bitVector.h"
#include "info.h"
#include "fileio.h"
#include "matchfiles.h"
#include "evalmatchfiles.h"


void
bl_matchfileCROSSERRGNUPLOT(void *space, matchfileindex_t *index) {
  char *name = "tempname.txt";
  double *y;
  Uint i;
 
  y = ALLOCMEMORY(space, NULL, double, (Uint)index->maxreadlen);

  for(i=0; i < index->maxreadlen; i++) {
    y[i] = index->P_ERR[i]/index->noofreads;
  }
 
  FILE *pipe = popen("gnuplot -persist 2>/dev/null","w");
  writeY(name, y, index->maxreadlen, 0, 0);
   
  fprintf(pipe, "set title 'read position dependent error rates'\n");
  fprintf(pipe, "set xlabel 'read position'\n");
  fprintf(pipe, "set ylabel 'error rate'\n"); 
  fprintf(pipe, "plot '%s' using 1:2 notitle w lines\n", name);

  pclose(pipe);
  FREEMEMORY(space, y);
}


void
bl_matchfileQERRGNUPLOT(void *space, matchfileindex_t *index) {
  char *name = "tempname.txt";
  char *name2 = "tempname2.txt";
  Uint i;
  double *arr, *arr2;

  arr = ALLOCMEMORY(space, NULL, double, QRNGE);
  arr2 = ALLOCMEMORY(space, NULL, double, QRNGE);

  for(i=0; i < QRNGE; i++) {
    arr[i] = (double)index->Q_ERR[i]/index->Q_N[i];  
    if(i==0 && !index->Q_N[i]) arr[i] = 1;
    arr2[i] = pow(10,(((double)i)/-10.0));
  }
 
  FILE *pipe = popen("gnuplot -persist 2>/dev/null","w");
  writeY(name, arr, QRNGE, 0, 0); 
  writeY(name2, arr2, QRNGE, 0, 0);

  fprintf(pipe, "set title 'quality dependent error rates'\n");
  fprintf(pipe, "set xlabel 'quality'\n");
  fprintf(pipe, "set ylabel 'errorrate'\n"); 
  fprintf(pipe, "plot '%s' using 1:2 w points notitle", name);
  fprintf(pipe, "   , '%s' using 1:2 smooth bezier title 'observed error' w lines", name);
  fprintf(pipe, "   , '%s' using 1:2 title 'expected error' w lines\n", name2);

  pclose(pipe);
  FREEMEMORY(space, arr);
  FREEMEMORY(space, arr2);
}

void
bl_matchfilePERRGNUPLOT(void *space, matchfileindex_t *index) {
  char *name = "tempname.txt";
  double *y;
  Uint i;
 
  y = ALLOCMEMORY(space, NULL, double, (Uint)index->maxreadlen);

  for(i=0; i < index->maxreadlen; i++) {
    y[i] = (double)index->P_ERR[i]/index->noofreads;
  }
 
  FILE *pipe = popen("gnuplot -persist 2>/dev/null","w");
  writeY(name, y, index->maxreadlen, 0, 0);
   
  fprintf(pipe, "set title 'read position dependent error rates'\n");
  fprintf(pipe, "set xlabel 'read position'\n");
  fprintf(pipe, "set ylabel 'error rate'\n"); 
  fprintf(pipe, "plot '%s' using 1:2 notitle w lines\n", name);

  pclose(pipe);
  FREEMEMORY(space, y);
}

void
bl_matchfileSUBGNUPLOT(void *space, matchfileindex_t *index) {
  char *name = "tempname.txt";
  double *y,*x,*z, sum=0;
  double cnt[4] = {0,0,0,0};

  Uint i, j, k, l, u=0;
 
  x = ALLOCMEMORY(space, NULL, double, 37);
  y = ALLOCMEMORY(space, NULL, double, 37);
  z = ALLOCMEMORY(space, NULL, double, 37);
 
  for(i=0; i < 5; i++) {
    for(j=0; j < 5; j++) { 
      for(k=0; k < QRNGE; k++) {
        for(l=0; l < index->maxreadlen; l++) {
          cnt[i] += MATRIX4D(index->submatrix, 6, 
              QRNGE, MAXREADLENGTH, i, j, k, l);
        }
      }
    }
  }

  for(i=0; i < 5; i++) {
    for(j=0; j < 5; j++) {
      sum = 0;
      for(k=0; k < QRNGE; k++) {
        for(l=0; l < index->maxreadlen; l++) {
          sum += MATRIX4D(index->submatrix, 6, 
              QRNGE, MAXREADLENGTH, i, j, k, l);
        }
      }
      fprintf(stderr, "i:%d, j:%d sum: %f\n", i, j, sum);
      x[u] = i;
      y[u] = j;
      z[u] = log10(sum/cnt[i]);
      u++;
    }
  }

  FILE *pipe = popen("gnuplot -persist 2>/dev/null","w");
  writeXYZ(name, x, y, z, 25);
   
  fprintf(pipe, "set title 'substitution rates (log10)'\n");
  fprintf(pipe, "set xlabel 'read'\n");
  fprintf(pipe, "set ylabel 'reference'\n"); 
  fprintf(pipe, "set xtics ('A' 0,'C' 1,'G' 2,'T' 3, '-' 4)\n");
  fprintf(pipe, "set ytics ('A' 0,'C' 1,'G' 2,'T' 3, '-' 4)\n"); 
  fprintf(pipe, "set label 1 'bla' at 1,1\n");
  fprintf(pipe, "plot '%s' using 1:2:3 w image\n", name);

  pclose(pipe);
  FREEMEMORY(space, y);
  FREEMEMORY(space, x);
  FREEMEMORY(space, z);

}


void
bl_matchfileCOVGNUPLOT(void *space, matchfileFrame_t *frame) {
  char *name = "tempname.txt";
  Uint i, *arr;

  arr = ALLOCMEMORY(space, NULL, Uint, frame->width);

  for(i=0; i < frame->width; i++) {
    arr[i] = frame->cs[i].len;
  }
 
  FILE *pipe = popen("gnuplot -persist 2>/dev/null","w");
  writeYUint(name, arr, frame->width, frame->start, 0);
 
  fprintf(pipe, "set title 'frame coverage %s[%d,%d]'\n", frame->chrname, 
      frame->start, frame->start+frame->width);
  fprintf(pipe, "set xlabel 'frame position'\n");
  fprintf(pipe, "set ylabel 'coverage'\n"); 
  fprintf(pipe, "plot '%s' using 1:2 notitle w lines\n", name);

  pclose(pipe);
}


void
bl_matchfileRSSGNUPLOT(void *space, matchfileFrame_t *frame,
    matchfileFrameStats_t *stats) {

  char *name  = "tempname.txt";
  char *name2 = "tempname2.txt";
  Uint *data, xmax, ymax, total, width;
  double avg, maxval;  
  FILE *pipe = popen("gnuplot -persist 2>/dev/null","w");

  data = stats->dist_rss;
  xmax = stats->dist_rss_xmax;
  ymax = stats->dist_rss_ymax;
  total = stats->rss;
  width = frame->width;
  avg = (double)total/width;

  maxval = MAX(poisson(avg, avg)*width, ymax);
  writeYUint(name, data, xmax, 0, 0);
  writeYUintNorm(name2, data, xmax, 0);
 
  fprintf(pipe, "set title 'read start site distribution %s[%d,%d]'\n", 
      frame->chrname, frame->start, frame->start+width);
  fprintf(pipe, "set xlabel 'number of read startsites'\n");
  fprintf(pipe, "set ylabel 'number of genomic loci'\n"); 
  fprintf(pipe, "poissondraw(x)= (x>=0) ? exp(-l) * l**(x) / gamma(x+1) : 1\n"); 
  fprintf(pipe, "normal(x) = (1/(sd*sqrt(2*pi)))*exp(-(x-mu)**2/(2*sd**2))\n");
  fprintf(pipe, "l=%f; mu=%f; sd=1\n", avg, avg); 

  if(xmax < 100) {
    fprintf(pipe, "fit normal(x) 'tempname2.txt' via mu,sd\n");
    fprintf(pipe, "fit poissondraw(x) 'tempname2.txt' via l\n");
  }
  
  fprintf(pipe, "set label 1 'lambda=%%g',l at %d,%d\n", (int)(xmax/2), (int)(maxval/2));
  fprintf(pipe, "set label 2 'mu=%%g',mu at %d,%d\n", (int)(xmax/2), (int)(maxval/3));
  fprintf(pipe, "set parametric\n");
  fprintf(pipe, "set trange [0:%d]\n", xmax+1);
  fprintf(pipe, "set xrange [0:%d]\n", xmax+1);
  fprintf(pipe, "set yrange [0:%d]\n", (int)(maxval+(0.1*maxval)));

  fprintf(pipe, "plot '%s' using 1:2 notitle w points,", name);
  fprintf(pipe, "     '%s' using 1:2 smooth csplines title 'read start sites spline' w lines,  ", name);
  fprintf(pipe, "     '%s' using 1:2 smooth bezier title 'read start sites bezier' w lines", name);
 
  if(xmax < 100) {
    fprintf(pipe, "    , t,poissondraw(t)*%d w lines title 'poisson'", width);
    fprintf(pipe, "    , t,normal(t)*%d  w lines    title 'gaussian'\n", width);
  } else {
    fprintf(pipe, "\n");
  }

  pclose(pipe);
}
