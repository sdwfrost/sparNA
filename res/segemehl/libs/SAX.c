
/*
 *  SAX.c
 *  Keoghs symbolic aggregate approximation
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 07/06/2011 09:19:20 PM CEST
 *  
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "basic-types.h"
#include "info.h"
#include "manopt.h"
#include "stringutils.h"
#include "mathematics.h"
#include "fileio.h"
#include "SAX.h"

/*---------------------------------- bl_SAX ----------------------------------
 *    
 * @brief piecewise aggregate approximation of a 2d timeseries of length n
 * into w windows and creation of its symbolic representation.
 * @author Steve Hoffmann 
 *   
 */

char *
bl_SAX (void *space, double *C, Uint n, Uint w, Uint b, Uint *P, Uint **SxP)
{

  Uint i, j;
  double *Cb;
  char *Sx;

  assert(b > 0 && b < 5);

  Cb = ALLOCMEMORY(space, NULL, double, w);
  Sx = ALLOCMEMORY(space, NULL, char, w+1);

  for(i=0; i < w; i++) {
    Cb[i] = 0;
    (*SxP)[i] = P[(n/w)*(i)+1];
    for(j=(n/w)*(i)+1; j < (n/w)*(i+1); j++) {
      Cb[i] += C[j];
    }
  //  fprintf(stderr,"%d -> %f/%f=%f\n", j, Cb[i], (double)w/n, Cb[i]/((double)w/n));
    Cb[i]/=(double)w/n;
  }


  for(i=0; i< w; i++) {
    Sx[i] = equichar[b-1][0];
    for(j=0; j < b; j++) {
      if(Cb[i] > equiprob[b-1][j]) {
        Sx[i] = equichar[b-1][j+1];
      }
    }
  //  fprintf(stderr,"Cb[%d]=%fi -> %c\n", i, Cb[i], Sx[i]);
  }

  FREEMEMORY(space, Cb);

  Sx[w] = 0;
  return Sx;
}


Uint
descendRSS(void *space, double *Q, Uint u, Uint k, Uint v) {
  Uint i, min_1= k, min_2 = k;
  double mu_1 = .0, mu_2 = .0;

  if(u >= k || k >= v) { 
    fprintf(stderr, "%d < %d < %d ?\n", u, k, v);
    return k;
    exit(-1);
  }
  assert( u < k && k < v );

  for(i=u; i < k; i++) {
    mu_1 += Q[i];
    if(Q[i] < Q[min_1] || Q[min_1] < 5.0) {
      min_1 = i;
    }
  }
  mu_1 /= k-u+1;

  for(i=k+1; i < v; i++) {
    mu_2 += Q[i];
    if(Q[i] < Q[min_2] && Q[min_2] > 5.0) {
      min_2 = i;
    }
  }
  mu_2 /= v-k+1;

  if(mu_1 < mu_2) {
    return min_1;
  } else {
    return min_2;
  }

  return k;
}


void
slidingdiffRSS(void *space, double *Q, Uint n, Uint w, Uint off) { 
 double f = 0.12;
  Uint i, j, u, minarg2, noofbreaks = (w/(f*w)), k=0, 
       next = 0, last = 0, pos = 0, *bl = NULL, *bl2=NULL, q=0;
  double epsilon = 0.05, sum = 0, t = w*0.10, lsum, rsum, *D;
  breakpoints_t *bp2;

  fprintf(stderr, "looking for %d breakpoints\n", noofbreaks);
  D = ALLOCMEMORY(space, NULL, double, w+1);

  for(i=0; i < n-w; i+=off) {
    for(j=0, sum=0; j< w; j++) {
      sum += Q[i+j];
      D[j] = Q[i+j+1] - Q[i+j];
    }

    if(sum > t) { 
         
      bp2 = bl_RSSmatrix(space, D, w, f*w, noofbreaks);
      minarg2 = 0;
      for(j=1; j < noofbreaks+1; j++) {
        if(bp2[j].RSS > .00000000000000000001 && 
            bp2[j].BIC < bp2[minarg2].BIC-(epsilon*bp2[minarg2].BIC)) {
          minarg2 = j;
        }
      }
 
 
      if(bp2[minarg2].RSS > 0.0) {


        last = 0;

        for(j=0; j < bp2[minarg2].noofbreaks; j++) {
          if(j < bp2[minarg2].noofbreaks-1) {
            next = bp2[minarg2].breaks[j+1];
          } else {
            next = w-1;
          }

          pos =  bp2[minarg2].breaks[j];
          
          for(u=last, lsum=0; u < pos; u++) {
            //lsum += D[u];
            lsum += Q[i+u];

          }
          for(u=pos, rsum=0; u < next; u++) {
            //rsum += D[u];
            rsum += Q[i+u];
          }
         fprintf(stderr, "entered minarg2 to register %d breaks (l:%f,r:%f)\n", 
             bp2[minarg2].noofbreaks, lsum, rsum);
          if(lsum/(pos-last) > 5.0 || rsum/(next-pos) > 5.0 ) { 
            if(q==0 || bl2[q-1] < i+pos-20) { 
              bl2 = ALLOCMEMORY(space, bl2, Uint, q+1);
              bl2[q] = i+pos;
              q++;
            }           
          }
          last = bp2[minarg2].breaks[j];
        }
      }


      for(j=1; j < noofbreaks+1; j++) {
          if(bp2[j].noofbreaks)
          FREEMEMORY(space, bp2[j].breaks);
      }
      FREEMEMORY(space, bp2);

    }
    if(k > 1000) break;
  }

  fprintf(stdout, "track name=events description=\"transcription events\" useScore=0\n");
  for(i=0; i < k; i++) { 
    fprintf(stdout, "%s\t%d\t%d\tevent\n", "chr15", bl[i], bl[i]);
  }


  fprintf(stdout, "track name=changes description=\"transcription changes\" useScore=0\n");
  for(i=0; i < q; i++) { 
    fprintf(stdout, "%s\t%d\t%d\tchange\n", "chr15", bl2[i], bl2[i]);
  }
}


void 
slidingRSS(void *space, double *Q, Uint n, Uint w, Uint off) { 
  double f = 0.12;
  Uint i, j, u, minarg, noofbreaks = (w/(f*w)), k=0, 
       next = 0, last = 0, pos = 0, *bl = NULL, *bl2=NULL, q=0;
  double epsilon = 0.05, sum = 0, t = w*0.10, lsum, rsum;
  Uint min_1, min_2, min;
  double mu_1 = .0, mu_2 = .0;

  breakpoints_t *bp;

  fprintf(stderr, "looking for %d breakpoints\n", noofbreaks);

  for(i=0; i < n-w; i+=off) {
    for(j=0, sum=0; j< w; j++) {
      sum += Q[i+j];
     }

    if(sum > t) { 
      bp = bl_RSSmatrix(space, &Q[i], w, f*w, noofbreaks);
   
      minarg = 0;
      for(j=1; j < noofbreaks+1; j++) {
        if(bp[j].RSS > .00000000000000000001 && 
            bp[j].BIC < bp[minarg].BIC-(epsilon*bp[minarg].BIC)) {
          minarg = j;
        }
      }


      if(bp[minarg].RSS > 0.0) { 
        for(j=0; j < bp[minarg].noofbreaks; j++) {
          bp[minarg].breaks[j] += i; 
        }

        /*
        //correction

        if(bp[minarg].noofbreaks) {
        u = k;
        while (u >= 1 && bp[minarg].breaks[0] <= bl[u-1]){
          u--;
        }

        if(u > 0)
          last = bl[u-1];
        else
          last = 0;
        }

        for(j=0; j < bp[minarg].noofbreaks; j++) {
          
          if(j < bp[minarg].noofbreaks-1) {
            next = bp[minarg].breaks[j+1];
          } else {
            next = i+w-1;
          }

          if(bp[minarg].breaks[j] > 50)
            last = MAX(bp[minarg].breaks[j]-50, last);
        
          //fprintf(stderr, "last:%d\n", last);
          pos = descendRSS(space, Q, last, bp[minarg].breaks[j], 
              MIN(next, bp[minarg].breaks[j]+50));
          
          bp[minarg].breaks[j] = pos;
          last = pos;
        } */

        //evaluate and register
        last = 0;
        for(j=0; j < bp[minarg].noofbreaks; j++) {
          if(j < bp[minarg].noofbreaks-1) {
            next = bp[minarg].breaks[j+1];
          } else {
            next = i+w-1;
          }

          pos =  bp[minarg].breaks[j];
         
          if(pos > 100) last = MAX(last, pos - 100);
          if(next > pos + 100) next = pos+100;

          min_1 = pos;
          for(u=last, lsum=0; u < pos; u++) {
            lsum += Q[u];    
            if(Q[u] < Q[min_1] || Q[u] < 5.0) {
              min_1 = u;
            }
          }

          mu_1 = (double)lsum/(pos-last+1);

          min_2 = pos;
          for(u=pos, rsum=0; u < next; u++) {
            rsum += Q[u];
            if(Q[u] < Q[min_2] && Q[min_2] > 5.0) {
              min_2 = u;
            }
          }

          mu_2 = (double)rsum/(next-pos+1);

          if(mu_1 < mu_2) {
            min = min_1;
          } else {
            min = min_2;
          }

          if(lsum/(pos-last) > 5.0 || rsum/(next-pos) > 5.0 ) { 
            if(k==0 || bl[k-1] < min-20) { 
              bl = ALLOCMEMORY(space, bl, Uint, k+1);
              bl[k] = min;
              fprintf(stderr, "%d\t", min);
              k++;
            } else { 
              fprintf(stderr, "%d*\t", min);
            }
          }
          last = bp[minarg].breaks[j];
        }

        fprintf(stderr, "\n");
      }

 

      for(j=1; j < noofbreaks+1; j++) {
        if(bp[j].noofbreaks)
          FREEMEMORY(space, bp[j].breaks);
      }
      FREEMEMORY(space, bp);

    }
    if(k > 1000) break;
  }

  fprintf(stdout, "track name=events description=\"transcription events\" useScore=0\n");
  for(i=0; i < k; i++) { 
    fprintf(stdout, "%s\t%d\t%d\tevent\n", "chr15", bl[i], bl[i]);
  }


  fprintf(stdout, "track name=changes description=\"transcription changes\" useScore=0\n");
  for(i=0; i < q; i++) { 
    fprintf(stdout, "%s\t%d\t%d\tchange\n", "chr15", bl2[i], bl2[i]);
  }


} 


int
main(int argc, char **argv) {
  stringset_t **csv;
  double *C, sum=0, mu, *Q;
  Uint *P, *SxP;
  char *Sx;
  Uint lines, i, k=0; 
  void *space = NULL;
  //breakpoints_t *bp; Uint j; Uint d; Uint noofbreaks=5;
  
  csv = readcsv(space, "startdist2.out", "\t", &lines);
//  C = ALLOCMEMORY(space, NULL, double, lines);
//  P = ALLOCMEMORY(space, NULL, Uint, lines);
  Q = ALLOCMEMORY(space, NULL, double, lines);

  /*
  for(i=0; i < lines-1; i++) {
    if((atoi(csv[i]->strings[0].str) > 5 || atoi(csv[i+1]->strings[0].str) > 5) &&
       ((d=atof(csv[i]->strings[0].str) - atof(csv[i+1]->strings[0].str)) > 0 || d < 0)){ 
      P[k] = i;
      C[k] = (d > 0) ? d-1 : d+1; 
      sum += C[k];
      k++;
    }
  }*/
  fprintf(stderr, "read csv.\n");

   for(i=0; i < lines; i++) {
    Q[i] = atof(csv[i]->strings[0].str);
  }
  
  fprintf(stderr, "converted csv.\n");
  fprintf(stderr, "starting segmentation.\n");
  
  slidingRSS(space, Q, lines-1, 1000, 500); 

/* 
  bp = bl_RSSmatrix (space, &Q[0], 100, 15, noofbreaks);
  for(i=0; i < noofbreaks+1; i++) {
    fprintf(stderr, "rss=%f\tll=%f\tBIC=%f\t", bp[i].RSS, bp[i].LL, bp[i].BIC);
    for(j=0; j < bp[i].noofbreaks; j++) {
      fprintf(stderr, "%d\t", bp[i].breaks[j]);
    }
    fprintf(stderr, "\n");
  }
 */
  exit(-1);

 
  for(i=0; i < lines; i++) {
    destructStringset(space, csv[i]);
  }
  
  FREEMEMORY(space, csv);

  mu = sum/k;

  //fprintf(stderr, "mean:%f\n", mu);
  for(i=0; i < k; i++) {
    /*center*/
    C[i]-=mu;

   // fprintf(stderr, "%f\n", C[i]);
    /*normalize*/
    C[i]/=sum;
  }
 
 // fprintf(stdout, "normalized:\n");
  for(i=0; i < k; i++) {
    fprintf(stdout, "%d\t%f\n", P[i], C[i]);
  }

  SxP = ALLOCMEMORY(space, NULL, Uint, (k/20)+1);
  Sx = bl_SAX(space, C, k, k/20, 4, P, &SxP);
  
  fprintf(stderr, "track\tname=SAX\tdescription=\"SAX of expression\"\n");

  for(i=0; i < k/20; i++) {
    fprintf(stderr, "chr15\t%d\t%d\t%c\n", SxP[i], SxP[i]+k/(k/20)-1, Sx[i]);
  }

//  fprintf(stderr, "%s\n", Sx);

  FREEMEMORY(space, Sx);
  FREEMEMORY(space, SxP);
  FREEMEMORY(space, C);
  FREEMEMORY(space, P);
  FREEMEMORY(space, Q);

  return 0;
}
