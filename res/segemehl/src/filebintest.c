
/*
 *  filebintest.c
 *  segemehl
 *
 *  Created by Steve Hoffmann on 09.02.10.
 *  Copyright 2010 University Leipzig. 
 *  All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <sys/types.h>
#include "basic-types.h"
#include "radixsort.h"
#include "fileBins.h"
#include "filebintest.h"

unsigned char mute=0;

LLint getLineColKey(char *src, void *nfo) {

  lineColKeyInfo_t *lcknfo = (lineColKeyInfo_t*) nfo;
  char c;
  char *s, *start;
  int sepcnt;
  LLint key = LLONG_MIN;

  start = src;
  s = src;
  c = lcknfo->sep;
  sepcnt = 0;

  if (!*s) return 0;

  do {
    if (*s == c) {
      sepcnt++; 

      if(sepcnt == lcknfo->col+1) break;
      if(sepcnt == lcknfo->col) {
        start = s+1;
      }
    }		
  } while (*++s);

  if(sepcnt-1 == lcknfo->col) {
    *s = '\0';
    key = atol(start);
    *s = c;
  }

  if(sepcnt == lcknfo->col && *s == 0 && start < s) {
    key = atol(start);
  }
  return key;

}


int main (int argc, char** argv) {
	bl_fileBins_t *myBins;
    bl_fileBin_t* fb;
	lineColKeyInfo_t nfo;
	int i, j, a, b, c, d; 
	char *line;
	time_t startsuf, endsuf;
	double difsuf;
    char *buffer[]={"string1","string2","string3","string4", "string5"};
    
    int range = 2000000, 
		 noofbins = 5, 
		 nooflines = 2000;
	
	unsigned int iseed = (unsigned int)time(NULL);
	line = malloc(1000);
	nfo.sep='\t';
	nfo.col=1;
	
    myBins = calloc(1, sizeof(bl_fileBins_t));
	bl_fileBinsAdd(NULL, myBins, noofbins, bl_fileBinCClassAssign, 
        buffer, NULL, "bla", 3);

	for(i=0; i < noofbins; i++) {
		fprintf(stderr, "file %d: name=%s\n", i, myBins->b[i].fname);
	}
	
	srand (iseed);
	
	for(i=0; i < noofbins; i++) {
      fb = bl_fileBinsFind(NULL, myBins, bl_fileBinsCClassSelect, buffer[i]);
      fprintf(stderr, "try to open %s (i:%d)\n", fb->fname, i);
      bl_fileBinsOpen(NULL, fb, "w");

		for(j=0; j < nooflines; j++) {
			a = rand()%range;
			b = rand()%range;
			c = rand()%range;
			d = rand()%range;
			
			sprintf(line, "%d \t %d \t %d \t %d%c", 
					a, b, c, d, '\0');
            bl_fileBinsWriteLn(NULL, fb, line);
//			fprintf(fp, "%s", line);
//            myBins->b[i].lines++;
			memset(line, 0, 1000);
		}
        //bl_fileBinsClose(NULL, fb);
	
/*		fp = fopen(fb->fname, "r");
		if (fp == NULL){
			fprintf(stderr, "main Opening of file %s failed. Exit forced.\n", 
                fb->fname);
			exit(EXIT_FAILURE);
		}
        fclose(fp);
        */
    }
    free(line);

    fprintf(stderr, "start to sort\n");
    time (&startsuf);
	bl_fileBinsSortLine(NULL, myBins, 1, "trash.txt", 1, getLineColKey, &nfo);	
    time (&endsuf);
    difsuf = difftime (endsuf, startsuf);
    fprintf(stderr, "sort has taken %f seconds.\n", difsuf);
    bl_fileBinsDestruct(NULL, myBins);
    free(myBins);

    return 0;
}

