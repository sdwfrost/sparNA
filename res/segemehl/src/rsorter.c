/*
 *  rsorter.c
 *  segemehl
 *
 *  Created by Steve Hoffmann on 08.02.10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <assert.h>
#include "radixsort.h"
#include "basic-types.h"


Uint
uintkey(void *key) {
	Uint* number = (Uint*) key;
	return *number;
}

int main(int argc, char** argv) {
	int i, len = 200000000, range = 20000000;
	unsigned int iseed = (unsigned int)time(NULL);	
	double difmatch;
	
	time_t startmatch, endmatch;
	unsigned int a[]={
		123,432,654,3123,654,2123,543,131,653,123,
		533,1141,532,213,2241,824,1124,42,134,411,
		491,341,1234,527,388,245,1992,654,243,987};
	
	Uint *lr, *lr2, *lr3;
	
	printf("Before radix sort:\n");
	for(i=0; i<sizeof a/sizeof(unsigned int); ++i) 
		printf(" %d", a[i]);
	putchar('\n');
	
	bl_radixSort(NULL, a, sizeof(Uint), sizeof a/sizeof(Uint), uintkey, 1);
	
	printf("After radix sort:\n");
	for(i=0; i<sizeof a/sizeof(unsigned int); ++i) 
		printf(" %d", a[i]);
	putchar('\n');
	
	srand (iseed);
	lr  = malloc(sizeof(Uint)*len);
	lr2 = malloc(sizeof(Uint)*len);
	lr3 = malloc(sizeof(Uint)*len);
	
	for(i=0; i < len; i++) {
		lr[i] = rand()%range;
		lr2[i] = lr[i];
		lr3[i] = lr[i];
	}
	
    time (&startmatch);
	bl_radixSort(NULL, lr, sizeof(Uint), len, uintkey, 16);
	time (&endmatch);	
	difmatch = difftime (endmatch, startmatch);
	printf("sorting took %f seconds\n", difmatch);
	
	time (&startmatch);
	bl_radixSortKeyFirst(NULL, lr2, sizeof(Uint), len, 16);
	time (&endmatch);	
	difmatch = difftime (endmatch, startmatch);
	printf("sorting took %f seconds\n", difmatch);
	
	time (&startmatch);
	bl_radixSortUint(NULL, lr3, len, 16);
	time (&endmatch);	
	difmatch = difftime (endmatch, startmatch);
	printf("sorting took %f seconds\n", difmatch);
	
	for(i=1; i < len; i++) {
		if(lr[i-1]>lr[i]) {
			printf("lr[%d]=%d > lr[%d]=%d", i-1, lr[i-1], i, lr[i]);
			return 0;
		}
		assert(lr[i] == lr2[i]);
		assert(lr[i] == lr3[i]);
	}
	
	return 0;
}

