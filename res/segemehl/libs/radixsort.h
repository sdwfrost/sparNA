#ifndef RADIXSORT_H
#define RADIXSORT_H

/*
 *  radixsort.h
 *  segemehl
 *
 *  Created by Steve Hoffmann on 08.02.10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include <basic-types.h>

void
bl_radixSort(void *space, void *toSort, 
			 size_t size, size_t nelem, 
			 Uint (*keyaccess)(void *), 
			 Uint bits);

void
bl_radixSortKeyFirst(void *space, void *toSort, 
					 size_t size, size_t nelem, 
					 Uint bits);

void
bl_radixSortUint(void *space, Uint *toSort, 
				 size_t nelem, 
				 Uint bits);

#endif
