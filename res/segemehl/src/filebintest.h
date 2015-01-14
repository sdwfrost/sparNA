#ifndef FILEBINTEST_H
#define FILEBINTEST_H
/*
 *  filebintest.h
 *  segemehl
 *
 *  Created by Steve Hoffmann on 09.02.10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


typedef struct lineColKeyInfo_s {
	Uint col;
	char sep;
} lineColKeyInfo_t;

LLint getLineColKey(char *s, void *nfo);


#endif