/**
 * list.h
 * implementation of a simple lineary linked list for object pointer
 *
 * @author Christian Otto
 * @email christian@bioinf.uni-leipzig.de
 * @company Bioinformatics, University of Leipzig
 * @date Wed Oct 15 11:39:42 CEST 2008
 */

/*
 * SVN
 * Revision of last commit: $Rev: 73 $
 * Author: $Author: steve $
 * Date: $Date: 2008-10-29 10:03:28 +0100 (Wed, 29 Oct 2008) $
 * Id: $Id$
 * Url: $URL$
 */

#ifndef LIST_H
#define LIST_H

#include <stdlib.h>
#include "basic-types.h"

#define LISTINC 1000
#ifndef BASEINC
#define BASEINC LISTINC
#endif

typedef struct {
  int next;
  int prev;
} Listelem;


typedef struct {
  Listelem *nodes;
  void *data;
  int first;
  int last;
  int nextfree;
  Uint numofelem;
  int allocelem;
  size_t sizeofelem;
} List;

int bl_listGetCur (List *l, Uint n);
void bl_listInit(List *l, int allocelem, size_t sizeofelem);
void bl_listDestruct(List *l, void (*rmv)(void*));
BOOL bl_listIsEmpty(List *l);
void bl_listResize(List *l);
void bl_listInsert(List *l, int cur, void *elem);
void* bl_listUnlink(List *l, Uint cur, void (*rmv)(void*));
void bl_listSweep(List *l);
Uint bl_listSize(List *l);
Uint bl_listBinarySearchInsert (List *l, void *elem, 
    Uint (*cmp)(Uint, void *, void*, void*), void *nfo);
void* bl_listGetElem (List *l, int cur);

#ifdef LISTTEST

typedef struct {
  Uint unsigned1;
  Uint unsigned2;
  int integer;
  char *string;

} listtestobj_t;

#endif

#endif /* LIST_H */
