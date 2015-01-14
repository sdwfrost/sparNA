/**
 * list.c
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

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <assert.h>
#include "debug.h"
#include "basic-types.h"
#include "sort.h"
#include "list.h"

/*------------------------------ bl_listInit -----------------------------------
 *    
 * @brief 	init list
 * @author 	Christian Otto
 *   
 */
void bl_listInit(List *l, int allocelem, size_t sizeofelem){
  if (allocelem <= 0){
    DBG("list.c: Attempt to initialize a list of size %d.\
        Exit forced.\n", allocelem);
    exit(-1);
  }
  if (sizeofelem <= 0){
    DBG("list.c: Attempt to initialize a list with sizeofelem %d.\
        Exit forced.\n", sizeofelem);
    exit(-1);
  }
  l->nodes = (Listelem *) malloc(allocelem * sizeof(Listelem));  
  if (l->nodes == NULL){
    DBG("list.c: Memory allocation for nodes failed. Exit forced.\n", NULL);
    exit(-1);
  }
  l->data = malloc(allocelem * sizeofelem);  
  if (l->nodes == NULL){
    DBG("list.c: Memory allocation for data failed. Exit forced.\n", NULL);
    exit(-1);
  }
  l->first = -1;
  l->last = -1;
  l->nextfree = 0;
  l->numofelem = 0;
  l->allocelem = allocelem;
  l->sizeofelem = sizeofelem;
}

/*----------------------------- bl_listDestruct --------------------------------
 *    
 * @brief 	destruct list,
 *              remove method for elems as parameter possible
 * @author 	Christian Otto
 *   
 */
void bl_listDestruct(List *l, void (*rmv)(void*)){
  Uint cur;
  char *p;
  if (rmv != NULL && l->numofelem > 0){
    p = l->data;
    for (cur = l->first; cur != -1; cur = l->nodes[cur].next){
      rmv(p + (cur * l->sizeofelem));
    }
  }
  free(l->nodes);
  free(l->data);
  l->first = 0;
  l->last = 0;
  l->nextfree = 0;
  l->numofelem = 0;
  l->allocelem = 0;
  l->sizeofelem = 0;
}

/*----------------------------- bl_listIsEmpty ---------------------------------
 *    
 * @brief 	returns if the container is empty
 * @author 	Christian Otto
 *   
 */
BOOL bl_listIsEmpty(List *l){
  return (l->numofelem == 0);
}

/*----------------------------- bl_listInsert ----------------------------------
 *    
 * @brief 	adds element after an given element in the list
 *              (at beginning for cur == -1, at end for cur == l->last)
 * @author 	Christian Otto
 *   
 */
void bl_listInsert(List *l, int cur, void *elem){
  char *p;
  if (cur > l->allocelem || (cur < 0 && cur != -1)){
    return;
  }
  /* reallocation */
  if (l->nextfree >= l->allocelem){
    l->nodes = (Listelem *) realloc(l->nodes, sizeof(Listelem) *
        (l->allocelem + BASEINC));
    if (l->nodes == NULL){
      DBG("list.c: Memory reallocation of nodes failed. Exit forced.\n", NULL);
      exit(-1);
    }
    l->data = realloc(l->data, l->sizeofelem * (l->allocelem + BASEINC));
    if (l->data == NULL){
      DBG("list.c: Memory reallocation of data failed. Exit forced.\n", NULL);
      exit(-1);
    }
    l->allocelem += BASEINC;
  }  
  p = (char *) l->data;
  /* insert data */
  memmove(p + (l->nextfree * l->sizeofelem), elem, l->sizeofelem);
  /* insert at begin (or in empty list) */
  if (cur == -1){
    l->nodes[l->nextfree].next = l->first;
    l->nodes[l->nextfree].prev = -1;
    if (l->first != -1){
      l->nodes[l->first].prev = l->nextfree;
    } else {      
      l->last = l->nextfree;
    }
    l->first = l->nextfree;
  }
  /* insert after elem cur */
  else {
    l->nodes[l->nextfree].prev = cur;
    l->nodes[l->nextfree].next = l->nodes[cur].next;
    /* new elem is last one */
    if (cur == l->last){
      l->last = l->nextfree;
    }
    /* otherwise */
    else {      
      l->nodes[l->nodes[l->nextfree].next].prev = l->nextfree;
    }
    l->nodes[cur].next = l->nextfree;
  }
  l->numofelem++;
  l->nextfree++;
}

/*-------------------------------- bl_listGetCur ------------------------------
 *    
 * @brief get nth element from list
 * @author Steve Hoffmann 
 *   
 */
 
int
bl_listGetCur (List *l, Uint n)
{
  Uint i=0; 
  int cur = l->first;

  while(i < n && cur != -1) {
    cur = l->nodes[cur].next;
    i++;
  }
	
  return cur;
}


/*------------------------------ bl_listGetElem ------------------------------
 *    
 * @brief getting the element at cursor position cur
 * @author Steve Hoffmann 
 *   
 */
 
void*
bl_listGetElem (List *l, int cur)
{
  char *p;

  if(cur > l->allocelem || cur < 0)
    return NULL;

  p = l->data;
  p += (cur * l->sizeofelem);
  
  return p;
}

/*------------------------------ bl_listUnlink ---------------------------------
 *    
 * @brief 	removes element from the list
 *              does not free
 * @author 	Christian Otto
 *   
 */

void* bl_listUnlink(List *l, Uint cur, void (*rmv)(void*)){

  char *p, *elem;
  p = (char *) l->data;
  if (cur > l->allocelem || cur < 0){
    return NULL;
  }
  
  if(cur != l->first && l->nodes[cur].next == -1 &&
      l->nodes[cur].next == l->nodes[cur].prev) {
   /*previously unlinked element*/
    return NULL;
  }

  elem = (char *) malloc(l->sizeofelem);
  memmove(elem, p + (cur * l->sizeofelem), l->sizeofelem);
  
  if (rmv != NULL){
    rmv(p + (cur * l->sizeofelem));
  }
  
  if (l->nodes[cur].prev != -1){
    l->nodes[l->nodes[cur].prev].next = l->nodes[cur].next;
  } else {
    l->first = l->nodes[cur].next;
  }

  if (l->nodes[cur].next != -1){
    l->nodes[l->nodes[cur].next].prev = l->nodes[cur].prev;
  } else {
    l->last = l->nodes[cur].prev;
  }

  l->nodes[cur].prev = -1;
  l->nodes[cur].next = -1;
  l->numofelem--;
  
  return elem;
}

/*------------------------------ bl_listSweep ----------------------------------
 *    
 * @brief 	cleans the list of all unlinked elements,
 *              implicitly sorts the nodes
 * @author 	Christian Otto
 *   
 */
void bl_listSweep(List *l){
  Uint cur, last = 0;
  Listelem *bufnodes;
  char *bufdata, *p;
  p = (char *) l->data;
  bufnodes = (Listelem *) malloc(sizeof(Listelem) * (l->numofelem + BASEINC));
  if (bufnodes == NULL){
    DBG("list.c: Memory allocation for nodes in sweep failed. Exit forced.\n",
        NULL);
    exit(-1);
  }
  bufdata = (char *) malloc(l->sizeofelem * (l->numofelem + BASEINC));
  if (bufdata == NULL){
    DBG("list.c: Memory allocation for data in sweep failed. Exit forced.\n",
        NULL);
    exit(-1);
  }
  for(cur = l->first; cur != -1; cur = l->nodes[cur].next){
    bufnodes[last].prev = last - 1;
    if (l->nodes[cur].next != -1){
      bufnodes[last].next = last + 1;
    }
    else {
      bufnodes[last].next = -1;
    }
    memmove(bufdata + (last * l->sizeofelem), p + (cur * l->sizeofelem),
        l->sizeofelem);
    last++;
  }
  free(l->nodes);
  free(l->data);
  l->nodes = bufnodes;
  l->data = bufdata;
  if(l->numofelem)
    l->first = 0;
  else
    l->first = -1;
  l->last = l->numofelem - 1;
  l->allocelem = l->numofelem + BASEINC;
  l->nextfree = l->numofelem;
}

/*------------------------------ bl_listSize -----------------------------------
 *    
 * @brief 	returns number of elements in the list
 * @author 	Christian Otto
 *   
 */
Uint bl_listSize(List *l){
  return l->numofelem;
}


/*--------------------------- bl_listSearchInsert ----------------------------
 *    
 * @brief searches the list position if srch(cur, elem, nfo) returns 0.
 * @author Steve Hoffmann 
 *   
 */
 
Uint
bl_listBinarySearchInsert (List *l, void *elem, 
    Uint (*cmp)(Uint, void *, void*, void*), void *nfo)
{

  Uint cur, left, i=0;

  if(l->numofelem ==0) {
    bl_listInsert(l, -1, elem);
    return 0;
  }

  left = binarySearch_left(l, l->numofelem, elem, cmp, nfo);
  left = (left > 0) ? left - 1: 0;
  cur = l->first;

  while(i < left && cur != l->last) {
    cur = l->nodes[cur].next;
    i++;
  }


  if(cmp(left, l, elem, nfo) < 2 && cur == l->first){
    bl_listInsert(l, -1, elem);
    return 0;
  }

  if(cmp(left, l, elem, nfo) == 2 && cur == l->last){
    bl_listInsert(l, l->last, elem);
    return l->last+1;
  }


  if(cmp(left, l, elem, nfo) == 2 && 
      (cmp(left+1, l, elem, nfo) == 1 || cmp(left+1, l, elem, nfo) == 0)){
    bl_listInsert(l, cur, elem);
    return cur+1;
  }

  
  return -1;
}

#ifdef LISTTEST

Uint
cmp_listtestobj(Uint no, void *list, void *elem, void *nfo) {
  List *l;
  listtestobj_t *a, *b;

  l = (List*) list;
  a = (listtestobj_t*) bl_listGetElem (l, bl_listGetCur(l,no));
  b = (listtestobj_t*) elem;

  if(a->unsigned1 > b->unsigned1) return 1;
  if(a->unsigned1 < b->unsigned1) return 2;
  if(a->unsigned2 > b->unsigned2) return 1;
  if(a->unsigned2 < b->unsigned2) return 2;


  return 0;
}

void
rmv_listtestobj(void *elem) {
  listtestobj_t *o = (listtestobj_t*) elem;
  free(o->string);
}

int
main(int argv, char **argc) {
  List l, sl;
  listtestobj_t obj, *ret, *ret2;
  int i;
  Uint last=0;

  srand(time(NULL)); 

  bl_listInit(&l, 1000, sizeof(listtestobj_t));
  obj.unsigned1 = 1; obj.unsigned2 = 1;
  bl_listInsert(&l, l.last, &obj);
  obj.unsigned1 = 1; obj.unsigned2 = 100;
  bl_listInsert(&l, l.last, &obj);
  obj.unsigned1 = 1; obj.unsigned2 = 1000;
  bl_listInsert(&l, l.last, &obj);
  obj.unsigned1 = 1; obj.unsigned2 = 10000;
  bl_listInsert(&l, l.last, &obj);
  obj.unsigned1 = 1; obj.unsigned2 = 100000;
  bl_listInsert(&l, l.last, &obj);
  obj.unsigned1 = 1; obj.unsigned2 = 100000;
  bl_listInsert(&l, l.last, &obj);


  for(i=0; i < l.numofelem; i++) {
    ret = (listtestobj_t*) bl_listGetElem (&l, bl_listGetCur(&l, i));
    printf("elem: %d, value:%d\n", i, ret->unsigned2);
  }

  bl_listInit(&sl, 1000, sizeof(listtestobj_t)); 
   
  obj.unsigned1 = 4; obj.unsigned2 = 4603315; obj.integer=1; 
    obj.string = malloc(sizeof(char)*10);
    bl_listBinarySearchInsert (&sl, &obj, cmp_listtestobj, NULL);

    obj.unsigned1 = 4; obj.unsigned2 = 4604837; obj.integer=2; 
    obj.string = malloc(sizeof(char)*10);
    bl_listBinarySearchInsert (&sl, &obj, cmp_listtestobj, NULL);

  obj.unsigned1 = 4; obj.unsigned2 = 1822274; obj.integer=3; 
    obj.string = malloc(sizeof(char)*10);
    bl_listBinarySearchInsert (&sl, &obj, cmp_listtestobj, NULL);


 

/*  for(i=0; i < 5000; i++) { 
    if(i % 1000 == 0) printf("%d\n",i);
    obj.unsigned1 = 0; obj.unsigned2 = rand() % 5000; obj.integer=i; 
    obj.string = malloc(sizeof(char)*10);
    bl_listBinarySearchInsert (&sl, &obj, cmp_listtestobj, NULL);
  } 
*/
  fprintf(stdout, "sl has %d elems\n", sl.numofelem);
/*
  for(i=0; i < 100; i++) {
    last = rand() % 3000;
    printf("delete elem i= %d\n", last);
    ret = bl_listUnlink(&sl, last, NULL);
    if(ret) { 
        free(ret->string);
        free(ret);
    }
  }

  last = 0;
  bl_listSweep(&sl);
  fprintf(stdout, "sl has %d elems\n", sl.numofelem);
*/

  for(i=0; i < sl.numofelem; i++) {
    ret = (listtestobj_t*) bl_listGetElem (&sl, bl_listGetCur(&sl, i));
    ret2 = (listtestobj_t*) bl_listGetElem (&sl, i);

    if(ret) { 
        assert(last <= ret->unsigned2);
        printf("elem: %d, value:%d, rank:%d\n", i, ret->unsigned2, ret->integer);
        printf("elem: %d, value:%d, rank:%d\n", i, ret2->unsigned2, ret2->integer);
        last = ret->unsigned2;
    } else {
      printf("elem %d is empty\n",i);
    }
  }

  bl_listDestruct(&sl, rmv_listtestobj);
  bl_listDestruct(&l, NULL);
}

#endif 
