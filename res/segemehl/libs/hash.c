/**
 * hash.c
 * implementation of simple hashing
 *
 * @author Christian Otto
 * @email christian@bioinf.uni-leipzig.de
 * @date Wed Nov 24 15:48:15 CET 2010
 *
 */

/*
 * SVN
 * Revision of last commit: $Rev: 284 $
 * Author: $Author: steve $
 * Date: $Date: 2011-05-03 07:41:30 -0400 (Tue, 03 May 2011) $
 * Id: $Id: hash.c 284 2011-05-03 11:41:30Z steve $
 * Url: $URL: http://www2.bioinf.uni-leipzig.de/svn5/segemehl/libs/hash.c $
 */

#include <stdlib.h>
#include <stdio.h>
#include "basic-types.h"
#include "debug.h"
#include "info.h"
#include "biofiles.h"
#include "hash.h"

#ifdef HASHING

hash_t DJBHash(char* msg){
  Uint i, len;
  hash_t hash = 5381;

  len = strlen(msg);
  for(i = 0; i < len; msg++, i++){
    hash = ((hash << 5) + hash) + (*msg);
  }

  return hash;
}

hash_t* DJBHash2(char* msg, Uint len, Uint masksize){
  Uint i, j;
  hash_t *hash;
  hash = malloc(sizeof(hash_t) * masksize);
  for (j = 0; j < masksize; j++){
    hash[j] = 5381;
  }
  
  for(i = 0; i < len; msg++, i++){
    for (j = 0; j < masksize; j++){
      if (j == 0 || i & (1 << (j-1))){
	hash[j] = ((hash[j] << 5) + hash[j]) + (*msg);
	//DBG("%d %d\n", i, j);
      }
    }
  }
  return hash;
}

void bl_hashInit(Hash *h, Uint hashbitsize, size_t sizeofelem,
		 hash_t (*hashfunc)(char *)){
  if (hashbitsize < 0 || hashbitsize >= PRIMES_SIZE){
    DBG("hash.c: Attempt to initialize a simple hash of bit size %d but \
must be between 1 and 32. Exit forced.\n", hashbitsize);
    exit(-1);
  }
  if (sizeofelem <= 0){
    DBG("hash.c: Attempt to initialize a simple hash with data of size %d.\
Exit forced.\n", sizeofelem);
    exit(-1);
  }
  h->allocelem = PRIMES[hashbitsize - 1];
  h->sizeofelem = sizeofelem;
  h->numofelem = 0;
  h->hashfunc = hashfunc;

  h->hashspace = malloc(h->allocelem * h->sizeofelem);
  h->flag = malloc(h->allocelem * sizeof(char));
  if (h->hashspace == NULL || h->flag == NULL){
    DBG("hash.c: Memory allocation failed. Exit forced.\n", NULL);
    exit(-1);
  }
  memset(h->hashspace, 0, h->allocelem * h->sizeofelem);
  memset(h->flag, EMPTY, h->allocelem * sizeof(char));
}

void bl_hashDestruct(Hash *h, void (*rmv)(void *)){
  hash_t i;
  char *p;
  if (rmv != NULL){
    p = (char *) h->hashspace;
    for (i = 0; i < h->numofelem; i++){
      rmv(p + (i * h->sizeofelem));
    }
  }
  free(h->hashspace);
  free(h->flag);
  h->allocelem = 0;
  h->sizeofelem = 0;
  h->numofelem = 0;
}

hash_t bl_hashGetHashval(Hash *h, char *key){
  assert(key != NULL);
  hash_t hashval;
  
  hashval = (h->hashfunc)(key) % h->allocelem;
  return hashval;
}

hash_t bl_hashGetHashinc(Hash *h, char *key){
  assert(key != NULL);
  hash_t hashval;
  
  hashval = 1 + ((h->hashfunc)(key) % (h->allocelem - 1));
  return hashval;
}

unsigned char bl_hashGetFlagFromKey(Hash *h, char *key){
  return bl_hashGetFlag(h, bl_hashGetHashval(h, key));
}

unsigned char bl_hashGetFlag(Hash *h, hash_t hashval){
  return h->flag[hashval];
}

void *bl_hashGetDataFromKey(Hash *h, char *key){
  return bl_hashGetData(h, bl_hashGetHashval(h, key));
}

void *bl_hashGetData(Hash *h, hash_t hashval){
  char *p;
  if (h->numofelem == 0 || bl_hashGetFlag(h, hashval) == EMPTY){
    return NULL;
  }
  p = (char *) h->hashspace;
  return(p + (hashval * h->sizeofelem));
}

unsigned char bl_hashInsertFromKey(Hash *h, char *key, void *data){
  return bl_hashInsert(h, bl_hashGetHashval(h, key), data);
}

unsigned char bl_hashInsert(Hash *h, hash_t hashval, void *data){
  char *p;  
  if (bl_hashGetFlag(h, hashval) != EMPTY){
    return 0;
  }
  p = (char *) h->hashspace;
  memmove(p + (hashval * h->sizeofelem), data, h->sizeofelem);
  h->flag[hashval] = OCCUPIED;
  h->numofelem++;
  return 1;
}

/* 
 * used for non-indexed fasta files:
 * 
 * --> very performant with few overhead (only hashtable)
 */
void bl_fastxGetTags(void *space, fasta_t *f){
  Uint i, hashbitsize = 28, readlen, *data;
  hash_t dupcnt = 0, collcnt = 0, cnt = 0;
  Hash *hash;
  hash_t hashval, hashinc;
  char *read;

  hash = malloc(sizeof(Hash));
  bl_hashInit(hash, hashbitsize, sizeof(Uint), DJBHash);
  
  for (i = 0; i < f->noofseqs; i++){
    readlen = bl_fastaGetSequenceLength(f, i);
    read = malloc(readlen + 1);
    memmove(read, bl_fastaGetSequence(f, i), readlen+1);
    hashval = bl_hashGetHashval(hash, read);
    hashinc = bl_hashGetHashinc(hash, read);

    /* collision handling */
    while(bl_hashGetFlag(hash, hashval) != EMPTY){
      /*
       * no double hashing if too full
       * --> notice user? abort?
       * but: too full is known before since each read is
       * inserted at most at one position
       * but: amount of double hashing is not known
       */
      collcnt++;
      data = (Uint *) bl_hashGetData(hash, hashval);
      //DBG("coll:%u\t%u\t%s\n%u\t%u\t%s\n", i, hashval, read, *data, bl_hashGetHashval(hash, bl_fastaGetSequence(f, *data)), bl_fastaGetSequence(f, *data));
      /* duplicate */
      if (strcmp(bl_fastaGetSequence(f, *data), read) == 0){
	break;
      }
      /* double hashing */
      if (hashval + hashinc >= hash->allocelem){
	hashval = hashval + hashinc - hash->allocelem;
      }
      else {
	hashval = hashval + hashinc;
      } 
    }
    /* empty position --> insert into hash */
    if (bl_hashGetFlag(hash, hashval) == EMPTY){
      if (!bl_hashInsert(hash, hashval, &i)){
	DBG("Insert in hash failed. Exit forced.\n", NULL);
	exit(-1);
      }
      cnt++;
      /* store that read is unique up to now */
      
    }
    /* duplicate handling */
    else {
      dupcnt++;
    }
    free(read);
  }
  NFO("%u/%u unique reads and %u duplicates and %u collisions\n", hash->numofelem, cnt, dupcnt, collcnt);
  bl_hashDestruct(hash, NULL);
  free(hash);
}

/* 
 * used for indexed fasta files:
 * stores collisions with read indices and one read sequence for post-processing,
 * afterwards sorting to access indexed fasta files only block-wise and
 * resolve collisions by one linear scan
 * NOTE: correct number of collisions may differ to non-indexed variant
 */
void bl_fastxGetTags3(void *space, fasta_t *f){
  Uint i, j, hashbitsize = 28, readlen, *data, lastcoll;
  hash_t dupcnt = 0, cnt = 0;
  Hash *hash;
  hash_t hashval;
  char *read;
  collision_t *coll;

  hash = malloc(sizeof(Hash));
  bl_hashInit(hash, hashbitsize, sizeof(Uint), DJBHash);
  lastcoll = 0;

  /* TODO: more efficiently */
  coll = malloc(sizeof(collision_t) * f->noofseqs);
  

  MSG("Hashing\n");
  for (i = 0; i < f->noofseqs; i++){
    //DBG("%u\t%u\n", i, bl_fastaGetQuantity(f,i));
    bl_fastaSetQuantity(f, i, 2);
    //DBG("%u\t%u\n", i, bl_fastaGetQuantity(f,i));
    readlen = bl_fastaGetSequenceLength(f, i);
    read = malloc(readlen + 1);
    memmove(read, bl_fastaGetSequence(f, i), readlen+1);
    hashval = bl_hashGetHashval(hash, read);

    /* collision handling */
    if(bl_hashGetFlag(hash, hashval) != EMPTY){
      data = (Uint *) bl_hashGetData(hash, hashval);
      coll[lastcoll].a = *data;
      coll[lastcoll].b = i;
      coll[lastcoll].readb = read;
      lastcoll++;
    }
    /* empty position --> insert into hash */
    else {
      cnt++;
      /* init entry object */
      if (!bl_hashInsert(hash, hashval, &i)){
	DBG("Insert in hash failed. Exit forced.\n", NULL);
	exit(-1);
      }
      free(read);
    }
  }
  /* sort collisions by a idx and readb sequence */
  MSG("Sorting\n");
  qsort(coll, lastcoll, sizeof(collision_t), cmp_collision_qsort);
  NFO("Resolving %u collisions\n", lastcoll);
  /* resolving collisions */
  for (i = 0; i < lastcoll; i++){
    read = bl_fastaGetSequence(f, coll[i].a);
    readlen = bl_fastaGetSequenceLength(f, coll[i].a);
    /* duplicate with element in hash table */
    if (strncmp(read, coll[i].readb, readlen) == 0){
      bl_fastaSetQuantity(f, coll[i].a, bl_fastaGetQuantity(f, coll[i].a) + 1);
      /* TODO: slow if this needs sequence data reloading */
      //bl_fastaSetQuantity(f, coll[i].b, 0);
      dupcnt++;
      free(coll[i].readb);
    }
    else {
      cnt++;
      /*
       * equal read sequences have
       * collisions with same element
       * in hash table but element has
       * different sequence, e.g.,
       * hashtable[i] = element a with hash(a)=i
       * hash(b) = i but a.sequence != b.sequence
       * hash(c) = i but a.sequence != c.sequence
       * BUT: b.sequence == c.sequence
       * --> process all such cases and jump over
       */
      j = 1;
      while(i + j < lastcoll && coll[i].a == coll[i+j].a &&
	    strcmp(coll[i].readb, coll[i+j].readb) == 0){
	//DBG("coll:%u\t%u\t%s\t%u\t%u\t%s\t%u\n%u\t%u\t%s\t%u\t%u\t%s\t%u\n", i, coll[i].a, bl_fastaGetSequence(f, coll[i].a),
	//    bl_hashGetHashval(hash, bl_fastaGetSequence(f, coll[i].a)), coll[i].b, coll[i].readb, bl_hashGetHashval(hash, coll[i].readb),
	//    i+j, coll[i+j].a, bl_fastaGetSequence(f, coll[i+j].a), bl_hashGetHashval(hash, bl_fastaGetSequence(f, coll[i+j].a)),
	//    coll[i+j].b, coll[i+j].readb, bl_hashGetHashval(hash, coll[i+j].readb));
	dupcnt++;
	/* TODO: slow since it requires sequence data reloading */
	//bl_fastaSetQuantity(f, coll[i].b, bl_fastaGetQuantity(f, coll[i].b) + 1);
	//bl_fastaSetQuantity(f, coll[i+j].b, 0);
	free(coll[i+j].readb);
	j++;
      }
      free(coll[i].readb);

      /* jump over already processed entries */
      i += j - 1;
    }
  }
  free(coll);

  NFO("%u unique reads and %u duplicates and %u collisions\n", cnt, dupcnt, lastcoll);
  bl_hashDestruct(hash, NULL);
  free(hash);
}

void bl_entryDestruct(void *data){
  entry_t *entry = (entry_t *) data;
  free(entry->read);
  entry->idx = 0;
}

int cmp_collision_qsort(const void *a, const void *b){
  collision_t *first = (collision_t *)  a;
  collision_t *second = (collision_t *)  b;
  if (first->a > second->a) return 1;
  if (first->a < second->a) return -1;
  return strcmp(first->readb, second->readb);
}

#endif /* only required if hashing is defined */
