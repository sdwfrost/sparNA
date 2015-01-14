#ifndef HASH_H
#define HASH_H

/**
 * hash.h
 * implementation of simple hashing
 *
 * @author Christian Otto
 * @email christian@bioinf.uni-leipzig.de
 * @date Wed Nov 24 15:48:15 CET 2010
 *
 */

/*
 * SVN
 * Revision of last commit: $Rev: 278 $
 * Author: $Author: steve $
 * Date: $Date: 2011-04-04 11:06:15 -0400 (Mon, 04 Apr 2011) $
 * Id: $Id: hash.h 278 2011-04-04 15:06:15Z steve $
 * Url: $URL: http://www2.bioinf.uni-leipzig.de/svn5/segemehl/libs/hash.h $
 */

#include <stdlib.h>
#include <stdio.h>
#include "basic-types.h"

typedef uint32_t hash_t;

#define MAXFILL 0.77

#define EMPTY    0
#define OCCUPIED 1
#define LOADED   2

#define PRIMES_SIZE 32
static const hash_t PRIMES[PRIMES_SIZE] = 
{
  0ul,          3ul,          11ul,         23ul,         53ul,
  97ul,         193ul,        389ul,        769ul,        1543ul,
  3079ul,       6151ul,       12289ul,      24593ul,      49157ul,
  98317ul,      196613ul,     393241ul,     786433ul,     1572869ul,
  3145739ul,    6291469ul,    12582917ul,   25165843ul,   50331653ul,
  100663319ul,  201326611ul,  402653189ul,  805306457ul,  1610612741ul,
  3221225473ul, 4294967291ul
};

typedef struct
{
  void *hashspace;
  unsigned char *flag;
  hash_t allocelem;
  hash_t numofelem;
  size_t sizeofelem;
  hash_t (*hashfunc) (char *);
} Hash;

typedef struct
{
  Uint a;
  Uint b;
  char *readb; 
} collision_t;

typedef struct
{
  Uint idx;
  char *read; 
} entry_t;

void bl_entryDestruct(void *data);
hash_t DJBHash(char* msg);
int cmp_collision_qsort(const void *a, const void *b);
void bl_hashInit(Hash *h, Uint hashbitsize, size_t sizeofelem, hash_t (*hashfunc)(char*));
void bl_hashDestruct(Hash *h, void (void *));
hash_t bl_hashGetHashval(Hash *h, char *key);
hash_t bl_hashGetHashinc(Hash *h, char *key);
unsigned char bl_hashGetFlagFromKey(Hash *h, char *key);
unsigned char bl_hashGetFlag(Hash *h, hash_t hashval);
void *bl_hashGetDataFromKey(Hash *h, char *key);
void *bl_hashGetData(Hash *h, hash_t hashval);
unsigned char bl_hashInsertFromKey(Hash *h, char *key, void *data);
unsigned char bl_hashInsert(Hash *h, hash_t hashval, void *data);
void bl_hashDestruct(Hash *h, void (*rmv)(void*));
void bl_fastxGetTags(void *space, fasta_t *f);
void bl_fastxGetTags3(void *space, fasta_t *f);

#endif /* HASH_H */
