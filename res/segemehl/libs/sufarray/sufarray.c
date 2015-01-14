
/*
 *  sufarray.c
 *  implementations for enhanced suffix arrays
 *  for large integer alphabets
 *
 *  @author Steve Hoffmann
 *  @email shoffmann@zbh.uni-hamburg.de
 *  @date 12/11/06 14:56:57 CET
 *  
 *  SVN
 *  Revision of last commit: $Rev: 74 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-10-29 15:03:04 +0100 (Wed, 29 Oct 2008) $
 *
 *  Id: $Id: sufarray.c 74 2008-10-29 14:03:04Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/sufarray/sufarray.c $
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/fcntl.h>
#include <sys/mman.h>
#include <unistd.h>
#include "basic-types.h"
#include "memory.h"
#include "mathematics.h"
#include "sufarray.h"
#include "charsequence.h"
#include "falphabet.h"
#include "stack.h"
#include "vstack.h"
#include "sort.h"
#include "container.h"
#include "vtprogressbar.h"
#include "aluruSort.h"
#include "vqueue.h"
#include "debug.h"
#include "md5.h"
#include "info.h"
#include "stringutils.h"
#include "bitArray.h"

unsigned char sl_diskacc = 0;

void
destructinterval(void *space, void *data) {
  FREEMEMORY(space, data);
}


/*------------------------------ checksuflinks -------------------------------
 *    
 * @brief integrity check for suflinks
 * @author Steve Hoffmann 
 *   
 */
 
void checksuflinks(Suffixarray *s, Uint i, Uint j){
  Uint k, childlcp, suflcp, *space = NULL;
  PairUint* child, childsuf;
  Container *children;
  // ignore singletons as initial input
  if (i == j){
    return;
  }
  children = getChildintervals(space, s, i, j, 0);
  for (k = 0; k < bl_containerSize(children); k++){
    child = (PairUint *) bl_containerGet(children, k);
    // exclude singletons
    if (child->a == child->b){
      return;
    }
    // check suflink of child
    childlcp = getlcpval(s, child->a, child->b);
    childsuf = getSuflink(s, child->a, child->b);
    suflcp = getlcpval(s, childsuf.a, childsuf.b);
    if (childlcp != suflcp + 1){
      DBG("suf[%u, %u, %u]=[%u, %u, %u]\n", child->a, child->b, childlcp,
	  childsuf.a, childsuf.b, suflcp);
    }
    // recursively check all children of child
    checksuflinks(s, child->a, child->b);
  }
  bl_containerDestruct(children, NULL);
  free(children);
}

/* ------------------------------ cmpCharSequence ----------------------------
 *    
 * function to compare CharSequences for mulitkey sort (sort.c)
 * 
 */

  Uint
cmpCharSequence (Uint a, Uint b, Uint depth, void *data, void *info)
{
  char *s = (char*) data;	
  Uint *end;

  /*quick fix to meet end of multiintsequence criterion*/
  if (info == NULL) {
    if(s[b] == (char) 127) {
      if (s[a+depth] == (char) 127) {
        return 0;
      }
      return 1;
    }
  } else {
    end = (Uint*) info;
    if (*end == b) {
      if (s[a+depth] == (char) 127) {
        return 0;
      }
      return 1;
    }
  }


  /*real comparison*/
  if (s[a+depth] > s[b+depth]) return 1;
  if (s[a+depth] < s[b+depth]) return 2;

  return 0;
}



/* ---------------------------- constructSufArr -----------------------------
 *    
 * constructs a suffix array from an (unsigned) integer sequence
 * should be working in O(n). It uses linear sorting method
 * introduced by Aluru et al.
 * 
 */

  Suffixarray*
constructSufArr(void *space, 
    CharSequence **s, 
    Uint len, 
    FAlphabet* alphabet,
    unsigned char silent)
{

  Uint i, numofsuffixes,
  *sorted, 
  *inv_suftab;
  MultiCharSeq *mseq; 
  Suffixarray *arr;
  unsigned char *temp,
                *mdfive=NULL;


  
  mseq = concatCharSequences(space, s, len, (char)126, (char)127);
  numofsuffixes = mseq->totallength;
  mdfive  = ALLOCMEMORY(space, NULL, char, 16);
  temp = MD5R((unsigned char*)mseq->sequences, numofsuffixes, NULL);
  
  
  memmove(mdfive, temp, 16);


  if(!silent) NFO("alphabet of size (%d): %s\n", mseq->mapsize, mseq->map);
  if(!silent) NFO("size of db sequence: %u\n", numofsuffixes);
  inv_suftab = ALLOCMEMORY(space, NULL, Uint , numofsuffixes);
  arr = ALLOCMEMORY(space, NULL, Suffixarray, 1);

  if(!silent) MSG("constructing suftab.\n");

#ifdef SUF_MKQUICKSORT
  sorted = quickSortMultikey (space, mseq->sequences, numofsuffixes, 
      cmpCharSequence, numofsuffixes-1, NULL);     
#else
  sorted = alurusort(space, mseq->sequences, 
      &(numofsuffixes));
#endif

  if (!silent)NFO("constructing inv_suftab (%u).\n", numofsuffixes);
  for (i=0; i < numofsuffixes; i++) {
    if (sorted[i] > numofsuffixes) fprintf(stderr, "construction error? %u: %u\n",i, sorted[i]);
    inv_suftab[sorted[i]]=i;
  }
  if (!silent)MSG("inv_suftab constructed.\n");

  arr->seq = mseq;
  arr->numofsuffixes = numofsuffixes;
  arr->suftab = sorted;
  arr->inv_suftab = inv_suftab;
  arr->mdfive = mdfive;
  arr->lcpctab = NULL;
  arr->llvtab = NULL;
  
  arr->id = NULL;
  arr->idvtab = NULL;
  arr->chldtab = NULL;
  arr->bcktab = NULL;

  arr->suflink = NULL;
  arr->suflink_l = NULL;
  arr->suflink_r = NULL;
  arr->llint = 1;
  return arr;
}


void
writeSuffixarray(Suffixarray *s, char *filename) {
  FILE *fp; 
  Uint   nmemb,
         idvmemb,
         llvmemb;
  unsigned char flags = 0;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    DBG("Couldn't open file %s. Exit forced.\n", filename);
    exit(-1);
  }

  if (s->lcpctab != NULL) {
    flags |= LCP_TAB_STORED;
  }

  if (s->chldtab != NULL) {
    flags |= CHLD_TAB_STORED; 
  }

  if(s->suflink != NULL) {
    flags |= SUFLINK_TAB_STORED;
    flags |= SUFLINK_COMPRESSED;
    if(s->llint) flags |= LINT_SUFLINKS;
  }

  if(s->suflink_l != NULL) {
    flags |= SUFLINK_TAB_STORED;
  }

  if(s->mdfive != NULL) {
    flags |= MD5_STORED;
  }

  nmemb = s->numofsuffixes;
  fwrite(&nmemb, sizeof(Uint), 1, fp);
  fwrite(s->suftab, sizeof(Uint), nmemb, fp);
  fwrite(&flags, sizeof(char), 1, fp);

  if (s->lcpctab != NULL) {
    fwrite(s->lcpctab, sizeof(char), nmemb, fp);
    llvmemb = (Uint) s->llvcnt;
    fwrite(&llvmemb, sizeof(Uint), 1, fp);
    fwrite(s->llvtab, 2*sizeof(Uint), llvmemb, fp);
  }

  if (s->chldtab != NULL) {
    fwrite(s->chldtab, sizeof(Uint), nmemb, fp);
  }
  
  if (s->suflink != NULL) {
    fwrite(s->suflink, sizeof(Uint), nmemb, fp);
    fwrite(s->id, sizeof(char), nmemb, fp);
    idvmemb = (Uint) s->idvcnt;
    fwrite(&idvmemb, sizeof(Uint), 1, fp);
    fwrite(s->idvtab, sizeof(PairLSint), idvmemb, fp);
  }

  if (s->mdfive != NULL) {
    fwrite(s->mdfive, sizeof(char), 16, fp);
  }

  fclose(fp);
}


Suffixarray *
readSuffixarray(void *space, 
    char *idxfilename, 
    CharSequence **seqs,
    Uint len,
    unsigned char silent) {
  FILE *fp; 
  Uint     nmemb = 0,
           idvmemb = 0,
           llvmemb = 0,
           numofsuffixes,
           *suftab = NULL,
           idvi =0;
  childtab *chldtab = NULL; 
  unsigned char flags=0,
                *lcpctab = NULL;
  unsigned char *mdfive=NULL,
                *check=NULL;
  PairUint *llvtab = NULL;
  PairLSint *idvtab = NULL;
  PairSint *idvutab = NULL;

  MultiCharSeq *mseq;
  Suffixarray *s;

#ifdef SUFLINK_MMAP
  int fd;
  signed char   *id = NULL;
  long curiopos, offset;
  struct stat sb;
  char *suflinkptr;
  int pagediff_id;
  int pagediff_sl;
#elif SUFLINK_DISKACC
  int fd;
  off_t off_sl;
  off_t off_id;
#else
  signed char   *id = NULL;
  Uint *suflink = NULL;
#endif
  
  mseq = concatCharSequences(space, seqs, len, (char)126, (char)127);
  numofsuffixes = mseq->totallength; 

  fp = fopen(idxfilename, "r");
  if (fp == NULL) {
    DBG("Couldn't open file '%s'. Exit forced.\n", idxfilename);
    exit(-1);
  }

  fread(&nmemb, sizeof(Uint), 1, fp);
  suftab = ALLOCMEMORY(NULL, NULL, Uint, nmemb);
  fread(suftab, sizeof(Uint), nmemb, fp);
  fread(&flags, sizeof(char), 1, fp);

  if (flags & LCP_TAB_STORED) {
    if (!silent) MSG("reading lcpc/vtab.\n");
    lcpctab = ALLOCMEMORY(space, NULL, unsigned char, nmemb);
    fread(lcpctab, sizeof(unsigned char), nmemb, fp);

    fread(&llvmemb, sizeof(Uint), 1, fp);
    llvtab = ALLOCMEMORY(space, NULL, PairUint, nmemb);
    fread(llvtab, sizeof(PairUint), llvmemb, fp);
  }

  if (flags & CHLD_TAB_STORED) {
    if(!silent) MSG("reading childtab.\n");
    chldtab = ALLOCMEMORY(space, NULL, childtab, nmemb);
    fread(chldtab, sizeof(childtab), nmemb, fp);
  }

  if ((flags & SUFLINK_TAB_STORED)) {
    if(!silent) MSG("reading suflinks.\n");

#ifdef SUFLINK_MMAP 
    curiopos = ftell(fp);
    fd = open(idxfilename, O_RDONLY);
    if (fd == -1) {
      perror("open");
      exit(EXIT_FAILURE);
    }

    if (fstat(fd, &sb) == -1) {       
      perror("fstat");
      exit(EXIT_FAILURE);
    }

    offset = curiopos & ~(sysconf(_SC_PAGE_SIZE) - 1);
    if (curiopos >= sb.st_size) {
      fprintf(stderr, "offset is past end of file\n");
      exit(EXIT_FAILURE);
    }
    
    pagediff_sl = curiopos - offset;   
    suflinkptr = mmap(0, nmemb*sizeof(Uint) + pagediff_sl, PROT_READ, MAP_SHARED, fd, offset);

    if (suflinkptr == MAP_FAILED) {
      perror("mmap");
      exit(EXIT_FAILURE);
    }
#elif SUFLINK_DISKACC 
    sl_diskacc = 1;
    off_sl = ftell(fp);
    fd = open(idxfilename, O_RDONLY);
#else
    suflink = ALLOCMEMORY(space, NULL, Uint, nmemb);
    fread(suflink, sizeof(Uint), nmemb, fp);
#endif

#ifdef SUFLINK_MMAP
    offset = (curiopos+(nmemb*sizeof(Uint))) & ~(sysconf(_SC_PAGE_SIZE) - 1);
    if (curiopos >= sb.st_size) {
      fprintf(stderr, "offset is past end of file\n");
      exit(EXIT_FAILURE);
    }
    
    pagediff_id = (curiopos+(nmemb*sizeof(Uint))) - offset;   
    id = mmap(0, nmemb*sizeof(signed char) + pagediff_id, PROT_READ, MAP_SHARED, fd, offset);

    if (id == MAP_FAILED) {
      perror("mmap");
      exit(EXIT_FAILURE);
    }
    fseek(fp, nmemb*(sizeof(Uint)+sizeof(signed char)), SEEK_CUR); 

#elif SUFLINK_DISKACC
    off_id = off_sl+(nmemb*sizeof(Uint));
    fseek(fp, nmemb*(sizeof(Uint)+sizeof(signed char)), SEEK_CUR);    
#else   
    id = ALLOCMEMORY(space, NULL, signed char, nmemb);
    fread(id, sizeof(signed char), nmemb, fp);
#endif

    fread(&idvmemb, sizeof(Uint), 1, fp);
    idvtab = ALLOCMEMORY(space, NULL, PairLSint, idvmemb);
    if ((flags & LINT_SUFLINKS)) {
      if(!silent) MSG("reading lsint id.\n");
      fread(idvtab, sizeof(PairLSint), idvmemb, fp);
    } else { 
      idvutab = ALLOCMEMORY(space, NULL, PairSint, idvmemb);
      if(!silent) MSG("reading uint id.\n");
      fread(idvutab, sizeof(PairUint), idvmemb, fp);
      for(idvi=0; idvi < idvmemb; idvi++) {
        idvtab[idvi].a = idvutab[idvi].a;
        idvtab[idvi].b = idvutab[idvi].b;
      }
      free(idvutab);
    }
  }

  if ((flags & MD5_STORED)) {
    mdfive = ALLOCMEMORY(space, NULL, unsigned char, 16);
    fread(mdfive, sizeof(unsigned char), 16, fp);
  }

  s = ALLOCMEMORY(space, NULL, Suffixarray, 1);
        
  if ((flags & LINT_SUFLINKS)) 
  s->llint = 1; else s->llint=0;
  s->suftab = suftab;
  s->seq = mseq;
  s->numofsuffixes = numofsuffixes;
  s->lcpctab = lcpctab;
  s->llvtab = llvtab;
  s->llvcnt = llvmemb;
  s->inv_suftab=NULL;
  s->chldtab = chldtab;

#ifdef SUFLINK_MMAP
  s->suflink = (Uint*) &suflinkptr[pagediff_sl];
  s->id = &id[pagediff_id];
  s->pagediff_id = pagediff_id;
  s->pagediff_sl = pagediff_sl; 
#elif SUFLINK_DISKACC
  s->fd = fd;
  s->off_sl = off_sl;
  s->off_id = off_id;
#else
  s->suflink = suflink;
  s->id = id;
#endif

  s->idvtab = idvtab;
  s->idvcnt = idvmemb;
  s->mdfive = mdfive;

  if (!silent) NFO("read suffix array '%s' with %u elements.\n", 
      idxfilename, nmemb);

  fclose(fp);

  check = MD5R((unsigned char*)mseq->sequences, numofsuffixes, NULL); 
  //Uint di;
  //fprintf(stderr, "fasta:");
  //for(di=0; di < 16; di++) {
  //  fprintf(stderr, "%02x", check[di]);
  //}
  //fprintf(stderr, "\nindex:");
  //for(di=0; di < 16; di++) {
  //  fprintf(stderr, "%02x", mdfive[di]);
  //}
  //fprintf(stderr, "\n");

  if (mdfive == NULL) {
    MSG("warning: index does not contain md5 key.\n");
  } else {
    if(checkmd5(check, mdfive) != 0) {
      MSG("error: db and idx MD5 mismatch. Wrong db?\n");
      char op = 0;
      while(op != 'i' && op != 'u' && op != 'a') {
	MSG("options: (i)gnore  (u)pdate index file  (a)bort: ");
	do {
	  op = fgetc(stdin);
	} while(ISWHITESPACE(op));
      }
      if (op == 'u'){
	NFO("updating suffix array '%s' on disk.\n", idxfilename);
	s->mdfive = check;
	writeSuffixarray(s, idxfilename);
      }
      if (op == 'a'){
	exit(-1);
      }
    } else {
      MSG("md5 keys of index and db match.\n");
    }
  }    
  FREEMEMORY(space, check);

  return s;
}



/*------------------------------ destructSufArr ------------------------------
 *    
 * destruct a suffix array.
 * 
 */
  void
destructSufArr (void *space, Suffixarray *arr)
{
  FREEMEMORY(space, arr->suftab);
  if (arr->lcpctab != NULL)
    FREEMEMORY(space, arr->lcpctab);
  if (arr->inv_suftab != NULL)
    FREEMEMORY(space, arr->inv_suftab);
  if (arr->seq != NULL)
    destructMultiCharSeq(space, arr->seq);
  if (arr->idvtab != NULL) 
    FREEMEMORY(space, arr->idvtab);
  if (arr->suflink != NULL)
#ifdef SUFLINK_MMAP
    munmap(arr->id, arr->numofsuffixes*sizeof(Uint) + arr->pagediff_sl);
#else 
    FREEMEMORY(space, arr->suflink);
#endif
  if (arr->chldtab != NULL)
    FREEMEMORY(space, arr->chldtab);
  if (arr->mdfive)
    FREEMEMORY(space, arr->mdfive);
  if(arr->llvtab != NULL)
    FREEMEMORY(space, arr->llvtab);
  if(arr->id != NULL)
#ifdef SUFLINK_MMAP
    munmap(arr->id, arr->numofsuffixes*sizeof(signed char) + arr->pagediff_id);
#else
    FREEMEMORY(space, arr->id);
#endif
  FREEMEMORY(space, arr);
#ifdef SUFLINK_DISKACC
  close(arr->fd);
#endif
  return ;
}
inline Uint
lcp(Suffixarray *s, Uint i) {
  PairUint *ret;
  Uint val;

  /*  return s->lcptab[i];*/

  if(s->lcpctab[i] < 254) {
    return (Uint) s->lcpctab[i];
  } else { 
    val = i;
    ret=bsearch(&val, s->llvtab, s->llvcnt, sizeof(PairUint), cmp_PairUint_bsearch);
  }
  if (ret == NULL) {
    DBG("lcp '%d' not found. Exit forced.\n", i);
    exit(-1);
  }

  return ret->b;
}


inline Uint
maxlcp(Suffixarray *s) {
  Uint i;
  Uint max = 0;

  for(i=0; i < s->numofsuffixes; i++) {
    if (lcp(s,i) > max) max = lcp(s,i);
  }
  return max;
}

/*------------------------------ computeLcpTab -------------------------------
 *    
 * computes the lcp tab from suftab and inv_suftab in O(n).
 * 
 */

  void
constructLcp (void *space, Suffixarray *arr)
{
  Uint i, j, k;
  Uint max = 0;
  Lint l=0;

  /*arr->lcptab = ALLOCMEMORY(space, NULL, Uint, arr->numofsuffixes);
    memset(arr->lcptab, 0, sizeof(Uint)*arr->numofsuffixes);*/

  arr->lcpctab = ALLOCMEMORY(space, NULL, unsigned char, arr->numofsuffixes);
  arr->llvcnt = 0;
  arr->llvtab = NULL;

  initProgressBarVT();
  for(i=0; i < arr->numofsuffixes; i++) {    
    j = arr->inv_suftab[i];

    if (j > 0) {
     k = arr->suftab[j-1];
      l=l-1;
      if (l < 0) l=0;

      while (arr->seq->sequences[i+l] == arr->seq->sequences[k+l] 
          //            && arr->seq->sequences[i+l] != (char)126 
          //            && arr->seq->sequences[k+l] != (char)126
          ){ 
        l++;
      }

      /*    arr->lcptab[j] = l;*/
      if (l > max) max = l;
      if (l < 254) {
        arr->lcpctab[j] = (char) l;
      } else {
        arr->lcpctab[j] = 254;
        arr->llvtab = ALLOCMEMORY(space, arr->llvtab, PairUint, arr->llvcnt+1);
        arr->llvtab[arr->llvcnt].a = j;
        arr->llvtab[arr->llvcnt].b = l;
        arr->llvcnt++;
      }
    }
  }

  qsort(arr->llvtab, arr->llvcnt, sizeof(PairUint), cmp_PairUint_qsort);
  arr->maxlcp=max;
  arr->lcpctab[0]=0;

  return;
}

inline Lint
id(Suffixarray *s, Uint i) {
  PairLSint *retl;
  signed char ch;
  Lint vall;
  ssize_t nbytes;

  if(sl_diskacc) {
    lseek(s->fd, s->off_id+(i*sizeof(signed char)), SEEK_SET);
    nbytes = read(s->fd, &ch, sizeof(signed char));

    if(nbytes == -1) {
      perror("suflink access failed");
      exit(EXIT_FAILURE);
    }
  } else {
    ch = s->id[i];
  }

  if(ch != (signed char)-128) {
    return (Lint) ch;
  } else {

      vall = (Lint) i;
      retl = bsearch(&vall, s->idvtab, s->idvcnt, sizeof(PairLSint), 
          cmp_PairLSint_bsearch);
      if (retl == NULL) {
        DBG("id '%d' not found. Exit forced.\n", i);
        exit(-1);
      }
      return retl->b;
  }
  exit(-1);
  return retl->b;
}


inline unsigned char
isnextlIndex(Suffixarray *s, Uint i) {
  return (lcp(s,s->chldtab[i].val) == lcp(s,i));
}

inline unsigned char
isdownIndex(Suffixarray *s, Uint i) {
  return (lcp(s,s->chldtab[i].val) > lcp(s,i));
}

inline unsigned char
isupIndex(Suffixarray *s, Uint i) {
  return (lcp(s,i) > lcp(s, i+1));
}

inline Uint
getfirstlindex(Suffixarray *s, Uint i, Uint j){
  Uint val=0; 

  if((i==0 && j == s->numofsuffixes-1) || i==j) return 0;

  if (j < s->numofsuffixes && isupIndex(s,j) 
      && i < s->chldtab[j].val && j >= s->chldtab[j].val) {
    val = s->chldtab[j].val;
  } else if (isdownIndex(s,i)){
    val = s->chldtab[i].val;
  }

  return val;
}

inline Uint
getlcpval(Suffixarray *s, Uint i, Uint j){
  Uint val=0;

  if((i==0 && j == s->numofsuffixes-1) || i==j) return 0;

  if (j < s->numofsuffixes && isupIndex(s,j) 
      && i < s->chldtab[j].val && j >= s->chldtab[j].val) {
    val = lcp(s, s->chldtab[j].val);
  } else if (isdownIndex(s,i)){
    val = lcp(s, s->chldtab[i].val);
  }

  return val;
}


inline void
addinterval(void *space, Container *c, Uint a, Uint b) {
  PairUint range;
  PairUint *check;

  range.a=a;
  range.b=b;
  bl_containerAdd(c, &range);
  if(bl_containerSize(c) > 0){
    check = (PairUint*) bl_containerGet(c, bl_containerSize(c) - 1);
    if(range.a < check->a) {
      printf("check->a: %d, range.a: %d\n", check->a, range.a);
    }
  }
  return;
}

inline Lint*
getChildintervalsArr(void *space, 
    Suffixarray *s, 
    Uint i, 
    Uint j, 
    Uint *noofintervals,
    BOOL checkdelim) { 

  Lint *c;
  Uint count=0;
  Uint     i1,
           i2,
           lcp =0;
  unsigned char child;

  /*ALERT -1*/
  child = (i > 0 || j < s->numofsuffixes-1);
  c = (Lint*) malloc(sizeof(Lint)*s->seq->mapsize*2+1);
  if(checkdelim) lcp = getlcpval(s, i, j);
  
  if(child) {
    if (i < s->chldtab[j].val && s->chldtab[j].val <=j) {
      i1 = s->chldtab[j].val;
    } else {
      i1 = s->chldtab[i].val;
    }
    if(!checkdelim || s->seq->sequences[s->suftab[i] + lcp] != s->seq->delim) {
      c[count*2] = i;
      c[count*2+1] = i1-1;
      count++;
    }

  } else {
    i1 = i;
  }

  while(i1 < s->numofsuffixes-1 && isnextlIndex(s,i1) && !isupIndex(s,i1) 
      && s->chldtab[i1].val != 0) {
    i2 = s->chldtab[i1].val;
    if(!checkdelim || s->seq->sequences[s->suftab[i1] + lcp] != s->seq->delim) {
      c[count*2] = i1;
      c[count*2+1] = i2-1;
      count++;
    }
    i1 = i2;
  }

  if(child && (!checkdelim || s->seq->sequences[s->suftab[i1] + lcp] != s->seq->delim)) {
    c[count*2] = i1;
    c[count*2+1] = j;
    count++;
  }

  *noofintervals = count;
  return c;
}

inline Container*
getChildintervals(void *space, 
    Suffixarray *s, 
    Uint i, 
    Uint j,
    BOOL checkdelim) { 

  Container *c;
  Uint     i1,
           i2,
           lcp = 0;
  unsigned char child;

  /*ALERT -1*/
  child = (i > 0 || j < s->numofsuffixes-1);
  c = (Container *) malloc(sizeof(Container));
  bl_containerInit(c, 10, sizeof(PairUint));
  if(checkdelim) lcp = getlcpval(s, i, j);

  if(child) {
    if (i < s->chldtab[j].val && s->chldtab[j].val <=j) {
      i1 = s->chldtab[j].val;
    } else {
      i1 = s->chldtab[i].val;
    }
    if(!checkdelim || s->seq->sequences[s->suftab[i]+lcp] != s->seq->delim) {
      addinterval(space, c, i, i1-1);
    }

  } else {
    i1 = i;
  }

  while(i1 < s->numofsuffixes-1 && isnextlIndex(s,i1) && !isupIndex(s,i1) 
      && s->chldtab[i1].val != 0) {
    i2 = s->chldtab[i1].val;
    if(!checkdelim || s->seq->sequences[s->suftab[i1]+lcp] != s->seq->delim) {
      addinterval(space, c, i1, i2-1);
    }
    i1 = i2;
  }

  if(child && (!checkdelim || s->seq->sequences[s->suftab[i1]+lcp] != s->seq->delim)) {
    addinterval(space, c, i1,j);
  }
  return c;
}



inline PairUint
getSuflink(Suffixarray *s, Uint i, Uint j) {
  Uint slidx, base;
  ssize_t nbytes;
  Lint off;
  PairUint link;
  Lint a, b;

  slidx = getfirstlindex(s, i, j);

  if(sl_diskacc) {
    lseek(s->fd, s->off_sl+(slidx*sizeof(Uint)), SEEK_SET);
    nbytes = read(s->fd, &base, sizeof(Uint));
    if(nbytes == -1) {
      perror("suflink access failed");
      exit(EXIT_FAILURE);
    }
  } else {
    base = s->suflink[slidx];
  }

  if ((off=id(s, base)) > 0) {
    a = base;
    b = off + base;
  } else {
    a = base+off;
    b = base;
  }

  link.a = a;
  link.b = b;
  return link;
}

inline PairUint
jumpkSuflinks(Suffixarray *s, Uint i, Uint j, Uint k) {
    Uint v;
    PairUint link;

    link.a = i;
    link.b = j;

    for(v=0; v < k && getlcpval(s, link.a, link.b) > 0; v++) {
        link = getSuflink(s, link.a, link.b);
    }

    return link;
}

/*---------------------------- constructchildtab -----------------------------
 *    
 * @brief performs bottom-up traversals to construct childtable 
 * @author Steve Hoffmann 
 *   
 */
 
void
constructchildtab(void *space, Suffixarray *s) {
  Uint i;
  Lint lastIndex = -1;
  Stack *stack;

  s->chldtab = ALLOCMEMORY(space, NULL, childtab, s->numofsuffixes+1);
  memset(s->chldtab, 0, s->numofsuffixes*sizeof(childtab));
  stack = ALLOCMEMORY(space, NULL, Stack, 1);
  bl_stackInit(stack, 100000);

  bl_stackPush(stack, 0);

  for(i = 0; i < s->numofsuffixes; i++) 
  {
    while(lcp(s,i) < lcp(s, bl_stackTop(stack))) {
      lastIndex = bl_stackPop(stack);
      if(lcp(s,i) <= lcp(s, bl_stackTop(stack)) && 
          lcp(s,bl_stackTop(stack)) != lcp(s,lastIndex))
      {
        s->chldtab[bl_stackTop(stack)].val  = lastIndex;
      }
    }
    if (lastIndex != -1) {
      s->chldtab[i-1].val = lastIndex;
      lastIndex = -1;
    }
    bl_stackPush(stack, i);
  }

  /*construction of nextlIndex value*/
  bl_stackDestruct(stack);
  bl_stackInit(stack, 10000);
  bl_stackPush(stack, 0);

  for(i = 1; i < s->numofsuffixes; i++) {
    while(lcp(s,i) < lcp(s, bl_stackTop(stack))) {
      bl_stackPop(stack);
    }
    if (lcp(s,i) == lcp(s, bl_stackTop(stack))) {
      lastIndex = bl_stackPop(stack);
      s->chldtab[lastIndex].val = i;
    }
    bl_stackPush(stack, i);
  }

  bl_stackDestruct(stack);
  FREEMEMORY(space, stack);
  return;
}



/*------------------------------- computeId ----------------------------------
 *    
 * @brief performs a top down traversal on the tree represented 
 * by the suffix array to compute unique ids for each interval
 * @author Steve Hoffmann 
 *   
 */

void
computeId (void *space, Suffixarray *s) {
  Uint i; 
  Lint l, 
       r;
  Container *c;
  VQueue vqueue;
  PairUint ival;

  bl_vqueueInit(&vqueue, 1000, sizeof(PairUint));
  s->id = ALLOCMEMORY(space, NULL, char, s->numofsuffixes+2);
  memset(s->id, 0, sizeof(char)*s->numofsuffixes+2);

  s->idvtab = ALLOCMEMORY(space, NULL, PairLSint, 1);
  s->idvcnt = 1;

  ival.a = 0;
  ival.b = s->numofsuffixes-1;

  s->id[0] = (signed char) -128;
  s->idvtab[0].a = 0;
  s->idvtab[0].b = (s->numofsuffixes-1);

  bl_vqueueEnqueue(&vqueue, &ival);
  while (!bl_vqueueIsEmpty(&vqueue)){
        
    PairUint *tmp = (PairUint *) bl_vqueueDequeue(&vqueue, NULL);
    memcpy(&ival, tmp, sizeof(PairUint));
    free(tmp);
 
    c = getChildintervals(space, s, ival.a, ival.b, 0);
    for(i=0; i < bl_containerSize(c); i++){

      l = ((PairUint*)bl_containerGet(c, i))->a;
      r = ((PairUint*)bl_containerGet(c, i))->b;         
        
      if (l < r) {
        if(s->id[l] == 0) {
          if (r-l <= 127){
            s->id[l] = (signed char) r-l;
          } else {
            s->id[l] = (signed char)-128;
            s->idvtab = ALLOCMEMORY(space, s->idvtab, PairLSint, s->idvcnt+1);
            s->idvtab[s->idvcnt].a = l;
            s->idvtab[s->idvcnt].b = r-l;
            s->idvcnt = s->idvcnt+1; 
          }

          //        id[l] = r;
        } else if(s->id[r] == 0) {
          
          if(l-r > -128) {
            s->id[r] = (signed char) l-r;
          } else {
            s->id[r] = (signed char)-128;
            s->idvtab = ALLOCMEMORY(space, s->idvtab, PairLSint, s->idvcnt+1);
            s->idvtab[s->idvcnt].a = r;
            s->idvtab[s->idvcnt].b = l-r;
            s->idvcnt = s->idvcnt+1;
          } 
          //        id[r] = -l;
        } else {
          DBG("ID failed id[l]:%d, id[r]:%d\n\n", s->id[l], s->id[r]);
          exit(-1);
        }
	ival.a = l;
	ival.b = r;
	bl_vqueueEnqueue(&vqueue, &ival);
      }
    }
    bl_containerDestruct(c, NULL);
    free(c);
  }
  qsort(s->idvtab, s->idvcnt, sizeof(PairLSint), cmp_PairLSint_qsort);
  bl_vqueueDestruct(&vqueue, NULL);
  return;
}


/*------------------------------- getsuffsucc --------------------------------
 *    
 * @brief performs a bottom-up traversal of the suffix array collecting
 *        cause(p)-successors for suffix link construction
 * @author Steve Hoffmann 
 *   
 */

Uint *
getsufsucc(void *space, Suffixarray *s){
    
  Lint  i,
       lb,
    //llcp,
       llb,
       min1,
       min2,
       m;

  Uint  *A;
  Stack *stack,
        *mstack;
  
  A = ALLOCMEMORY(space, NULL, Uint, s->numofsuffixes+2);
  memset(A, 255, sizeof(Uint)*s->numofsuffixes+2);

  stack = ALLOCMEMORY(space, NULL, Stack, 1);
  mstack = ALLOCMEMORY(space, NULL, Stack, 1);
  
  bl_stackInit(stack, 100000);
  bl_stackInit(mstack, 100000);

  /*push lcp and lbound*/
  bl_stackPush(stack, 0);
  bl_stackPush(stack, 0);


  bl_stackPush(mstack, s->suftab[0]);
  bl_stackPush(mstack, 0);
  
  for(i = 1; i < s->numofsuffixes; i++) {
    lb = i-1;
    
    bl_stackPush(mstack, s->suftab[i-1]);
    bl_stackPush(mstack, i-1);
    
    while (lcp(s,i) < bl_stackTop(stack)) {    
      bl_stackPop(stack);
      //not used: llcp = bl_stackPop(stack);
        llb = bl_stackPop(stack);
        
        /*child interval is given by llcp-[llb, i-1]*/
        /*cycle children here*/
        
        min1 = s->numofsuffixes+1;
        min2 = s->numofsuffixes+1;
        while (!bl_stackIsEmpty(mstack) && llb <= bl_stackTop(mstack) && 
                bl_stackTop(mstack) <= i-1) {
            bl_stackPop(mstack);
            m = bl_stackPop(mstack);
            
            if (m < min1) {
                min2 = min1;
                min1 = m;
            } else {
              if (m < min2 && m != min1)
                min2 = m;
            } 
        }        
        lb = llb;

        bl_stackPush(mstack, min1);
        bl_stackPush(mstack, lb);
        if (id(s, lb) + lb == i-1) {
            A[min2+1] = lb;
        } else {
            A[min2+1] = i-1;
        }
    }

    if(lcp(s,i) > bl_stackTop(stack)){
        bl_stackPush(stack, lb);
        bl_stackPush(stack, lcp(s,i));
    }
  }

  bl_stackDestruct(stack);
  bl_stackDestruct(mstack);
  FREEMEMORY(space, stack);
  FREEMEMORY(space, mstack);
  return A;
}


/*---------------------------- constructsuflinks ----------------------------
 *    
 * @brief performs a top down traversal on the tree represented 
 * by the suffix array
 * @author Steve Hoffmann 
 *   
 */

void
constructsuflinks (void *space, Suffixarray *s, Uint *succ) {
  Uint i, 
       a, b, 
       l, r,
       max=0, pushes=0;
  PairUint ab, lr, *tmp;
  Lint  u, v; 
  Uint d, 
       slidx,
       lidx,
       *B; 
  Container *c;
  VStack *vstack;
  
#ifdef EXPLICITSUFLINKS

  PairUint suflink;
  Lint m, n;
  
  s->suflink_l = ALLOCMEMORY(space, NULL, Uint, s->numofsuffixes+1);
  s->suflink_r = ALLOCMEMORY(space, NULL, Uint, s->numofsuffixes+1);
  memset(s->suflink_l, 0, s->numofsuffixes*sizeof(Uint));
  memset(s->suflink_r, 0, s->numofsuffixes*sizeof(Uint));
#else
  s->suflink = ALLOCMEMORY(space, NULL, Uint, s->numofsuffixes+1); 
  memset(s->suflink, 0, s->numofsuffixes*sizeof(Uint));
#endif

  vstack = (VStack *) malloc(sizeof(VStack));
  bl_vstackInit(vstack, 100000, sizeof(PairUint));
 
  B = ALLOCMEMORY(space, NULL, Uint, s->maxlcp+1);

  ab.a = 0;
  ab.b = s->numofsuffixes-1;
  bl_vstackPush(vstack, &ab);

  while(!bl_vstackIsEmpty(vstack)){
    if(max < bl_vstackSize(vstack)) max = bl_vstackSize(vstack);
    tmp = (PairUint *) bl_vstackPop(vstack, NULL);
    a = tmp->a;
    b = tmp->b;
    free(tmp);

    c = getChildintervals(space, s, a, b, 0);
    d = getlcpval(s, a, b);

    if (id(s, a)+a == b) {
      B[d] = a;
    } else if (a == b+id(s, b)){
      B[d] = b;
    } else {
      DBG("Id failed. id[a]: %d\n", id(s,a));
    }

    for (i=0; i < bl_containerSize(c); i++) {

      lr = *((PairUint *) bl_containerGet(c,i));
      l = lr.a;
      r = lr.b;

      if(l < r) {
        bl_vstackPush(vstack, &lr);
        pushes++;
      } else {

        lidx = s->suftab[l];

        if((succ[lidx] > 0 && succ[lidx] < s->numofsuffixes-1) &&
            (abs(id(s, succ[lidx])) < s->numofsuffixes-1)) {  

          if (id(s, succ[lidx]) > 0) {
            u = succ[lidx];
            v = id(s, succ[lidx]) + succ[lidx];
          } else {
            u = succ[lidx] + id(s, succ[lidx]);
            v = succ[lidx];
          }

          d = getlcpval(s, u, v);
          slidx = getfirstlindex(s, u, v);

#ifdef EXPLICITSUFLINKS 

          if (id(s, B[d-1]) > 0) {
            m = B[d-1];
            n = id(s, B[d-1])+B[d-1];
          } else {
            m = B[d-1] + id(s, B[d-1]);
            n = B[d-1];
          }

          s->suflink_l[slidx] = m;
          s->suflink_r[slidx] = n;
          suflink = getSuflink(s, u, v);
#else
          d = MAX(1,d);
          s->suflink[slidx] = B[d-1];
#endif

#ifdef SUFLINKDEBUG        

          fprintf(stderr, "linking [%d,%d] -> [%d,%d] {%d,%d}\n", 
              u,v,m,n, s->inv_suftab[s->suftab[u]+1], 
              s->inv_suftab[s->suftab[v]+1]); 


          { Lint w;
            for(w=0; w < getlcpval(s,u,v)-1; w++) {
              if(s->seq->sequences[s->suftab[u] + w+1]!= s->seq->sequences[s->suftab[m]+ w] 
                  || getlcpval(s, u, v) != getlcpval(s, m, n)+1){
                DBG("Suffixlink construction failed with %d-[%d,%d] -> %d-[%d,%d]\n", getlcpval(s, u, v), u, v, getlcpval(s, m, n),m, n);
                exit(-1);
              }
              if(s->seq->sequences[s->suftab[v]+w+1]!= s->seq->sequences[s->suftab[n] + w] 
                  || getlcpval(s, u, v) != getlcpval(s, m, n)+1){
                DBG("Suffixlink construction failed with %d-[%d,%d] -> %d-[%d,%d]\n", getlcpval(s, u, v), u, v, getlcpval(s, m, n), m, n);
                exit(-1);
              }
            }
          }
#endif
        }
      }
    }
    bl_containerDestruct(c, NULL);
    free(c);
  }
  bl_vstackDestruct(vstack, NULL);
  free(vstack);
  FREEMEMORY(space, B);
  fprintf(stderr, "suflink construction. pushes: %d, maxstack: %d\n", pushes, max);
  return;
}



Uint
childCount(void *space, 
    Suffixarray *s, 
    Uint i, 
    Uint j) { 

  Uint     childcount=0; 
  Uint     i1,
           i2;
  unsigned char child;

  /*ALERT -1*/
  child = (i > 0 || j < s->numofsuffixes-1);

  if(child) {
    if (i < s->chldtab[j].val && s->chldtab[j].val <=j) {
      i1 = s->chldtab[j].val;
    } else {
      i1 = s->chldtab[i].val;
    }
    childcount++;

  } else {
    i1 = i;
  }

  while(isnextlIndex(s,i1) && !isupIndex(s,i1) 
      && s->chldtab[i1].val != 0) {
    i2 = s->chldtab[i1].val;
    childcount++;
    i1 = i2;
  }

  if(child) {
    childcount++;
  }

  return childcount;
}



/*----------------------------- getCharInterval ------------------------------
 *    
 * @brief  gets a child interval starting with character ch
 *         pos deprecated (write what you like)
 * @return returns an interval [l,r], empty interval l > r
 * @author Steve Hoffmann 
 *   
 */


inline PairUint
getCharIntervalArr(void *space,
    Suffixarray *s,
    Uint i,
    Uint j,
    Uint pos,
    char ch) 
{
  Lint *c;
  Uint count=0;
  Uint lcp=0;
  PairUint lr;

  lr.a = 1;
  lr.b = 0;

  if(i==j) return lr;

  c = getChildintervalsArr(space,s, i, j, &count, 1);
  lcp = getlcpval(s, i, j);

  for(i=0; i < count; i++) {
    if(s->seq->sequences[ s->suftab[c[i*2]] + lcp] == ch){
      lr.a = c[i*2];       
      lr.b = c[i*2+1];
      break;
    }
  }
  
  free(c);
  return lr;
}



/*----------------------------- getCharInterval ------------------------------
 *    
 * @brief  gets a child interval starting with character ch
 *         pos deprecated (write what you like)
 * @return returns an interval [l,r], empty interval l > r
 * @author Steve Hoffmann 
 *   
 */

inline PairUint
getCharInterval(void *space,
    Suffixarray *s,
    Uint i,
    Uint j,
    Uint pos,
    char ch) 
{
  Container *c;
  Uint lcp=0;
  PairUint lr;

  lr.a = 1;
  lr.b = 0;

  if(i==j) return lr;

  c = getChildintervals(space,s, i, j, 1);
  lcp = getlcpval(s, i, j);

  for(i=0; i < bl_containerSize(c); i++) {

    if(s->seq->sequences[ s->suftab[((PairUint*)bl_containerGet(c, i))->a] + lcp] == ch) {
      lr.a = ((PairUint*)bl_containerGet(c, i))->a;       
      lr.b = ((PairUint*)bl_containerGet(c, i))->b;

      break;
    }
  }
  bl_containerDestruct(c, NULL);
  free(c);
  return lr;
}



/*-------------------------------- dumpSufArr --------------------------------
 *    
 * dumps a suffix array to a screen
 * 
 */

  void
dumpSufArr (Suffixarray *arr)
{
  Uint i;

  for(i=0; i < arr->numofsuffixes; i++) {
    printf("%d \t %d \t %d \t %d \t %d \t %s\n", i, 
        arr->suftab[i], 
        lcp(arr,i),
        arr->inv_suftab[i], 
        arr->seq->sequences[arr->suftab[i]],
        &arr->seq->sequences[arr->suftab[i]]);
  }

  return;
}

void
dumplcps(Suffixarray *arr) {
  Uint i, j, s, t;


  for(i=0; i < arr->numofsuffixes; i++) {
    if (lcp(arr,i) > 0) {
      s = &(arr->seq->sequences[arr->suftab[i-1]])-arr->seq->sequences;
      t = &(arr->seq->sequences[arr->suftab[i]])-arr->seq->sequences;
      printf("lcp of suffix %d and %d has length %d\t:\n", i-1, i, lcp(arr,i));
      for(j=0; j <= lcp(arr,i); j++) printf(" %d ", arr->seq->sequences[s+j]);
      printf("\n");
      for(j=0; j <= lcp(arr,i); j++) printf(" %d ", arr->seq->sequences[t+j]);
      printf("\n");
    }
  }
}

void
dumplcptab(Suffixarray *s) {

  Uint i;

  for(i=0; i < s->numofsuffixes; i++) {
    printf("i:%d lcp:%d\n", 
        i, lcp(s,i));
  }

  printf("\n");

}


void
dumpchildtab(Suffixarray *s) {
  Uint i;

  for(i=0; i < s->numofsuffixes; i++) {
    printf("i:%d up:%d, down:%d, nextlIndex:%d := %d\n", 
        i, isnextlIndex(s,i), isdownIndex(s,i), isnextlIndex(s,i), s->chldtab[i].val);
  }

  printf("\n");
}




/*------------------------------- searchSuffix -------------------------------
 *    
 * @brief looking up a suffix
 * @author Steve Hoffmann 
 *   
 */

  PairUint
searchSuffix (void *space, Suffixarray *arr, char *p, Uint len)
{

  PairUint res, cur;
  Uint i = 0, ell;
  char *suf, *q, *qend, *sufend;

  res.a = 1;
  res.b = 0;
  cur.a = 0;
  cur.b = arr->numofsuffixes-1;

  q = p;
  qend = &q[len-1];

  do {

    cur = getCharInterval(space, arr, cur.a, cur.b, 0, *q);
    if (cur.a > cur.b) return res;

    if(cur.a < cur.b) {
      ell = getlcpval(arr, cur.a, cur.b);
    } else {
      ell = len;
    }

    suf = &arr->seq->sequences[arr->suftab[cur.a]];
    sufend = &suf[MIN(ell, len)-1];
    suf = &suf[i];

    while(*q && suf <= sufend && q <= qend && *suf == *q) {
      ++q;
      ++suf;
      ++i;
    }

    if(*q && suf <= sufend && q <= qend && *suf != *q) return res;

  } while (i < len);

  return cur;
}

