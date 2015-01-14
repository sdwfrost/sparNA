 #ifndef INTSEQUENCE_H
 #define INTSEQUENCE_H

/*
 * charsequence.h
 * declaration of char sequence
 * and functions working on it
 *
 * @author Steve Hoffmann
 * @date Mon 27 Nov 2006
 *
 *  SVN
 *  Revision of last commit: $Rev: 87 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-11-20 11:24:26 +0100 (Thu, 20 Nov 2008) $
 *
 *  Id: $Id: charsequence.h 87 2008-11-20 10:24:26Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/sufarray/charsequence.h $
 */

 #include "basic-types.h"
 #include <stdio.h>
 #include <stdlib.h>
 #include "stringutils.h"

 typedef struct {
	Uint descrlen;			
    Uint namelen;
	Uint urllen;
    Uint noofinfo;
	Uint *infolength;
	
    char *description; 		/*a description*/
    char *alphabetname;		/*the name of the corresponding alphabet*/
	char *url;				/*the name of the sequences url*/	
	char *sequence;			/*the sequence itself*/
    char *info;		        /*additional information*/
	Uint length;

    Uint clip5[2];          /*clipping information*/
    Uint clip3[2];          
    
   #ifdef HASHING
   Uint quantity;          /* quantity of equal sequences */
   #endif
    Uint *map;
    Uint mapsize;
	
 } CharSequence;
 
 void destructSequence(void *, CharSequence *);
 CharSequence* initSequence(void *);
 void resetSequence(CharSequence *);
 char* printSequence(void *, CharSequence *, Uint);
 void dumpSequence(CharSequence *s);
 void saveSequence (CharSequence *s, char *filename); 
 CharSequence* loadSequence (void *space, char *filename);
 char * printAlignment (void *, int *, Uint, CharSequence *, CharSequence *, 
	 Uint);
 CharSequence** createSequenceHash(void *, Uint);

 static inline char* charDNAcomplement(void *space, char *s, Uint len) {
    Uint i,k=0;
    char* buffer;
    
    buffer = ALLOCMEMORY(space, NULL, char, len+1);
    for(i=len; i > 0; i--) {
        switch(s[i-1]) {
          case 'a':
            buffer[k] = 't';
            break;
          case 't':
            buffer[k] = 'a';
            break;
          case 'c':
            buffer[k] = 'g';
            break;
          case 'g':
            buffer[k] = 'c';
            break;
          case 'n':
            buffer[k] = 'n';
            break;
          case 'A':
            buffer[k] = 'T';
            break;
          case 'T':
            buffer[k] = 'A';
            break;
          case 'C':
            buffer[k] = 'G';
            break;
          case 'G':
            buffer[k] = 'C';
            break;
          case 'N':
            buffer[k] = 'N';
            break;
          default:
            buffer[k] = s[i-1]; 
            break;
        }
     k++;
    }
    buffer[k] = '\0';
    return buffer;
 }


 static inline char charComplementChar(char ch) {
 
   switch(ch) {
          case 'a':
            return 't';
            break;
          case 't':
            return 'a';
            break;
          case 'c':
            return 'g';
            break;
          case 'g':
            return 'c';
            break;
          case 'n':
            return 'n';
            break;
          case 'A':
            return 'T';
            break;
          case 'T':
            return 'A';
            break;
          case 'C':
            return 'G';
            break;
          case 'G':
            return 'C';
            break;
          case 'N':
            return 'N';
            break;
          default:
            return ch; 
            break;
   }

   return ch;
 }

static inline char* charIUPACcomplement(void *space, char *s, Uint len) {
  Uint i,k=0;
  char* buffer;

  buffer = ALLOCMEMORY(space, NULL, char, len+1);
  for(i=len; i > 0; i--) {
    switch(s[i-1]) {
      case 'a':
        buffer[k] = 't';
        break;
      case 't':
        buffer[k] = 'a';
        break;
      case 'c':
        buffer[k] = 'g';
        break;
      case 'g':
        buffer[k] = 'c';
        break;
      case 'r':
        buffer[k] = 'y';
        break;
      case 'y':
        buffer[k] = 'r';
        break;
      case 's':
        buffer[k] = 's';
        break;
      case 'w':
        buffer[k] = 'w';
        break;
      case 'k':
        buffer[k] = 'm';
        break;
      case 'm':
        buffer[k] = 'k';
        break;
      case 'b':
        buffer[k] = 'v';
        break;
      case 'd':
        buffer[k] = 'h';
        break;
      case 'h':
        buffer[k] = 'd';
        break;
      case 'v':
        buffer[k] = 'b';
        break;
      case 'n':
        buffer[k] = 'n';
        break;
      case 'A':
        buffer[k] = 'T';
        break;
      case 'T':
        buffer[k] = 'A';
        break;
      case 'C':
        buffer[k] = 'G';
        break;
      case 'G':
        buffer[k] = 'C';
        break;
      case 'R':
        buffer[k] = 'Y';
        break;
      case 'Y':
        buffer[k] = 'R';
        break;
      case 'S':
        buffer[k] = 'S';
        break;
      case 'W':
        buffer[k] = 'W';
        break;
      case 'K':
        buffer[k] = 'M';
        break;
      case 'M':
        buffer[k] = 'K';
        break;
      case 'B':
        buffer[k] = 'V';
        break;
      case 'D':
        buffer[k] = 'H';
        break;
      case 'H':
        buffer[k] = 'D';
        break;
      case 'V':
        buffer[k] = 'B';
        break;
      case 'N':
        buffer[k] = 'N';
        break;
      default:
        buffer[k] = s[i-1]; 
        break;
    }
    k++;
  }
  buffer[k] = '\0';
  return buffer;
}

static inline void bl_convertBisulfite(char *seq, Uint len, Uint bisulfite, Uint seed) {
  if (seed){
    /* bisulfite or PARCLIP in run 1 */
    if (bisulfite <=4 && bisulfite % 2 == 1){
      //fprintf(stderr,"seed conv of reads: C --> T\n");
      strconvert(seq, len, 'C', 'T');
    }
    /* bisulfite or PARCLIP in run 2 */
    if (bisulfite <=4 && bisulfite % 2 == 0){
      //fprintf(stderr,"seed conv of reads: G --> A\n");
      strconvert(seq, len, 'G', 'A');
    }
  }
  else {  
    /* bisulfite or PARCLIP with 4SG in run 1 */
    if (bisulfite == 1){
      //fprintf(stderr,"align conv of reads: T --> Y\n");
      strconvert(seq, len, 'T', 'Y');
    }
    /* bisulfite or PARCLIP with 4SG in run 2 */
    if (bisulfite == 2){
      //fprintf(stderr,"align conv of reads: A --> R\n");
      strconvert(seq, len, 'A', 'R');
    }
    /* PARCLIP with 4SU in run 1 */
    if (bisulfite == 3 || bisulfite == 5){
      //fprintf(stderr,"align conv of reads: C --> Y\n");
      strconvert(seq, len, 'C', 'Y');
    }  
    /* PARCLIP with 4SU in run 2 */
    if (bisulfite == 4 || bisulfite == 6){
      //fprintf(stderr,"align conv of reads: G --> R\n");
      strconvert(seq, len, 'G', 'R');
    }  
  }
}

static inline void bl_reconvertBisulfite(char *seq, Uint len, Uint bisulfite) { 
  /* 
   * restoring original state is only possible in case of seed == 0
   * (collapsing the alphabet is unrecoverable)
   * => hence seed parameter omitted here
   */
  
  /* bisulfite or PARCLIP with 4SG in run 1 */
  if (bisulfite == 1){
    strconvert(seq, len, 'Y', 'T');
  }
  /* bisulfite or PARCLIP with 4SG in run 2 */
  if (bisulfite == 2){
    strconvert(seq, len, 'R', 'A');
    }
  /* PARCLIP with 4SU in run 1 */
  if (bisulfite == 3 || bisulfite == 5){
    strconvert(seq, len, 'Y', 'C');
  }  
  /* PARCLIP with 4SU in run 2 */
  if (bisulfite == 4 || bisulfite == 6){
    strconvert(seq, len, 'R', 'G');
  }
}
 
 #endif
