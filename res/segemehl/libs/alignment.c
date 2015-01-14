
/*
 *  alignment.c
 *  implementation to handle alignments
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 02/03/2009 02:50:06 PM CET
 *  
 */
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <math.h>
#include "mathematics.h"
#include "basic-types.h"
#include "alignment.h"
#include "iupac.h"
#include "debug.h"
#include "charsequence.h"

const char decodeEop[] = {'R','D','I'};
const char ntdecode[]  = {'A', 'C', 'G', 'T', '-', 'N' };

char*
getNTcodekey(void *space) {
  char *code;
  code = ALLOCMEMORY(space, NULL, char, 256);
  memset(code, 5, sizeof(char)*256);

  code['A'] = 0;
  code['a'] = 0;
  code['C'] = 1;
  code['c'] = 1;
  code['G'] = 2;
  code['g'] = 2;
  code['T'] = 3;
  code['t'] = 3;
  code['-'] = 4;
  
  return code;
}

void 
initAlignment(Alignment *al, 
    char *u, Uint ulen, Uint uoff, 
    char *v, Uint vlen, Uint voff) {
  
  assert(uoff < ulen && voff < vlen);
  
  al->u = u;
  al->v = v;
  al->ulen = ulen;
  al->vlen = vlen;
  al->uoff = uoff;
  al->voff = voff;
  al->numofmeops = 0;
  al->meops = malloc(sizeof(Multieop)*(ulen+vlen));
  memset(al->meops, 0, sizeof(Multieop)*(ulen+vlen));
}

void 
wrapAlignment(Alignment *al) {
  free(al->meops);
  al->numofmeops = 0;
  al->u = NULL;
  al->v = NULL;
  al->vlen = 0;
  al->ulen = 0;
  al->uoff = 0;
  al->voff = 0;
}

void
copyAlignment(Alignment *to, Alignment *from) {
  
  to->meops = malloc(sizeof(Multieop)*(from->ulen+from->vlen));
  memmove(to->meops, from->meops, sizeof(Multieop)*(from->ulen+from->vlen));
  to->numofmeops = from->numofmeops;
  to->u = from->u;
  to->v = from->v;
  to->vlen = from->vlen;
  to->ulen = from->ulen;
  to->uoff = from->uoff;
  to->voff = from->voff;
}

void
countEops(Alignment *al, Uint *mat, Uint *mis, Uint *ins, Uint *del) {
  Uint i,j,k=0,l=0;
  
  *mat = 0;
  *mis = 0;
  *ins = 0;
  *del = 0;

  for(i=0; i < al->numofmeops; i++) {
    if(al->meops[i].eop != Replacement) {
      if (al->meops[i].eop == Deletion) {
        *del += al->meops[i].steps;
        l += al->meops[i].steps;
      } else {
        *ins += al->meops[i].steps;
        k += al->meops[i].steps;
      }
    } else {
      for(j=0; j < al->meops[i].steps; j++) {
        if(!matchIUPAC(al->u[k+al->uoff], al->v[l+al->voff]))
          *mis += 1; 
        else 
          *mat += 1;
        
        k++; l++;
      }
    }
  }
  return;
}

Uint
getEdist(Alignment *al) {
  Uint i,j,k=0,l=0;
  Uint edist=0;

  for(i=0; i < al->numofmeops; i++) {
    if(al->meops[i].eop != Replacement) {
      edist += al->meops[i].steps;
      if (al->meops[i].eop == Deletion) {
        l+= al->meops[i].steps;
      } else {
        k+= al->meops[i].steps;
      }
    } else {
      for(j=0; j < al->meops[i].steps; j++) {
        if(!matchIUPAC(al->u[k+al->uoff], al->v[l+al->voff])) edist++;
        k++; l++;
      }
    }
  }
  return edist;
}

Uint
getBisulfiteMismatches(Alignment *al, Uint bisulfite){
  Uint i,j,k=0,l=0;
  Uint mis=0;
  char *seq;

  seq = malloc(al->ulen);
  memmove(seq, al->u, al->ulen);
  bl_reconvertBisulfite(seq, al->ulen, bisulfite);

  for(i=0; i < al->numofmeops; i++) {
    if(al->meops[i].eop != Replacement) {
      if (al->meops[i].eop == Deletion) {
        l+= al->meops[i].steps;
      } else {
        k+= al->meops[i].steps;
      }
    } else {
      for(j=0; j < al->meops[i].steps; j++) {
        /* bisulfite mismatch: IUPAC match on converted query but mismatch on recoverted one */
        if(matchIUPAC(al->u[k+al->uoff], al->v[l+al->voff]) &&
           seq[k+al->uoff] != al->v[l+al->voff]) mis++;
        k++; l++;
      }
    }
  }
  free(seq);
  return mis;
}

Uint getWrongStrandBisulfiteMismatches(Alignment *al, Uint bisulfite){
  Uint i,j,k=0,l=0;
  Uint mis=0;
  char *seq;

  seq = malloc(al->ulen);
  memmove(seq, al->u, al->ulen);
  /* get other bisulfite run by altering bisulfite (e.g. 1=>2, 2=>1) */
  bl_convertBisulfite(seq, al->ulen, (bisulfite % 2) ? bisulfite + 1 : bisulfite - 1, 0);

  for(i=0; i < al->numofmeops; i++) {
    if(al->meops[i].eop != Replacement) {
      if (al->meops[i].eop == Deletion) {
        l+= al->meops[i].steps;
      } else {
        k+= al->meops[i].steps;
      }
    } else {
      for(j=0; j < al->meops[i].steps; j++) {
        /* 
         * wrong strand bisulfite mismatch: IUPAC mismatch on converted query
         * but match on coverted one in other bisulfite run
         */
        if(!matchIUPAC(al->u[k+al->uoff], al->v[l+al->voff]) &&
           matchIUPAC(seq[k+al->uoff], al->v[l+al->voff])) mis++;
        k++; l++;
      }
    }
  }
  free(seq);
  return mis;
}

void
getSoftClipScores(Alignment *al, int polyAlen, int *scores, int indel, 
    int *pAscr, int *adscr, int *adlen) {

  Uint i,j,k=0,l=0;
  int polyAscore=0, adapterScore=0, adapterLen=0;

  for(i=0; i < al->numofmeops; i++) {

    if(k+al->uoff < polyAlen) { 
      if(al->meops[i].eop != Replacement) {
        polyAscore += indel * al->meops[i].steps;
        if (al->meops[i].eop == Deletion) {
          l+= al->meops[i].steps;
        } else {
          k+= al->meops[i].steps;
        }
      } else {
        for(j=0; j < al->meops[i].steps; j++) {
          if(al->u[k+al->uoff] != 'N' && al->v[k+al->voff] != 'N' 
              && !matchIUPAC(al->u[k+al->uoff], al->v[l+al->voff])) {
            polyAscore += scores[1];
          } else {
            polyAscore += scores[0];
          }

          k++; l++;
        }
      }
    } else {
      adapterLen++;
      if(al->meops[i].eop != Replacement) {
        adapterScore += indel * al->meops[i].steps;
        if (al->meops[i].eop == Deletion) {
          l+= al->meops[i].steps;
        } else {
          k+= al->meops[i].steps;
        }
      } else {
        for(j=0; j < al->meops[i].steps; j++) {
          if(!matchIUPAC(al->u[k+al->uoff], al->v[l+al->voff])) {
            adapterScore += scores[1];
          } else {
            adapterScore += scores[0];
          }

          k++; l++;
        }
      }
    }

  }

  *pAscr = polyAscore;
  *adscr = adapterScore;
  *adlen = adapterLen;

  return;
}

int
getSWScore(Alignment *al, int *scores, int indel) {
  Uint i,j,k=0,l=0;
  int score=0;

  for(i=0; i < al->numofmeops; i++) {

    if(al->meops[i].eop != Replacement) {
      score += indel * al->meops[i].steps;
      if (al->meops[i].eop == Deletion) {
        l+= al->meops[i].steps;
      } else {
        k+= al->meops[i].steps;
      }
    } else {
      for(j=0; j < al->meops[i].steps; j++) {
        if(!matchIUPAC(al->u[k+al->uoff], al->v[l+al->voff])) {
          score += scores[1];
        } else {
          score += scores[0];
        }

        k++; l++;
      }
    }
  }
  return score;
}

int
getAlignScore(Alignment *al, int *scores, int indel) {
  Uint i,j,k=0,l=0;
  int score=0;

  for(i=0; i < al->numofmeops; i++) {

    if(al->meops[i].eop != Replacement) {
      score += indel * al->meops[i].steps;
      if (al->meops[i].eop == Deletion) {
        l+= al->meops[i].steps;
      } else {
        k+= al->meops[i].steps;
      }
    } else {
      for(j=0; j < al->meops[i].steps; j++) {
        if(!matchIUPAC(al->u[k+al->uoff], al->v[l+al->voff])) {
          score += scores[1];
        } else {
          score += scores[0];
        }

        k++; l++;
      }
    }
  }
  return score;
}



int
getSubstringEdist(Alignment *al, Uint u, Uint v) {

  Uint i,j,k=0,l=0;
  Uint edist=0, scr=0;
  for(i=0; i < al->numofmeops; i++) {
    if(al->meops[i].eop != Replacement) {
      if(k >= u && k < v)
        edist += al->meops[i].steps;
      if (al->meops[i].eop == Deletion) {
        l+= al->meops[i].steps;
      } else {
        k+= al->meops[i].steps;
      }
    } else {
      for(j=0; j < al->meops[i].steps; j++) {
        if(k >=u && k < v && !matchIUPAC(al->u[k+al->uoff], al->v[l+al->voff]))
	  edist++;
        else if (k >=u && k<v) scr += 1;
        k++; l++;
      }
    }
  }
  return scr;

}

//dumps visual representation of alignments and should be shorter!
void 
showmultieoplist(FILE *dev, Alignment *al) {

  Uint i=0;
  fprintf(dev, "[");
  if(al->numofmeops) {
    for(i=0; i < al->numofmeops; i++) {
      fprintf(dev, "%c %d, ", decodeEop[al->meops[i].eop], al->meops[i].steps);
    }
  fprintf(dev, "%c %d",decodeEop[al->meops[i].eop], al->meops[i].steps);    
  }
  fprintf(dev, "]\n");
}

//dumps visual representation of alignments and should be shorter!
char *
multieopstring(Alignment *al, Uint leftclip, Uint rightclip, unsigned char rev) {
  Uint i, j, k, q=0, p=0, cur=0, strsize, steps, msteps, ssteps;
  char *meopstr;
  char eopc=0;

  meopstr = (char*) malloc(sizeof(char)*(3*(al->vlen+al->ulen+leftclip+rightclip)+1));

  if(leftclip || (rightclip && rev)) {
    steps = (rev) ? rightclip : leftclip;
    //strsize = floor(log(steps)/log(10))+3;
    strsize = snprintf(NULL, 0, "%d", steps)+2;
    meopstr[cur] = 'C';
    sprintf(&meopstr[cur+1], "%d;", steps);
    cur+=strsize;
  }

  for(k=0; k < al->numofmeops; k++) {
    i = (rev) ? al->numofmeops - k -1 : k;
    //if Replacement occured
    steps=0;
    if (al->meops[i].eop == Replacement) {
      //iter over all steps
      msteps=0;
      ssteps=0;
      for (j=0; j < al->meops[i].steps; j++) {
        if (!matchIUPAC(al->u[j+p+al->uoff], al->v[j+q+al->voff])) {
          if (j==0 || eopc == 'S') {
            ssteps++;
          } else {
            //strsize = floor(log(msteps)/log(10))+3;
            strsize = snprintf(NULL, 0, "%d", msteps)+2;
            meopstr[cur] = eopc;
            sprintf(&meopstr[cur+1], "%d;", msteps);
            cur+=strsize;
            msteps=0;
            ssteps=1;
          }
          eopc = 'S';
        } else {
          if (j==0 || eopc == 'M') {
            msteps++;
          } else {
            //strsize = floor(log(ssteps)/log(10))+3;
            strsize = snprintf(NULL, 0, "%d", ssteps)+2;
            meopstr[cur] = eopc;
            sprintf(&meopstr[cur+1], "%d;", ssteps);
            cur+=strsize;
            msteps=1;
            ssteps=0;
          }
          eopc = 'M';
        }
      }
      steps = msteps + ssteps;
      assert(msteps == 0 || ssteps == 0);
      //set string ptrs
      p+=j;
      q+=j;
    } 
    //if deletion occured
    if (al->meops[i].eop == Deletion) {
      eopc = 'D';
      //set ptrs
      steps = al->meops[i].steps;
      q+=steps;
    }
    //if insertions occured
    if(al->meops[i].eop == Insertion)  {
      eopc = 'I';  
      steps = al->meops[i].steps;
      p+=steps;
    }

    //strsize = floor(log(steps)/log(10))+3;
    strsize = snprintf(NULL, 0, "%d", steps)+2;
    meopstr[cur] = eopc;
    sprintf(&meopstr[cur+1], "%d;", steps);
    cur+=strsize;
  }

  if(rightclip || (leftclip && rev)) {
    steps = (rev) ? leftclip : rightclip;
    //strsize = floor(log(steps)/log(10))+3;
    strsize = snprintf(NULL, 0, "%d", steps)+2;
    meopstr[cur] = 'C';
    sprintf(&meopstr[cur+1], "%d;", steps);
    cur+=strsize;
  }
  return meopstr;
}

char*
mdstring(Alignment *al, unsigned char rev) {
  Uint i, j, k, q=0, p=0, cur=0, strsize, steps, msteps=0, ssteps=0;
  char *mdstr;
  char eopc=0;

  mdstr = (char*) malloc(sizeof(char)*(3*(al->vlen+al->ulen)+1));
  memset(mdstr, 0, sizeof(char)*(3*(al->vlen+al->ulen)+1));

    msteps=0;
    ssteps=0;
  for(k=0; k < al->numofmeops; k++) {
    i = (rev) ? al->numofmeops - k - 1 : k;
    //if Replacement occured
    steps=0;
    
    if (al->meops[i].eop == Replacement) {
      //iter over all steps

      for (j=0; j < al->meops[i].steps; j++) {  
        if (!matchIUPAC(al->u[j+p+al->uoff], al->v[j+q+al->voff])) {
          if(msteps) {
            //strsize = floor(log(msteps)/log(10))+1;
            strsize = snprintf(NULL, 0, "%d", msteps);
            sprintf(&mdstr[cur], "%d", msteps);
            cur += strsize;
            msteps = 0;
          }
          if(eopc != 'M') {
            sprintf(&mdstr[cur], "0");
            cur += 1;
          }
          sprintf(&mdstr[cur], "%c", al->v[j+q+al->voff]);
          cur += 1;
          ssteps++;
          eopc = 'S';
        } else {
          if (msteps) {
            msteps++;
          } else {
            msteps=1;
            ssteps=0;
          }
          eopc = 'M';
        }
      }
      steps = msteps + ssteps;
      assert(msteps == 0 || ssteps == 0);
      //set string ptrs
      p+=j;
      q+=j;
    } 

    //if deletion occured
    if (al->meops[i].eop == Deletion) {
      if (msteps) {
        //strsize = floor(log(msteps)/log(10))+1;
        strsize = snprintf(NULL, 0, "%d", msteps);
        sprintf(&mdstr[cur], "%d", msteps);
        cur += strsize;
        msteps = 0;
      } else { 
        sprintf(&mdstr[cur], "0");
        cur += 1;
      }
      
      eopc = 'D';
      //set ptrs
      steps = al->meops[i].steps;
      sprintf(&mdstr[cur], "^");
      cur+=1;
      for(j=0; j < steps; j++) {
        sprintf(&mdstr[cur], "%c", al->v[j+q+al->voff]);
        cur+=1;
      }
      q+=steps;
    }

    //if insertions occured
    if(al->meops[i].eop == Insertion)  {
      //eopc = 'I'; 
      steps = al->meops[i].steps;
      p+=steps;
    }
  }

  if(eopc != 'M') {
    sprintf(&mdstr[cur], "0");
    cur += 1;
  }
  
  if (msteps) {
    //strsize = floor(log(msteps)/log(10))+1;
    strsize = snprintf(NULL, 0, "%d", msteps);
    sprintf(&mdstr[cur], "%d", msteps);
    cur += strsize;
    msteps = 0;
  } 

  return mdstr;
}


/*-------------------------- bl_cigarGetAlignString --------------------------
 *    
 * @brief decode a cigar string
 * @author Steve Hoffmann 
 *   
 */
 
char*
bl_cigarGetAlignString(char *cigar) {
  Uint i, len, allen=0, cur=0;
  char *buffer, *string = NULL;

  len = strlen(cigar);
  buffer = calloc(len, sizeof(char));

  for(i=0; i < len; i++) {
    switch (cigar[i]) {
      case 'S':
        string = realloc(string, allen+atoi(buffer)+1);
        memset(&string[allen], 'S', atoi(buffer));
        allen += atoi(buffer);
        string[allen] = 0;
        memset (buffer, 0, len);
        cur = 0;
        break;
      case 'M':
      case 'N':
        string = realloc(string, allen+atoi(buffer)+1);
        memset(&string[allen], 'M', atoi(buffer));
        allen += atoi(buffer);
        string[allen] = 0;
        memset (buffer, 0, len);
        cur = 0;
        break;
      case 'D':         
        string = realloc(string, allen+atoi(buffer)+1);
        memset(&string[allen], 'D', atoi(buffer));
        allen += atoi(buffer);
        string[allen] = 0;
        memset (buffer, 0, len);
        cur = 0;
        break;
      case 'I':
        string = realloc(string, allen+atoi(buffer)+1);
        memset(&string[allen], 'I', atoi(buffer));
        allen += atoi(buffer);
        string[allen] = 0;
        memset (buffer, 0, len);
        cur =0;
        break;
      default :
        buffer[cur++] = cigar[i];
    }
  }

  free(buffer);
  return string;
}


/*--------------------------- bl_cigarGetAlignLen ----------------------------
 *    
 * @brief get alignment length from a cigar string
 * @author Steve Hoffmann 
 *   
 */
 
Uint
bl_cigarGetAlignLen(char *cigar) {
  Uint i, len, allen=0, cur=0;
  char *buffer;

  len = strlen(cigar);
  buffer = calloc(len, sizeof(char));

  for(i=0; i < len; i++) {
    switch (cigar[i]) {
      case 'M':
      case 'X':
      case '=':
        allen += atoi(buffer);
        memset (buffer, 0, len);
        cur = 0;
        break;
      case 'I': 
        memset (buffer, 0, len);
        cur = 0;
        break;
      case 'D':
        allen += atoi(buffer);
        memset (buffer, 0, len);
        cur =0;
        break;
      case 'N':
        memset (buffer, 0, len);
        cur = 0;
        break;
      case 'S':
        memset (buffer, 0, len);
        cur = 0;
        break;
      case 'H':
        memset (buffer, 0, len);
        cur = 0;
        break;
      case 'P':
        memset (buffer, 0, len);
        cur = 0;
        break;
      default :
        buffer[cur++] = cigar[i];
    }
  }

  free(buffer);
  return allen;
}

char*
bl_mdGetDiffString(char *MD) {
  Uint i, k=0, MDlen, allen=0, nops=0, buffersize=100;
  char chr;
  unsigned char del = 0;
  char *buffer = NULL;
  char *alignstr = NULL;

  buffer = ALLOCMEMORY(space, NULL, char, buffersize);
  memset(buffer, 0, sizeof(char)*buffersize);
  MDlen = strlen(MD); 

  for(i=0; i < MDlen; i++) {
    
    if(!isalpha((int)MD[i]) && MD[i] != '^') {
      if(k >= buffersize-2) {
        buffersize += 100;
        buffer = ALLOCMEMORY(space, buffer, char, buffersize);
      }
      buffer[k] = MD[i];
      buffer[k+1] = 0;
      k++;
      del = 0;
    } else {
      
      nops = atoi(buffer);
      if(nops) {
        alignstr = ALLOCMEMORY(space, alignstr, char, allen+nops+1);
        memset(&alignstr[allen], 'M', nops);
        alignstr[allen+nops] = 0;
        allen+=nops;
       }
       memset(buffer, 0, sizeof(char)*buffersize);
       k=0;

      if(MD[i] == '^') {
        del = 1;
      } else {
        chr = (del) ? 'D' : MD[i];
        alignstr = ALLOCMEMORY(space, alignstr, char, allen+2);
        alignstr[allen] = chr;
        alignstr[allen+1] = 0;
        allen += 1;
      }
    }
  }

  nops = atoi(buffer);
  if(nops) {
    alignstr = ALLOCMEMORY(space, alignstr, char, allen+nops+1);
    memset(&alignstr[allen], 'M', nops);
    alignstr[allen+nops] = 0;
    allen+=nops;
  }

  FREEMEMORY(space, buffer);
  return alignstr;
}

char*
cigarstring(Alignment *al, Uint leftclip, Uint rightclip, char clipch, unsigned char rev) {
  Uint i, j, k, q=0, p=0, cur=0, strsize, steps, msteps, ssteps;
  char *meopstr;
  char eopc=0;
  meopstr = (char*) malloc(sizeof(char)*(3*(al->vlen+al->ulen+leftclip+rightclip)+1));

  if(leftclip || (rightclip && rev)) {
    steps = (rev) ? rightclip : leftclip;
//    strsize = floor(log(steps)/log(10))+2;
    strsize = snprintf(NULL, 0, "%d", steps)+1;
    sprintf(&meopstr[cur], "%d%c", steps, clipch);
    cur+=strsize;
  }

  for(k=0; k < al->numofmeops; k++) {
    i = (rev) ? al->numofmeops - k - 1 : k;
    //if Replacement occured
    steps=0;
    if (al->meops[i].eop == Replacement) {
      //iter over all steps
      msteps=0;
      ssteps=0;
      for (j=0; j < al->meops[i].steps; j++) {
        if (j==0 || eopc == 'M') {
          ssteps++;
        } else {
          //strsize = floor(log(msteps)/log(10))+2;
          strsize = snprintf(NULL, 0, "%d", msteps)+1;
          sprintf(&meopstr[cur], "%d%c", msteps, eopc);
          cur+=strsize;
          msteps=0;
          ssteps=1;
        }
        eopc = 'M';
      }
      steps = msteps + ssteps;
      assert(msteps == 0 || ssteps == 0);
      //set string ptrs
      p+=j;
      q+=j;
    } 
    //if deletion occured
    if (al->meops[i].eop == Deletion) {
      eopc = 'D';
      //set ptrs
      steps = al->meops[i].steps;
      q+=steps;
    }
    //if insertions occured
    if(al->meops[i].eop == Insertion)  {
      eopc = 'I';  
      steps = al->meops[i].steps;
      p+=steps;
    }

    //strsize = floor(log(steps)/log(10))+2;
    strsize = snprintf(NULL, 0, "%d", steps)+1;
    sprintf(&meopstr[cur], "%d%c", steps, eopc);
    cur+=strsize;
  }

  if(rightclip || (leftclip && rev)) {
    steps = (rev) ? leftclip : rightclip;
    //strsize = floor(log(steps)/log(10))+2;
    strsize = snprintf(NULL, 0, "%d", steps)+1;
    sprintf(&meopstr[cur], "%d%c", steps, clipch);
    cur+=strsize;
  }
  return meopstr;
}

//shows muliteoplist of all Alignments in *al
void 
showDynmultieoplist(Alignment* al, int size) {

  int i;
  for (i=0; i < size; i++) {
    showmultieoplist(stdout, &al[i]);
  }
}

//dumps visual representation of alignments and should be shorter!
void 
showAlign(Alignment* al, FILE *dev) {
  int i, j , k, nlines, len;

  Uint p = 0, q = 0, r = 0;
  //output strings
  char* utemp = (char*) malloc(sizeof(char)*(al->vlen+al->ulen));
  char* vtemp = (char*) malloc(sizeof(char)*(al->vlen+al->ulen));
  char* comp = (char*) malloc(sizeof(char)*(al->vlen+al->ulen));

  //iter over all multieops
  for (i=0; i < al->numofmeops; i++) {
    //if Replacement occured
    if (al->meops[i].eop == Replacement) {
      //iter over all steps
      for (j=0; j < al->meops[i].steps; j++) {
        utemp[j+r] = al->u[j+p+al->uoff];
        vtemp[j+r] = al->v[j+q+al->voff];
        //real Replacement?
        if (!matchIUPAC(utemp[j+r], vtemp[j+r]))
          comp[j+r] = ' ';
        else
          comp[j+r] = '|';
      }
      //set string ptrs
      p+=j;
      q+=j;
      r+=j;
    }
    //if deletion occured
    if (al->meops[i].eop == Deletion) {
      //iter over all steps
      for (j=0; j < al->meops[i].steps; j++) {
        utemp[j+r] = '-';
        vtemp[j+r] = al->v[j+q+al->voff];
        comp[j+r] = ' ';
      }
      //set ptrs
      r+=j;
      q+=j;
    }
    //if insertions occured
    if (al->meops[i].eop == Insertion) {
      for (j=0; j < al->meops[i].steps; j++) {
        utemp[j+r] = al->u[j+p+al->uoff];
        vtemp[j+r] = '-';
        comp[j+r] = ' ';
      }
      r+=j;
      p+=j;
    }
    if(i == al->numofmeops-1) {
      //terminate strings
      utemp[r]='\0';
      vtemp[r]='\0';
      comp[r] ='\0';
      
      nlines = r/60;
      nlines += (r % 60) ? 1 : 0;
      //dump strings
      for(k=0; k < nlines; k++) {
        len = (k*60 > r) ? r % 60 : 60;
        fprintf(dev, "%.*s\n", len, &utemp[k*60]);
        fprintf(dev, "%.*s\n", len, &comp[k*60]);
        fprintf(dev, "%.*s\n", len, &vtemp[k*60]);
      }
      fprintf(dev, "\n");
      memset(utemp, 0, sizeof(char)*(al->ulen+al->vlen));
      memset(comp, 0, sizeof(char)*(al->ulen+al->vlen));
      memset(vtemp, 0, sizeof(char)*(al->ulen+al->vlen));
    }
  }

  free(utemp);
  free(comp);
  free(vtemp);
}

//dumps visual representation of alignments and should be shorter!
void 
showAlignLF(Alignment* al, FILE *dev, char lf) {
  int i, j , k, nlines, len;

  Uint p = 0, q = 0, r = 0;
  //output strings
  char* utemp = (char*) malloc(sizeof(char)*(al->vlen+al->ulen));
  char* vtemp = (char*) malloc(sizeof(char)*(al->vlen+al->ulen));
  char* comp = (char*) malloc(sizeof(char)*(al->vlen+al->ulen));

  //iter over all multieops
  for (i=0; i < al->numofmeops; i++) {
    //if Replacement occured
    if (al->meops[i].eop == Replacement) {
      //iter over all steps
      for (j=0; j < al->meops[i].steps; j++) {
        utemp[j+r] = al->u[j+p+al->uoff];
        vtemp[j+r] = al->v[j+q+al->voff];
        //real Replacement?
        if (!matchIUPAC(utemp[j+r], vtemp[j+r]))
          comp[j+r] = ' ';
        else
	     comp[j+r] = '|';
      }
      //set string ptrs
      p+=j;
      q+=j;
      r+=j;
    }
    //if deletion occured
    if (al->meops[i].eop == Deletion) {
      //iter over all steps
      for (j=0; j < al->meops[i].steps; j++) {
        utemp[j+r] = '-';
        vtemp[j+r] = al->v[j+q+al->voff];
        comp[j+r] = ' ';
      }
      //set ptrs
      r+=j;
      q+=j;
    }
    //if insertions occured
    if (al->meops[i].eop == Insertion) {
      for (j=0; j < al->meops[i].steps; j++) {
        utemp[j+r] = al->u[j+p+al->uoff];
        vtemp[j+r] = '-';
        comp[j+r] = ' ';
      }
      r+=j;
      p+=j;
    }
    if(i == al->numofmeops-1) {
      //terminate strings
      utemp[r]='\0';
      vtemp[r]='\0';
      comp[r] ='\0';
      
      nlines = r/60;
      nlines += (r % 60) ? 1 : 0;
      //dump strings
      for(k=0; k < nlines; k++) {
        len = (k*60 > r) ? r % 60 : 60;
        fprintf(dev, "%.*s%c", len, &utemp[k*60], lf);
        fprintf(dev, "%.*s%c", len, &comp[k*60], lf);
        if(k < nlines-1)
        fprintf(dev, "%.*s%c", len, &vtemp[k*60], lf);
        else
        fprintf(dev, "%.*s", len, &vtemp[k*60]);
      }
    
      memset(utemp, 0, sizeof(char)*(al->ulen+al->vlen));
      memset(comp, 0, sizeof(char)*(al->ulen+al->vlen));
      memset(vtemp, 0, sizeof(char)*(al->ulen+al->vlen));
    }
  }

  free(utemp);
  free(comp);
  free(vtemp);
}

void
insertEop(Alignment *al, Eoptype eop) {

  //if previous multieops have been set up
  if (al->numofmeops > 0)  { 
    //inc steps if curr eop matches prev eops
    if (al->meops[al->numofmeops-1].eop == eop) {
      al->meops[al->numofmeops-1].steps++;
      //set up new multieop otherwise
    } else {
      al->numofmeops++;
      al->meops[al->numofmeops-1].eop =  eop;
      al->meops[al->numofmeops-1].steps = 1;
    }
    //set up first multieop
  } else {
    al->numofmeops = 1;
    al->meops[0].eop = eop;
    al->meops[0].steps = 1;
  }
}


void
revMeops(Alignment *al) {
  Uint start = 0;
  Uint end = al->numofmeops-1;
  Multieop *meops = al->meops; 

  if (al->numofmeops == 0) return;

  while (start<end) {
    meops[start].eop ^= meops[end].eop;
    meops[start].steps ^= meops[end].steps;
    meops[end].eop ^= meops[start].eop;
    meops[end].steps ^= meops[start].steps;
    meops[start].eop ^= meops[end].eop;
    meops[start].steps ^= meops[end].steps;

    start++;
    end--;
  } 
}


Uint
getValignlen(Alignment *al) {
  Uint i, vallen=0, steps;
  Eoptype eop;

  for(i=0; i < al->numofmeops; i++) {
    eop = al->meops[i].eop;
    steps = al->meops[i].steps;
    switch(eop) {
      case Replacement:
        vallen += steps;
        break;
      case Deletion:
        vallen += steps;
        break;
      case Insertion:
        break;
    }
  }

  return vallen;
}

Uint
getUalignlen(Alignment *al) {
  Uint i, uallen=0, steps;
  Eoptype eop;

  for(i=0; i < al->numofmeops; i++) {
    eop = al->meops[i].eop;
    steps = al->meops[i].steps;
    switch(eop) {
      case Replacement:
        uallen += steps;
        break;
      case Deletion:
        break;
      case Insertion:
        uallen += steps;
        break;
    }
  }

  return uallen;
}



