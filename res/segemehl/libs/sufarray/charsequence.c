
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "memory.h"
#include "stringutils.h"
#include "charsequence.h"

/*  
 *  SVN
 *  Revision of last commit: $Rev: 87 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-11-20 11:24:26 +0100 (Thu, 20 Nov 2008) $
 *
 *  Id: $Id: charsequence.c 87 2008-11-20 10:24:26Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/sufarray/charsequence.c $
 */

/* ----------------------------- printSequence -------------------------------
 *    
 * prints a beautiful CharSequence to a buffer 
 * 
 */
 
char * 
printSequence (void *space, CharSequence *s, Uint cols)
{
    Uint i, c, k, pos=0, l=0, width=0;
	char *buf, *buf2;
    stringset_t *entries;

	entries=initStringset(space);	
	for (i=0; i < s->length; i++)
	{	    
	  	/*buf = ALLOCMEMORY(space, NULL, char, 32);
	    buf = my_itoa(s->sequence[i], buf, 10);*/	
		addString(space, entries, &s->sequence[i], strlen(&s->sequence[i]));
		if(SETSTRLEN(entries, i) > width) 
			width = SETSTRLEN(entries, i);		
	} 
	
	/*add spacer*/
	width++;
	c = cols / (width+1);
	l = (s->descrlen+2)+(s->namelen+2)+(s->length*(width+1))+((entries->noofstrings/c)*(5+1+1));
	
	buf = ALLOCMEMORY(space, NULL, char, l);
	memset(&buf[pos++], '>', 1);
	memmove(&buf[pos], s->url, s->urllen);
	pos+=s->urllen;
	memset(&buf[pos++], '\n', 1);
	memset(&buf[pos++], '>', 1);
	memmove(&buf[pos], s->alphabetname, s->namelen);
	pos+=s->namelen;
	memset(&buf[pos++], '\n', 1);
	for (i=0; i < entries->noofstrings; i++) {
		if((i%c) == 0) {
			memset(&buf[pos++], '\n', 1);	
	    	buf2 = ALLOCMEMORY(space, NULL, char, 5);
			buf2 = my_itoa(i, buf2, 10);
			memset(&buf[pos], ' ', 5-strlen(buf2));
			pos += 5-strlen(buf2);
			memmove(&buf[pos], buf2, strlen(buf2));
			pos += strlen(buf2);
			memset(&buf[pos++], '\t', 1);
			FREEMEMORY(space, buf2);
		}
	    k = (width-SETSTRLEN(entries,i));
		memset(&buf[pos], ' ', k);
		pos += k;

	  	memmove(&buf[pos], SETSTR(entries, i), SETSTRLEN(entries, i));
		pos += SETSTRLEN(entries, i);
	}
	
	buf[pos]='\0';
		
	destructStringset(space, entries);
	return buf;
}

/* ----------------------------- printAlignment -------------------------------
 *    
 * prints a beautiful alignment to a buffer 
 * 
 */
 
char * 
printAlignment (void *space, int *align, Uint size, 
				CharSequence *a, CharSequence *b, Uint cols)
{
    Uint i, c, k, pos1=0, pos2=0, pos3=0, l=0, width=0, m=0, d=size;
	char *buf, *bufa, *bufb, *nobuf;
    stringset_t *first, *second;

	first=initStringset(space);	
	for (i=0; i < a->length; i++)
	{	    
	  	/*buf = ALLOCMEMORY(space, NULL, char, 32);
	    buf = my_itoa(a->sequence[i], buf, 10);*/	
		addString(space, first, &a->sequence[i], strlen(&a->sequence[i]));
		if(SETSTRLEN(first, i) > width) 
			width = SETSTRLEN(first, i);		 
	}

	
	second = initStringset(space);
	for (i=0; i < b->length; i++)
	{	    
	  	/*buf = ALLOCMEMORY(space, NULL, char, 32);
	    buf = my_itoa(b->sequence[i], buf, 10);*/	
		addString(space, second, &b->sequence[i], strlen(&b->sequence[i]));
		if(SETSTRLEN(second, i) > width) 
			width = SETSTRLEN(second, i);		
	} 

	
	/*add spacer*/
	width++;
	c = cols / (width+1);
	l = (a->descrlen+2)+(a->namelen+2)
	  +(a->length*(width+1)) +((first->noofstrings/c)*(5+1+1));
	m= (b->length*(width+1)) +((second->noofstrings/c)*(5+1+1));
	
	bufa  = ALLOCMEMORY(space, NULL, char, l);
	bufb = ALLOCMEMORY(space, NULL, char, m);
	buf = ALLOCMEMORY(space, NULL, char, (l+m)*2);
	
	memset(bufa, 0, l);
	memset(bufb, 0, m);
	memset(buf, 0, (l+m)*2);
	
	memset(&bufa[pos1++], '>', 1);
	memmove(&bufa[pos1], a->url, a->urllen);
	pos1+=a->urllen;
	memset(&bufa[pos1++], '\n', 1);
	memset(&bufa[pos1++], '>', 1);
	memmove(&bufa[pos1], a->alphabetname, a->namelen);
	pos1+=a->namelen;
	memset(&bufa[pos1++], '\n', 1);
	
	for (i=0; i < first->noofstrings; i++) {
		if((i%c) == 0) {
			memset(&bufa[pos1++], '\n', 1);	
	    	memset(&bufb[pos2++], '\n', 1);
			
			memmove(&buf[pos3], bufa, pos1);
			pos3+=pos1;
			memmove(&buf[pos3], bufb, pos2);
			pos3+=pos2;
			memset(&buf[pos3++], '\n', 1);

			pos1 =0;
			pos2 =0;
			
			nobuf = ALLOCMEMORY(space, NULL, char, 5);
			nobuf = my_itoa(i, nobuf, 10);
			memset(&bufa[pos1], ' ', 5-strlen(nobuf));
			pos1 += 5-strlen(nobuf);		
			memmove(&bufa[pos1], nobuf, strlen(nobuf));
			pos1 += strlen(nobuf);
			
			if (d > 1)
			  nobuf = my_itoa(align[d-1], nobuf, 10);
			else
			  nobuf = my_itoa(align[1], nobuf, 10);
			  
			memset(&bufb[pos2], ' ', 5-strlen(nobuf));
			pos2 += 5-strlen(nobuf);	
			memmove(&bufb[pos2], nobuf, strlen(nobuf));
			pos2 += strlen(nobuf);
			
			memset(&bufa[pos1++], '\t', 1);
			memset(&bufb[pos2++], '\t', 1);	
		
			FREEMEMORY(space, nobuf);	
		}

	    k = (width-SETSTRLEN(first,i));
		memset(&bufa[pos1], ' ', k);
		pos1 += k;
	  	memmove(&bufa[pos1], SETSTR(first, i), SETSTRLEN(first, i));
		pos1 += SETSTRLEN(first, i);
		
		if (d > 1 && align[d-2]-1 == i) {
			k = (width-SETSTRLEN(second, align[d-1]-1));
			memset(&bufb[pos2], ' ', k);
		    pos2 += k;
	  		memmove(&bufb[pos2], SETSTR(second, align[d-1]-1), 
				SETSTRLEN(second, align[d-1]-1));
			pos2 += SETSTRLEN(second, align[d-1]-1);
			d-=2;
		} else {
		  	k = width-1;
			memset(&bufb[pos2], ' ', k);
		    pos2 += k;
	  		memset(&bufb[pos2], '-', 1);
			pos2++; 
		}		
	}
	
	memset(&bufa[pos1++], '\n', 1);	
	memset(&bufb[pos2++], '\n', 1);
			
	memmove(&buf[pos3], bufa, pos1);
	pos3+=pos1;
	memmove(&buf[pos3], bufb, pos2);
	pos3+=pos2;
	memset(&buf[pos3++], '\n', 1);


	buf[pos3]='\0';
		
	destructStringset(space, first);
	destructStringset(space, second);
	FREEMEMORY(space, bufa);
	FREEMEMORY(space, bufb);

	return buf;
}




/* ------------------------------ dumpSequence -------------------------------
 *  dumps the sequence to the screen.  
 * 
 */

void
dumpSequence (CharSequence *s)
{
  Uint i;
  printf("sequence:\n");
  for (i=0; i < s->length; i++) { 
	printf("%c", s->sequence[i]);
	if (i!=(s->length-1)) printf("-");
  }
  printf("\n");
/*  printf("info:\n");
  for (i=0; i < s->length; i++) { 
	printf("%c", s->info[i]);
	if (i!=(s->length-1)) printf("-"); 
  }
  printf("\n");
*/
}

/*------------------------------- loadSequence -------------------------------
 *  loads a sequence from a file.   
 * 
 */
 
CharSequence*
loadSequence (void *space, char *filename)
{	  	
  //long size;
  FILE *infile;	
  char *sequence;
  /*Uint *info;*/
  CharSequence *s;
  
  infile	= fopen( filename, "r" );
  if ( infile == NULL )
  {
	fprintf ( stderr, "couldn't open file '%s'; %s\n",
		filename, strerror(errno) );
	exit (EXIT_FAILURE);
  }
  fseek(infile, 0, SEEK_END);
  //not used: size = ftell(infile);
  rewind(infile);
   
  s = initSequence(space);
  fread(s, sizeof(CharSequence), 1, infile);
  
  s->description = ALLOCMEMORY(space, NULL, char, s->descrlen+1);
  s->alphabetname = ALLOCMEMORY(space, NULL, char, s->namelen+1);
  s->url = ALLOCMEMORY(space, NULL, char, s->urllen+1);
  sequence = ALLOCMEMORY(space, NULL, char, s->length);
  /*info = ALLOCMEMORY(space, NULL, Uint, s->length);*/
 
  fread(s->description, sizeof(char), s->descrlen+1, infile);
  fread(s->alphabetname, sizeof(char), s->namelen+1, infile);
  fread(s->url, sizeof(char), s->urllen+1, infile);
  fread(sequence, sizeof(char), s->length, infile);
  /*fread(info, sizeof(Uint), s->length, infile);*/
 
  s->sequence = sequence;
  /*s->info = info;*/
  s->noofinfo = 0;
  s->info = NULL;

  if( fclose(infile) == EOF)			/* close input file */
  {
	fprintf ( stderr, "couldn't close file '%s'; %s\n",
		filename, strerror(errno) );
	exit (EXIT_FAILURE);
  }
  
  return s;
}

/*------------------------------- saveSequence -------------------------------
 *
 *  saves the sequences to a file
 *  
 */

void
saveSequence (CharSequence *s, char *filename)
{   
  FILE	*outfile;
 
  outfile = fopen (filename, "w");
  if (outfile == NULL)
  {
	fprintf ( stderr, "couldn't open file '%s'; %s\n",
		filename, strerror(errno) );
	exit (EXIT_FAILURE);
  }
  
  fwrite(s, sizeof(CharSequence), 1, outfile);
  fwrite(s->description, sizeof(char), s->descrlen+1, outfile);
  fwrite(s->alphabetname, sizeof(char), s->namelen+1, outfile);
  fwrite(s->url, sizeof(char), s->urllen+1, outfile);
  fwrite(s->sequence, sizeof(char), s->length, outfile);
  fwrite(s->info, sizeof(Uint), s->length, outfile);
  
  if( fclose(outfile) == EOF )		
  {
	fprintf ( stderr, "couldn't close file '%s'; %s\n",
		filename, strerror(errno) );
	exit (EXIT_FAILURE);
  }
	return ;
}



/*---------------------------- createSequenceHash ----------------------------
 *    
 * creates a hash table (array) for CharSequences
 * 
 */
 
CharSequence**
createSequenceHash (void *space, Uint hashsize)
{
    CharSequence** hashTable;
    hashTable = ALLOCMEMORY(space, NULL, CharSequence*, hashsize);
	memset(hashTable, 0, sizeof(CharSequence*)*hashsize);
	return hashTable;
}



/*----------------------------- findSequenceHash -----------------------------
 *    
 * find a slot in the sequence hash
 * 
 */
 
void
findSequenceHash (  )
{
	return ;
}


/*---------------------------- lookupSequenceHash ----------------------------
 *    
 * lookup if a key is already present in a hash table
 * 
 */
 
void
lookupSequenceHash (  )
{
	return ;
}






/*------------------------------ resetSequence -------------------------------
 *    
 * @brief resets (initalizes) a Sequence
 * @author Steve Hoffmann 
 *   
 */
 
void
resetSequence (CharSequence* s)
{
	s->description=NULL;
	s->descrlen=0;
	s->alphabetname=NULL;
	s->namelen=0;
	s->sequence=NULL;
	s->info=NULL;
	s->length=0;
	s->urllen =0;
	
	return ;
}

/*------------------------------- initSequence -------------------------------
 *    
 * allocates and initializes a sequence struct
 * 
 */
 
CharSequence*
initSequence (void *space)
{
	CharSequence *s;

	s=ALLOCMEMORY(space, NULL, CharSequence, 1);
	s->description=NULL;
	s->descrlen=0;
	s->alphabetname=NULL;
	s->namelen=0;
	s->sequence=NULL;
	s->info=NULL;
	s->length=0;
	s->urllen =0;
	s->url = NULL;
  	s->map=NULL;
    s->mapsize=0;
    s->clip3[0] = 0;
    s->clip3[1] = 0;
    s->clip5[0] = 0;
    s->clip5[1] = 0;
    #ifdef HASHING
    s->quantity = 1;
    #endif

    return s;
}

void destructSequence(void *space, CharSequence *sequence) {
  if (sequence->sequence != NULL)
  	FREEMEMORY(space, sequence->sequence);
	FREEMEMORY(space, sequence->description);
	if (sequence->alphabetname != NULL)
    FREEMEMORY(space, sequence->alphabetname);
    if (sequence->url != NULL)
	FREEMEMORY(space, sequence->url);
    if (sequence->info != NULL)
	FREEMEMORY(space, sequence->info);

	FREEMEMORY(space, sequence);
	return;
 }

