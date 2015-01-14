#ifndef MANOUTFORMATS_H
#define MANOUTFORMATS_H

/*
 * outformats.h
 * definition of used symbols and output formats
 *
 * @author Christian Otto
 * @email christian@bioinf.uni-leipzig.de
 * @date Thu Oct  2 09:59:57 CEST 2008
 */

#define SEPARATOR "\t"
#define OUTLENGTH 56


/* definition of used symbols*/

#define QRY_LEN 0    /* length of query */

#define SCR 1        /* score of match */
#define EVALUE 2     /* evalue of match */
#define QRY_S 3      /* start position on query */
#define QRY_E 4      /* end position on query */

#define SEQ_S 5      /* start position on sequence (absolute) */
#define SEQ_E 6      /* end position on sequence (absolute) */
#define MAT 7        /* number of matching symbols */
#define MIS 8        /* number of mismatching symbols */
#define INS 9        /* number */
#define DEL 10       /* number of deletions */
#define EDIST 11     /* alignment edist */
#define REF_SEQ 12       /* sequence of match */
#define QRY_SEQ 13       /* the query sequence*/
#define STRAND 14    /* strand of match */

#define QRY_DESC 15  /* description of query */
#define NOOFMATCHES 16 /* no of matches */
#define PAIR_STATUS 17 /* pair status */
#define SEQ_DESC 18  /* sequence description in multi fasta file */
#define MEOP_STR 19 /* meop string*/

/*mate information*/

#define MATE_LEN 20
#define MATE_SCR 21
#define MATE_EVALUE 22
#define MATE_QRY_S 23
#define MATE_QRY_E 24

#define MATE_SEQ_S 25
#define MATE_SEQ_E 26
#define MATE_MAT 27
#define MATE_MIS 28
#define MATE_INS 29
#define MATE_DEL 30
#define MATE_EDIST 31 
#define MATE_REF_SEQ 32
#define MATE_QRY_SEQ 33
#define MATE_STRAND 34
#define MATE_SEQ_DESC 35
#define MATE_MEOP 36
#define MATE_NOOFMATCHES 37
#define MATE_DESC 38 
#define MATE_QUAL 39

#define QUAL 41 /*quality string*/
#define SAM_QRY 42 /*SAM query name*/
#define SAM_FLAG 43  /*SAM FLAGS*/
#define SAM_MAPQ 44 /*SAM reference name*/
#define SAM_CIGAR 45 /*SAM CIGAR*/
#define SAM_QRY_REF 46 /*SAM mate query*/
#define SAM_ISIZE 47 /*insert size*/

#define SAM_MATE_FLAG 48  /*SAM FLAGS*/
#define SAM_MATE_MAPQ 49 /*SAM reference name*/
#define SAM_MATE_CIGAR 50 /*SAM CIGAR*/
#define SAM_MATE_REF 51 /*SAM reference of mate*/
#define SAM_MATE_ISIZE 52
#define TAG 53
#define MATE_TAG 54
#define SAM_MATE_QRY 55

 const char EMPTY[] = " ";
 const char HEAD0[] = "#descr;score;Evalue;qstart;qend;matches;mismatches;insertions;deletions;strand;sstart;send;sequence;sequence descr\n";
 const int FORMAT0[] = {QRY_DESC, SCR, EVALUE, QRY_S, QRY_E, MAT, MIS, INS, DEL, STRAND, SEQ_S, SEQ_E, REF_SEQ, SEQ_DESC, -1};
 const char HEAD1[] = "#descr;score;qstart;qend;matches;mismatches;insertions;deletions;strand;sstart;send;sequence\n";
 const int FORMAT1[] = {QRY_DESC, SCR, QRY_S, QRY_E, MAT, MIS, INS, DEL, STRAND, SEQ_S, SEQ_E, REF_SEQ, -1};
 const char HEAD2[] = "#gff-format\n";
 const int FORMAT2[] = {SEQ_S, SEQ_E, SCR, STRAND, QRY_S, QRY_E, MAT, MIS, INS, DEL, QRY_DESC, QRY_LEN, REF_SEQ, -1};
 const char HEAD3[] = "#descr;score;Evalue;qstart;qend;matches;mismatches;insertions;deletions;strand;sstart;send;sequence descr";
 const int FORMAT3[] = {QRY_DESC, SCR, EVALUE, QRY_S, QRY_E, MAT, MIS, INS, DEL, STRAND, SEQ_S, SEQ_E, SEQ_DESC, -1};
 const char HEAD4[] = "#descr;full alignment edist;fragment score;fragment Evalue;fragment qstart;fragment qend;fragment matches;fragment mismatches;fragment insertions;fragment deletions;strand;sstart;send;sequence descr";
 const int FORMAT4[] = {QRY_DESC, EDIST, SCR, EVALUE, QRY_S, QRY_E, MAT, MIS, INS, DEL, STRAND, SEQ_S, SEQ_E, SEQ_DESC, -1};
 const char HEAD5[] = "#descr;sstart;send;strand;edist;sequence descr";
 const int FORMAT5[] = {QRY_DESC, SEQ_S, SEQ_E, STRAND, EDIST, SEQ_DESC, -1};
 const char SORTBIN5[] = "-k2,2n";
 const char SORT5[] = "-k5,5 -k2,2n";

 const char HEAD6[] = "#descr;sstart;send;strand;edist;sequence descr\n";
 const int FORMAT6[] = {QRY_DESC, MATE_SEQ_S, MATE_SEQ_E, MATE_STRAND, MATE_EDIST, MATE_SEQ_DESC, -1};

// const char HEAD6[] = "#descr;full alignment edist;fragment score;fragment Evalue;fragment qstart;fragment qend;fragment matches;fragment mismatches;fragment insertions;fragment deletions;strand;sstart;send;subject;query;sequence descr\n";
// const int FORMAT6[] = {QRY_DESC, EDIST, SCR, EVALUE, QRY_S, QRY_E, MAT, MIS, INS, DEL, STRAND, SEQ_S, SEQ_E, REF_SEQ, QRY_SEQ, SEQ_DESC, -1};
 const char HEAD7[] = "#descr;seed score;seed Evalue;seed qstart;seed qend;semi global alignment matches;semi global alignment mismatches;semi global alignment insertions;semi global alginment deletions;strand;start of semi global alignment in subject(reference) sequence;end of semi global alignment in subject sequence;sequence descr;meop string";
 const int FORMAT7[] = {QRY_DESC, SCR, EVALUE, QRY_S, QRY_E, MAT, MIS, INS, DEL, STRAND, SEQ_S, SEQ_E, SEQ_DESC, MEOP_STR, -1};
 const char HEAD8[] = "#pair status;descr;semi global alignment distance;seed score;seed Evalue;seed qstart;seed qend;semi global alignment matches;semi global alignment mismatches;semi global alignment insertions;semi global alginment deletions;strand;start of semi global alignment in subject(reference) sequence;end of semi global alignment in subject sequence;sequence descr;meop string;number of matches";
 const int FORMAT8[] = {PAIR_STATUS, QRY_DESC, EDIST, SCR, EVALUE, QRY_S, QRY_E, MAT, MIS, INS, DEL, STRAND, SEQ_S, SEQ_E, SEQ_DESC, MEOP_STR, NOOFMATCHES, -1};
 const char HEAD9[] = "#descr;semi global alignment distance;seed score;seed Evalue;seed qstart;seed qend;semi global alignment matches;semi global alignment mismatches;semi global alignment insertions;semi global alginment deletions;strand;start of semi global alignment in subject(reference) sequence;end of semi global alignment in subject sequence;sequence descr;meop stringi;query";
 const int FORMAT9[] = {QRY_DESC, EDIST, SCR, EVALUE, QRY_S, QRY_E, MAT, MIS, INS, DEL, STRAND, SEQ_S, SEQ_E, SEQ_DESC, MEOP_STR, QRY_SEQ, -1};

 const char HEAD10[] = "#descr;semi global alignment distance;seed score;seed Evalue;seed qstart;seed qend;semi global alignment matches;semi global alignment mismatches;semi global alignment insertions;semi global alginment deletions;strand;start of semi global alignment in subject(reference) sequence;end of semi global alignment in subject sequence;sequence descr;meop stringi;query";
 const int FORMAT10[] = {QRY_DESC, EDIST, SCR, EVALUE, QRY_S, QRY_E, MAT, MIS, INS, DEL, STRAND, SEQ_S, SEQ_E, SEQ_DESC, MEOP_STR, QRY_SEQ, -1};

/*query and mate*/
 const char HEAD11[] = "#pair status;descr;semi global alignment distance;seed score;seed Evalue;seed qstart;seed qend;semi global alignment matches;semi global alignment mismatches;semi global alignment insertions;semi global alginment deletions;strand;start of semi global alignment in subject(reference) sequence;end of semi global alignment in subject sequence;sequence descr;meop string;number of matches";
 const int FORMAT11[] = {PAIR_STATUS, QRY_DESC, EDIST, QRY_S, QRY_E, MAT, MIS, INS, DEL, STRAND, SEQ_S, SEQ_E, SEQ_DESC, MEOP_STR, MATE_EDIST, MATE_STRAND, MATE_SEQ_S, MATE_SEQ_E, NOOFMATCHES, -1};

/*unpaired query or mate*/
 const char HEAD12[] = "#pair status;descr;semi global alignment distance;seed score;seed qstart;seed qend;semi global alignment matches;semi global alignment mismatches;semi global alignment insertions;semi global alginment deletions;strand;start of semi global alignment in subject(reference) sequence;end of semi global alignment in subject sequence;sequence descr;meop string;number of matches;number of mate matches";
 const int FORMAT12[] = {PAIR_STATUS, QRY_DESC, EDIST, QRY_S, QRY_E, MAT, MIS, INS, DEL, STRAND, SEQ_S, SEQ_E, SEQ_DESC, MEOP_STR, NOOFMATCHES, -1};

 const char SORTBIN12[] = "-k11,11n";
 const char SORT12[] = "-k13,13 -k11,11n";

/*mate and query*/
 const char HEAD13[] = "#pair status;descr;semi global alignment distance;seed score;seed Evalue;seed qstart;seed qend;semi global alignment matches;semi global alignment mismatches;semi global alignment insertions;semi global alignment deletions;strand;start of semi global alignment in subject(reference) sequence;end of semi global alignment in subject sequence;sequence descr;meop string;number of matches;number of mate matches";
 const int FORMAT13[] = {PAIR_STATUS, MATE_DESC, MATE_EDIST, MATE_QRY_S, MATE_QRY_E, MATE_MAT, MATE_MIS, MATE_INS, MATE_DEL, MATE_STRAND, MATE_SEQ_S, MATE_SEQ_E, MATE_SEQ_DESC, MATE_MEOP, EDIST, STRAND, SEQ_S, SEQ_E, NOOFMATCHES, -1};

 const char HEAD14[] = "#descr;semi global alignment distance;seed score;seed Evalue;seed qstart;seed qend;semi global alignment matches;semi global alignment mismatches;semi global alignment insertions;semi global alignment deletions;strand;start of semi global alignment in subject(reference) sequence;end of semi global alignment in subject sequence;sequence descr;meop string;query";
 const int FORMAT14[] = {QRY_DESC, EDIST, SCR, EVALUE, QRY_S, QRY_E, MAT, MIS, INS, DEL, STRAND, SEQ_S, SEQ_E, SEQ_DESC, MEOP_STR, QRY_SEQ, MATE_EDIST, MATE_MAT, MATE_MIS, MATE_INS, MATE_DEL, MATE_STRAND, MATE_SEQ_S, MATE_SEQ_E, MATE_MEOP, MATE_QRY_SEQ, -1};

/*SAM*/

 const char HEAD15[] = "SAM";
 const int FORMAT15[] = {SAM_QRY, SAM_FLAG, SEQ_DESC, SEQ_S, SAM_MAPQ, SAM_CIGAR, SAM_MATE_REF, MATE_SEQ_S, SAM_ISIZE, QRY_SEQ, QUAL, TAG, -1};
 const char SORTBIN15[] = "-k4,4n";
 const char SORT15[] = "-k3,3 -k4,4n";

 const char HEAD16[] = "SAM";
 const int FORMAT16[] = {SAM_MATE_QRY, SAM_MATE_FLAG, MATE_SEQ_DESC, MATE_SEQ_S, SAM_MATE_MAPQ, SAM_MATE_CIGAR, SAM_QRY_REF, SEQ_S, SAM_MATE_ISIZE, MATE_QRY_SEQ, MATE_QUAL, MATE_TAG, -1};

/* definition of constant arrays */
 const char* HEAD[] = {HEAD0, HEAD1, HEAD2, HEAD3, HEAD4, HEAD5, HEAD6, HEAD7, HEAD8, HEAD9, HEAD10, HEAD11, HEAD12, HEAD13, HEAD14, HEAD15};
 const int* FORMAT[] = {FORMAT0, FORMAT1, FORMAT2, FORMAT3, FORMAT4, FORMAT5, FORMAT6, FORMAT7, FORMAT8, FORMAT9, FORMAT10, FORMAT11, FORMAT12, FORMAT13, FORMAT14, FORMAT15, FORMAT16};

 const char* SORT[] = {EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, SORT5, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, SORT12, EMPTY, EMPTY, SORT15};
 const char* SORTBIN[] = {EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, SORTBIN5, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, SORTBIN12, EMPTY, EMPTY, SORTBIN15};
 const char SORTDELIM = 9;

 

#endif
