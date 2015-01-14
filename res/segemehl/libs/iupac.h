#ifndef IUPAC_H
#define IUPAC_H

/**
 * iupac.h
 * declarations for IUPAC nucleotide code
 *
 * @author Christian Otto
 * @email christian@bioinf.uni-leipzig.de
 * @company Bioinformatics, University of Leipzig
 * @date Fri Jul 23 15:03:08 CEST 2010
 */

/*
 *  SVN
 *  Revision of last commit: $Rev: 144 $
 *  Author: $Author: steve $
 *  Date: $Date: 2010-09-02 05:58:04 -0400 (Thu, 02 Sep 2010) $
 *  Id: $Id: iupac.h 144 2010-09-02 09:58:04Z steve $
 *  Url: $URL: http://www2.bioinf.uni-leipzig.de/svn5/segemehl/libs/iupac.h $
 */

#include "basic-types.h"

void initIUPAC(Uint maxqryamb, Uint maxseqamb);
BOOL isallowedIUPAC();
BOOL matchIUPAC(char qrych, char seqch);
BOOL couldMatchIUPAC(char qrych);
Uint countAmbChars(char *seq, Uint len);
Uint countNonMatchingChars(char *seq, Uint len);
double minshannonentropy(char *seq, Uint len);

#endif /* IUPAC_H */
 
