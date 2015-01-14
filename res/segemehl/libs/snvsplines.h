
/*
 *
 *	snvsplines.h
 *  
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 06/18/14 14:37:44 CEST  
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include "sort.h"
#include "alignment.h"
#include "stringutils.h"
#include "basic-types.h"
#include "mathematics.h"
#include "matfile.h"
#include "bitVector.h"
#include "info.h"
#include "vtprogressbar.h"
#include "fileio.h"
#include "matchfilesfields.h"
#include "matchfiles.h"
#include "debug.h"
#include "evalmatchfiles.h"
#include "biofiles.h"
#include "splicesites.h"


/*-------------------------------- getcutoff ---------------------------------
 *    
 * @brief calculate the cutoff for the calculation
 * @author Steve Hoffmann 
 *   
 */

void
getcutoff (matchfileSampleStats_t *stats, char *histofilename, 
    char *scorefilename, char *cutfilename, char* splinefilename, 
    char* estimatefilename);

