#ifndef MMCHAR_H
#define MMCHAR_H

/*
 *
 *	mmchar.h
 *  declaration for searches manber myers style
 *  on enhanced suffix arrays
 * 
 *  @author Steve Hoffmann, shoffmann@zbh.uni-hamburg.de
 *  @company Center for Bioinformatics, Hamburg 
 *  @date 12/22/06 19:11:16 CET  
 *
 *  SVN
 *  Revision of last commit: $Rev: 19 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-05-14 15:43:29 +0200 (Wed, 14 May 2008) $
 *
 *  Id: $Id: mmchar.h 19 2008-05-14 13:43:29Z steve $
 *  Url: $URL: file:///homes/bierdepot/steve/svn/segemehl/trunk/libs/sufarray/mmchar.h $
 */
 
 #include "basic-types.h"
 #include "sufarray.h"

 int mmleft(Suffixarray *, char*, Uint, int, int, int);
 int mmright(Suffixarray *, char*, Uint, int, int, int);
 PairSint mmcompare(Suffixarray *, char*, Uint, int, int);
 PairSint mmsearch(Suffixarray *, char*, Uint, Uint, Uint, Uint);

#endif
