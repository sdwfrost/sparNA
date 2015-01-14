#ifndef SAX_H
#define SAX_H

/*
 *
 *	SAX.h
 *  Keoghs symbolic aggregate approximation
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 07/06/2011 09:38:28 PM CEST  
 *
 */


const double equiprob[4][4] = { {.0,.0,.0,.0}, {-.43,.43,.0,.0}, {-.67,.0,.67,.0}, {-.84,-.25,.25,.84}};
const char equichar[4][5] = {{'a','b','-','-','-'}, {'a','b','c','-','-'}, {'a','b','c','d','-'}, {'a','b','c','d','e'}};

char* bl_SAX (void *space, double *C, Uint n, Uint w, Uint b, Uint *P, Uint **SxP);

#endif
