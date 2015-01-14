
/*
 *
 *	ncursesext.h
 *  extensions to ncurses
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 10/16/2010 03:13:29 PM CEST  
 *
 */

#include <ncurses.h>
#include "basic-types.h"

typedef struct {
  WINDOW *main;
  WINDOW *shadow;
} SHADOWEDWINDOW;


void shadowedwbkgd(SHADOWEDWINDOW *W, int maincol, int shadowcol);
void shadowedwrefresh(SHADOWEDWINDOW *w);
SHADOWEDWINDOW* newshadowedwin(Uint nlines, Uint ncols, Uint xpos, Uint ypos);
WINDOW* dershadowedwin(SHADOWEDWINDOW*, Uint nlines, Uint ncols, Uint xpos, Uint ypos);
void delshadowedwin(SHADOWEDWINDOW *w);
int mvwaddchattr(WINDOW *win, int y, int x, int attr, chtype ch);
int mvwprintwattr(WINDOW *win, int y, int x, int attr, char *fmt, ...);
