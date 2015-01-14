
/*
 *  ncursesext.c
 *  extensions for ncurses
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 10/16/2010 03:14:16 PM CEST
 *  
 */
#include <stdlib.h>

#include <ncurses.h>
#include "ncursesext.h"
#include "basic-types.h"

SHADOWEDWINDOW*
newshadowedwin(Uint nlines, Uint ncols, Uint xpos, Uint ypos) {
  SHADOWEDWINDOW *w;

  w = calloc(1, sizeof(SHADOWEDWINDOW));
  w->main = newwin(nlines, ncols, xpos, ypos);
  w->shadow = newwin(nlines, ncols, xpos+1, ypos+2);
  return w;
}

void 
shadowedwbkgd(SHADOWEDWINDOW *w, int maincol, int shadowcol) {
  if(!w) return;
  wbkgd(w->shadow, shadowcol);
  wbkgd(w->main, maincol);
}

void
shadowedwrefresh(SHADOWEDWINDOW *w){
  if(!w) return;
  wrefresh(w->shadow);
  wrefresh(w->main);
}

void
delshadowedwin(SHADOWEDWINDOW *w) {
  if(!w) return;
  delwin(w->shadow);
  delwin(w->main);
  free(w);
}


WINDOW* 
dershadowedwin(SHADOWEDWINDOW* w, Uint nlines, Uint ncols, Uint xpos, Uint ypos) {
  WINDOW *dw;

  if(!w) return NULL;
  dw = derwin(w->main, nlines, ncols, xpos, ypos);
  return dw;
}

int
mvwaddchattr(WINDOW *win, int y, int x, int attr, chtype ch) {

  wattrset(win, attr);  
  return mvwaddch(win, y, x, ch);
}

int
mvwprintwattr(WINDOW *win, int y, int x, int attr, char *fmt, ...) {
  va_list ap;
  int ret;

  va_start(ap, fmt);
  wattrset(win, attr);  
  ret = mvwprintw(win, y, x, fmt, ap);
  va_end(ap);

  return ret;
}

