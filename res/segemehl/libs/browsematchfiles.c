
/*
 *  browsematchfiles.c
 *  a small browser
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 06.10.2010 01:39:45 CEST
 *  
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "alignment.h"
#include "stringutils.h"
#include "basic-types.h"
#include "mathematics.h"
#include "matfile.h"
#include "bitVector.h"
#include "matchfiles.h"
#include "browsematchfiles.h"
#include "plotmatchfiles.h"
#include "info.h"
#include <form.h>
#include <menu.h>
#include <curses.h>
#include "ncursesext.h"

#define MAXCHROMSIZE 1500000000


/*--------------------------- bl_matchfileGetInfo ----------------------------
 *    
 * @brief display window with stats and cs info
 * @author Steve Hoffmann 
 *   
 */
 

SHADOWEDWINDOW*
bl_matchfileGetInfo(matchfileView_t *view) {
  int j,k;
  Uint *ntcounts, curptr=0 ;  
  double *ntprobs, *ntvarpos, *ntreadpos, *ntedist, *ntredund;
  char NT[]={'A', 'C', 'G', 'T', '-'};
  WINDOW *win;
  SHADOWEDWINDOW *shadowed;
  matchfileCross_t *cs;

  shadowed = newshadowedwin(36, 70, 9, 12);
  win = dershadowedwin(shadowed, 34, 67, 1, 1);
  shadowedwbkgd(shadowed, COLOR_PAIR(2), COLOR_PAIR(8));
  box(win,0,0);
  wbkgd(win, COLOR_PAIR(2));
  mvwprintw(win, 22, 2, "(Press 'q' to exit)");
  
  if(view->curframe && view->curframe->cs) { 
    cs = view->curframe->cs;

    if(view->imap[view->offset-1]) { 
      curptr = view->imap[view->offset-1];
    } else { 
      k = view->offset-1;
      while(view->imap[k] == 0 && k > 0) k--;
      curptr = view->imap[k];
    }

   
    bl_matchfileTest(NULL, 0, 0, view->curframe->start+curptr, 
        &view->curframe->cs[curptr], view->curframe->ref[curptr], view->idx, 1, NULL);

    wattron(win, A_BOLD);
    mvwprintw(win, 1, 2, "cross section %d:", 
        view->curframe->start+curptr);
    mvwprintw(win, 1, 40, "r: %c <> %c :c", 
        cs[curptr].ref, cs[curptr].cons);
    wattroff(win, A_BOLD);
    
    ntredund = bl_matchfileGetNTRedundancy(view->curframe, curptr);
    ntcounts = bl_matchfileGetNTCounts(&view->curframe->cs[curptr]);
    ntprobs = bl_matchfileGetNTError(view->curframe, curptr);
    ntedist = bl_matchfileGetNTEdist(view->curframe, curptr);
    ntreadpos = bl_matchfileGetNTReadPos(view->curframe, curptr);
    ntvarpos = bl_matchfileGetNTReadPosVar(&view->curframe->cs[curptr]);
 
    mvwprintw(win, 3, 2, "coverage:  %d ",
        cs[curptr].len);
    mvwprintw(win, 3, 21, "5'-ends: %d", 
        cs[curptr].starts);
    mvwprintw(win, 3, 38, "3'-ends: %d", 
        cs[curptr].ends);

    wattron(win, A_UNDERLINE);
    mvwprintw(win, 4, 20, " ");
    mvwprintw(win, 4, 26, "    ");
    mvwprintw(win, 4, 35, "     ");
    mvwprintw(win, 4, 44, "    ");
    wattroff(win,A_UNDERLINE);
    mvwprintw(win, 5, 2, 
        "      #       log(E)    edst     multi    rpos    s(rpos) ");
    mvwprintw(win, 6, 2, 
        "     ---------------------------------------------------- ");

    for(j=0,k=0; j < 5; j++) {

      if(ntcounts[(int)NT[j]]) {
        if(NT[j] == cs[curptr].ref || 
           NT[j] == cs[curptr].cons)
          wattron(win, A_BOLD);
        mvwprintw(win, 7+k, 2,  "%c:", NT[j]);
        if(NT[j] == cs[curptr].ref ||
           NT[j] == cs[curptr].ref)
          wattroff(win, A_BOLD);

        mvwprintw(win, 7+k, 8,  "%d", ntcounts[(int)NT[j]]);
        mvwprintw(win, 7+k, 16, "%.3f", ntprobs[(int)NT[j]]);
        mvwprintw(win, 7+k, 26, "%.3f", ntedist[(int)NT[j]]);
        mvwprintw(win, 7+k, 35, "%.1f", ntredund[(int)NT[j]]);
        mvwprintw(win, 7+k, 44, "%.1f", ntreadpos[(int)NT[j]]);
        mvwprintw(win, 7+k, 52, "%.3f", ntvarpos[(int)NT[j]]);
        
        mvwprintw(win,  7+k+1, 26, "%.1f", 
            log(((double)view->stats->RR[(int)MIN(ntedist[(int)NT[j]],6)] + 1.0)/
            ((double)view->stats->RR_N + 1.0)));
        
        mvwprintw(win, 7+k+1, 35, "%.1f", 
            log(((double)view->stats->MM[(int)MIN(ntredund[(int)NT[j]],50)] + 1.0)/
            ((double)view->stats->MM_N + 1.0)));
           
        mvwprintw(win, 7+k+1, 52, "%.1f",
            log(univarnormcdf(ntvarpos[(int)NT[j]], view->stats->V_mu, view->stats->V_sd)));
        
        k+=2;
      }
    }
 
    wattron(win, A_BOLD);
    mvwprintw(win, 7+k+1, 2, "frame statistics");
    wattroff(win, A_BOLD);

    wattron(win, A_UNDERLINE);
    mvwprintw(win, 7+k+2, 2, "        ");
    wattroff(win, A_UNDERLINE);
    mvwprintw(win, 7+k+3, 2, "coverage: %.2f", 
        view->curframestats->mean_cov);

    for(j=0; j < 5; j++) {

      mvwprintw(win, 7+k+5, 2,  "%c:", NT[j]);
      mvwprintw(win, 7+k+5, 8,  "%d",   
          view->curframestats->ntcnt[(int)NT[j]]);
      mvwprintw(win, 7+k+5, 16, "%.3f", 
          view->curframestats->mean_err[(int)NT[j]]);
      mvwprintw(win, 7+k+5, 26, "%.3f",  
          view->curframestats->mean_dis[(int)NT[j]]);
      mvwprintw(win, 7+k+5, 35, "%.1f",  
          view->curframestats->mean_mul[(int)NT[j]]);
      mvwprintw(win, 7+k+5, 44, "%.1f",  
          view->curframestats->mean_pos[(int)NT[j]]);
      mvwprintw(win, 7+k+5, 52, "%.1f",  
          view->curframestats->mean_sde[(int)NT[j]]);
      k++;
    }


    mvwprintw(win, 7+k+6, 2, "P(V)=%.2f, P(NV)=%.2f, P(Entropy=%.3f)=%.2f",  
        log((double)view->stats->X/view->stats->N),
        log((double)view->stats->P/view->stats->N),
        cs[curptr].entropy, cs[curptr].pentropy);

    mvwprintw(win, 7+k+7, 2, "P(V)*P(D|V)=%.2f,%.2f P(NV)*P(D|NV)=%.2f,%.2f S=%.2f",
        cs[curptr].p_refx, cs[curptr].s_refx, cs[curptr].p_ref, cs[curptr].s_ref, logadd(cs[curptr].p_ref, 
          cs[curptr].p_refx));

    mvwprintw(win, 7+k+8, 2, "P(V)*P(D|V)=%.2f,%.2f P(NV)*P(D|NV)=%.2f,%.2f S=%.2f",
        cs[curptr].p_consx, cs[curptr].s_consx, cs[curptr].p_cons, cs[curptr].s_cons, logadd(cs[curptr].p_cons, 
          cs[curptr].p_consx));
    
//    mvwprintw(win, 7+k+9, 2, "Phom=%.2f, ee_cons=%.2f, ee_consx=%.2f, ee_ref=%.2f, ee_refx=%.2f", cs[curptr].p_hom, 
//        cs[curptr].ee_cons, cs[curptr].ee_consx, cs[curptr].ee_ref, cs[curptr].ee_refx);



    FREEMEMORY(space, ntcounts);
    FREEMEMORY(space, ntprobs);
    FREEMEMORY(space, ntedist);
    FREEMEMORY(space, ntredund);
    FREEMEMORY(space, ntreadpos);
    FREEMEMORY(space, ntvarpos);

  } else {
 
    mvwprintw(win, 1, 2, "no cross section at position %d available", 
        view->curframe->start+curptr);
  }

  return shadowed;
}


/*------------------------ bl_matchfileSelectChrMenu -------------------------
 *    
 * @brief display menu for chromosome selection
 * @author Steve Hoffmann 
 *   
 */
 

char*
bl_matchfileSelectChrMenu(fasta_t *set) {
  ITEM   **it;
  MENU   *me;
  WINDOW *shadoww,*mainwin, *win;
  char **mi, *newchrname = NULL;
  int ch, i, sel, nchr;
  
  nchr = set->noofseqs;
  sel= nchr+1;

  it = (ITEM **)calloc(nchr+2, sizeof(ITEM *));
  mi = ALLOCMEMORY(space, NULL, char*, nchr);

  for(i=0;i < nchr; i++) {
    mi[i] = ALLOCMEMORY(space, NULL, char, log10(nchr)+3);
    snprintf(mi[i], log10(nchr)+2, "%d", i);
    it[i] = new_item(mi[i], set->seqs[i]->description);
  }

  it[i] = new_item("E","Exit");
  it[i+1] = NULL; 
  me = new_menu(it);

  shadoww = newwin(13, 63, 11, 12);
  mainwin = newwin(13, 63, 10, 10);
  win = derwin(mainwin, 11, 60, 1, 1);
  set_menu_win (me, win);
  set_menu_sub (me, derwin(win, 6, 38, 3, 1));
  set_menu_format(me, 5, 1);
  set_menu_mark(me, " * ");

  box(win, 0, 0);  
  mvwaddstr(win, 1, 2, "Select Reference");
  set_menu_fore(me, COLOR_PAIR(5)|A_REVERSE);
  set_menu_back(me, COLOR_PAIR(6));
  set_menu_grey(me, COLOR_PAIR(7));
  wbkgd(shadoww, COLOR_PAIR(8));
  wbkgd(mainwin, COLOR_PAIR(2));
  wbkgd(win, COLOR_PAIR(2));
  post_menu(me);

  attron(COLOR_PAIR(1));
  mvwprintw(win, 9, 2, "Use PageUp, PageDown and Arrows to scroll (F1 exits).");
  attroff(COLOR_PAIR(1));

  refresh();
  wrefresh(shadoww);
  wrefresh(mainwin);

  while(sel > nchr && (ch=getch()) != KEY_F(1))
  {
    switch(ch)
    {
      case KEY_DOWN:
        menu_driver(me, REQ_DOWN_ITEM);
        break;
      case KEY_UP:
        menu_driver(me, REQ_UP_ITEM);
        break;
      case KEY_NPAGE:
        menu_driver(me, REQ_SCR_DPAGE);
        break;
      case KEY_PPAGE:
        menu_driver(me, REQ_SCR_UPAGE);
        break;
      case 0xA: /* Return- bzw. Enter-Taste -> ASCII-Code */
        sel = item_index(current_item(me));
    }

    wrefresh(shadoww);
    wrefresh(mainwin);
  } 

  unpost_menu(me);
  free_menu(me);
                  
  for(i = 0; i < nchr; ++i) {
    free(mi[i]);
    free_item(it[i]);
  }

  free(it);
  free(mi);

  delwin(win);
  delwin(shadoww);
  delwin(mainwin);

  if(sel < nchr) {
    newchrname = ALLOCMEMORY(space, NULL, char, 
        strlen(set->seqs[sel]->description)+1);
    memset(newchrname, 0, 
        strlen(set->seqs[sel]->description)+1);
    memmove(newchrname, set->seqs[sel]->description, 
        strlen(set->seqs[sel]->description));
  }

  return newchrname;
}

/*------------------------ bl_matchfileJumpToInitForm ------------------------
 *    
 * @brief initialize the form to jump to position
 * @author Steve Hoffmann 
 *   
 */

FORM*
bl_matchfileJumpToInitForm(FIELD ***flds) {
  FIELD **fi;
  FORM *fo;
  int i;

  fi = (FIELD **)calloc(3, sizeof(FIELD *));
  fi[0] = new_field(1, 15, 0, COLS-20, 0, 0);
  fi[1] = new_field(1, 15, 0, COLS-45, 0, 0);
  fi[2] = 0;

  set_field_type(fi[0], TYPE_INTEGER, 0, 1, 999999999);
  set_field_type(fi[1], TYPE_REGEXP, "^.*$");

  for(i=0; i < 2; i++) {
    set_field_fore(fi[i], COLOR_PAIR(2));
    set_field_back(fi[i], COLOR_PAIR(2));
    field_opts_on(fi[i], O_EDIT);
  }
  field_opts_off(fi[1], O_AUTOSKIP);
  field_opts_on(fi[1], O_STATIC);


  fo = new_form(fi);
  post_form(fo);

  mvaddstr(0, COLS-25, "pos: ");     
  mvaddstr(0, COLS-50, "chr: ");

  *flds = fi;
  return fo;
}


/*----------------------- bl_matchfileJumpToWrapField ------------------------
 *    
 * @brief destruct the form to jump
 * @author Steve Hoffmann 
 *   
 */

void
bl_matchfileJumpToWrapField(FIELD **fi) {
  free_field(fi[0]);
  free_field(fi[1]);
  free(fi);
  curs_set(0);
}


/*--------------------------- bl_machfileViewQuit ----------------------------
 *    
 * @brief Adios!
 * @author Steve Hoffmann 
 *   
 */

void 
bl_matchfileViewQuit() {
  endwin();
}

/*-------------------------- bl_matchfileGetTrackName ---------------------------
 *    
 * @brief get the track names
 * @author Steve Hoffmann 
 *   
 */

char*
bl_matchfileItemDisplayName (annotationtrack_t *track, Uint k)
{

  Uint i, len, nlen, dlen;
  char *name=NULL, GFFname = 0;
  annotationitem_t *items;

  
  items = track->items; 
  nlen = strlen(items[k].name);

  if(items[k].type == GFFITEM) { 
    for(i=0; i < items[k].noofattributes; i++) {

      if(items[k].attributelen[i] > 5 && 
          (strncmp(items[k].attributes[i], "Name=", 5) == 0 ||
           strncmp(items[k].attributes[i], "name=", 5) == 0)) {             

        if(strlen(items[k].name)) { 
          dlen = strlen(&items[k].attributes[i][5]);
          len = nlen + dlen + 1;
          name = ALLOCMEMORY(space, NULL, char, len+2);
          memset(name, ' ', len);
          memmove(name, items[k].name, nlen);
          memmove(&name[nlen+1], &items[k].attributes[i][5], dlen);
          name[len] = 0;
        }

        GFFname = 1;
        break;
      }
    }
  } 

  if(!GFFname) {  
    name = ALLOCMEMORY(space, NULL, char, nlen+1);
    memmove(name, items[k].name, nlen);
    name[nlen] = 0;
  }

  return name;
}


/*--------------------- bl_matchfileDrawAnnotationTrack ----------------------
 *    
 * @brief draw the annotation track
 * @author Steve Hoffmann 
 *   
 */

  void
bl_matchfileDrawAnnotationTrack (matchfileView_t *view)
{ 
  Uint i=0, k, pos=0, noofitems = 0,
    off; //, start; 
  int  dpos=0, dend=0, dstart=0;
  Uint *imap, curptr, endpos, itemnamelen, p, startitem=0;
  annotationitem_t *items = NULL;
  char *itemname, trackchar;
  int attr, istart, iend, fstart, fend;
  Uint xoff = 3, yoff=0;
  WINDOW *pad;

  attr = A_BOLD | COLOR_PAIR(4) ;

  if(!view->annotation) return;

  if(view->annotationpad) {
    delwin(view->annotationpad);
  }
  
  pad = newpad(MAXPADLINES, COLS);
  off = view->offset;
  // not used: start = view->curframe->start;
  imap = view->imap;

  if(imap[view->offset-1]) { 
    curptr = view->curframe->start + imap[off-1];
  } else { 
    k = view->offset-1;
    while(imap[k] == 0 && k > 0) k--;
    curptr = view->curframe->start + imap[k];
  }

  view->annotationoffset = 0;
  noofitems = view->annotation->noofitems;
  items = view->annotation->items;

  for(k=0; k < noofitems; k++) { 
    if(!strcmp(items[k].chromname, view->curframe->chrname) && 
        items[k].end >= curptr) break;
  }

  startitem = k;

  for(i=0; i < COLS-xoff; i++) {
    if(i+off-1 == 0 || imap[i+off-1] > 0) {  
      pos = view->curframe->start+imap[i+off-1];

      for(k=startitem; k < noofitems; k++) { 
        if(!strcmp(items[k].chromname, view->curframe->chrname)) { 
          if(pos >= items[k].start && pos <= items[k].end) {

            if(items[k].strand == '+') {
              trackchar = '>';
            } else {
              trackchar = '<';
            }

            if(items[k].start == pos) {
              wattrset(pad, COLOR_PAIR(BLUEONYELLOW)); 
              mvwprintw(pad, yoff+items[k].level, i, "|");
            } else if (items[k].end == pos) {
              wattrset(pad, A_BOLD | COLOR_PAIR(REDONYELLOW)); 
              mvwprintw(pad, yoff+items[k].level, i, "|");
            } else { 
              wattrset(pad, A_BOLD | COLOR_PAIR(BLUEONYELLOW)); 
              if(i%2) mvwprintw(pad, yoff+items[k].level, i, "%c", trackchar);
              else mvwprintw(pad, yoff+items[k].level, i, " ");
            }
          }
        }
        if(pos+COLS+10 < items[k].start) {
          break;
        }
      }
    }
  }

  wattrset(pad, COLOR_PAIR(BLUEONYELLOW)); 
  endpos = view->curframe->start+view->map[off-1+COLS-1-(2*xoff)];

  for(k=startitem; k < noofitems; k++) { 

    if(!strcmp(items[k].chromname, view->curframe->chrname)) { 
      if(OVERLAP(items[k].start, items[k].end, curptr, endpos)) {

        itemname = bl_matchfileItemDisplayName(view->annotation, k);
        itemnamelen = strlen(itemname);

        fstart = MAX(0, (int)items[k].start-(int)view->curframe->start);
        fend = items[k].end - view->curframe->start;

        istart = view->map[fstart];
        iend = view->map[fend];

        dstart = MAX((int)istart-(int)off+2, 0); //+ xoff;             
        dend = MIN(COLS-1, iend-(off-1)); //- (2*xoff);

        if(dend - dstart > itemnamelen+2) { 
          dpos = dstart+(((dend-dstart)-itemnamelen)/2);
        } else {  
          dpos = dstart;
        } 

        if(itemnamelen+2 > (dend  - dpos)) {
          for(p=dpos; p < dend && p-dstart < itemnamelen ; p++) {
            mvwaddch(pad, yoff+items[k].level, p, itemname[p-dstart]);
          }
        } else { 
          mvwprintw(pad, yoff+items[k].level, dpos, " %s ", itemname, dstart, dend);
        }

        FREEMEMORY(space, itemname);
      }
    }
  }

  view->annotationpad = pad;
  attrset(A_NORMAL | attr);

  return ;
}

/*-------------------------- bl_matchfileDrawRuler ---------------------------
 *    
 * @brief draw the ruler
 * @author Steve Hoffmann 
 *   
 */


void
bl_matchfileDrawRuler (matchfileView_t *view) {
  Uint i=0, k, pos=0,
    e=0, z=0, h=0, dz=100, dh=1000, off; //, start; 
  int  attr=0;
  Uint *imap, curptr;
  Uint xoff = 3;

  off = view->offset;
  // not used: start = view->curframe->start;
  imap = view->imap;

  if(imap[view->offset-1]) { 
    curptr = view->curframe->start + imap[off-1];
  } else { 
    k = view->offset-1;
    while(imap[k] == 0 && k > 0) k--;
    curptr = view->curframe->start + imap[k];
  }

  attr = A_BOLD | COLOR_PAIR(4) ;
  attrset(attr);

  mvprintw(0 ,3+xoff,"%d (chr: '%s')",  
      curptr, view->curframe->chrname);

  mvaddch (0, 0+xoff, ACS_ULCORNER);

  for(i=0; i < COLS-(2*xoff); i++) {
    if(i+off-1 == 0 || imap[i+off-1] > 0) {  
      pos = view->curframe->start+imap[i+off-1];
      /*ruler*/
      e = (pos)%10;
      attrset(A_UNDERLINE | attr);
      mvprintw(3, i+xoff, "%d", e);
      attrset(A_NORMAL | attr);

      z = (pos+1-e)%100; 
      if(dz != z || !i)
        mvprintw(2, i+xoff, "%d", z/10);
      else
        mvprintw(2, i+xoff, " ");
      dz = z;

      h = (pos+1-e-z)%1000;
      if(dh != h || !i)
        mvprintw(1, i+xoff, "%d", h/100);
      else
        mvprintw(1, i+xoff, " ");
      dh = h;
    }
  }

  attrset(A_NORMAL | attr);
  return;
}



/*------------------------ bl_matchfileViewUpdateFrame -------------------------
 *    
 * @brief load a frame
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchfileViewUpdateFrame(void *space, matchfileView_t *vw, char *chrom, 
    Uint start, Uint width)
{ 
  char *chrname;
  Uint fs=1, fw = 2*width; //, k, i;

  if(vw->curframe && vw->curframe->cs) {  
    bl_matchfileDestructView(space, vw);
  }
    
  if(!chrom) {
    chrname = vw->file->index->chromnames[0]; 
  } else {
    chrname = chrom;
  }

  if(start <= width) {
    vw->offset = start; 
  } else {
    fs = start-width;
    vw->offset = width; 
  }
  /* not used 
  k = bl_matchfileGetChromIndexNumber(vw->file->index, chrname);
  i = bl_fastxFindIDIdx(chrname, vw->set); 
  */
  vw->curframe = bl_matchfileGetFrame(space, vw->file, chrname, fs, 
      fw, vw->set, 20000, NULL);
  
  bl_matchfileGetConsensus(vw->curframe);
  vw->curframestats = bl_matchfileFrameStats(space, vw->curframe);
  
  return ;
}

/*----------------------- bl_matchfileViewUpdateFrame ------------------------
 *    
 * @brief draw the frame to the view
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchfileViewDrawFrame(void *space, matchfileView_t *view, 
    matchfilePanel_t *panel, Uint frameno) {

  Uint i=0, v, j, k, l, m, starts, cover;
  Uint *imap, *map, maxdel=0;
  char *dels;
  matchfileCross_t *cs, *xs;
  int col, curcol=0, errcol=0;

  delwin(view->pad);   
  view->pad = newpad(MAXPADLINES, view->curframe->width+2);
  wclear(view->pad);
  assert(view->pad);


  cs = view->curframe->cs;

  map = calloc(view->curframe->width*10, sizeof(Uint));
  imap = calloc(view->curframe->width*10, sizeof(Uint));

  for(k=0, i=0; i < view->curframe->width; i++) {
    map[i] = k-maxdel;
    imap[k] = i;

    if(view->curframe->ref) { 

      bl_matchfileTest(space, 0, 0, view->curframe->start+i, &view->curframe->cs[i], 
          view->curframe->ref[i], view->idx, 0, NULL);
      /*
         if(!isinf(cs[i].p_hom)) { 
         errcol = COLOR_PAIR(REDONWHITE);
         curcol = COLOR_PAIR(BLUEONWHITE);
         } else { 
         if (cs[i].s_cons >= cs[i].s_consx && cs[i].s_ref >= cs[i].s_refx) { 
         curcol = COLOR_PAIR(WHITEONBLACK);
         errcol = COLOR_PAIR(REDONBLACK);
         } else { 
         if(cs[i].cons == view->curframe->ref[i] || 
         (cs[i].s_ref <  cs[i].s_refx && cs[i].s_cons < cs[i].s_consx)) { 
         errcol = COLOR_PAIR(REDONYELLOW);
         curcol = COLOR_PAIR(BLUEONYELLOW);
         } else {
         errcol = COLOR_PAIR(REDONBLUE); 
         curcol = COLOR_PAIR(WHITEONBLUE);

         }
         }
         }
         */ 
      if( (cs[i].p_consx != log(0) && 
            cs[i].p_consx > cs[i].p_cons) || 
          (cs[i].p_refx != log(0) && 
           cs[i].p_refx > cs[i].p_ref) ||
          !isinf(cs[i].p_hom)) {

        if(!isinf(cs[i].p_hom)) {
          errcol = COLOR_PAIR(REDONWHITE);
          curcol = COLOR_PAIR(BLUEONWHITE);
        } else if(cs[i].p_consx > cs[i].p_cons && 
            cs[i].p_refx > cs[i].p_ref) { 
          errcol = COLOR_PAIR(REDONYELLOW);
          curcol = COLOR_PAIR(BLUEONYELLOW);
        } else { 
          errcol = COLOR_PAIR(REDONBLUE); 
          curcol = COLOR_PAIR(WHITEONBLUE);
        }

      } else {
        curcol = COLOR_PAIR(WHITEONBLACK);
        errcol = COLOR_PAIR(REDONBLACK);
      }
    }

    if(view->curframe->ref){
      mvwaddchattr(view->pad, 0, k, A_DIM | curcol, 
          view->curframe->ref[i]);
    }

    col = (view->curframe->ref && cs[i].cons != view->curframe->ref[i]) 
      ? errcol : curcol;
    mvwaddchattr(view->pad, 1, k, A_UNDERLINE | col, cs[i].cons);

    starts = (cs[i].starts/10 > 9) ? 9 : cs[i].starts/10;
    cover =  (cs[i].len/10 > 9) ? 9 : cs[i].len/10;

    wattrset(view->pad, curcol);
    mvwprintw(view->pad, 2, k, "%d", cover);  
    wattrset(view->pad, curcol | A_UNDERLINE);
    mvwprintw(view->pad, 3, k, "%d", starts);  
    dels = calloc((cs[i].maxrow*2)+1, sizeof(char));

    /*determine max deletion accros all views*/
    for(v=0, maxdel=0; v < panel->noofviews; v++) {
      xs = panel->views[v]->curframe->cs;
      for(l=0; l < xs[i].noofdels; l++) {
        maxdel = MAX(maxdel, xs[i].dels[l].len);
      }
    }

    for(l=0; l < cs[i].noofdels; l++) {
      dels[cs[i].dels[l].row] = 1;
    }


    for(j=0; j < cs[i].len; j++) {

      if(cs[i].row[j] < MAXPADLINES){
        col = (view->curframe->ref && 
            cs[i].chars[j] != view->curframe->ref[i]) ? errcol : curcol;
        mvwaddchattr(view->pad, cs[i].row[j]+4, k, col, cs[i].chars[j]);  
      }

      assert(cs[i].row[j] <= cs[i].maxrow);
      wattrset(view->pad, A_BOLD | COLOR_PAIR(4));

      if(dels && dels[cs[i].row[j]]) {
        for(l=0; l < cs[i].noofdels; l++) {
          if(cs[i].dels[l].row == cs[i].row[j]) {
            break;
          }
        }

        for(m=0; m < cs[i].dels[l].len; m++) {
          mvwaddch(view->pad, cs[i].row[j]+4, k+m+1, 
              cs[i].dels[l].string[m]);  
        }

        for(;m < maxdel;m++) {
          mvwaddch(view->pad, cs[i].row[j]+4, k+m+1, '^');  
        }

      } else {

        for(m=0; m < maxdel; m++) {
          if(cs[i].feat[j] != '$'){ 
            mvwaddch(view->pad, cs[i].row[j]+4, k+m+1, '^');  
          }
        }
      }
    }


    FREEMEMORY(space, dels);
    k+=maxdel+1;
  }

  view->offset = map[view->offset];
  view->map = map;
  view->imap = imap;
}

/*----------------------- bl_matchfileViewDrawPanel ------------------------
 *    
 * @brief update the panel
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchfileViewDrawPanel (void *space, matchfilePanel_t *panel)
{

  Uint i;

  for(i=0; i < panel->noofviews; i++) { 
    bl_matchfileViewDrawFrame(space, panel->views[i], panel, i);
  }

  return ;
}

/*----------------------- bl_matchfileViewUpdatePanel ------------------------
 *    
 * @brief update the panel
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchfileViewUpdatePanel (void *space, matchfilePanel_t *panel, 
    char *chrom, Uint start, Uint width)
{

  Uint i;

  for(i=0; i < panel->noofviews; i++) { 
    bl_matchfileViewUpdateFrame(space, panel->views[i], chrom, start, width); 
  }



  return ;
}

/*----------------------- bl_matchfileViewRefreshPanel -----------------------
 *    
 * @brief panel refresh
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchfileRefreshViewPanel (void *space, matchfilePanel_t *panel)
{
  Uint i, lpp = MAX(10, ((LINES-1-8)/(panel->noofviews)));
  Uint al = 0;
  Uint xoff = 3;


  if(panel->views[0]->annotation) {
    al = 5; 
  }
  
  refresh();


  for(i=0; i < panel->noofviews; i++) { 
    
    if(i+1 == panel->activeview){ 
      delwin(panel->activebox);
      
      panel->activebox = newwin(lpp, COLS-2,
          4+al+(i*lpp)+1, panel->smincol[i]+1);

      wbkgd(panel->activebox, COLOR_PAIR(1));
      box(panel->activebox, ACS_VLINE, ACS_HLINE);  
    }
    
    wrefresh(panel->activebox);
    
      /*   
        prefresh(panel->views[i]-header, 5, panel->views[i]->offset-1, 
        4+al+(i*lpp)+3, panel->smincol[i]+3, 
        4+al+((i+1)*lpp)-3, COLS-1-3);
       
       */

    prefresh(panel->views[i]->pad, panel->views[i]->scroll, 
        panel->views[i]->offset-1, 
        4+al+(i*lpp)+2, panel->smincol[i]+3, 
        4+al+((i+1)*lpp)-3, COLS-1-3);

    refresh();
  }

  if(panel->views[0]->annotation) {  
    if(panel->activeview == 0){ 
      delwin(panel->activebox);
      panel->activebox = newwin(7, COLS-2, 4, 1);

      wbkgd(panel->activebox, COLOR_PAIR(1));
      box(panel->activebox, ACS_VLINE, ACS_HLINE);  
    }

    wrefresh(panel->activebox);
    prefresh(panel->views[0]->annotationpad, 
        panel->views[0]->annotationscroll, 0, 5, 3, 9, COLS-1-xoff);
    refresh();
  }

  return ;
}

/*---------------------- bl_matchfileControlScreenSize -----------------------
 *    
 * @brief control the size of the screen
 * @author Steve Hoffmann 
 *   
 */
 
int
bl_matchfileControlScreenSize ( )
{

  if (COLS-1 < 10) {
    mvprintw(0,0,"screensize!");
    return 0;
  }

  return 1;
}


/*-------------------------- bl_matchfileViewUpdate --------------------------
 *    
 * @brief user interface 
 * @author Steve Hoffmann 
 *   
 */

void
bl_matchfileViewUpdateScreen(void *space, matchfilePanel_t *panel, 
    Uint width, fasta_t *set) {
  int ch, 
      nch;
  int curexit = 'q', redraw;
  Uint chrnamelen, pos, k, i1, margin, refend, vwidth, vstart, j, chrlen;
  FIELD **fi;
  FORM *fo;
  SHADOWEDWINDOW *infowin = NULL;
  unsigned char showcsinfo=0;
  char *chrname = NULL, 
       *newchrname;
  matchfileView_t **views;


  bl_matchfileDrawRuler(panel->views[0]);
  bl_matchfileDrawAnnotationTrack (panel->views[0]);

  bl_matchfileRefreshViewPanel (space, panel);
  views = panel->views;

  while((ch=getch()) != curexit) {

    if(!bl_matchfileControlScreenSize()) continue;

    /*
     *  key up jump +100 nucleotides
     *
     * */

    if(ch == 'a' || ch == KEY_LEFT || ch == KEY_UP) {
      clear();

      redraw = 0;
      for(j=0; j < panel->noofviews; j++) { 
        if(ch != KEY_UP) {
          views[j]->offset++;
        } else {  
          views[j]->offset += 100;
        }

        vstart = views[j]->curframe->start;
        vwidth = views[j]->curframe->width;
        margin = views[j]->offset + COLS + 100;
        refend = vstart + margin;
        chrlen = views[j]->curframe->chrlen;

        if(margin >= vwidth &&  refend < MIN(chrlen, MAXCHROMSIZE)) {        
          
          k = views[j]->offset;
          
          if(views[j]->imap[k]) { 
            pos = vstart + views[j]->imap[k];
          } else { 
            while(views[j]->imap[k] == 0 && k > 0) k--;
            pos = vstart + views[j]->imap[k];
          }

          bl_matchfileViewUpdateFrame(space, views[j], chrname, pos, width); 
          redraw = 1;
          clear();
        }
      }

      if(redraw) bl_matchfileViewDrawPanel (space, panel);
    }

    /*
     *  key down jump -100 nucleotides
     *
     * */

    if(ch == 'd' || ch == KEY_RIGHT || ch == KEY_DOWN) {
      clear();

      redraw = 0;
      for(j=0; j < panel->noofviews; j++) { 

        if(ch != KEY_DOWN) {
          views[j]->offset--;
        } else {
          if (views[j]->offset < 100) {
            views[j]->offset = 1;
          } else {
            views[j]->offset -= 100;
          }
        }
 
        vstart = views[j]->curframe->start;

        if(views[j]->offset <= 1 && vstart > 1 
            && views[j]->offset + vstart < MAXCHROMSIZE) {

          bl_matchfileViewUpdateFrame(space, views[j], chrname, vstart, width);
          redraw = 1;
          clear();
        }
        if (views[j]->offset < 1) views[j]->offset = 1;
      }

      if(redraw) bl_matchfileViewDrawPanel (space, panel);
      bl_matchfileDrawRuler(views[0]);
      bl_matchfileDrawAnnotationTrack (views[0]);
    }

    /*
     *  menu to change chromosomes
     *
     **/

    if(ch == 'm') {
      newchrname = bl_matchfileSelectChrMenu(set);
      if(newchrname) {
        if(chrname) {
          FREEMEMORY(space, chrname);
        }
        chrname = newchrname;

        for(j=0; j < panel->noofviews; j++) {
          bl_matchfileViewUpdateFrame(space, views[j], chrname, 100, width);
        }
        clear();
        bl_matchfileViewDrawPanel (space, panel);
      }   
    }

    if(ch == KEY_NPAGE) {      
      if(panel->activeview > 0 && panel->views[panel->activeview-1]->scroll+1 < MAXPADLINES-5)
        panel->views[panel->activeview-1]->scroll++;

      if(!panel->activeview && panel->views[0]->annotationscroll+1 < MAXPADLINES-5) {
        panel->views[0]->annotationscroll++;
      }
    }
  
    if(ch == KEY_PPAGE) {
      if(panel->activeview > 0 && panel->views[panel->activeview-1]->scroll > 0)
        panel->views[panel->activeview-1]->scroll--;

      if(!panel->activeview && panel->views[0]->annotationscroll > 0) {
        panel->views[0]->annotationscroll--;
      }
    }

    /*
     *  change active panel
     *
     **/

    if (ch == 'n') {
      if(panel->activeview+1 <= panel->noofviews)
        panel->activeview++;
      clear();
    }

    if (ch == 'p') {
      if(panel->activeview >= 1)
        panel->activeview--;
      clear();
    }


    if(ch == 'P') {
      bl_matchfilePERRGNUPLOT(space, views[0]->file->index);
    }

    if(ch == 'E') {
      bl_matchfileQERRGNUPLOT(space, views[0]->file->index);
    }
    
    if(ch == 'C') {
      bl_matchfileCOVGNUPLOT(space, views[0]->curframe);
    }

    if(ch == 'W'){
      bl_matchfileSUBGNUPLOT(space, views[0]->file->index);
    }

    if(ch == 'S') {
      bl_matchfileRSSGNUPLOT(space, views[0]->curframe, 
          views[0]->curframestats);
    }
    
    if(ch == 'i' || (showcsinfo && ch == 'q')) {
      if(showcsinfo) {
        delshadowedwin(infowin);
        infowin = NULL;
        showcsinfo = 0;
        curexit = 'q';
      } else {
        showcsinfo = 1;
        curexit = 0;
      }
    }

    if(ch == ':') {

      clear();
      fo = bl_matchfileJumpToInitForm(&fi);   
      bl_matchfileDrawRuler(views[0]);    
      bl_matchfileDrawAnnotationTrack (views[0]);
      bl_matchfileRefreshViewPanel (space, panel);

      curs_set(1);
      move(0,COLS-20);

      while((nch=wgetch(stdscr)) != KEY_F(1)) {
        switch(nch) {
          case KEY_BTAB:
            form_driver(fo, REQ_END_LINE);
            form_driver(fo, REQ_PREV_FIELD);
            form_driver(fo, REQ_END_LINE);
            break;
          case 9:
            form_driver(fo, REQ_END_LINE);
            form_driver(fo, REQ_NEXT_FIELD);
            form_driver(fo, REQ_END_LINE);
            break;
          case KEY_LEFT:
          case '\b':
          case KEY_BACKSPACE:
            form_driver(fo, REQ_DEL_PREV);
            break;
          case '\n':
            form_driver(fo, REQ_END_LINE);
            form_driver(fo, REQ_CLR_FIELD);
            break;      
          default:
            form_driver(fo, nch);
        } 
        bl_matchfileRefreshViewPanel (space, panel);
        if(nch == '\n') break;
      }

      curs_set(0);
      i1 = atoi(field_buffer(fi[0], 0));
      chrnamelen = strlen(field_buffer(fi[1],0));
      newchrname = strtrim(NULL, field_buffer(fi[1],0), &chrnamelen);

      if(newchrname) {
        if(chrname) {
          FREEMEMORY(space, chrname);
        }
        chrname = newchrname;

        for(j=0; j < panel->noofviews; j++) {

          bl_matchfileViewUpdateFrame(space, views[j], chrname, 100, width);
        }
        bl_matchfileViewDrawPanel (space, panel);
      }

      unpost_form(fo);
      free_form(fo);
      bl_matchfileJumpToWrapField(fi);     

      if(i1 && i1 >= 0 && i1 < MAXCHROMSIZE) {
        for(j=0; j < panel->noofviews; j++) {
          bl_matchfileViewUpdateFrame(space, views[j], chrname, i1, width);
        }      
        bl_matchfileViewDrawPanel (space, panel);
      }

      for(j=0; j < panel->noofviews; j++) {
        if(views[j]->offset < 1) views[j]->offset = 1; 
      }
    }

    bl_matchfileDrawRuler(views[0]);
    bl_matchfileDrawAnnotationTrack (views[0]);
    bl_matchfileRefreshViewPanel (space, panel);

    if(showcsinfo) {
      delshadowedwin(infowin);
      infowin = bl_matchfileGetInfo(views[0]);
      refresh();
      shadowedwrefresh(infowin);
    }
  }

  if(chrname) {
    FREEMEMORY(space, chrname);
  }
}


/*--------------------------- bl_matchfileViewInit ---------------------------
 *    
 * @brief initialize the match file viewer
 * @author Steve Hoffmann 
 *   
 */

void
bl_matchfileViewInit() {
  initscr();
  atexit(bl_matchfileViewQuit);
  clear();
  noecho();
  curs_set(0);
  cbreak();
  keypad(stdscr, TRUE);
  start_color();

  init_pair(1, COLOR_WHITE, COLOR_BLACK);
  init_pair(2, COLOR_BLACK, COLOR_WHITE); 
  init_pair(3, COLOR_RED, COLOR_BLACK);
  init_pair(4, COLOR_BLUE, COLOR_BLACK);
  init_pair(5, COLOR_WHITE, COLOR_BLUE);
  
  init_pair(6, COLOR_BLUE, COLOR_GREEN);

  init_pair(7, COLOR_BLACK, COLOR_BLUE); 
  init_pair(8, COLOR_BLACK, COLOR_BLACK);
  init_pair(9, COLOR_BLUE, COLOR_YELLOW);
  init_pair(10, COLOR_RED, COLOR_YELLOW);
  init_pair(11, COLOR_BLUE, COLOR_MAGENTA);

  init_pair(12, COLOR_WHITE, COLOR_GREEN);
  init_pair(13, COLOR_RED, COLOR_GREEN);
  init_pair(14, COLOR_WHITE, COLOR_YELLOW);

  init_pair(15, COLOR_RED, COLOR_WHITE);
  init_pair(16, COLOR_BLUE, COLOR_WHITE);
  init_pair(17, COLOR_GREEN, COLOR_WHITE);
  init_pair(18, COLOR_MAGENTA, COLOR_WHITE);
  init_pair(19, COLOR_CYAN, COLOR_WHITE);
  
  init_pair(20, COLOR_GREEN, COLOR_YELLOW);
  init_pair(21, COLOR_RED, COLOR_BLUE);

}



/*------------------------- bl_matchfileDestructView -------------------------
 *    
 * @brief remove the view from the heap
 * @author Steve Hoffmann 
 *   
 */
 

void
bl_matchfileDestructView(void *space, matchfileView_t *view) {

  bl_matchfileDestructCross(space, view->curframe->cs, view->curframe->width); 
  bl_matchfileDestructFrameStats(space, view->curframestats);
  
  FREEMEMORY(space, view->curframe->cs);
  FREEMEMORY(space, view->curframestats);
  FREEMEMORY(space, view->curframe);
  FREEMEMORY(space, view->imap);
  FREEMEMORY(space, view->map);


}

/*-------------------------- bl_matchfileInitView ----------------------------
 *    
 * @brief init the views
 * @author Steve Hoffmann 
 *   
 */
 
matchfileView_t*
bl_matchfileInitView(void *space, matchfile_t *file, 
    matchfileindex_t *idx, fasta_t *set, Uint fw)
{
  matchfileView_t *view;
  matchfileSampleStats_t *stats = idx->stats;

  assert(file->index);
  assert(file->index->chromnames);
  assert(file->index->chromnames[0]);
  
  view = ALLOCMEMORY(space, NULL, matchfileView_t, 1);
  view->file = file;
  view->stats = stats;
  view->offset = 0;
  view->set = set;
  view->map = NULL;
  view->imap = NULL;
  view->curframe = NULL;
  view->curframestats = NULL;	
  view->scroll = 0;
  view->pad = NULL;
  view->annotation = NULL;
  view->annotationoffset = 0;
  view->scroll = 0;
  view->annotationscroll = 0;
  view->pad = NULL;
  view->annotationpad = NULL;
  return view;
}

/*------------------------ bl_matchfileInitViewPanel -------------------------
 *    
 * @brief init the view panel
 * @author Steve Hoffmann 
 *   
 */
 
matchfilePanel_t*
bl_matchfileInitViewPanel(void *space, matchfile_t **files, 
    Uint nooffiles, fasta_t *set, annotationtrack_t *track, Uint fw)
{ 

  Uint i, nrow = 0, rowsum = 0;;
  matchfilePanel_t *panel;
  matchfileView_t **views;


  panel = ALLOCMEMORY(space, NULL, matchfilePanel_t, 1);
  views = ALLOCMEMORY(space, NULL, matchfileView_t*, nooffiles);

  panel->activeview = 0;
  panel->activebox = NULL;
  panel->pminrow = ALLOCMEMORY(space, NULL, int, nooffiles);
  panel->pmincol = ALLOCMEMORY(space, NULL, int, nooffiles);
  panel->sminrow = ALLOCMEMORY(space, NULL, int, nooffiles);
  panel->smincol = ALLOCMEMORY(space, NULL, int, nooffiles);
  panel->smaxrow = ALLOCMEMORY(space, NULL, int, nooffiles);
  panel->smaxcol = ALLOCMEMORY(space, NULL, int, nooffiles);
 
  nrow = MIN(20, (LINES-1-4/nooffiles));

  for(i=0; i < nooffiles; i++) {
    views[i] = bl_matchfileInitView(space, files[i], 
        files[i]->index, set, fw);

    views[i]->annotation = track;

    views[i]->pad = newpad(MAXPADLINES, fw);
    assert(views[i]->pad);

    panel->pminrow[i] = 5;
    panel->pmincol[i] = 0;
    panel->sminrow[i] = rowsum+4;
    panel->smincol[i] = 0;
    panel->smaxrow[i] = rowsum+4+nrow;
    panel->smaxcol[i] = COLS-1;
    rowsum += nrow;
  }

  panel->noofviews = nooffiles;
  panel->views = views;
  
  return panel;
}


/*---------------------- bl_matchfileDestructViewPanel -----------------------
 *    
 * @brief destruct the view panel
 * @author Steve Hoffmann 
 *   
 */
 
void
bl_matchfileDestructViewPanel(void *space, matchfilePanel_t *panel)
{
  Uint i;

  for(i=0; i < panel->noofviews; i++) { 
    bl_matchfileDestructView(space, panel->views[i]);
    FREEMEMORY(space, panel->views[i]);
  }

  FREEMEMORY(space, panel->pminrow);
  FREEMEMORY(space, panel->pmincol);
  FREEMEMORY(space, panel->sminrow);
  FREEMEMORY(space, panel->smincol);
  FREEMEMORY(space, panel->smaxrow);
  FREEMEMORY(space, panel->smaxcol);
  FREEMEMORY(space, panel->views);
  return ;
}

/*----------------------------- matchfileViewer ------------------------------
 *    
 * @brief start the browser
 * @author Steve Hoffmann 
 *   
 */

void
bl_matchfileViewer(void *space, matchfile_t **files, Uint nooffiles, 
    fasta_t *set, annotationtrack_t *track, Uint start, Uint width) {

  matchfilePanel_t *panel;
  Uint fw;
    

  bl_matchfileViewInit();
  while(!bl_matchfileControlScreenSize());

  fw = (start<= width) ? width : 2*width;
  panel = bl_matchfileInitViewPanel(space, files, nooffiles, set, track, fw);

  bl_matchfileViewUpdatePanel(space, panel, NULL, start, width);
  bl_matchfileViewDrawPanel (space, panel);
 // bl_matchfileDumpSampleStats (panel->views[0]->stats);
  bl_matchfileViewUpdateScreen(space, panel, width, set);
  
  bl_matchfileDestructViewPanel(space, panel);
  FREEMEMORY(space, panel);
}

