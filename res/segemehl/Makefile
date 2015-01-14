  CC=gcc
  LD=${CC} 
  CFLAGS= -Wall -pedantic -std=c99 -g -O3 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DDBGNFO -DSHOWALIGN -DDBGLEVEL=0 -DPROGNFO -Isrc -Ilibs -Ilibs/sufarray -Lsrc
  LDFLAGS= -lm -lpthread -lz -lncurses 
  CTAGS=ctags > tags
  LIBS=-lob -lm -lpthread 



  GENFASTAOBJ = libs/vtprogressbar.o\
                libs/debug.o\
				libs/info.o\
				libs/seqclip.o\
				libs/stringutils.o\
				libs/mathematics.o\
				libs/fileio.o\
				libs/zran.o\
				libs/biofiles.o\
				libs/randseqs.o\
				libs/memory.o\
				libs/alignment.o\
				libs/sw.o\
				libs/sufarray/charsequence.o\
				libs/manopt.o\
				libs/iupac.o\
				libs/bitVector.o\
				src/genfasta.o

  SEGEMEHLOBJ =	libs/debug.o\
				libs/info.o\
				libs/md5.o\
				libs/stringutils.o\
				libs/708.o\
				libs/mathematics.o\
                libs/memory.o\
				libs/zran.o\
				libs/fileio.o\
				libs/biofiles.o\
				libs/karlin.o\
				libs/kdseed.o\
				libs/queue.o\
				libs/sufarray/sufarray.o\
				libs/sufarray/charsequence.o\
				libs/sufarray/multicharseq.o\
			    libs/radixsort.o\
			    libs/fileBins.o\
			 	libs/stack.o\
				libs/aluruSort.o\
				libs/sort.o\
				libs/vtprogressbar.o\
				libs/sufarray/mmchar.o\
                libs/bitArray.o\
				libs/bitVector.o\
				libs/bitvectoralg.o\
                libs/manout.o\
				libs/kdchain.o\
				libs/manopt.o\
				libs/container.o\
				libs/vstack.o\
				libs/vqueue.o\
				libs/alignment.o\
				libs/sw.o\
				libs/seqclip.o\
				libs/iupac.o\
				libs/hash.o\
				libs/matchfilesfields.o\
				libs/merge.o\
				src/kdmatch.o\
				src/segemehl.o


  TESTMATOBJ =  libs/vtprogressbar.o\
				libs/debug.o\
				libs/info.o\
				libs/bitVector.o\
				libs/seqclip.o\
				libs/stringutils.o\
				libs/mathematics.o\
				libs/splines.o\
				libs/snvsplines.o\
				libs/evalmatchfileshelper.o\
				libs/708.o\
				libs/fileio.o\
				libs/zran.o\
				libs/biofiles.o\
				libs/randseqs.o\
				libs/sufarray/charsequence.o\
				libs/memory.o\
				libs/iupac.o\
				libs/alignment.o\
				libs/sw.o\
				libs/nw.o\
				libs/manopt.o\
				libs/container.o\
				libs/vstack.o\
				libs/vqueue.o\
			    libs/sort.o\
				libs/list.o\
				libs/matepairs.o\
				libs/matchfilesfields.o\
				libs/matchfiles.o\
				libs/evalmatchfiles.o\
				libs/evalmethylmatchfiles.o\
				libs/browsematchfiles.o\
				libs/plotmatchfiles.o\
				libs/ncursesext.o\
				libs/splicesites.o\
				libs/startsites.o\
				libs/matfile.o

  TESTREALIGNOBJ =  libs/fileBins.o\
				libs/manout.o\
				libs/evalmatchfileshelper.o\
				libs/708.o\
				libs/debug.o\
				libs/info.o\
				libs/bitVector.o\
				libs/seqclip.o\
				libs/stringutils.o\
				libs/mathematics.o\
				libs/fileio.o\
				libs/zran.o\
				libs/biofiles.o\
				libs/randseqs.o\
				libs/sufarray/charsequence.o\
				libs/sufarray/multicharseq.o\
				libs/memory.o\
				libs/iupac.o\
				libs/alignment.o\
				libs/sw.o\
				libs/nw.o\
				libs/manopt.o\
				libs/container.o\
				libs/vstack.o\
				libs/vqueue.o\
			    libs/radixsort.o\
				libs/sort.o\
				libs/list.o\
				libs/fqueue.o\
				libs/matepairs.o\
				libs/matchfilesfields.o\
				libs/matchfiles.o\
				libs/vtprogressbar.o\
				libs/evalmatchfiles.o\
				libs/splicesites.o\
				libs/startsites.o\
				libs/realign.o

  TESTREMAPPINGOBJ =  libs/vtprogressbar.o\
			    libs/fileBins.o\
				libs/manout.o\
				libs/debug.o\
				libs/info.o\
				libs/bitVector.o\
				libs/seqclip.o\
				libs/stringutils.o\
				libs/mathematics.o\
				libs/evalmatchfileshelper.o\
				libs/708.o\
				libs/fileio.o\
				libs/zran.o\
				libs/biofiles.o\
				libs/randseqs.o\
				libs/sufarray/charsequence.o\
				libs/sufarray/multicharseq.o\
				libs/memory.o\
				libs/iupac.o\
				libs/alignment.o\
				libs/sw.o\
				libs/nw.o\
				libs/manopt.o\
				libs/container.o\
				libs/vstack.o\
				libs/vqueue.o\
			    	libs/radixsort.o\
				libs/sort.o\
				libs/list.o\
				libs/fqueue.o\
				libs/matepairs.o\
				libs/matchfilesfields.o\
				libs/matchfiles.o\
				libs/evalmatchfiles.o\
				libs/ncursesext.o\
				libs/splicesites.o\
				libs/startsites.o\
				libs/realign.o\
				libs/remapping.o

all: segemehl.x lack.x testrealign.x

haarz.x: ${TESTMATOBJ}
					gcc $(CFLAGS) ${TESTMATOBJ} -o $@ $(LDFLAGS) -lform -lmenu -lgsl -lgslcblas


segemehl.x: ${SEGEMEHLOBJ}
					echo " " >> src/version.h
					gcc $(CFLAGS) ${SEGEMEHLOBJ} -o $@ $(LDFLAGS)

genfasta.x: ${GENFASTAOBJ}
					gcc $(CFLAGS) ${GENFASTAOBJ} -o $@ $(LDFLAGS)


testrealign.x: ${TESTREALIGNOBJ}
					gcc $(CFLAGS) -DREALIGNTEST -c -o libs/realign.o libs/realign.c
					gcc $(CFLAGS) ${TESTREALIGNOBJ} -o $@ $(LDFLAGS)

lack.x: ${TESTREMAPPINGOBJ}					
					gcc $(CFLAGS) -c -o libs/realign.o libs/realign.c
					gcc $(CFLAGS) -DREMAPPINGTEST -c -o libs/remapping.o libs/remapping.c
					gcc $(CFLAGS) ${TESTREMAPPINGOBJ} -o $@ $(LDFLAGS)

clean: 
	rm -f ${TESTMATOBJ} 
	rm -f ${SEGEMEHLOBJ} ${GENFASTAOBJ}
	rm -f ${TESTREALIGNOBJ} 
	rm -f *~
	rm -f libs/*~
	rm -f libs/sufarray/*~
	rm -f src/*~



