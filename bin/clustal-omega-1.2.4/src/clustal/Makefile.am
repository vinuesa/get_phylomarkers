AM_CFLAGS = @OPENMP_CFLAGS@ @AM_CFLAGS@
AM_LDFLAGS = @AM_LDFLAGS@

noinst_LTLIBRARIES = libclustal.la

libclustal_la_SOURCES = \
        hhalign_wrapper.c hhalign_wrapper.h \
		ktuple_pair.c ktuple_pair.h \
        list.c list.h \
		log.c log.h \
        muscle_upgma.c muscle_upgma.h \
        muscle_tree.c muscle_tree.h \
        mbed.c mbed.h \
        pair_dist.h pair_dist.c \
        progress.c progress.h \
        queue.h \
        seq.c seq.h \
        symmatrix.c symmatrix.h \
        tree.c tree.h \
        util.c util.h \
        weights.c weights.h

library_includedir=$(includedir)/clustalo/clustal

library_include_HEADERS = \
	hhalign_wrapper.h \
	ktuple_pair.h \
	log.h \
	mbed.h \
	muscle_tree.h \
	pair_dist.h \
	progress.h \
	seq.h \
	symmatrix.h \
	tree.h \
	util.h \
	weights.h

