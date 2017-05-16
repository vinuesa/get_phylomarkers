# Replaces squids own automake files
#
# Some useful automake documentation:
# http://www.openismus.com/documents/linux/automake/automake.shtml
# http://www.bioinf.uni-freiburg.de/~mmann/HowTo/automake.html

#SUBDIRS =
EXTRA_DIST = COPYRIGHT LICENSE 00README clustalo.README squid-1.9g.tar.gz
#LDADD =

# need DCLUSTALO here to enable Clustal specific squid patches
# CLUSTALO is also defined in ../config.h which is included by squid.h
# but some files do not include squid.h
# not using -DSRE_STRICT_ANSI allows to use gzipped files as input and lists user *and* system time during timings
# not sure if there are any downsides
# AM_CFLAGS = -DCLUSTALO -DSRE_STRICT_ANSI @AM_CFLAGS@ 
AM_CFLAGS = -DCLUSTALO @AM_CFLAGS@ 

noinst_LTLIBRARIES = libsquid.la

libsquid_la_SOURCES =  a2m.c aligneval.c alignio.c \
	clustal.c cluster.c \
	dayhoff.c \
	eps.c \
	file.c \
	getopt.c gki.c gki.h gsi.c gsi.h gsi64.c gsi64.h \
	hsregex.c \
	iupac.c \
	msa.c msa.h msf.c \
	phylip.c \
	revcomp.c rk.c rk.h \
	selex.c seqencode.c shuffle.c sqerror.c sqfuncs.h sqio.c squid.h squidcore.c \
	sre_ctype.c sre_math.c sre_random.c sre_random.h sre_string.c \
	ssi.c ssi.h stack.c stockholm.c stockholm.h stopwatch.c stopwatch.h \
	translate.c types.c vectorops.c \
	vectorops.h version.h \
	weight.c

library_includedir=$(includedir)/clustalo/squid

library_include_HEADERS = gki.h gsi64.h msa.h rk.h \
	squid.h stopwatch.h sqfuncs.h \
	sre_random.h ssi.h stockholm.h \
	vectorops.h version.h


