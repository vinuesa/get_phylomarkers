# Some useful automake documentation:
# http://www.openismus.com/documents/linux/automake/automake.shtml
# http://www.bioinf.uni-freiburg.de/~mmann/HowTo/automake.html

#SUBDIRS = 
#EXTRA_DIST =
# LDADD = foofiles/libfoo.a

AM_CXXFLAGS = -fno-strict-aliasing -DCLUSTALO -DCLUSTALO_NOFILE -DDEFAULT_FILTER=90 @OPENMP_CXXFLAGS@ @AM_CXXFLAGS@ 

noinst_LTLIBRARIES = libhhalign.la

libhhalign_la_SOURCES =  hhalign.cpp hhalign.h \
	hash-C.h hash.h hhdecl-C.h hhfunc-C.h hhfunc.h \
	hhmatrices-C.h hhutil-C.h \
	hhalignment-C.h hhalignment.h \
	hhfullalignment-C.h hhfullalignment.h \
	hhhalfalignment-C.h hhhalfalignment.h \
	hhhit-C.h hhhit.h hhhitlist-C.h hhhitlist.h \
	hhhmm-C.h hhhmm.h \
	util-C.h list-C.h list.h \
	general.h

library_includedir=$(includedir)/clustalo/hhalign

library_include_HEADERS = general.h hhfunc.h
