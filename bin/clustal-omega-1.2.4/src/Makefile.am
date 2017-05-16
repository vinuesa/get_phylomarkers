EXTRA_DIST = clustalo-api-test.c

SUBDIRS = hhalign squid kmpp clustal

AM_LDFLAGS = @AM_LDFLAGS@
AM_CXXFLAGS = @AM_CXXFLAGS@ @OPENMP_CXXFLAGS@
AM_CFLAGS = @AM_CFLAGS@ @OPENMP_CFLAGS@


# the following is the clustalo convenience library which wraps all
# the other noinst libraries

lib_LTLIBRARIES = libclustalo.la

libclustalo_la_SOURCES = clustal-omega.c
libclustalo_la_LIBADD = clustal/libclustal.la \
	hhalign/libhhalign.la \
	squid/libsquid.la	\
	kmpp/libkmpp.la 
libclustalo_la_LDFLAGS = -static
library_includedir=$(includedir)/clustalo/
library_include_HEADERS = clustal-omega-config.h clustal-omega.h



bin_PROGRAMS = clustalo

clustalo_SOURCES =  main.cpp \
	mymain.c mymain.h

clustalo_LDADD = libclustalo.la
