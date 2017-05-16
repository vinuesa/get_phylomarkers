# Replaces squids own automake files
#
# Some useful automake documentation:
# http://www.openismus.com/documents/linux/automake/automake.shtml
# http://www.bioinf.uni-freiburg.de/~mmann/HowTo/automake.html

#SUBDIRS =
#EXTRA_DIST = 
#LDADD =

# need DCLUSTALO here to enable Clustal specific squid patches
AM_CXXFLAGS = -DCLUSTALO @AM_CXXFLAGS@ 

noinst_LTLIBRARIES = libkmpp.la

libkmpp_la_SOURCES =  KMeans.cpp KMeans.h \
	KmTree.cpp KmTree.h \
	KmUtils.cpp KmUtils.h
