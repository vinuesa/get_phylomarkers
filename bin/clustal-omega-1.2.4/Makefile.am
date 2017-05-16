# Some useful automake documentation:
# http://www.openismus.com/documents/linux/automake/automake.shtml
# http://www.bioinf.uni-freiburg.de/~mmann/HowTo/automake.html

# we do not use GNU-style files: NEWS README AUTHORS ChangeLog
#AUTOMAKE_OPTIONS = foreign

ACLOCAL_AMFLAGS=-I m4
SUBDIRS = src
# SUBDIRS = src tests

dist_doc_DATA = README

EXTRA_DIST = Doxyfile

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = clustalo.pc


bug-tests: all
	cd bug-tests && $(SHELL) bug-tests.sh;

#dist-hook:
#	svnversion . > $(distdir)/src/svnversion.h 2>/dev/null

