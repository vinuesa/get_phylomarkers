# configure.ac for Clustal Omega
#
# RCS $Id: configure.ac 317 2016-12-16 17:44:50Z fabian $


# release
#
#AC_INIT([Clustal Omega], [0.1.1], [clustalw@ucd.ie])
#PACKAGE_CODENAME="OppidumLeonis"
#r238
#AC_INIT([Clustal Omega], [0.2.0], [clustalw@ucd.ie])
#PACKAGE_CODENAME="IbiMeVidebunt"
#r250
#AC_INIT([Clustal Omega], [1.0.1], [clustalw@ucd.ie])
#PACKAGE_CODENAME="DeshilHollesEamus"
#r255
#AC_INIT([Clustal Omega], [1.0.2], [clustalw@ucd.ie])
#PACKAGE_CODENAME="DeshilHollesEamus"
#AC_INIT([Clustal Omega], [1.0.3], [clustalw@ucd.ie])
#PACKAGE_CODENAME="DeshilHollesEamus"
#AC_INIT([Clustal Omega], [1.0.4], [clustalw@ucd.ie])
#PACKAGE_CODENAME="DeshilHollesEamus"
#AC_INIT([Clustal Omega], [1.1.0], [clustalw@ucd.ie])
#PACKAGE_CODENAME="FilumVitae"
#AC_INIT([Clustal Omega], [1.1.1], [clustalw@ucd.ie])
#PACKAGE_CODENAME="FilumVitae"
#AC_INIT([Clustal Omega], [1.2.0], [clustalw@ucd.ie])
#PACKAGE_CODENAME="AndreaGiacomo"
#AC_INIT([Clustal Omega], [1.2.1], [clustalw@ucd.ie])
#PACKAGE_CODENAME="AndreaGiacomo"
# r310
#AC_INIT([Clustal Omega], [1.2.2], [clustalw@ucd.ie])
#PACKAGE_CODENAME="AndreaGiacomo"
# r313
#AC_INIT([Clustal Omega], [1.2.3], [clustalw@ucd.ie])
#PACKAGE_CODENAME="AndreaGiacomo"
# r317
AC_INIT([Clustal Omega], [1.2.4], [clustalw@ucd.ie])
PACKAGE_CODENAME="AndreaGiacomo"

# The AC_INIT macro can take any source file as an argument. It just
# checks that the file is there, which should, in turn, mean that the
# source directory is there. AC_INIT generates the following m4
# macros, output variables and preprocessor symbols:
AC_SUBST(PACKAGE_NAME)
AC_SUBST(PACKAGE_TARNAME)
AC_SUBST(PACKAGE_VERSION)
AC_SUBST(PACKAGE_STRING)
AC_SUBST(PACKAGE_BUGREPORT)
AC_SUBST(PACKAGE_CODENAME)

AC_MSG_NOTICE([

    Welcome to $PACKAGE_NAME - version $PACKAGE_VERSION ($PACKAGE_CODENAME)
                                               
                       +NMMMMMMMMMS=             
                    MMMMM?      :MMMMM8          
                  IMMMMS           MMMMM,        
                ~MMMMN              ,MMMMM       
               =MMMMM                 MMMMM      
              ?MMMMM                  +MMMMM     
              MMMMMS                   MMMMM     
              MMMMM                    MMMMM:    
              MMMMM                    MMMMM+    
              MMMMM                    MMMMM~    
              MMMMM                    MMMMM     
              MMMMMZ                   MMMMM     
               MMMMM                  ,MMMMZ     
                MMMM8                 MMMM+      
                 MMMMS               MMMM,       
                   MMMMO           MMMM          
            ,N       SMMM        MMM8       +8   
            ,M7         M+       M+        :M8   
            ,MMMMMMMMMMMMI       MMMMMMMMMMMM8   
            ,MMMMMMMMMMMM8       MMMMMMMMMMMM8   
            ,MMMMMMMMMMMMN       MMMMMMMMMMMM8   
                                               
              ])
#
# Created with http://www.glassgiant.com/ascii/ and Wikipedia's Omega
# Replaced all dollars with S's.
#
# 
#                 ,MMMMMMMMMD:        MMMMM     
#               NMMMMMMMMMMMMMM+      SMMMM     
#             DMMMMMMMMMMMMMMMMMM,    +MMMM     
#            MMMMMMMMM8ZSO8MMMMMMMM   +MMMM     
#           MMMMM7             SMMMM  +MMMM     
#          MMMO                   MM? +MMMM     
#         ?MMI                     NMIIMMMM     
#         MMM                       MMMMMMM     
#         MM                                    
#         MN                                    
#         MD                                    
#         MM                                    
#         MM?                       ,~+ISON     
#         =MMI                     MMSOMMMM     
#          MMMZ                   MMS +MMMM     
#           MMMM8:             :NMMM  +MMMM     
#            MMMMMMMMMDOODMMMMMMMMM   +MMMM     
#             MMMMMMMMMMMMMMMMMMMM    +MMMM     
#               OMMMMMMMMMMMMMMZ      SMMMM     
#                  MMMMMMMMMM+        MMMMM     
#                                                


# AC_MSG_NOTICE([Configuring for your system...])

AC_PREREQ(2.6) # 2.60 needed for AC_USE_SYSTEM_EXTENSIONS
# consider 2.64 which allows to define a URL as well

# The AM_INIT_AUTOMAKE line adds several standard checks. 
# This macro is always needed for automake
# Obsolete: It takes the program name and version number as arguments.
AM_INIT_AUTOMAKE

# config header: provides a way for people to customise the
# configuration for their platform, via #defines. config.h.in can be
# generated automatically with the autoheader tool. However, you need
# a stamp-h file in your project to ensure that automake regenerates
# config.h from config.h.in. Type 'touch stamp-h' to add this file to
# your project.
#
AC_CONFIG_HEADERS([src/config.h]) 

# the following allows to install a copy config.h for libclustalo in a
# namespace like way, i.e. no name clash with your own config.h
AX_PREFIX_CONFIG_H(src/clustal-omega-config.h)

# This macro was introduced in Autoconf 2.60. If possible, enable
# extensions to C or Posix on hosts that normally disable the
# extensions.
AC_USE_SYSTEM_EXTENSIONS


# for cross compilation
#
AC_CANONICAL_HOST
## AC_CANONICAL_BUILD
## AC_CANONICAL_TARGET

# Set the @SET_MAKE@ variable=make if $(MAKE) not set
AC_PROG_MAKE_SET

# AC_PROG_RANLIB will allow you to build code in sub-directories into
# temporary libraries, which make will then link in with the rest of
# the code.
# Rendered obsolte by LT_INIT: AC_PROG_RANLIB
#AC_PROG_LIBTOOL
LT_INIT

# AC_PROG_INSTALL will generate an install target so that users may
# just type 'make install' to install the software.
AC_PROG_INSTALL

# check for C++ preprocessor and compiler and the library compiler
# (might change the compiler flags)
#
# Determine a C compiler to use. If using the GNU C compiler, set
# shell variable GCC to `yes'. If output variable CFLAGS was not
# already set, set it to `-g -O2' for the GNU C compiler (`-O2' on
# systems where GCC does not accept `-g'), or `-g' for other
# compilers.
AC_PROG_CC
# Set output variable CPP to a command that runs the C preprocessor. 
AC_PROG_CPP
#
# Determine a C++ compiler to use. If output variable CXXFLAGS was not
# already set, set it to `-g -O2' for the GNU C++ compiler (`-O2' on
# systems where G++ does not accept `-g'), or `-g' for other
# compilers.
AC_PROG_CXX
# Set CXXCPP to a command that runs the C++ preprocessor
AC_PROG_CXXCPP
#
# reset compiler flags to initial flags
#CXXFLAGS=$OLD_CXXFLAGS




# safety check: AC_CONFIG_SRCDIR([src/clustal.h])

# Let source files know in which context they are being used
AC_DEFINE([CLUSTALO], [1],
          [Let all files know they are being used inside Clustal Omega])

AC_DEFINE_UNQUOTED(PACKAGE_CODENAME, "$PACKAGE_CODENAME", [The package code name])


# Get svn revision
#AC_PATH_PROG(svnversioncommand, svnversion)
#if test "X$svnversioncommand" = "X" || test `$svnversioncommand -n '.'` = "exported"; then
#AC_DEFINE_UNQUOTED(SVN_VERSION, ["`cat saved_svn_version.txt`"], [repository svn version])
#else
#AC_DEFINE_UNQUOTED(SVN_VERSION, ["`svnversion -n`"], [repository svn version])
#`svnversion -n > saved_svn_version.txt`
#fi




# AW: we might need this in the future
#
case "$host_os" in
linux*)
#      AC_DEFINE(LINUX,1,[This is a Linux system])
     AC_DEFINE_UNQUOTED(LINUX,1,[This is a LINUX system])
#     AM_CONDITIONAL(HAVE_GNU_LD, true)    
#     AM_CONDITIONAL(SOMEBSD, false)
#     XTRA_CPPLIBS=-lstdc++
#     LIBEXT=.so
     ;;
freebsd*)
     AC_DEFINE_UNQUOTED(SOMEBSD,1,[This is a BSD system])
#     AC_CHECK_LIB(c_r, pthread_create)
#     AM_CONDITIONAL(HAVE_GNU_LD, true)
#     AM_CONDITIONAL(SOMEBSD, true)
#     XTRA_CPPLIBS=-lstdc++
#     LIBEXT=.so
     ;;
openbsd*)
     AC_DEFINE_UNQUOTED(SOMEBSD,1,[This is a BSD system])
#     AC_CHECK_LIB(c_r, pthread_create)
#     AM_CONDITIONAL(HAVE_GNU_LD, true)
#     AM_CONDITIONAL(SOMEBSD, true)
#     XTRA_CPPLIBS=-lstdc++
#     LIBEXT=.so
     ;;
netbsd*)
     AC_DEFINE_UNQUOTED(SOMEBSD,1,[This is a BSD system])
#     AC_CHECK_LIB(c_r, pthread_create)
#     AM_CONDITIONAL(HAVE_GNU_LD, true)
#     AM_CONDITIONAL(SOMEBSD, true)
#     XTRA_CPPLIBS=-lstdc++
#     LIBEXT=.so
     ;;
*solaris*)
     AC_DEFINE_UNQUOTED(SOLARIS,1,[This is a Solaris system])
#     AC_CHECK_LIB(resolv, res_init)
#     XTRA_CPPLIBS=-lstdc++
#     AM_CONDITIONAL(HAVE_GNU_LD, false)
#     AM_CONDITIONAL(SOMEBSD, false)
#     CFLAGS="-D_POSIX_PTHREAD_SEMANTICS $CFLAGS"
#     LIBEXT=.so
     ;;
darwin*)
     AC_DEFINE_UNQUOTED(DARWIN,1,[This is a Darwin system])
#     AM_CONDITIONAL(HAVE_GNU_LD, false)
#     AM_CONDITIONAL(SOMEBSD, false)
#     CFLAGS="-fno-common $CFLAGS"
#     LIBEXT=.so
     ;;
cygwin*)
     AC_DEFINE_UNQUOTED(CYGWIN,1,[This is a CYGWIN system])
#     LDFLAGS="$LDFLAGS -no-undefined"
#     AM_CONDITIONAL(SOMEBSD, false)
#     AM_CONDITIONAL(HAVE_GNU_LD, false)
#     LIBEXT=.dll
     ;;
mingw*)
     AC_DEFINE_UNQUOTED(MINGW,1,[This is a MinGW system])
#     AC_DEFINE_UNQUOTED(WINDOWS,1,[This is a Windows system])
#     AC_CHECK_LIB(intl, gettext)
#     AC_CHECK_LIB(ws2_32, _head_libws2_32_a,,[AC_MSG_ERROR([libextractor requires Winsock2])])
#     AC_CHECK_LIB(plibc, plibc_init,,[AC_MSG_ERROR([libextractor requires PlibC])])
#     LDFLAGS="$LDFLAGS -no-undefined -Wl,--export-all-symbols"
#		 LIBSOLD=$LIBS
#     AM_CONDITIONAL(HAVE_GNU_LD, true)
#     AM_CONDITIONAL(SOMEBSD, false)
#     LIBEXT=.dll
     ;;
*)
#     AC_MSG_RESULT(Unrecognised OS $host_os)
     AC_DEFINE_UNQUOTED(OTHEROS,1,[Some strange OS])
#     AC_MSG_RESULT(otheros)
#     AM_CONDITIONAL(HAVE_GNU_LD, false)
#     AM_CONDITIONAL(SOMEBSD, false)
#     LIBEXT=.so
     ;;
esac



# Checks for standard header files.
# obsolescent according to http://www.gnu.org/software/hello/manual/autoconf/Particular-Headers.html
AC_HEADER_STDC
# Guess the next one is obsolete as well then?
AC_CHECK_HEADERS(assert.h ctype.h float.h limits.h math.h memory.h stdarg.h stdio.h stdlib.h string.h time.h unistd.h)
#  awk  '/^#include *</ {print $2}' $(find src/ -name \*.[ch] -or -name \*.cpp) | grep -v argtable2.h  | sort -u | tr -d '<>'

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T
AC_C_INLINE

# Checks for library functions.
# AC_CHECK_FUNCS(strdup strstr strchr erand48)
# AW is there an automatic way to extract library functions used in the program?


# explicit libm check
AC_CHECK_LIB(m, log,, AC_MSG_ERROR([Could not find libm]))
AC_CHECK_LIB(m, log2, AC_DEFINE(HAVE_LOG2, 1, [Has log2()]))


# ---    argtable header and library check
#
#
AC_CHECK_HEADERS([argtable2.h], [],
                 AC_MSG_ERROR([Could not find argtable2.h. Try $ ./configure CFLAGS='-Iyour-argtable2-include-path]))
AC_CHECK_LIB(argtable2, arg_parse, [],
             AC_MSG_ERROR([Could not find libargtable2. Try $ ./configure LDFLAGS="-Lyour-argtable2-lib-path']))


# enable python if requested (--with-python)
#
#AC_ARG_WITH([python],
#        [AS_HELP_STRING([--with-python],
#        [Enable Python Interface])],
#        [with_python=$withval],
#        [with_python=no])
#if test "$with_python" = "yes"; then
#   AC_MSG_CHECKING([for python-config to get compiler flags])
#   AC_CHECK_PROG(python_config, 'python-config', yes)
#   if test "${python_config}" != "yes"; then
#      AC_MSG_ERROR("Could not run python-config.")
#   fi
#   PYTHON_CFLAGS=`python-config --cflags`
#   PYTHON_LDFLAGS=`python-config --ldflags`
#   AC_DEFINE([PYTHON_ENABLED], [1], [define for python usage])
#fi
#AC_SUBST(PYTHON_CFLAGS)
#AC_SUBST(PYTHON_LDFLAGS)


# ---   SSE support
# "One instruction to do the same operation on 4 packed elements simultaneously."
#
# http://www.codeproject.com/KB/recipes/sseintro.aspx?msg=568491
# http://www.drdobbs.com/go-parallel/blog/archives/2010/03/getting_amazing.html
# http://software.intel.com/file/1000
#
# AC_MSG_CHECKING(for SSE in current arch/CFLAGS)
# AC_LINK_IFELSE([
# AC_LANG_PROGRAM([[
# #include <xmmintrin.h>
# __m128 testfunc(float *a, float *b) {
#   return _mm_add_ps(_mm_loadu_ps(a), _mm_loadu_ps(b));
# }
# ]])],
# [has_sse=yes], [has_sse=no]
# )
# AC_MSG_RESULT($has_sse)  
#  
# AC_ARG_ENABLE([sse],
#         [AS_HELP_STRING([--enable-sse],
#                         [Enable SSE support])],
#         [has_sse="$enableval"],
#         [has_sse=no])
#  
# if test x"$has_sse" = x"yes"; then
#    AC_DEFINE([USE_SSE], , [Enable SSE support])
#    CFLAGS="${CFLAGS} -msse"
# fi


# --- openmp
#
AC_CONFIG_MACRO_DIR([m4])
m4_include([m4/ax_openmp.m4])
AC_ARG_WITH([openmp],
    AS_HELP_STRING([--without-openmp], [Do not use OpenMP]))
AS_IF([test "x$with_openmp" = "xno"], 
            [], 
            [AX_OPENMP(
                AC_DEFINE([HAVE_OPENMP],[1],[Defined if OpenMP should and can be used])
             )]
     )
AC_SUBST(OPENMP_CFLAGS)
# according to http://beagle-lib.googlecode.com/svn-history/r824/trunk/configure.ac:
# "OpenMP checker only defines for C when compiling both C and C++"
OPENMP_CXXFLAGS=$OPENMP_CFLAGS
AC_SUBST(OPENMP_CXXFLAGS)

# version 2: enable openmp only if explicitely asked for it
#AC_ARG_ENABLE(openmp, AS_HELP_STRING([--enable-openmp], [Multi-threaded parallel version using OpenMP]))
#if test -n "${enable_openmp}" -a x"${enable_openmp}" != x"no"; then
#   m4_include([m4/ax_openmp.m4])
#   AX_OPENMP(AC_DEFINE(HAVE_OPENMP, 1, [enable OpenMP support]),
#             AC_MSG_ERROR([Compiler does not seem to support OpenMP]))
#   AC_SUBST(OPENMP_CFLAGS)
#   AC_SUBST(OPENMP_CXXFLAGS)
#fi

# ---   squid
#
#
# Extra checks for squid; see their configure, config.log, squidconf.h
# etc. For example, I'm not sure about USE_HOST_BYTESWAP_FUNCTIONS,
# SIZEOF_UNSIGNED_LONG, HAVE_STRTOUL and HAVE_STRTOULL
#
# Might not just have to set compiler flags but create use squidconf.h.in
#
# The following hopefully takes about byte ordering (big-endian...).
# My hope is that it will set WORDS_BIGENDIAN and
# USE_HOST_BYTESWAP_FUNCTIONS (as used by squid) if necessary:
#
AC_C_BIGENDIAN(AC_MSG_WARN([$PACKAGE_NAME hasn't been tested on Big-Endian machines...expecting trouble!]))
AC_CHECK_FUNCS([ntohs ntohl ntons ntonl])
#
AC_CHECK_SIZEOF(unsigned short)
AC_CHECK_SIZEOF(unsigned int)
AC_CHECK_SIZEOF(unsigned long)
AC_CHECK_SIZEOF(unsigned long long)
#
AC_CHECK_FUNCS([stroul strtoull])
AC_CHECK_FUNCS([ftello fseeko])
AC_CHECK_FUNCS([ftello64 fseeko64])
AC_CHECK_FUNCS([ftell64 fseek64])
#
AC_TYPE_OFF_T
AC_CHECK_SIZEOF(fpos_t)
#
# struct stat64 and stat64 test:
# version 1:
#AC_CHECK_TYPE(struct stat64,
#               [AC_CHECK_FUNC([stat64],
#                             AC_DEFINE([HAVE_STAT64], [1], [Define if you have the stat64 function]))])
# version 2:
#AC_MSG_CHECKING(for struct stat64)
#AC_CACHE_VAL(sc_cv_type_stat64,
#[AC_TRY_COMPILE([#include <sys/stat.h>],[struct stat64 s;],
#[sc_cv_type_stat64=yes],
#[sc_cv_type_stat64=no])])
#if test $sc_cv_type_stat64 = yes; then
#   AC_MSG_RESULT($sc_cv_type_stat64)
#   AC_CHECK_FUNC([stat64], AC_DEFINE([HAVE_STAT64], [1], [Define if you have the stat64 function]))
#   #AC_DEFINE(HAVE_STAT64, [1], [Define if you have the stat64 function])
#fi
#
#
# Above fails on myosin because storage size of struct 64 is unknown:
#  ssi.c:480: error: storage size of ‘s1’ isn’t known 
#
# On my MacBook OS X 10.6 I get:
#  ssi.c:481: warning: ‘stat64’ is deprecated (declared at /usr/include/sys/stat.h:465)
#
# So forget about it: we won't use stat64 or define HAVE_STAT64 and
# consequently do no tests here
#



# Add support for the Dmalloc package. If the user runs configure with
# --with-dmalloc, then define WITH_DMALLOC and add -ldmalloc to LIBS.
AM_WITH_DMALLOC


# Check if svnversion is installed
#AC_PATH_PROG([SVNVERSION_BIN], [svnversionFIXME], [no])
## checking if we build in a subversion WC
#wc_found="no";
#if test "$SVNVERSION_BIN" != "no"; then
#   SVNVERSION=`$SVNVERSION_BIN -n $srcdir`
#   if test $SVNVERSION != "exported";
#   then
#       AC_MSG_NOTICE($SVNVERSION);
#   fi
#fi


# ---   debugging
# 
# 
#
AC_MSG_CHECKING([whether to build with debug information])
AC_ARG_ENABLE([debug],
              [AS_HELP_STRING([--enable-debug],
                              [Enable debugging (default=no)])],
              [debugit="$enableval"],
              [debugit=no])
AC_MSG_RESULT([$debugit])

if test x"$debugit" = x"yes"; then
    # User variables CFLAGS/CXXFLAGS are used last during compilation and
    # can therefore overwrite system settings.
    #
    # if debugging was requested, add debugging info at the end of
    # C[XX]FLAGS, which makes sure previous optimisation flags are
    # overwritten (normally C[XX]FLAGS should not be modified as they
    # are user variables):
    # switch off assert etc
    # AC_DEFINE([DEBUG],[],[Debug Mode])
    CFLAGS="${CFLAGS} -g -ansi -Wall -O0"
    CXXFLAGS="${CXXFLAGS} -g -ansi -Wall -O0"
else
    # switch off assert etc
    AC_DEFINE([NDEBUG],[],[No-debug Mode])
    # if no debugging is requested, we make optimisation explicit
    AM_CXXFLAGS="${AM_CXXFLAGS} -O3"
    AM_CFLAGS="${AM_CFLAGS} -O3"
fi



AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_CFLAGS])
AC_SUBST([AM_LDFLAGS])


AC_CONFIG_FILES(Makefile
 src/Makefile
 clustalo.pc
 Doxyfile
 src/clustal/Makefile
 src/hhalign/Makefile
 src/kmpp/Makefile
 src/squid/Makefile)
AC_OUTPUT
