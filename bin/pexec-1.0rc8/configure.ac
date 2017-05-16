# Initialize autoconf:
AC_INIT(pexec, 1.0rc8, apal@szofi.net)

# our default $(CFLAGS) are:
CFLAGS="-Wall -pedantic -O3 -D_GNU_SOURCE"

# check compiler:
AC_PROG_CC

# check used headers:
AC_CHECK_HEADERS(ctype.h)
AC_CHECK_HEADERS(errno.h)
AC_CHECK_HEADERS(signal.h)
AC_CHECK_HEADERS(stdarg.h)
AC_CHECK_HEADERS(time.h)
AC_CHECK_HEADERS(sys/resource.h)
AC_CHECK_HEADERS(sys/select.h)
AC_CHECK_HEADERS(sys/socket.h)
AC_CHECK_HEADERS(sys/time.h)
AC_CHECK_HEADERS(sys/un.h)
AC_CHECK_HEADERS(sys/wait.h)
AC_CHECK_HEADERS(fcntl.h)
AC_CHECK_HEADERS(netdb.h)
AC_CHECK_HEADERS(netinet/in.h)
AC_CHECK_HEADERS(sys/loadavg.h)
AC_CHECK_HEADERS(libintl.h)

# Some OS-dependent functions:
AC_CHECK_LIB(c,setenv)
AC_CHECK_LIB(c,putenv)
AC_CHECK_LIB(c,socket)
AC_CHECK_LIB(socket,socket)
AC_CHECK_LIB(c,gethostbyname)
AC_CHECK_LIB(nsl,gethostbyname)
AC_CHECK_LIB(c,getloadavg)
AC_CHECK_LIB(c,gettext)

# Debian archiver:
AC_MSG_CHECKING([for a debian archiver])
if dpkg-deb --version >/dev/null 2>&1 ; then
	AC_MSG_RESULT([dpkg-deb])
	ac_cv_debian_archiver="yes"
	DPKG_DEB=dpkg-deb
else
	AC_MSG_RESULT([no])
	ac_cv_debian_archiver="no"
	DPKG_DEB=false
fi

# help2man utility:
AC_MSG_CHECKING([help2man])
if help2man --version >/dev/null 2>&1 ; then
	AC_MSG_RESULT([help2man])
	ac_cv_help2man="yes"
	HELP2MAN=help2man
else
	AC_MSG_RESULT([no])
	ac_cv_help2man="no"
	HELP2MAN=true
fi

# Groff:
AC_MSG_CHECKING([for a roff formatter])
if groff --version >/dev/null 2>&1 ; then
	AC_MSG_RESULT([groff])
	ac_cv_roff_groff="yes"
	ROFF=groff
else
	AC_MSG_RESULT([no])
	ac_cv_roff_groff="no"
	ROFF=false
fi

# Check paper size:
AC_MSG_CHECKING([papersize])
if test -x /usr/bin/paperconf ; then 
	PAPERSIZE=`/usr/bin/paperconf`
	PAPER_WIDTH=`/usr/bin/paperconf -w`
	PAPER_HEIGHT=`/usr/bin/paperconf -h`
elif test -f /etc/papersize ; then 
	PAPERSIZE=`grep -v "^\#" /etc/papersize | head -1`
	PAPER_WIDTH=0
	PAPER_HEIGHT=0
else
	PAPERSIZE=letter
	PAPER_WIDTH=612 
	PAPER_HEIGHT=792
fi
AC_MSG_RESULT([$PAPERSIZE: ${PAPER_WIDTH}x${PAPER_HEIGHT}])

xlibs=""

if test "$ac_cv_lib_c_socket" == no && test "$ac_cv_lib_socket_socket" == no ; then
	AC_MSG_ERROR([no proper socket library is available])
elif test "$ac_cv_lib_socket_socket" == yes ; then
	AC_MSG_NOTICE([using -lsocket])
	xlibs="${xlibs} -lsocket"
fi

if test "$ac_cv_lib_c_gethostbyname" == no && test "$ac_cv_lib_nsl_gethostbyname" == no ; then
	AC_MSG_ERROR([error: no proper gethostbyname function is available])
elif test "$ac_cv_lib_c_gethostbyname" == yes ; then
	true
elif test "$ac_cv_lib_nsl_gethostbyname" == yes ; then
	AC_MSG_NOTICE([using -lnsl])
	xlibs="${xlibs} -lnsl"
fi

if test "$ac_cv_lib_c_setenv" == yes ; then
	CFLAGS="${CFLAGS} -DHAVE_SETENV"
fi

if test "$ac_cv_header_libintl_h" == yes && test "$ac_cv_lib_c_gettext" == yes ; then
	CFLAGS="${CFLAGS} -DHAVE_LIBINTL_H"
	AC_MSG_NOTICE([with internationalization support])
else
	AC_MSG_NOTICE([no support for internationalization])
fi

if test "$ac_cv_header_sys_loadavg_h" == yes ; then
	CFLAGS="${CFLAGS} -DHAVE_SYS_LOADAVG_H"
fi

AC_MSG_NOTICE([using CFLAGS: ${CFLAGS}])

DEB_ARCH=`uname -m | sed -e 's/x86_64/amd64/' -e 's/i486/i386/' -e 's/i586/i386/' -e 's/i686/i386/'`
AC_MSG_NOTICE([using Debian architecture identifier: ${DEB_ARCH}])
DEB_VERSION=$PACKAGE_VERSION

AC_SUBST(XLIBS, ${xlibs})
AC_SUBST(DEB_VERSION, ${DEB_VERSION})
AC_SUBST(DEB_ARCH, ${DEB_ARCH})
AC_SUBST(DPKG_DEB, ${DPKG_DEB})
AC_SUBST(HELP2MAN, ${HELP2MAN})
AC_SUBST(ROFF, ${ROFF})
AC_SUBST(PAPERSIZE, ${PAPERSIZE})
AC_SUBST(PAPER_WIDTH, ${PAPER_WIDTH})
AC_SUBST(PAPER_HEIGHT, ${PAPER_HEIGHT})

AC_OUTPUT([Makefile src/Makefile doc/Makefile])
