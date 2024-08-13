# ------------------------------------------------------------------------------
# Package setup
# ------------------------------------------------------------------------------

AC_INIT([MDaq], [2.1.0], [shebli.anvar@cea.fr;patrick.sizun@cea.fr])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR([m4])
AC_CONFIG_MACRO_DIR([m4])

DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(OFF)
DX_PS_FEATURE(OFF)
DX_INIT_DOXYGEN([$PACKAGE_NAME], [doc/doxyfile])

AM_INIT_AUTOMAKE([-Wall -Werror -Wno-portability-recursive])
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

INSTALL="${INSTALL} -C"

# ------------------------------------------------------------------------------
# Default flags reinitialization
# ------------------------------------------------------------------------------
CXXFLAGS="-Wall -Wextra ${CXXFLAGS}"
CFLAGS="-Wall -Wextra ${CFLAGS}"

# ------------------------------------------------------------------------------
# Check for tool support
# ------------------------------------------------------------------------------
AC_PROG_CXX
AC_PROG_LIBTOOL

# ------------------------------------------------------------------------------
# Check for C++ support
# ------------------------------------------------------------------------------
AC_LANG([C++])
AC_CHECK_HEADERS([algorithm bitset cassert cctype cmath cstdlib cstring \
   deque exception fstream functional iomanip ios iterator limits list locale \
   map memory ostream set sstream string typeinfo vector], [],
   [AC_MSG_ERROR([missing standard C++ header])]) 

# ------------------------------------------------------------------------------
# Enable silent build rules by default
# ------------------------------------------------------------------------------
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])],[AC_SUBST([AM_DEFAULT_VERBOSITY], [1])])

# ------------------------------------------------------------------------------
# Detect compilation mode: debug or release
# ------------------------------------------------------------------------------
AX_DEBUG_OR_RELEASE_BUILD
AC_SUBST([DEBUG_SUFFIX])

# ------------------------------------------------------------------------------
# Define log level
# ------------------------------------------------------------------------------
AX_UTL_LOG

# ------------------------------------------------------------------------------
# Enable VxWorks portability
# ------------------------------------------------------------------------------
AX_VXWORKS_PORTABILITY

# ------------------------------------------------------------------------------
# Enable WINDOWS portability
# ------------------------------------------------------------------------------
AX_WINDOWS_PORTABILITY

# ------------------------------------------------------------------------------
# Define variable for Apple hosts to disable features not yet ported
# ------------------------------------------------------------------------------
AM_CONDITIONAL([APPLE], [test "$host_vendor" = "apple"])

# ------------------------------------------------------------------------------
# Establish endianness
# ------------------------------------------------------------------------------
AC_C_BIGENDIAN()

# ------------------------------------------------------------------------------
# Option to disable compilation of command line client
# ------------------------------------------------------------------------------
AC_ARG_ENABLE([cli], [AS_HELP_STRING([--disable-cli],[Disable command line client])], [], [])
AS_IF([test "$host_os" = "vxworks"], [enable_cli=no])
AM_CONDITIONAL([ENABLE_CLI], [test "x$enable_cli" != "xno"], [Whether to build the command line client.])

# ------------------------------------------------------------------------------
# Check for BOOST libraries
# ------------------------------------------------------------------------------
boost_version_min_req=1.40.0
AX_BOOST_BASE([$boost_version_min_req],[],
	[AC_MSG_ERROR([Boost headers are needed to build this library.
      Please go to www.boost.org to download them.])]
)

BOOST_CPPFLAGS="-D__STDC_CONSTANT_MACROS -DBOOST_NO_CWCHAR ${BOOST_CPPFLAGS}"
if test "$host_os" = "vxworks"
then
	# Hack to get boost/iostreams 1.33 to work on VxWorks 6.3 (with wchar.hx header)
	BOOST_CPPFLAGS="-DBOOST_IOSTREAMS_NO_WIDE_STREAMS -DBOOST_NO_INTRINSIC_WCHAR_T ${BOOST_CPPFLAGS}"
fi

if test -n "$ac_boost_path"
then
	DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-boost=$ac_boost_path"
fi

# ------------------------------------------------------------------------------
# Check for specific Boost libraries
# ------------------------------------------------------------------------------
if test "$host_os" != "vxworks"
then
  AX_BOOST_SYSTEM
  AX_BOOST_THREAD
  if test "$enable_ecc" != "no"
  then
	  AX_BOOST_FILESYSTEM
	  AX_BOOST_PROGRAM_OPTIONS
  fi
fi
# ------------------------------------------------------------------------------
# Enable Ice implementation
# ------------------------------------------------------------------------------
AX_LIB_ICE_OR_ICEE

# ------------------------------------------------------------------------------
# Check for IceStorm library
# ------------------------------------------------------------------------------
AX_LIB_ICESTORM([have_icestorm=yes],
	[have_icestorm=no
	AC_MSG_WARN([Remote logging features will not be built!])])
AM_CONDITIONAL([HAVE_ICESTORM], [test "x$have_icestorm" = "xyes"])

# ------------------------------------------------------------------------------
# Check for log4cxx library
# ------------------------------------------------------------------------------
MORD_LIB_LOG4CXX([have_log4cxx=yes],
	[have_log4cxx=no])
AM_CONDITIONAL([HAVE_LOG4CXX], [test "x$have_log4cxx" = "xyes"])

# ------------------------------------------------------------------------------
# Check for readline library
# ------------------------------------------------------------------------------
# Do no let AX_LIB_READLINE pollute the LIBS variable
saved_LIBS="$LIBS"
AC_ARG_WITH([readline],
  [AS_HELP_STRING([--with-readline=DIR], [root directory of the readline library])],
  [readline_includedir="$withval/include"
   readline_libdir="$withval/lib"
  ]
)
AC_ARG_WITH([readline-includedir],
  [AS_HELP_STRING([--with-readline-includedir=DIR], [directory for the readline headers])],
  [readline_includedir="$withval"]
)
AC_ARG_WITH([readline-libdir],
  [AS_HELP_STRING([--with-readline-libdir=DIR], [directory for the readline library])],
  [readline_libdir="$withval"]
)
AS_IF([test "x$readline_includedir" != "x" ], [READLINE_CPPFLAGS="-I$readline_includedir"])
AS_IF([test "x$readline_libdir" != "x" ], [READLINE_LDFLAGS="-L$readline_libdir"])

if test "x$enable_cli" != "xno"
then
	saved_CPPFLAGS="$CPPFLAGS"
    CPPFLAGS="$CPPFLAGS $READLINE_CPPFLAGS"
    saved_LDFLAGS="$LDFLAGS"
    LDFLAGS="$LDFLAGS $READLINE_LDFLAGS"
	AX_LIB_READLINE
	AS_IF([test "$ac_cv_header_readline_history_h" = "yes"],
		[AC_SUBST([READLINE_LIB], [$ax_cv_lib_readline])],
		[AC_MSG_ERROR([Could not find GNU readline and history libraries development files. Install them or use the --with-readline* options.])]
	)
fi

# ------------------------------------------------------------------------------
# Check for ncurses library
# ------------------------------------------------------------------------------
if test "x$enable_cli" != "xno"
then
	AX_WITH_CURSES
	if test "x$ax_cv_ncursesw" != xyes && test "x$ax_cv_ncurses" != xyes; then
		AC_MSG_ERROR([requires ncurses library!])
	fi
	AC_SUBST([CURSES_LIB])
fi

# Restore initial flags before check for readline and ncurses
CPPFLAGS=$saved_CPPFLAGS"
LDFLAGS=$saved_LDFLAGS"
LIBS="$saved_LIBS"
AC_SUBST([READLINE_CPPFLAGS])
AC_SUBST([READLINE_LDFLAGS])

# ------------------------------------------------------------------------------
# Check for Utilities library
# ------------------------------------------------------------------------------
PKG_CHECK_MODULES([UTILS], [Utilities$DEBUG_SUFFIX >= 1.2.2])
AC_SUBST([UTILS_CFLAGS])
AC_SUBST([UTILS_LIBS])

# ------------------------------------------------------------------------------
# Check for CompoundConfig library
# ------------------------------------------------------------------------------
PKG_CHECK_MODULES([COMPOUND_CONFIG], [CompoundConfig$DEBUG_SUFFIX >= 1.0.0])
AC_SUBST([COMPOUND_CONFIG_CFLAGS])
AC_SUBST([COMPOUND_CONFIG_LIBS])

# ------------------------------------------------------------------------------
# Check for TinyXml library
# ------------------------------------------------------------------------------
PKG_CHECK_MODULES([TINY_XML], [tinyxml$DEBUG_SUFFIX >= 2.5.2])
AC_SUBST([TINY_XML_CFLAGS])
AC_SUBST([TINY_XML_LIBS])

# ------------------------------------------------------------------------------
# Declare output files
# ------------------------------------------------------------------------------
AC_CONFIG_FILES([
	Makefile
    mdaq-interfaces${DEBUG_SUFFIX}.pc:mdaq-interfaces.pc.in
    mdaq-mt${DEBUG_SUFFIX}.pc:mdaq-mt.pc.in
    mdaq-cli${DEBUG_SUFFIX}.pc:mdaq-cli.pc.in
    mdaq-daq${DEBUG_SUFFIX}.pc:mdaq-daq.pc.in
    mdaq-utl${DEBUG_SUFFIX}.pc:mdaq-utl.pc.in
    mdaq-hwclient${DEBUG_SUFFIX}.pc:mdaq-hwclient.pc.in
    mdaq-hwserver${DEBUG_SUFFIX}.pc:mdaq-hwserver.pc.in
    mdaq-eccserver${DEBUG_SUFFIX}.pc:mdaq-eccserver.pc.in
])
AC_OUTPUT
