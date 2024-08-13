dnl Process this file with autoconf to produce a configure script.

AC_INIT([CoBoFrameViewer], [1.1.0], [patrick.sizun@cea.fr], [CoBoFrameViewer])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR([m4])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([-Wall -Werror -Wno-portability-recursive])
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

INSTALL="${INSTALL} -C"

# ------------------------------------------------------------------------------
# Default flag reinitialization
# ------------------------------------------------------------------------------
CXXFLAGS="-Wall -Wno-deprecated-declarations -Wextra ${CXXFLAGS}"
CFLAGS="-Wall -Wextra ${CFLAGS}"

# ------------------------------------------------------------------------------
# Host setup
# ------------------------------------------------------------------------------
if test -n "$host"
then
	DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --host=$host_alias"
fi
case "$host_os" in
mingw*)
	CXXFLAGS="-Wno-strict-aliasing $CXXFLAGS"
esac

# ------------------------------------------------------------------------------
# Check for tool support
# ------------------------------------------------------------------------------
AC_PROG_CXX
AC_PROG_LIBTOOL

# ------------------------------------------------------------------------------
# Check for C++ support
# ------------------------------------------------------------------------------
AC_LANG([C++])
AC_CHECK_HEADERS([exception iostream list memory set sstream string utility vector], [],
	[AC_MSG_ERROR([missing standard C++ header])])

# ------------------------------------------------------------------------------
# Check for BOOST libraries
# ------------------------------------------------------------------------------
# Version 3 of the Boost  Filesystem library is only available starting with version 1.44 of Boost.
AX_BOOST_BASE([],[],
	[AC_MSG_ERROR([Boost headers are needed to build this library.
      Please go to www.boost.org to download them.])]
)

BOOST_CPPFLAGS="-D__STDC_CONSTANT_MACROS -DBOOST_NO_CWCHAR ${BOOST_CPPFLAGS}"
if test -n "$ac_boost_path"
then
	DISTCHECK_CONFIGURE_FLAGS="--with-boost=$ac_boost_path $DISTCHECK_CONFIGURE_FLAGS"
fi

AX_BOOST_PROGRAM_OPTIONS
AX_BOOST_FILESYSTEM
AX_BOOST_SYSTEM

# ------------------------------------------------------------------------------
# Enable silent build rules by default
# ------------------------------------------------------------------------------
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])],[AC_SUBST([AM_DEFAULT_VERBOSITY], [0])])

# ------------------------------------------------------------------------------
# Enable WINDOWS portability
# ------------------------------------------------------------------------------
AX_WINDOWS_PORTABILITY

# ------------------------------------------------------------------------------
# Detect compilation mode: debug or release
# ------------------------------------------------------------------------------
AX_DEBUG_OR_RELEASE_BUILD
if test "$enable_debug" = "yes"
then
	QT_CONFIG="+debug -debug_and_release"
else
	QT_CONFIG="+release -debug_and_release"
fi
AC_SUBST([DEBUG_SUFFIX])

# ------------------------------------------------------------------------------
# Define log level
# ------------------------------------------------------------------------------
AX_UTL_LOG

# ------------------------------------------------------------------------------
# GProf profiler
# ------------------------------------------------------------------------------
AC_ARG_ENABLE([gprof],
	[AS_HELP_STRING([--enable-gprof], [Enable GProf GNU profiler])],
	[],
	[])
AS_IF([test "x$enable_gprof" = "xyes"], [
	AC_MSG_NOTICE([Enabling gprof profiling data (to gmon.out).])
	CFLAGS="$CFLAGS -pg"
	LDFLAGS="$LDFLAGS -pg"
])

# ------------------------------------------------------------------------------
# Construct library name
# ------------------------------------------------------------------------------
MFM_LIB_NAME="MultiFrame"
AC_SUBST([MFM_LIB_NAME])

# ------------------------------------------------------------------------------
# Check for Qt
# ------------------------------------------------------------------------------
AT_WITH_QT([], [$QT_CONFIG], [QT += opengl widgets svg printsupport])

AT_REQUIRE_QT_VERSION([4.3])
if test -n "$QT_PATH"
then
	DISTCHECK_CONFIGURE_FLAGS="--with-qt=$QT_PATH $DISTCHECK_CONFIGURE_FLAGS"
fi

# ------------------------------------------------------------------------------
# Check for Qwt
# ------------------------------------------------------------------------------
AX_LIB_QWT
if test "$have_qwt" = "no"
then
	AC_MSG_ERROR([Could not find the Qwt library. Use the --with-qwt5 option.])
fi
# Check version
if test "$QWT_VERSION_MAJOR" -lt 6 ; then
	AC_MSG_WARN([Found version $QWT_VERSION_STR of the Qwt library. While this should be Ok, you may benefit from upgrading to a more recent version.])
fi
# Hack to remove strip pedantic flag from CXXFLAGS since Qwt has many extra ';' errors
CXXFLAGS=`echo $CXXFLAGS | sed 's,-pedantic-errors,,' | sed 's,-pedantic,,'`

# ------------------------------------------------------------------------------
# Check for QwtPlot3d
# ------------------------------------------------------------------------------
AX_LIB_QWTPLOT3D([have_qwt_plot_3d=yes],
	[have_qwt_plot_3d=no
	AC_MSG_WARN([To enable the 3D mode of the frame viewer, you need to install the QwtPlot3d library!])])
AM_CONDITIONAL([BUILD_3DVIEWER], [test "x$have_qwt_plot_3d" = "xyes"])

# ------------------------------------------------------------------------------
# Check for MFM library
# ------------------------------------------------------------------------------
PKG_CHECK_MODULES([MFM], [$MFM_LIB_NAME >= 2.8.1])
AC_SUBST([MFM_CFLAGS])
AC_SUBST([MFM_LIBS])

# ------------------------------------------------------------------------------
# Check for specific version of Utilities library
# ------------------------------------------------------------------------------
PKG_CHECK_MODULES([UTILS], [Utilities >= 1.2.1])

# ------------------------------------------------------------------------------
# Check for ROOT framework
# ------------------------------------------------------------------------------
AC_ARG_ENABLE([root],
	[AS_HELP_STRING([--disable-root], [Disable ROOT-related features])],
	[],
	[enable_root=yes])

AS_IF([test "x$enable_root" = "xyes"], [
dnl Check for ROOT
	ROOT_PATH([], [have_root=yes], [have_root=no])
dnl Check whether we found ROOT5 or ROOT6
	AS_IF([test "x$have_root" = "xyes"], [
		AC_MSG_CHECKING([for ROOT major version])
		ROOT_VERSION_MAJOR=`expr $ROOTVERSION : '\([[0-9]]*\)'`
		AC_MSG_RESULT([$ROOT_VERSION_MAJOR])
dnl If ROOT6 was found, check for and enable C++11 support.
		AS_IF([test "x$ROOT_VERSION_MAJOR" = "x6"], [AX_CXX_COMPILE_STDCXX_11], [])
	], [])
], [])

AS_IF([test "x$have_root" = "xyes"], [AC_DEFINE([HAVE_ROOT], [1], [Whether the ROOT framework has been found])])
AS_IF([test "x$have_root" = "xno"], [AC_MSG_WARN([To build the .graw to ROOT converter, you need to install the ROOT framework!])])
AM_CONDITIONAL([BUILD_GRAW2ROOT], [test "x$have_root" = "xyes"])
AM_CONDITIONAL([HAVE_ROOT6], [test "x$ROOT_VERSION_MAJOR" = "x6"])

# ------------------------------------------------------------------------------
# Check for CCfits library
# ------------------------------------------------------------------------------
PKG_CHECK_MODULES([CCFITS], [cfitsio CCfits >= 2.5],
	[AC_DEFINE([HAVE_CCFITS], [1], [Defined if the CCfits library is available])
	 have_ccfits=yes],
	[AC_MSG_WARN([To install the raw to FITS data converter, install the CCfits library])])
AM_CONDITIONAL([HAVE_CCFITS], [test "x$have_ccfits" = "xyes"]) 

# ------------------------------------------------------------------------------
# Declare output files
# ------------------------------------------------------------------------------
AC_CONFIG_FILES([Makefile
    src/Makefile
    src/get/Makefile
    src/graw2gmat/Makefile
    src/root/Makefile
    src/graw2root/Makefile
    src/root2disp/Makefile
    src/merge/Makefile
    src/graw2fits/Makefile
    data/Makefile
])
AC_OUTPUT
