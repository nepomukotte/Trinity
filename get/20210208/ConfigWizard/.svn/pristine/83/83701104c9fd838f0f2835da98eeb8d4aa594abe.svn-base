dnl Process this file with autoconf to produce a configure script.

AC_INIT([GetConfigWizard], [1.0.2], [patrick.sizun@cea.fr], [get-config-wizard])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR([m4])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([-Wall -Werror])
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

# ------------------------------------------------------------------------------
# Default flag reinitialization
# ------------------------------------------------------------------------------
CXXFLAGS="-Wall -Wextra ${CXXFLAGS}"
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
AC_PROG_LN_S
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

CPPFLAGS="-DBOOST_NO_CWCHAR ${CPPFLAGS}"
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
# Let users choose log level
# ------------------------------------------------------------------------------
AX_UTL_LOG

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
# Check for Utilities library
# ------------------------------------------------------------------------------
PKG_CHECK_MODULES([UTL], [Utilities${DEBUG_SUFFIX}])
AC_SUBST([UTL_CFLAGS])
AC_SUBST([UTL_LIBS])

# ------------------------------------------------------------------------------
# Check for ConfigWidgets library
# ------------------------------------------------------------------------------
PKG_CHECK_MODULES([CONFIG_WIDGETS], [ConfigWidgets${DEBUG_SUFFIX} >= 1.2.0])
AC_SUBST([CONFIG_WIDGETS_CFLAGS])
AC_SUBST([CONFIG_WIDGETS_LIBS])

# ------------------------------------------------------------------------------
# Check for Qt
# ------------------------------------------------------------------------------
AT_WITH_QT([+svg], $QT_CONFIG)
AT_REQUIRE_QT_VERSION([4.3])
if test -n "$QT_PATH"
then
	DISTCHECK_CONFIGURE_FLAGS="--with-qt=$QT_PATH $DISTCHECK_CONFIGURE_FLAGS"
fi

# ------------------------------------------------------------------------------
# Check for ROOT framework
# ------------------------------------------------------------------------------
ROOT_PATH([],[],[AC_MSG_ERROR([This package requires the ROOT framework!])])
ROOT_FEATURE([qt],
	[have_qtroot=yes
	 CPPFLAGS="$CPPFLAGS -DHAVE_QTROOT"
	 QT_CPPFLAGS="$QT_CPPFLAGS -DHAVE_QTROOT"],
	[AC_MSG_WARN([It seems ROOT was installed without the qt feature! You should re-install ROOT with the --enable-qt option.])])
AM_CONDITIONAL([HAVE_QTROOT], [test "x$have_qtroot" = "xyes"])

dnl Check whether we found ROOT5 or ROOT6
AC_MSG_CHECKING([for ROOT major version])
ROOT_VERSION_MAJOR=`expr $ROOTVERSION : '\([[0-9]]*\)'`
AC_MSG_RESULT([$ROOT_VERSION_MAJOR])
dnl If ROOT6 was found, check for and enable C++11 support.
AS_IF([test "x$ROOT_VERSION_MAJOR" = "x6"], [AX_CXX_COMPILE_STDCXX_11], [])

# ------------------------------------------------------------------------------
# Declare output files
# ------------------------------------------------------------------------------
AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
