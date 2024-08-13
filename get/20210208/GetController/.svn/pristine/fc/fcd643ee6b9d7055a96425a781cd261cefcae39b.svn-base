# ------------------------------------------------------------------------------
# Package setup
# ------------------------------------------------------------------------------

AC_INIT([GetController], [2.0.1], [patrick.sizun@cea.fr], [GetController])
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
DX_INIT_DOXYGEN([$PACKAGE_NAME], [doc/Doxyfile], [doc])

AM_INIT_AUTOMAKE([-Wall -Werror])
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

INSTALL="${INSTALL} -C"

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
AC_PROG_CXX
AC_PROG_LIBTOOL

# ------------------------------------------------------------------------------
# Check for C++ support
# ------------------------------------------------------------------------------
AC_LANG([C++])
AC_CHECK_HEADERS([exception iostream list memory set sstream string utility vector], [],
	[AC_MSG_ERROR([missing standard C++ header])])

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
# Choose log level
# ------------------------------------------------------------------------------
AX_UTL_LOG

# ------------------------------------------------------------------------------
# Enable Ice implementation
# ------------------------------------------------------------------------------
ACX_PTHREAD
AX_LIB_ICE_OR_ICEE

# ------------------------------------------------------------------------------
# Check for IceStorm library
# ------------------------------------------------------------------------------
AX_LIB_ICESTORM([have_icestorm=yes],
	[have_icestorm=no
	AC_MSG_WARN([The ECC logger cannot be subscribed to without IceStorm!])])
AM_CONDITIONAL([HAVE_ICESTORM], [test "x$have_icestorm" = "xyes"])

# ------------------------------------------------------------------------------
# Check for Qt
# ------------------------------------------------------------------------------
AT_WITH_QT([], [$QT_CONFIG], [QT += widgets xml network svg printsupport])
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
	AC_MSG_ERROR([Could not find Qwt library. Use --with-qwt5 option.])
fi
# Check version
if test "$QWT_VERSION_MAJOR" -lt 6 ; then
	AC_MSG_WARN([Found version $QWT_VERSION_STR of the Qwt library. It might work but a more recent version would be better.])
fi

# ------------------------------------------------------------------------------
# Check for CompoundConfig library
# ------------------------------------------------------------------------------
PKG_CHECK_MODULES([COMPOUND_CONFIG], [CompoundConfig$DEBUG_SUFFIX >= 1.0.0])
AC_SUBST([COMPOUND_CONFIG_CFLAGS])
AC_SUBST([COMPOUND_CONFIG_LIBS])

# ------------------------------------------------------------------------------
# Check for CompoundConfig library
# ------------------------------------------------------------------------------
CCFG_LIB_NAME="CompoundConfig${DEBUG_SUFFIX}"
AC_SUBST([CCFG_LIB_NAME])

PKG_CHECK_MODULES(COMPOUND_CONFIG, $CCFG_LIB_NAME >= 2.2.0)
AC_SUBST(COMPOUND_CONFIG_CFLAGS)
AC_SUBST(COMPOUND_CONFIG_LIBS)

# ------------------------------------------------------------------------------
# Check for ConfigWidgets library
# ------------------------------------------------------------------------------
CONFIG_WIDGETS_LIB_NAME="ConfigWidgets${DEBUG_SUFFIX}"
AC_SUBST([CONFIG_WIDGETS_LIB_NAME])

PKG_CHECK_MODULES([CONFIG_WIDGETS], [$CONFIG_WIDGETS_LIB_NAME >= 1.2.0])
AC_SUBST([CONFIG_WIDGETS_CFLAGS])
AC_SUBST([CONFIG_WIDGETS_LIBS])

# ------------------------------------------------------------------------------
# Check for BOOST libraries
# ------------------------------------------------------------------------------
AX_BOOST_BASE([1.40.0],[],
	[AC_MSG_ERROR([Boost headers are needed to build this library.
      Please go to www.boost.org to download them.])]
)
BOOST_CPPFLAGS="-D__STDC_CONSTANT_MACROS -DBOOST_NO_CWCHAR ${BOOST_CPPFLAGS}"
if test -n "$ac_boost_path"
then
	DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-boost=$ac_boost_path"
fi

AX_BOOST_PROGRAM_OPTIONS

if test "$enable_threadsafe_logging" = "yes"
then
	AX_BOOST_SYSTEM
	AX_BOOST_THREAD
fi

# ------------------------------------------------------------------------------
# Check for GetBench library
# ------------------------------------------------------------------------------

PKG_CHECK_MODULES([GET_BENCH], [getbench-daq${DEBUG_SUFFIX} >= 2.2.1 getbench-cfg${DEBUG_SUFFIX} >= 2.2.1 getbench-eccclient${DEBUG_SUFFIX} >= 2.2.1])
AC_SUBST([GET_BENCH_CFLAGS])
AC_SUBST([GET_BENCH_LIBS])

# ------------------------------------------------------------------------------
# Check for MultiFrame library
# ------------------------------------------------------------------------------
if test "$host_os" != "vxworks"
then
	PKG_CHECK_MODULES([MFM], [MultiFrame${DEBUG_SUFFIX} >= 2.6.0])
	AC_SUBST([MFM_CFLAGS])
	AC_SUBST([MFM_LIBS])
fi

# ------------------------------------------------------------------------------
# Declare output files
# ------------------------------------------------------------------------------
AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
