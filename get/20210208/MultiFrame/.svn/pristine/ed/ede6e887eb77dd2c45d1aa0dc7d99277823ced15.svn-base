# ------------------------------------------------------------------------------
# Package setup
# ------------------------------------------------------------------------------

AC_INIT([MultiFrame], [2.8.1], [shebli.anvar@cea.fr;patrick.sizun@cea.fr])
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

AM_INIT_AUTOMAKE([-Wall -Werror])
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

# ------------------------------------------------------------------------------
# Enable silent build rules by default
# ------------------------------------------------------------------------------
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])],[AC_SUBST([AM_DEFAULT_VERBOSITY], [1])])

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
# Detect compilation mode: debug or release
# ------------------------------------------------------------------------------
AX_DEBUG_OR_RELEASE_BUILD
AC_SUBST([DEBUG_SUFFIX])

# ------------------------------------------------------------------------------
# Enable VxWorks portability
# ------------------------------------------------------------------------------
AX_VXWORKS_PORTABILITY

# ------------------------------------------------------------------------------
# Enable WINDOWS portability
# ------------------------------------------------------------------------------
AX_WINDOWS_PORTABILITY

# ------------------------------------------------------------------------------
# Define log level
# ------------------------------------------------------------------------------
AX_UTL_LOG([info])

# ------------------------------------------------------------------------------
# Establish endianness
# ------------------------------------------------------------------------------
AC_C_BIGENDIAN()

# ------------------------------------------------------------------------------
# Check for BOOST headers
boost_version_min_req=1.33.1
AX_BOOST_BASE([$boost_version_min_req], [], [AC_MSG_ERROR([Boost headers are needed to build this library.])])

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
# Check for CompoundConfig library
CCFG_LIB_NAME="CompoundConfig${DEBUG_SUFFIX}"
PKG_CHECK_MODULES(CCFG, $CCFG_LIB_NAME >= 2.1.0)
AC_SUBST([CCFG_CFLAGS])
AC_SUBST([CCFG_LIBS])

# ------------------------------------------------------------------------------
# Check for Utilities library
# ------------------------------------------------------------------------------
UTILS_LIB_NAME="Utilities${DEBUG_SUFFIX}"
PKG_CHECK_MODULES(UTILS, $UTILS_LIB_NAME >= 1.0.0)
AC_SUBST([UTILS_CFLAGS])
AC_SUBST([UTILS_LIBS])

# ------------------------------------------------------------------------------
# Construct library name
# ------------------------------------------------------------------------------
MFM_LIB_NAME="${PACKAGE_NAME}${DEBUG_SUFFIX}"
AC_SUBST([MFM_LIB_NAME])

LIBRARY_VERSION=1:0:0
AC_SUBST([LIBRARY_VERSION])

# ------------------------------------------------------------------------------
# Declare output files
# ------------------------------------------------------------------------------
AC_CONFIG_FILES([Makefile
	src/Makefile
	test/Makefile
    ${MFM_LIB_NAME}.pc:MultiFrame.pc.in
])
AC_OUTPUT
