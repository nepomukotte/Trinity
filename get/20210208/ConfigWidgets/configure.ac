# ------------------------------------------------------------------------------
# Package setup
# ------------------------------------------------------------------------------

AC_INIT([ConfigWidgets], [2.0.3], [frederic.chateau@cea.fr])
AC_CONFIG_HEADERS([config.h])
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
# Default flag reinitialization
# ------------------------------------------------------------------------------
CXXFLAGS="-Wall -Wextra ${CXXFLAGS}"
CFLAGS="-Wall -Wextra ${CFLAGS}"
CXXFLAGS_CCFG_WIDGETS="-Wno-deprecated-declarations"
AC_SUBST([CXXFLAGS_CCFG_WIDGETS])

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

AX_UTL_LOG

# ------------------------------------------------------------------------------
# Construct library name
# ------------------------------------------------------------------------------
CONFIG_WIDGETS_LIB_NAME="${PACKAGE_NAME}${DEBUG_SUFFIX}"
CCFG_LIB_NAME="CompoundConfig${DEBUG_SUFFIX}"
AC_SUBST([CONFIG_WIDGETS_LIB_NAME])
AC_SUBST([CCFG_LIB_NAME])

# ------------------------------------------------------------------------------
# Include Qt stuff
# ------------------------------------------------------------------------------
AT_WITH_QT([], [$QT_CONFIG], [QT += widgets])
AT_REQUIRE_QT_VERSION([4.3])
if test -n "$QT_PATH"
then
	DISTCHECK_CONFIGURE_FLAGS="--with-qt=$QT_PATH $DISTCHECK_CONFIGURE_FLAGS"
fi

# ------------------------------------------------------------------------------
# Check for CompoundConfig library
# ------------------------------------------------------------------------------
PKG_CHECK_MODULES(COMPOUND_CONFIG, $CCFG_LIB_NAME >= 2.1.0)
AC_SUBST(COMPOUND_CONFIG_CFLAGS)
AC_SUBST(COMPOUND_CONFIG_LIBS)

# ------------------------------------------------------------------------------
# Check for Utilities library
# ------------------------------------------------------------------------------
PKG_CHECK_MODULES(UTILS, Utilities$DEBUG_SUFFIX >= 1.1.0)
AC_SUBST([UTILS_CFLAGS])
AC_SUBST([UTILS_LIBS])

# ------------------------------------------------------------------------------
# Declare output files
# ------------------------------------------------------------------------------
AC_CONFIG_FILES([Makefile
    src/Makefile
    ${CONFIG_WIDGETS_LIB_NAME}.pc:ConfigWidgets.pc.in
])
AC_OUTPUT
