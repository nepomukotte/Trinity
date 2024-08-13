# ------------------------------------------------------------------------------
# Package setup
# ------------------------------------------------------------------------------

AC_INIT([DatabaseWidgets], [1.3.0], [frederic.chateau@cea.fr], [DatabaseWidgets])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([-Wall -Werror -Wno-portability-recursive])
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

INSTALL="${INSTALL} -C"

# ------------------------------------------------------------------------------
# Default flag reinitialization
# ------------------------------------------------------------------------------
CXXFLAGS="-Wall -Wextra -Wno-deprecated-declarations ${CXXFLAGS}"
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
AC_CHECK_HEADERS([cassert map memory sstream stdexcept string vector], [],
	[AC_MSG_ERROR([missing standard C++ header])])

# ------------------------------------------------------------------------------
# Enable silent build rules by default
# ------------------------------------------------------------------------------
m4_ifdef([AM_SILENT_RULES],
	[AM_SILENT_RULES([yes])],
	[AC_SUBST([AM_DEFAULT_VERBOSITY], [1])]
)

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
# Construct library name
# ------------------------------------------------------------------------------
DATABASE_WIDGETS_LIB_NAME="$PACKAGE_NAME${DEBUG_SUFFIX}"
CONFIG_WIDGETS_LIB_NAME="ConfigWidgets${DEBUG_SUFFIX}"
CONFIG_DATABASE_LIB_NAME="ConfigDatabase${DEBUG_SUFFIX}"
AC_SUBST([DATABASE_WIDGETS_LIB_NAME])
AC_SUBST([CONFIG_WIDGETS_LIB_NAME])
AC_SUBST([CONFIG_DATABASE_LIB_NAME])

# ------------------------------------------------------------------------------
# Include Qt stuff
# ------------------------------------------------------------------------------
AT_WITH_QT([], [$QT_CONFIG], [QT += widgets])
AT_REQUIRE_QT_VERSION([4.3])

# ------------------------------------------------------------------------------
# Check for ConfigDatabase and ConfigWidgets library
# ------------------------------------------------------------------------------
PKG_CHECK_MODULES(CONFIG_DB_AND_WIDGETS, $CONFIG_DATABASE_LIB_NAME >= 2.2.0 $CONFIG_WIDGETS_LIB_NAME >= 2.0.3)
AC_SUBST(CONFIG_DB_AND_WIDGETS_CFLAGS)
AC_SUBST(CONFIG_DB_AND_WIDGETS_LIBS)

# ------------------------------------------------------------------------------
# Declare output files
# ------------------------------------------------------------------------------
AC_CONFIG_FILES([Makefile
    src/Makefile
    ${DATABASE_WIDGETS_LIB_NAME}.pc:DatabaseWidgets.pc.in
])
AC_OUTPUT
