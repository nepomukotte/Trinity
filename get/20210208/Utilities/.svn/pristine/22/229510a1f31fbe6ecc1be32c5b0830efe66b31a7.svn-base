# ------------------------------------------------------------------------------
# Package setup
# ------------------------------------------------------------------------------

AC_INIT([Utilities], [1.3.1], [frederic.chateau@cea.fr;patrick.sizun@cea.fr])
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
# Default flag reinitialization
# ------------------------------------------------------------------------------
CXXFLAGS="-Wall -Wextra ${CXXFLAGS}"
CFLAGS="-Wall -Wextra ${CFLAGS}"
CXXFLAGS_UTILS="-std=c++98 -DBOOST_NO_CXX11_TEMPLATE_ALIASES"
AC_SUBST([CXXFLAGS_UTILS])

# ------------------------------------------------------------------------------
# Check for tool support
# ------------------------------------------------------------------------------
AC_PROG_CXX
AC_PROG_LIBTOOL

# ------------------------------------------------------------------------------
# Enable VxWorks portability
# ------------------------------------------------------------------------------
AX_VXWORKS_PORTABILITY

# ------------------------------------------------------------------------------
# Enable RTEMS portability
# ------------------------------------------------------------------------------
AX_RTEMS_PORTABILITY

# ------------------------------------------------------------------------------
# Enable WINDOWS portability
# ------------------------------------------------------------------------------
AX_WINDOWS_PORTABILITY

# ------------------------------------------------------------------------------
# Check for C++ support
# ------------------------------------------------------------------------------
AC_LANG([C++])
AC_CHECK_HEADERS([algorithm bitset cassert cctype cmath cstdlib cstring \
   deque exception fstream functional iomanip ios iterator limits list locale \
   map memory ostream set sstream string typeinfo vector], [],
   [AC_MSG_ERROR([missing standard C++ header])])


# ------------------------------------------------------------------------------
# Detect endianness
# ------------------------------------------------------------------------------
AC_C_BIGENDIAN(
	[AC_DEFINE([WORDS_BIG_ENDIAN],[1],[Platform is big-endian])],
	[AC_DEFINE([WORDS_LITTLE_ENDIAN],[1],[Platform is little-endian])],
	[AC_MSG_ERROR([Endianness cannot be defined!])]
	)

# ------------------------------------------------------------------------------
# Detect compilation mode: debug or release
# ------------------------------------------------------------------------------
AX_DEBUG_OR_RELEASE_BUILD

# ------------------------------------------------------------------------------
# Check for BOOST libraries
# -----------------------------------------------------------------------------
AX_BOOST_BASE([1.40.0],[],
	[AC_MSG_ERROR([Boost headers are needed to build this library.
      Please go to www.boost.org to download them.])]
)

if test -n "$ac_boost_path"
then
        DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --with-boost=$ac_boost_path"
fi

# ------------------------------------------------------------------------------
# Detect log level in Utilities Logger
# ------------------------------------------------------------------------------
AX_UTL_LOG

# ------------------------------------------------------------------------------
# Enable silent build rules by default
# ------------------------------------------------------------------------------
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])],[AC_SUBST([AM_DEFAULT_VERBOSITY], [1])])

# ------------------------------------------------------------------------------
# Define library to link with for clock and timing functions
# ------------------------------------------------------------------------------
AS_CASE(
	[$host_os],
	[*linux*], [CLOCK_LIBS="-lrt"],
	[CLOCK_LIBS=""])
AC_SUBST([CLOCK_LIBS])

# ------------------------------------------------------------------------------
# Construct library name
# ------------------------------------------------------------------------------
UTILS_LIB_NAME="${PACKAGE_NAME}${DEBUG_SUFFIX}"
AC_SUBST([UTILS_LIB_NAME])

UTILS_LIB_VERSION=0:0:0
AC_SUBST([UTILS_LIB_VERSION])

# ------------------------------------------------------------------------------
# Hands on configure flags to "make distcheck"
# ------------------------------------------------------------------------------
if test "$enable_shared" = "no"
then
	DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --disable-shared"
fi
AC_SUBST([DISTCHECK_CONFIGURE_FLAGS])

# ------------------------------------------------------------------------------
# Declare output files
# ------------------------------------------------------------------------------
AC_CONFIG_FILES([Makefile
    ${UTILS_LIB_NAME}.pc:Utilities.pc.in
])
AC_OUTPUT
