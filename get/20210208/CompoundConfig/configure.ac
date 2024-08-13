# ------------------------------------------------------------------------------
# Package setup
# ------------------------------------------------------------------------------

AC_INIT([CompoundConfig], [2.2.0], [frederic.chateau@cea.fr])
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

AM_INIT_AUTOMAKE([-Wall -Werror])
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

INSTALL="${INSTALL} -C"

# ------------------------------------------------------------------------------
# Default flag reinitialization
# ------------------------------------------------------------------------------
CXXFLAGS="-Wall -Wextra ${CXXFLAGS}"
CFLAGS="-Wall -Wextra ${CFLAGS}"
CXXFLAGS_CCFG="-std=c++98"
AC_SUBST([CXXFLAGS_CCFG])

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
# Detect endianness
# ------------------------------------------------------------------------------
CCFG_ENDIANNESS_FLAGS=""
AC_C_BIGENDIAN(
	[AC_DEFINE([WORDS_BIG_ENDIAN],[1],[Platform is big-endian])
	 CCFG_ENDIANNESS_FLAGS="-DWORDS_BIG_ENDIAN=1"],
	[AC_DEFINE([WORDS_LITTLE_ENDIAN],[1],[Platform is little-endian])
	 CCFG_ENDIANNESS_FLAGS="-DWORDS_LITTLE_ENDIAN=1"],
	[AC_MSG_ERROR([Endianness cannot be defined!])]
	)
AC_SUBST([CCFG_ENDIANNESS_FLAGS])

# ------------------------------------------------------------------------------
# Check for BOOST libraries
# ------------------------------------------------------------------------------
AX_BOOST_BASE([1.40.0],[],
	[AC_MSG_ERROR([Boost headers are needed to build this library.
      Please go to www.boost.org to download them.])]
)

if test -n "$ac_boost_path"
then
	DISTCHECK_CONFIGURE_FLAGS="--with-boost=$ac_boost_path $DISTCHECK_CONFIGURE_FLAGS"
fi

# ------------------------------------------------------------------------------
# Enable silent build rules by default
# ------------------------------------------------------------------------------
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])],[AC_SUBST([AM_DEFAULT_VERBOSITY], [1])])

# ------------------------------------------------------------------------------
# Detect compilation mode: debug or release
# ------------------------------------------------------------------------------
AX_DEBUG_OR_RELEASE_BUILD

# ------------------------------------------------------------------------------
# Enable VxWorks portability
# ------------------------------------------------------------------------------
AX_VXWORKS_PORTABILITY

# ------------------------------------------------------------------------------
# Enable WINDOWS portability
# ------------------------------------------------------------------------------
AX_WINDOWS_PORTABILITY

# ------------------------------------------------------------------------------
# Construct library name
# ------------------------------------------------------------------------------
CCFG_LIB_NAME="${PACKAGE_NAME}${DEBUG_SUFFIX}";
TINYXML_LIB_NAME="tinyxml${DEBUG_SUFFIX}"
AC_SUBST([CCFG_LIB_NAME])
AC_SUBST([TINYXML_LIB_NAME])

# ------------------------------------------------------------------------------
# Enable or disable XML format (enabled by default)
# ------------------------------------------------------------------------------
BUILD_CCFG_FMT_XML=1
AC_ARG_ENABLE([xml-format],
	[AS_HELP_STRING([--disable-xml-format],[disable CompoundConfig XML format])],
	[want_xml_format=$enableval],[want_xml_format=yes]
)
AC_MSG_CHECKING([whether to build CompoundConfig XML format])
AC_MSG_RESULT([$want_xml_format])
AM_CONDITIONAL([BUILD_CCFG_FMT_XML], [test "$want_xml_format" = yes])
AS_IF([test "$want_xml_format" = yes],
[AC_DEFINE([BUILD_CCFG_FMT_XML], [1], [Build CompoundConfig XML format])],
[])


# ------------------------------------------------------------------------------
# Enable or disable binary format (enabled by default)
# ------------------------------------------------------------------------------
BUILD_CCFG_FMT_BINARY=1
AC_ARG_ENABLE([binary-format],
	[AS_HELP_STRING([--disable-binary-format],[disable CompoundConfig binary format])],
	[want_binary_format=$enableval],[want_binary_format=yes]
)
AC_MSG_CHECKING([whether to build CompoundConfig binary format])
AC_MSG_RESULT([$want_binary_format])
AM_CONDITIONAL([BUILD_CCFG_FMT_BINARY], [test "$want_binary_format" = yes])
AS_IF([test "$want_binary_format" = yes],
	[AC_DEFINE([BUILD_CCFG_FMT_BINARY], [1], [Build CompoundConfig binary format])],
	[])

# ------------------------------------------------------------------------------
# Check for TinyXml library
# ------------------------------------------------------------------------------
TINYXML_LIB_MIN_VERSION="2.5.2"
AS_IF([test "$want_xml_format" = yes],
	[PKG_CHECK_MODULES([TINYXML], [$TINYXML_LIB_NAME >= $TINYXML_LIB_MIN_VERSION])
	 TINYXML_REQ="$TINYXML_LIB_NAME >= $TINYXML_LIB_MIN_VERSION"],
	[TINYXML_REQ=""
	 TINYXML_CFLAGS=""
	 TINYXML_LIBS=""])
AC_SUBST([TINYXML_CFLAGS])
AC_SUBST([TINYXML_LIBS])
AC_SUBST([TINYXML_REQ])
AC_MSG_RESULT(["tinyxmllib flags = " $TINYXML_LIBS])

# ------------------------------------------------------------------------------
# Declare output files
# ------------------------------------------------------------------------------
AC_CONFIG_FILES([Makefile
    src/Makefile
    test/Makefile
    ${CCFG_LIB_NAME}.pc:CompoundConfig.pc.in
])
AC_OUTPUT
