# ------------------------------------------------------------------------------
# Package setup
# ------------------------------------------------------------------------------

AC_INIT([TinyXml], [2.5.2], [http://sourceforge.net/projects/tinyxml])
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
CXXFLAGS_TINYXML="-std=c++98"
AC_SUBST([CXXFLAGS_TINYXML])

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
# Detect compilation mode: debug or release
# ------------------------------------------------------------------------------
AX_DEBUG_OR_RELEASE_BUILD

# ------------------------------------------------------------------------------
# Detect use of STL
# ------------------------------------------------------------------------------
AC_ARG_ENABLE([stl],
	AS_HELP_STRING([--enable-stl@<:@=yes|no@:>@], [Turn on STL support for strings. Default is yes]),
	[case "${enableval}" in
		yes) use_stl=true ;;
		no)  use_stl=false ;;
		*) AC_MSG_ERROR([bad value ${enableval} for --enable-stl]) ;;
	esac],
	[use_stl=true]
)
if test "$use_stl" = "true"
then
	DEFINE_FLAGS="-DTIXML_USE_STL"
#	AC_CHECK_HEADERS([iostream sstream string], [], AC_MSG_ERROR([missing standard C++ headers]))
else
	FULL_SUFFIX="${FULL_SUFFIX}l"
fi
AC_SUBST([DEFINE_FLAGS])

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
# Construct library name
# ------------------------------------------------------------------------------
TINYXML_LIB_NAME=${PACKAGE_TARNAME}
if test -n "${FULL_SUFFIX}"
then
	TINYXML_LIB_NAME="${TINYXML_LIB_NAME}-${FULL_SUFFIX}"
fi
AC_SUBST([TINYXML_LIB_NAME])

# ------------------------------------------------------------------------------
# Declare output files
# ------------------------------------------------------------------------------
AC_CONFIG_FILES([
	Makefile
	${TINYXML_LIB_NAME}.pc:tinyxml.pc.in
])
AC_OUTPUT
