AC_INIT([ConfigServer], [1.1.0], [frederic.chateau@cea.fr])
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

# ------------------------------------------------------------------------------
# Enable silent build rules by default
# ------------------------------------------------------------------------------
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])],[AC_SUBST([AM_DEFAULT_VERBOSITY], [1])])

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
AX_PKG_SWIG([1.3.17], [], [AC_MSG_ERROR([SWIG is required to build.])])

# ------------------------------------------------------------------------------
# Check for Java Native Interface (JNI) header
# ------------------------------------------------------------------------------
AX_JNI_INCLUDE_DIR
JNI_CPPFLAGS=""
for JNI_INCLUDE_DIR in $JNI_INCLUDE_DIRS
do
	JNI_CPPFLAGS="$JNI_CPPFLAGS -I$JNI_INCLUDE_DIR"
done
AC_SUBST([JNI_CPPFLAGS])

# ------------------------------------------------------------------------------
# Check for C++ support
# ------------------------------------------------------------------------------
AC_LANG([C++])

# ------------------------------------------------------------------------------
# Detect compilation mode: debug or release
# ------------------------------------------------------------------------------
AX_DEBUG_OR_RELEASE_BUILD
AC_SUBST([DEBUG_SUFFIX])

# ------------------------------------------------------------------------------
# Check for BOOST libraries
# ------------------------------------------------------------------------------
AX_BOOST_BASE([1.33.1])
if test "$have_boost" = "no" ; then
   AC_MSG_ERROR([Boost headers are needed to build this library.
      Please go to www.boost.org to download them.])
fi

CPPFLAGS="-DBOOST_NO_CWCHAR ${CPPFLAGS}"
if test -n "$ac_boost_path"
then
	DISTCHECK_CONFIGURE_FLAGS="--with-boost=$ac_boost_path $DISTCHECK_CONFIGURE_FLAGS"
fi

# ------------------------------------------------------------------------------
# Enable Ice implementation
# ------------------------------------------------------------------------------
ACX_PTHREAD
AX_LIB_ICE("default-disabled")
AX_LIB_ICEE("default-disabled")
if test "$have_ice" = "yes" -o "$have_icee" = "yes"
then
	use_ice="yes"
	ICE_DEFINE="-DICE_EDITION=$ICE_EDITION -D'ICE_HEADER(header)=<ICE_EDITION/header>'"
	AC_SUBST([ICE_DEFINE])
fi
AM_CONDITIONAL([ICE], [test "$use_ice" = "yes"])

# ------------------------------------------------------------------------------
# Check for CompoundConfig library
# ------------------------------------------------------------------------------
PKG_CHECK_MODULES(COMPOUND_CONFIG, CompoundConfig >= 2.1.0)
AC_SUBST([COMPOUND_CONFIG_CFLAGS])
AC_SUBST([COMPOUND_CONFIG_LIBS])

# ------------------------------------------------------------------------------
# Declare output files
# ------------------------------------------------------------------------------
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
