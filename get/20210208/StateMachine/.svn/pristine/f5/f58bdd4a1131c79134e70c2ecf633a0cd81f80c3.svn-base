# ------------------------------------------------------------------------------
# Package setup
# ------------------------------------------------------------------------------

AC_INIT([dhsm], [1.1.0], [frederic.chateau@cea.fr])
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
CXXFLAGS_DHSM="-std=c++98"
AC_SUBST([CXXFLAGS_DHSM])

# ------------------------------------------------------------------------------
# Check for tool support
# ------------------------------------------------------------------------------
AC_PROG_CXX
AC_PROG_LIBTOOL

# ------------------------------------------------------------------------------
# Check for C++ support
# ------------------------------------------------------------------------------
AC_LANG([C++])

# ------------------------------------------------------------
# Enable silent build rules by default
# ------------------------------------------------------------------------------
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])],[AC_SUBST([AM_DEFAULT_VERBOSITY], [1])])

# ------------------------------------------------------------------------------
# Detect compilation mode: debug or release, and detect log level
# ------------------------------------------------------------------------------
AX_DEBUG_OR_RELEASE_BUILD
AC_SUBST([DEBUG_SUFFIX])
AX_UTL_LOG

# ------------------------------------------------------------------------------
# Enable VxWorks portability
# ------------------------------------------------------------------------------
AX_VXWORKS_PORTABILITY

# ------------------------------------------------------------------------------
# Enable RTEMS portability
# ------------------------------------------------------------------------------
AX_RTEMS_PORTABILITY

# ------------------------------------------------------------------------------
# Enable Ice implementation
# ------------------------------------------------------------------------------
AX_LIB_ICE_OR_ICEE

# ------------------------------------------------------------------------------
# Enable dhsm-master library build
# ------------------------------------------------------------------------------
AC_ARG_ENABLE([distributed-sm-libs],
	AS_HELP_STRING([--enable-distributed-sm-libs=(no|none|servant|master|all|yes)],
		[selects which libraries for distributed state machine to build.
		"master", "all" and "yes" have the same effect, "no" and "none" also. default is "no".]),
	[case "${enableval}" in
		yes)	enabled_dsm_libs=all ;;
		all)	enabled_dsm_libs=all ;;
		master)	enabled_dsm_libs=all ;;
		servant)enabled_dsm_libs=servant ;;
		none)	enabled_dsm_libs=none ;;
		no)		enabled_dsm_libs=none ;;
		*)		AC_MSG_ERROR([bad value ${enableval} for --enable-distributed-sm-libs])
	esac],
	[enabled_dsm_libs=none]
)
if test "$enabled_dsm_libs" = "all" -a "$have_ice" != "yes"
then
	AC_MSG_ERROR([The master library of the distributed state machine requires Ice (not IceE) to work.])
fi
AM_CONDITIONAL([MASTER], [test "$enabled_dsm_libs" = "all"])

if test "$enabled_dsm_libs" != "none" -a "$have_ice" != "yes" -a "$have_icee" != "yes"
then
	AC_MSG_ERROR([The servant library of the distributed state machine requires Ice or IceE to work.])
fi
AM_CONDITIONAL([SERVANT], [test "$enabled_dsm_libs" != "none"])

# ------------------------------------------------------------------------------
# Check for Utilities library
# ------------------------------------------------------------------------------
PKG_CHECK_MODULES(UTILS, Utilities$DEBUG_SUFFIX >= 1.2.2)
AC_SUBST([UTILS_CFLAGS])
AC_SUBST([UTILS_LIBS])

# ------------------------------------------------------------------------------
# Hands on configure flags to "make distcheck"
# ------------------------------------------------------------------------------
if test "$enabled_dsm_libs" != "none"
then
	DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS --enable-distributed-sm-libs=$enabled_dsm_libs"
fi
AC_SUBST([DISTCHECK_CONFIGURE_FLAGS])

# ------------------------------------------------------------------------------
# Declare output files
# ------------------------------------------------------------------------------
AC_CONFIG_FILES([test/test_dhsm.sh], [chmod +x test/test_dhsm.sh])
AC_CONFIG_FILES([
	Makefile
    src/Makefile
    test/Makefile
    dhsm${DEBUG_SUFFIX}.pc:dhsm.pc.in
])
if test "$enabled_dsm_libs" = "all"
then
	AC_CONFIG_FILES([dhsm-master${DEBUG_SUFFIX}.pc:dhsm-master.pc.in])
fi
if test "$enabled_dsm_libs" != "none"
then
	AC_CONFIG_FILES([dhsm-servant${DEBUG_SUFFIX}.pc:dhsm-servant.pc.in])
fi
AC_OUTPUT
