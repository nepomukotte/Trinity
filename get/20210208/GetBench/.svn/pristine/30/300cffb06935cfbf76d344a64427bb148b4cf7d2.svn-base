# ------------------------------------------------------------------------------
# Package setup
# ------------------------------------------------------------------------------

AC_INIT([gettestbench], [2.3.0], [shebli.anvar@cea.fr;patrick.sizun@cea.fr])
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
DX_INIT_DOXYGEN([$PACKAGE_NAME], [doc/doxyfile], [doc])

AM_INIT_AUTOMAKE([-Wall -Werror])
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
AC_PROG_LN_S
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
# Enable VxWorks portability
# ------------------------------------------------------------------------------
AX_VXWORKS_PORTABILITY

# ------------------------------------------------------------------------------
# Check for BOOST libraries
# ------------------------------------------------------------------------------
boost_version_min_req=1.40.0
AX_BOOST_BASE([$boost_version_min_req],[],
	[AC_MSG_ERROR([Boost headers are needed to build this library.
      Please go to www.boost.org to download them.])]
)

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
# Check for specific Boost libraries
# ------------------------------------------------------------------------------
AX_BOOST_SYSTEM
if test "$host_os" != "vxworks"
then
  AX_BOOST_THREAD
  if test "$enable_ecc" != "no"
  then
      AX_BOOST_THREAD
	  AX_BOOST_FILESYSTEM
	  AX_BOOST_PROGRAM_OPTIONS
  fi
fi
# ------------------------------------------------------------------------------
# Enable silent build rules by default
# ------------------------------------------------------------------------------
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])],[AC_SUBST([AM_DEFAULT_VERBOSITY], [1])])

# ------------------------------------------------------------------------------
# Let users choose the log level
# ------------------------------------------------------------------------------
AX_UTL_LOG

# ------------------------------------------------------------------------------
# Detect compilation mode: debug or release
# ------------------------------------------------------------------------------
AX_DEBUG_OR_RELEASE_BUILD
AC_SUBST([DEBUG_SUFFIX])

# ------------------------------------------------------------------------------
# Enable Windows portability
# ------------------------------------------------------------------------------
AX_WINDOWS_PORTABILITY

# ------------------------------------------------------------------------------
# Check for pkg-config
# ------------------------------------------------------------------------------
PKG_PROG_PKG_CONFIG

# ------------------------------------------------------------------------------
# Define variable for Apple hosts to disable features not yet ported
# ------------------------------------------------------------------------------
AM_CONDITIONAL([APPLE], [test "$host_vendor" = "apple"])

# ------------------------------------------------------------------------------
# Define variable for ARM hosts to disable features not yet ported
# ------------------------------------------------------------------------------
AM_CONDITIONAL([ARM], [test "$host_cpu" = "arm"])
AS_IF([test "$host_cpu" = "arm"], [CPPFLAGS="${CPPFLAGS} -DZEDBOARD=1"])

# ------------------------------------------------------------------------------
# Options to enable build of target software, ECC server, ECC client
# Better than relying on host platform to decide what to build
# ------------------------------------------------------------------------------
AC_ARG_ENABLE([cobo],
	AS_HELP_STRING([--enable-cobo], [enables build of (CoBo) embedded target server])]
)
AM_CONDITIONAL([BUILD_TARGET_SERVER], [test "x$enable_cobo" = "xyes"])

AC_ARG_ENABLE([ecc],
	AS_HELP_STRING([--disable-ecc], [disables build of ECC server and client])
)
AS_IF([test "x$enable_ecc" != "xno"],
	[AS_IF([test "$host_os" = "vxworks" -o "$host_cpu" = "arm"],
		[AC_MSG_ERROR([The ECC server and client cannot be build on VxWorks nor ARM. Use --disable-ecc.])])]
)
AM_CONDITIONAL([BUILD_ECC], [test "x$enable_ecc" != "xno"])

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
	AC_MSG_WARN([Remote logging features will not be built!])])
AM_CONDITIONAL([HAVE_ICESTORM], [test "x$have_icestorm" = "xyes"])

# ------------------------------------------------------------------------------
# Check for log4cxx library
# ------------------------------------------------------------------------------
MORD_LIB_LOG4CXX([have_log4cxx=yes],
	[have_log4cxx=no])
AM_CONDITIONAL([HAVE_LOG4CXX], [test "x$have_log4cxx" = "xyes"])

# ------------------------------------------------------------------------------
# Check for ConfigDatabase if required (disabled by default)
# ------------------------------------------------------------------------------
db_support=false
if test "$enable_ecc" != "no"
then
AC_ARG_ENABLE([database-support],
	[AS_HELP_STRING([--enable-database-support],[enable ConfigDatabase support])],
	[use_config_database=$enableval],[use_config_database=no]
)
AS_IF([test "x$use_config_database" = "xyes"],
	[PKG_CHECK_MODULES([CONFIG_DATABASE], [ConfigDatabase$DEBUG_SUFFIX],
		[db_support=true])
		CPPFLAGS="$CPPFLAGS -DCFG_DB_SUPPORT_ENABLED=1"],
	[AC_MSG_NOTICE([support for configuration database disabled on user request])
	 db_support=false]
)
fi
AM_CONDITIONAL([CFG_DB_SUPPORT_ENABLED], [test x$db_support = xtrue])


# ------------------------------------------------------------------------------
# Check for MDaq libraries
# ------------------------------------------------------------------------------
PKG_CHECK_MODULES([MDAQ_UTL], [mdaq-utl$DEBUG_SUFFIX >= 2.1.0])
AC_SUBST([MDAQ_UTL_CFLAGS])
AC_SUBST([MDAQ_UTL_LIBS])

if test "$enable_ecc" != "no"
then
	PKG_CHECK_MODULES([MDAQ_ECCSERVER], [mdaq-eccserver$DEBUG_SUFFIX >= 2.1.0])
	AC_SUBST([MDAQ_ECCSERVER_CFLAGS])
	AC_SUBST([MDAQ_ECCSERVER_LIBS])
	
	PKG_CHECK_MODULES([MDAQ_CLI], [mdaq-cli$DEBUG_SUFFIX])
	AC_SUBST([MDAQ_CLI_CFLAGS])
	AC_SUBST([MDAQ_CLI_LIBS])
	
	PKG_CHECK_MODULES([MDAQ_IF], [mdaq-interfaces$DEBUG_SUFFIX])
	AC_SUBST([MDAQ_IF_CFLAGS])
	AC_SUBST([MDAQ_IF_LIBS])
fi

PKG_CHECK_MODULES([MDAQ_DAQ], [mdaq-daq$DEBUG_SUFFIX])
AC_SUBST([MDAQ_DAQ_CFLAGS])
AC_SUBST([MDAQ_DAQ_LIBS])

PKG_CHECK_MODULES([MDAQ_HWSERVER], [mdaq-hwserver$DEBUG_SUFFIX])
AC_SUBST([MDAQ_HWSERVER_CFLAGS])
AC_SUBST([MDAQ_HWSERVER_LIBS])

PKG_CHECK_MODULES([MDAQ_HWCLIENT], [mdaq-hwclient$DEBUG_SUFFIX])
AC_SUBST([MDAQ_HWCLIENT_CFLAGS])
AC_SUBST([MDAQ_HWCLIENT_LIBS])

# ------------------------------------------------------------------------------
# Check for location of MDaq slice definitions
# ------------------------------------------------------------------------------
AC_MSG_CHECKING([for MDaq slice directory])
MDAQ_SLICEDIR=`$PKG_CONFIG --variable=slicedir mdaq-interfaces$DEBUG_SUFFIX`
if test -z $MDAQ_SLICEDIR; then
  # In the long term, this should become a fatal error
  AC_MSG_RESULT([not found])
  MDAQ_SLICEDIR="."
  AC_MSG_WARN([You may need to upgrade the MDaq library!])
else
  AC_MSG_RESULT([found])
fi
AC_SUBST([MDAQ_SLICEDIR])

# ------------------------------------------------------------------------------
# Check for MFM library (only on Linux)
# ------------------------------------------------------------------------------
AC_ARG_WITH([mfm-support],
	[AS_HELP_STRING([--with-mfm-support],[enable MFM-related features])],
	[],
	[with_mfm_support=no]
)
AS_IF([test "$enable_ecc" != "no" -o "x$with_mfm_support" != "xno"],
	[PKG_CHECK_MODULES([MFM], [MultiFrame$DEBUG_SUFFIX >= 2.7.0])
	 AC_SUBST([MFM_CFLAGS])
	 AC_SUBST([MFM_LIBS])
	],
	[have_mfm=no]
)
AM_CONDITIONAL([HAVE_MFM], [test "x$have_mfm" != "xno"], [Whether the MFM library has been searched and found.])
AS_IF([test "x$have_mfm" != "xno"], AC_SUBST([HAVE_MFM], [1]), [])

# ------------------------------------------------------------------------------
# Check for Utilities library
# ------------------------------------------------------------------------------
PKG_CHECK_MODULES(UTILS, Utilities$DEBUG_SUFFIX >= 1.2.1)
AC_SUBST([UTILS_CFLAGS])
AC_SUBST([UTILS_LIBS])

# ------------------------------------------------------------------------------
# Check for CompoundConfig library
# ------------------------------------------------------------------------------
PKG_CHECK_MODULES([COMPOUND_CONFIG], [CompoundConfig$DEBUG_SUFFIX >= 2.1.2])
AC_SUBST([COMPOUND_CONFIG_CFLAGS])
AC_SUBST([COMPOUND_CONFIG_LIBS])

# ------------------------------------------------------------------------------
# Check for Mutant library
# ------------------------------------------------------------------------------
dnl PKG_CHECK_MODULES([MUTANT], [mutant$DEBUG_SUFFIX >= 1.0.0],[],[AC_MSG_ERROR([You need to install the Mutant software library.])])
dnl AC_SUBST([MUTANT_CFLAGS])
dnl AC_SUBST([MUTANT_LIBS])

# ------------------------------------------------------------------------------
# Check for DHSM library (except on VxWorks and Zed board)
# ------------------------------------------------------------------------------
if test "$enable_ecc" != "no"
then
	PKG_CHECK_MODULES([DHSM], [dhsm$DEBUG_SUFFIX >= 1.0.0],[],[
		AC_MSG_ERROR([You need to install the Distributed Hierarchical State Machine (dhsm) library.
See https://dsm-trac.cea.fr/mordicus/wiki])])
	AC_SUBST([DHSM_CFLAGS])
	AC_SUBST([DHSM_LIBS])
fi

# ------------------------------------------------------------------------------
# Check for gSOAP tools
# ------------------------------------------------------------------------------
have_gsoap=no

if test "$enable_ecc" != "no"
then
	AX_LIB_GSOAP([have_gsoap=yes], [AC_MSG_WARN([Could not find the gSOAP toolkit])])
fi

if test "x$have_gsoap" = "xyes"
then
dnl Starting from 2.8.13, old-style C++ service proxies and objects are deprecated and -z1 option must be used to generate them
	AC_MSG_CHECKING([whether gSOAP version is older than 20813])
	SOAP_CPP2_FLAGS=""
	if test $GSOAP_VERSION -ge 20813
	then
		AC_MSG_RESULT([no])
		SOAP_CPP2_FLAGS="-z1"
	else
		AC_MSG_RESULT([yes])
	fi
	AC_SUBST([SOAP_CPP2_FLAGS])
fi
AM_CONDITIONAL([HAVE_SOAP], [test "$have_gsoap" = "yes"])

# ------------------------------------------------------------------------------
# Check if compiler support flag to silence warning about gsoap using variable length arrays
# ------------------------------------------------------------------------------
saved_CPPFLAGS="$CPPFLAGS"
VLA_CPPFLAGS="-Wno-vla"
CPPFLAGS="$VLA_CPPFLAGS"
AC_MSG_CHECKING([whether compiler supports -Wno-vla])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([])],
    [AC_MSG_RESULT([yes])],
    [AC_MSG_RESULT([no])]
    [VLA_CPPFLAGS=""]
)
CPPFLAGS="$saved_CPPFLAGS"
AC_SUBST([VLA_CPPFLAGS])

# ------------------------------------------------------------------------------
# Check for BSP XGPIO headers needed for CoBo/Mutant link
# ------------------------------------------------------------------------------
AC_ARG_WITH([bsp-drv-xsrc],
  [AS_HELP_STRING([--with-bsp-drv-xsrc=DIR], [Xilinx BSP device driver source directory])],
  [XGPIO_CFLAGS="-I$withval -I$withval/../.."], [XGPIO_CFLAGS=""]
)
saved_CPPFLAGS=$CPPFLAGS
CPPFLAGS="$CPPFLAGS $XGPIO_CFLAGS"
AC_CHECK_HEADER([xgpio.h], [], [have_xgpio_headers=no])
AC_CHECK_HEADER([xparameters.h], [], [have_xgpio_headers=no])
CPPFLAGS="$saved_CPPFLAGS"
AS_IF([test x"$have_xgpio_headers" != xno],
 [XGPIO_CFLAGS="-DHAVE_XGPIO $XGPIO_CFLAGS"
  AC_SUBST([XGPIO_CFLAGS])], [])
AS_IF([test "$have_xgpio_headers" != "no" -o "$host_os" != "vxworks"], [], [AC_MSG_WARN([CoBo's Mutant link manager will not be built!])])
AM_CONDITIONAL([HAVE_XGPIO], [test x$have_xgpio_headers != xno])

# ------------------------------------------------------------------------------
# Declare output files
# ------------------------------------------------------------------------------
AC_CONFIG_FILES([
	Makefile
	data/Makefile
	getbench-daq${DEBUG_SUFFIX}.pc:getbench-daq.pc.in
	getbench-cfg${DEBUG_SUFFIX}.pc:getbench-cfg.pc.in
    getbench-eccclient${DEBUG_SUFFIX}.pc:getbench-eccclient.pc.in
    getbench-eccserver${DEBUG_SUFFIX}.pc:getbench-eccserver.pc.in
])
AC_OUTPUT
