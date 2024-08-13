# ------------------------------------------------------------------------------
# Package setup
# ------------------------------------------------------------------------------

AC_INIT([ConfigDatabase], [2.2.0], [frederic.chateau@cea.fr], [ConfigDatabase])
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
# Check for tool support
# ------------------------------------------------------------------------------
AC_PROG_CXX
AC_PROG_LIBTOOL

# ------------------------------------------------------------------------------
# Check for C++ support
# ------------------------------------------------------------------------------
AC_LANG([C++])
AC_CHECK_HEADERS([algorithm cassert cctype ctime exception functional iomanip \
    iostream istream map memory ostream sstream stdexcept string utility vector],
    [], [AC_MSG_ERROR([missing standard C++ header])])

# ------------------------------------------------------------------------------
# Check for database client libraries
# ------------------------------------------------------------------------------
AX_LIB_MYSQL([5.0.0])
if test "$found_mysql" != "yes" -o "$mysql_version_check" != "1"
then
	AC_MSG_ERROR([missing required Mysql client SDK]) 
fi

CCFG_BACKENDS_FLAGS="-DCCFG_USE_MYSQL=1"

dnl Disable Oracle backend by default
AC_ARG_ENABLE([oracle],
	[AS_HELP_STRING([--enable-oracle],[enable Oracle backend])],
	[use_oracle=$enableval],[use_oracle=no]
)
AS_IF([test "x$use_oracle" = "xyes"],
	[AX_LIB_ORACLE_OCI([9.0])
	 AS_IF([test "$oci_header_found" != "yes" -o "$oci_lib_found" != "yes" -o "$oracle_version_check" != "1"],
		[AC_MSG_ERROR([missing required Oracle instantclient SDK])],
		[AM_CONDITIONAL([CCFG_USE_ORACLE],true)
		 CCFG_BACKENDS_FLAGS+=" -DCCFG_USE_ORACLE=1"
		]
	 )
	],
	[AC_MSG_NOTICE([Oracle backend disabled on user request])
	 AM_CONDITIONAL([CCFG_USE_ORACLE],false)]
)

dnl Disable PostgreSql backend by default
USE_POSTGRESQL="0"
AC_ARG_ENABLE([postgresql],
	[AS_HELP_STRING([--enable-postgresql],[enable PostgreSQL backend])],
	[use_postgresql=$enableval],[use_postgresql=no]
)
AS_IF([test "x$use_postgresql" = "xyes"],
	[AX_LIB_POSTGRESQL([7.3])
	 AS_IF([test "$found_postgresql" != "yes" -o "$postgresql_version_check" != "1"],
		[AC_MSG_ERROR([missing required PostgreSQL client SDK])],
		[AM_CONDITIONAL([CCFG_USE_POSTGRESQL],true)
		 CCFG_BACKENDS_FLAGS+=" -DCCFG_USE_POSTGRESQL=1"
		]
	 )
	],
	[AC_MSG_NOTICE([PostgreSQL backend disabled on user request])
	 AM_CONDITIONAL([CCFG_USE_POSTGRESQL],false)]
)

dnl Disable SQLite3 backend by default
USE_SQLITE3="0"
AC_ARG_ENABLE([sqlite3],
	[AS_HELP_STRING([--enable-sqlite3],[enable SQLite3 backend])],
	[use_sqlite3=$enableval],[use_sqlite3=no]
)
AS_IF([test "x$use_sqlite3" = "xyes"],
	[AX_LIB_SQLITE3([3.1.0])
	 AS_IF([test "$found_sqlite3" != "yes" -o "$sqlite3_version_check" != "1"],
		[AC_MSG_ERROR([missing required SQLite3 client SDK])],
		[AM_CONDITIONAL([CCFG_USE_SQLITE3],true)
		 CCFG_BACKENDS_FLAGS+=" -DCCFG_USE_SQLITE3=1"
		]
	 )
	],
	[AC_MSG_NOTICE([SQLite3 backend disabled on user request])
	 AM_CONDITIONAL([CCFG_USE_SQLITE3],false)]
)

AC_SUBST([CCFG_BACKENDS_FLAGS])


# ------------------------------------------------------------------------------
# Enable silent build rules by default
# ------------------------------------------------------------------------------
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])],[AC_SUBST([AM_DEFAULT_VERBOSITY], [1])])

# ------------------------------------------------------------------------------
# Hands on configure flags to "make distcheck"
# ------------------------------------------------------------------------------
if test -n "$MYSQL_CONFIG"
then
	DISTCHECK_CONFIGURE_FLAGS+=" --with-mysql=$MYSQL_CONFIG"
fi
if test -n "$oracle_home_dir"
then
	DISTCHECK_CONFIGURE_FLAGS+=" --with-oracle=$oracle_home_dir"
fi
if test -n "$oracle_home_include_dir"
then
	DISTCHECK_CONFIGURE_FLAGS+=" --with-oracle-include=$oracle_home_include_dir"
fi
if test -n "$oracle_home_lib_dir"
then
	DISTCHECK_CONFIGURE_FLAGS+=" --with-oracle-lib=$oracle_home_lib_dir"
fi
if test -n "$soci_dir"
then
	DISTCHECK_CONFIGURE_FLAGS+=" --with-soci=$soci_dir"
fi
AC_SUBST(DISTCHECK_CONFIGURE_FLAGS)

# ------------------------------------------------------------------------------
# Detect compilation mode: debug or release
# ------------------------------------------------------------------------------
AX_DEBUG_OR_RELEASE_BUILD
if test "$enable_debug" = "yes"
then
	SOCI_BUILD="gcc-g"
else
	SOCI_BUILD="gcc"
fi

# ------------------------------------------------------------------------------
# Construct library name
# ------------------------------------------------------------------------------
CONFIG_DATABASE_LIB_NAME="${PACKAGE_NAME}${DEBUG_SUFFIX}"
CCFG_LIB_NAME="CompoundConfig${DEBUG_SUFFIX}"
AC_SUBST([CONFIG_DATABASE_LIB_NAME])
AC_SUBST([CCFG_LIB_NAME])

# ------------------------------------------------------------------------------
# Check for SOCI library
# ------------------------------------------------------------------------------
AX_LIB_SOCI([2_2], [$SOCI_BUILD])
AC_SUBST([SOCI_CPPFLAGS])
AC_SUBST([SOCI_LDFLAGS])

# ------------------------------------------------------------------------------
# Check for CompoundConfig library
# ------------------------------------------------------------------------------
PKG_CHECK_MODULES(COMPOUND_CONFIG, $CCFG_LIB_NAME >= 2.0.0)
AC_SUBST(COMPOUND_CONFIG_CFLAGS)
AC_SUBST(COMPOUND_CONFIG_LIBS)

# ------------------------------------------------------------------------------
# Declare output files
# ------------------------------------------------------------------------------
AC_CONFIG_FILES([Makefile
    src/Makefile
    test/Makefile
    ${CONFIG_DATABASE_LIB_NAME}.pc:ConfigDatabase.pc.in
])
AC_OUTPUT
