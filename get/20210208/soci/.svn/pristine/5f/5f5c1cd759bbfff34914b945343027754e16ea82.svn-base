# $Id: Makefile.am,v 1.3 2006/11/27 09:59:32 mloskot Exp $
#
# Author: Mateusz Loskot <mateusz@loskot.net>
#
INCLUDES = -I$(top_srcdir)/src/core
AM_CXXFLAGS = -fPIC
AM_CPPFLAGS = @POSTGRESQL_CFLAGS@

if POSTGRESQL
lib_LTLIBRARIES	= \
	lib@SOCI_POSTGRESQL_LIBRARY_NAME@.la

lib@SOCI_POSTGRESQL_LIBRARY_NAME@_la_LDFLAGS = \
	-release @LIBTOOL_RELEASE@ -no-undefined @POSTGRESQL_LDFLAGS@

lib@SOCI_POSTGRESQL_LIBRARY_NAME@_la_LIBADD = \
	@top_builddir@/src/core/lib@SOCI_CORE_LIBRARY_NAME@.la

lib@SOCI_POSTGRESQL_LIBRARY_NAME@_la_SOURCES = \
	blob.cpp \
	common.cpp \
	factory.cpp \
	row-id.cpp \
	session.cpp \
	standard-into-type.cpp \
	standard-use-type.cpp \
	statement.cpp \
	vector-into-type.cpp \
	vector-use-type.cpp

soci_postgresqldir = $(includedir)/soci/postgresql

soci_postgresql_HEADERS = \
	common.h \
	soci-postgresql.h
endif

EXTRA_DIST = Makefile.basic

MAINTAINERCLEANLOCAL = Makefile.in

maintainer-clean-local:
	-test -z "$(MAINTAINERCLEANLOCAL)" || rm -rf $(MAINTAINERCLEANLOCAL)
