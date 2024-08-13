include $(top_srcdir)/am/cppcheck.am

AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4

AM_CPPFLAGS = -I@top_srcdir@/src $(DEFINE_FLAGS)
AM_CXXFLAGS = $(CXXFLAGS_TINYXML)

headers =  \
	src/tinyxml.h

sources = \
	src/tinyxml.cpp \
	src/tinyxmlerror.cpp \
	src/tinyxmlparser.cpp

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = @TINYXML_LIB_NAME@.pc

library_includedir = $(includedir)
library_include_HEADERS = $(headers)

lib_LTLIBRARIES = lib@TINYXML_LIB_NAME@.la
lib@TINYXML_LIB_NAME@_la_SOURCES = $(sources)
lib@TINYXML_LIB_NAME@_la_LDFLAGS = -release $(PACKAGE_VERSION) -no-undefined

check_PROGRAMS = xmltest
xmltest_SOURCES = test/xmltest.cpp
xmltest_LDADD = lib@TINYXML_LIB_NAME@.la

EXTRA_DIST = autogen.sh test/utf8test.gif test/utf8test.xml test/utf8testverify.xml
CLEANFILES = utf8test.xml utf8testout.xml utf8testverify.xml demotest.xml \
	test5.xml test6.xml test7.xml textfile.txt

check-am: check-before
check-before:
	cp @srcdir@/test/utf8test.xml @builddir@
	cp @srcdir@/test/utf8testverify.xml @builddir@

include $(srcdir)/Makefile.vxworks.am
include $(top_srcdir)/am/pkg-linux.am
include $(top_srcdir)/am/pkg-win.am
