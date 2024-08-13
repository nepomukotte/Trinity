AUTOMAKE_OPTIONS = subdir-objects

headers = CCfg/Io/Db/ConfigDaoTest.hpp \
	CCfg/Io/Db/Sleep.h \
	CCfg/Io/Db/TableBuilderTest.h \
	TestParameters.h \
	tut.h \
	tut_reporter.h

sources = CCfg/Io/DatabaseIoFactoryTest.cpp \
	CCfg/Io/Db/MySqlTableBuilderTest.cpp \
	CCfg/Io/Db/MySqlConfigDaoTest.cpp \
	CCfg/Io/Db/TableBuilderTest.cpp \
	main.cpp \
	TestParameters.cpp

if CCFG_USE_ORACLE
sources += 	CCfg/Io/Db/OracleTableBuilderTest.cpp \
	CCfg/Io/Db/OracleConfigDaoTest.cpp
endif

if CCFG_USE_POSTGRESQL
sources += 	CCfg/Io/Db/PostgreSqlTableBuilderTest.cpp \
	CCfg/Io/Db/PostgreSqlConfigDaoTest.cpp
endif

if CCFG_USE_SQLITE3
sources += 	CCfg/Io/Db/Sqlite3TableBuilderTest.cpp \
	CCfg/Io/Db/Sqlite3ConfigDaoTest.cpp
endif

AM_CPPFLAGS = $(CCFG_BACKENDS_FLAGS) $(COMPOUND_CONFIG_CFLAGS) -I@top_srcdir@/src

#TESTS = ConfigDatabaseTest
noinst_PROGRAMS = ConfigDatabaseTest
ConfigDatabaseTest_SOURCES = $(sources)
ConfigDatabaseTest_LDADD = $(COMPOUND_CONFIG_LIBS) @top_builddir@/src/lib@CONFIG_DATABASE_LIB_NAME@.la

EXTRA_DIST = testconfig1.xml testconfig2.xml $(headers)
CLEANFILES = testconfig1.xml testconfig2.xml

cppcheck:

#check-am: check-before
#check-before:
#	cp @srcdir@/*.xml @builddir@
#	cp @srcdir@/connections.txt @builddir@
