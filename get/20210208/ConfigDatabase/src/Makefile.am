include $(top_srcdir)/am/cppcheck.am

AUTOMAKE_OPTIONS = subdir-objects

headers = CCfg/Io/Db/Config.h \
	CCfg/Io/Db/AbstractEntity.h \
	CCfg/Io/Db/ConfigEntity.h \
	CCfg/Io/Db/ConfigDatabase.h \
	CCfg/Io/Db/ConfigAttribute.h \
	CCfg/Io/Db/ConfigStructureVersion.h \
	CCfg/Io/Db/ConfigVersion.h \
	CCfg/Io/Db/ConfigClass.h \
	CCfg/Io/Db/ConfigClassInfo.h \
	CCfg/Io/Db/Site.h \
	CCfg/Io/Db/Run.h \
	CCfg/Io/Db/RunCondition.h \
	CCfg/Io/Db/EntitySets.h \
	CCfg/Io/Db/ConfigDao.h \
	CCfg/Io/Db/AbstractEntityBatch.h \
	CCfg/Io/Db/Transaction.h \
	CCfg/Io/Db/MySqlTableBuilder.h \
	CCfg/Io/Db/TableBuilder.h \
	CCfg/Io/Db/MySqlConfigDao.h \
	CCfg/Io/Db/ConfigBatch.h \
	CCfg/Io/Db/ConfigStructureVersionBatch.h \
	CCfg/Io/Db/ConfigVersionBatch.h \
	CCfg/Io/Db/ConfigClassBatch.h \
	CCfg/Io/Db/ConfigClassInfoBatch.h \
	CCfg/Io/Db/ConfigAttributeBatch.h \
	CCfg/Io/Db/SiteBatch.h \
	CCfg/Io/Db/RunBatch.h \
	CCfg/Io/Db/RunConditionBatch.h \
	CCfg/Io/Db/ExtractionProcess.h \
	CCfg/Io/Db/AbstractUpdateController.h \
	CCfg/Io/Db/InsertionProcess.h \
	CCfg/Io/Db/PrimaryKeyGenerator.h \
	CCfg/Io/Db/Exception.h \
	CCfg/Io/Db/ObjRelMapping.hpp \
	CCfg/Io/AbstractLoginCallback.h \
	CCfg/Io/DatabaseFormat.h \
	CCfg/Io/DatabaseStorage.h \
	CCfg/Io/DatabaseFormatFactory.h \
	CCfg/Io/DatabaseStorageFactory.h \
	Utils/Date/StdTmHelpers.h \
	Utils/Stream/Require.h

sources = CCfg/Io/Db/MySqlConfigDao.cpp \
	CCfg/Io/Db/AbstractEntity.cpp \
	CCfg/Io/Db/ConfigEntity.cpp \
	CCfg/Io/Db/ConfigDatabase.cpp \
	CCfg/Io/Db/ConfigStructureVersion.cpp \
	CCfg/Io/Db/ConfigVersion.cpp \
	CCfg/Io/Db/Config.cpp \
	CCfg/Io/Db/ConfigAttribute.cpp \
	CCfg/Io/Db/ConfigClass.cpp \
	CCfg/Io/Db/ConfigClassInfo.cpp \
	CCfg/Io/Db/Site.cpp \
	CCfg/Io/Db/Run.cpp \
	CCfg/Io/Db/RunCondition.cpp \
	CCfg/Io/Db/InsertionProcess.cpp \
	CCfg/Io/Db/ConfigDao.cpp \
	CCfg/Io/Db/Transaction.cpp \
	CCfg/Io/Db/TableBuilder.cpp \
	CCfg/Io/Db/MySqlTableBuilder.cpp \
	CCfg/Io/Db/AbstractEntityBatch.cpp \
	CCfg/Io/Db/ConfigBatch.cpp \
	CCfg/Io/Db/ConfigStructureVersionBatch.cpp \
	CCfg/Io/Db/ConfigVersionBatch.cpp \
	CCfg/Io/Db/ConfigClassBatch.cpp \
	CCfg/Io/Db/ConfigClassInfoBatch.cpp \
	CCfg/Io/Db/ConfigAttributeBatch.cpp \
	CCfg/Io/Db/SiteBatch.cpp \
	CCfg/Io/Db/RunBatch.cpp \
	CCfg/Io/Db/RunConditionBatch.cpp \
	CCfg/Io/Db/AbstractUpdateController.cpp \
	CCfg/Io/Db/Exception.cpp \
	CCfg/Io/Db/PrimaryKeyGenerator.cpp \
	CCfg/Io/Db/ExtractionProcess.cpp \
	CCfg/Io/DatabaseFormat.cpp \
	CCfg/Io/DatabaseStorage.cpp \
	CCfg/Io/AbstractLoginCallback.cpp \
	CCfg/Io/DatabaseFormatFactory.cpp \
	CCfg/Io/DatabaseStorageFactory.cpp \
	Utils/Date/StdTmHelpers.cpp

if CCFG_USE_ORACLE
headers += 	CCfg/Io/Db/OracleTableBuilder.h \
	CCfg/Io/Db/OracleConfigDao.h
sources += 	CCfg/Io/Db/OracleTableBuilder.cpp \
	CCfg/Io/Db/OracleConfigDao.cpp
endif

if CCFG_USE_POSTGRESQL
headers += CCfg/Io/Db/PostgreSqlTableBuilder.h \
	CCfg/Io/Db/PostgreSqlConfigDao.h
sources += CCfg/Io/Db/PostgreSqlTableBuilder.cpp \
	CCfg/Io/Db/PostgreSqlConfigDao.cpp
endif

if CCFG_USE_SQLITE3
headers += CCfg/Io/Db/Sqlite3TableBuilder.h \
	CCfg/Io/Db/Sqlite3ConfigDao.h
sources += CCfg/Io/Db/Sqlite3TableBuilder.cpp \
	CCfg/Io/Db/Sqlite3ConfigDao.cpp
endif

library_includedir = $(includedir)
nobase_library_include_HEADERS = $(headers)

AM_CPPFLAGS = $(CCFG_BACKENDS_FLAGS) $(COMPOUND_CONFIG_CFLAGS) $(SOCI_CPPFLAGS)
AM_CPPFLAGS += $(MYSQL_CFLAGS) $(ORACLE_OCI_CFLAGS) $(POSTGRESQL_CFLAGS) $(SQLITE3_CFLAGS)

lib_LTLIBRARIES = lib@CONFIG_DATABASE_LIB_NAME@.la
lib@CONFIG_DATABASE_LIB_NAME@_la_SOURCES = $(sources)
lib@CONFIG_DATABASE_LIB_NAME@_la_LDFLAGS = -version-info 0:0:0 \
	-no-undefined $(COMPOUND_CONFIG_LIBS) $(SOCI_LDFLAGS) \
	$(MYSQL_LDFLAGS) $(ORACLE_OCI_LDFLAGS) $(POSTGRESQL_LDFLAGS) $(SQLITE3_LDFLAGS)
