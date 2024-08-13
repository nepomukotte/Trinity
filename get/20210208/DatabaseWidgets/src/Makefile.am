include $(top_srcdir)/am/autotroll.am
include $(top_srcdir)/am/cppcheck.am

AUTOMAKE_OPTIONS = subdir-objects

headers = CCfg/Model/DatabaseModel.h \
	CCfg/Model/DatabaseFilterProxyModel.h \
	CCfg/Model/ConnectionRootItem.h \
	CCfg/Model/AbstractDatabaseItem.h \
	CCfg/Model/ConnectionItem.h \
	CCfg/Model/SchemaItem.h \
	CCfg/Model/ConfigItem.h \
	CCfg/Model/VersionItem.h \
	CCfg/Model/SiteItem.h \
	CCfg/Model/RunItem.h \
	CCfg/Model/ConditionItem.h \
	CCfg/Model/ConfigModel2.h \
	CCfg/Ui/DbConfigChooserDialog.h \
	CCfg/Ui/DbRunEditorDialog.h \
	CCfg/Ui/DatabaseBrowserTree.h \
	CCfg/Ui/EditorMainWindow2.h \
	CCfg/Ui/LoginDialog.h \
	CCfg/Ui/InteractiveLoginCallback.h \
	CCfg/Ui/UpdateControllerDialog.h \
	CCfg/Ui/InteractiveUpdateController.h \
	CCfg/Ui/ConnectionPropertiesDialog.h \
	CCfg/Ui/SchemaPropertiesDialog.h \
	CCfg/Ui/UserAccessDialog.h \
	CCfg/Control/Workspace2.h

sources = CCfg/Model/DatabaseModel.cpp \
	CCfg/Model/DatabaseFilterProxyModel.cpp \
	CCfg/Model/ConnectionRootItem.cpp \
	CCfg/Model/AbstractDatabaseItem.cpp \
	CCfg/Model/ConnectionItem.cpp \
	CCfg/Model/SchemaItem.cpp \
	CCfg/Model/ConfigItem.cpp \
	CCfg/Model/VersionItem.cpp \
	CCfg/Model/SiteItem.cpp \
	CCfg/Model/RunItem.cpp \
	CCfg/Model/ConditionItem.cpp \
	CCfg/Model/ConfigModel2.cpp \
	CCfg/Ui/DbConfigChooserDialog.cpp \
	CCfg/Ui/DbRunEditorDialog.cpp \
	CCfg/Ui/DatabaseBrowserTree.cpp \
	CCfg/Ui/EditorMainWindow2.cpp \
	CCfg/Ui/LoginDialog.cpp \
	CCfg/Ui/InteractiveLoginCallback.cpp \
	CCfg/Ui/UpdateControllerDialog.cpp \
	CCfg/Ui/InteractiveUpdateController.cpp \
	CCfg/Ui/ConnectionPropertiesDialog.cpp \
	CCfg/Ui/SchemaPropertiesDialog.cpp \
	CCfg/Ui/UserAccessDialog.cpp \
	CCfg/Control/Workspace2.cpp

moc_sources = CCfg/Ui/DatabaseBrowserTree.moc.cpp \
	CCfg/Ui/DbConfigChooserDialog.moc.cpp \
	CCfg/Ui/DbRunEditorDialog.moc.cpp \
	CCfg/Ui/EditorMainWindow2.moc.cpp \
	CCfg/Ui/ConnectionPropertiesDialog.moc.cpp \
	CCfg/Ui/UserAccessDialog.moc.cpp \
	CCfg/Control/Workspace2.moc.cpp

qrc_files = DatabaseWidgets.qrc
qrc_sources = DatabaseWidgets.qrc.cpp

ui_files = CCfg/Ui/ConnectionPropertiesDialog.ui \
	CCfg/Ui/DbConfigChooserDialog.ui \
	CCfg/Ui/DbRunEditorDialog.ui \
	CCfg/Ui/LoginDialog.ui \
	CCfg/Ui/SchemaPropertiesDialog.ui \
	CCfg/Ui/UpdateControllerDialog.ui \
	CCfg/Ui/UserAccessDialog.ui

ui_headers = CCfg/Ui/ConnectionPropertiesDialog.ui.h \
	CCfg/Ui/DbConfigChooserDialog.ui.h \
	CCfg/Ui/DbRunEditorDialog.ui.h \
	CCfg/Ui/LoginDialog.ui.h \
	CCfg/Ui/SchemaPropertiesDialog.ui.h \
	CCfg/Ui/UpdateControllerDialog.ui.h \
	CCfg/Ui/UserAccessDialog.ui.h

gen_sources = $(moc_sources) $(qrc_sources) $(ui_headers)

library_includedir = $(includedir)
nobase_library_include_HEADERS = $(headers)

AM_CPPFLAGS = -I@top_srcdir@/src $(QT_CPPFLAGS) $(CONFIG_DB_AND_WIDGETS_CFLAGS)
AM_CXXFLAGS = $(QT_CXXFLAGS)

lib_LTLIBRARIES = lib@DATABASE_WIDGETS_LIB_NAME@.la
lib@DATABASE_WIDGETS_LIB_NAME@_la_SOURCES = $(sources)
nodist_lib@DATABASE_WIDGETS_LIB_NAME@_la_SOURCES = $(gen_sources)
lib@DATABASE_WIDGETS_LIB_NAME@_la_LDFLAGS = -version-info 0:0:0 $(QT_LDFLAGS) -no-undefined
lib@DATABASE_WIDGETS_LIB_NAME@_la_LIBADD = $(CONFIG_DB_AND_WIDGETS_LIBS) $(QT_LIBS)

BUILT_SOURCES = $(gen_sources)

EXTRA_DIST = $(qrc_files) $(ui_files)


bin_PROGRAMS = dbconfiged@DEBUG_SUFFIX@
dbconfiged@DEBUG_SUFFIX@_SOURCES = main.cpp
dbconfiged@DEBUG_SUFFIX@_LDADD = lib@DATABASE_WIDGETS_LIB_NAME@.la $(CONFIG_DB_AND_WIDGETS_LIBS) $(QT_LIBS)
