include $(top_srcdir)/am/autotroll.am
include $(top_srcdir)/am/cppcheck.am

AUTOMAKE_OPTIONS = subdir-objects

headers = CCfg/Control/AbstractDelegate.h \
	CCfg/Control/BaseConfigEditorController.h \
	CCfg/Control/ConfigEditorController.h \
	CCfg/Control/Delegate/ConfigEditorDelegate.h \
	CCfg/Control/ReadOnlyConfigEditorController.h \
	CCfg/Control/Workspace.h \
	CCfg/Model/AbstractColumn.h \
	CCfg/Model/AbstractCommand.h \
	CCfg/Model/AbstractConfigModelAccessor.h \
	CCfg/Model/AbstractRow.h \
	CCfg/Model/AttributeBuilderModel.h \
	CCfg/Model/Column/AttributeBacktrack.h \
	CCfg/Model/Column/AttributeInclude.h \
	CCfg/Model/Column/AttributeName.h \
	CCfg/Model/Column/AttributeRange.h \
	CCfg/Model/Column/AttributeUnit.h \
	CCfg/Model/Column/AttributeValue.h \
	CCfg/Model/Column/ConfigTreeColumn.h \
	CCfg/Model/Column/ConfigTreeRow.h \
	CCfg/Model/Command/AbstractConfigModelCommand.h \
	CCfg/Model/Command/Commands.h \
	CCfg/Model/Command/DeleteConfigNode.h \
	CCfg/Model/Command/ImportConfig.h \
	CCfg/Model/Command/IncludeConfig.h \
	CCfg/Model/Command/InsertConfigNode.h \
	CCfg/Model/Command/OverrideConfigNode.h \
	CCfg/Model/Command/RenameConfigNode.h \
	CCfg/Model/Command/SetConfigNodeIncludeUrl.h \
	CCfg/Model/Command/SetConfigNodeRange.h \
	CCfg/Model/Command/SetConfigNodeUnit.h \
	CCfg/Model/Command/SetConfigNodeValue.h \
	CCfg/Model/ConfigModelAccessor.h \
	CCfg/Model/ConfigModel.h \
	CCfg/Model/ConfigNode.h \
	CCfg/Model/ConfigSetModel.h \
	CCfg/Model/EditorConfiguration.h \
	CCfg/Model/UndoableConfigModelAccessor.h \
	CCfg/Ui/AttributeBuilderWidget.h \
	CCfg/Ui/AttributeFinderDialog.h \
	CCfg/Ui/ClosableTabWidget.h \
	CCfg/Ui/CommonParentSelectionModel.h \
	CCfg/Ui/ConfigEditorTree.h \
	CCfg/Ui/ConfigSetWidget.h \
	CCfg/Ui/ConfigSetWorkspaceWidget.h \
	CCfg/Ui/EditorMainWindow.h \
	CCfg/Ui/EditConfigDialog.h \
	CCfg/Ui/Enums.h \
	CCfg/Ui/ErrorGuiPrinter.h \
	CCfg/Ui/MyComboBox.h \
	CCfg/Ui/TypeDisplayInfo.h \
	CCfg/Ui/ValidationReportDialog.h \
	CCfg/Ui/ValueEditionWidget.h \
	Utils/Error/AbstractFormatter.h \
	Utils/Error/AbstractPrinter.h \
	Utils/Error/ExcepHandling.h \
	Utils/Error/Exception.h \
	Utils/Error/Reporter.h \
	Utils/Ui/DesignedWidget.hpp \
	Utils/Ui/fromArgs.h \
	Utils/Ui/FileChooser.h

sources = CCfg/Control/AbstractDelegate.cpp \
	CCfg/Control/BaseConfigEditorController.cpp \
	CCfg/Control/ConfigEditorController.cpp \
	CCfg/Control/Delegate/ConfigEditorDelegate.cpp \
	CCfg/Control/ReadOnlyConfigEditorController.cpp \
	CCfg/Control/Workspace.cpp \
	CCfg/Model/AbstractColumn.cpp \
	CCfg/Model/AbstractCommand.cpp \
	CCfg/Model/AbstractConfigModelAccessor.cpp \
	CCfg/Model/AbstractRow.cpp \
	CCfg/Model/AttributeBuilderModel.cpp \
	CCfg/Model/Column/AttributeBacktrack.cpp \
	CCfg/Model/Column/AttributeInclude.cpp \
	CCfg/Model/Column/AttributeName.cpp \
	CCfg/Model/Column/AttributeRange.cpp \
	CCfg/Model/Column/AttributeUnit.cpp \
	CCfg/Model/Column/AttributeValue.cpp \
	CCfg/Model/Column/ConfigTreeColumn.cpp \
	CCfg/Model/Column/ConfigTreeRow.cpp \
	CCfg/Model/Command/AbstractConfigModelCommand.cpp \
	CCfg/Model/Command/DeleteConfigNode.cpp \
	CCfg/Model/Command/ImportConfig.cpp \
	CCfg/Model/Command/IncludeConfig.cpp \
	CCfg/Model/Command/InsertConfigNode.cpp \
	CCfg/Model/Command/OverrideConfigNode.cpp \
	CCfg/Model/Command/RenameConfigNode.cpp \
	CCfg/Model/Command/SetConfigNodeIncludeUrl.cpp \
	CCfg/Model/Command/SetConfigNodeRange.cpp \
	CCfg/Model/Command/SetConfigNodeUnit.cpp \
	CCfg/Model/Command/SetConfigNodeValue.cpp \
	CCfg/Model/ConfigModelAccessor.cpp \
	CCfg/Model/ConfigModel.cpp \
	CCfg/Model/ConfigNode.cpp \
	CCfg/Model/ConfigSetModel.cpp \
	CCfg/Model/EditorConfiguration.cpp \
	CCfg/Model/UndoableConfigModelAccessor.cpp \
	CCfg/Ui/AttributeBuilderWidget.cpp \
	CCfg/Ui/AttributeFinderDialog.cpp \
	CCfg/Ui/ClosableTabWidget.cpp \
	CCfg/Ui/CommonParentSelectionModel.cpp \
	CCfg/Ui/ConfigEditorTree.cpp \
	CCfg/Ui/ConfigSetWidget.cpp \
	CCfg/Ui/ConfigSetWorkspaceWidget.cpp \
	CCfg/Ui/EditorMainWindow.cpp \
	CCfg/Ui/EditConfigDialog.cpp \
	CCfg/Ui/ErrorGuiPrinter.cpp \
	CCfg/Ui/MyComboBox.cpp \
	CCfg/Ui/TypeDisplayInfo.cpp \
	CCfg/Ui/ValidationReportDialog.cpp \
	CCfg/Ui/ValueEditionWidget.cpp \
	Utils/Error/AbstractFormatter.cpp \
	Utils/Error/AbstractPrinter.cpp \
	Utils/Error/Exception.cpp \
	Utils/Error/Reporter.cpp \
	Utils/Ui/fromArgs.cpp \
	Utils/Ui/FileChooser.cpp

moc_sources = CCfg/Control/BaseConfigEditorController.moc.cpp \
	CCfg/Control/ConfigEditorController.moc.cpp \
	CCfg/Control/ReadOnlyConfigEditorController.moc.cpp \
	CCfg/Control/Workspace.moc.cpp \
	CCfg/Model/AttributeBuilderModel.moc.cpp \
	CCfg/Model/ConfigModelAccessor.moc.cpp \
	CCfg/Model/ConfigModel.moc.cpp \
	CCfg/Model/ConfigSetModel.moc.cpp \
	CCfg/Model/UndoableConfigModelAccessor.moc.cpp \
	CCfg/Ui/AttributeBuilderWidget.moc.cpp \
	CCfg/Ui/AttributeFinderDialog.moc.cpp \
	CCfg/Ui/ClosableTabWidget.moc.cpp \
	CCfg/Ui/ConfigEditorTree.moc.cpp \
	CCfg/Ui/ConfigSetWidget.moc.cpp \
	CCfg/Ui/ConfigSetWorkspaceWidget.moc.cpp \
	CCfg/Ui/EditorMainWindow.moc.cpp \
	CCfg/Ui/EditConfigDialog.moc.cpp \
	CCfg/Ui/MyComboBox.moc.cpp \
	CCfg/Ui/ValidationReportDialog.moc.cpp \
	CCfg/Ui/ValueEditionWidget.moc.cpp \
	Utils/Ui/FileChooser.moc.cpp

qrc_files = ConfigWidgets.qrc
qrc_sources = ConfigWidgets.qrc.cpp

ui_files = CCfg/Ui/AttributeBuilder.ui \
	CCfg/Ui/AttributeFinder.ui \
	CCfg/Ui/ConfigSetWidget.ui \
	CCfg/Ui/ConfigSetWorkspaceWidget.ui \
	CCfg/Ui/EditConfigDialog.ui \
	CCfg/Ui/MainWindow.ui \
	CCfg/Ui/ValidationReportDialog.ui

ui_headers = CCfg/Ui/AttributeBuilder.ui.h \
	CCfg/Ui/AttributeFinder.ui.h \
	CCfg/Ui/ConfigSetWidget.ui.h \
	CCfg/Ui/ConfigSetWorkspaceWidget.ui.h \
	CCfg/Ui/EditConfigDialog.ui.h \
	CCfg/Ui/MainWindow.ui.h \
	CCfg/Ui/ValidationReportDialog.ui.h

gen_sources = $(moc_sources) $(qrc_sources) $(ui_headers)

library_includedir = $(includedir)
nobase_library_include_HEADERS = $(headers)

AM_CPPFLAGS = -I@srcdir@ $(QT_CPPFLAGS) $(COMPOUND_CONFIG_CFLAGS) $(UTILS_CFLAGS)
AM_CXXFLAGS = $(CXXFLAGS_CCFG_WIDGETS) -Wno-long-long $(QT_CXXFLAGS)

lib_LTLIBRARIES = lib@CONFIG_WIDGETS_LIB_NAME@.la
lib@CONFIG_WIDGETS_LIB_NAME@_la_SOURCES = $(sources)
nodist_lib@CONFIG_WIDGETS_LIB_NAME@_la_SOURCES = $(gen_sources)
lib@CONFIG_WIDGETS_LIB_NAME@_la_LDFLAGS = -no-undefined -version-info 1:0:0 $(QT_LDFLAGS)
lib@CONFIG_WIDGETS_LIB_NAME@_la_LIBADD = $(COMPOUND_CONFIG_LIBS) $(QT_LIBS) $(UTILS_LIBS)

BUILT_SOURCES = $(gen_sources)

EXTRA_DIST = $(qrc_files) $(ui_files)


bin_PROGRAMS = cconfiged@DEBUG_SUFFIX@
cconfiged@DEBUG_SUFFIX@_SOURCES = main.cpp
cconfiged@DEBUG_SUFFIX@_LDFLAGS = $(QT_LDFLAGS)
cconfiged@DEBUG_SUFFIX@_LDADD = lib@CONFIG_WIDGETS_LIB_NAME@.la $(QT_LIBS) $(COMPOUND_CONFIG_LIBS)

