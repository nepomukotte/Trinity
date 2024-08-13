include $(top_srcdir)/am/cppcheck.am

AUTOMAKE_OPTIONS = subdir-objects

headers = CCfg/Attribute.h \
	CCfg/CConfig.h \
	CCfg/ClassInfo.h \
	CCfg/Comment.h \
	CCfg/Enum.h \
	CCfg/Exception.h \
	CCfg/Interval.h \
	CCfg/Object.h \
	CCfg/Path.h \
	CCfg/PathNode.h \
	CCfg/Range.h \
	CCfg/ScalarClassInfo.h \
	CCfg/Scalar.h \
	CCfg/StructureNode.h \
	CCfg/TypeInfo.h \
	CCfg/UnsupportedOperation.h \
	CCfg/ValidationVisitor.h \
	CCfg/Value.h \
	CCfg/Visitor.h \
	CCfg/Io/AbstractFormatFactory.h \
	CCfg/Io/AbstractStorageFactory.h \
	CCfg/Io/BufferStorage.h \
	CCfg/Io/ConfigToTextVisitor.h \
	CCfg/Io/Document.h \
	CCfg/Io/FileStorageFactory.h \
	CCfg/Io/FileStorage.h \
	CCfg/Io/FormatFactory.h \
	CCfg/Io/Format.h \
	CCfg/Io/IncludeHandler.h \
	CCfg/Io/StorageFactory.h \
	CCfg/Io/Storage.h \
	CCfg/Io/StreamStorage.h \
	CCfg/Io/TextFormat.h \
	CCfg/Io/UnhandledUrl.h \
	CCfg/Io/Url.h \
	CCfg/View/ClassIterator.h \
	CCfg/View/InheritanceIterator.h \
	CCfg/View/InheritedIterator.h \
	CCfg/View/InheritersIterator.h \
	CCfg/View/Result.h \
	CCfg/View/Search.h \
	CCfg/View/Search.hpp \
	Utils/Io/PointerValueStreamer.hpp \
	Utils/Math/float.h \
	Utils/String/Comparison.hpp

noinst_headers = CCfg/Grammar/EnumGrammar.h \
	CCfg/Grammar/EnumGrammar.hpp \
	CCfg/Grammar/IntervalGrammar.h \
	CCfg/Grammar/IntervalGrammar.hpp \
	CCfg/Grammar/RangeGrammar.h \
	CCfg/Grammar/RangeGrammar.hpp \
	CCfg/Grammar/ValueGrammar.h \
	CCfg/Grammar/ValueGrammar.hpp

sources = CCfg/Attribute.cpp \
	CCfg/CConfig.cpp \
	CCfg/ClassInfo.cpp \
	CCfg/Comment.cpp \
	CCfg/Enum.cpp \
	CCfg/Exception.cpp \
	CCfg/Interval.cpp \
	CCfg/Object.cpp \
	CCfg/Path.cpp \
	CCfg/PathNode.cpp \
	CCfg/Range.cpp \
	CCfg/ScalarClassInfo.cpp \
	CCfg/Scalar.cpp \
	CCfg/StructureNode.cpp \
	CCfg/TypeInfo.cpp \
	CCfg/UnsupportedOperation.cpp \
	CCfg/ValidationVisitor.cpp \
	CCfg/Value.cpp \
	CCfg/Visitor.cpp \
	CCfg/Grammar/EnumGrammar.cpp \
	CCfg/Grammar/IntervalGrammar.cpp \
	CCfg/Grammar/RangeGrammar.cpp \
	CCfg/Grammar/ValueGrammar.cpp \
	CCfg/Io/AbstractFormatFactory.cpp \
	CCfg/Io/AbstractStorageFactory.cpp \
	CCfg/Io/BufferStorage.cpp \
	CCfg/Io/ConfigToTextVisitor.cpp \
	CCfg/Io/Document.cpp \
	CCfg/Io/FileStorage.cpp \
	CCfg/Io/FileStorageFactory.cpp \
	CCfg/Io/Format.cpp \
	CCfg/Io/FormatFactory.cpp \
	CCfg/Io/IncludeHandler.cpp \
	CCfg/Io/Storage.cpp \
	CCfg/Io/StorageFactory.cpp \
	CCfg/Io/StreamStorage.cpp \
	CCfg/Io/TextFormat.cpp \
	CCfg/Io/UnhandledUrl.cpp \
	CCfg/Io/Url.cpp \
	CCfg/View/ClassIterator.cpp \
	CCfg/View/InheritanceIterator.cpp \
	CCfg/View/InheritedIterator.cpp \
	CCfg/View/InheritersIterator.cpp \
	CCfg/View/Result.cpp \
	CCfg/View/Search.cpp \
	Utils/Math/float.cpp

AM_CXXFLAGS = $(CXXFLAGS_CCFG)
AM_CPPFLAGS = $(CCFG_ENDIANNESS_FLAGS) $(BOOST_CPPFLAGS)

if BUILD_CCFG_FMT_XML
    headers += \
		CCfg/Io/XmlFormatFactory.h \
		CCfg/Io/XmlFormat.h \
		CCfg/Io/XmlToConfigVisitor.h \
		CCfg/Io/ConfigToXmlVisitor.h
    sources += \
		CCfg/Io/XmlFormat.cpp \
		CCfg/Io/XmlFormatFactory.cpp \
		CCfg/Io/XmlToConfigVisitor.cpp \
		CCfg/Io/ConfigToXmlVisitor.cpp
    AM_CPPFLAGS += $(TINYXML_CFLAGS) -DBUILD_CCFG_FMT_XML=1
endif

if BUILD_CCFG_FMT_BINARY
    headers += \
		CCfg/Io/BinaryFormat.h \
		CCfg/Io/BinaryFormatFactory.h
    noinst_headers += \
		CCfg/ValueRange.h \
		CCfg/AttributeClass.h
    sources += \
		CCfg/ValueRange.cpp \
		CCfg/AttributeClass.cpp \
		CCfg/Io/BinaryFormat.cpp \
		CCfg/Io/BinaryFormatFactory.cpp
    AM_CPPFLAGS += -DBUILD_CCFG_FMT_BINARY=1
endif

library_includedir = $(includedir)
nobase_library_include_HEADERS = $(headers)
noinst_HEADERS = $(noinst_headers)

lib_LTLIBRARIES = lib@CCFG_LIB_NAME@.la
lib@CCFG_LIB_NAME@_la_SOURCES = $(sources)
lib@CCFG_LIB_NAME@_la_LDFLAGS = -version-info 1:0:0 -no-undefined
lib@CCFG_LIB_NAME@_la_LIBADD = $(TINYXML_LIBS)
