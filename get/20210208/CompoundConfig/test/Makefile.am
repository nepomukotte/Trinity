AUTOMAKE_OPTIONS = subdir-objects

headers = CCfg/Io/Config.h \
	CCfg/View/Configs.h

sources = CCfg/AttributeTest.cpp \
	CCfg/CConfigTest.cpp \
	CCfg/EnumTest.cpp \
	CCfg/IntervalTest.cpp \
	CCfg/Io/Config.cpp \
	CCfg/Io/UrlTest.cpp \
	CCfg/ObjectTest.cpp \
	CCfg/PathNode.cpp \
	CCfg/PathTest.cpp \
	CCfg/RangeTest.cpp \
	CCfg/ScalarTest.cpp \
	CCfg/StructureNodeTest.cpp \
	CCfg/ValidationVisitorTest.cpp \
	CCfg/ValueTest.cpp \
	CCfg/View/ClassIteratorTest.cpp \
	CCfg/View/Configs.cpp \
	CCfg/View/InheritedIteratorTest.cpp \
	CCfg/View/InheritersIteratorTest.cpp \
	CCfg/VisitorTest.cpp \
	main.cpp

test_files = data/daq_ill.xml \
	data/daq_inc_mod.xml \
	data/daq_inc.xml \
	data/daq_mod.xml \
	data/daq_rec.xml \
	data/daq.xml \
	data/daq.bcfg \
	data/daq2_ill.bcfg \
	data/daq2.bcfg \
	data/daq2.ccfg \
	data/default_channel_mod.xml \
	data/default_channel_rec.xml \
	data/channels/default_channel.xml

if BUILD_CCFG_FMT_XML
    sources += \
		CCfg/Io/XmlFormatTest.cpp \
		CCfg/Io/ConfigToXmlVisitorTest.cpp \
		CCfg/Io/XmlToConfigVisitorTest.cpp \
		CCfg/Io/DocumentTest.cpp
endif
	
if BUILD_CCFG_FMT_BINARY
    sources += \
		CCfg/Io/BinaryFormatTest.cpp \
		CCfg/Io/DocumentBinaryTest.cpp
endif

if BUILD_CCFG_FMT_XML
if BUILD_CCFG_FMT_BINARY
    sources += \
		CCfg/Io/XmlBinaryFormatsTest.cpp
endif
endif

AM_CXXFLAGS = $(CXXFLAGS_CCFG)
AM_CPPFLAGS = $(CCFG_ENDIANNESS_FLAGS) -I@srcdir@/../src $(TINYXML_CFLAGS) $(BOOST_CPPFLAGS)

TESTS = CompoundConfigTest
noinst_PROGRAMS = CompoundConfigTest
CompoundConfigTest_SOURCES = $(sources)
CompoundConfigTest_LDADD = @builddir@/../src/lib@CCFG_LIB_NAME@.la $(TINYXML_LIBS)

EXTRA_DIST = $(headers) $(test_files)
DISTCLEANFILES = $(test_files)

check-am: check-before
check-before:
	(cd $(srcdir) && rsync -aR $(test_files) $(abs_builddir) && chmod -R u+w $(abs_builddir)/data)

cppcheck:

include $(srcdir)/Makefile.vxworks.am
