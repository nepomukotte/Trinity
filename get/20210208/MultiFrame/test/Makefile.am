AUTOMAKE_OPTIONS = subdir-objects

headers = \
	tut.h \
	tut_reporter.h

sources = \
	mfm/PrimaryHeaderTest.cpp \
	mfm/BlobHeaderTest.cpp \
	mfm/StandardHeaderTest.cpp \
	mfm/FrameTest.cpp \
	mfm/FieldTest.cpp \
	mfm/BitFieldTest.cpp \
	mfm/FrameFormatTest.cpp \
	mfm/FrameDictionaryTest.cpp \
	mfm/SignedTest.cpp \
	mfm/BlobFrameTest.cpp \
	mfm/LayeredFrameTest.cpp \
	mfm/FrameBuilderTest.cpp \
	main.cpp
	
test_files = \
	mfm/sample_formats.xml

AM_CPPFLAGS = -I@top_srcdir@/src $(UTILS_CFLAGS) $(CCFG_CFLAGS) $(BOOST_CPPFLAGS)

TESTS = MultiFrameTest
noinst_PROGRAMS = MultiFrameTest
MultiFrameTest_SOURCES = $(sources)
MultiFrameTest_LDADD = @top_builddir@/src/lib@MFM_LIB_NAME@.la

EXTRA_DIST = $(headers) $(test_files)
CLEANFILES = $(test_files)

cppcheck:

check-am: check-before
check-before:
	(cd $(srcdir) && rsync -R $(test_files) $(abs_builddir))

include Makefile.VxWorks.am
