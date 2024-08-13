AUTOMAKE_OPTIONS = subdir-objects
include $(top_srcdir)/am/cppcheck.am

headers =				         \
	mfm/Common.h		         \
	mfm/Buffer.h				\
	mfm/Exception.h              \
	mfm/PrimaryHeader.h	         \
	mfm/Header.h		         \
	mfm/BlobHeader.h 	         \
	mfm/StandardHeader.h         \
	mfm/Serializer.h	         \
	mfm/Frame.h			         \
	mfm/Item.h 			         \
	mfm/AbstractFieldContainer.h \
	mfm/Field.h			         \
	mfm/Field.hpp		         \
	mfm/FrameFormat.h            \
	mfm/FrameDictionary.h        \
	mfm/BitField.h               \
	mfm/BitField.hpp             \
	mfm/FrameBuilder.h           \
	mfm/ItemIterator.h

sources =				           \
	mfm/Buffer.cpp					\
    mfm/Exception.cpp              \
	mfm/PrimaryHeader.cpp	       \
	mfm/Header.cpp			       \
	mfm/BlobHeader.cpp 		       \
	mfm/StandardHeader.cpp	       \
	mfm/Serializer.cpp 		       \
	mfm/Frame.cpp			       \
	mfm/Item.cpp			       \
	mfm/Field.cpp			       \
	mfm/FrameFormat.cpp  	       \
	mfm/FrameDictionary.cpp        \
	mfm/BitField.cpp               \
	mfm/FrameBuilder.cpp           \
	mfm/ItemIterator.cpp

AM_CPPFLAGS = -I@top_srcdir@/src $(UTILS_CFLAGS) $(CCFG_CFLAGS) $(BOOST_CPPFLAGS)

library_includedir = $(includedir)
nobase_library_include_HEADERS = $(headers)

lib_LTLIBRARIES = lib@MFM_LIB_NAME@.la
lib@MFM_LIB_NAME@_la_SOURCES = $(sources)
lib@MFM_LIB_NAME@_la_LDFLAGS = -version-info $(LIBRARY_VERSION) -no-undefined
lib@MFM_LIB_NAME@_la_LIBADD = $(UTILS_LIBS) $(CCFG_LIBS)
