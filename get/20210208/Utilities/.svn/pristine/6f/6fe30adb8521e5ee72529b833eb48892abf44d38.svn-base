include $(top_srcdir)/am/buildinfo.am
include $(top_srcdir)/am/cppcheck.am

AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4

AM_CXXFLAGS = $(CXXFLAGS_UTILS)

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = @UTILS_LIB_NAME@.pc

dist_noinst_SCRIPTS = autogen.sh cleanup.sh

headers =                           \
	src/utl/BinIO.hpp               \
	src/utl/BitFieldHelper.hpp      \
	src/utl/BitFieldHelperBitSet.hpp\
	src/utl/BitFieldHelperVector.hpp\
	src/utl/WordIterator.hpp        \
	src/utl/RtBitSet.hpp            \
	src/utl/BitSet.hpp              \
	src/utl/BitFieldRegister.hpp    \
	src/utl/BuildInfo.h             \
	src/utl/Configurator.h          \
	src/utl/Construct.hpp           \
	src/utl/EnumStringMapper.hpp    \
	src/utl/EnumStringMapper.ipp    \
	src/utl/Exception.h             \
	src/utl/Factory.hpp             \
	src/utl/Factory.ipp             \
	src/utl/ForwardCtor.hpp	        \
	src/utl/Logging.h               \
	src/utl/MemOperator.hpp         \
	src/utl/Rounding.hpp            \
	src/utl/Singleton.hpp           \
	src/utl/Singleton.ipp

compat_headers =						\
	src/utl/compat/ArrayHelper.hpp

io_headers =                        \
	src/utl/io/Buffer.hpp           \
	src/utl/io/InplaceParser.hpp    \
	src/utl/io/InplaceParser.ipp

net_headers =                       \
	src/utl/net/IpAddress.h         \
	src/utl/net/SocketAddress.h

stack_headers =                     \
	src/utl/stack/LineResolver.h    \
	src/utl/stack/LineResolver.ipp  \
	src/utl/stack/StackFrame.h      \
	src/utl/stack/StackFrame.ipp    \
	src/utl/stack/StackWalker.h     \
	src/utl/stack/SymbolResolver.h  \
	src/utl/stack/SymbolResolver.ipp

sources =                           \
	src/utl/BinIO.cpp               \
	src/utl/BuildInfo.cpp           \
	src/utl/Configurator.cpp        \
	src/utl/EnumStringMapper.cpp    \
	src/utl/Exception.cpp           \
	src/utl/Logging.cpp             \
	src/utl/net/IpAddress.cpp       \
	src/utl/net/SocketAddress.cpp   \
	src/utl/stack/StackWalker.cpp

test_sources =                      \
	test/main.cpp                   \
	test/utl/BitFieldHelperTest.cpp \
	test/utl/ConfiguratorTest.cpp   \
	test/utl/FactoryTest.cpp        \
	test/utl/ForwardCtorTest.cpp    \
	test/utl/RoundingTest.cpp       \
	test/utl/BitSetTest.cpp
	
utl_includedir = $(includedir)/utl
utl_include_HEADERS = $(headers) $(top_builddir)/config.h

utl_compat_includedir = $(utl_includedir)/compat
utl_compat_include_HEADERS = $(compat_headers)

utl_io_includedir = $(utl_includedir)/io
utl_io_include_HEADERS = $(io_headers)

utl_net_includedir = $(utl_includedir)/net
utl_net_include_HEADERS = $(net_headers)

utl_stack_includedir = $(utl_includedir)/stack
utl_stack_include_HEADERS = $(stack_headers)

lib_LTLIBRARIES = lib@UTILS_LIB_NAME@.la
lib@UTILS_LIB_NAME@_la_SOURCES = $(sources)
lib@UTILS_LIB_NAME@_la_LDFLAGS = -version-info $(UTILS_LIB_VERSION) -no-undefined

if WINDOWS
lib@UTILS_LIB_NAME@_la_LDFLAGS += -ldbghelp -lws2_32
endif

TESTS = UtilitiesTest BitSetStandaloneTest
noinst_PROGRAMS = UtilitiesTest
UtilitiesTest_SOURCES = $(test_sources)
UtilitiesTest_LDADD = @builddir@/lib@UTILS_LIB_NAME@.la

noinst_PROGRAMS += BitSetStandaloneTest 
BitSetStandaloneTest_SOURCES = test/utl/BitSetStandaloneTest.cpp
BitSetStandaloneTest_LDADD = @builddir@/lib@UTILS_LIB_NAME@.la

AM_CPPFLAGS = $(BOOST_CPPFLAGS) -I@srcdir@/src

BUILT_SOURCES = genbuildinfo.h
CLEANFILES = genbuildinfo.h

include Makefile.VxWorks.am

include $(top_srcdir)/am/pkg-linux.am
include $(top_srcdir)/am/pkg-win.am
