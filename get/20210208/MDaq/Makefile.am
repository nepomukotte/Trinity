include $(top_srcdir)/am/doxygen.am
include $(top_srcdir)/am/ice.am

ACLOCAL_AMFLAGS = -I m4
AUTOMAKE_OPTIONS = subdir-objects

dist_noinst_SCRIPTS = autogen.sh cleanup.sh

# Initialisation
BUILT_SOURCES	= 
CLEANFILES		= 
EXTRA_DIST      =
lib_LTLIBRARIES = 
noinst_LTLIBRARIES = 
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = 

# slice2cpp options
SPPFLAGS      = -I@top_srcdir@/src
slicedir = $(datadir)/slice

# flags
AM_CPPFLAGS =                       \
	-I@top_srcdir@/src              \
	-I@top_srcdir@/demo             \
	-I@top_builddir@/src            \
	$(UTILS_CFLAGS)                 \
	$(COMPOUND_CONFIG_CFLAGS)       \
	$(ICE_CPPFLAGS) $(ICE_DEFINE)   \
	$(ICESTORM_CPPFLAGS)            \
	$(LOG4CXX_CPPFLAGS)             \
	$(READLINE_CPPFLAGS)            \
	$(UTL_LOG_THREADSAFE_CPPFLAGS)

AM_CXXFLAGS = $(PTHREAD_CFLAGS)
AM_LDFLAGS = $(ICE_LDFLAGS) $(READLINE_LDFLAGS)
LIBS += $(UTL_LOG_THREADSAFE_LIBS)

####################################################################################################
# 'mdaq-interfaces' library
####################################################################################################
ice_mdaq_sources     =      \
	src/mdaq/Ecc.ice        \
	src/mdaq/DaqControl.ice \
	src/mdaq/hw/HardwareTypes.ice

ice_mdaq_gen_sources =      \
	src/mdaq/Ecc.cpp        \
	src/mdaq/DaqControl.cpp \
	src/mdaq/hw/HardwareTypes.cpp

ice_mdaq_gen_headers =      \
	src/mdaq/Ecc.h          \
	src/mdaq/DaqControl.h   \
	src/mdaq/hw/HardwareTypes.h

# Library
pkgconfig_DATA += mdaq-interfaces@DEBUG_SUFFIX@.pc
lib_LTLIBRARIES += libmdaq-interfaces@DEBUG_SUFFIX@.la
nodist_libmdaq_interfaces@DEBUG_SUFFIX@_la_SOURCES = $(ice_mdaq_gen_sources)
libmdaq_interfaces@DEBUG_SUFFIX@_la_LIBADD  = $(ICE_LIBS)

# Installation of Slice interface
mdaq_slicedir = $(slicedir)/mdaq
mdaq_slice_DATA = $(ice_mdaq_sources)

mdaq_hw_slicedir = $(datadir)/slice/mdaq/hw
mdaq_hw_slice_DATA = src/mdaq/hw/HardwareTypes.ice

# Installation of header files
mdaq_hw_includedir = $(mdaq_includedir)/hw
mdaq_hw_include_HEADERS = src/mdaq/hw/HardwareTypes.h

mdaq_includedir = $(includedir)/mdaq
mdaq_include_HEADERS =         \
	$(ice_mdaq_gen_headers)    \
	src/mdaq/DefaultPortNums.h \
	src/mdaq/EccCmd.h          \
	src/mdaq/EccCmd.hpp        \
	src/mdaq/EccBackEnd.h      \
	src/mdaq/EccImpl.h

BUILT_SOURCES	+= $(ice_mdaq_gen_headers) $(ice_mdaq_gen_sources)
EXTRA_DIST      += $(ice_mdaq_sources)

####################################################################################################
# 'mdaq-utl' library
####################################################################################################
ice_mdaq_utl_sources  =      \
	src/mdaq/utl/LoggingService.ice \
	src/mdaq/utl/CmdException.ice

ice_mdaq_utl_gen_sources =   \
	src/mdaq/utl/LoggingService.cpp \
	src/mdaq/utl/CmdException.cpp

ice_mdaq_utl_gen_headers =   \
	src/mdaq/utl/LoggingService.h \
	src/mdaq/utl/CmdException.h

# Installation of Slice interface
mdaq_utl_slicedir = $(datadir)/slice/mdaq/utl
mdaq_utl_slice_DATA = $(ice_mdaq_utl_sources)

# Library
mdaq_utl_headers =                            \
	$(ice_mdaq_utl_gen_headers)               \
	src/mdaq/utl/CmdLineArgs.h                \
	src/mdaq/utl/Server.h                     \
	src/mdaq/utl/ConsoleLoggingBackend.h      \
	src/mdaq/utl/VersaPointer.hpp

mdaq_utl_sources =                            \
	src/mdaq/utl/CmdLineArgs.cpp              \
	src/mdaq/utl/Server.cpp                   \
	src/mdaq/utl/ConsoleLoggingBackend.cpp	  \
	src/mdaq/utl/net/IpDetect.cpp             \
	src/mdaq/utl/net/Socket.cpp               \
	src/mdaq/utl/net/TcpAcceptor.cpp          \
	src/mdaq/utl/net/TcpConnection.cpp        \
	src/mdaq/utl/net/UdpTransceiver.cpp

if HAVE_ICESTORM
mdaq_utl_headers += \
	src/mdaq/utl/IceStormLoggingBackend.h  \
	src/mdaq/utl/IceStormLoggingSubscriber.h
mdaq_utl_sources += \
	src/mdaq/utl/IceStormLoggingBackend.cpp  \
	src/mdaq/utl/IceStormLoggingSubscriber.cpp
endif

if HAVE_LOG4CXX
mdaq_utl_headers += \
	src/mdaq/utl/Log4cxxLoggingBackend.h
mdaq_utl_sources += \
	src/mdaq/utl/Log4cxxLoggingBackend.cpp
endif

if !VXWORKS
mdaq_utl_sources +=                           \
	src/mdaq/utl/net/SocketAddress.cpp
endif

if VXWORKS
mdaq_utl_sources += src/mdaq/utl/net/ZBufConnection.cpp
endif

pkgconfig_DATA += mdaq-utl@DEBUG_SUFFIX@.pc

lib_LTLIBRARIES += libmdaq-utl@DEBUG_SUFFIX@.la
nodist_libmdaq_utl@DEBUG_SUFFIX@_la_SOURCES = $(ice_mdaq_utl_gen_sources)
libmdaq_utl@DEBUG_SUFFIX@_la_SOURCES = $(mdaq_utl_sources)
libmdaq_utl@DEBUG_SUFFIX@_la_LDFLAGS = $(LOG4CXX_LDFLAGS)
libmdaq_utl@DEBUG_SUFFIX@_la_LIBADD  = $(ICE_LIBS) $(ICESTORM_LIBS)
libmdaq_utl@DEBUG_SUFFIX@_la_LIBADD += $(LOG4CXX_LIBS)
libmdaq_utl@DEBUG_SUFFIX@_la_LIBADD += $(COMPOUND_CONFIG_LIBS) $(UTILS_LIBS)

# Installation of header files
mdaq_utl_includedir  = $(includedir)/mdaq/utl
mdaq_utl_include_HEADERS = $(mdaq_utl_headers)

mdaq_utl_net_includedir  = $(includedir)/mdaq/utl/net
mdaq_utl_net_include_HEADERS = 			\
	src/mdaq/utl/net/Buffer.hpp			\
	src/mdaq/utl/net/IpDetect.h			\
	src/mdaq/utl/net/Socket.h			\
	src/mdaq/utl/net/TcpAcceptor.h		\
	src/mdaq/utl/net/TcpConnection.h    \
	src/mdaq/utl/net/UdpTransceiver.h


if !VXWORKS
mdaq_utl_net_include_HEADERS += src/mdaq/utl/net/SocketAddress.h
else
mdaq_utl_net_include_HEADERS += src/mdaq/utl/net/ZBufConnection.h
endif

BUILT_SOURCES	+= $(ice_mdaq_utl_gen_headers) $(ice_mdaq_utl_gen_sources)
EXTRA_DIST      += $(ice_mdaq_utl_sources)

####################################################################################################
# 'mdaq-cli' library
####################################################################################################

mdaq_cli_sources = \
	src/mdaq/cli/UserQuitting.cpp

mdaq_cli_headers = \
	src/mdaq/cli/ClientCmd.hpp \
	src/mdaq/cli/UserQuitting.h
if ENABLE_CLI
if !ICEE
mdaq_cli_headers += src/mdaq/cli/Cmd.h
mdaq_cli_sources += src/mdaq/cli/Cmd.cpp
endif #!ICEE
endif #ENABLE_CLI

pkgconfig_DATA += mdaq-cli@DEBUG_SUFFIX@.pc

lib_LTLIBRARIES += libmdaq-cli@DEBUG_SUFFIX@.la
nodist_libmdaq_cli@DEBUG_SUFFIX@_la_SOURCES = 
libmdaq_cli@DEBUG_SUFFIX@_la_SOURCES = $(mdaq_cli_sources)
libmdaq_cli@DEBUG_SUFFIX@_la_LDFLAGS = $(LOG4CXX_LDFLAGS)
libmdaq_cli@DEBUG_SUFFIX@_la_LIBADD   = libmdaq-utl@DEBUG_SUFFIX@.la
libmdaq_cli@DEBUG_SUFFIX@_la_LIBADD  += $(ICE_LIBS) $(ICESTORM_LIBS)
libmdaq_cli@DEBUG_SUFFIX@_la_LIBADD  += $(LOG4CXX_LIBS)
libmdaq_cli@DEBUG_SUFFIX@_la_LIBADD  += $(READLINE_LIB) $(CURSES_LIB)
libmdaq_cli@DEBUG_SUFFIX@_la_LIBADD  += $(COMPOUND_CONFIG_LIBS) $(UTILS_LIBS)

# Installation of header files
mdaq_cli_includedir  = $(includedir)/mdaq/cli
mdaq_cli_include_HEADERS = $(mdaq_cli_headers)

####################################################################################################
# 'mdaq-mt' convenience library (part of 'mdaq' library)
####################################################################################################
mdaq_mt_sources = \
	src/mdaq/mt/EventConsumer.cpp            \
	src/mdaq/mt/EventProducer.cpp            \
	src/mdaq/mt/EventMonitor.cpp             \
	src/mdaq/mt/PollingProducer.cpp          \
	src/mdaq/mt/InterruptProducer.cpp        \
	src/mdaq/mt/DefaultForTemplate.cpp       \
	src/mdaq/mt/EventConsumerAdapter.cpp     \
	src/mdaq/mt/PollingProducerAdapter.cpp   \
	src/mdaq/mt/InterruptProducerAdapter.cpp \
	src/mdaq/mt/PollingMonitorAdapter.cpp    \
	src/mdaq/mt/InterruptMonitorAdapter.cpp  \
	src/mdaq/mt/Semaphore.cpp                \
	src/mdaq/mt/SysClock.cpp
if !APPLE
mdaq_mt_sources += \
	src/mdaq/mt/Timer.cpp \
	src/mdaq/mt/SigAction.cpp
endif

pkgconfig_DATA += mdaq-mt@DEBUG_SUFFIX@.pc

lib_LTLIBRARIES += libmdaq-mt@DEBUG_SUFFIX@.la
nodist_libmdaq_mt@DEBUG_SUFFIX@_la_SOURCES = 
libmdaq_mt@DEBUG_SUFFIX@_la_SOURCES = $(mdaq_mt_sources)
libmdaq_mt@DEBUG_SUFFIX@_la_LDFLAGS = $(BOOST_LDFLAGS)
libmdaq_mt@DEBUG_SUFFIX@_la_LIBADD   = 

mdaq_mt_includedir = $(mdaq_includedir)/mt
mdaq_mt_include_HEADERS =			\
	src/mdaq/mt/EventConsumer.h     \
	src/mdaq/mt/EventProducer.h     \
	src/mdaq/mt/PollingProducer.h   \
	src/mdaq/mt/InterruptProducer.h \
	src/mdaq/mt/EventMonitor.h      \
	src/mdaq/mt/DefaultForTemplate.h   \
	src/mdaq/mt/EventConsumerAdapter.hpp \
	src/mdaq/mt/PollingProducerAdapter.hpp \
	src/mdaq/mt/InterruptProducerAdapter.hpp \
	src/mdaq/mt/PollingMonitorAdapter.hpp \
	src/mdaq/mt/InterruptMonitorAdapter.hpp \
	src/mdaq/mt/Semaphore.h			\
	src/mdaq/mt/SysClock.h			\
	src/mdaq/mt/SigEvent.hpp
if !APPLE
mdaq_mt_include_HEADERS += \
	src/mdaq/mt/Timer.h \
	src/mdaq/mt/SigAction.h
endif

####################################################################################################
# 'mdaq-daq' library
####################################################################################################
ice_mdaq_daq_sources     = src/mdaq/daq/DataProcessor.ice src/mdaq/daq/RunProcessor.ice
ice_mdaq_daq_gen_sources = src/mdaq/daq/DataProcessor.cpp src/mdaq/daq/RunProcessor.cpp
ice_mdaq_daq_gen_headers = src/mdaq/daq/DataProcessor.h   src/mdaq/daq/RunProcessor.h

# Installation of Slice interface
mdaq_daq_slicedir = $(datadir)/slice/mdaq/daq
mdaq_daq_slice_DATA = $(ice_mdaq_daq_sources)

# Library
mdaq_daq_headers =                      \
	$(ice_mdaq_daq_gen_headers)         \
	src/mdaq/daq/CircularBuffer.h       \
	src/mdaq/daq/DataSender.h           \
	src/mdaq/daq/TcpDataSender.hpp      \
	src/mdaq/daq/TcpDataTransmit.h      \
	src/mdaq/daq/IceDataSender.h		\
	src/mdaq/daq/FdtDataSender.hpp

if !VXWORKS
mdaq_daq_headers += src/mdaq/daq/DataProcessorCore.h
mdaq_daq_headers += src/mdaq/daq/DataReceiver.h
mdaq_daq_headers += src/mdaq/daq/TcpDataReceiver.h
mdaq_daq_headers += src/mdaq/daq/IceDataReceiver.h
if !APPLE
mdaq_daq_headers += src/mdaq/daq/NetPerfProcessor.h
endif
endif

mdaq_daq_sources =				        \
	src/mdaq/daq/CircularBuffer.cpp     \
	src/mdaq/daq/DataSender.cpp         \
	src/mdaq/daq/IceDataSender.cpp
	
if !VXWORKS
mdaq_daq_sources += src/mdaq/daq/DataProcessorCore.cpp
mdaq_daq_sources += src/mdaq/daq/DataReceiver.cpp
mdaq_daq_sources += src/mdaq/daq/TcpDataReceiver.cpp
mdaq_daq_sources += src/mdaq/daq/IceDataReceiver.cpp
if !APPLE
mdaq_daq_sources += src/mdaq/daq/NetPerfProcessor.cpp
endif
endif

pkgconfig_DATA += mdaq-daq@DEBUG_SUFFIX@.pc

lib_LTLIBRARIES += libmdaq-daq@DEBUG_SUFFIX@.la
nodist_libmdaq_daq@DEBUG_SUFFIX@_la_SOURCES = $(ice_mdaq_daq_gen_sources)
libmdaq_daq@DEBUG_SUFFIX@_la_SOURCES = $(mdaq_daq_sources)
libmdaq_daq@DEBUG_SUFFIX@_la_LIBADD  = libmdaq-mt@DEBUG_SUFFIX@.la libmdaq-utl@DEBUG_SUFFIX@.la
libmdaq_daq@DEBUG_SUFFIX@_la_LIBADD += $(ICE_LIBS)
libmdaq_daq@DEBUG_SUFFIX@_la_LIBADD += $(COMPOUND_CONFIG_LIBS) $(UTILS_LIBS)

# Installation of header files
mdaq_daq_includedir  = $(includedir)/mdaq/daq
mdaq_daq_include_HEADERS = $(mdaq_daq_headers)

BUILT_SOURCES	+= $(ice_mdaq_daq_gen_headers) $(ice_mdaq_daq_gen_sources)
EXTRA_DIST      += $(ice_mdaq_daq_sources)

####################################################################################################
# Tests
####################################################################################################

noinst_PROGRAMS = mdaq_test mdaq_net_test mdaq_net_test_server

noinst_HEADERS =						\
	test/mdaq/utl/net/NetTestClient.h	\
	test/mdaq/utl/net/NetTestData.h

mdaq_test_SOURCES =	\
	test/main.cpp   \
	test/mdaq/hw/StdRegisterTest.cpp

if !VXWORKS
mdaq_test_SOURCES += test/mdaq/hw/RegisterProcessorTest.cpp
endif

mdaq_test_LDADD = 
mdaq_test_LDADD += $(ICE_LIBS)
mdaq_test_LDADD += $(UTILS_LIBS)
mdaq_test_LDADD += $(COMPOUND_CONFIG_LIBS)
mdaq_test_LDADD += libmdaq-utl@DEBUG_SUFFIX@.la
mdaq_test_LDADD += libmdaq-hwclient@DEBUG_SUFFIX@.la
if !VXWORKS
mdaq_test_LDADD += $(BOOST_LDFLAGS)
endif

noinst_PROGRAMS += threadsafe_logger_test
threadsafe_logger_test_SOURCES = test/mdaq/utl/ThreadSafeLoggerTest.cpp
threadsafe_logger_test_LDADD  = libmdaq-utl@DEBUG_SUFFIX@.la
threadsafe_logger_test_LDADD += libmdaq-mt@DEBUG_SUFFIX@.la
threadsafe_logger_test_LDADD += $(UTILS_LIBS)
if !VXWORKS
threadsafe_logger_test_LDADD += $(BOOST_LDFLAGS) $(BOOST_SYSTEM_LIB) $(BOOST_THREAD_LIB)
endif

mdaq_net_test_SOURCES =						\
	test/main.cpp							\
	test/mdaq/utl/net/NetTestClient.cpp		\
	test/mdaq/utl/net/NetTestData.cpp		\
	test/mdaq/utl/net/TcpConnectionTest.cpp	\
	test/mdaq/utl/net/UdpTransceiverTest.cpp
mdaq_net_test_LDADD  = libmdaq-utl@DEBUG_SUFFIX@.la
mdaq_net_test_LDADD += libmdaq-mt@DEBUG_SUFFIX@.la
mdaq_net_test_LDADD += $(ICE_LIBS)
mdaq_net_test_LDADD += $(UTILS_LIBS)
mdaq_net_test_LDADD += $(COMPOUND_CONFIG_LIBS)
if !VXWORKS
mdaq_net_test_LDADD += $(BOOST_LDFLAGS)
endif

mdaq_net_test_server_SOURCES = 				\
	test/mdaq/utl/net/NetTestData.cpp		\
	test/mdaq/utl/net/TestServer.cpp
mdaq_net_test_server_LDADD  = libmdaq-utl@DEBUG_SUFFIX@.la
mdaq_net_test_server_LDADD += libmdaq-mt@DEBUG_SUFFIX@.la
mdaq_net_test_server_LDADD += $(ICE_LIBS)
mdaq_net_test_server_LDADD += $(UTILS_LIBS) 
mdaq_net_test_server_LDADD += $(COMPOUND_CONFIG_LIBS)
mdaq_net_test_server_LDADD += $(BOOST_LDFLAGS)

if HAVE_ICESTORM

noinst_PROGRAMS += icestorm_logger_test_publisher icestorm_logger_test_subscriber

icestorm_logger_test_publisher_SOURCES = \
	test/mdaq/utl/IceStormLoggingBackendTest.cpp
icestorm_logger_test_publisher_LDADD  = libmdaq-utl@DEBUG_SUFFIX@.la libmdaq-mt@DEBUG_SUFFIX@.la
icestorm_logger_test_publisher_LDADD += $(ICE_LIBS)
icestorm_logger_test_publisher_LDADD += $(ICESTORM_LIBS)
icestorm_logger_test_publisher_LDADD += $(UTILS_LIBS)
icestorm_logger_test_publisher_LDADD += $(COMPOUND_CONFIG_LIBS)
icestorm_logger_test_publisher_LDADD += $(BOOST_LDFLAGS)

icestorm_logger_test_subscriber_SOURCES = \
	test/mdaq/utl/IceStormLoggingSubscriberTest.cpp
icestorm_logger_test_subscriber_LDADD  = libmdaq-utl@DEBUG_SUFFIX@.la
icestorm_logger_test_subscriber_LDADD += $(ICE_LIBS)
icestorm_logger_test_subscriber_LDADD += $(ICESTORM_LIBS)
icestorm_logger_test_subscriber_LDADD += $(UTILS_LIBS)
icestorm_logger_test_subscriber_LDADD += $(COMPOUND_CONFIG_LIBS)
icestorm_logger_test_subscriber_LDADD += $(BOOST_LDFLAGS)

endif

if HAVE_LOG4CXX
noinst_PROGRAMS += log4cxx_logger_test

log4cxx_logger_test_SOURCES = \
	test/mdaq/utl/Log4cxxLoggingBackendTest.cpp
log4cxx_logger_test_LDFLAGS = $(LOG4CXX_LDFLAGS)
log4cxx_logger_test_LDADD  = libmdaq-utl@DEBUG_SUFFIX@.la
log4cxx_logger_test_LDADD += libmdaq-mt@DEBUG_SUFFIX@.la
log4cxx_logger_test_LDADD += $(LOG4CXX_LIBS)
log4cxx_logger_test_LDADD += $(ICE_LIBS)
log4cxx_logger_test_LDADD += $(COMPOUND_CONFIG_LIBS) $(UTILS_LIBS)
endif

if !APPLE
noinst_PROGRAMS += perfsend
perfsend_SOURCES  = test/mdaq/daq/PerfSend.cpp
perfsend_LDADD  = libmdaq-daq@DEBUG_SUFFIX@.la
perfsend_LDADD += libmdaq-mt@DEBUG_SUFFIX@.la
perfsend_LDADD += libmdaq-utl@DEBUG_SUFFIX@.la
perfsend_LDADD += $(ICE_LIBS)
perfsend_LDADD += $(UTILS_LIBS)
perfsend_LDADD += $(COMPOUND_CONFIG_LIBS)
perfsend_LDADD += $(BOOST_LDFLAGS)

noinst_PROGRAMS += fileSend
fileSend_SOURCES  = test/mdaq/daq/FileSend.cpp
fileSend_LDADD  = libmdaq-daq@DEBUG_SUFFIX@.la
fileSend_LDADD += libmdaq-mt@DEBUG_SUFFIX@.la
fileSend_LDADD += libmdaq-utl@DEBUG_SUFFIX@.la
fileSend_LDADD += $(ICE_LIBS)
fileSend_LDADD += $(UTILS_LIBS)
fileSend_LDADD += $(COMPOUND_CONFIG_LIBS)
fileSend_LDADD += $(BOOST_LDFLAGS)
endif

if !VXWORKS
if !APPLE
noinst_PROGRAMS += perfrcve
perfrcve_SOURCES  = test/mdaq/daq/PerfRcve.cpp
perfrcve_LDFLAGS = $(ICE_LDFLAGS)
perfrcve_LDADD  = libmdaq-daq@DEBUG_SUFFIX@.la
perfrcve_LDADD += libmdaq-mt@DEBUG_SUFFIX@.la
perfrcve_LDADD += libmdaq-utl@DEBUG_SUFFIX@.la
perfrcve_LDADD += $(ICE_LIBS)
perfrcve_LDADD += $(UTILS_LIBS)
perfrcve_LDADD += $(COMPOUND_CONFIG_LIBS)
perfrcve_LDADD += $(BOOST_LDFLAGS)
endif
endif

if !VXWORKS
TESTS = mdaq_test test/net_test.sh
endif

####################################################################################################
# 'mdaq-hwclient' library
####################################################################################################
# Plain sources
mdaq_hw_sources =                           \
	src/mdaq/hw/RegisterProcessor.cpp       \
	src/mdaq/cfg/NodeAtTop.cpp

if !VXWORKS
mdaq_hw_sources += src/mdaq/hw/RemoteRegisterProcessor.cpp
endif

# Ice source and generated files
ice_mdaq_hw_sources      = src/mdaq/hw/Control.ice
ice_mdaq_hw_gen_sources  = src/mdaq/hw/Control.cpp
ice_mdaq_hw_gen_headers  = src/mdaq/hw/Control.h

# Installation of Slice interface
#mdaq_hw_slicedir = $(mdaq_slicedir)/hw
mdaq_hw_slice_DATA += $(ice_mdaq_hw_sources)

# Library
lib_LTLIBRARIES += libmdaq-hwclient@DEBUG_SUFFIX@.la
libmdaq_hwclient@DEBUG_SUFFIX@_la_SOURCES = $(mdaq_hw_sources)
 # These generated sources are needed but should not be distributed
nodist_libmdaq_hwclient@DEBUG_SUFFIX@_la_SOURCES  = $(ice_mdaq_hw_gen_sources)
 # These libraries are needed
libmdaq_hwclient@DEBUG_SUFFIX@_la_LDFLAGS = -version-info 0:0:0 -no-undefined
libmdaq_hwclient@DEBUG_SUFFIX@_la_LIBADD = libmdaq-interfaces@DEBUG_SUFFIX@.la
libmdaq_hwclient@DEBUG_SUFFIX@_la_LIBADD += libmdaq-utl@DEBUG_SUFFIX@.la
libmdaq_hwclient@DEBUG_SUFFIX@_la_LIBADD += $(COMPOUND_CONFIG_LIBS)
libmdaq_hwclient@DEBUG_SUFFIX@_la_LIBADD += $(TINY_XML_LIBS)
libmdaq_hwclient@DEBUG_SUFFIX@_la_LIBADD += $(ICE_LIBS)

pkgconfig_DATA += mdaq-hwclient@DEBUG_SUFFIX@.pc

# Installation of header files
mdaq_cfg_includedir = $(mdaq_includedir)/cfg
mdaq_cfg_include_HEADERS =						\
	src/mdaq/cfg/NodeAtTop.h

#mdaq_hw_includedir = $(mdaq_includedir)/hw
mdaq_hw_include_HEADERS +=                      \
	$(ice_mdaq_hw_gen_headers)                  \
	src/mdaq/hw/AbstractControllerDevice.hpp    \
	src/mdaq/hw/AbstractRegister.hpp            \
	src/mdaq/hw/ControllerDeviceSimulator.hpp   \
	src/mdaq/hw/DeviceStoragePolicy.hpp         \
	src/mdaq/hw/HwControllerException.h         \
	src/mdaq/hw/MemStoragePolicy.hpp            \
	src/mdaq/hw/ObservableValue.hpp             \
	src/mdaq/hw/RegisterProcessor.h             \
	src/mdaq/hw/RemoteRegisterProcessor.h       \
	src/mdaq/hw/SmartRegister.hpp				\
	src/mdaq/hw/StdRegister.hpp                 \
	src/mdaq/hw/NestedStoragePolicy.hpp

BUILT_SOURCES	+= $(ice_mdaq_hw_gen_headers)  $(ice_mdaq_hw_gen_sources)
EXTRA_DIST      += $(ice_mdaq_hw_sources)

####################################################################################################
# 'mdaq-hwserver' library
####################################################################################################
# Source files
hwserver_sources =								\
	src/mdaq/hw/server/BaseDynDevice.cpp		\
	src/mdaq/hw/server/BaseDynRegister.cpp		\
	src/mdaq/hw/server/DeviceCreationParams.cpp \
	src/mdaq/hw/server/DeviceFactory.cpp		\
	src/mdaq/hw/server/DeviceImpl.cpp			\
	src/mdaq/hw/server/HwNodeImpl.cpp

# Library
lib_LTLIBRARIES += libmdaq-hwserver@DEBUG_SUFFIX@.la
libmdaq_hwserver@DEBUG_SUFFIX@_la_SOURCES = $(hwserver_sources) $(ice_mdaq_hw_gen_sources)
libmdaq_hwserver@DEBUG_SUFFIX@_la_LDFLAGS = -version-info 0:0:0 -no-undefined
libmdaq_hwserver@DEBUG_SUFFIX@_la_LIBADD = libmdaq-utl@DEBUG_SUFFIX@.la
libmdaq_hwserver@DEBUG_SUFFIX@_la_LIBADD += libmdaq-hwclient@DEBUG_SUFFIX@.la
libmdaq_hwserver@DEBUG_SUFFIX@_la_LIBADD += libmdaq-interfaces@DEBUG_SUFFIX@.la
libmdaq_hwserver@DEBUG_SUFFIX@_la_LIBADD += $(ICE_LIBS) $(LOG4CXX_LIBS)
libmdaq_hwserver@DEBUG_SUFFIX@_la_LIBADD += $(COMPOUND_CONFIG_LIBS) $(TINY_XML_LIBS)

# Installation of header files
mdaq_hw_server_includedir = $(mdaq_hw_includedir)/server
mdaq_hw_server_include_HEADERS =               \
	src/mdaq/hw/server/BaseDynDevice.h         \
	src/mdaq/hw/server/BaseDynRegister.h       \
	src/mdaq/hw/server/DeviceCreationParams.h  \
	src/mdaq/hw/server/DeviceFactory.h         \
	src/mdaq/hw/server/DeviceImpl.h            \
	src/mdaq/hw/server/DynDevice.hpp           \
	src/mdaq/hw/server/DynRegister.hpp         \
	src/mdaq/hw/server/HwNodeImpl.h

pkgconfig_DATA += mdaq-hwserver@DEBUG_SUFFIX@.pc

####################################################################################################
# 'mdaq-eccserver' library
####################################################################################################
if !VXWORKS
lib_LTLIBRARIES += libmdaq-eccserver@DEBUG_SUFFIX@.la
libmdaq_eccserver@DEBUG_SUFFIX@_la_SOURCES = src/mdaq/EccBackEnd.cpp src/mdaq/rc/HardwareNode.cpp src/mdaq/EccImpl.cpp

 # These libraries are needed
libmdaq_eccserver@DEBUG_SUFFIX@_la_LDFLAGS = -version-info 0:0:0 -no-undefined
libmdaq_eccserver@DEBUG_SUFFIX@_la_LIBADD  = libmdaq-hwclient@DEBUG_SUFFIX@.la
libmdaq_eccserver@DEBUG_SUFFIX@_la_LIBADD += libmdaq-interfaces@DEBUG_SUFFIX@.la
#libmdaq_eccserver@DEBUG_SUFFIX@_la_LIBADD += 
libmdaq_eccserver@DEBUG_SUFFIX@_la_LIBADD += $(ICE_LIBS) $(COMPOUND_CONFIG_LIBS) $(TINY_XML_LIBS) $(UTILS_LIBS)

pkgconfig_DATA += mdaq-eccserver@DEBUG_SUFFIX@.pc
endif #!VXWORKS

mdaq_rc_includedir = $(mdaq_includedir)/rc
mdaq_rc_include_HEADERS =	     \
	src/mdaq/rc/HardwareNode.h

####################################################################################################
# 'HwServerTest_local'
####################################################################################################
noinst_PROGRAMS += HwServerTest_local
HwServerTest_local_SOURCES =                    \
	test/main.cpp								\
	test/mdaq/hw/server/DeviceFactoryTest.cpp	\
	test/mdaq/hw/server/DynDeviceTest.cpp		\
	test/mdaq/hw/server/DynRegisterTest.cpp

HwServerTest_local_LDADD  = libmdaq-hwserver@DEBUG_SUFFIX@.la 
HwServerTest_local_LDADD += libmdaq-hwclient@DEBUG_SUFFIX@.la 
HwServerTest_local_LDADD += libmdaq-utl@DEBUG_SUFFIX@.la 
HwServerTest_local_LDADD += libmdaq-daq@DEBUG_SUFFIX@.la
HwServerTest_local_LDADD += libmdaq-interfaces@DEBUG_SUFFIX@.la 
HwServerTest_local_LDADD += libmdaq-mt@DEBUG_SUFFIX@.la
HwServerTest_local_LDADD += $(ICE_LIBS)
HwServerTest_local_LDADD += $(UTILS_LIBS)
HwServerTest_local_LDADD += $(COMPOUND_CONFIG_LIBS)
HwServerTest_local_LDADD += $(BOOST_LDFLAGS)

####################################################################################################
# 'HwServerTest'
####################################################################################################
noinst_PROGRAMS += HwServerTest_server
HwServerTest_server_SOURCES = 					\
	test/mdaq/hw/server/HwServerMain.cpp		\
	test/mdaq/hw/server/ServerDevices.cpp
	
# These libraries are needed
HwServerTest_server_LDADD  = libmdaq-hwserver@DEBUG_SUFFIX@.la
HwServerTest_server_LDADD += libmdaq-hwclient@DEBUG_SUFFIX@.la
HwServerTest_server_LDADD += libmdaq-utl@DEBUG_SUFFIX@.la
HwServerTest_server_LDADD += libmdaq-daq@DEBUG_SUFFIX@.la
HwServerTest_server_LDADD += libmdaq-interfaces@DEBUG_SUFFIX@.la
HwServerTest_server_LDADD += libmdaq-mt@DEBUG_SUFFIX@.la
HwServerTest_server_LDADD += $(ICE_LIBS) 
HwServerTest_server_LDADD += $(COMPOUND_CONFIG_LIBS)
HwServerTest_server_LDADD += $(UTILS_LIBS)
HwServerTest_server_LDADD += $(BOOST_LDFLAGS)

####################################################################################################
# 'HwServerTest_client'
####################################################################################################
if !VXWORKS
noinst_PROGRAMS += HwServerTest_client
HwServerTest_client_SOURCES = 					\
	test/main.cpp								\
	test/mdaq/hw/server/HwNodeTest.cpp


HwServerTest_client_LDADD  = libmdaq-hwclient@DEBUG_SUFFIX@.la
HwServerTest_client_LDADD += libmdaq-utl@DEBUG_SUFFIX@.la
HwServerTest_client_LDADD += libmdaq-mt@DEBUG_SUFFIX@.la
HwServerTest_client_LDADD += libmdaq-interfaces@DEBUG_SUFFIX@.la
HwServerTest_client_LDFLAGS = $(COMPOUND_CONFIG_LIBS) $(TINY_XML_LIBS) $(UTILS_LIBS)
HwServerTest_client_LDADD += $(ICE_LIBS)
endif # !VXWORKS

####################################################################################################
# 'hwServer'
####################################################################################################
bin_PROGRAMS = hwServer 
hwServer_SOURCES  = src/mdaq/HwServer.cpp		
if VXWORKS
hwServer_SOURCES +=	src/ServerDevices.cpp
else
hwServer_SOURCES +=	src/SimServerDevices.cpp
endif

# These libraries are needed 
hwServer_LDADD  = libmdaq-hwclient@DEBUG_SUFFIX@.la 
hwServer_LDADD += libmdaq-hwserver@DEBUG_SUFFIX@.la
hwServer_LDADD += libmdaq-utl@DEBUG_SUFFIX@.la
hwServer_LDADD += libmdaq-interfaces@DEBUG_SUFFIX@.la
hwServer_LDADD += $(ICE_LIBS)
hwServer_LDADD += $(UTILS_LIBS)
hwServer_LDADD += $(COMPOUND_CONFIG_LIBS)
hwServer_LDADD += $(BOOST_LDFLAGS)

####################################################################################################
# 'hwLoad'
####################################################################################################
if !VXWORKS
bin_PROGRAMS += hwLoad
hwLoad_SOURCES =		\
	demo/mdaq/hw/hwLoad.cpp      

# These libraries are needed 
hwLoad_LDADD  = libmdaq-interfaces@DEBUG_SUFFIX@.la
hwLoad_LDADD += libmdaq-hwclient@DEBUG_SUFFIX@.la 
hwLoad_LDADD += libmdaq-utl@DEBUG_SUFFIX@.la
hwLoad_LDADD += $(ICE_LIBS)
hwLoad_LDADD += $(UTILS_LIBS) 
hwLoad_LDADD += $(COMPOUND_CONFIG_LIBS)
endif # !VXWORKS

####################################################################################################
# 'hwFieldWrite'
####################################################################################################
if !VXWORKS
bin_PROGRAMS += hwFieldWrite
hwFieldWrite_SOURCES =		\
	demo/mdaq/hw/hwFieldWrite.cpp      	

# These libraries are needed
hwFieldWrite_LDADD  = libmdaq-interfaces@DEBUG_SUFFIX@.la
hwFieldWrite_LDADD += libmdaq-hwclient@DEBUG_SUFFIX@.la 
hwFieldWrite_LDADD += libmdaq-utl@DEBUG_SUFFIX@.la
hwFieldWrite_LDADD +=  $(ICE_LIBS)
hwFieldWrite_LDADD += $(COMPOUND_CONFIG_LIBS) $(UTILS_LIBS)
endif # !VXWORKS

####################################################################################################
# 'hwFieldRead'
####################################################################################################
if !VXWORKS
bin_PROGRAMS += hwFieldRead
hwFieldRead_SOURCES =		\
	demo/mdaq/hw/hwFieldRead.cpp      	

# These libraries are needed
hwFieldRead_LDADD  = libmdaq-utl@DEBUG_SUFFIX@.la
hwFieldRead_LDADD += libmdaq-interfaces@DEBUG_SUFFIX@.la
hwFieldRead_LDADD += libmdaq-hwclient@DEBUG_SUFFIX@.la
hwFieldRead_LDADD += $(COMPOUND_CONFIG_LIBS) $(ICE_LIBS)
endif # !VXWORKS

####################################################################################################
## Programs for the use of mdaq/ecc
####################################################################################################

if !VXWORKS

bin_PROGRAMS += eccServer@DEBUG_SUFFIX@

eccServer@DEBUG_SUFFIX@_SOURCES =		\
	demo/eccServer.cpp
eccServer@DEBUG_SUFFIX@_LDADD  = libmdaq-interfaces@DEBUG_SUFFIX@.la
eccServer@DEBUG_SUFFIX@_LDADD += libmdaq-eccserver@DEBUG_SUFFIX@.la
eccServer@DEBUG_SUFFIX@_LDADD += libmdaq-hwclient@DEBUG_SUFFIX@.la
eccServer@DEBUG_SUFFIX@_LDADD += libmdaq-utl@DEBUG_SUFFIX@.la
eccServer@DEBUG_SUFFIX@_LDADD += $(ICE_LIBS) $(ICESTORM_LIBS)
eccServer@DEBUG_SUFFIX@_LDADD += $(UTILS_LIBS)
eccServer@DEBUG_SUFFIX@_LDADD += $(COMPOUND_CONFIG_LIBS)

endif #!VXWORKS

if ENABLE_CLI
if !ICEE

bin_PROGRAMS += eccClient@DEBUG_SUFFIX@

eccClient@DEBUG_SUFFIX@_SOURCES =		\
	demo/eccClient.cpp
eccClient@DEBUG_SUFFIX@_LDADD  = libmdaq-cli@DEBUG_SUFFIX@.la
eccClient@DEBUG_SUFFIX@_LDADD += libmdaq-interfaces@DEBUG_SUFFIX@.la
eccClient@DEBUG_SUFFIX@_LDADD += libmdaq-utl@DEBUG_SUFFIX@.la
eccClient@DEBUG_SUFFIX@_LDADD += $(ICE_LIBS)
eccClient@DEBUG_SUFFIX@_LDADD += $(ICESTORM_LIBS)
eccClient@DEBUG_SUFFIX@_LDADD += $(UTILS_LIBS)
eccClient@DEBUG_SUFFIX@_LDADD += $(COMPOUND_CONFIG_LIBS)

endif #!ICEE
endif #ENABLE_CLI

####################################################################################################
## Demos for the use of mdaq (mt)
####################################################################################################

if VXWORKS

bin_PROGRAMS += demoMonitor
demoMonitor_SOURCES =				\
	demo/mdaq/mt/demoMonitor.cpp  	
demoMonitor_LDADD  = libmdaq-hwclient@DEBUG_SUFFIX@.la
demoMonitor_LDADD += $(UTILS_LIBS)
demoMonitor_LDADD += $(BOOST_LDFLAGS)

endif # VXWORKS

####################################################################################################

include Makefile.VxWorks.am
include Makefile.packages.am

CLEANFILES		+= $(BUILT_SOURCES)
DISTCLEANFILES   = $(CLEANFILES)
