include $(top_srcdir)/am/doxygen-ext.am
include $(top_srcdir)/am/ice.am

ACLOCAL_AMFLAGS     = -I m4
AUTOMAKE_OPTIONS    = subdir-objects

dist_noinst_SCRIPTS = autogen.sh cleanup.sh
pkgconfigdir        = $(libdir)/pkgconfig

pkgconfig_DATA  =
lib_LTLIBRARIES = 
bin_PROGRAMS    =

AM_CPPFLAGS         =               \
	-I@top_srcdir@/src              \
	-I@top_srcdir@/src/get          \
	-I@top_builddir@/src            \
	-I@top_builddir@/src/get        \
	-I@top_builddir@/src/get/rc     \
	-I@top_builddir@/src/get/mt     \
	-I@top_builddir@/src/get/hw     \
	-I@top_builddir@/src/get/daq    \
	-I@top_builddir@/src/get/soap   \
	-I@top_builddir@/src/get/cobo   \
	-I@top_builddir@/src/get/mutant \
	$(ICE_CPPFLAGS)                 \
	$(MDAQ_UTL_CFLAGS)              \
	$(MDAQ_HWSERVER_CFLAGS)         \
	$(MDAQ_CLI_CFLAGS)              \
	$(MDAQ_DAQ_CFLAGS)              \
	$(UTILS_CFLAGS)                 \
	$(MFM_CFLAGS)                   \
	$(GSOAP_CFLAGS)                 \
	$(DHSM_CFLAGS)                  \
	$(ICESTORM_CPPFLAGS)            \
	$(BOOST_CPPFLAGS)               \
	$(LOG4CXX_CPPFLAGS)             \
	$(XGPIO_CFLAGS)                 \
	$(UTL_LOG_THREADSAFE_CPPFLAGS)

AM_CXXFLAGS = $(PTHREAD_CFLAGS)
AM_LDFLAGS = $(ICE_LDFLAGS) $(READLINE_LDFLAGS) $(ICE_LDFLAGS)
LIBS += $(UTL_LOG_THREADSAFE_LIBS)

################################################################################
# slice 
################################################################################

#####################################
# slice2cpp options
#####################################

SPPFLAGS = -I. -I$(MDAQ_SLICEDIR) -I@datadir@/slice -I$(top_srcdir)/src

#####################################
# slice sources
#####################################

ice_ecc_sources            = src/get/GetEcc.ice
ice_ecc_gen_sources        = src/get/GetEcc.cpp
noinst_ice_ecc_gen_headers = src/get/GetEcc.h

slicedir = $(datadir)/slice
get_slicedir = $(slicedir)/get
get_slice_DATA = $(ice_ecc_sources)

get_daq_slicedir = $(get_slicedir)/daq
get_daq_slice_DATA = $(ice_dataSender_sources)

################################################################################
# installation of include files
################################################################################

get_daq_includedir = $(includedir)/get/daq
get_daq_include_HEADERS =                   \
	src/get/daq/SyncedCircularBuffer.h      \
	src/get/daq/MemRead.h            \
	$(noinst_ice_dataSender_gen_headers)       \
	src/get/daq/DaqCtrlNodeI.h  \
	src/get/daq/FakeDataSender.h \
	src/get/daq/FdtDataReceiver.h
if HAVE_MFM
get_daq_include_HEADERS += src/get/daq/FakeFrameSender.h
endif

if BUILD_ECC
get_daq_include_HEADERS +=           \
	src/get/daq/DataRouter.h         \
	src/get/daq/FrameStorage.h		 \
	src/get/daq/ByteStorage.h
if !APPLE
get_daq_include_HEADERS += 	src/get/daq/FrameCounter.h
endif #!APPLE
endif #BUILD_ECC

get_asad_includedir = $(includedir)/get/asad
get_asad_include_HEADERS = src/get/asad/utl.h 

get_includedir = $(includedir)/get
get_include_HEADERS =                \
	src/get/CoBoCircularBufferRegisters.h    \
	$(noinst_ice_ecc_gen_headers)    \
	src/get/GetEccImpl.h \
	src/get/DefaultPortNums.h

get_rc_includedir = $(includedir)/get/rc
get_rc_include_HEADERS =                \
	src/get/rc/SubConfigIdSet.h		\
	src/get/rc/ConfigId.h	\
	src/get/rc/ConfigIdSet.h	\
	src/get/rc/DataSenderId.h \
	src/get/rc/DataSenderSet.h \
	src/get/rc/DataRouterId.h \
	src/get/rc/DataLink.h \
	src/get/rc/DataLinkSet.h \
	src/get/rc/CoBoPLLManager.h \
	src/get/rc/BackEnd.h	\
	src/get/rc/ConfigManager.h	\
	src/get/rc/ConfigFileManager.h	\
	src/get/rc/SM.h	\
	src/get/rc/NodeId.h \
	src/get/rc/Node.h \
	src/get/rc/CoBoNode.h \
	src/get/rc/CoBoConfigParser.h \
	src/get/rc/MutantConfigParser.h \
	src/get/rc/MutantNode.h \
	src/get/rc/SystemManager.h \
	src/get/rc/AmcSlotManager.h
if CFG_DB_SUPPORT_ENABLED
get_rc_include_HEADERS += src/get/rc/ConfigDatabaseManager.h
endif

get_drv_includedir = $(get_includedir)/drv
get_drv_include_HEADERS =                \
	src/get/drv/AGetControllerDevice.hpp \
	src/get/drv/AGetControllerDevice.h   \
	src/get/drv/AsAdControllerDevice.hpp \
	src/get/drv/AsAdControllerDevice.h

get_zboard_hw_includedir = $(get_includedir)/zboard/hw
get_zboard_hw_include_HEADERS = \
	src/get/zboard/hw/ControllerDeviceAdc.hpp \
	src/get/zboard/hw/ControllerDeviceAget.hpp \
	src/get/zboard/hw/ControllerDeviceDac.hpp \
	src/get/zboard/hw/ControllerDeviceIsp.hpp
	
get_utl_includedir = $(get_includedir)/utl
get_utl_include_HEADERS  =

if BUILD_ECC
get_utl_include_HEADERS += src/get/utl/Time.h 
endif #BUILD_ECC

#####################################
# CoBo-Mutant link convenience library
#####################################

ice_mutant_link_sources = src/get/cobo/MutantLinkManager.ice
ice_mutant_link_gen_sources = src/get/cobo/MutantLinkManager.cpp
ice_mutant_link_gen_headers = src/get/cobo/MutantLinkManager.h

cobo_mutant_link_includedir = $(includedir)/get/cobo
cobo_mutant_link_include_HEADERS = $(ice_mutant_link_gen_headers) 

noinst_LTLIBRARIES = libcobo-mutant-link@DEBUG_SUFFIX@.la
nodist_libcobo_mutant_link@DEBUG_SUFFIX@_la_SOURCES = $(ice_mutant_link_gen_sources) $(ice_mutant_link_gen_headers)
libcobo_mutant_link@DEBUG_SUFFIX@_la_LIBADD = $(ICE_LIBS)

#####################################
# CoBo CtrlNode interface convenience library
#####################################

ice_cobo_sources = src/get/cobo/CtrlNode.ice
ice_cobo_gen_sources = src/get/cobo/CtrlNode.cpp
ice_cobo_gen_headers = src/get/cobo/CtrlNode.h

noinst_LTLIBRARIES += libcobo-ctrl-node@DEBUG_SUFFIX@.la
nodist_libcobo_ctrl_node@DEBUG_SUFFIX@_la_SOURCES = src/get/cobo/CtrlNode.cpp src/get/cobo/CtrlNode.h
libcobo_ctrl_node@DEBUG_SUFFIX@_la_LIBADD = $(ICE_LIBS)

#####################################
# Alarm interface convenience library
#####################################

ice_mt_alarm_sources = src/get/mt/AlarmService.ice
ice_mt_alarm_gen_sources = src/get/mt/AlarmService.cpp
ice_mt_alarm_gen_headers = src/get/mt/AlarmService.h

get_mt_alarm_includedir = $(includedir)/get/mt
get_mt_alarm_include_HEADERS = $(ice_mt_alarm_gen_headers) src/get/mt/AlarmLogger.h

noinst_LTLIBRARIES += libmt-alarm@DEBUG_SUFFIX@.la
nodist_libmt_alarm@DEBUG_SUFFIX@_la_SOURCES = $(ice_mt_alarm_gen_sources) $(ice_mt_alarm_gen_headers)
libmt_alarm@DEBUG_SUFFIX@_la_LIBADD = $(ICE_LIBS)

#####################################
# CoBo LEDs convenience library
#####################################

ice_cobo_sources += src/get/cobo/LedManager.ice
ice_cobo_gen_sources += src/get/cobo/LedManager.cpp
ice_cobo_gen_headers += src/get/cobo/LedManager.h

get_cobo_includedir = $(includedir)/get/cobo
get_cobo_include_HEADERS = $(ice_cobo_gen_headers) 

noinst_LTLIBRARIES += libcobo-led-mgr@DEBUG_SUFFIX@.la
nodist_libcobo_led_mgr@DEBUG_SUFFIX@_la_SOURCES = src/get/cobo/LedManager.cpp src/get/cobo/LedManager.h
libcobo_led_mgr@DEBUG_SUFFIX@_la_LIBADD = $(ICE_LIBS)

#####################################
# Monitoring convenience library
#####################################
if BUILD_ECC

noinst_LTLIBRARIES += libgetbench-mt@DEBUG_SUFFIX@.la
libgetbench_mt@DEBUG_SUFFIX@_la_SOURCES = src/get/mt/AlarmLogger.cpp
libgetbench_mt@DEBUG_SUFFIX@_la_LIBADD  = libmt-alarm@DEBUG_SUFFIX@.la
libgetbench_mt@DEBUG_SUFFIX@_la_LIBADD  += $(UTILS_LIBS)
libgetbench_mt@DEBUG_SUFFIX@_la_LIBADD  += $(ICE_LIBS)

endif #BUILD_ECC

######################################
# Mutant interface convenience library
######################################

if BUILD_ECC

ice_mutant_sources     = src/get/mutant/AlignmentServer.ice
ice_mutant_gen_sources = src/get/mutant/AlignmentServer.cpp
ice_mutant_gen_headers = src/get/mutant/AlignmentServer.h

get_mutant_includedir = $(includedir)/get/mutant
get_mutant_include_HEADERS = $(ice_mutant_gen_headers)

noinst_LTLIBRARIES += libgetbench-mutant@DEBUG_SUFFIX@.la
nodist_libgetbench_mutant@DEBUG_SUFFIX@_la_SOURCES = $(ice_mutant_gen_sources) $(ice_mutant_gen_headers)
libgetbench_mutant@DEBUG_SUFFIX@_la_LIBADD = $(ICE_LIBS)

endif #BUILD_ECC

#__________________________________________________________________________________________________
# Library libgetbench-eccclient

eccclient_sources =         \
	$(ice_ecc_gen_sources)	\
	src/get/rc/SM.cpp       \
	src/get/asad/utl.cpp

if BUILD_ECC

pkgconfig_DATA += getbench-eccclient@DEBUG_SUFFIX@.pc

lib_LTLIBRARIES += libgetbench-eccclient@DEBUG_SUFFIX@.la
libgetbench_eccclient@DEBUG_SUFFIX@_la_SOURCES = $(eccclient_sources) 
libgetbench_eccclient@DEBUG_SUFFIX@_la_LDFLAGS = $(ICE_LDFLAGS) -version-info 1:0:0 -no-undefined
libgetbench_eccclient@DEBUG_SUFFIX@_la_LIBADD  = libcobo-mutant-link@DEBUG_SUFFIX@.la
libgetbench_eccclient@DEBUG_SUFFIX@_la_LIBADD  += libgetbench-mutant@DEBUG_SUFFIX@.la
libgetbench_eccclient@DEBUG_SUFFIX@_la_LIBADD += libgetbench-mt@DEBUG_SUFFIX@.la
libgetbench_eccclient@DEBUG_SUFFIX@_la_LIBADD += libcobo-led-mgr@DEBUG_SUFFIX@.la
libgetbench_eccclient@DEBUG_SUFFIX@_la_LIBADD += $(MDAQ_UTL_LIBS) $(MDAQ_DAQ_LIBS) $(MDAQ_IF_LIBS) $(MDAQ_CLI_LIBS)
libgetbench_eccclient@DEBUG_SUFFIX@_la_LIBADD += $(ICE_LIBS)

#__________________________________________________________________________________________________
# getEccClient

bin_PROGRAMS += getEccClient@DEBUG_SUFFIX@

getEccClient@DEBUG_SUFFIX@_SOURCES = \
	src/main/getEccClient.cpp \
	src/get/GetEccCmd.cpp \
	src/get/GetEccCmd.h

getEccClient@DEBUG_SUFFIX@_LDADD  = libgetbench-eccclient@DEBUG_SUFFIX@.la
getEccClient@DEBUG_SUFFIX@_LDADD += $(MDAQ_UTL_LIBS) $(MDAQ_CLI_LIBS) $(MDAQ_HWCLIENT_LIBS) $(MDAQ_IF_LIBS) $(MFM_LIBS)
getEccClient@DEBUG_SUFFIX@_LDADD += $(ICE_LIBS) $(ICESTORM_LIBS)

#__________________________________________________________________________________________________
# Library libgetbench-cfg
pkgconfig_DATA += getbench-cfg@DEBUG_SUFFIX@.pc
lib_LTLIBRARIES += libgetbench-cfg@DEBUG_SUFFIX@.la
libgetbench_cfg@DEBUG_SUFFIX@_la_SOURCES = \
	src/get/rc/CoBoConfigParser.cpp	src/get/rc/MutantConfigParser.cpp
libgetbench_cfg@DEBUG_SUFFIX@_la_LDFLAGS = -version-info 1:0:0 -no-undefined
libgetbench_cfg@DEBUG_SUFFIX@_la_LIBADD  = $(UTILS_LIBS) $(COMPOUND_CONFIG_LIBS) $(MDAQ_IF_LIBS) $(MDAQ_UTL_LIBS)

endif #BUILD_ECC

#__________________________________________________________________________________________________
# Library libgetbench-eccserver
eccserver_sources =                         \
	src/get/GetEccImpl.cpp                  \
	src/get/rc/SubConfigIdSet.cpp           \
	src/get/rc/ConfigId.cpp	                \
	src/get/rc/ConfigIdSet.cpp              \
	src/get/rc/DataSenderId.cpp             \
	src/get/rc/DataSenderSet.cpp            \
	src/get/rc/DataRouterId.cpp             \
	src/get/rc/DataLink.cpp                 \
	src/get/rc/DataLinkSet.cpp              \
	src/get/rc/CoBoPLLManager.cpp           \
	src/get/rc/BackEnd.cpp                  \
	src/get/rc/ConfigManager.cpp            \
	src/get/rc/ConfigFileManager.cpp        \
	src/get/rc/NodeId.cpp                   \
	src/get/rc/Node.cpp                     \
	src/get/rc/CoBoNode.cpp                 \
	src/get/rc/MutantNode.cpp               \
	src/get/rc/SystemManager.cpp            \
	src/get/rc/AmcSlotManager.cpp
if CFG_DB_SUPPORT_ENABLED
eccserver_sources += src/get/rc/ConfigDatabaseManager.cpp
endif

if BUILD_ECC

pkgconfig_DATA += getbench-eccserver@DEBUG_SUFFIX@.pc

lib_LTLIBRARIES += libgetbench-eccserver@DEBUG_SUFFIX@.la
libgetbench_eccserver@DEBUG_SUFFIX@_la_SOURCES = $(eccserver_sources)
libgetbench_eccserver@DEBUG_SUFFIX@_CPPFLAGS   = 
libgetbench_eccserver@DEBUG_SUFFIX@_la_LDFLAGS = $(MDAQ_ECCSERVER_LDFLAGS) $(ICE_LDFLAGS) -version-info 1:0:0 -no-undefined
libgetbench_eccserver@DEBUG_SUFFIX@_la_LIBADD  = libcobo-ctrl-node@DEBUG_SUFFIX@.la
libgetbench_eccserver@DEBUG_SUFFIX@_la_LIBADD += libgetbench-cfg@DEBUG_SUFFIX@.la
libgetbench_eccserver@DEBUG_SUFFIX@_la_LIBADD += libgetbench-eccclient@DEBUG_SUFFIX@.la
libgetbench_eccserver@DEBUG_SUFFIX@_la_LIBADD += $(MDAQ_ECCSERVER_LIBS)
libgetbench_eccserver@DEBUG_SUFFIX@_la_LIBADD += $(COMPOUND_CONFIG_LIBS)
if CFG_DB_SUPPORT_ENABLED
libgetbench_eccserver@DEBUG_SUFFIX@_la_LIBADD += $(CONFIG_DATABASE_LIBS)
endif
libgetbench_eccserver@DEBUG_SUFFIX@_la_LIBADD += $(UTILS_LIBS)
libgetbench_eccserver@DEBUG_SUFFIX@_la_LIBADD += $(ICE_LIBS)
libgetbench_eccserver@DEBUG_SUFFIX@_la_LIBADD += $(BOOST_FILESYSTEM_LIB) $(BOOST_SYSTEM_LIB) $(BOOST_THREAD_LIB)
libgetbench_eccserver@DEBUG_SUFFIX@_la_LIBADD += $(TINY_XML_LIBS)
libgetbench_eccserver@DEBUG_SUFFIX@_la_LIBADD += $(DHSM_LIBS)

# $(BOOST_LDFLAGS) -lboost_timer
endif #BUILD_ECC

#####################################
# getEccServer
#####################################

if BUILD_ECC

bin_PROGRAMS += getEccServer@DEBUG_SUFFIX@

getEccServer@DEBUG_SUFFIX@_SOURCES =             \
	$(ice_dataSender_gen_sources)     \
	src/main/getEccServer.cpp
getEccServer@DEBUG_SUFFIX@_LDADD  = libgetbench-cfg@DEBUG_SUFFIX@.la
getEccServer@DEBUG_SUFFIX@_LDADD += libgetbench-eccserver@DEBUG_SUFFIX@.la
getEccServer@DEBUG_SUFFIX@_LDADD += libgetbench-eccclient@DEBUG_SUFFIX@.la
getEccServer@DEBUG_SUFFIX@_LDADD += libcobo-ctrl-node@DEBUG_SUFFIX@.la
getEccServer@DEBUG_SUFFIX@_LDADD += $(COMPOUND_CONFIG_LIBS)
getEccServer@DEBUG_SUFFIX@_LDADD += $(MDAQ_ECCSERVER_LIBS)
getEccServer@DEBUG_SUFFIX@_LDADD += $(MDAQ_UTL_LIBS)
getEccServer@DEBUG_SUFFIX@_LDADD += $(BOOST_FILESYSTEM_LIB) $(BOOST_SYSTEM_LIB) $(BOOST_PROGRAM_OPTIONS_LIB)
getEccServer@DEBUG_SUFFIX@_LDADD += $(ICE_LIBS) $(ICESTORM_LIBS)
if HAVE_LOG4CXX
getEccServer@DEBUG_SUFFIX@_LDADD += $(LOG4CXX_LIBS)
endif
endif #BUILD_ECC

#####################################
# getEccSoapServer
#####################################
if HAVE_SOAP

soap_eccServer_gen_sources = \
	src/get/soap/eccC.cpp \
	src/get/soap/eccServer.cpp

soap_eccServer_gen_extra = \
	src/get/soap/eccClient.cpp \
	src/get/soap/ecceccObject.h \
	src/get/soap/ecceccProxy.h \
	src/get/soap/eccH.h \
	src/get/soap/ecc.nsmap \
	src/get/soap/eccStub.h \
	src/get/soap/ecc.wsdl \
	src/get/soap/ecc.xsd

soap_sources = $(srcdir)/src/get/soap/ecc.hh
	
$(soap_eccServer_gen_sources) $(soap_eccServer_gen_extra) : $(soap_sources)
	$(SOAP_CPP2) $(SOAP_CPP2_FLAGS) -L -x -pecc -n -d@top_builddir@/src/get/soap $(soap_sources)

bin_PROGRAMS += getEccSoapServer@DEBUG_SUFFIX@

getEccSoapServer@DEBUG_SUFFIX@_SOURCES = src/get/soap/getEccSoapServer.cpp src/get/soap/EccSoapServer.cpp src/get/soap/EccSoapServer.h
nodist_getEccSoapServer@DEBUG_SUFFIX@_SOURCES = $(soap_eccServer_gen_sources)  
getEccSoapServer@DEBUG_SUFFIX@_LDADD  = libgetbench-cfg@DEBUG_SUFFIX@.la
getEccSoapServer@DEBUG_SUFFIX@_LDADD += libgetbench-eccclient@DEBUG_SUFFIX@.la
getEccSoapServer@DEBUG_SUFFIX@_LDADD += libgetbench-eccserver@DEBUG_SUFFIX@.la
getEccSoapServer@DEBUG_SUFFIX@_LDADD += libcobo-ctrl-node@DEBUG_SUFFIX@.la
getEccSoapServer@DEBUG_SUFFIX@_LDADD += $(MDAQ_ECCSERVER_LIBS)
getEccSoapServer@DEBUG_SUFFIX@_LDADD += $(MDAQ_UTL_LIBS)
getEccSoapServer@DEBUG_SUFFIX@_LDADD += $(COMPOUND_CONFIG_LIBS)
if CFG_DB_SUPPORT_ENABLED
getEccSoapServer@DEBUG_SUFFIX@_LDADD += $(CONFIG_DATABASE_LIBS)
endif
if HAVE_LOG4CXX
getEccSoapServer@DEBUG_SUFFIX@_LDADD += $(LOG4CXX_LIBS)
endif
getEccSoapServer@DEBUG_SUFFIX@_LDADD += $(BOOST_FILESYSTEM_LIB)  $(BOOST_SYSTEM_LIB) $(BOOST_PROGRAM_OPTIONS_LIB)
getEccSoapServer@DEBUG_SUFFIX@_LDADD += $(ICE_LIBS)
getEccSoapServer@DEBUG_SUFFIX@_LDADD += $(GSOAP_LIBS) -lpthread
getEccSoapServer@DEBUG_SUFFIX@_CPPFLAGS = $(AM_CPPFLAGS) $(VLA_CPPFLAGS)

#####################################
# getEccSoapClient
#####################################

bin_PROGRAMS += getEccSoapClient@DEBUG_SUFFIX@

getEccSoapClient@DEBUG_SUFFIX@_SOURCES = src/get/soap/getEccSoapClient.cpp
nodist_getEccSoapClient@DEBUG_SUFFIX@_SOURCES = \
	src/get/soap/eccC.cpp \
	src/get/soap/eccClient.cpp
	
getEccSoapClient@DEBUG_SUFFIX@_LDADD  = libgetbench-eccclient@DEBUG_SUFFIX@.la
getEccSoapClient@DEBUG_SUFFIX@_LDADD += $(MDAQ_UTL_LIBS) $(MDAQ_CLI_LIBS) $(MDAQ_IF_LIBS) $(MFM_LIBS)
getEccSoapClient@DEBUG_SUFFIX@_LDADD += $(ICE_LIBS)
getEccSoapClient@DEBUG_SUFFIX@_LDADD += $(GSOAP_LIBS) -lpthread
getEccSoapClient@DEBUG_SUFFIX@_CPPFLAGS = $(AM_CPPFLAGS) $(VLA_CPPFLAGS)

endif
	
################################################################################
# DAQ targets
################################################################################	
	
#####################################
# getbench-daq
#####################################

if BUILD_ECC
pkgconfig_DATA += getbench-daq@DEBUG_SUFFIX@.pc
getbench_daq_sources =               \
	src/get/daq/DataRouter.cpp       \
	src/get/daq/FrameStorage.cpp    \
	src/get/daq/ByteStorage.cpp 	\
	src/get/daq/FdtDataReceiver.cpp \
	src/get/utl/Time.cpp
if !APPLE
getbench_daq_sources += \
	src/get/daq/FrameCounter.cpp
endif #!APPLE

lib_LTLIBRARIES += libgetbench-daq@DEBUG_SUFFIX@.la
libgetbench_daq@DEBUG_SUFFIX@_la_SOURCES = $(getbench_daq_sources)
libgetbench_daq@DEBUG_SUFFIX@_la_LDFLAGS = $(ICE_LDFLAGS) $(UTILS_LDFLAGS) $(MFM_LDFLAGS) -version-info 1:0:0 -no-undefined

libgetbench_daq@DEBUG_SUFFIX@_la_LIBADD  = $(MDAQ_DAQ_LIBS) $(MDAQ_UTL_LIBS) $(MDAQ_HWSERVER_LIBS) $(MFM_LIBS)
libgetbench_daq@DEBUG_SUFFIX@_la_LIBADD += $(ICE_LIBS) $(BOOST_FILESYSTEM_LIB)

#####################################
# dataRouter
#####################################

bin_PROGRAMS += dataRouter

dataRouter_SOURCES = src/main/dataRouter.cpp

dataRouter_LDADD  = libgetbench-daq@DEBUG_SUFFIX@.la
dataRouter_LDADD += $(MDAQ_ECCSERVER_LIBS) $(MDAQ_UTL_LIBS) $(MDAQ_DAQ_LIBS)
dataRouter_LDADD += $(MFM_LIBS)
dataRouter_LDADD += $(ICE_LIBS)
dataRouter_LDADD += $(BOOST_FILESYSTEM_LIB) $(BOOST_SYSTEM_LIB)
 
dataRouter_LDFLAGS  = $(UTILS_LDFLAGS)
dataRouter_LDFLAGS += $(MFM_LDFLAGS)
dataRouter_LDFLAGS += $(ICE_LDFLAGS)

# Temporary hack for backward compatibility
install-exec-hook:
	cd $(DESTDIR)$(bindir) && \
	$(LN_S) -f dataRouter$(EXEEXT) daqServer$(EXEEXT)

endif #!BUILD_ECC

if BUILD_TARGET_SERVER

if !VXWORKS
if !ARM

#####################################
# getHwServer
#####################################

bin_PROGRAMS += getHwServer 

getHwServer_includedir =  $(includedir)
getHwServer_SOURCES =                      \
	src/GetHwServer.cpp	                   \
	src/SimServerDevices.cpp \
	src/get/cobo/MutantLinkManagerI.cpp src/get/cobo/MutantLinkManagerI.h \
	src/get/utl/PeriodicTask.cpp src/get/utl/PeriodicTask.h \
	src/get/cobo/AsAdPulserMgrI.cpp src/get/cobo/AsAdPulserMgrI.h \
	src/get/cobo/LedManagerI.cpp src/get/cobo/LedManagerI.h \
	src/get/cobo/AsAdAlarmWatcher.cpp src/get/cobo/AsAdAlarmWatcher.h \
	src/get/mt/AlarmServiceI.cpp src/get/mt/AlarmServiceI.h \
	src/get/cobo/CtrlNodeI.cpp src/get/cobo/CtrlNodeI.h

getHwServer_LDFLAGS = $(MDAQ_HWSERVER_LDFLAGS) $(UTILS_LDFLAGS) $(ICE_LDFLAGS)   

getHwServer_LDADD  = libcobo-ctrl-node@DEBUG_SUFFIX@.la
getHwServer_LDADD += libcobo-mutant-link@DEBUG_SUFFIX@.la
getHwServer_LDADD += libmt-alarm@DEBUG_SUFFIX@.la
getHwServer_LDADD += libcobo-led-mgr@DEBUG_SUFFIX@.la
getHwServer_LDADD += libgetbench-eccserver@DEBUG_SUFFIX@.la
getHwServer_LDADD += $(MDAQ_UTL_LIBS) $(MDAQ_DAQ_LIBS) $(MDAQ_HWSERVER_LIBS)
getHwServer_LDADD += $(MFM_LIBS) $(UTILS_LIBS)
getHwServer_LDADD += $(PTHREAD_LIBS)
getHwServer_LDADD += $(ICE_LIBS)
getHwServer_LDADD += $(BOOST_SYSTEM_LIB)

endif #!ARM
endif #!VXWORKS

#####################################
# getHwServer (VxWorks --> includes all the functions, as there is no main in VxWorks)
#####################################

if VXWORKS

bin_PROGRAMS += getHwServer 

getHwServer_SOURCES =                \
	$(ice_dataSender_gen_sources)       \
	src/get/daq/SyncedCircularBuffer.cpp \
	src/get/daq/MemRead.cpp   \
	src/get/daq/DaqCtrlNodeI.cpp \
	src/get/daq/FakeDataSender.cpp \
	src/get/cobo/MutantLinkManagerI.cpp src/get/cobo/MutantLinkManagerI.h \
	src/get/utl/PeriodicTask.cpp src/get/utl/PeriodicTask.h \
	src/get/cobo/AsAdPulserMgrI.cpp src/get/cobo/AsAdPulserMgrI.h \
	src/get/cobo/LedManagerI.cpp src/get/cobo/LedManagerI.h \
	src/get/cobo/AsAdAlarmWatcher.cpp src/get/cobo/AsAdAlarmWatcher.h \
	src/get/mt/AlarmServiceI.cpp src/get/mt/AlarmServiceI.h \
	src/get/cobo/CtrlNodeI.cpp src/get/cobo/CtrlNodeI.h
if HAVE_MFM
getHwServer_SOURCES += src/get/daq/FakeFrameSender.cpp
endif

# HwServer
getHwServer_SOURCES +=               \
	src/GetHwServer.cpp	             \
	src/ServerDevices.cpp

getHwServer_LDFLAGS = $(UTILS_LDFLAGS) $(ICE_LDFLAGS) $(BOOST_LDFLAGS)
getHwServer_LDADD  = libcobo-mutant-link@DEBUG_SUFFIX@.la
getHwServer_LDADD += libcobo-ctrl-node@DEBUG_SUFFIX@.la
getHwServer_LDADD += libmt-alarm@DEBUG_SUFFIX@.la
getHwServer_LDADD += libcobo-led-mgr@DEBUG_SUFFIX@.la
getHwServer_LDADD += $(MFM_LIBS)
getHwServer_LDADD += $(MDAQ_UTL_LIBS) $(MDAQ_DAQ_LIBS) $(MDAQ_HWSERVER_LIBS)
getHwServer_LDADD += $(UTILS_LIBS)
getHwServer_LDADD += $(PTHREAD_LIBS)
getHwServer_LDADD += $(ICE_LIBS)
getHwServer_LDADD += $(BOOST_SYSTEM_LIB)

endif #VXWORKS

if ARM

bin_PROGRAMS += getHwServer 

getHwServer_SOURCES =                \
	$(ice_dataSender_gen_sources)       \
	src/get/daq/SyncedCircularBuffer.cpp \
	src/get/daq/MemRead.cpp   \
	src/get/daq/FakeDataSender.cpp \
	src/get/cobo/MutantLinkManagerI.cpp src/get/cobo/MutantLinkManagerI.h \
	src/get/utl/PeriodicTask.cpp src/get/utl/PeriodicTask.h \
	src/get/cobo/AsAdPulserMgrI.cpp src/get/cobo/AsAdPulserMgrI.h \
	src/get/cobo/LedManagerI.cpp src/get/cobo/LedManagerI.h \
	src/get/cobo/AsAdAlarmWatcher.cpp src/get/cobo/AsAdAlarmWatcher.h \
	src/get/mt/AlarmServiceI.cpp src/get/mt/AlarmServiceI.h \
	src/get/cobo/CtrlNodeI.cpp src/get/cobo/CtrlNodeI.h
if HAVE_MFM
getHwServer_SOURCES += src/get/daq/FakeFrameSender.cpp
endif

# HwServer
getHwServer_SOURCES +=               \
	src/GetHwServer.cpp	             \
	src/ServerDevices.cpp

getHwServer_LDFLAGS = $(UTILS_LDFLAGS) $(ICE_LDFLAGS) $(BOOST_LDFLAGS)
getHwServer_LDADD  = libcobo-mutant-link@DEBUG_SUFFIX@.la
getHwServer_LDADD += libcobo-ctrl-node@DEBUG_SUFFIX@.la
getHwServer_LDADD += libmt-alarm@DEBUG_SUFFIX@.la
getHwServer_LDADD += libcobo-led-mgr@DEBUG_SUFFIX@.la
getHwServer_LDADD += $(MFM_LIBS)
getHwServer_LDADD += $(MDAQ_UTL_LIBS) $(MDAQ_DAQ_LIBS) $(MDAQ_HWSERVER_LIBS)
getHwServer_LDADD += $(UTILS_LIBS)
getHwServer_LDADD += $(PTHREAD_LIBS)
getHwServer_LDADD += $(ICE_LIBS)
getHwServer_LDADD += $(BOOST_SYSTEM_LIB)

endif #ARM


endif #BUILD_TARGET_SERVER

####################################################################################################
# Tests
####################################################################################################

noinst_PROGRAMS =

if BUILD_ECC
noinst_PROGRAMS += get_bench_test
test_files = test/data/mutant-config.xcfg
get_bench_test_SOURCES = \
	test/main.cpp \
	test/get/rc/SubConfigIdSetTest.cpp \
	test/get/rc/ConfigIdTest.cpp \
	test/get/rc/ConfigIdSetTest.cpp \
	test/get/rc/DataLinkTest.cpp \
	test/get/rc/MutantConfigParserTest.cpp \
	test/get/rc/AgetLsbThresholdTest.cpp \
	test/get/rc/AgetPeakingTimeTest.cpp \
	test/get/asad/Lcm1Test.cpp
get_bench_test_LDADD  = libgetbench-cfg@DEBUG_SUFFIX@.la
get_bench_test_LDADD += libcobo-ctrl-node@DEBUG_SUFFIX@.la
get_bench_test_LDADD += libgetbench-eccclient@DEBUG_SUFFIX@.la
get_bench_test_LDADD += libgetbench-eccserver@DEBUG_SUFFIX@.la
get_bench_test_LDADD +=	$(UTILS_LIBS)
get_bench_test_LDADD += $(COMPOUND_CONFIG_LIBS)
get_bench_test_LDADD += $(MDAQ_ECCSERVER_LIBS)
get_bench_test_LDADD += $(MDAQ_UTL_LIBS)
get_bench_test_LDADD += $(BOOST_FILESYSTEM_LIB) $(BOOST_SYSTEM_LIB)
get_bench_test_LDADD += $(ICE_LIBS) $(ICESTORM_LIBS)
endif #BUILD_ECC

if HAVE_MFM
if !APPLE
noinst_PROGRAMS += frameSender
frameSender_SOURCES  = test/get/daq/FrameSender.cpp
frameSender_LDFLAGS = $(ICE_LDFLAGS)
frameSender_LDADD   = $(MDAQ_UTL_LIBS)
frameSender_LDADD  += $(MDAQ_DAQ_LIBS)
frameSender_LDADD  += $(MFM_LIBS)
frameSender_LDADD  += $(ICE_LIBS)
endif
endif

if !VXWORKS
if !ARM
TESTS = get_bench_test
endif #!ARM
endif #!VXWORKS

################################################################################
# General
################################################################################
check-am: check-before
check-before:
	(cd $(srcdir) && rsync -aR $(test_files) $(abs_builddir) && chmod -R u+w $(abs_builddir)/test)

BUILT_SOURCES       = $(ice_dataSender_gen_sources) $(noinst_ice_dataSender_gen_headers)  
BUILT_SOURCES      += $(ice_dataRouter_gen_sources) $(noinst_ice_dataRouter_gen_headers)
BUILT_SOURCES      += $(ice_ecc_gen_sources) $(noinst_ice_ecc_gen_headers)
BUILT_SOURCES      += $(ice_mutant_link_gen_sources) $(ice_mutant_link_gen_headers)
BUILT_SOURCES      += $(ice_cobo_gen_sources) $(ice_cobo_gen_headers)
BUILT_SOURCES 	   += $(ice_mt_alarm_gen_sources) $(ice_mt_alarm_gen_headers)
BUILT_SOURCES 	   += $(ice_mutant_gen_sources) $(ice_mutant_gen_headers)
BUILT_SOURCES 	   += $(soap_eccServer_gen_sources) $(soap_eccServer_gen_extra)
CLEANFILES          = $(BUILT_SOURCES)
DISTCLEANFILES      = $(CLEANFILES) $(test_files)

EXTRA_DIST  = doc/doxyfile
EXTRA_DIST += $(ice_dataSender_sources) $(ice_dataRouter_sources) $(ice_ecc_sources)
EXTRA_DIST += $(ice_cobo_sources) $(ice_mt_alarm_sources) $(ice_mutant_link_sources) $(ice_mutant_sources)
EXTRA_DIST += $(soap_sources) $(test_files)

SUBDIRS = data

dist-hook:
	mkdir -p $(distdir)/doc
	cp -p $(srcdir)/doc/*.dox $(srcdir)/doc/*.png $(srcdir)/doc/*.svg $(distdir)/doc/

include $(top_srcdir)/Makefile.VxWorks.am
include $(top_srcdir)/Makefile.packages.am
