include $(top_srcdir)/am/ice.am
include $(top_srcdir)/am/doxygen.am

ACLOCAL_AMFLAGS = -I m4
AUTOMAKE_OPTIONS = subdir-objects

dist_noinst_SCRIPTS = autogen.sh cleanup.sh

AM_CPPFLAGS = -I@builddir@/src/CCfg 		\
	$(COMPOUND_CONFIG_CFLAGS)				\
	$(ICE_CPPFLAGS) $(ICE_DEFINE)			\
	$(BOOST_CPPFLAGS)						\
	$(JNI_CPPFLAGS)
# See http://www.swig.org/Doc2.0/Java.html#Java_compiling_dynamic :
# If you are going to use optimisations turned on with gcc, ensure you also compile with -fno-strict-aliasing.
# The GCC optimisations have become more aggressive from gcc-4.0 onwards and will result in code that fails with strict aliasing optimisations turned on.
AM_CXXFLAGS = $(PTHREAD_CFLAGS) -fno-strict-aliasing
	
ice_sources = src/CCfg/ConfigServer.ice
ice_gen_sources =
ice_gen_headers =
if ICE
ice_gen_sources += src/CCfg/ConfigServer.cpp
ice_gen_headers += src/CCfg/ConfigServer.h
endif

jar_file = @top_builddir@/jccfg.jar
JAVAROOT = @top_builddir@/src/

java_gen_sources =								\
	src/CCfg/CCfg.java							\
	src/CCfg/CCfgJNI.java						\
	src/CCfg/SWIGTYPE_p_CCfg__CConfig.java 		\
	src/CCfg/CConfigFacade.java					\
	src/CCfg/DocumentFacade.java				\
	src/CCfg/StringVector.java

# classes are defined explicitly, and we omit the src/ for
# the first argument to cope with a bug in jar -C
class_files2 = 									\
	CCfg/CCfg.class								\
	src/CCfg/CCfgJNI.class						\
	src/CCfg/SWIGTYPE_p_CCfg__CConfig.class 	\
	src/CCfg/CConfigFacade.class				\
	src/CCfg/DocumentFacade.class				\
	src/CCfg/StringVector.class

nodist_noinst_JAVA = $(java_sources) $(java_gen_sources)
swig_gen_cpp_sources = src/CCfg/SwigFacade.cxx

all_java_sources = $(java_sources) $(java_gen_sources)
class_files = ${all_java_sources:%.java=%.class}

BUILT_SOURCES =			\
	$(ice_gen_sources)	\
	$(ice_gen_headers)

ccfg_includedir = $(includedir)
ccfg_include_HEADERS =			\
	src/CCfg/CConfigFacade.h	\
	src/CCfg/ConfigServerI.h	\
	src/CCfg/DocumentFacade.h	\
	src/CCfg/RCConfigI.h		\
	src/CCfg/RDocumentI.h

noinst_LTLIBRARIES = libccfgfacade.la
libccfgfacade_la_SOURCES = 			\
	src/CCfg/CConfigFacade.cpp		\
	src/CCfg/DocumentFacade.cpp

lib_LTLIBRARIES = libccfg.la
libccfg_la_SOURCES = 
nodist_libccfg_la_SOURCES = $(swig_gen_cpp_sources)
libccfg_la_LDFLAGS = -no-undefined -version-info 0:0:0
libccfg_la_LIBADD = libccfgfacade.la $(COMPOUND_CONFIG_LIBS)

bin_PROGRAMS = 
if ICE
bin_PROGRAMS += client server testClient
client_SOURCES =		\
	$(ice_gen_headers)	\
	$(ice_gen_sources)	\
	src/Client.cpp
client_LDFLAGS = $(ICE_LDFLAGS)
client_LDADD = $(ICE_LIBS)

server_SOURCES =				\
	$(ice_gen_headers)			\
	$(ice_gen_sources)			\
	src/CCfg/ConfigServerI.cpp	\
	src/CCfg/RCConfigI.cpp		\
	src/CCfg/RDocumentI.cpp		\
	src/Server.cpp
server_LDFLAGS = $(ICE_LDFLAGS)
server_LDADD = libccfgfacade.la $(COMPOUND_CONFIG_LIBS) $(ICE_LIBS)
endif

testClient_SOURCES =			\
	$(ice_gen_headers)			\
	$(ice_gen_sources)			\
	test/RDocumentTest.cpp		\
	test/ConfigServerTest.cpp	\
	test/RCConfigTest.cpp 		\
	test/main.cpp
testClient_LDFLAGS = $(ICE_LDFLAGS)
testClient_LDADD = $(ICE_LIBS)
EXTRA_DIST =		\
	doc/doxyfile	\
	$(ice_sources)	\
	$(swig_sources)

CLEANFILES = $(BUILT_SOURCES)	\
	$(swig_gen_cpp_sources)		\
	$(java_gen_sources)			\
	$(jar_file)					\
	${class_files}


###############
# Rules added #
###############

${jar_file}: ${class_files}
	jar vcf $@ -C src ${class_files2}	

define SWIG_RUN
	$(SWIG) -Wextra -java -c++ -package CCfg -o $(swig_gen_cpp_sources) @top_srcdir@/src/CCfg/SwigFacade.i
endef

src/CCfg/SwigFacade.cxx: @top_srcdir@/src/CCfg/SwigFacade.i
	$(SWIG_RUN)


############
# PKG INIT #
############
#Init Inno Setup Compiler
WIN_BUILD_DIR=$(PWD)/_pkgbuild
WIN_INSTALL_DIR=$(PWD)/_pkginst

#The Needed Files for the Packages Creation
DEB_FOLDER=$(srcdir)/packages/debian
WIN_FOLDER=$(srcdir)/packages/win/$(PACKAGE_TARNAME).iss

#The location of the Packages after their creation
DEB_CREATE_DIR=$(builddir)/pkg
WIN_CREATE_DIR=$(builddir)/pkg

#Extra files (not present in archive with make dist) needed to create the Deb file
EXTRA_FILES= 							\
$(srcdir)/src/CCfg/SwigFacade.i 		\
$(srcdir)/src/CCfg/CatchExceptions.h	

#Add the jar to the packages - Comment the Jar path depending on platform.
inst_DATA = ${jar_file}
#Jar path for Linux Packages
#instdir = /usr/lib
#Jar path for Windows Package
instdir = $(top_builddir)/../_pkginst/lib

#######
# DEB #
#######
deb:
	mkdir -p $(DEB_CREATE_DIR);
	$(MAKE) dist;
	tar -xvzf $(PACKAGE_TARNAME)-$(PACKAGE_VERSION).tar.gz;
	cp -R $(DEB_FOLDER) $(PACKAGE_TARNAME)-$(PACKAGE_VERSION)/;
	cp $(EXTRA_FILES) $(builddir)/$(PACKAGE_TARNAME)-$(PACKAGE_VERSION)/src/CCfg/;
	cd $(PACKAGE_TARNAME)-$(PACKAGE_VERSION); 									\
	chmod 755 debian/rules; 													\
	dpkg-buildpackage -rfakeroot -us -uc -B;
	mv $(PACKAGE_TARNAME)_$(PACKAGE_VERSION)-*.deb 								\
	$(DEB_CREATE_DIR)/$(PACKAGE_TARNAME)_$(PACKAGE_VERSION)_$(host_cpu).deb;
	rm -rf $(PACKAGE_TARNAME)-$(PACKAGE_VERSION).tar.gz 						\
	$(PACKAGE_TARNAME)-$(PACKAGE_VERSION) $(PACKAGE_TARNAME)*.changes


#######
# RPM # - We use Alien to manage easily the add of some files (jccfg, etc) during the Rpm creation 
#######
rpm:
	mkdir -p $(DEB_CREATE_DIR);
	$(MAKE) dist;
	tar -xvzf $(PACKAGE_TARNAME)-$(PACKAGE_VERSION).tar.gz;
	cp -R $(DEB_FOLDER) $(PACKAGE_TARNAME)-$(PACKAGE_VERSION)/;
	cp $(EXTRA_FILES) $(builddir)/$(PACKAGE_TARNAME)-$(PACKAGE_VERSION)/src/CCfg/;
	cd $(PACKAGE_TARNAME)-$(PACKAGE_VERSION); 								\
	chmod 755 debian/rules; 												\
	dpkg-buildpackage -rfakeroot -us -uc -B;
	mv $(PACKAGE_TARNAME)_$(PACKAGE_VERSION)-*.deb 							\
	$(DEB_CREATE_DIR)/$(PACKAGE_TARNAME)_$(PACKAGE_VERSION)_$(host_cpu)_rpm.deb;
	cd $(DEB_CREATE_DIR)/;										\
	fakeroot alien -r -k $(PACKAGE_TARNAME)_$(PACKAGE_VERSION)_$(host_cpu)_rpm.deb;
	mv $(DEB_CREATE_DIR)/$(PACKAGE_TARNAME)-$(PACKAGE_VERSION)-*.rpm 				\
	$(DEB_CREATE_DIR)/$(PACKAGE_TARNAME)-$(PACKAGE_VERSION)_$(host_cpu).rpm;
	rm -rf $(PACKAGE_TARNAME)-$(PACKAGE_VERSION).tar.gz 						\
	$(PACKAGE_TARNAME)-$(PACKAGE_VERSION) $(PACKAGE_TARNAME)*.changes			\
	$(DEB_CREATE_DIR)/$(PACKAGE_TARNAME)_$(PACKAGE_VERSION)_$(host_cpu)_rpm.deb;


# Windows # - Cross compile and create the Windows installer using InnoSetupCompiler
win:
	mkdir -p $(WIN_BUILD_DIR)
	mkdir -p $(WIN_INSTALL_DIR)
	cd $(WIN_BUILD_DIR); cp $(EXTRA_FILES) .; \
	$(top_srcdir)/configure --prefix=$(WIN_INSTALL_DIR) $(DISTCHECK_CONFIGURE_FLAGS); \
	make install
	cp $(WIN_FOLDER) .
	$(SETUP_EXE) /Q /O$(WIN_CREATE_DIR) $(PACKAGE_TARNAME).iss
	mv $(WIN_CREATE_DIR)/setup.exe $(WIN_CREATE_DIR)/$(PACKAGE_TARNAME)-$(PACKAGE_VERSION)_$(host_cpu).exe
	rm -rf $(PACKAGE_TARNAME).iss $(WIN_INSTALL_DIR) $(WIN_BUILD_DIR)
	