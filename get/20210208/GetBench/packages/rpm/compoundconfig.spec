%define name    compoundconfig
%define version 2.1.0 
%define release noarch

Name:          %{name} 
Summary:       compoundconfig library
Version:       %{version}
Release:       %{release}
Source:        %{name}-%{version}.tar.gz 
URL:           http://www.cea.fr/ 
Group:         Libraries
BuildRoot:     %{_tmppath}/%{name}-buildroot
License:       lherve CEA
Requires:      tinyxml >= 2.5.2

%description
The CompoundConfig library defines API to manipulate tree-like configurations.

%prep 
%setup

%build 
./configure --prefix=/usr
make

%install
rm -rf $RPM_BUILD_ROOT
%makeinstall

%clean 
rm -rf $RPM_BUILD_ROOT 

%files
%defattr(-,root,root,0755) 
%{_includedir}/Utils/
%{_includedir}/Utils/Math/
%{_includedir}/Utils/Math/float.h
%{_includedir}/Utils/Io/
%{_includedir}/Utils/Io/PointerValueStreamer.hpp
%{_includedir}/Utils/String/
%{_includedir}/Utils/String/Comparison.hpp
%{_includedir}/CCfg/
%{_includedir}/CCfg/Visitor.h
%{_includedir}/CCfg/StructureNode.h
%{_includedir}/CCfg/Range.h
%{_includedir}/CCfg/Scalar.h
%{_includedir}/CCfg/ClassInfo.h
%{_includedir}/CCfg/Value.h
%{_includedir}/CCfg/CConfig.h
%{_includedir}/CCfg/Interval.h
%{_includedir}/CCfg/Enum.h
%{_includedir}/CCfg/Exception.h
%{_includedir}/CCfg/View/
%{_includedir}/CCfg/View/Search.hpp
%{_includedir}/CCfg/View/Result.h
%{_includedir}/CCfg/View/InheritersIterator.h
%{_includedir}/CCfg/View/InheritanceIterator.h
%{_includedir}/CCfg/View/InheritedIterator.h
%{_includedir}/CCfg/View/Search.h
%{_includedir}/CCfg/View/ClassIterator.h
%{_includedir}/CCfg/Attribute.h
%{_includedir}/CCfg/Path.h
%{_includedir}/CCfg/ScalarClassInfo.h
%{_includedir}/CCfg/ValidationVisitor.h
%{_includedir}/CCfg/UnsupportedOperation.h
%{_includedir}/CCfg/PathNode.h
%{_includedir}/CCfg/Io/
%{_includedir}/CCfg/Io/ConfigToTextVisitor.h
%{_includedir}/CCfg/Io/AbstractStorageFactory.h
%{_includedir}/CCfg/Io/XmlToConfigVisitor.h
%{_includedir}/CCfg/Io/UnhandledUrl.h
%{_includedir}/CCfg/Io/Storage.h
%{_includedir}/CCfg/Io/FormatFactory.h
%{_includedir}/CCfg/Io/BufferStorage.h
%{_includedir}/CCfg/Io/AbstractFormatFactory.h
%{_includedir}/CCfg/Io/Format.h
%{_includedir}/CCfg/Io/TextFormat.h
%{_includedir}/CCfg/Io/Document.h
%{_includedir}/CCfg/Io/XmlFormatFactory.h
%{_includedir}/CCfg/Io/IncludeHandler.h
%{_includedir}/CCfg/Io/FileStorage.h
%{_includedir}/CCfg/Io/XmlFormat.h
%{_includedir}/CCfg/Io/StorageFactory.h
%{_includedir}/CCfg/Io/ConfigToXmlVisitor.h
%{_includedir}/CCfg/Io/FileStorageFactory.h
%{_includedir}/CCfg/Io/StreamStorage.h
%{_includedir}/CCfg/Io/Url.h
%{_includedir}/CCfg/Comment.h
%{_includedir}/CCfg/Object.h
%{_includedir}/CCfg/TypeInfo.h
%{_libdir}/pkgconfig/
%{_libdir}/pkgconfig/CompoundConfig.pc
%{_libdir}/libCompoundConfig.a
%{_libdir}/libCompoundConfig.so.1.0.0
%{_libdir}/libCompoundConfig.la
%{_libdir}/libCompoundConfig.so
%{_libdir}/libCompoundConfig.so.1
