%define name    tinyxml
%define version 2.5.2
%define release 1%{?dist}

Name:          %{name}
Summary:       A simple, small, C++ XML parser
Version:       %{version}
Release:       %{release}
Source:        %{name}-%{version}.tar.gz
URL:           http://www.grinninglizard.com/tinyxml/
License:       zlib
Group:         Development/Libraries

%description
TinyXML is a simple, small, C++ XML parser that can be easily integrating
into other programs. Have you ever found yourself writing a text file parser
every time you needed to save human readable data or serialize objects?
TinyXML solves the text I/O file once and for all.
(Or, as a friend said, ends the Just Another Text File Parser problem.)

%prep
%setup

%build
%configure
make CXXFLAGS="-g0 -O3"

%install
rm -rf %{buildroot}
%makeinstall

%clean
rm -rf %{buildroot}

%files
%defattr(-,root,root)
%{_libdir}/
%{_includedir}/
