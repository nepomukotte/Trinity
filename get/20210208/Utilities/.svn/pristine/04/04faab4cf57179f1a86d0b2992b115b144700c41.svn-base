%define name    utilities
%define version 1.3.0
%define release 1%{?dist}

Name:          %{name}
Summary:       Mordicus utilities library
Version:       %{version}
Release:       %{release}
Source:        %{name}-%{version}.tar.gz
Group:         Development/Libraries
License:       CeCILL
Requires:      boost-devel >= 1.33.1

%description
This library contains general-purpose utility classes for C++ application development.

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
