%define name    compoundconfig
%define version 2.1.2
%define release 1%{?dist}

Name:          %{name}
Summary:       compoundconfig library
Version:       %{version}
Release:       %{release}
Source:        %{name}-%{version}.tar.gz
Group:         Development/Libraries
License:       CeCILL
Requires:      tinyxml = 2.5.2, boost-devel >= 1.33.1

%description
The CompoundConfig library defines API to manipulate tree-like configurations.

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
