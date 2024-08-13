%define name    configwidgets
%define version 2.0.3
%define release 1%{?dist}

Name:          %{name}
Summary:       configwidgets library
Version:       %{version}
Release:       %{release}
Source:        %{name}-%{version}.tar.gz
Group:         Development/Libraries
License:       CeCILL
Requires:      compoundconfig >= 2.1.0, utilities >= 1.1.0, qt >= 4.3.0
BuildRequires: qt-devel >= 4.3.0

%description
This is a Configuration Widget used with the CompoundConfig Library.

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
%doc README
%{_bindir}/
%{_includedir}/
%{_libdir}/
