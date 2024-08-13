# GET

For installing ROOT on Debian:
1.	sudo apt install dpkg-dev cmake g++ gcc binutils libx11-dev libxpm-dev libxft-dev libxext-dev python libssl-dev python-dev gfortran libpcre3-dev xlibmesa-glu-dev libglew1.5-dev libftgl-dev libmysqlclient-dev libfftw3-dev libcfitsio-dev graphviz-dev libavahi-compat-libdnssd-dev libldap2-dev libxml2-dev libkrb5-dev libgsl0-dev libqt4-dev

For installing GET software on Debian:
1.	sudo apt-get install g++ automake autoconf libtool pkg-config libboost-all-dev libreadline-dev libncurses-dev subversion qtbase5-dev libqt5opengl5 libqwt-qt5-dev libqwtplot3d-qt5-dev liblog4cxx-dev zeroc-ice-all-dev gsoap libmariadb-dev-compat libmariadb-dev libpq-dev libsqlite3-dev alien
2.	Download these three software from below link:
	oracle-instantclient-basic.rpm, oracle-instantclient-sqlplus.rpm, oracle-instantclient-devel.rpm
	https://www.oracle.com/database/technologies/instant-client/linux-x86-64-downloads.html
3.	Then convert the rpm files to deb files using the command below
	sudo alien -d oracle-instantclient-basic.rpm oracle-instantclient-sqlplus.rpm oracle-instantclient-devel.rpm
4.	sudo apt install libaio1
5.	sudo dpkg -i oracle-instantclient-basic.deb oracle-instantclient-sqlplus.deb oracle-instantclient-devel.deb
6.	add this line to .bashrc file: export LD_LIBRARY_PATH=/usr/lib/oracle/11.2/client64/lib:$LD_LIBRARY_PATH
7.	At the end, go to the GET/20210208/GetBench/src/get/soap/getEccSoapClient.cpp at line 174 and replace this line “client.endpoint = endpoint.c_str();” with this one: “client.soap_endpoint_url = endpoint.c_str();”