#!/bin/bash
# $Id$
#
# Example:
# SRC=/tmp/src PREFIX=/usr/local source ./scripts/install_get_components.sh

#set -x

# Default arguments
BUILD_CONFIG="${pentium-linux-release}"
PREFIX=${PREFIX:-/usr/local)}
SRC=${SRC:-$PWD}
INSTALL_CMD="sudo -u root"

# install function
function install_component()
{
	component=$1
	prefix=$2
	shift
	shift
	src_dir=${SRC}/${component}
	options="$@"

	echo "===== Installing '${component}' into '${prefix}' with options '${options}' ..." >&2
	top_dir=${PWD}
	cd $src_dir && \
	./autogen.sh  > /dev/null && \
	rm -rf build/${BUILD_CONFIG} && \
	mkdir -p build/${BUILD_CONFIG} && \
	cd build/${BUILD_CONFIG} && \
	cmd="${src_dir}/configure --prefix=$prefix ${options} > /dev/null" && echo $cmd  && eval $cmd && \
	make  > /dev/null && \
	${INSTALL_CMD} make install  > /dev/null && \
	cd ${top_dir} && \
	echo "===== Successfully installed package '${component}'" >&2
}

# Searches recursively for existing parent directory  
function get_existing_parent()
{
	dir=$1
	echo "Searching for existing parent directory of $dir" >&2
	while test -n "$dir"
	do
		if [[ $dir == "." ]] || [[ $dir == "/" ]]; then
			break
		elif test -d $dir; then
			break
		fi
		dir=`dirname $dir`
	done
	echo $dir
}

# Installation directory exists but user does not have write permission
if test -d "$PREFIX" && ! test -w "$PREFIX"; then
	echo "Installation in $PREFIX requires administrative privileges" >&2
	sudo echo "Caching sudoer rights"
	INSTALL_CMD="sudo -u root"
# Installation directory does not exist
elif ! test -d "$PREFIX"; then
	topdir=$(get_existing_parent $PREFIX)
	# User does not have write permission
	if ! test -w "$topdir"; then
		echo "Installation in $PREFIX requires administrative privileges" >&2
		sudo mkdir -p $topdir
		INSTALL_CMD="sudo -u root"
	else
		INSTALL_CMD=""
	fi
	${INSTALL_CMD} mkdir -p $topdir
# Installation directory exists and user has write permission
else
	INSTALL_CMD=""
fi

# Confirm
read -p  "Ready to install all components into '$PREFIX' (yes/no)? " ans
if test "$ans" != "yes" ; then
	echo "Aborting"
	return 1
fi

export PKG_CONFIG_PATH=${PKG_CONFIG_PATH:+:}"${PREFIX}/lib/pkgconfig"
echo "PKG_CONFIG_PATH=$PKG_CONFIG_PATH"

# Install
install_component TinyXml ${PREFIX} && \
install_component CompoundConfig ${PREFIX} && \
install_component Utilities ${PREFIX}  && \
install_component ConfigWidgets ${PREFIX} && \
#install_component soci ${PREFIX} --enable-backend-mysql && \
#install_component ConfigDatabase ${PREFIX} --with-soci=${PREFIX} && \
#install_component DatabaseWidgets ${PREFIX} && \
# ConfigServer requires swig and is used by the Ganil configuration GUI
# install_component ConfigServer ${PREFIX} && \
install_component MDaq ${PREFIX} && \
install_component StateMachine ${PREFIX} --enable-distributed-sm-libs=none --with-log-level=info && \
install_component MultiFrame ${PREFIX} && \
install_component GetBench ${PREFIX}  && \
install_component GetController ${PREFIX} #&& \
install_component CoBoFrameViewer ${PREFIX}
# ConfigWizard requires QtRoot and on many distributions root-config says qt feature is installed while they are not
# So we shall let users decide to install ConfigWizard only when they decide to
#install_component ConfigWizard ${PREFIX}
