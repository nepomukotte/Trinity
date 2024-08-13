#!/bin/bash
# $Id$
#
# Examples:
#
# RELEASE=20131202 DEST=/tmp/src source ./scripts/download_get_components.sh
# TAG=HEAD DEST=/tmp/src source ./scripts/download_components.sh
#
# Default arguments
TAG=${TAG:-release-${RELEASE:-$(date +%Y%m%d)}}
DEBUG=${DEBUG:-no}
DEST=${DEST:-$PWD}
SVN_USER=${SVN_USER:-}
SVN_PASSWORD=${SVN_PASSWORD:-}

# Function to download a software tag
function download_component()
{
	repo=$1
	component=$2
	tag=$3
	subdir=$4
	
	repourl="https://dsm-trac.cea.fr/svn/"${repo}
	
	src_url=${repourl}/${component}/tags/${component}${subdir:+-${subdir}}-${tag}
	if test "x$tag" = "xHEAD"; then
		src_url=${repourl}/${component}/trunk${subdir:+/${subdir}}
	fi
	dest_url=${DEST}/${component}

	SVN_CRED=""
	if test "x$SVN_USER" != "x" ; then SVN_CRED="--username $SVN_USER"; fi
	if test "x$SVN_PASSWORD" != "x"; then SVN_CRED="$SVN_CRED --password $SVN_PASSWORD"; fi

	cmd="svn checkout -q ${SVN_CRED} ${src_url} ${dest_url}"
	if test "$DEBUG" != "no"; then
		echo $cmd
	else
		mkdir -p ${DEST} && rm -rf ${dest_url} && echo "Downloading ${src_url} into ${dest_url}" && eval $cmd
	fi
}

# Confirm
read -p  "Ready to download tag '${TAG}' of all components into ${DEST}? (yes/no) " ans
if test "$ans" != "yes" ; then
	echo "Aborting"
	exit 1
fi

# Download all packages
download_component cconfig TinyXml $TAG && \
download_component cconfig CompoundConfig $TAG && \
download_component cconfig ConfigWidgets $TAG && \
download_component cconfig soci $TAG && \
download_component cconfig ConfigDatabase $TAG && \
download_component cconfig DatabaseWidgets $TAG && \
download_component cconfig ConfigServer $TAG && \
download_component mordicus Utilities $TAG && \
download_component mordicus MDaq $TAG && \
download_component mordicus StateMachine $TAG && \
download_component mordicus MultiFrame $TAG && \
download_component get GetBench $TAG software && \
download_component get ConfigWizard $TAG && \
download_component get GetController $TAG && \
download_component get CoBoFrameViewer $TAG
