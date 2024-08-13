#! /bin/sh
svn status --no-ignore | grep '^I' |  sed "s/^[I] //" | xargs -I{} rm -rf "{}"
