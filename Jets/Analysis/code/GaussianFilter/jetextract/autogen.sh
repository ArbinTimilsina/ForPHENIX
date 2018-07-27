#!/bin/sh
srcdir=`dirname $0`
test -z "$srcdir" && srcdir=.

(cd $srcdir; aclocal -I ${OFFLINE_MAIN}/share;\
libtoolize --copy --force;
autoheader;
automake --copy --add-missing;
autoconf)

$srcdir/configure "$@"

