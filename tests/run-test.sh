#!/bin/sh

topdir="$(readlink -f ..)"
for f in `find . -name '*.rb'`; do
	if test "$f" = "./sf/make_test.rb"; then continue; fi
	if test "$f" = "./speed/bench.rb"; then continue; fi
	if test "$f" = "./tensor.rb"; then continue; fi
	pushd `dirname $f` > /dev/null
	ruby -w -I "$topdir/lib" -I "$topdir/ext" `basename $f`
#	ruby -w `basename $f`
	result=$?
	popd > /dev/null
	if test "$result" != "0"; then
		echo $f failed
		exit $result
	fi
done
