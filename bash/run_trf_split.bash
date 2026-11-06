#!/bin/bash -e

i=$1

if [ -z "$i" ] ; then
	echo "Missing split index" >> /dev/stderr
	exit 1
fi

cd ../big_splits

# -rw-r--r-- 1 robedgar ablab 281M Nov  5 23:34 all.fa.2.7.7.80.10.50.500.dat
# -rw-r--r-- 1 robedgar ablab 1.0G Nov  5 23:34 all.fa.2.7.7.80.10.50.500.mask

../bin/trf409.linux64 split$i.fa 2 7 7 80 10 50 500 -m -h \
  > /dev/null \
  2> /dev/null
