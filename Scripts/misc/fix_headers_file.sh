#!/bin/bash

DIR=`dirname $0`
FILE=/tmp/headers-`perl -e 'print rand();'`

( ls *.h ; cd $PUBLIC ; find . -name '*.h' ) | sed 's@\./@@' > $FILE;
echo $* | xargs perl $DIR/fix_headers.pl $FILE 2>&1

rm $FILE
