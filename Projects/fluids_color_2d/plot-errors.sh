#!/bin/bash

base_name="$1"
type="$2"

shift 2

M=`mktemp`

(
    for d in "$@" ; do
        cmd=`grep command "$d/common/log.txt"`
        base_res=`echo "$cmd" | sed 's/.*-resolution *\([0-9][0-9]*\) .*/\1/';`
        refine=`echo "$cmd" | sed 's/.*-refine *\([0-9][0-9]*\) .*/\1/';`
        echo -n $((base_res*refine))
        grep max_error "$d/common/log.txt" | tail -n 1
    done
) | sort -n > $M

CMD=""
if [ "x$type" = "xeps" ] ; then
    CMD="set terminal postscript eps color"
elif [ "x$type" = "xpng" ] ; then
    CMD="set terminal png"
else
    echo "unrecognized image type: $type"
    exit 1
fi

gnuplot -p -e "$CMD ; set output '$base_name-u.$type' ; plot '$M' u (log10(\$1)):(log10(\$3)) title 'L-inf error' , '$M' u (log10(\$1)):(log10(\$4)) title 'L-2 error' , -2*x , -2*x-1 , -x-2;"
gnuplot -p -e "$CMD ; set output '$base_name-p.$type' ; plot '$M' u (log10(\$1)):(log10(\$8)) title 'L-inf error' , '$M' u (log10(\$1)):(log10(\$9)) title 'L-2 error' , -2*x , -2*x-1 , -x-2;"
rm $M
