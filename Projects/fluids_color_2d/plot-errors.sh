#!/bin/bash

echo "$@" >> args.txt

base_name="$1"
shift

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

echo T >> args.txt
CMD=""
UFILE=""
PFILE=""
if `echo "$base_name" | grep -q '\.eps$'` ; then
    CMD="set terminal postscript eps color"
    UFILE="${base_name/.eps/-u.eps}"
    PFILE="${base_name/.eps/-p.eps}"
elif `echo $base_name | grep -q '\.png$'` ; then
    CMD="set terminal png"
    UFILE="${base_name/.png/-u.png}"
    PFILE="${base_name/.png/-p.png}"
else
    echo "unrecognized image extension: $base_name" >> args.txt
    exit 1
fi
echo $UFILE $PFILE >> args.txt
gnuplot -p -e "$CMD ; set output '$UFILE' ; plot '$M' u (log10(\$1)):(log10(\$3)) title 'L-inf error' , '$M' u (log10(\$1)):(log10(\$4)) title 'L-2 error' , -2*x , -2*x-1 , -x-2;"
gnuplot -p -e "$CMD ; set output '$PFILE' ; plot '$M' u (log10(\$1)):(log10(\$8)) title 'L-inf error' , '$M' u (log10(\$1)):(log10(\$9)) title 'L-2 error' , -2*x , -2*x-1 , -x-2;"
rm $M
