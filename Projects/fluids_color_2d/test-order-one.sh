#!/bin/bash
MASTER="$PHYSBAM/Tools/batch/master"
SLAVE="$PHYSBAM/Tools/batch/slave"

M=0
if [ "x$1" != "x-b" ] ; then
    M=1
    $MASTER &
else
    shift
fi

KP=0
if [ "x$1" == "x-k" ] ; then
    KP=1
    shift
fi

min="$1"
max="$2"
out="$3"
shift 3
L=""
J=""
OO=""
for r in `seq $min $max` ; do
    T=`mktemp`
    O=`mktemp -d`
    echo "$((8*$r))" > $T
    [ $KP = 1 ] && echo OUTPUT $r $O
    K=`$SLAVE -a -o $T -p $r -- "$@" -refine $r -o $O`
    J="$J -d $K"
    OO="$OO $O"
    L="$L $T"
done
PPJ=`$SLAVE $J -p 100 -- /bin/bash ./post-process.sh "$out" $L`
[ $KP = 1 ] || $SLAVE -d $PPJ -p 100 -- /bin/bash -c "rm -r $OO" >/dev/null

if [ $M = 1 ] ; then
    $SLAVE -k
    wait
fi
