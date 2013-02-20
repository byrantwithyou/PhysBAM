#!/bin/bash
MASTER="$PHYSBAM/Tools/batch/master"
SLAVE="$PHYSBAM/Tools/batch/slave"

$MASTER &

function emit_test()
{
    out="$1"
    shift
    L=""
    J=""
    OO=""
    for r in {2..9} ; do
        T=`mktemp`
        O=`mktemp -d`
        echo "$((8*$r))" > $T
        K=`$SLAVE -a -o $T -p $r -- nice ./fluids_color_2d -resolution 8 -last_frame 1 -refine $r -s 1.3 -m .8 -kg 1.2 -o $O "$@"`
        J="$J -d $K"
        OO="$OO $O"
        L="$L $T"
    done
    PPJ=`$SLAVE $J -p 10 -- /bin/bash ./post-process.sh "$out" $L`
    $SLAVE -d $PPJ -p 10 -- /bin/bash -c "echo rm -r $OO >> rm-list.txt" >/dev/null
}

rm -rf new_test_order
mkdir new_test_order
for t in 02 03 04 05 06 10 11 16 17 18 19 24 ; do
    for b in d n s ; do
        emit_test  new_test_order/conv-$t-$b.png -bc_$b -dt .05 $t
    done
done
for t in 00 01 07 08 09 12 13 14 20 21 ; do
    emit_test new_test_order/conv-$t-x.png -dt .05 $t
done

$SLAVE -k
wait
