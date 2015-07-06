#!/bin/bash
MASTER="$PHYSBAM/Tools/batch/master"
SLAVE="$PHYSBAM/Tools/batch/slave"

max=9
if [ "x$1" = "x-m" ] ; then
    max=$2
    shift 2
fi

$MASTER &

rm -rf new_test_order
mkdir new_test_order
for t in 02 03 04 05 06 10 11 16 17 18 19 24 28 ; do
    for b in d n s ; do
        ./test-order-one.sh -b 2 $max new_test_order/conv-$t-$b.png nice ./fluids_color -resolution 8 -last_frame 1 -s 1.3 -m .8 -kg 1.2 -bc_$b -dt .05 $t
    done
done
for t in 00 01 08 09 12 13 14 20 21  250 251 252 253 254 255 256 257 258 259  ; do
    ./test-order-one.sh -b 2 $max new_test_order/conv-$t.png nice ./fluids_color -resolution 8 -last_frame 1 -s 1.3 -m .8 -kg 1.2 -dt .05 $t
done
for t in 250 251 252 253 254 255 256 257 258 259 260 ; do
    ./test-order-one.sh -b 2 $max new_test_order/conv-$t.png nice ./fluids_color -resolution 8 -last_frame 1 -s 1.3 -m .8 -kg 1.2 -dt .05 $t
done
U=()
U[0]="u=1;v=1;"
U[1]="u=t+1;v=t+1;"
U[2]="u=exp(t)+1;v=exp(-t)+1;"
U[3]="pi=3.14159265358979;x=sin(x*2*pi);y=cos(y*2*pi);u=-y;v=x;"
U[4]="pi=3.14159265358979;x=sin(x*2*pi);y=cos(y*2*pi);u=-y*exp(-t)+exp(t);v=x*exp(-t)+exp(t);"
U[5]="pi=3.14159265358979;u=sin(x*2*pi)*sin(y*2*pi);v=cos(x*2*pi)*cos(y*2*pi);"
U[6]="pi=3.14159265358979;T=exp(t)+tan(t+1);u=sin(x*2*pi)*sin(y*2*pi)*T;v=cos(x*2*pi)*cos(y*2*pi)*T;"

P=()
P[0]="p=0;"

S=()
S[0]="s00=1;s10=0;s11=1;"
S[1]="pi=3.14159265358979;x=sin(x*2*pi);y=cos(y*2*pi);s00=x;s10=0;s11=1;"
S[2]="pi=3.14159265358979;x=sin(x*2*pi);y=cos(y*2*pi);s00=x;s10=0;s11=x;"
S[3]="s00=t;s10=0;s11=t;"
S[4]="s00=1;s10=0;s11=t;"
S[5]="pi=3.14159265358979;x=sin(x*2*pi);y=cos(y*2*pi);s00=x*x;s10=x*y;s11=y*y;"
S[6]="pi=3.14159265358979;x=sin(x*2*pi);y=cos(y*2*pi);s00=x*x*exp(t);s10=x*y*exp(2*t);s11=y*y;"
for u in `seq 0 $((${#U[@]}-1))` ; do
    for p in `seq 0 $((${#P[@]}-1))` ; do
        for s in `seq 0 $((${#S[@]}-1))` ; do
            ./test-order-one.sh -b 2 $max new_test_order/conv-custom-$u-$p-$s.png nice ./fluids_color -resolution 8 -last_frame 1 -s 1.3 -m .8 -kg 1.2 -dt .05 250 -set_u0 "${U[$u]}" -set_S0 "${S[$s]}" -set_p0 "${P[$p]}"
        done
    done
done

$SLAVE -k
wait
