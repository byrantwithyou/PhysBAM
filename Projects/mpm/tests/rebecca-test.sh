#!/bin/bash

CORES=`cat /proc/cpuinfo  | grep processor | wc -l`
TF=('-affine' '-no_affine' '-no_affine -flip 1')
TN=('-apic' '-pic' '-flip')

FF=('' '-use_exp_F')
FN=('' '-F')

EF=('' '-symplectic_euler')
EN=('-imp' '-exp')

CF=('' '-I 1')
CN=('' '-tait')

for K in 1 10 100 1000 ; do
    for dt in 0.03 0.01 0.003 0.001 0.0003 ; do
        for o in 2 3 ; do
            for T in 0 1 2 ; do
                for F in 0 1 ; do
                    for E in 0 1 ; do
                        for C in 0 1 ; do
                            name="block${EN[$E]}${TN[$T]}${FN[$F]}${CN[$C]}-$K-$o-$dt"
                            echo "../mpm 71 -resolution 100 -cfl 0.2 -scale_E $K ${FF[$F]} ${TF[$T]} ${CF[$C]} -max_dt $dt -order $o -framerate 20 -last_frame 100 ${EF[$E]} -o $name >& $name.txt"
                        done
                    done
                done
            done
        done
    done
done | shuf | xargs -P $CORES -n 1 -d '\n' bash -c
