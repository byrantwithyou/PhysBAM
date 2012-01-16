#!/bin/bash

# for m in '' -use_corotated -use_constant_ife ; do
#     for e in 1 "7 -dampen 10 -last_frame 500" 8 16 "17 -dampen 10 -last_frame 1000" "18 -dampen 10 -last_frame 500" 24 25 26 27 ; do
#         nice extrap_3d $e $m -o out$m-`echo $e | sed 's/ .*//'` >&/dev/null &
#     done
# done

for m in "$@" ; do
    for e in 1 "7 -dampen 10 -last_frame 500" 8 16 "17 -dampen 10 -last_frame 1000" "18 -dampen 10 -last_frame 500" 24 25 26 27 ; do
        nice extrap_3d $e $m -o "out`echo $m | sed 's/ .*//'`-`echo $e | sed 's/ .*//'`" >&/dev/null &
    done
done

