#!/bin/bash

test -e analytic-output || mkdir analytic-output
for fi in '' '-fully_implicit'
  do for tr in '' '-use_trapezoid'
    do for ratio in 1 10
      do for fr in 100 1000 10000
        do ./asynchronous $fi $tr -framerate $fr -substeps_per_frame $ratio 19 | grep errors > analytic-output/out-19$fi$tr-$ratio-$fr.txt
        for p in 1 2
        do ./asynchronous $fi $tr -framerate $fr -substeps_per_frame $ratio -parameter $p 20 | grep errors > analytic-output/out-20-$p$fi$tr-$ratio-$fr.txt
          for d in 0 1
            do for ov in '' '-orthogonal_velocity'
              do  ./asynchronous $fi $tr -framerate $fr -substeps_per_frame $ratio -parameter $p -edge_damping $d $ov 21 | grep errors > analytic-output/out-21-$p-$d$fi$tr-$ratio-$fr$ov.txt
            done
          done
        done
      done
    done
  done
done


